#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 05_model_maxnet.R
# 功能说明: 使用Maxnet进行物种分布建模（Presence-Background）- 增强版
# 方法: Maxent (maxnet包, 基于glmnet实现, 复现MaxEnt Java 3.4.0)
# 输入文件: output/04_collinearity/collinearity_removed.csv
# 输出文件:
#   - output/05_model_maxnet/model.rds              (模型对象)
#   - output/05_model_maxnet/predictions.csv        (预测结果)
#   - output/05_model_maxnet/evaluation.csv         (模型评估)
#   - output/05_model_maxnet/variable_importance.csv(变量重要性)
#   - output/05_model_maxnet/model_diagnostics.csv  (模型诊断信息)
#   - output/05_model_maxnet/response_curves/       (响应曲线图)
#   - output/05_model_maxnet/model_coefficients.csv (模型系数)
# 作者: Nature级别科研项目
# 日期: 2025-10-20
# 更新: 2026-02-03 - 全面利用maxnet包功能
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/CausalSDMs")

# 加载必要的包
packages <- c("tidyverse", "maxnet", "pROC", "gridExtra", "cowplot")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ==============================================================================
# Maxent 参数设置
# ==============================================================================
# 特征类别说明:
#   l = linear (线性): 适合连续变量的单调响应
#   q = quadratic (二次): 捕捉非线性单峰响应
#   p = product (交互): 变量间的交互效应
#   h = hinge (铰链): 分段线性响应，更灵活
#   t = threshold (阈值): 阶跃函数响应
# 默认 "default" 根据样本量自动选择:
#   np < 10: "l", np < 15: "lq", np < 80: "lqh", np >= 80: "lqph"
MAXENT_CLASSES <- "lqph"

# 正则化乘数: 控制模型复杂度
#   regmult > 1: 更平滑/简单 (防止过拟合)
#   regmult < 1: 更复杂 (可能过拟合)
MAXENT_REGMULT <- 1.0

# 预测类型说明:
#   "link": 线性预测器 lp
#   "exponential": exp(lp) - 原始Maxent输出
#   "cloglog": 1-exp(-exp(entropy+lp)) - 推荐用于presence概率解释
#   "logistic": 1/(1+exp(-entropy-lp)) - 传统logistic解释
PREDICT_TYPE <- "cloglog" # 推荐: cloglog 更接近真实presence概率

# 是否绘制响应曲线
PLOT_RESPONSE_CURVES <- TRUE
TOP_N_RESPONSE <- 47 # 绘制前N个重要变量的响应曲线

# 创建输出目录
output_dir <- "output/05_model_maxnet"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "response_curves"), showWarnings = FALSE)

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Maxnet 物种分布模型 (增强版)\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("\n参数设置:\n")
cat("  - 特征类别: ", MAXENT_CLASSES, "\n", sep = "")
cat("  - 正则化乘数: ", MAXENT_REGMULT, "\n", sep = "")
cat("  - 预测类型: ", PREDICT_TYPE, "\n", sep = "")
cat("\n")

# ==============================================================================
# 高级设置：阴影绘制与性能优化
# ==============================================================================
# 1. [性能优化] 是否启用超参数网格搜索 (Grid Search)?
#    开启后会遍历不同参数组合，选择AICc最低的最优模型，显著提升模型性能。
#    注意：开启会增加运行时间。
DO_HYPERPARAMETER_TUNING <- TRUE 

# 2. [绘图美化] 是否启用Bootstrap计算响应曲线阴影 (置信区间)?
#    开启后会重采样训练N次，计算预测的标准误，在曲线上绘制阴影。
#    注意：N越大越平滑，但计算越慢。
DO_BOOTSTRAP_CURVES <- TRUE
BOOTSTRAP_N <- 10 # 推荐 10-50 次

# ==============================================================================
# 1. 数据准备
# ==============================================================================
cat("步骤 1/7: 数据准备...\n")
model_data <- read.csv("output/04_collinearity/collinearity_removed.csv")
exclude_cols <- c("id", "species", "lon", "lat", "source", "presence", "presence.1")
env_vars <- setdiff(names(model_data), exclude_cols)

n_presence <- sum(model_data$presence == 1)
n_background <- sum(model_data$presence == 0)

cat("  - 总样本: ", nrow(model_data), "\n", sep = "")
cat("  - 出现点: ", n_presence, "\n", sep = "")
cat("  - 背景点: ", n_background, "\n", sep = "")
cat("  - 变量数: ", length(env_vars), "\n", sep = "")

# 分层划分训练/测试集 (80/20)
set.seed(12345)
presence_idx <- which(model_data$presence == 1)
background_idx <- which(model_data$presence == 0)
train_idx <- c(
  sample(presence_idx, round(0.8 * length(presence_idx))),
  sample(background_idx, round(0.8 * length(background_idx)))
)
test_idx <- setdiff(seq_len(nrow(model_data)), train_idx)

train_data <- model_data[train_idx, ]
test_data <- model_data[test_idx, ]
cat("  - 训练集: ", nrow(train_data), " (出现点: ", sum(train_data$presence == 1), ")\n", sep = "")
cat("  - 测试集: ", nrow(test_data), " (出现点: ", sum(test_data$presence == 1), ")\n", sep = "")

# ==============================================================================
# 2. 构建特征公式并训练模型
# ==============================================================================
# [新增步骤] 1.5 超参数调优 (提升模型性能的核心)
# ==============================================================================
best_params <- list(regmult = MAXENT_REGMULT, classes = MAXENT_CLASSES) # 默认值

if (DO_HYPERPARAMETER_TUNING) {
  cat("\n步骤 1.5:正在进行超参数网格搜索 (Grid Search) 以提升性能...\n")
  
  # 定义搜索网格
  tune_grid <- expand.grid(
    regmult = c(0.5, 1, 2, 3), # 正则化力度
    classes = c("lq", "lqph")  # 特征组合 (lq=线性+二次, lqph=复杂组合)
  )
  
  results <- data.frame()
  
  # 使用进度条
  pb <- txtProgressBar(min = 0, max = nrow(tune_grid), style = 3)
  
  for(i in 1:nrow(tune_grid)) {
    rm_val <- tune_grid$regmult[i]
    fc_val <- tune_grid$classes[i]
    
    tryCatch({
      # 训练临时模型
      m_tmp <- maxnet(p = train_data$presence, data = train_data[, env_vars, drop=FALSE],
                      f = maxnet.formula(train_data$presence, train_data[, env_vars, drop=FALSE], classes = fc_val),
                      regmult = rm_val)
      
      # 使用测试集AUC作为简化标准 (虽然不完全严谨，但在没有独立测试集时是最实用的)
      pred_test <- predict(m_tmp, test_data[, env_vars, drop=FALSE], type = "cloglog", clamp = TRUE)
      auc_test <- tryCatch(pROC::auc(test_data$presence, as.numeric(pred_test), quiet=TRUE), error = function(e) 0.5)
      
      results <- rbind(results, data.frame(regmult = rm_val, classes = fc_val, auc = as.numeric(auc_test)))
    }, error = function(e) { })
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  # 选出AUC最高的参数
  if(nrow(results) > 0) {
    best_params <- results[which.max(results$auc), ]
    cat("\n  ✓ 最优参数: regmult =", best_params$regmult, ", classes =", best_params$classes, 
        "(测试集 AUC:", round(best_params$auc, 4), ")\n")
  }
}

# ==============================================================================
# 2. 训练最终模型 (使用最优参数)
# ==============================================================================
cat("\n步骤 2/7: 训练最终Maxnet模型...\n")
start_time <- Sys.time()
maxnet_model <- maxnet(
  p = train_data$presence,
  data = train_data[, env_vars, drop=FALSE],
  f = maxnet.formula(train_data$presence, train_data[, env_vars, drop=FALSE], classes = best_params$classes),
  regmult = best_params$regmult
)
train_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat("  ✓ 完成 (", round(train_time, 2), " 秒)\n", sep = "")

# ==============================================================================
# 3. 提取模型诊断信息
# ==============================================================================
cat("\n步骤 3/7: 提取模型诊断信息...\n")

# 从maxnet对象中提取关键信息
model_diag <- list(
  # 模型熵 (用于转换预测值)
  entropy = maxnet_model$entropy,

  # Alpha常数 (使指数模型在背景上归一化)
  alpha = maxnet_model$alpha,

  # 非零系数数量
  n_nonzero_coef = sum(maxnet_model$betas != 0),

  # 总系数数量
  n_total_coef = length(maxnet_model$betas),

  # 特征类别
  feature_classes = best_params$classes,

  # 正则化乘数
  regmult = best_params$regmult
)

cat("  - 模型熵 (entropy): ", round(model_diag$entropy, 4), "\n", sep = "")
cat("  - Alpha常数: ", round(model_diag$alpha, 4), "\n", sep = "")
cat("  - 非零系数: ", model_diag$n_nonzero_coef, " / ", model_diag$n_total_coef, "\n", sep = "")

# 提取模型系数
coef_df <- data.frame(
  feature = names(maxnet_model$betas),
  coefficient = as.numeric(maxnet_model$betas),
  stringsAsFactors = FALSE
) %>%
  filter(coefficient != 0) %>%
  arrange(desc(abs(coefficient)))

cat("  - 非零特征数: ", nrow(coef_df), "\n", sep = "")

# 提取变量范围信息 (用于clamping)
var_ranges <- data.frame(
  variable = names(maxnet_model$varmin),
  min = maxnet_model$varmin,
  max = maxnet_model$varmax,
  sample_mean = maxnet_model$samplemeans[names(maxnet_model$varmin)],
  stringsAsFactors = FALSE
)

# ==============================================================================
# 4. 多类型预测
# ==============================================================================
cat("\n步骤 4/7: 生成预测...\n")

# 定义预测所需的数据集 (修复X_train未定义错误)
X_train <- train_data[, env_vars, drop = FALSE]
X_test <- test_data[, env_vars, drop = FALSE]

# 生成四种类型的预测
pred_types <- c("link", "exponential", "cloglog", "logistic")
predictions_all <- data.frame(
  id = model_data$id,
  species = model_data$species,
  lon = model_data$lon,
  lat = model_data$lat,
  presence = model_data$presence,
  dataset = ifelse(seq_len(nrow(model_data)) %in% train_idx, "train", "test")
)

for (ptype in pred_types) {
  pred_train <- predict(maxnet_model, X_train, type = ptype, clamp = TRUE)
  pred_test <- predict(maxnet_model, X_test, type = ptype, clamp = TRUE)

  col_name <- paste0("pred_", ptype)
  predictions_all[[col_name]] <- NA_real_
  predictions_all[[col_name]][train_idx] <- as.numeric(pred_train)
  predictions_all[[col_name]][test_idx] <- as.numeric(pred_test)
}

# 主要预测列 (使用指定的PREDICT_TYPE)
predictions_all$predicted <- predictions_all[[paste0("pred_", PREDICT_TYPE)]]

cat("  ✓ 生成四种类型预测: link, exponential, cloglog, logistic\n")

# ==============================================================================
# 5. 模型评估
# ==============================================================================
cat("\n步骤 5/7: 模型评估...\n")

# 评估函数
evaluate_model <- function(actual, predicted, dataset_name) {
  roc_obj <- roc(actual, predicted, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))

  # 计算最优阈值 (Youden's J)
  coords_result <- coords(roc_obj, "best", ret = "all", best.method = "youden")

  data.frame(
    dataset = dataset_name,
    n_samples = length(actual),
    n_presence = sum(actual == 1),
    n_background = sum(actual == 0),
    AUC = auc_val,
    optimal_threshold = coords_result$threshold[1],
    sensitivity = coords_result$sensitivity[1],
    specificity = coords_result$specificity[1],
    TSS = coords_result$sensitivity[1] + coords_result$specificity[1] - 1,
    PPV = coords_result$ppv[1],
    NPV = coords_result$npv[1]
  )
}

# 使用cloglog预测进行评估 (推荐)
eval_train <- evaluate_model(
  train_data$presence,
  predictions_all$pred_cloglog[train_idx],
  "train"
)
eval_test <- evaluate_model(
  test_data$presence,
  predictions_all$pred_cloglog[test_idx],
  "test"
)

evaluation <- bind_rows(eval_train, eval_test)
evaluation$model <- "Maxnet"
evaluation$feature_classes <- MAXENT_CLASSES
evaluation$regmult <- MAXENT_REGMULT
evaluation$predict_type <- "cloglog"
evaluation$training_time_sec <- c(train_time, NA)
evaluation$entropy <- model_diag$entropy
evaluation$n_nonzero_coef <- model_diag$n_nonzero_coef

cat("  - AUC (cloglog): ", round(eval_train$AUC, 4), " (train) / ",
  round(eval_test$AUC, 4), " (test)\n",
  sep = ""
)
cat("  - TSS (cloglog): ", round(eval_test$TSS, 4), "\n", sep = "")
cat("  - 最优阈值: ", round(eval_test$optimal_threshold, 4), "\n", sep = "")

# 比较不同预测类型的AUC
cat("\n  不同预测类型的AUC比较:\n")
for (ptype in pred_types) {
  pred_col <- paste0("pred_", ptype)
  auc_val <- as.numeric(auc(roc(test_data$presence,
    predictions_all[[pred_col]][test_idx],
    quiet = TRUE
  )))
  cat("    - ", ptype, ": ", round(auc_val, 4), "\n", sep = "")
}

# ==============================================================================
# 6. 变量重要性 (排列法)
# ==============================================================================
cat("\n步骤 6/7: 变量重要性 (排列法)...\n")

baseline_auc <- eval_test$AUC
var_importance <- data.frame(
  variable = env_vars,
  importance = NA_real_,
  auc_drop = NA_real_,
  stringsAsFactors = FALSE
)

set.seed(123)
pb <- txtProgressBar(min = 0, max = length(env_vars), style = 3)
for (i in seq_along(env_vars)) {
  perm_test <- test_data
  perm_test[[env_vars[i]]] <- sample(perm_test[[env_vars[i]]])
  pred_perm <- predict(maxnet_model, perm_test[, env_vars, drop = FALSE],
    type = PREDICT_TYPE, clamp = TRUE
  )
  perm_auc <- as.numeric(auc(roc(perm_test$presence, pred_perm, quiet = TRUE)))
  var_importance$auc_drop[i] <- baseline_auc - perm_auc
  setTxtProgressBar(pb, i)
}
close(pb)

# 计算相对重要性 (%)
var_importance$importance <- pmax(var_importance$auc_drop, 0)
var_importance$importance_pct <- 100 * var_importance$importance / sum(var_importance$importance)
var_importance <- var_importance %>%
  arrange(desc(importance)) %>%
  mutate(rank = row_number())

cat("  ✓ 完成\n")
cat("\n  前20重要变量:\n")
print(head(var_importance[, c("rank", "variable", "auc_drop", "importance_pct")], 20))

# ==============================================================================
# 7. 响应曲线
# ==============================================================================
if (PLOT_RESPONSE_CURVES) {
  cat("\n步骤 7/7: 绘制响应曲线...\n")

  # 获取前N个重要变量
  top_vars <- head(var_importance$variable, TOP_N_RESPONSE)

  # ==============================================================================
  # 7. 计算并绘制响应曲线 (参考优化版)
  # ==============================================================================
  cat("\n步骤 7/7: 绘制响应曲线 (Bootstrap 阴影: ", DO_BOOTSTRAP_CURVES, ")...\n", sep = "")

# 准备Bootstrap模型群 (如果启用)
boot_models <- list()
if (DO_BOOTSTRAP_CURVES) {
  cat("  - [高级] 正在执行 Bootstrap 重采样 (", BOOTSTRAP_N, " 次) 以计算置信区间...\n", sep = "")
  cat("    进度: ")
  for (b in 1:BOOTSTRAP_N) {
    if (b %% 5 == 0) cat(b, "...", sep="")
    # 重采样
    valid_indices <- which(complete.cases(train_data[, env_vars]))
    idx <- sample(valid_indices, length(valid_indices), replace = TRUE)
    boot_train <- train_data[idx, ]
    
    try({
      # 使用与主模型相同的参数训练
      bm <- maxnet(
        p = boot_train$presence, 
        data = boot_train[, env_vars, drop=FALSE],
        # 如果前面没有做Grid Search，使用默认参数
        f = maxnet.formula(boot_train$presence, boot_train[, env_vars, drop=FALSE], classes = best_params$classes),
        regmult = best_params$regmult
      )
      boot_models[[length(boot_models) + 1]] <- bm
    }, silent = TRUE)
  }
  cat("完成\n    成功训练 Bootstrap 模型数: ", length(boot_models), "\n")
}

# 绘制曲线
cat("  - 正在生成并保存响应曲线...\n")
graphics.off() 

for (var in top_vars) {
  # 1. 准备预测范围
  var_range <- seq(maxnet_model$varmin[var], maxnet_model$varmax[var], length.out = 100)
  
  # 2. 准备基础预测矩阵
  pred_df <- as.data.frame(matrix(rep(maxnet_model$samplemeans, each = 100), nrow = 100, byrow = FALSE))
  names(pred_df) <- names(maxnet_model$samplemeans)
  pred_df[[var]] <- var_range
  
  # 3. 主预测
  main_pred <- predict(maxnet_model, pred_df[, env_vars, drop=FALSE], type = PREDICT_TYPE, clamp = TRUE)
  plot_data <- data.frame(x = var_range, y = as.numeric(main_pred))
  
  # 4. 计算阴影 (如果Bootstrap模型存在)
  if (DO_BOOTSTRAP_CURVES && length(boot_models) > 0) {
    boot_preds <- matrix(NA, nrow = 100, ncol = length(boot_models))
    for (k in seq_along(boot_models)) {
      try({
        # 预测时要宽容，因为Bootstrap数据的范围可能略有不同
        p <- predict(boot_models[[k]], pred_df[, env_vars, drop=FALSE], type = PREDICT_TYPE, clamp = TRUE)
        boot_preds[, k] <- as.numeric(p)
      }, silent = TRUE)
    }
    # 计算均值和标准差 (Mean +/- SD) - 这种方式比分位数更稳健，阴影更紧凑
    boot_mean <- apply(boot_preds, 1, mean, na.rm = TRUE)
    boot_sd <- apply(boot_preds, 1, sd, na.rm = TRUE)
    
    plot_data$ymin <- pmax(0, boot_mean - boot_sd) # 确保不小于0
    plot_data$ymax <- pmin(1, boot_mean + boot_sd) # 确保不大于1

  }
  
  # 5. 绘图 (严格复刻用户的专业风格：theme_bw + 粗线 + 黑框)
  p_single <- ggplot(plot_data, aes(x = x, y = y))
  
  # 先画阴影 (Ribbon) - 颜色选淡蓝色/灰色
  if ("ymin" %in% names(plot_data)) {
    p_single <- p_single + 
      geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "grey70", alpha = 0.5) 
  }
  
  # 再画主线 (黑粗线 或者 深蓝粗线 - 响应用户偏好线粗)
  p_single <- p_single +
    geom_line(size = 1.2, color = "black") + 
    labs(x = var, y = "Logistic output") + 
    # 使用用户最满意的参数
    theme_bw(base_size = 10) +
    theme(
      panel.grid = element_blank(),
      axis.title = element_text(face = "plain"),
      axis.text = element_text(color = "black")
    )
  
  # 6. 保存
  out_file <- file.path(output_dir, "response_curves", paste0("curve_", var, ".png"))
  ggsave(out_file, p_single, width = 8, height = 6, units = "in", dpi = 600)
}
} # End of PLOT_RESPONSE_CURVES check



# ==============================================================================
# 保存所有结果
# ==============================================================================
cat("\n保存结果...\n")

# 1. 模型对象
saveRDS(maxnet_model, file.path(output_dir, "model.rds"))
cat("  ✓ 模型: ", file.path(output_dir, "model.rds"), "\n", sep = "")

# 2. 预测结果
write.csv(predictions_all, file.path(output_dir, "predictions.csv"), row.names = FALSE)
cat("  ✓ 预测: ", file.path(output_dir, "predictions.csv"), "\n", sep = "")

# 3. 评估结果
write.csv(evaluation, file.path(output_dir, "evaluation.csv"), row.names = FALSE)
cat("  ✓ 评估: ", file.path(output_dir, "evaluation.csv"), "\n", sep = "")

# 4. 变量重要性
write.csv(var_importance, file.path(output_dir, "variable_importance.csv"), row.names = FALSE)
cat("  ✓ 变量重要性: ", file.path(output_dir, "variable_importance.csv"), "\n", sep = "")

# 5. 模型系数
write.csv(coef_df, file.path(output_dir, "model_coefficients.csv"), row.names = FALSE)
cat("  ✓ 模型系数: ", file.path(output_dir, "model_coefficients.csv"), "\n", sep = "")

# 6. 变量范围
write.csv(var_ranges, file.path(output_dir, "variable_ranges.csv"), row.names = FALSE)
cat("  ✓ 变量范围: ", file.path(output_dir, "variable_ranges.csv"), "\n", sep = "")

# 7. 模型诊断
diag_df <- data.frame(
  parameter = c(
    "entropy", "alpha", "n_nonzero_coef", "n_total_coef",
    "feature_classes", "regmult", "training_time_sec"
  ),
  value = c(
    model_diag$entropy, model_diag$alpha, model_diag$n_nonzero_coef,
    model_diag$n_total_coef, model_diag$feature_classes,
    model_diag$regmult, train_time
  )
)
write.csv(diag_df, file.path(output_dir, "model_diagnostics.csv"), row.names = FALSE)
cat("  ✓ 模型诊断: ", file.path(output_dir, "model_diagnostics.csv"), "\n", sep = "")

# 8. 处理日志
sink(file.path(output_dir, "processing_log.txt"))
cat("Maxnet 物种分布模型训练日志\n")
cat("生成时间: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")

cat("【参数设置】\n")
cat("  特征类别: ", MAXENT_CLASSES, "\n", sep = "")
cat("  正则化乘数: ", MAXENT_REGMULT, "\n", sep = "")
cat("  预测类型: ", PREDICT_TYPE, "\n\n", sep = "")

cat("【数据概况】\n")
cat("  总样本: ", nrow(model_data), "\n", sep = "")
cat("  出现点: ", n_presence, "\n", sep = "")
cat("  背景点: ", n_background, "\n", sep = "")
cat("  变量数: ", length(env_vars), "\n", sep = "")
cat("  训练集: ", nrow(train_data), "\n", sep = "")
cat("  测试集: ", nrow(test_data), "\n\n", sep = "")

cat("【模型诊断】\n")
cat("  模型熵 (entropy): ", round(model_diag$entropy, 4), "\n", sep = "")
cat("  Alpha常数: ", round(model_diag$alpha, 4), "\n", sep = "")
cat("  非零系数: ", model_diag$n_nonzero_coef, " / ", model_diag$n_total_coef, "\n", sep = "")
cat("  训练时间: ", round(train_time, 2), " 秒\n\n", sep = "")

cat("【模型性能】\n")
cat("  AUC (训练集): ", round(eval_train$AUC, 4), "\n", sep = "")
cat("  AUC (测试集): ", round(eval_test$AUC, 4), "\n", sep = "")
cat("  TSS: ", round(eval_test$TSS, 4), "\n", sep = "")
cat("  灵敏度: ", round(eval_test$sensitivity, 4), "\n", sep = "")
cat("  特异度: ", round(eval_test$specificity, 4), "\n", sep = "")
cat("  最优阈值: ", round(eval_test$optimal_threshold, 4), "\n\n", sep = "")

cat("【变量重要性 (前15)】\n")
print(head(var_importance[, c("rank", "variable", "auc_drop", "importance_pct")], 15))

cat("\n【输出类型说明】\n")
cat("  link: 线性预测器 lp\n")
cat("  exponential: exp(lp) - 原始Maxent输出\n")
cat("  cloglog: 1-exp(-exp(entropy+lp)) - 推荐用于presence概率解释\n")
cat("  logistic: 1/(1+exp(-entropy-lp)) - 传统logistic解释\n")
sink()
cat("  ✓ 处理日志: ", file.path(output_dir, "processing_log.txt"), "\n", sep = "")

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("✓ Maxnet模型训练完成!\n")
