#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 07_model_gam.R
# 功能说明: 使用广义加性模型(GAM)进行物种分布建模
# 方法: GAM with smooth terms
# 输入文件: ../output/04_collinearity_removed.csv
# 输出文件: ../output/07_gam_model.rds
#          ../output/07_gam_predictions.csv
#          ../output/07_gam_evaluation.csv
#          ../output/07_gam_variable_importance.csv
# 作者: SDM Analysis Pipeline
# 日期: 2025-10-02
# ==============================================================================

# 清空环境
rm(list = ls())
gc()

# 设置工作目录
setwd("E:/SDM01")

# 加载必要的包
packages <- c("tidyverse", "mgcv", "caret", "pROC")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 创建输出目录
if(!dir.exists("output")) dir.create("output")
if(!dir.exists("output/07_model_gam")) dir.create("output/07_model_gam", recursive = TRUE)

cat("======================================\n")
cat("GAM (广义加性模型) 训练\n")
cat("======================================\n\n")

# ------------------------------------------------------------------------------
# 1. 读取数据
# ------------------------------------------------------------------------------
cat("步骤 1/8: 读取建模数据...\n")

model_data <- read.csv("output/04_collinearity/collinearity_removed.csv")
cat("  - 总记录数: ", nrow(model_data), "\n", sep = "")
cat("  - 出现点: ", sum(model_data$presence == 1), "\n", sep = "")
cat("  - 背景点: ", sum(model_data$presence == 0), "\n", sep = "")

# 提取环境变量（排除ID、坐标、presence列）
# 数据驱动：只保留有效的环境变量列
exclude_cols <- c("id", "species", "lon", "lat", "source", "presence", "presence.1")
all_cols <- names(model_data)
env_vars <- setdiff(all_cols, exclude_cols)

cat("  - 原始环境变量数: ", length(env_vars), "\n", sep = "")
cat("  - 数据完整性: ", sum(complete.cases(model_data)), " / ", nrow(model_data), "\n", sep = "")

# ------------------------------------------------------------------------------
# 2. 数据划分: 训练集和测试集
# ------------------------------------------------------------------------------
cat("\n步骤 2/8: 划分训练集和测试集...\n")

# 设置随机种子
set.seed(12345)

# 分别对出现点和背景点进行分层抽样
presence_indices <- which(model_data$presence == 1)
background_indices <- which(model_data$presence == 0)

# 80%训练, 20%测试
train_presence <- sample(presence_indices, 
                         size = round(0.8 * length(presence_indices)))
train_background <- sample(background_indices, 
                           size = round(0.8 * length(background_indices)))

train_indices <- c(train_presence, train_background)
test_indices <- setdiff(1:nrow(model_data), train_indices)

train_data <- model_data[train_indices, ]
test_data <- model_data[test_indices, ]

cat("  - 训练集样本数: ", nrow(train_data), 
    " (出现:", sum(train_data$presence == 1), 
    ", 背景:", sum(train_data$presence == 0), ")\n", sep = "")
cat("  - 测试集样本数: ", nrow(test_data), 
    " (出现:", sum(test_data$presence == 1), 
    ", 背景:", sum(test_data$presence == 0), ")\n", sep = "")

# ------------------------------------------------------------------------------
# 3. 构建GAM模型公式（数据驱动的k值选择）
# ------------------------------------------------------------------------------
cat("\n步骤 3/8: 构建GAM模型公式...\n")

# 数据驱动策略：检查每个变量的唯一值数量，过滤掉不适合平滑的变量
valid_vars <- c()
k_values <- c()

for(v in env_vars) {
  n_unique <- length(unique(train_data[[v]]))
  
  # 至少需要10个唯一值才能使用平滑项
  if(n_unique >= 10) {
    valid_vars <- c(valid_vars, v)
    # 提升灵活度：k 随唯一值数量而增，最小5，最大15（避免过拟合）
    k_val <- min(max(5, floor(n_unique / 4)), 15)
    k_values <- c(k_values, k_val)
  } else {
    cat("    警告: 变量", v, "唯一值太少(", n_unique, ")，跳过\n")
  }
}

cat("  有效变量数: ", length(valid_vars), "/", length(env_vars), "\n", sep = "")
cat("  变量自由度设置: k范围 =", min(k_values), "-", max(k_values), "\n")

# 构建平滑项
smooth_terms <- paste0("s(", valid_vars, ", k=", k_values, ")", collapse = " + ")

# 添加空间平滑项（控制空间自相关）——提升 k 以提高空间分辨率
formula_str <- paste("presence ~", smooth_terms, "+ s(lon, lat, k=80)")
gam_formula <- as.formula(formula_str)

cat("  - GAM公式构建完成\n")
cat("    * 环境变量平滑项: ", length(valid_vars), " 个\n", sep = "")
cat("    * 空间平滑项: s(lon, lat, k=80)\n")

# ------------------------------------------------------------------------------
# 4. 训练GAM模型
# ------------------------------------------------------------------------------
cat("\n步骤 4/8: 训练GAM模型...\n")
cat("  注意: 这可能需要较长时间...\n")

# 为类别不平衡设置权重
n_presence_train <- sum(train_data$presence == 1)
n_background_train <- sum(train_data$presence == 0)
weight_ratio <- n_background_train / n_presence_train

train_data$weight <- ifelse(train_data$presence == 1, 
                             weight_ratio, 
                             1)

start_time <- Sys.time()

# 使用binomial family (logistic regression)
# method = "REML" (Restricted Maximum Likelihood)
gam_model <- tryCatch({
  bam(  # 使用bam代替gam,处理大数据集更快
    formula = gam_formula,
    family = binomial(link = "logit"),
    data = train_data,
    weights = weight,
    method = "fREML",      # fast REML
    select = TRUE,         # 额外惩罚进行变量选择（防止过拟合）
    discrete = TRUE,       # 加速
    gamma = 1.2            # 轻微提高惩罚强度，提升稳健性
  )
}, error = function(e) {
  cat("  - bam失败,尝试使用标准gam...\n")
  gam(
    formula = gam_formula,
    family = binomial(link = "logit"),
    data = train_data,
    weights = weight,
    method = "REML",
    select = TRUE,
    gamma = 1.2
  )
})

end_time <- Sys.time()
training_time <- as.numeric(difftime(end_time, start_time, units = "secs"))

cat("  - 模型训练完成 (耗时: ", round(training_time, 2), " 秒)\n", sep = "")
cat("  - Deviance explained: ", 
    round(summary(gam_model)$dev.expl * 100, 2), "%\n", sep = "")

# 诊断：并发曲率与 k-index（保存日志，便于后续增强）
try({
  conc <- mgcv::concurvity(gam_model, full = TRUE)
  sink("output/07_model_gam/diagnostics.txt")
  cat("Concurvity summary (diagonal/overall)\n\n")
  print(lapply(conc, function(m) round(apply(m, 2, max, na.rm=TRUE), 3)))
  cat("\nGAM check (k-index)\n\n")
  print(utils::capture.output(mgcv::gam.check(gam_model)) )
  sink()
}, silent = TRUE)

# ------------------------------------------------------------------------------
# 5. 模型预测
# ------------------------------------------------------------------------------
cat("\n步骤 5/8: 进行预测...\n")

# 训练集预测(使用response类型得到概率)
train_pred <- predict(gam_model, 
                      newdata = train_data, 
                      type = "response")

# 测试集预测
test_pred <- predict(gam_model, 
                     newdata = test_data, 
                     type = "response")

cat("  - 训练集预测完成\n")
cat("  - 测试集预测完成\n")

# ------------------------------------------------------------------------------
# 6. 模型评估
# ------------------------------------------------------------------------------
cat("\n步骤 6/8: 模型评估...\n")

# 计算AUC
train_roc <- roc(train_data$presence, train_pred, quiet = TRUE)
test_roc <- roc(test_data$presence, test_pred, quiet = TRUE)

train_auc <- auc(train_roc)
test_auc <- auc(test_roc)

cat("  - 训练集 AUC: ", round(train_auc, 4), "\n", sep = "")
cat("  - 测试集 AUC: ", round(test_auc, 4), "\n", sep = "")

# 使用最大化TSS的阈值
coords_result <- coords(test_roc, "best", ret = "all", 
                        best.method = "youden")
optimal_threshold <- coords_result$threshold
test_sensitivity <- coords_result$sensitivity
test_specificity <- coords_result$specificity
test_tss <- test_sensitivity + test_specificity - 1

cat("  - 最优阈值: ", round(optimal_threshold, 4), "\n", sep = "")
cat("  - TSS: ", round(test_tss, 4), "\n", sep = "")

# 计算Kappa
test_pred_binary <- ifelse(test_pred >= optimal_threshold, 1, 0)
confusion_matrix <- table(Predicted = test_pred_binary, 
                          Observed = test_data$presence)

if(nrow(confusion_matrix) == 2 && ncol(confusion_matrix) == 2) {
  # 计算Kappa
  observed_accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  expected_accuracy <- sum(rowSums(confusion_matrix) * 
                           colSums(confusion_matrix)) / sum(confusion_matrix)^2
  kappa <- (observed_accuracy - expected_accuracy) / (1 - expected_accuracy)
  
  cat("  - Kappa: ", round(kappa, 4), "\n", sep = "")
  
  # 打印混淆矩阵
  cat("\n  混淆矩阵 (测试集):\n")
  print(confusion_matrix)
} else {
  kappa <- NA
  cat("  - Kappa: 无法计算\n")
}

# ------------------------------------------------------------------------------
# 7. 变量重要性
# ------------------------------------------------------------------------------
cat("\n步骤 7/8: 计算变量重要性...\n")

# GAM的变量重要性基于统计显著性和效应大小
gam_summary <- summary(gam_model)

# 初始化变量重要性数据框
var_importance <- data.frame(
  variable = character(),
  type = character(),
  statistic = numeric(),
  p_value = numeric(),
  importance = numeric(),
  stringsAsFactors = FALSE
)

# 1. 提取平滑项的统计量
if("s.table" %in% names(gam_summary) && nrow(gam_summary$s.table) > 0) {
  smooth_stats <- gam_summary$s.table
  
  # 检查列名（可能是F或Chi.sq）
  stat_col <- if("F" %in% colnames(smooth_stats)) {
    "F"
  } else if("Chi.sq" %in% colnames(smooth_stats)) {
    "Chi.sq"
  } else {
    NULL
  }
  
  if(!is.null(stat_col)) {
    # 提取平滑项变量名（从行名中提取）
    smooth_vars <- gsub("s\\((.+?)\\)", "\\1", rownames(smooth_stats))
    smooth_vars <- gsub(",k=.*", "", smooth_vars)  # 移除k参数
    
    smooth_importance <- data.frame(
      variable = smooth_vars,
      type = "smooth",
      statistic = smooth_stats[, stat_col],
      p_value = smooth_stats[, "p-value"],
      stringsAsFactors = FALSE
    )
    
    # 计算重要性分数
    smooth_importance$importance <- smooth_importance$statistic * 
                                    (1 - pmin(smooth_importance$p_value, 0.999))
    
    var_importance <- rbind(var_importance, smooth_importance)
  }
}

# 2. 提取线性项（参数项）的统计量
if("p.table" %in% names(gam_summary) && nrow(gam_summary$p.table) > 1) {
  param_stats <- gam_summary$p.table[-1, , drop = FALSE]  # 移除截距
  
  if(nrow(param_stats) > 0) {
    param_vars <- rownames(param_stats)
    
    param_importance <- data.frame(
      variable = param_vars,
      type = "linear",
      statistic = abs(param_stats[, "t value"]),
      p_value = param_stats[, "Pr(>|t|)"],
      stringsAsFactors = FALSE
    )
    
    # 计算重要性分数
    param_importance$importance <- param_importance$statistic * 
                                   (1 - pmin(param_importance$p_value, 0.999))
    
    var_importance <- rbind(var_importance, param_importance)
  }
}

# 按重要性排序
if(nrow(var_importance) > 0) {
  var_importance <- var_importance[order(var_importance$importance, 
                                         decreasing = TRUE), ]
  
  
  # 打印前10个最重要变量
  cat("\n  前10个最重要变量:\n")
  print(head(var_importance[, c("variable", "type", "importance")], 10))
  
} else {
  # 如果无法提取统计表,使用排列重要性
  cat("  - 警告: 无法从模型中提取统计表,使用排列重要性方法...\n")
  
  var_importance <- data.frame(
    variable = valid_vars,
    type = "permutation",
    statistic = NA,
    p_value = NA,
    importance = 0,
    stringsAsFactors = FALSE
  )
  
  baseline_auc <- test_auc
  
  for(i in 1:length(valid_vars)) {
    var_name <- valid_vars[i]
    test_permuted <- test_data
    test_permuted[, var_name] <- sample(test_permuted[, var_name])
    
    pred_permuted <- predict(gam_model, 
                             newdata = test_permuted, 
                             type = "response")
    roc_permuted <- roc(test_permuted$presence, pred_permuted, quiet = TRUE)
    auc_permuted <- auc(roc_permuted)
    
    var_importance$importance[i] <- baseline_auc - auc_permuted
    
    if(i %% 5 == 0) {
      cat("    - 已处理 ", i, "/", length(valid_vars), " 个变量\n", sep = "")
    }
  }
  
  var_importance <- var_importance[order(var_importance$importance, 
                                         decreasing = TRUE), ]
  
  # 打印前10个最重要变量
  cat("\n  前10个最重要变量:\n")
  print(head(var_importance[, c("variable", "importance")], 10))
}

# 保存变量重要性
write.csv(var_importance,
          "output/07_model_gam/variable_importance.csv",
          row.names = FALSE)
cat("\n  - 已保存: output/07_gam_variable_importance.csv\n")

# ------------------------------------------------------------------------------
# 8. 保存模型和结果
# ------------------------------------------------------------------------------
cat("\n步骤 8/8: 保存模型和结果...\n")

# 保存模型对象
saveRDS(gam_model, "output/07_model_gam/model.rds")
cat("  - 已保存模型: output/07_gam_model.rds\n")

# 保存预测结果
predictions_df <- data.frame(
  id = model_data$id,
  species = model_data$species,
  lon = model_data$lon,
  lat = model_data$lat,
  presence = model_data$presence,
  dataset = ifelse(1:nrow(model_data) %in% train_indices, "train", "test"),
  predicted = NA
)
predictions_df$predicted[train_indices] <- train_pred
predictions_df$predicted[test_indices] <- test_pred

write.csv(predictions_df,
          "output/07_model_gam/predictions.csv",
          row.names = FALSE)
cat("  - 已保存预测结果: output/07_gam_predictions.csv\n")

# 保存评估指标
evaluation_df <- data.frame(
  model = "GAM",
  dataset = c("train", "test"),
  n_samples = c(nrow(train_data), nrow(test_data)),
  n_presence = c(sum(train_data$presence == 1), sum(test_data$presence == 1)),
  AUC = c(train_auc, test_auc),
  TSS = c(NA, test_tss),
  Sensitivity = c(NA, test_sensitivity),
  Specificity = c(NA, test_specificity),
  Kappa = c(NA, kappa),
  optimal_threshold = c(NA, optimal_threshold),
  training_time_sec = c(training_time, NA),
  deviance_explained = c(summary(gam_model)$dev.expl, NA)
)

write.csv(evaluation_df,
          "output/07_model_gam/evaluation.csv",
          row.names = FALSE)
cat("  - 已保存评估结果: output/07_gam_evaluation.csv\n")

# ------------------------------------------------------------------------------
# 9. 数据摘要
# ------------------------------------------------------------------------------
cat("\n======================================\n")
cat("GAM 模型训练完成 - 摘要\n")
cat("======================================\n")
cat("模型类型: GAM (Generalized Additive Model)\n")
cat("训练样本数: ", nrow(train_data), "\n", sep = "")
cat("测试样本数: ", nrow(test_data), "\n", sep = "")
cat("环境变量数: ", length(valid_vars), "\n", sep = "")
cat("平滑项数: ", length(valid_vars), "\n", sep = "")
cat("\n模型性能 (测试集):\n")
cat("  AUC: ", round(test_auc, 4), "\n", sep = "")
cat("  TSS: ", round(test_tss, 4), "\n", sep = "")
cat("  Sensitivity: ", round(test_sensitivity, 4), "\n", sep = "")
cat("  Specificity: ", round(test_specificity, 4), "\n", sep = "")
if(!is.na(kappa)) {
  cat("  Kappa: ", round(kappa, 4), "\n", sep = "")
}
cat("  Deviance Explained: ", 
    round(summary(gam_model)$dev.expl * 100, 2), "%\n", sep = "")

# 保存日志
sink("output/07_model_gam/processing_log.txt")
cat("GAM 模型训练日志\n")
cat("处理时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("模型参数:\n")
cat("  family: binomial (logit link)\n")
cat("  平滑项k值: 5\n")
cat("  method: fREML/REML\n")
cat("  使用权重平衡类别\n")
cat("\n训练信息:\n")
cat("  训练集样本:", nrow(train_data), "\n")
cat("  测试集样本:", nrow(test_data), "\n")
cat("  环境变量数:", length(valid_vars), "\n")
cat("  训练时间:", round(training_time, 2), "秒\n")
cat("\n模型性能:\n")
print(evaluation_df)
cat("\n变量重要性 (前10):\n")
print(head(var_importance, 10))
sink()

cat("\n脚本执行完成!\n")

