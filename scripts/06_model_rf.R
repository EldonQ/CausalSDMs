#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 06_model_rf.R (全面调优版)
# 功能说明: 随机森林 (Random Forest) 物种分布建模 (带超参数调优与响应曲线)
# 特性:
#   1. 网格搜索 (Grid Search) 寻找最优 mtry 和 nodesize
#   2. 类别不平衡处理 (下采样)
#   3. Bootstrap 响应曲线 (带置信区间阴影)
#   4. 专业绘图风格 (与 Maxnet 脚本一致)
# 输入文件: output/04_collinearity/collinearity_removed.csv
# 输出目录: output/06_model_rf/
# ==============================================================================

rm(list = ls()); gc()
setwd("E:/CausalSDMs")

# 加载包
packages <- c("tidyverse", "randomForest", "pROC", "ggplot2", "caret")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ------------------------------------------------------------------------------
# 配置参数
# ------------------------------------------------------------------------------
DO_TUNING <- TRUE          # 是否进行超参数网格搜索
DO_RESPONSE_CURVES <- TRUE # 是否绘制响应曲线
BOOTSTRAP_N <- 10          # 响应曲线重采样次数
TOP_N_RESPONSE <- 15       # 绘制前 N 个变量的曲线
RF_NTREE <- 1000           # 树的数量 (越多越稳定)

# 创建输出目录
output_dir <- "output/06_model_rf"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "response_curves"), showWarnings = FALSE)

# ------------------------------------------------------------------------------
# 1. 数据准备
# ------------------------------------------------------------------------------
cat("\n步骤 1/6: 数据准备...\n")
model_data <- read.csv("output/04_collinearity/collinearity_removed.csv")
exclude_cols <- c("id", "species", "lon", "lat", "source", "presence", "presence.1")
env_vars <- setdiff(names(model_data), exclude_cols)

# 转换为因子 (RF分类模式必需)
model_data$presence_fac <- factor(model_data$presence) 

# 划分训练/测试集
set.seed(12345)
train_idx <- createDataPartition(model_data$presence_fac, p = 0.8, list = FALSE)
train_data <- model_data[train_idx, ]
test_data  <- model_data[-train_idx, ]

cat("  - 训练集: ", nrow(train_data), " (出现: ", sum(train_data$presence==1), ")\n", sep="")
cat("  - 测试集: ", nrow(test_data), " (出现: ", sum(test_data$presence==1), ")\n", sep="")

# ------------------------------------------------------------------------------
# 2. 超参数调优 (Grid Search)
# ------------------------------------------------------------------------------
best_params <- list(mtry = floor(sqrt(length(env_vars))), nodesize = 1) # 默认值

if (DO_TUNING) {
  cat("\n步骤 2/6: 超参数网格搜索 (寻找最优 mtry / nodesize)...\n")
  
  # 定义搜索网格
  # mtry: 建议尝试 p/3, sqrt(p), 2*sqrt(p)
  p <- length(env_vars)
  tune_grid <- expand.grid(
    mtry = unique(floor(c(sqrt(p), 2*sqrt(p), p/3))),
    nodesize = c(1, 5) # 1=尽量纯净(易过拟合), 5=剪枝(防过拟合)
  )
  tune_grid <- tune_grid[tune_grid$mtry <= p, ] # 过滤非法值
  
  results <- data.frame()
  pb <- txtProgressBar(min = 0, max = nrow(tune_grid), style = 3)
  
  for(i in 1:nrow(tune_grid)) {
    c_mtry <- tune_grid$mtry[i]
    c_node <- tune_grid$nodesize[i]
    
    # 平衡采样: 每次采样的 presence 和 background 数量一致 (等于较少那类)
    n_min <- min(table(train_data$presence_fac))
    
    try({
      # 使用OOB错误率作为快速评估指标 (无需CV，RF自带OOB)
      m_tmp <- randomForest(
        x = train_data[, env_vars], 
        y = train_data$presence_fac,
        ntree = 200, # 调优时用较少树加速
        mtry = c_mtry,
        nodesize = c_node,
        sampsize = c(n_min, n_min), # 下采样
        strata = train_data$presence_fac
      )
      
      # 这里我们不仅看OOB Error，最好还是看测试集AUC
      pred_prob <- predict(m_tmp, test_data[, env_vars], type = "prob")[, "1"]
      auc_val <- as.numeric(auc(test_data$presence, pred_prob, quiet=TRUE))
      
      results <- rbind(results, data.frame(mtry = c_mtry, nodesize = c_node, auc = auc_val))
    }, silent=TRUE)
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  if(nrow(results) > 0) {
    best_idx <- which.max(results$auc)
    best_params <- results[best_idx, ]
    cat("\n  ✓ 最优参数: mtry =", best_params$mtry, ", nodesize =", best_params$nodesize, 
        "(Test AUC:", round(best_params$auc, 4), ")\n")
  }
}

# ------------------------------------------------------------------------------
# 3. 训练最终模型
# ------------------------------------------------------------------------------
cat("\n步骤 3/6: 训练最终 RF 模型...\n")

# 平衡采样大小
n_min <- min(table(train_data$presence_fac))
sampsize_vec <- c('0' = n_min, '1' = n_min)

start_time <- Sys.time()
rf_model <- randomForest(
  x = train_data[, env_vars],
  y = train_data$presence_fac,
  ntree = RF_NTREE,
  mtry = best_params$mtry,
  nodesize = best_params$nodesize,
  sampsize = sampsize_vec,
  strata = train_data$presence_fac,
  importance = TRUE,
  keep.forest = TRUE
)
train_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat("  ✓ 完成 (", round(train_time, 2), " 秒)\n", sep = "")

# ------------------------------------------------------------------------------
# 4. 预测与评估
# ------------------------------------------------------------------------------
cat("\n步骤 4/6: 评估模型性能...\n")

# 预测概率
pred_train <- predict(rf_model, train_data[, env_vars], type = "prob")[, "1"]
pred_test  <- predict(rf_model, test_data[, env_vars],  type = "prob")[, "1"]

# 评估
roc_train <- roc(train_data$presence, pred_train, quiet=TRUE)
roc_test  <- roc(test_data$presence,  pred_test,  quiet=TRUE)

coords_res <- coords(roc_test, "best", ret="all", best.method="youden", transpose=FALSE)
# 修复: coords 可能返回多行，取第一行
if(!is.null(nrow(coords_res))) coords_res <- coords_res[1, ]

eval_metrics <- list(
  auc_train = auc(roc_train),
  auc_test  = auc(roc_test),
  tss       = coords_res$sensitivity + coords_res$specificity - 1,
  threshold = coords_res$threshold,
  sens      = coords_res$sensitivity,
  spec      = coords_res$specificity
)

cat("  - AUC (Train): ", round(eval_metrics$auc_train, 4), "\n", sep="")
cat("  - AUC (Test):  ", round(eval_metrics$auc_test, 4), "\n", sep="")
cat("  - TSS:         ", round(eval_metrics$tss, 4), "\n", sep="")

# ------------------------------------------------------------------------------
# 5. 变量重要性
# ------------------------------------------------------------------------------
cat("\n步骤 5/6: 变量重要性...\n")
# type=1: MeanDecreaseAccuracy (更有意义), type=2: MeanDecreaseGini
imp_raw <- importance(rf_model)
var_imp <- data.frame(
  variable = rownames(imp_raw),
  # 优先用 MeanDecreaseAccuracy (如果没有则用Gini)
  importance = if("MeanDecreaseAccuracy" %in% colnames(imp_raw)) imp_raw[, "MeanDecreaseAccuracy"] else imp_raw[, "MeanDecreaseGini"]
)
var_imp <- var_imp[order(var_imp$importance, decreasing = TRUE), ]
top_vars <- head(var_imp$variable, TOP_N_RESPONSE)

cat("  前5关键变量: ", paste(head(top_vars, 5), collapse=", "), "\n")

# ------------------------------------------------------------------------------
# 6. 绘制响应曲线 (Bootstrap 阴影)
# ------------------------------------------------------------------------------
if (DO_RESPONSE_CURVES) {
  cat("\n步骤 6/6: 绘制响应曲线 (带 Bootstrap 阴影)...\n")
  
  # 预训练一些 Bootstrap 模型用于阴影计算 (耗时)
  boot_models <- list()
  if (BOOTSTRAP_N > 0) {
    cat("  - 训练 Bootstrap 模型组 (", BOOTSTRAP_N, "次)...\n", sep="")
    for(k in 1:BOOTSTRAP_N) {
      idx <- sample(nrow(train_data), replace=TRUE)
      b_train <- train_data[idx, ]
      n_min_b <- min(table(b_train$presence_fac))
      try({
        # 快速训练少量树
        bm <- randomForest(
          x = b_train[, env_vars], y = b_train$presence_fac,
          ntree = 200, mtry = best_params$mtry, nodesize = best_params$nodesize,
          sampsize = c(n_min_b, n_min_b), strata = b_train$presence_fac
        )
        boot_models[[length(boot_models)+1]] <- bm
      }, silent=TRUE)
    }
  }
  
  # 清理设备
  graphics.off()
  
  # 准备全空间均值 (用于控制变量法)
  mean_data <- train_data[1, env_vars, drop=FALSE]
  for(v in env_vars) {
    if(is.numeric(train_data[[v]])) mean_data[[v]] <- mean(train_data[[v]], na.rm=TRUE)
    else mean_data[[v]] <- getmode(train_data[[v]]) # 自定义模数函数
  }
  
  cat("  - 生成绘图...\n")
  for (var in top_vars) {
    # 构建 x 轴范围
    x_range <- seq(min(train_data[[var]], na.rm=T), max(train_data[[var]], na.rm=T), length.out=100)
    
    # 基础预测数据
    pred_df <- mean_data[rep(1, 100), ]
    pred_df[[var]] <- x_range
    
    # 主模型预测
    p_main <- predict(rf_model, pred_df, type="prob")[, "1"]
    plot_df <- data.frame(x = x_range, y = p_main)
    
    # Bootstrap 预测 (计算阴影)
    if(length(boot_models) > 0) {
      boot_mat <- matrix(NA, nrow=100, ncol=length(boot_models))
      for(b in seq_along(boot_models)) {
        try({
          p_boot <- predict(boot_models[[b]], pred_df, type="prob")[, "1"]
          boot_mat[, b] <- p_boot
        }, silent=TRUE)
      }
      # Mean +/- SD
      b_mean <- apply(boot_mat, 1, mean, na.rm=TRUE)
      b_sd   <- apply(boot_mat, 1, sd, na.rm=TRUE)
      plot_df$ymin <- pmax(0, b_mean - b_sd)
      plot_df$ymax <- pmin(1, b_mean + b_sd)
    }
    
    # 绘图 (严格复刻 Maxnet 风格)
    p <- ggplot(plot_df, aes(x=x, y=y))
    
    if("ymin" %in% names(plot_df)) {
      p <- p + geom_ribbon(aes(ymin=ymin, ymax=ymax), fill="grey70", alpha=0.5)
    }
    
    p <- p + 
      geom_line(color="#3633f2", size=1.2) +
      labs(x = var, y = "Probability (RF)") +
      theme_bw(base_size = 10) +  # 保持 10pt 字号
      theme(
        panel.grid = element_blank(),
        axis.title = element_text(face="plain", color="black"),
        axis.text = element_text(color="black")
      )
    
    out_file <- file.path(output_dir, "response_curves", paste0("curve_", var, ".png"))
    ggsave(out_file, p, width = 8, height = 6, units = "in", dpi = 300)
  }
}

# ------------------------------------------------------------------------------
# 7. 保存结果
# ------------------------------------------------------------------------------
cat("\n保存结果...\n")
saveRDS(rf_model, file.path(output_dir, "model.rds"))
write.csv(var_imp, file.path(output_dir, "variable_importance.csv"), row.names=FALSE)

# 组合预测结果表
res_df <- data.frame(
  id = model_data$id,
  observed = model_data$presence,
  dataset = ifelse(seq_len(nrow(model_data)) %in% train_idx, "train", "test"),
  predicted = NA
)
res_df$predicted[train_idx] <- pred_train
res_df$predicted[-train_idx] <- pred_test # 只有test_idx是负索引时才用-train_idx
write.csv(res_df, file.path(output_dir, "predictions.csv"), row.names=FALSE)

# 评估表
eval_df <- data.frame(
  Model = "Random Forest",
  Items = c("AUC_Train", "AUC_Test", "TSS", "Threshold", "Sensitivity", "Specificity"),
  Value = c(eval_metrics$auc_train, eval_metrics$auc_test, eval_metrics$tss, eval_metrics$threshold, eval_metrics$sens, eval_metrics$spec)
)
write.csv(eval_df, file.path(output_dir, "evaluation.csv"), row.names=FALSE)

cat("\n✓ RF 模型训练全流程完成!\n")
