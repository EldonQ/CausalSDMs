#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 05_model_maxnet.R
# 功能说明: 使用Maxnet进行物种分布建模（Presence-Background）
# 方法: Maxent (maxnet包, 特征类别 lqph)
# 输入文件: output/04_collinearity/collinearity_removed.csv
# 输出文件: output/05_model_maxnet/model.rds, predictions.csv, evaluation.csv, variable_importance.csv
# 作者: Nature级别科研项目
# 日期: 2025-10-20
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/SDM01")

# 加载必要的包
packages <- c("tidyverse", "maxnet", "pROC")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

dir.create("output/05_model_maxnet", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("Maxnet 模型训练 (83变量)\n")
cat("======================================\n\n")

# 1. 读取数据并划分训练/测试集
cat("步骤 1/5: 数据准备...\n")
model_data <- read.csv("output/04_collinearity/collinearity_removed.csv")
exclude_cols <- c("id", "species", "lon", "lat", "source", "presence", "presence.1")
env_vars <- setdiff(names(model_data), exclude_cols)

cat("  - 样本: ", nrow(model_data), " (出现点: ", sum(model_data$presence == 1), 
    ", 背景点: ", sum(model_data$presence == 0), ")\n", sep = "")
cat("  - 变量数: ", length(env_vars), "\n", sep = "")

# 分层划分
set.seed(12345)
presence_idx <- which(model_data$presence == 1)
background_idx <- which(model_data$presence == 0)
train_idx <- c(sample(presence_idx, round(0.8 * length(presence_idx))),
               sample(background_idx, round(0.8 * length(background_idx))))
test_idx <- setdiff(seq_len(nrow(model_data)), train_idx)

train_data <- model_data[train_idx, ]
test_data  <- model_data[test_idx, ]
cat("  - 训练/测试: ", nrow(train_data), " / ", nrow(test_data), "\n", sep = "")

# 2. 训练Maxnet模型
cat("\n步骤 2/5: 训练Maxnet模型...\n")
X_train <- train_data[, env_vars, drop = FALSE]  # 保持数据框格式
P_train <- as.numeric(train_data$presence == 1)
fml <- maxnet.formula(P_train, X_train, classes = "lqph")

start_time <- Sys.time()
maxnet_model <- maxnet(p = P_train, data = X_train, f = fml)
train_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat("  ✓ 完成 (", round(train_time, 2), " 秒)\n", sep = "")

# 3. 预测与评估
cat("\n步骤 3/5: 预测与评估...\n")
pred_train <- predict(maxnet_model, X_train, type = "logistic")
pred_test  <- predict(maxnet_model, test_data[, env_vars, drop = FALSE], type = "logistic")

roc_train <- roc(train_data$presence, pred_train, quiet = TRUE)
roc_test  <- roc(test_data$presence, pred_test, quiet = TRUE)
auc_train <- as.numeric(auc(roc_train))
auc_test  <- as.numeric(auc(roc_test))

coords_result <- coords(roc_test, "best", ret = "all", best.method = "youden")
thr <- coords_result$threshold
tss <- coords_result$sensitivity + coords_result$specificity - 1

cat("  - AUC: ", round(auc_train, 4), " (train) / ", round(auc_test, 4), " (test)\n", sep = "")
cat("  - TSS: ", round(tss, 4), ", 阈值: ", round(thr, 4), "\n", sep = "")

# 4. 变量重要性（排列法）
cat("\n步骤 4/5: 变量重要性 (排列法, 83变量)...\n")
baseline_auc <- auc_test
var_importance <- data.frame(variable = env_vars, importance = NA_real_)

set.seed(123)
for(i in seq_along(env_vars)) {
  perm_test <- test_data
  perm_test[[env_vars[i]]] <- sample(perm_test[[env_vars[i]]])
  pred_perm <- predict(maxnet_model, perm_test[, env_vars, drop = FALSE], type = "logistic")
  var_importance$importance[i] <- baseline_auc - as.numeric(auc(roc(perm_test$presence, pred_perm, quiet = TRUE)))
  if(i %% 20 == 0) cat("  - 进度: ", i, " / ", length(env_vars), "\n", sep = "")
}
var_importance <- var_importance[order(var_importance$importance, decreasing = TRUE), ]

# 5. 保存结果
cat("\n步骤 5/5: 保存结果...\n")
saveRDS(maxnet_model, "output/05_model_maxnet/model.rds")

predictions <- data.frame(
  id = model_data$id, species = model_data$species,
  lon = model_data$lon, lat = model_data$lat, presence = model_data$presence,
  dataset = ifelse(seq_len(nrow(model_data)) %in% train_idx, "train", "test"),
  predicted = NA_real_
)
predictions$predicted[train_idx] <- pred_train
predictions$predicted[test_idx]  <- pred_test
write.csv(predictions, "output/05_model_maxnet/predictions.csv", row.names = FALSE)

evaluation <- data.frame(
  model = "Maxnet", dataset = c("train", "test"),
  n_samples = c(nrow(train_data), nrow(test_data)),
  n_presence = c(sum(train_data$presence == 1), sum(test_data$presence == 1)),
  AUC = c(auc_train, auc_test), TSS = c(NA, tss),
  Sensitivity = c(NA, coords_result$sensitivity),
  Specificity = c(NA, coords_result$specificity),
  optimal_threshold = c(NA, thr), training_time_sec = c(train_time, NA)
)
write.csv(evaluation, "output/05_model_maxnet/evaluation.csv", row.names = FALSE)
write.csv(var_importance, "output/05_model_maxnet/variable_importance.csv", row.names = FALSE)

# 日志
sink("output/05_model_maxnet/processing_log.txt")
cat("Maxnet 模型训练日志\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
cat("训练/测试: ", nrow(train_data), " / ", nrow(test_data), "\n", sep = "")
cat("AUC: ", round(auc_train, 4), " / ", round(auc_test, 4), "\n", sep = "")
cat("TSS: ", round(tss, 4), ", 阈值: ", round(thr, 4), "\n\n", sep = "")
cat("变量重要性(前15):\n")
print(head(var_importance, 15))
sink()

cat("  ✓ 模型: output/05_model_maxnet/model.rds\n")
cat("  ✓ 预测: output/05_model_maxnet/predictions.csv\n")
cat("  ✓ 评估: output/05_model_maxnet/evaluation.csv\n")
cat("  ✓ 变量重要性: output/05_model_maxnet/variable_importance.csv\n")
cat("\n✓ 脚本执行完成!\n\n")


