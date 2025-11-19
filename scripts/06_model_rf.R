#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 06_model_rf.R
# 功能说明: 使用随机森林(Random Forest)进行物种分布建模
# 方法: randomForest (分类), 类别不平衡通过采样权重/分层划分缓解
# 输入文件: output/04_collinearity/collinearity_removed.csv
# 输出文件: output/06_model_rf/model.rds, predictions.csv, evaluation.csv, variable_importance.csv
# 作者: Nature级别科研项目
# 日期: 2025-10-20
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/SDM01")

packages <- c("tidyverse", "randomForest", "pROC")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

dir.create("output/06_model_rf", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("随机森林 模型训练 (83变量)\n")
cat("======================================\n\n")

# 1. 数据准备
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

train <- model_data[train_idx, ]
test  <- model_data[test_idx, ]
cat("  - 训练/测试: ", nrow(train), " / ", nrow(test), "\n", sep = "")

# 2. 训练随机森林
cat("\n步骤 2/5: 训练随机森林...\n")
X_train <- train[, env_vars, drop = FALSE]
y_train <- factor(train$presence)

# 平衡类别采样
class_weights <- table(y_train)
sampl_size <- c(`0` = min(class_weights), `1` = min(class_weights))
mtry_val <- max(1, floor(sqrt(length(env_vars))))  # 83变量 -> mtry=9

cat("  - ntree: 500, mtry: ", mtry_val, "\n", sep = "")

start_time <- Sys.time()
rf_model <- randomForest(x = X_train, y = y_train,
                         ntree = 500, mtry = mtry_val,
                         strata = y_train, sampsize = sampl_size,
                         importance = TRUE, na.action = na.omit)
train_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat("  ✓ 完成 (", round(train_time, 2), " 秒)\n", sep = "")

# 3. 预测与评估
cat("\n步骤 3/5: 预测与评估...\n")
pred_train <- predict(rf_model, newdata = train[, env_vars, drop = FALSE], type = "prob")[, "1"]
pred_test  <- predict(rf_model, newdata = test[, env_vars, drop = FALSE], type = "prob")[, "1"]

roc_train <- roc(train$presence, pred_train, quiet = TRUE)
roc_test  <- roc(test$presence, pred_test, quiet = TRUE)
auc_train <- as.numeric(auc(roc_train))
auc_test  <- as.numeric(auc(roc_test))

coords_result <- coords(roc_test, "best", ret = "all", best.method = "youden")
thr <- coords_result$threshold
tss <- coords_result$sensitivity + coords_result$specificity - 1

cat("  - AUC: ", round(auc_train, 4), " (train) / ", round(auc_test, 4), " (test)\n", sep = "")
cat("  - TSS: ", round(tss, 4), ", 阈值: ", round(thr, 4), "\n", sep = "")

# 4. 变量重要性
cat("\n步骤 4/5: 变量重要性 (Gini)...\n")
imp <- importance(rf_model, type = 2)
imp_df <- data.frame(variable = rownames(imp), importance = imp[, 1], row.names = NULL)
imp_df <- imp_df[order(imp_df$importance, decreasing = TRUE), ]

# 5. 保存结果
cat("\n步骤 5/5: 保存结果...\n")
saveRDS(rf_model, "output/06_model_rf/model.rds")

predictions <- data.frame(
  id = model_data$id, species = model_data$species,
  lon = model_data$lon, lat = model_data$lat, presence = model_data$presence,
  dataset = ifelse(seq_len(nrow(model_data)) %in% train_idx, "train", "test"),
  predicted = NA_real_
)
predictions$predicted[train_idx] <- pred_train
predictions$predicted[test_idx]  <- pred_test
write.csv(predictions, "output/06_model_rf/predictions.csv", row.names = FALSE)

evaluation <- data.frame(
  model = "RF", dataset = c("train", "test"),
  n_samples = c(nrow(train), nrow(test)),
  n_presence = c(sum(train$presence == 1), sum(test$presence == 1)),
  AUC = c(auc_train, auc_test), TSS = c(NA, tss),
  Sensitivity = c(NA, coords_result$sensitivity),
  Specificity = c(NA, coords_result$specificity),
  optimal_threshold = c(NA, thr), training_time_sec = c(train_time, NA)
)
write.csv(evaluation, "output/06_model_rf/evaluation.csv", row.names = FALSE)
write.csv(imp_df, "output/06_model_rf/variable_importance.csv", row.names = FALSE)

# 日志
sink("output/06_model_rf/processing_log.txt")
cat("RF 模型训练日志\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
cat("训练/测试: ", nrow(train), " / ", nrow(test), "\n", sep = "")
cat("ntree: 500, mtry: ", mtry_val, "\n", sep = "")
cat("AUC: ", round(auc_train, 4), " / ", round(auc_test, 4), "\n", sep = "")
cat("TSS: ", round(tss, 4), ", 阈值: ", round(thr, 4), "\n\n", sep = "")
cat("变量重要性(前15):\n")
print(head(imp_df, 15))
sink()

cat("  ✓ 模型: output/06_model_rf/model.rds\n")
cat("  ✓ 预测: output/06_model_rf/predictions.csv\n")
cat("  ✓ 评估: output/06_model_rf/evaluation.csv\n")
cat("  ✓ 变量重要性: output/06_model_rf/variable_importance.csv\n")
cat("\n✓ 脚本执行完成!\n\n")


