#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 05b_model_nn.R
# 功能说明: 使用多层感知机(MLP)进行物种分布建模
# 方法: nnet包 (单隐藏层), 概率输出
# 输入文件: output/04_collinearity/collinearity_removed.csv
# 输出文件: output/05b_model_nn/model.rds, predictions.csv, evaluation.csv
# 作者: Nature级别科研项目
# 日期: 2025-10-20
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/SDM01")

packages <- c("tidyverse", "nnet", "pROC")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

dir.create("output/05b_model_nn", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("神经网络(MLP) 模型训练 (83变量)\n")
cat("======================================\n\n")

# 1. 数据准备与标准化
cat("步骤 1/4: 数据准备与标准化...\n")
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

# 标准化（基于训练集统计量）
scale_mean <- sapply(train[, env_vars, drop = FALSE], mean, na.rm = TRUE)
scale_sd   <- sapply(train[, env_vars, drop = FALSE], sd, na.rm = TRUE)
scale_sd[scale_sd == 0 | is.na(scale_sd)] <- 1

scale_apply <- function(df) {
  vals <- sweep(as.matrix(df[, env_vars, drop = FALSE]), 2, scale_mean, "-")
  vals <- sweep(vals, 2, scale_sd, "/")
  as.data.frame(vals)
}

X_train <- scale_apply(train)
X_test  <- scale_apply(test)
y_train <- train$presence
y_test  <- test$presence

# 2. 训练神经网络
cat("\n步骤 2/4: 训练神经网络...\n")
set.seed(123)
size_hidden <- max(3, floor(length(env_vars) / 5))  # 83变量 -> 16个隐藏单元
# 权重数 = (83+1)*16 + (16+1)*1 = 1361，需要设置MaxNWts
max_weights <- (length(env_vars) + 1) * size_hidden + (size_hidden + 1) * 1 + 100
cat("  - 隐藏单元数: ", size_hidden, ", 权重数: ", max_weights - 100, "\n", sep = "")

start_time <- Sys.time()
nn_model <- nnet(x = X_train, y = y_train, size = size_hidden,
                 linout = FALSE, rang = 0.1, decay = 5e-4, maxit = 500,
                 MaxNWts = max_weights, trace = FALSE)
train_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
cat("  ✓ 完成 (", round(train_time, 2), " 秒)\n", sep = "")

# 3. 预测与评估
cat("\n步骤 3/4: 预测与评估...\n")
pred_train <- predict(nn_model, X_train, type = "raw")
pred_test  <- predict(nn_model, X_test, type = "raw")

roc_train <- roc(y_train, pred_train, quiet = TRUE)
roc_test  <- roc(y_test, pred_test, quiet = TRUE)
auc_train <- as.numeric(auc(roc_train))
auc_test  <- as.numeric(auc(roc_test))

coords_result <- coords(roc_test, "best", ret = "all", best.method = "youden")
thr <- coords_result$threshold
tss <- coords_result$sensitivity + coords_result$specificity - 1

cat("  - AUC: ", round(auc_train, 4), " (train) / ", round(auc_test, 4), " (test)\n", sep = "")
cat("  - TSS: ", round(tss, 4), ", 阈值: ", round(thr, 4), "\n", sep = "")

# 4. 保存结果
cat("\n步骤 4/4: 保存结果...\n")
saveRDS(list(model = nn_model, mean = scale_mean, sd = scale_sd, vars = env_vars),
        "output/05b_model_nn/model.rds")

predictions <- data.frame(
  id = model_data$id, species = model_data$species,
  lon = model_data$lon, lat = model_data$lat, presence = model_data$presence,
  dataset = ifelse(seq_len(nrow(model_data)) %in% train_idx, "train", "test"),
  predicted = NA_real_
)
predictions$predicted[train_idx] <- pred_train
predictions$predicted[test_idx]  <- pred_test
write.csv(predictions, "output/05b_model_nn/predictions.csv", row.names = FALSE)

evaluation <- data.frame(
  model = "NN", dataset = c("train", "test"),
  n_samples = c(length(y_train), length(y_test)),
  n_presence = c(sum(y_train == 1), sum(y_test == 1)),
  AUC = c(auc_train, auc_test), TSS = c(NA, tss),
  Sensitivity = c(NA, coords_result$sensitivity),
  Specificity = c(NA, coords_result$specificity),
  optimal_threshold = c(NA, thr), training_time_sec = c(train_time, NA)
)
write.csv(evaluation, "output/05b_model_nn/evaluation.csv", row.names = FALSE)

# 日志
sink("output/05b_model_nn/processing_log.txt")
cat("NN 模型训练日志\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
cat("训练/测试: ", nrow(train), " / ", nrow(test), "\n", sep = "")
cat("隐藏单元: ", size_hidden, "\n", sep = "")
cat("AUC: ", round(auc_train, 4), " / ", round(auc_test, 4), "\n", sep = "")
cat("TSS: ", round(tss, 4), ", 阈值: ", round(thr, 4), "\n", sep = "")
sink()

cat("  ✓ 模型: output/05b_model_nn/model.rds\n")
cat("  ✓ 预测: output/05b_model_nn/predictions.csv\n")
cat("  ✓ 评估: output/05b_model_nn/evaluation.csv\n")
cat("\n✓ 脚本执行完成!\n\n")


