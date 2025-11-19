#!/usr/bin/env Rscript
# 计算GAM模型的平均预测概率

# 读取数据
gam_pred <- read.csv("output/07_model_gam/predictions.csv")
data <- read.csv("output/04_collinearity/collinearity_removed.csv")

# 计算统计量
cat("\n=== 数据集基本信息 ===\n")
cat("总样本数:", nrow(data), "\n")
cat("存在点数:", sum(data$presence == 1), "\n")
cat("缺失点数:", sum(data$presence == 0), "\n")
cat("存在比例:", round(mean(data$presence), 4), "\n")

cat("\n=== GAM 模型预测统计 ===\n")
cat("平均预测概率:", round(mean(gam_pred$predicted), 6), "\n")
cat("中位数预测概率:", round(median(gam_pred$predicted), 6), "\n")
cat("预测概率范围:", round(min(gam_pred$predicted), 6), "到", round(max(gam_pred$predicted), 6), "\n")
cat("标准差:", round(sd(gam_pred$predicted), 6), "\n")

# 按 presence 分组统计
cat("\n=== 按实际标签分组 ===\n")
presence_1 <- gam_pred[gam_pred$presence == 1, ]
presence_0 <- gam_pred[gam_pred$presence == 0, ]
cat("存在点的平均预测概率:", round(mean(presence_1$predicted), 6), "\n")
cat("缺失点的平均预测概率:", round(mean(presence_0$predicted), 6), "\n")

cat("\n=== ALE 基线值验证 ===\n")
cat("从 ALE CSV 读取的基线值约为: 0.254\n")
cat("实际计算的平均预测概率:", round(mean(gam_pred$predicted), 6), "\n")
cat("差异:", round(abs(mean(gam_pred$predicted) - 0.254), 6), "\n")










