#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 08_model_evaluation.R
# 功能说明: 综合评估所有模型性能（Maxnet, NN, RF, GAM）
# 方法: ROC曲线对比、性能指标对比、校准曲线
# 输入文件: output/05_model_maxnet/evaluation.csv, predictions.csv
#          output/05b_model_nn/evaluation.csv, predictions.csv
#          output/06_model_rf/evaluation.csv, predictions.csv
#          output/07_model_gam/evaluation.csv, predictions.csv
# 输出文件: output/08_model_evaluation/evaluation_summary.csv
#          figures/08_model_evaluation/roc_curves.png
#          figures/08_model_evaluation/performance_comparison.png
# 作者: Nature级别科研项目
# 日期: 2025-10-20
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/SDM01")

# 加载必要的包
packages <- c("tidyverse", "pROC", "ggplot2", "gridExtra", "scales")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

dir.create("output/08_model_evaluation", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/08_model_evaluation", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("综合模型评估 (4个模型)\n")
cat("======================================\n\n")

# 1. 读取所有模型的评估结果和预测
cat("步骤 1/4: 读取所有模型结果...\n")

models <- c("Maxnet", "NN", "RF", "GAM")
model_dirs <- c("05_model_maxnet", "05b_model_nn", "06_model_rf", "07_model_gam")

# 读取评估结果
eval_list <- list()
pred_list <- list()

for(i in seq_along(models)) {
  eval_file <- paste0("output/", model_dirs[i], "/evaluation.csv")
  pred_file <- paste0("output/", model_dirs[i], "/predictions.csv")
  
  if(file.exists(eval_file)) {
    eval_list[[models[i]]] <- read.csv(eval_file)
    cat("  ✓ ", models[i], " evaluation\n", sep = "")
  } else {
    cat("  ✗ ", models[i], " evaluation 文件不存在\n", sep = "")
  }
  
  if(file.exists(pred_file)) {
    pred_list[[models[i]]] <- read.csv(pred_file)
  }
}

# 合并评估结果
all_eval <- bind_rows(eval_list)
test_eval <- all_eval %>% filter(dataset == "test")

cat("  - 成功读取 ", length(eval_list), " 个模型\n", sep = "")

# 2. 绘制ROC曲线对比
cat("\n步骤 2/4: 绘制ROC曲线对比...\n")

# 定义颜色方案（Nature配色）
model_colors <- c("Maxnet" = "#E41A1C", "NN" = "#377EB8", "RF" = "#4DAF4A", "GAM" = "#984EA3")

# 计算ROC曲线
roc_list <- list()
for(model in names(pred_list)) {
  pred_data <- pred_list[[model]] %>% filter(dataset == "test")
  if(nrow(pred_data) > 0) {
    roc_list[[model]] <- roc(pred_data$presence, pred_data$predicted, quiet = TRUE)
  }
}

# PNG版本
png("figures/08_model_evaluation/roc_curves.png",
    width = 3000, height = 3000, res = 1200, family = "Arial")
par(mar = c(4, 4, 2, 1))
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "False Positive Rate (1 - Specificity)", 
     ylab = "True Positive Rate (Sensitivity)",
     main = "ROC Curves Comparison", cex.main = 0.9, cex.lab = 0.8, cex.axis = 0.7)
abline(a = 0, b = 1, lty = 3, col = "gray50", lwd = 1)

legend_text <- c()
legend_cols <- c()
for(model in names(roc_list)) {
  plot(roc_list[[model]], col = model_colors[model], lwd = 2, add = TRUE)
  auc_val <- as.numeric(auc(roc_list[[model]]))
  legend_text <- c(legend_text, sprintf("%s (AUC = %.3f)", model, auc_val))
  legend_cols <- c(legend_cols, model_colors[model])
}

legend("bottomright", legend = legend_text, col = legend_cols,
       lwd = 2, cex = 0.7, bg = "white")
dev.off()

cat("  ✓ ROC曲线: figures/08_model_evaluation/roc_curves.png\n")

# 3. 绘制性能指标对比
cat("\n步骤 3/4: 绘制性能指标对比...\n")

# 准备数据
metrics_data <- test_eval %>%
  select(model, AUC, TSS, Sensitivity, Specificity) %>%
  pivot_longer(cols = c(AUC, TSS, Sensitivity, Specificity),
               names_to = "Metric", values_to = "Value")

# PNG版本
png("figures/08_model_evaluation/performance_comparison.png",
    width = 3600, height = 2400, res = 1200, family = "Arial")

p <- ggplot(metrics_data, aes(x = Metric, y = Value, fill = model)) +
  geom_bar(stat = "identity", position = "dodge", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = model_colors) +
  labs(title = "Model Performance Comparison",
       x = "Metric", y = "Value", fill = "Model") +
  ylim(0, 1) +
  theme_minimal(base_size = 8) +
  theme(
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 7, face = "bold"),
    axis.text = element_text(size = 6),
    legend.position = "bottom",
    legend.title = element_text(size = 6, face = "bold"),
    legend.text = element_text(size = 6),
    panel.grid.minor = element_blank()
  )
print(p)
dev.off()

cat("  ✓ 性能对比: figures/08_model_evaluation/performance_comparison.png\n")

# 4. 保存综合评估结果
cat("\n步骤 4/4: 保存综合评估结果...\n")

summary_table <- test_eval %>%
  select(model, n_samples, n_presence, AUC, TSS, Sensitivity, Specificity, optimal_threshold) %>%
  arrange(desc(AUC))

write.csv(summary_table, "output/08_model_evaluation/evaluation_summary.csv", row.names = FALSE)
write.csv(all_eval, "output/08_model_evaluation/evaluation_all.csv", row.names = FALSE)

# 日志
sink("output/08_model_evaluation/processing_log.txt")
cat("综合模型评估日志\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
cat("模型数量: ", length(models), "\n\n", sep = "")
cat("测试集性能排名 (按AUC):\n")
print(summary_table)
sink()

cat("  ✓ 评估摘要: output/08_model_evaluation/evaluation_summary.csv\n")

# 输出摘要
cat("\n======================================\n")
cat("模型评估完成\n")
cat("======================================\n\n")

cat("测试集性能排名 (按AUC):\n")
for(i in 1:nrow(summary_table)) {
  cat(sprintf("  %d. %-8s AUC=%.4f, TSS=%.4f\n", 
              i, summary_table$model[i], summary_table$AUC[i], summary_table$TSS[i]))
}

cat("\n✓ 脚本执行完成!\n\n")
