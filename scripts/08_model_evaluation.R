#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 08_model_evaluation.R (Advanced Edition)
# 功能说明: Nature级模型综合评估系统 (Maxnet vs RF)
# 特性:
#   1. 多维度指标: AUC, TSS, Kappa, F1-score, Boyce Index(Approx)
#   2. 统计检验: DeLong Test (AUC差异显著性)
#   3. 高级绘图: ROC, 密度分布图(Density), 校准曲线(Calibration)
#   4. 输出管理: 自动归档至 MaxnetRF 子目录
# ==============================================================================

# 初始化
rm(list = ls()); gc()
setwd("E:/CausalSDMs")

# 加载包
packages <- c("tidyverse", "pROC", "ggplot2", "gridExtra", "scales", "caret")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 定义输出目录
out_dir <- "output/08_model_evaluation/MaxnetRF"
fig_dir <- "figures/08_model_evaluation/MaxnetRF"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n==================================================\n")
cat("      高级模型评估系统 (Maxnet vs RF)      \n")
cat("==================================================\n\n")

# ------------------------------------------------------------------------------
# 辅助函数: 计算 Boyce Index (基于分箱的简易版)
# ------------------------------------------------------------------------------
calc_boyce_binned <- function(fit, obs, n_bins = 10) {
  # fit: 预测值, obs: 观测值 (1/0)
  # 仅使用 Presence 点的预测值分布与背景分布对比
  # 这里简化为: 预测值分箱后，计算 P/E 比率的 Spearman 相关系数
  
  if(sum(obs==1) == 0) return(NA)
  
  # 过滤
  pres_pred <- fit[obs == 1]
  bg_pred <- fit[obs == 0] # 或者是全部点，视数据类型而定，这里假设0是背景
  
  # 分箱
  bins <- quantile(fit, probs = seq(0, 1, length.out = n_bins + 1), na.rm=TRUE)
  # 避免重复断点
  if(any(duplicated(bins))) bins <- seq(min(fit), max(fit), length.out = n_bins + 1)
  
  # 计算每个箱内的点数
  p_counts <- hist(pres_pred, breaks = bins, plot = FALSE)$counts
  b_counts <- hist(bg_pred, breaks = bins, plot = FALSE)$counts
  
  # 频率
  p_freq <- p_counts / sum(p_counts)
  b_freq <- b_counts / sum(b_counts)
  
  # P/E ratio
  pe_ratio <- p_freq / (b_freq + 1e-10) # 避免除0
  
  # 取箱的中点作为预测值等级
  bin_mids <- (bins[-1] + bins[-(n_bins+1)]) / 2
  
  # 计算 Spearman 相关
  cor_val <- cor(bin_mids[b_counts > 0], pe_ratio[b_counts > 0], method = "spearman")
  return(cor_val)
}

# ------------------------------------------------------------------------------
# 1. 数据读取与预处理
# ------------------------------------------------------------------------------
cat("步骤 1/5: 读取并标准化数据...\n")

models <- c("Maxnet", "RF")
model_dirs <- c("05_model_maxnet", "06_model_rf")
pred_list <- list()
eval_df_list <- list()

for (i in seq_along(models)) {
  m_name <- models[i]
  p_file <- paste0("output/", model_dirs[i], "/predictions.csv")
  e_file <- paste0("output/", model_dirs[i], "/evaluation.csv")
  
  # 读取预测 (关键)
  if (file.exists(p_file)) {
    tmp <- read.csv(p_file) %>% filter(dataset == "test") # 只看测试集
    
    # 列名标准化
    if ("observed" %in% names(tmp)) tmp <- rename(tmp, presence = observed)
    
    # 确保 predicted 列存在且有效
    if ("predicted" %in% names(tmp)) {
      tmp <- tmp %>% filter(!is.na(predicted), !is.na(presence))
      pred_list[[m_name]] <- tmp
      cat("  ✓ ", m_name, ": 读取测试样本 ", nrow(tmp), " 条\n", sep="")
    }
  }
  
  # 读取预测 (关键)
  if (file.exists(p_file)) {
    tmp <- read.csv(p_file) %>% filter(dataset == "test") # 只看测试集
    
    # 列名标准化
    if ("observed" %in% names(tmp)) tmp <- rename(tmp, presence = observed)
    
    # 确保 predicted 列存在且有效
    if ("predicted" %in% names(tmp)) {
      tmp <- tmp %>% filter(!is.na(predicted), !is.na(presence))
      pred_list[[m_name]] <- tmp
      cat("  ✓ ", m_name, ": 读取测试样本 ", nrow(tmp), " 条\n", sep="")
    }
  }
}

# ------------------------------------------------------------------------------
# 2. 高级指标计算与统计检验
# ------------------------------------------------------------------------------
cat("\n步骤 2/5: 计算高级指标与统计检验...\n")

advanced_metrics <- data.frame()
roc_objects <- list()

for (m_name in names(pred_list)) {
  dat <- pred_list[[m_name]]
  
  # 1. 基础 ROC 对象
  roc_obj <- roc(dat$presence, dat$predicted, quiet = TRUE)
  roc_objects[[m_name]] <- roc_obj
  
  # 2. 确定最优阈值 (Youden Index)
  coords_best <- coords(roc_obj, "best", best.method = "youden", ret = c("threshold", "specificity", "sensitivity"))
  # 处理可能返回多个阈值的情况
  if(is.matrix(coords_best) || is.data.frame(coords_best)) coords_best <- coords_best[1, ]
  thresh <- coords_best$threshold
  
  # 3. 基于阈值的混淆矩阵指标
  pred_class <- ifelse(dat$predicted >= thresh, 1, 0)
  cm <- confusionMatrix(factor(pred_class, levels=c(0,1)), factor(dat$presence, levels=c(0,1)))
  
  # 4. Boyce Index
  boyce <- calc_boyce_binned(dat$predicted, dat$presence)
  
  # 汇总
  metrics_row <- data.frame(
    Model = m_name,
    AUC = as.numeric(auc(roc_obj)),
    TSS = as.numeric(coords_best$sensitivity + coords_best$specificity - 1),
    Kappa = as.numeric(cm$overall["Kappa"]),
    F1 = as.numeric(cm$byClass["F1"]),
    Accuracy = as.numeric(cm$overall["Accuracy"]),
    Boyce = boyce,
    Optimal_Threshold = thresh
  )
  advanced_metrics <- rbind(advanced_metrics, metrics_row)
}

# DeLong Test (两两比较)
cat("  - 执行 DeLong Test (AUC 差异显著性检验)...\n")
delong_res <- "N/A"
if (length(roc_objects) == 2) {
  test_res <- roc.test(roc_objects[[1]], roc_objects[[2]], method = "delong")
  p_val <- test_res$p.value
  signif_symbol <- ifelse(p_val < 0.001, "***", ifelse(p_val < 0.01, "**", ifelse(p_val < 0.05, "*", "ns")))
  delong_res <- sprintf("p = %.4f (%s)", p_val, signif_symbol)
  cat(sprintf("    对比 %s vs %s: %s\n", names(roc_objects)[1], names(roc_objects)[2], delong_res))
  
  # 添加到结果表
  advanced_metrics$Difference_Significance <- c(delong_res, "")
}

write.csv(advanced_metrics, file.path(out_dir, "advanced_evaluation_metrics.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 3. 绘图 A: ROC 曲线 (含 DeLong P值标注)
# ------------------------------------------------------------------------------
cat("\n步骤 3/5: 绘制 ROC 曲线 (修复截断问题)...\n")

png(file.path(fig_dir, "roc_curves_advanced.png"), width = 2400, height = 2400, res = 300, bg = "transparent")
# 简洁风格：减小顶部边距，移除标题 (符合奥卡姆剃刀原理)
par(mar = c(5, 5, 2, 2), family = "serif", cex.lab=1.2, cex.axis=1.1)

col_map <- c("Maxnet" = "#E64B35", "RF" = "#4DBBD5")

# 空画布 - ylim 设为 1.02 (略高于1) 且使用 yaxs="i" 刚性约束，保证 1.0 处的线条完整且不留过多空白
plot(roc_objects[[1]], type = "n", legacy.axes = TRUE,
     xlab = "False Positive Rate (1 - Specificity)",
     ylab = "True Positive Rate (Sensitivity)",
     main = "",
     ylim = c(0, 1.02), xlim = c(1, 0), xaxs="i", yaxs="i")
grid(col = "grey90")

# 移除所有标题函数 (title/mtext)，保持图形纯净

# 参考线
abline(a = 0, b = 1, lty = 2, col = "grey60")

# 循环画线
for (m in names(roc_objects)) {
  plot(roc_objects[[m]], add = TRUE, col = col_map[m], lwd = 3, legacy.axes=TRUE)
}

# 图例 (inset参数微调位置: 第一个数是左右偏移，第二个数是上下偏移)
legend_txt <- paste0(advanced_metrics$Model, " (AUC = ", round(advanced_metrics$AUC, 3), ")")
# 上移图例: inset = c(0.05, 0.08)
legend("bottomright", legend = legend_txt, col = col_map[names(roc_objects)], 
       lwd = 3, bty = "n", cex = 1.1, inset = c(0.05, 0.08))

# 标注 P 值
if (delong_res != "N/A") {
  text(0.5, 0.15, paste("DeLong Test:\n", delong_res), adj = 0, col = "black", font = 3, cex = 1)
}

dev.off()

# ------------------------------------------------------------------------------
# 4. 绘图 B: SDM 专属评估图 (Boyce & TSS 曲线)
# ------------------------------------------------------------------------------
cat("步骤 4/5: 绘制 SDM 专属评估图 (Boyce, TSS)...\n")

# --- 4.1 Boyce Index P/E 曲线 ---
# 这是一个展示模型校准度的黄金指标图
png(file.path(fig_dir, "boyce_pe_curves.png"), width = 2400, height = 1800, res = 300, bg = "white")
par(mar = c(5, 5, 4, 2), family = "serif")
plot(0, 0, type="n", xlim=c(0, 1), ylim=c(0, 5), # 通常P/E比率在0-5之间
     xlab="Predicted Probability", ylab="Predicted-to-Expected Ratio (P/E)",
     main="Boyce Index (Model Calibration)")
grid()
abline(h=1, col="grey", lty=2) # 参考线: P/E=1 表示随机预测

for (m in names(pred_list)) {
  dat <- pred_list[[m]]
  # 计算分箱 P/E
  obs <- dat$presence
  fit <- dat$predicted
  # 简易分箱
  bins <- quantile(fit, probs = seq(0, 1, length.out = 11), na.rm=T)
  digits <- cut(fit, unique(bins), include.lowest=T)
  
  # 计算每个箱的 P/E
  bin_stats <- data.frame(fit=fit, obs=obs, bin=digits) %>%
    group_by(bin) %>%
    summarise(
      mean_pred = mean(fit),
      n_pres = sum(obs==1),
      n_total = n(),
      .groups = 'drop'
    ) %>%
    mutate(
      f_pres = n_pres / sum(n_pres),     # 出现频率
      f_total = n_total / sum(n_total),  # 期望频率
      pe = f_pres / f_total              # P/E Ratio
    ) %>%
    filter(!is.na(pe) & is.finite(pe))
  
  # 画点和线 (平滑)
  points(bin_stats$mean_pred, bin_stats$pe, col=alpha(col_map[m], 0.4), pch=16)
  lines(lowess(bin_stats$mean_pred, bin_stats$pe, f=0.8), col=col_map[m], lwd=3)
}
legend("topleft", legend=names(pred_list), col=col_map, lwd=3, bty="n", cex=1.2)
dev.off()


# --- 4.2 TSS vs Threshold 曲线 ---
# 展示 TSS 随阈值变化的趋势，直观显示最优阈值
tss_data <- data.frame()

for (m in names(roc_objects)) {
  roc_obj <- roc_objects[[m]]
  # 提取所有阈值的指标
  # coords 返回矩阵: threshold, specificity, sensitivity
  res <- coords(roc_obj, "all", ret=c("threshold", "specificity", "sensitivity"), transpose=FALSE)
  res$TSS <- res$sensitivity + res$specificity - 1
  res$Model <- m
  tss_data <- rbind(tss_data, res)
}

p_tss <- ggplot(tss_data, aes(x=threshold, y=TSS, color=Model)) +
  geom_line(linewidth=1.2) +
  scale_color_manual(values=col_map) +
  theme_bw(base_size = 14) +
  labs(title="TSS Sensitivity to Threshold", x="Threshold", y="TSS") +
  theme(panel.grid.minor = element_blank(), legend.position=c(0.85, 0.85))

ggsave(file.path(fig_dir, "tss_threshold_curves.png"), p_tss, width = 6, height = 5, dpi = 300)


# ------------------------------------------------------------------------------
# 5. 绘图 C: 模型性能雷达图/柱状图 (多指标)
# ------------------------------------------------------------------------------
cat("步骤 5/5: 绘制综合指标对比图...\n")

metrics_long <- advanced_metrics %>%
  select(Model, AUC, TSS, Kappa, F1, Accuracy) %>%
  pivot_longer(cols = -Model, names_to = "Metric", values_to = "Value")

p_metrics <- ggplot(metrics_long, aes(x = Metric, y = Value, fill = Model)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6, color="black", size=0.3) +
  scale_fill_manual(values = col_map) +
  # 添加数值标签
  geom_text(aes(label = sprintf("%.2f", Value)), 
            position = position_dodge(width = 0.7), vjust = -0.5, size = 3.5) +
  ylim(0, 1.1) +
  labs(title = "Comprehensive Performance Metrics", x = "", y = "Value") +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    panel.grid.major.x = element_blank(),
    axis.text = element_text(color="black")
  )

ggsave(file.path(fig_dir, "metrics_barplot.png"), p_metrics, width = 8, height = 6, dpi = 300)

# ==============================================================================
# 完成
# ==============================================================================
sink(file.path(out_dir, "evaluation_report.txt"))
cat("高级模型评估报告\n")
cat("Generated at:", format(Sys.time()), "\n\n")
cat("1. 详细指标:\n")
print(advanced_metrics)
cat("\n2. 统计检验 (DeLong):\n")
cat("   Maxnet vs RF P-value:", delong_res, "\n")
sink()

cat("\n所有评估结果已保存至:", out_dir, "\n")
cat("✓ 脚本执行完成!\n")
