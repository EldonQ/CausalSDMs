#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 09_variable_importance_viz.R (Advanced Aesthetic)
# 功能说明: 绘制 Nature 级变量重要性图 (紫色渐变 + 数值标注 + 极简风格)
# 输入文件: output/05_model_maxnet/variable_importance.csv
#          output/06_model_rf/variable_importance.csv
# 输出: figures/09_variable_importance/MaxnetRF/feature_importance.png
# ==============================================================================

# 初始化
rm(list = ls()); gc()
setwd("E:/CausalSDMs")

packages <- c("tidyverse", "ggplot2", "scales", "gridExtra")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 目录设置
out_dir <- "output/09_variable_importance/MaxnetRF"
fig_dir <- "figures/09_variable_importance/MaxnetRF"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n==================================================\n")
cat("      变量重要性可视化 (Nature Aesthetic)      \n")
cat("==================================================\n\n")

# 1. 数据读取与预处理
# ------------------------------------------------------------------------------
cat("步骤 1/3: 读取与标准化数据...\n")

# 读取 Maxnet (Permutation Importance)
imp_maxnet <- read.csv("output/05_model_maxnet/variable_importance.csv") %>%
  dplyr::select(variable, importance = importance_pct) %>% # 假设有一列百分比，或者用原始值归一化
  mutate(Model = "Maxnet") %>%
  arrange(desc(importance)) %>%
  mutate(importance = importance / max(importance) * 100) # 归一化到 0-100 相对值

# 读取 RF (MeanDecreaseAccuracy or Gini)
imp_rf <- read.csv("output/06_model_rf/variable_importance.csv") %>%
  mutate(Model = "Random Forest") %>%
  arrange(desc(importance)) %>%
  mutate(importance = importance / max(importance) * 100) # 归一化到 0-100 相对值

# 合并
plot_data <- bind_rows(imp_maxnet, imp_rf)

# 2. 绘图函数 (极简大气风)
# ------------------------------------------------------------------------------
create_importance_plot <- function(data, model_name, top_n = 20) {
  
  # 筛选前N个变量，并重新排序因子
  df_sub <- data %>%
    filter(Model == model_name) %>%
    arrange(desc(importance)) %>%
    head(top_n) %>%
    mutate(variable = factor(variable, levels = rev(variable))) # 反转因子顺序使最高的在上面
  
  p <- ggplot(df_sub, aes(x = importance, y = variable, fill = importance)) +
    geom_col(width = 0.8) +
    # 紫色渐变风格 (参考用户附图)
    scale_fill_gradient(low = "#D8BFD8", high = "#800080") + 
    # 添加数值标签 (白色，加粗，位于条形图内部右侧)
    geom_text(aes(label = sprintf("%.1f", importance)), 
              hjust = 1.1, color = "white", size = 3.5, fontface = "bold") +
    scale_x_continuous(expand = c(0, 0)) +
    labs(title = paste0("Feature Importance (", model_name, ")"),
         x = "Relative Importance (%)", y = NULL) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.text.y = element_text(color = "black", size = 12),
      axis.text.x = element_text(color = "black", size = 10),
      axis.line.y = element_blank(), # 移除y轴线
      axis.ticks.y = element_blank(), # 移除y轴刻度
      legend.position = "none", # 移除图例 (颜色本身已代表数值)
      panel.grid.major.x = element_line(color = "grey90", linetype = "dotted"), # 仅保留横向网格
      plot.margin = margin(10, 15, 10, 5)
    )
  
  return(p)
}

# 3. 生成并保存图表
# ------------------------------------------------------------------------------
cat("\n步骤 2/3: 生成图表...\n")

p_maxnet <- create_importance_plot(plot_data, "Maxnet", top_n = 25)
p_rf <- create_importance_plot(plot_data, "Random Forest", top_n = 25)

# 单独保存
ggsave(file.path(fig_dir, "importance_maxnet.png"), p_maxnet, width = 6, height = 8, dpi = 300)
ggsave(file.path(fig_dir, "importance_rf.png"), p_rf, width = 6, height = 8, dpi = 300)

# 组合保存 (并排)
p_combined <- grid.arrange(p_maxnet, p_rf, ncol = 2)
ggsave(file.path(fig_dir, "importance_combined.png"), p_combined, width = 12, height = 8, dpi = 300)

cat("  ✓ Maxnet图: figures/09_variable_importance/MaxnetRF/importance_maxnet.png\n")
cat("  ✓ RF图:     figures/09_variable_importance/MaxnetRF/importance_rf.png\n")
cat("  ✓ 组合图:   figures/09_variable_importance/MaxnetRF/importance_combined.png\n")

# 保存汇总数据
write.csv(plot_data, file.path(out_dir, "importance_summary_all.csv"), row.names = FALSE)

cat("\n✓ 脚本执行完成!\n")
