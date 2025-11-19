#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 17_local_sensitivity_analysis.R
# 功能说明: GAM模型的局部敏感度分析和偏依赖图
# 方法: 1) 数值梯度法计算敏感度 2) 偏依赖图展示边际效应
# 输入文件: output/07_model_gam/model.rds
#          output/04_collinearity/collinearity_removed.csv
#          output/09_variable_importance/importance_summary.csv
# 输出文件: figures/17_local_sensitivity/sensitivity_violin.png
#          figures/17_local_sensitivity/partial_dependence.png
#          output/17_local_sensitivity/sensitivity_summary.csv
# 作者: Nature级别科研项目
# 日期: 2025-10-20
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/SDM01")

# 加载必要的包
packages <- c("tidyverse", "mgcv", "ggplot2", "viridis", "patchwork")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

dir.create("output/17_local_sensitivity", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/17_local_sensitivity", showWarnings = FALSE, recursive = TRUE)

# 统一绘图工具（Arial、1200dpi、Nature风格）
source("scripts/visualization/viz_utils.R")

cat("\n======================================\n")
cat("局部敏感度分析与偏依赖图\n")
cat("======================================\n\n")

# 1. 读取模型和数据
cat("步骤 1/5: 读取GAM模型和数据...\n")

gam_model <- readRDS("output/07_model_gam/model.rds")
model_data <- read.csv("output/04_collinearity/collinearity_removed.csv")
var_importance <- read.csv("output/09_variable_importance/importance_summary.csv") %>%
  dplyr::filter(model == "GAM", variable != "lon,lat")

exclude_cols <- c("id", "species", "lon", "lat", "source", "presence", "presence.1")
env_vars <- setdiff(names(model_data), exclude_cols)

# 选择Top 10最重要的变量
top_vars <- var_importance %>%
  dplyr::arrange(desc(importance_normalized)) %>%
  head(10) %>%
  dplyr::pull(variable)

# 仅保留在建模数据中存在的变量，避免列名不一致导致预测失败
top_vars <- intersect(top_vars, env_vars)
if(length(top_vars) == 0) {
  stop("未找到可用于敏感度/偏依赖的变量（请检查 importance_summary 与建模数据是否一致）")
}

cat("  ✓ GAM模型已加载\n")
cat("  ✓ 环境变量: ", length(env_vars), " 个, 选择Top 10进行分析\n", sep = "")

# 2. 定义空间分组
cat("\n步骤 2/5: 定义空间分组...\n")

analysis_data <- model_data %>%
  dplyr::mutate(
    lat_zone = cut(lat, 
                   breaks = quantile(lat, probs = c(0, 1/3, 2/3, 1)),
                   labels = c("South", "Central", "North"),
                   include.lowest = TRUE)
  ) %>%
  dplyr::select(tidyselect::all_of(c("lon", "lat", "lat_zone", env_vars)))  # 包含所有环境变量，GAM需要

# 清除用于预测的数据中的缺失，确保 mgcv::predict 正常
analysis_data <- analysis_data[stats::complete.cases(analysis_data[, c("lon","lat", env_vars)]), ]

cat("  - 纬度带分布: 南 ", sum(analysis_data$lat_zone == "South"),
    ", 中 ", sum(analysis_data$lat_zone == "Central"),
    ", 北 ", sum(analysis_data$lat_zone == "North"), "\n", sep = "")

# 3. 计算局部敏感度
cat("\n步骤 3/5: 计算局部敏感度...\n")

# 敏感度计算函数
compute_sensitivity <- function(data, model, var_name, delta = 0.01) {
  data_plus <- data
  data_plus[[var_name]] <- data[[var_name]] * (1 + delta)
  
  pred_original <- predict(model, newdata = data, type = "response")
  pred_perturbed <- predict(model, newdata = data_plus, type = "response")
  
  sensitivity <- (pred_perturbed - pred_original) / (data[[var_name]] * delta + 1e-10)
  return(sensitivity)
}

# 对每个变量计算敏感度
sensitivity_results <- list()

for(var in top_vars) {
  cat("  - ", var, "\n", sep = "")
  
  tryCatch({
    sens <- compute_sensitivity(analysis_data, gam_model, var)
    
    sensitivity_results[[var]] <- data.frame(
      variable = var,
      sensitivity = sens,
      lat_zone = analysis_data$lat_zone
    )
  }, error = function(e) {
    cat("    ✗ 失败\n")
  })
}

all_sensitivity <- dplyr::bind_rows(sensitivity_results)

# 检查是否有有效结果
if(nrow(all_sensitivity) == 0) {
  cat("  ✗ 警告: 没有成功计算的敏感度，跳过后续分析\n")
  cat("\n脚本终止: 敏感度计算失败\n\n")
  quit(status = 1)
}

# 汇总统计
sensitivity_summary <- all_sensitivity %>%
  dplyr::group_by(variable, lat_zone) %>%
  dplyr::summarise(
    mean_sensitivity = mean(sensitivity, na.rm = TRUE),
    sd_sensitivity = sd(sensitivity, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

write.csv(sensitivity_summary, "output/17_local_sensitivity/sensitivity_summary.csv", 
          row.names = FALSE)
cat("  ✓ 已保存: output/17_local_sensitivity/sensitivity_summary.csv\n")

# 4. 计算偏依赖
cat("\n步骤 4/5: 计算偏依赖（边际效应）...\n")

# 偏依赖计算函数
compute_partial_dependence <- function(data, model, var_name, n_points = 100) {
  # 创建变量范围
  var_range <- seq(
    from = quantile(data[[var_name]], 0.01, na.rm = TRUE),
    to = quantile(data[[var_name]], 0.99, na.rm = TRUE),
    length.out = n_points
  )
  
  # 对每个值计算平均预测
  pd_values <- numeric(n_points)
  
  for(i in seq_along(var_range)) {
    pred_data <- data
    pred_data[[var_name]] <- var_range[i]
    preds <- predict(model, newdata = pred_data, type = "response")
    pd_values[i] <- mean(preds, na.rm = TRUE)
  }
  
  return(data.frame(
    variable = var_name,
    x = var_range,
    y = pd_values
  ))
}

# 对每个变量计算偏依赖
pd_results <- list()

for(var in top_vars) {
  cat("  - ", var, "\n", sep = "")
  
  tryCatch({
    pd_results[[var]] <- compute_partial_dependence(analysis_data, gam_model, var)
  }, error = function(e) {
    cat("    ✗ 失败\n")
  })
}

all_pd <- dplyr::bind_rows(pd_results)

write.csv(all_pd, "output/17_local_sensitivity/partial_dependence_data.csv", 
          row.names = FALSE)
cat("  ✓ 已保存: output/17_local_sensitivity/partial_dependence_data.csv\n")

# 5. 可视化
cat("\n步骤 5/5: 生成可视化图表...\n")

# 5a. 敏感度小提琴图
cat("  (1) 局部敏感度小提琴图...\n")

p_violin <- ggplot(all_sensitivity, aes(x = variable, y = sensitivity, fill = lat_zone)) +
  geom_violin(position = position_dodge(0.8), alpha = 0.7, scale = "width") +
  geom_boxplot(position = position_dodge(0.8), width = 0.15, outlier.size = 0.3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  scale_fill_viridis(discrete = TRUE, option = "viridis", name = "Latitude Zone") +
  labs(title = "Local Sensitivity Analysis by Latitude Zone",
       subtitle = "GAM Model - Top 10 Variables",
       x = "Environmental Variable",
       y = expression("Sensitivity ("*partialdiff*"P/"*partialdiff*"X)")) +
  viz_theme_nature(base_size = 8, title_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 6))

ggsave("figures/17_local_sensitivity/sensitivity_violin.png",
       plot = p_violin, width = 4, height = 3, dpi = 1200)
ggsave("figures/17_local_sensitivity/sensitivity_violin.svg",
       plot = p_violin, width = 4, height = 3)

cat("      ✓ 小提琴图已保存\n")

# 5b. 偏依赖图
cat("  (2) 偏依赖图...\n")

p_pd <- ggplot(all_pd, aes(x = x, y = y)) +
  geom_line(color = "#2166AC", linewidth = 0.8) +
  geom_rug(sides = "b", alpha = 0.3, color = "gray50") +
  facet_wrap(~ variable, scales = "free", ncol = 5) +
  labs(title = "Partial Dependence Plots",
       subtitle = "Marginal Effect of Each Variable on Occurrence Probability",
       x = "Variable Value",
       y = "Occurrence Probability") +
  viz_theme_nature(base_size = 7, title_size = 9) +
  theme(axis.text = element_text(size = 5),
        axis.title = element_text(size = 6),
        strip.text = element_text(size = 6, face = "bold"))

ggsave("figures/17_local_sensitivity/partial_dependence.png",
       plot = p_pd, width = 5, height = 2.5, dpi = 1200)
ggsave("figures/17_local_sensitivity/partial_dependence.svg",
       plot = p_pd, width = 5, height = 2.5)

cat("      ✓ 偏依赖图已保存\n")

# 5c. 敏感度热图
cat("  (3) 敏感度热图...\n")

p_heatmap <- ggplot(sensitivity_summary, 
                    aes(x = lat_zone, y = variable, fill = mean_sensitivity)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.3f", mean_sensitivity)), 
            color = "white", size = 2.5, fontface = "bold") +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#D73027",
                       midpoint = 0, name = "Mean\nSensitivity") +
  labs(title = "Mean Sensitivity by Latitude Zone",
       x = "Latitude Zone", y = "Variable") +
  viz_theme_nature(base_size = 8, title_size = 9) +
  theme(panel.grid = element_blank())

ggsave("figures/17_local_sensitivity/sensitivity_heatmap.png",
       plot = p_heatmap, width = 3, height = 3, dpi = 1200)
ggsave("figures/17_local_sensitivity/sensitivity_heatmap.svg",
       plot = p_heatmap, width = 3, height = 3)

cat("      ✓ 热图已保存\n")

# 日志
sink("output/17_local_sensitivity/processing_log.txt")
cat("局部敏感度分析与偏依赖图日志\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
cat("分析变量数: ", length(top_vars), "\n", sep = "")
cat("样本数: ", nrow(analysis_data), "\n\n", sep = "")
cat("敏感度汇总:\n")
print(sensitivity_summary)
cat("\n最敏感变量 (按平均绝对敏感度):\n")
var_summary <- all_sensitivity %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(mean_abs_sens = mean(abs(sensitivity), na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(desc(mean_abs_sens))
print(var_summary)
cat("\n偏依赖分析:\n")
cat("计算了", length(unique(all_pd$variable)), "个变量的偏依赖\n")
cat("每个变量使用100个点进行插值\n")
sink()

cat("\n======================================\n")
cat("局部敏感度与偏依赖分析完成\n")
cat("======================================\n")
cat("输出文件:\n")
cat("  数据:\n")
cat("    - output/17_local_sensitivity/sensitivity_summary.csv\n")
cat("    - output/17_local_sensitivity/partial_dependence_data.csv\n")
cat("  图表:\n")
cat("    - figures/17_local_sensitivity/sensitivity_violin.png\n")
cat("    - figures/17_local_sensitivity/partial_dependence.png\n")
cat("    - figures/17_local_sensitivity/sensitivity_heatmap.png\n\n")

cat("✓ 脚本执行完成!\n\n")
