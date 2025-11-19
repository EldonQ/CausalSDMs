#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 04_collinearity_analysis.R
# 功能说明: 环境变量共线性分析，筛选非共线性变量子集
# 方法: Pearson相关系数 + VIF方差膨胀因子
# 输入文件: output/03_background_points/combined_presence_absence.csv
# 输出文件: output/04_collinearity/collinearity_removed.csv
#          output/04_collinearity/selected_variables.csv
# 作者: Nature级别科研项目
# 日期: 2025-10-20
# ==============================================================================

# 清空环境
rm(list = ls())
gc()

# 设置工作目录
setwd("E:/SDM01")

# 加载必要的包
packages <- c("tidyverse", "corrplot", "usdm")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 创建输出目录
if(!dir.exists("output/04_collinearity")) {
  dir.create("output/04_collinearity", recursive = TRUE)
}
if(!dir.exists("figures/04_collinearity")) {
  dir.create("figures/04_collinearity", recursive = TRUE)
}

cat("======================================\n")
cat("环境变量共线性分析\n")
cat("======================================\n\n")

# ------------------------------------------------------------------------------
# 1. 读取数据
# ------------------------------------------------------------------------------
cat("步骤 1/7: 读取组合数据...\n")

model_data <- read.csv("output/03_background_points/combined_presence_absence.csv")
cat("  - 总记录数: ", nrow(model_data), "\n", sep = "")
cat("  - 出现点: ", sum(model_data$presence == 1), "\n", sep = "")
cat("  - 背景点: ", sum(model_data$presence == 0), "\n", sep = "")

# 提取环境变量
env_vars <- names(model_data)[6:(ncol(model_data)-1)]
env_data <- model_data[, env_vars]

cat("  - 环境变量数: ", length(env_vars), "\n", sep = "")

# ------------------------------------------------------------------------------
# 严格变量模式：若检测到 scripts/variables_selected_47.csv，则直接保留该白名单
# ------------------------------------------------------------------------------
override_path <- "scripts/variables_selected_47.csv"
if (file.exists(override_path)) {
  cat("\n检测到变量白名单: ", override_path, "\n", sep = "")
  cat("启用严格变量模式：跳过VIF与相关性剔除，保留白名单全部变量。\n")
  whitelist_df <- read.csv(override_path, stringsAsFactors = FALSE)
  whitelist <- whitelist_df$variable
  # 按白名单顺序保留（而非 env_vars 顺序），确保四组变量的既定排列
  selected_vars_cor <- whitelist[whitelist %in% env_vars]
  cat("  - 白名单变量数: ", length(selected_vars_cor), "\n", sep = "")

  # 相关性矩阵（用于保存图表与表格，满足期刊留痕）
  cor_matrix <- cor(model_data[, selected_vars_cor], use = "complete.obs", method = "pearson")
  if(!dir.exists("output/04_collinearity")) dir.create("output/04_collinearity", recursive = TRUE)
  if(!dir.exists("figures/04_collinearity")) dir.create("figures/04_collinearity", recursive = TRUE)
  write.csv(cor_matrix, "output/04_collinearity/correlation_matrix.csv", row.names = TRUE)

  # 保存选中变量列表（保持与下游兼容，仅含variable列；index用于可读性）
  selected_var_list <- data.frame(index = 1:length(selected_vars_cor), variable = selected_vars_cor, stringsAsFactors = FALSE)
  write.csv(selected_var_list, "output/04_collinearity/selected_variables.csv", row.names = FALSE)
  # 额外保存变量-分组映射（便于后续健壮性与分组可视化）
  sel_map <- whitelist_df[whitelist_df$variable %in% selected_vars_cor, c("variable","category","file","band")]
  # 按 selected_vars_cor 的顺序排列
  sel_map <- sel_map[match(selected_vars_cor, sel_map$variable), ]
  write.csv(sel_map, "output/04_collinearity/selected_variables_with_group.csv", row.names = FALSE)

  # 创建最终建模数据集（仅保留白名单变量）
  model_data_final <- model_data[, c("id", "species", "lon", "lat", "source", selected_vars_cor, "presence")]
  write.csv(model_data_final, "output/04_collinearity/collinearity_removed.csv", row.names = FALSE)

  # 可视化：原始（即白名单）相关性热图
  png("figures/04_collinearity/correlation_heatmap_original.png",
      width = 4000, height = 4000, res = 1200, family = "Arial")
  par(mar = c(1, 1, 2, 1))
  corrplot(cor_matrix,
           method = "color",
           type = "upper",
           tl.col = "black",
           tl.srt = 45,
           tl.cex = 0.3,
           col = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200),
           title = "Correlation Matrix (Selected 47 Variables)",
           mar = c(0, 0, 1, 0))
  dev.off()
  cat("  ✓ 已保存: figures/04_collinearity/correlation_heatmap_original.png\n")

  # 结束并摘要
  cat("\n严格变量模式完成：已将47变量直接用于后续建模。\n")
  cat("  - 最终变量数: ", length(selected_vars_cor), "\n", sep = "")
  STRICT_MODE_DONE <- TRUE
} else {

# ------------------------------------------------------------------------------
# 2. 方差筛选
# ------------------------------------------------------------------------------
cat("\n步骤 2/7: 方差筛选（移除常数和低方差变量）...\n")

var_sd <- apply(env_data, 2, sd, na.rm = TRUE)

# 识别常数或低方差变量
zero_var <- names(var_sd)[var_sd == 0 | is.na(var_sd)]
low_var <- names(var_sd)[var_sd > 0 & var_sd < 0.01]

cat("  - 常数变量数量: ", length(zero_var), "\n", sep = "")
cat("  - 低方差变量数量 (SD < 0.01): ", length(low_var), "\n", sep = "")

# 移除常数和低方差变量
vars_to_remove <- c(zero_var, low_var)
if(length(vars_to_remove) > 0) {
  cat("  - 移除 ", length(vars_to_remove), " 个低方差变量\n", sep = "")
  env_vars_filtered <- env_vars[!env_vars %in% vars_to_remove]
  env_data_filtered <- env_data[, env_vars_filtered]
} else {
  env_vars_filtered <- env_vars
  env_data_filtered <- env_data
  cat("  - 所有变量均有足够方差\n")
}

cat("  - 方差筛选后保留: ", length(env_vars_filtered), " 个变量\n", sep = "")

# ------------------------------------------------------------------------------
# 3. 计算相关系数矩阵
# ------------------------------------------------------------------------------
cat("\n步骤 3/7: 计算Pearson相关系数矩阵...\n")

cor_matrix <- cor(env_data_filtered, use = "complete.obs", method = "pearson")

# 统计高度相关的变量对
high_cor_threshold <- 0.7
high_cor <- which(abs(cor_matrix) > high_cor_threshold & 
                  abs(cor_matrix) < 1, arr.ind = TRUE)

if(length(high_cor) > 0) {
  high_cor_pairs <- data.frame(
    var1 = rownames(cor_matrix)[high_cor[, 1]],
    var2 = colnames(cor_matrix)[high_cor[, 2]],
    correlation = cor_matrix[high_cor],
    stringsAsFactors = FALSE
  )
  
  # 移除重复对
  high_cor_pairs <- high_cor_pairs[high_cor_pairs$var1 < high_cor_pairs$var2, ]
  high_cor_pairs <- high_cor_pairs[order(abs(high_cor_pairs$correlation), 
                                         decreasing = TRUE), ]
  
  cat("  - 高相关变量对数量 (|r| > 0.7): ", nrow(high_cor_pairs), "\n", sep = "")
} else {
  cat("  - 无高相关变量对 (|r| > 0.7)\n")
}

write.csv(cor_matrix, "output/04_collinearity/correlation_matrix.csv", 
          row.names = TRUE)

# ------------------------------------------------------------------------------
# 4. 使用VIF进行共线性诊断
# ------------------------------------------------------------------------------
cat("\n步骤 4/7: 计算方差膨胀因子(VIF)...\n")
cat("  注意: 这可能需要较长时间...\n")

vif_results <- tryCatch({
  vifstep(env_data_filtered, th = 10)
}, error = function(e) {
  cat("  - VIF计算出错，将使用相关系数法\n")
  return(NULL)
})

if(!is.null(vif_results)) {
  selected_vars_vif <- vif_results@results$Variables
  excluded_vars_vif <- vif_results@excluded
  
  cat("  - VIF分析完成\n")
  cat("  - 保留的变量数: ", length(selected_vars_vif), "\n", sep = "")
  cat("  - 移除的变量数: ", length(excluded_vars_vif), "\n", sep = "")
  
  selected_vars <- selected_vars_vif
} else {
  cat("  - 使用备选方法: 基于相关系数的变量选择\n")
  selected_vars <- env_vars_filtered
}

# ------------------------------------------------------------------------------
# 5. 基于相关系数的额外筛选
# ------------------------------------------------------------------------------
cat("\n步骤 5/7: 基于相关系数的额外筛选 (|r| < 0.8)...\n")

selected_vars_cor <- selected_vars
cor_threshold <- 0.8
removed_count <- 0

repeat {
  if(length(selected_vars_cor) < 2) break
  
  cor_current <- cor(env_data_filtered[, selected_vars_cor], use = "complete.obs")
  
  # 找出最高相关对
  cor_upper <- cor_current
  cor_upper[lower.tri(cor_upper, diag = TRUE)] <- 0
  max_cor <- max(abs(cor_upper), na.rm = TRUE)
  
  if(is.na(max_cor) || is.infinite(max_cor) || max_cor < cor_threshold) break
  
  # 找到相关性最高的变量对
  max_pos <- which(abs(cor_upper) == max_cor, arr.ind = TRUE)[1, ]
  var1 <- selected_vars_cor[max_pos[1]]
  var2 <- selected_vars_cor[max_pos[2]]
  
  # 计算每个变量与其他变量的平均相关性
  cor_row1 <- abs(cor_current[max_pos[1], ])
  cor_row1[max_pos[1]] <- NA
  mean_cor1 <- mean(cor_row1, na.rm = TRUE)
  
  cor_row2 <- abs(cor_current[max_pos[2], ])
  cor_row2[max_pos[2]] <- NA
  mean_cor2 <- mean(cor_row2, na.rm = TRUE)
  
  # 移除平均相关性更高的变量
  if(is.na(mean_cor1) || is.na(mean_cor2)) {
    selected_vars_cor <- selected_vars_cor[selected_vars_cor != var1]
  } else if(mean_cor1 > mean_cor2) {
    selected_vars_cor <- selected_vars_cor[selected_vars_cor != var1]
  } else {
    selected_vars_cor <- selected_vars_cor[selected_vars_cor != var2]
  }
  
  removed_count <- removed_count + 1
  
  # 安全机制
  if(removed_count > 100) {
    cat("  - 警告: 已移除100个变量，停止进一步筛选\n")
    break
  }
}

cat("  - 共移除 ", removed_count, " 个高相关变量\n", sep = "")
cat("  - 最终保留的变量数: ", length(selected_vars_cor), "\n", sep = "")

# ------------------------------------------------------------------------------
# 6. 保存筛选结果
# ------------------------------------------------------------------------------
cat("\n步骤 6/7: 保存筛选结果...\n")

# 保存选中的变量列表
selected_var_list <- data.frame(
  index = 1:length(selected_vars_cor),
  variable = selected_vars_cor,
  stringsAsFactors = FALSE
)

write.csv(selected_var_list,
          "output/04_collinearity/selected_variables.csv",
          row.names = FALSE)
cat("  ✓ 已保存: output/04_collinearity/selected_variables.csv\n")

# 创建最终建模数据集
model_data_final <- model_data[, c("id", "species", "lon", "lat", 
                                   "source", selected_vars_cor, "presence")]

write.csv(model_data_final,
          "output/04_collinearity/collinearity_removed.csv",
          row.names = FALSE)
cat("  ✓ 已保存: output/04_collinearity/collinearity_removed.csv\n")

# ------------------------------------------------------------------------------
# 7. 生成可视化
# ------------------------------------------------------------------------------
cat("\n步骤 7/7: 生成可视化...\n")

# 绘制原始相关性热图
png("figures/04_collinearity/correlation_heatmap_original.png",
    width = 4000, height = 4000, res = 1200, family = "Arial")

par(mar = c(1, 1, 2, 1))
corrplot(cor_matrix,
         method = "color",
         type = "upper",
         tl.col = "black",
         tl.srt = 45,
         tl.cex = 0.3,
         col = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200),
         title = "Correlation Matrix (Original Variables)",
         mar = c(0, 0, 1, 0))

dev.off()
cat("  ✓ 已保存: figures/04_collinearity/correlation_heatmap_original.png\n")

# 绘制筛选后的相关性热图
cor_final <- cor(model_data_final[, selected_vars_cor], use = "complete.obs")

png("figures/04_collinearity/correlation_heatmap_final.png",
    width = 3000, height = 3000, res = 1200, family = "Arial")

if(length(selected_vars_cor) <= 20) {
  corrplot(cor_final,
           method = "color",
           type = "upper",
           tl.col = "black",
           tl.srt = 45,
           tl.cex = 0.6,
           col = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200),
           addCoef.col = "black",
           number.cex = 0.5,
           title = "Correlation Matrix (Selected Variables)",
           mar = c(0, 0, 1, 0))
} else {
  corrplot(cor_final,
           method = "color",
           type = "upper",
           tl.col = "black",
           tl.srt = 45,
           tl.cex = 0.4,
           col = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200),
           title = "Correlation Matrix (Selected Variables)",
           mar = c(0, 0, 1, 0))
}

dev.off()
cat("  ✓ 已保存: figures/04_collinearity/correlation_heatmap_final.png\n")

# 变量选择对比图
var_counts <- data.frame(
  stage = c("Original", "After Collinearity Removal"),
  n_vars = c(length(env_vars), length(selected_vars_cor))
)

png("figures/04_collinearity/variable_selection_summary.png",
    width = 2400, height = 1800, res = 1200, family = "Arial")

par(mar = c(4, 4, 2, 1))
barplot(var_counts$n_vars,
        names.arg = var_counts$stage,
        col = c("#377EB8", "#4DAF4A"),
        ylab = "Number of Variables",
        main = "Variable Selection Summary",
        ylim = c(0, max(var_counts$n_vars) * 1.2),
        cex.names = 0.8,
        cex.axis = 0.8,
        cex.lab = 0.9)
text(x = 1:2, y = var_counts$n_vars + max(var_counts$n_vars) * 0.05,
     labels = var_counts$n_vars, cex = 0.8)

dev.off()
cat("  ✓ 已保存: figures/04_collinearity/variable_selection_summary.png\n")

# ------------------------------------------------------------------------------
# 摘要
# ------------------------------------------------------------------------------
cat("\n======================================\n")
cat("共线性分析完成\n")
cat("======================================\n")
cat("原始变量数: ", length(env_vars), "\n", sep = "")
cat("筛选后变量数: ", length(selected_vars_cor), "\n", sep = "")
cat("移除变量数: ", length(env_vars) - length(selected_vars_cor), 
    " (", round((1 - length(selected_vars_cor) / length(env_vars)) * 100, 1), "%)\n", sep = "")

# 相关性统计
cor_abs <- abs(cor_final[upper.tri(cor_final)])
cor_abs <- cor_abs[!is.na(cor_abs)]

cat("\n筛选后相关性统计:\n")
cat("  平均相关系数: ", round(mean(cor_abs), 3), "\n", sep = "")
cat("  最大相关系数: ", round(max(cor_abs), 3), "\n", sep = "")
cat("  |r| > 0.7 的变量对数: ", sum(cor_abs > 0.7), "\n", sep = "")

# 保存处理日志
sink("output/04_collinearity/processing_log.txt")
cat("共线性分析日志\n")
cat("处理时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("=== 数据统计 ===\n")
cat("总记录数: ", nrow(model_data_final), "\n", sep = "")
cat("出现点: ", sum(model_data_final$presence == 1), "\n", sep = "")
cat("背景点: ", sum(model_data_final$presence == 0), "\n\n", sep = "")

cat("=== 变量筛选 ===\n")
cat("原始变量数: ", length(env_vars), "\n", sep = "")
cat("筛选后变量数: ", length(selected_vars_cor), "\n", sep = "")
cat("移除变量数: ", length(env_vars) - length(selected_vars_cor), "\n\n", sep = "")

cat("=== 相关性统计 ===\n")
cat("平均相关系数: ", round(mean(cor_abs, na.rm = TRUE), 3), "\n", sep = "")
cat("最大相关系数: ", round(max(cor_abs, na.rm = TRUE), 3), "\n", sep = "")
cat("|r| > 0.7 的变量对数: ", sum(cor_abs > 0.7, na.rm = TRUE), "\n\n", sep = "")

cat("=== 筛选后的变量列表 ===\n")
print(selected_var_list)

sink()

cat("  ✓ 已保存: output/04_collinearity/processing_log.txt\n")
cat("\n✓ 数据已准备就绪，可用于建模\n\n")

cat("脚本执行完成！\n")
}
