#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 04_uncertainty_map_clean.R
# 功能说明: 重新绘制12_uncertainty_map.R的所有结果（高质量版：2400dpi、透明背景、专业科研配色）
#           分析4个模型预测的不确定性（标准差）与一致性
# 参考脚本: scripts/12_uncertainty_map.R（完整功能复刻，优化绘图质量与配色）
# 输入文件: output/11_prediction_maps/rasters/pred_*_river.tif
# 输出文件: Drawing/output/uncertainty/uncertainty_map.png (2400 dpi)
#          Drawing/output/uncertainty/model_agreement.png (2400 dpi)
#          Drawing/output/uncertainty/*.tif (GeoTIFF)
#          Drawing/output/uncertainty/uncertainty_summary.csv
# 作者: Nature级别科研项目
# 日期: 2025-11-09
# ==============================================================================

# -----------------------------
# 0. 初始化环境
# -----------------------------
rm(list = ls())
gc()
setwd("E:/SDM01")

# -----------------------------
# 1. 加载必要的包（缺则自动安装）
# -----------------------------
packages <- c(
  "tidyverse",     # 数据整理
  "sf",            # 矢量数据
  "ggplot2",       # 绘图
  "raster",        # 栅格处理（兼容性）
  "terra",         # 现代栅格处理
  "viridis",       # 配色
  "rnaturalearth", # 边界
  "scico",         # 科学配色包（比viridis更专业）
  "patchwork"      # 图形组合
)
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 统一可视化工具（Nature风格、Arial、高dpi）
source("scripts/visualization/viz_utils.R")

# -----------------------------
# 2. 目录准备
# -----------------------------
dir.create("Drawing/output/uncertainty", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("模型不确定性分析（高质量重绘版）\n")
cat("======================================\n\n")

# -----------------------------
# 3. 读取河网预测栅格
# -----------------------------
cat("步骤 1/4: 读取河网预测栅格...\n")

ras_paths <- list(
  Maxnet = "output/11_prediction_maps/rasters/pred_maxnet_river.tif",
  NN = "output/11_prediction_maps/rasters/pred_nn_river.tif",
  RF = "output/11_prediction_maps/rasters/pred_rf_river.tif",
  GAM = "output/11_prediction_maps/rasters/pred_gam_river.tif"
)

missing <- names(ras_paths)[!file.exists(unlist(ras_paths))]
if (length(missing) > 0) {
  stop(paste0("缺少预测栅格: ", paste(missing, collapse = ", "),
              "。请先运行 11_current_prediction_maps.R 生成河网预测。"))
}

stk <- raster::stack(unlist(ras_paths))
names(stk) <- names(ras_paths)

cat("  ✓ 输入层数: ", raster::nlayers(stk), "\n", sep = "")

# -----------------------------
# 4. 计算不确定性指标（逐像元）
# -----------------------------
cat("\n步骤 2/4: 计算不确定性指标（逐像元）...\n")

mean_r <- raster::calc(stk, fun = function(x) { 
  if (all(is.na(x))) NA else mean(x, na.rm = TRUE) 
})
sd_r <- raster::calc(stk, fun = function(x) { 
  if (all(is.na(x))) NA else sd(x, na.rm = TRUE) 
})
min_r <- raster::calc(stk, fun = function(x) { 
  if (all(is.na(x))) NA else min(x, na.rm = TRUE) 
})
max_r <- raster::calc(stk, fun = function(x) { 
  if (all(is.na(x))) NA else max(x, na.rm = TRUE) 
})
range_r <- max_r - min_r
agreement_r <- 1 - range_r  # 一致性：范围越小，一致性越高

# 保存GeoTIFF
raster::writeRaster(mean_r, "Drawing/output/uncertainty/mean_prediction_river.tif", overwrite = TRUE)
raster::writeRaster(sd_r, "Drawing/output/uncertainty/sd_prediction_river.tif", overwrite = TRUE)
raster::writeRaster(agreement_r, "Drawing/output/uncertainty/agreement_river.tif", overwrite = TRUE)
raster::writeRaster(range_r, "Drawing/output/uncertainty/range_river.tif", overwrite = TRUE)

cat("  ✓ 均值栅格: Drawing/output/uncertainty/mean_prediction_river.tif\n")
cat("  ✓ 标准差栅格: Drawing/output/uncertainty/sd_prediction_river.tif\n")
cat("  ✓ 一致性栅格: Drawing/output/uncertainty/agreement_river.tif\n")

# -----------------------------
# 5. 绘制高质量热图（2400 dpi，透明背景，专业配色）
# -----------------------------
cat("\n步骤 3/4: 绘制高质量热图...\n")

# 读取中国边界
china <- viz_read_china("earthenvstreams_china/china_boundary.shp")
if (is.null(china)) {
  china <- rnaturalearth::ne_countries(country = "China", scale = "medium", returnclass = "sf")
}

# 转为数据框用于 ggplot2 绘图
sd_df <- as.data.frame(raster::rasterToPoints(sd_r))
colnames(sd_df) <- c("lon", "lat", "value")
agreement_df <- as.data.frame(raster::rasterToPoints(agreement_r))
colnames(agreement_df) <- c("lon", "lat", "value")

# 提升对比度：使用 1%-99% 分位数裁剪极端值
sd_limits <- quantile(sd_df$value, probs = c(0.01, 0.99), na.rm = TRUE)
agreement_limits <- quantile(agreement_df$value, probs = c(0.01, 0.99), na.rm = TRUE)

cat("  - 不确定性（SD）范围: [", sprintf("%.4f", sd_limits[1]), ", ", sprintf("%.4f", sd_limits[2]), "]\n", sep = "")
cat("  - 一致性范围: [", sprintf("%.4f", agreement_limits[1]), ", ", sprintf("%.4f", agreement_limits[2]), "]\n", sep = "")

# 绘图1: 不确定性（标准差）热图
# 专业科研配色：scico "vik"（发散型，适合不确定性）或 "lajolla"（连续型，黄-红）
cat("\n  -> 绘制不确定性热图（专业配色：scico 'lajolla'）...\n")

p_uncertainty <- ggplot() +
  geom_raster(data = sd_df, aes(x = lon, y = lat, fill = value)) +
  # 专业科研热图配色：lajolla（黄-橙-红，适合不确定性/热度）
  scale_fill_scico(
    name = "SD",
    palette = "lajolla",
    limits = sd_limits,
    oob = scales::squish,
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 12,
      barheight = 0.5
    )
  ) +
  geom_sf(data = china, fill = NA, color = "black", linewidth = 0.3) +
  labs(
    title = "Prediction Uncertainty (Standard Deviation)",
    x = "Longitude (°E)",
    y = "Latitude (°N)"
  ) +
  coord_sf() +
  theme_minimal(base_size = 10, base_family = "Arial") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", color = NA),
    axis.title = element_text(size = 10, family = "Arial"),
    plot.title = element_text(size = 12, family = "Arial", face = "bold", hjust = 0.5)
  )

# 保存 PNG (2400 dpi) + SVG
ggsave(
  "Drawing/output/uncertainty/uncertainty_map.png",
  plot = p_uncertainty, width = 8, height = 6, dpi = 2400, bg = "transparent"
)
ggsave(
  "Drawing/output/uncertainty/uncertainty_map.svg",
  plot = p_uncertainty, width = 8, height = 6, bg = "transparent"
)

cat("    ✓ 不确定性热图: Drawing/output/uncertainty/uncertainty_map.png\n")

# 绘图2: 模型一致性热图
# 专业科研配色：scico "bamako"（绿-黄-棕，适合一致性）或 viridis
cat("\n  -> 绘制模型一致性热图（专业配色：scico 'bamako'）...\n")

p_agreement <- ggplot() +
  geom_raster(data = agreement_df, aes(x = lon, y = lat, fill = value)) +
  # 专业科研热图配色：bamako（绿-黄-棕，适合一致性/质量）
  scale_fill_scico(
    name = "Agreement",
    palette = "bamako",
    limits = agreement_limits,
    oob = scales::squish,
    direction = 1,
    guide = guide_colorbar(
      title.position = "top",
      title.hjust = 0.5,
      barwidth = 12,
      barheight = 0.5
    )
  ) +
  geom_sf(data = china, fill = NA, color = "black", linewidth = 0.3) +
  labs(
    title = "Model Agreement (1 - Range)",
    x = "Longitude (°E)",
    y = "Latitude (°N)"
  ) +
  coord_sf() +
  theme_minimal(base_size = 10, base_family = "Arial") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    legend.position = "bottom",
    legend.background = element_rect(fill = "transparent", color = NA),
    axis.title = element_text(size = 10, family = "Arial"),
    plot.title = element_text(size = 12, family = "Arial", face = "bold", hjust = 0.5)
  )

# 保存 PNG (2400 dpi) + SVG
ggsave(
  "Drawing/output/uncertainty/model_agreement.png",
  plot = p_agreement, width = 8, height = 6, dpi = 2400, bg = "transparent"
)
ggsave(
  "Drawing/output/uncertainty/model_agreement.svg",
  plot = p_agreement, width = 8, height = 6, bg = "transparent"
)

cat("    ✓ 一致性热图: Drawing/output/uncertainty/model_agreement.png\n")

# -----------------------------
# 6. 输出统计摘要
# -----------------------------
cat("\n步骤 4/4: 输出统计摘要...\n")

vals_sd <- raster::getValues(sd_r)
vals_sd <- vals_sd[!is.na(vals_sd)]
vals_ag <- raster::getValues(agreement_r)
vals_ag <- vals_ag[!is.na(vals_ag)]
vals_mean <- raster::getValues(mean_r)
vals_mean <- vals_mean[!is.na(vals_mean)]
vals_range <- raster::getValues(range_r)
vals_range <- vals_range[!is.na(vals_range)]

uncertainty_summary <- data.frame(
  metric = c("mean_prediction", "sd_prediction", "agreement", "range"),
  mean = c(mean(vals_mean), mean(vals_sd), mean(vals_ag), mean(vals_range)),
  sd = c(sd(vals_mean), sd(vals_sd), sd(vals_ag), sd(vals_range)),
  min = c(min(vals_mean), min(vals_sd), min(vals_ag), min(vals_range)),
  max = c(max(vals_mean), max(vals_sd), max(vals_ag), max(vals_range)),
  p10 = c(
    as.numeric(quantile(vals_mean, 0.1)),
    as.numeric(quantile(vals_sd, 0.1)),
    as.numeric(quantile(vals_ag, 0.1)),
    as.numeric(quantile(vals_range, 0.1))
  ),
  p50 = c(
    as.numeric(quantile(vals_mean, 0.5)),
    as.numeric(quantile(vals_sd, 0.5)),
    as.numeric(quantile(vals_ag, 0.5)),
    as.numeric(quantile(vals_range, 0.5))
  ),
  p90 = c(
    as.numeric(quantile(vals_mean, 0.9)),
    as.numeric(quantile(vals_sd, 0.9)),
    as.numeric(quantile(vals_ag, 0.9)),
    as.numeric(quantile(vals_range, 0.9))
  )
)

write.csv(uncertainty_summary, "Drawing/output/uncertainty/uncertainty_summary.csv", row.names = FALSE)

cat("  ✓ 统计摘要: Drawing/output/uncertainty/uncertainty_summary.csv\n")

# 日志
sink("Drawing/output/uncertainty/processing_log.txt")
cat("不确定性分析日志（高质量重绘版）\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
cat("绘图设置:\n")
cat("  - 输出分辨率: 2400 dpi\n")
cat("  - 画布背景: 透明\n")
cat("  - 不确定性配色: scico 'lajolla' (黄-橙-红，科研热图)\n")
cat("  - 一致性配色: scico 'bamako' (绿-黄-棕，科研热图)\n")
cat("  - 对比度增强: 1%-99% 分位数裁剪\n\n")
cat("输入预测栅格: 4 (Maxnet, NN, RF, GAM)\n\n")
cat("统计摘要（河网像元）:\n")
print(uncertainty_summary)
sink()

cat("\n======================================\n")
cat("不确定性分析完成（高质量版）\n")
cat("======================================\n\n")

cat(sprintf("✓ 不确定性热图: Drawing/output/uncertainty/uncertainty_map.png\n"))
cat(sprintf("✓ 一致性热图: Drawing/output/uncertainty/model_agreement.png\n"))
cat(sprintf("✓ GeoTIFF栅格: Drawing/output/uncertainty/*.tif\n"))
cat(sprintf("✓ 统计摘要: Drawing/output/uncertainty/uncertainty_summary.csv\n"))
cat(sprintf("✓ 分辨率: 2400 dpi（超高质量）\n"))
cat(sprintf("✓ 配色: 专业科研热图（scico palette）\n"))

cat("\n✓ 脚本执行完成!\n\n")

