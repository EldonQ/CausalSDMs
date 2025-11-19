#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 13_study_area_and_points.R
# 功能说明: 绘制研究区域图、物种分布点图和背景点图
# 输入文件: output/01_data_preparation/species_occurrence_cleaned.csv
#          output/03_background_points/background_points.csv
# 输出文件: figures/13_study_area_and_points/study_area.png
#          figures/13_study_area_and_points/species_points.png
#          figures/13_study_area_and_points/background_points.png
#          figures/13_study_area_and_points/combined_map.png
# 作者: Nature级别科研项目
# 日期: 2025-10-20
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/SDM01")

# 加载必要的包
packages <- c("tidyverse", "sf", "ggplot2", "rnaturalearth", "viridis", "patchwork", "raster", "terra")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 统一可视化工具（Nature风格、Arial、1200dpi、viridis配色）
source("scripts/visualization/viz_utils.R")

dir.create("output/13_study_area_and_points", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/13_study_area_and_points", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("研究区域与分布点可视化\n")
cat("======================================\n\n")

# 1. 读取数据
cat("步骤 1/4: 读取数据...\n")

species_data <- read.csv("output/01_data_preparation/species_occurrence_cleaned.csv")
background_data <- read.csv("output/03_background_points/background_points.csv")

# 优先使用项目内中国边界（与其他图一致），失败则回退到 rnaturalearth
china <- viz_read_china("earthenvstreams_china/china_boundary.shp")
if (is.null(china)) {
  china <- ne_countries(country = "China", scale = "medium", returnclass = "sf")
}

# 读取河网强度图层（flow_acc 对数增强，归一化到[0,1]），用于优雅底图渲染（viridis透明叠加）
cat("  ✓ 读取并处理河网强度(flow_acc)用于底图...\n")
fa_spat <- terra::rast("earthenvstreams_china/flow_acc.tif")[[2]]
fa_spat[fa_spat <= 0] <- NA
fa_log <- log1p(fa_spat)
mx <- suppressWarnings(as.numeric(terra::global(fa_log, "max", na.rm = TRUE)[1,1]))
if (is.finite(mx) && mx > 0) fa_log <- fa_log / mx
fa_df <- as.data.frame(fa_log, xy = TRUE, na.rm = TRUE)
colnames(fa_df) <- c("lon","lat","val")
if(nrow(fa_df) > 120000) {
  set.seed(1)
  fa_df <- fa_df[sample(seq_len(nrow(fa_df)), 120000), ]
}

cat("  ✓ 物种分布点: ", nrow(species_data), "\n", sep = "")
cat("  ✓ 背景点: ", nrow(background_data), "\n", sep = "")

# 2. 绘制研究区域图
cat("\n步骤 2/4: 绘制研究区域图...\n")

p_area <- ggplot() +
  # 中国底图内部填充为RGB(245,245,245)，边界线为黑色
  geom_sf(data = china, fill = "#DCDCDC", color = "black", linewidth = 0.4) +
  labs(title = "Study Area: China",
       x = "Longitude (\u00B0E)", y = "Latitude (\u00B0N)") +
  viz_theme_nature(base_size = 9) +
  # 去掉坐标轴与背景网格线（仅保留图形主体与标题）
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 画布与绘图区背景透明，确保导出为透明背景
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

ggsave("figures/13_study_area_and_points/study_area.png",
       plot = p_area, width = 4, height = 3.2, dpi = 1200, bg = "transparent")
ggsave("figures/13_study_area_and_points/study_area.svg",
       plot = p_area, width = 4, height = 3.2, bg = "transparent")

cat("  ✓ 研究区域图\n")

# 3. 绘制物种分布点图
cat("\n步骤 3/4: 绘制物种分布点图...\n")

p_species <- ggplot() +
  # 底图先填充中国区域为RGB(245,245,245)，仅填充不绘制边界
  geom_sf(data = china, fill = "#DCDCDC", color = NA) +
  geom_raster(data = fa_df, aes(x = lon, y = lat, fill = val), alpha = 0.35) +
  scale_fill_viridis(name = NULL, option = "C", guide = "none") +
  # 边界线单独置于最上层，保证边界清晰
  geom_sf(data = china, fill = NA, color = "black", linewidth = 0.3) +
  # 物种发生点单独图：使用小黑点，适度缩小点大小
  geom_point(data = species_data, aes(x = lon, y = lat),
             color = "black", size = 0.4, alpha = 0.8) +
  labs(title = "Species Occurrence Points",
       x = "Longitude (\u00B0E)", y = "Latitude (\u00B0N)") +
  viz_theme_nature(base_size = 9) +
  # 去掉坐标轴与背景网格线
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 画布与绘图区背景透明
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

ggsave("figures/13_study_area_and_points/species_points.png",
       plot = p_species, width = 4.2, height = 3.2, dpi = 1200, bg = "transparent")
ggsave("figures/13_study_area_and_points/species_points.svg",
       plot = p_species, width = 4.2, height = 3.2, bg = "transparent")

cat("  ✓ 物种分布点图\n")

# 绘制背景点图
p_background <- ggplot() +
  # 底图先填充中国区域为RGB(245,245,245)
  geom_sf(data = china, fill = "#DCDCDC", color = NA) +
  geom_raster(data = fa_df, aes(x = lon, y = lat, fill = val), alpha = 0.35) +
  scale_fill_viridis(name = NULL, option = "C", guide = "none") +
  geom_sf(data = china, fill = NA, color = "black", linewidth = 0.3) +
  # 背景点适度缩小点大小
  geom_point(data = background_data, aes(x = lon, y = lat),
             color = viridis::viridis(1, option = "D"), size = 0.25, alpha = 0.7) +
  labs(title = "Background Points",
       x = "Longitude (\u00B0E)", y = "Latitude (\u00B0N)") +
  viz_theme_nature(base_size = 9) +
  # 去掉坐标轴与背景网格线
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 画布与绘图区背景透明
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

ggsave("figures/13_study_area_and_points/background_points.png",
       plot = p_background, width = 4.2, height = 3.2, dpi = 1200, bg = "transparent")
ggsave("figures/13_study_area_and_points/background_points.svg",
       plot = p_background, width = 4.2, height = 3.2, bg = "transparent")

cat("  ✓ 背景点图\n")

# 4. 绘制组合图
cat("\n步骤 4/4: 绘制组合图...\n")

p_combined <- ggplot() +
  # 底图先填充中国区域为RGB(245,245,245)
  geom_sf(data = china, fill = "#DCDCDC", color = NA) +
  geom_raster(data = fa_df, aes(x = lon, y = lat, fill = val), alpha = 0.35) +
  scale_fill_viridis(name = NULL, option = "C", guide = "none") +
  geom_sf(data = china, fill = NA, color = "black", linewidth = 0.3) +
  # 组合图：缩小背景点与发生点
  geom_point(data = background_data, aes(x = lon, y = lat),
             color = viridis::viridis(1, option = "D"), size = 0.25, alpha = 0.55) +
  geom_point(data = species_data, aes(x = lon, y = lat),
             color = "#E41A1C", size = 0.4, alpha = 0.85) +
  labs(title = "Species Presence and Background Points",
       x = "Longitude (\u00B0E)", y = "Latitude (\u00B0N)") +
  viz_theme_nature(base_size = 9) +
  # 去掉坐标轴与背景网格线
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # 画布与绘图区背景透明
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

ggsave("figures/13_study_area_and_points/combined_map.png",
       plot = p_combined, width = 4.5, height = 3.4, dpi = 1200, bg = "transparent")
ggsave("figures/13_study_area_and_points/combined_map.svg",
       plot = p_combined, width = 4.5, height = 3.4, bg = "transparent")

cat("  ✓ 组合图\n")

# 保存统计信息
summary_data <- data.frame(
  type = c("Species Points", "Background Points", "Total"),
  count = c(nrow(species_data), nrow(background_data), 
            nrow(species_data) + nrow(background_data)),
  lon_range = c(paste(range(species_data$lon), collapse = " to "),
                paste(range(background_data$lon), collapse = " to "),
                paste(range(c(species_data$lon, background_data$lon)), collapse = " to ")),
  lat_range = c(paste(range(species_data$lat), collapse = " to "),
                paste(range(background_data$lat), collapse = " to "),
                paste(range(c(species_data$lat, background_data$lat)), collapse = " to "))
)

write.csv(summary_data, "output/13_study_area_and_points/points_summary.csv", row.names = FALSE)

# 日志
sink("output/13_study_area_and_points/processing_log.txt")
cat("研究区域与分布点可视化日志\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
print(summary_data)
sink()

cat("\n======================================\n")
cat("可视化完成\n")
cat("======================================\n\n")

cat(sprintf("物种分布点: %d\n", nrow(species_data)))
cat(sprintf("背景点: %d\n", nrow(background_data)))
cat(sprintf("背景点/物种点比例: %.2f:1\n", nrow(background_data)/nrow(species_data)))

cat("\n✓ 脚本执行完成!\n\n")

