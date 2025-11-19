#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 02_grid_colored_regions.R
# 功能说明: 在中国区域创建精细方形网格（5km），随机均匀选择部分网格上色（黑色、蓝色）
#           用于Nature级别期刊配图（英文标注、Arial、2400dpi、透明画布）
# 输入文件: earthenvstreams_china/china_boundary.shp
# 输出文件: Drawing/output/grid_colored_regions.png
#          Drawing/output/grid_colored_regions.svg
#          Drawing/output/grid_statistics.csv
#          Drawing/output/grid_log.txt
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
  "terra",         # 栅格处理
  "rnaturalearth"  # 边界回退
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
dir.create("Drawing/output", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("中国区域网格化与随机上色\n")
cat("======================================\n\n")

# -----------------------------
# 3. 读取中国边界
# -----------------------------
cat("步骤 1/4: 读取中国边界...\n")

viz_china <- viz_read_china("earthenvstreams_china/china_boundary.shp")
if (is.null(viz_china)) {
  viz_china <- rnaturalearth::ne_countries(country = "China", scale = "medium", returnclass = "sf")
}

# 转为投影坐标系（等积投影，便于精确创建5km网格）
# 使用 Albers Equal Area Conic 投影（中国标准）
china_proj <- sf::st_transform(viz_china, 
                               crs = "+proj=aea +lat_1=25 +lat_2=47 +lat_0=35 +lon_0=105 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs")

cat("  ✓ 中国边界已加载并投影至等积坐标系（单位：米）\n")

# -----------------------------
# 4. 创建精细方形网格（5km）
# -----------------------------
cat("\n步骤 2/4: 创建 5km 方形网格...\n")

# 获取中国边界的范围
bbox <- sf::st_bbox(china_proj)

# 创建 5km × 5km 网格（5000米）
grid_size <- 5000  # 米
grid <- sf::st_make_grid(
  china_proj,
  cellsize = c(grid_size, grid_size),
  what = "polygons",
  square = TRUE
)

# 转为 sf 对象并添加 ID
grid_sf <- sf::st_sf(grid_id = seq_along(grid), geometry = grid)

# 仅保留与中国边界相交的网格
grid_china <- sf::st_intersection(grid_sf, china_proj)

cat("  ✓ 网格生成完成，共 ", nrow(grid_china), " 个有效网格单元（5km × 5km）\n", sep = "")

# -----------------------------
# 5. 随机均匀选择网格并上色（黑色、蓝色）
# -----------------------------
cat("\n步骤 3/4: 随机选择网格并分配颜色...\n")

set.seed(42)  # 保证可重复性

# 设置上色比例（可调整）
prop_colored <- 0.15  # 总共15%的网格上色
prop_black <- 0.5     # 上色网格中，50%为黑色，50%为蓝色

n_total <- nrow(grid_china)
n_colored <- round(n_total * prop_colored)
n_black <- round(n_colored * prop_black)
n_blue <- n_colored - n_black

# 随机选择网格索引（均匀分布）
colored_idx <- sample(1:n_total, n_colored, replace = FALSE)
black_idx <- sample(colored_idx, n_black, replace = FALSE)
blue_idx <- setdiff(colored_idx, black_idx)

# 为网格添加颜色标签
grid_china$color <- "none"
grid_china$color[black_idx] <- "black"
grid_china$color[blue_idx] <- "blue"

cat("  ✓ 上色网格统计:\n")
cat("    - 总网格数: ", n_total, "\n", sep = "")
cat("    - 黑色网格: ", n_black, " (", round(n_black/n_total*100, 2), "%)\n", sep = "")
cat("    - 蓝色网格: ", n_blue, " (", round(n_blue/n_total*100, 2), "%)\n", sep = "")
cat("    - 无色网格: ", n_total - n_colored, " (", round((n_total-n_colored)/n_total*100, 2), "%)\n", sep = "")

# -----------------------------
# 6. 绘图
# -----------------------------
cat("\n步骤 4/4: 绘制网格图...\n")

# 转回经纬度坐标系以便绘图
grid_china_lonlat <- sf::st_transform(grid_china, crs = 4326)
china_lonlat <- sf::st_transform(china_proj, crs = 4326)

# 分离不同颜色的网格
grid_black <- grid_china_lonlat[grid_china_lonlat$color == "black", ]
grid_blue <- grid_china_lonlat[grid_china_lonlat$color == "blue", ]

p <- ggplot() +
  # 中国边界底图（浅灰色填充）
  geom_sf(data = china_lonlat, fill = "#f5f5f5", color = "black", linewidth = 0.4) +
  # 黑色网格
  geom_sf(data = grid_black, fill = "black", color = NA, alpha = 0.85) +
  # 蓝色网格
  geom_sf(data = grid_blue, fill = "#1f78b4", color = NA, alpha = 0.85) +
  # 中国边界线（前景，突出显示）
  geom_sf(data = china_lonlat, fill = NA, color = "black", linewidth = 0.4) +
  labs(title = "Grid-based Colored Regions (5km × 5km)",
       subtitle = sprintf("Black: %d grids, Blue: %d grids", n_black, n_blue),
       x = "Longitude (°E)", y = "Latitude (°N)") +
  viz_theme_nature(base_size = 9) +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# 最高质量输出：PNG 2400 dpi + SVG 矢量
ggsave("Drawing/output/grid_colored_regions.png",
       plot = p, width = 5.0, height = 4.0, dpi = 2400, bg = "transparent")
ggsave("Drawing/output/grid_colored_regions.svg",
       plot = p, width = 5.0, height = 4.0, bg = "transparent")

cat("  ✓ 图形已输出（PNG 2400 dpi + SVG，透明背景）\n")

# -----------------------------
# 7. 输出统计表与日志
# -----------------------------
cat("\n步骤 5/5: 输出统计信息...\n")

grid_stats <- data.frame(
  category = c("Total Grids", "Black Grids", "Blue Grids", "Uncolored Grids"),
  count = c(n_total, n_black, n_blue, n_total - n_colored),
  percentage = c(100, 
                 round(n_black/n_total*100, 2),
                 round(n_blue/n_total*100, 2),
                 round((n_total-n_colored)/n_total*100, 2))
)

write.csv(grid_stats, "Drawing/output/grid_statistics.csv", row.names = FALSE)

sink("Drawing/output/grid_log.txt")
cat("中国区域网格化与随机上色日志\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
cat("网格参数:\n")
cat("  - 网格大小: 5km × 5km\n")
cat("  - 投影坐标系: Albers Equal Area Conic (中国标准)\n")
cat("  - 随机种子: 42\n\n")
cat("网格统计:\n")
print(grid_stats)
sink()

cat("\n======================================\n")
cat("绘图完成\n")
cat("======================================\n\n")

cat(sprintf("总网格数: %d\n", n_total))
cat(sprintf("黑色网格: %d (%.2f%%)\n", n_black, n_black/n_total*100))
cat(sprintf("蓝色网格: %d (%.2f%%)\n", n_blue, n_blue/n_total*100))

cat("\n✓ 脚本执行完成!\n\n")

