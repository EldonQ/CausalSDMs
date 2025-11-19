#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 01_occurrence_background_terrain.R
# 功能说明: 绘制单张中国区域底图（地形阴影+高程着色），叠加物种发生点（黑色）
#           与背景点（白色），用于Nature级别期刊配图（英文标注、Arial、1200dpi、透明画布）
# 参考脚本: scripts/13_study_area_and_points.R
# 输入文件: output/01_data_preparation/species_occurrence_cleaned.csv
#          output/03_background_points/background_points.csv
#          earthenvstreams_china/china_boundary.shp
#          Drawing/data/china_dem_30s.tif（若不存在则自动下载 geodata::elevation_30s 中国高程数据）
# 输出文件: Drawing/output/occurrence_background_terrain.png
#          Drawing/output/occurrence_background_terrain.svg
#          Drawing/output/points_summary.csv
#          Drawing/output/elevation_stats.csv
#          Drawing/output/processing_log.txt
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
  "rnaturalearth", # 边界回退
  "ggnewscale",    # 多重fill映射
  "colorspace",    # 科研配色方案（连续/发散/定性，色盲友好）
  "geodata"        # 下载高分辨率高程数据
)
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 统一可视化工具（Nature风格、Arial、1200dpi）
source("scripts/visualization/viz_utils.R")

# -----------------------------
# 2. 目录准备
# -----------------------------
dir.create("Drawing", showWarnings = FALSE, recursive = TRUE)
dir.create("Drawing/output", showWarnings = FALSE, recursive = TRUE)
dir.create("Drawing/data", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("Single Figure: Terrain + Elevation + Points\n")
cat("======================================\n\n")

# -----------------------------
# 3. 读取数据
# -----------------------------
cat("步骤 1/3: 读取数据与底图...\n")

# 物种与背景点数据（假定包含列: lon, lat）
species_data_raw <- read.csv("output/01_data_preparation/species_occurrence_cleaned.csv")
background_data_raw <- read.csv("output/03_background_points/background_points.csv")

cat("  - 原始物种点: ", nrow(species_data_raw), "\n", sep = "")
cat("  - 原始背景点: ", nrow(background_data_raw), "\n", sep = "")

# 中国边界：优先项目内边界，失败回退至 rnaturalearth
viz_china <- viz_read_china("earthenvstreams_china/china_boundary.shp")
if (is.null(viz_china)) {
  viz_china <- rnaturalearth::ne_countries(country = "China", scale = "medium", returnclass = "sf")
}

# -----------------------------
# 3.0 空间采样减少点数量（保持均匀分布）
# -----------------------------
cat("\n  -> 空间采样减少点数量（保持均匀分布）...\n")

# 空间分层采样函数：将区域分为网格，每个网格随机选取若干点（保持均匀分布）
spatial_sample <- function(data, target_n, grid_size = 0.5) {
  # 若目标数量大于等于原数量，直接返回
  if (target_n >= nrow(data)) {
    return(data)
  }
  
  # 创建空间网格（经纬度）
  grid_lon <- floor(data$lon / grid_size)
  grid_lat <- floor(data$lat / grid_size)
  grid_id <- paste(grid_lon, grid_lat, sep = "_")
  
  # 计算每个网格应采样的点数
  n_grids <- length(unique(grid_id))
  points_per_grid <- ceiling(target_n / n_grids)
  
  # 从每个网格中随机采样
  set.seed(42)  # 保证可重复性
  sampled_list <- list()
  for (gid in unique(grid_id)) {
    grid_points <- data[grid_id == gid, ]
    n_sample <- min(points_per_grid, nrow(grid_points))
    if (n_sample > 0) {
      sampled_idx <- sample(seq_len(nrow(grid_points)), n_sample)
      sampled_list[[length(sampled_list) + 1]] <- grid_points[sampled_idx, ]
    }
  }
  sampled <- do.call(rbind, sampled_list)
  
  # 若采样数超过目标，随机选择目标数量
  if (nrow(sampled) > target_n) {
    sampled <- sampled[sample(seq_len(nrow(sampled)), target_n), ]
  }
  
  return(sampled)
}

# 减少背景点（保留30%，约均匀分布）
target_background <- ceiling(nrow(background_data_raw) * 0.3)
background_data <- spatial_sample(background_data_raw, target_background, grid_size = 0.5)

# 减少物种发生点（保留50%，保留更多信息）
target_species <- ceiling(nrow(species_data_raw) * 0.5)
species_data <- spatial_sample(species_data_raw, target_species, grid_size = 0.3)

cat("  ✓ 采样后背景点: ", nrow(background_data), " (", 
    round(nrow(background_data)/nrow(background_data_raw)*100, 1), "%)\n", sep = "")
cat("  ✓ 采样后物种点: ", nrow(species_data), " (", 
    round(nrow(species_data)/nrow(species_data_raw)*100, 1), "%)\n", sep = "")

# 读取高程，构建地形阴影（hillshade）
# 优先使用更高分辨率的 SRTM 90m 数据，失败回退至 30s，最后回退至项目数据
dem_path_srtm <- file.path("Drawing", "data", "china_dem_srtm_90m.tif")
dem_path_30s <- file.path("Drawing", "data", "china_dem_30s.tif")
dem_source_used <- NA_character_
elev <- NULL

# 尝试1: SRTM 30m（最高质量，若 elevation_3s 支持更高分辨率则使用）
if (!file.exists(dem_path_srtm)) {
  cat("  - 尝试下载最高分辨率高程 (geodata::elevation_3s, 尽可能精细) ...\n")
  elev_try <- tryCatch({
    # elevation_3s 会根据区域大小自动选择最优分辨率（SRTM 90m 或更高）
    geodata::elevation_3s(country = "China", path = "Drawing/data", mask = TRUE)
  }, error = function(e) {
    cat("  ! 高分辨率高程下载失败: ", e$message, "\n", sep = "")
    NULL
  })
  if (!is.null(elev_try)) {
    cat("  ✓ 高分辨率高程下载完成，写入本地...\n")
    terra::writeRaster(elev_try, dem_path_srtm, overwrite = TRUE)
    elev <- elev_try
    dem_source_used <- "geodata::elevation_3s (highest available resolution, ~90m)"
  }
}
if (is.null(elev) && file.exists(dem_path_srtm)) {
  cat("  ✓ 加载本地高分辨率高程: ", dem_path_srtm, "\n", sep = "")
  elev <- terra::rast(dem_path_srtm)
  dem_source_used <- "geodata::elevation_3s (highest available resolution, ~90m)"
}

# 尝试2: 30s 数据（中等质量）
if (is.null(elev) && !file.exists(dem_path_30s)) {
  cat("  - SRTM 不可用，回退至 30s 数据 (geodata::elevation_30s) ...\n")
  elev_try <- tryCatch({
    geodata::elevation_30s(country = "China", path = "Drawing/data", mask = TRUE, keepzip = FALSE)
  }, error = function(e) {
    cat("  ! 30s 数据下载失败: ", e$message, "\n", sep = "")
    NULL
  })
  if (!is.null(elev_try)) {
    cat("  ✓ 30s 数据下载完成，写入本地...\n")
    terra::writeRaster(elev_try, dem_path_30s, overwrite = TRUE)
    elev <- elev_try
    dem_source_used <- "geodata::elevation_30s (30 arc-second, ~1km)"
  }
}
if (is.null(elev) && file.exists(dem_path_30s)) {
  cat("  ✓ 加载本地 30s 高程: ", dem_path_30s, "\n", sep = "")
  elev <- terra::rast(dem_path_30s)
  dem_source_used <- "geodata::elevation_30s (30 arc-second, ~1km)"
}

# 备用: 项目内高程
if (is.null(elev)) {
  cat("  ! 所有高分辨率数据不可用，使用备用: earthenvstreams_china/elevation.tif\n")
  elev <- terra::rast("earthenvstreams_china/elevation.tif")
  dem_source_used <- "earthenvstreams_china/elevation.tif (project backup)"
}
china_vect <- terra::vect(viz_china)
# 投影统一
if (!terra::same.crs(elev, china_vect)) {
  # 若坐标参考系不同，则将中国边界重投影至高程栅格的CRS
  china_vect <- terra::project(china_vect, terra::crs(elev, proj = TRUE))
}
# 裁剪并掩膜至中国区域
elev_china <- terra::mask(terra::crop(elev, china_vect), china_vect)
# 若为多波段高程，先合成为单层（取平均），以满足后续地形分析函数需求
elev_single <- if (terra::nlyr(elev_china) > 1) {
  suppressWarnings(terra::mean(elev_china, na.rm = TRUE))
} else {
  elev_china
}
# 完全不降采样，保持最高质量（性能强劲无需聚合）
cat("  ✓ 保持原始高程分辨率（无降采样），最大化细节...\n")
elev_single_ag <- elev_single

# 计算坡度与坡向（弧度），生成多方位阴影（增强渲染）
slope <- terra::terrain(elev_single_ag, v = "slope", unit = "radians")
aspect <- terra::terrain(elev_single_ag, v = "aspect", unit = "radians")
azimuths <- c(45, 90, 135, 225, 315)
hill_list <- lapply(azimuths, function(az) terra::shade(slope, aspect, 45, az))
hill <- Reduce("+", hill_list) / length(hill_list)

# 转为数据框（经纬度 + 值）
hill_df <- as.data.frame(hill, xy = TRUE, na.rm = TRUE)
colnames(hill_df) <- c("lon", "lat", "hill")
# 归一化到[0,1]
if (nrow(hill_df) > 0) {
  mx_h <- suppressWarnings(max(hill_df$hill, na.rm = TRUE))
  if (is.finite(mx_h) && mx_h > 0) hill_df$hill <- hill_df$hill / mx_h
}
# 完全不采样，保留所有栅格点（最高质量）
cat("  ✓ 地形阴影栅格: ", nrow(hill_df), " 点（全部保留）\n", sep = "")

elev_df <- as.data.frame(elev_single_ag, xy = TRUE, na.rm = TRUE)
colnames(elev_df) <- c("lon", "lat", "elev")
if (nrow(elev_df) > 0) {
  mn_e <- suppressWarnings(min(elev_df$elev, na.rm = TRUE))
  mx_e <- suppressWarnings(max(elev_df$elev, na.rm = TRUE))
  if (is.finite(mn_e) && is.finite(mx_e) && mx_e > mn_e) {
    elev_df$elev_norm <- (elev_df$elev - mn_e) / (mx_e - mn_e)
  } else {
    elev_df$elev_norm <- 0
  }
}
# 完全不采样，保留所有栅格点（最高质量）
cat("  ✓ 高程栅格: ", nrow(elev_df), " 点（全部保留）\n", sep = "")

# -----------------------------
# 3.1 河网数据处理（HydroRIVERS矢量数据，高质量分级渲染）
# -----------------------------
cat("\n步骤 1b: 处理河网数据（HydroRIVERS亚洲）...\n")

# 读取HydroRIVERS亚洲河网矢量数据
hydrorivers_path <- "E:/HydroRIVERS_v10_as_shp/HydroRIVERS_v10_as_shp"
shp_files <- list.files(hydrorivers_path, pattern = "*.shp$", full.names = TRUE)

if (length(shp_files) == 0) {
  stop("未找到HydroRIVERS shp文件，请检查路径: ", hydrorivers_path)
}

# 读取第一个shp文件（通常为主河网文件）
rivers_raw <- sf::st_read(shp_files[1], quiet = TRUE)
cat("  ✓ 读取河网矢量: ", nrow(rivers_raw), " 条河流\n", sep = "")

# 转换坐标系为中国边界相同的CRS
rivers_crs <- sf::st_transform(rivers_raw, sf::st_crs(viz_china))

# 裁剪至中国范围（保留与中国相交的河流）
rivers_china <- sf::st_intersection(rivers_crs, viz_china)
cat("  ✓ 裁剪至中国范围: ", nrow(rivers_china), " 条河流\n", sep = "")

# 使用 Strahler 等级（ORD_STRA 字段）或其他等级字段进行分级
# HydroRIVERS 常用字段: ORD_STRA (Strahler), ORD_CLAS (分类), UPLAND_SKM (上游面积)
if ("ORD_STRA" %in% names(rivers_china)) {
  rivers_china$river_order <- rivers_china$ORD_STRA
  order_field <- "ORD_STRA"
} else if ("ORD_CLAS" %in% names(rivers_china)) {
  rivers_china$river_order <- rivers_china$ORD_CLAS
  order_field <- "ORD_CLAS"
} else if ("UPLAND_SKM" %in% names(rivers_china)) {
  # 使用上游面积分级（对数变换）
  rivers_china$river_order <- cut(
    log10(rivers_china$UPLAND_SKM + 1),
    breaks = 5,
    labels = 1:5
  )
  order_field <- "UPLAND_SKM (derived)"
} else {
  # 备用：使用行索引简单分级
  rivers_china$river_order <- 3
  order_field <- "default"
}

# 转换为因子并统计
rivers_china$river_order <- as.factor(rivers_china$river_order)
order_levels <- levels(rivers_china$river_order)
n_classes <- length(order_levels)

# 分配蓝色层次（浅→深）
all_colors <- c("#B3D9FF", "#66B3FF", "#1A8CFF", "#0066CC", "#003D7A")
if (n_classes <= 5) {
  river_colors <- all_colors[1:n_classes]
} else {
  # 若超过5级，动态生成渐变色
  river_colors <- colorRampPalette(c("#B3D9FF", "#003D7A"))(n_classes)
}
names(river_colors) <- order_levels

cat("  ✓ 河流分级字段: ", order_field, "\n", sep = "")
cat("  ✓ 河流等级数: ", n_classes, "\n", sep = "")
for (i in seq_along(order_levels)) {
  ord <- order_levels[i]
  cnt <- sum(rivers_china$river_order == ord, na.rm = TRUE)
  cat("    - 等级 ", ord, ": ", cnt, " 条 (", river_colors[i], ")\n", sep = "")
}

# -----------------------------
# 3.2 补充流域累积栅格（填补HydroRIVERS未覆盖区域）
# -----------------------------
cat("\n步骤 1c: 补充栅格河网（填补数据空白）...\n")

# 读取 flow_acc 作为补充数据（用于覆盖HydroRIVERS缺失区域）
flow_acc_raw <- terra::rast("earthenvstreams_china/flow_acc.tif")
if (terra::nlyr(flow_acc_raw) > 1) {
  flow_acc_raw <- flow_acc_raw[[2]]
}
if (!terra::same.crs(flow_acc_raw, china_vect)) {
  flow_acc_raw <- terra::project(flow_acc_raw, terra::crs(elev_single_ag))
}
flow_acc <- terra::mask(terra::crop(flow_acc_raw, china_vect), china_vect)

# 仅保留中高累积量河网（避免过于密集）
flow_acc[flow_acc <= 200] <- NA

# 对数变换
flow_acc_log <- log10(flow_acc + 1)

# 转为数据框
flow_df <- as.data.frame(flow_acc_log, xy = TRUE, na.rm = TRUE)
colnames(flow_df) <- c("lon", "lat", "flow")

if (nrow(flow_df) > 0) {
  # 归一化
  mn_f <- suppressWarnings(min(flow_df$flow, na.rm = TRUE))
  mx_f <- suppressWarnings(max(flow_df$flow, na.rm = TRUE))
  if (is.finite(mn_f) && is.finite(mx_f) && mx_f > mn_f) {
    flow_df$flow_norm <- (flow_df$flow - mn_f) / (mx_f - mn_f)
  } else {
    flow_df$flow_norm <- 0.5
  }
  
  # 统一蓝色（浅蓝色填充，与矢量河网协调）
  flow_df$river_color <- "#80B3D9"  # 中等蓝色，介于1-3级之间
  
  cat("  ✓ 补充栅格河网: ", nrow(flow_df), " 点（填补空白区域）\n", sep = "")
} else {
  flow_df <- NULL
  cat("  ! 无补充栅格数据\n")
}

cat("  ✓ 物种分布点: ", nrow(species_data), "\n", sep = "")
cat("  ✓ 背景点: ", nrow(background_data), "\n", sep = "")
cat("  ✓ 地形网格(阴影): ", nrow(hill_df), "\n", sep = "")
cat("  ✓ 高程网格: ", nrow(elev_df), "\n", sep = "")
cat("  ✓ 高程数据来源: ", dem_source_used, "\n", sep = "")

# -----------------------------
# 4. 绘图（单张）
#    要求：
#    - 同时展示发生点（黑）与背景点（白）
#    - 中国区域背景渲染地形阴影 + 高程着色
#    - 英文标注、Arial、透明画布、1200dpi，单图输出
# -----------------------------
cat("\n步骤 2/3: 绘制单图...\n")

# 统一点大小（背景点与发生点一致）
point_size <- 0.3

# 专业地形配色：基于 ETOPO1 / Terrain 经典配色方案（蓝-绿-黄-棕-白，适合高程渲染）
# 低海拔→高海拔：深绿 → 浅绿 → 黄 → 橙棕 → 棕红 → 灰白
pal_elev <- grDevices::colorRampPalette(c(
  "#2A6F3F",  # 低海拔：深绿（平原、盆地）
  "#4FA15C",  # 浅绿（丘陵）
  "#8FBC8F",  # 黄绿（低山）
  "#C8B56B",  # 黄（中山）
  "#D4955A",  # 橙（高山）
  "#B8704F",  # 棕（高原）
  "#A0524F",  # 棕红（极高山）
  "#D9D9D9"   # 高海拔：灰白（雪线、冰川）
))(256)

# 构建基础图层
p <- ggplot() +
  # 第一层：地形阴影（灰度，增强对比度）
  geom_raster(data = hill_df, aes(x = lon, y = lat, fill = hill)) +
  scale_fill_gradient(name = NULL, low = "#1a1a1a", high = "#fafafa", guide = "none") +
  ggnewscale::new_scale_fill() +
  # 第二层：高程着色（半透明叠加，增加透明度以更好融合阴影）
  geom_raster(data = elev_df, aes(x = lon, y = lat, fill = elev_norm), alpha = 0.5) +
  scale_fill_gradientn(name = NULL, colours = pal_elev, guide = "none")

# 第三层A：补充栅格河网（底层，填补HydroRIVERS未覆盖区域）
if (!is.null(flow_df) && nrow(flow_df) > 0) {
  p <- p + geom_raster(
    data = flow_df,
    aes(x = lon, y = lat),
    fill = "#80B3D9",  # 中等蓝色
    alpha = 0.5        # 半透明，作为底图
  )
}

# 第三层B：河网矢量分级渲染（上层，蓝色层次，按等级从小到大绘制）
# 使用 geom_sf 绘制矢量河流，按等级分配颜色与线宽
# 先绘制小河流，后绘制大河流，确保大河流覆盖在上层
linewidth_values <- seq(0.15, 0.6, length.out = n_classes)
for (i in seq_along(order_levels)) {
  ord <- order_levels[i]
  river_subset <- rivers_china[rivers_china$river_order == ord, ]
  if (nrow(river_subset) > 0) {
    p <- p + geom_sf(
      data = river_subset,
      color = river_colors[i],
      linewidth = linewidth_values[i],
      alpha = 0.85
    )
  }
}

# 继续添加其余图层
p <- p +
  # 中国边界线
  geom_sf(data = viz_china, fill = NA, color = "black", linewidth = 0.3) +
  # 背景点：两层叠加实现白色+黑边（严格统一大小）
  # 第1层：黑色底（作为极细边框，仅比白色大 0.05）
  geom_point(data = background_data, aes(x = lon, y = lat),
             shape = 19, color = "black", size = point_size + 0.05, alpha = 1.0) +
  # 第2层：白色填充（与发生点完全相同的 size）
  geom_point(data = background_data, aes(x = lon, y = lat),
             shape = 19, color = "white", size = point_size, alpha = 1.0) +
  # 发生点：黑色实心圆（shape=19），与背景点白色层大小严格一致
  geom_point(data = species_data, aes(x = lon, y = lat),
             shape = 19, color = "black", size = point_size, alpha = 1.0) +
  labs(title = "Species Presence and Background Points (Terrain, Elevation, and River Network)",
       x = "Longitude (\u00B0E)", y = "Latitude (\u00B0N)") +
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

# 最高质量输出：PNG 2400 dpi（超高分辨率，适合印刷）+ SVG 矢量
ggsave("Drawing/output/occurrence_background_terrain.png",
       plot = p, width = 4.8, height = 3.6, dpi = 2400, bg = "transparent")
ggsave("Drawing/output/occurrence_background_terrain.svg",
       plot = p, width = 4.8, height = 3.6, bg = "transparent")
cat("  ✓ 输出分辨率: 2400 dpi (超高质量印刷级)\n")

cat("  ✓ 单图已输出（PNG+SVG，透明背景）\n")

# -----------------------------
# 5. 输出统计表与日志
# -----------------------------
cat("\n步骤 3/3: 输出统计信息与日志...\n")

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

elev_stats <- data.frame(
  metric = c("source", "elev_min", "elev_max", "elev_median"),
  value = c(
    dem_source_used,
    sprintf("%.4f", suppressWarnings(min(elev_df$elev, na.rm = TRUE))),
    sprintf("%.4f", suppressWarnings(max(elev_df$elev, na.rm = TRUE))),
    sprintf("%.4f", suppressWarnings(stats::median(elev_df$elev, na.rm = TRUE)))
  )
)

write.csv(summary_data, "Drawing/output/points_summary.csv", row.names = FALSE)
write.csv(elev_stats, "Drawing/output/elevation_stats.csv", row.names = FALSE)

sink("Drawing/output/processing_log.txt")
cat("Single Figure (Terrain + Elevation + Points) Log\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
cat("Species and background points counts:\n")
print(summary_data)
cat("\nElevation statistics:\n")
print(elev_stats)
sink()

cat("\n======================================\n")
cat("绘图完成\n")
cat("======================================\n\n")

cat(sprintf("物种分布点: %d\n", nrow(species_data)))
cat(sprintf("背景点: %d\n", nrow(background_data)))
cat(sprintf("背景点/物种点比例: %.2f:1\n", nrow(background_data)/nrow(species_data)))

cat("\n✓ 脚本执行完成!\n\n")


