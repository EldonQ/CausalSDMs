################################################################################
# 脚本名称：01_data_preparation_NEW.R
# 功能描述：基于 binomial/long/lat 的物种出现数据清洗（简化版）
# 输入文件：Carassius_auratusOCC/Carassius_auratus.csv (binomial,long,lat)
# 输出文件：output/01_data_preparation/species_occurrence_cleaned.csv
#          output/01_data_preparation/species_occurrence_cleaned.shp
# 作者：Nature级别科研项目
# 日期：2025-10-19
# 研究区域：中国水系（全国）
################################################################################

# 清空环境
rm(list = ls())
gc()

# 设置工作目录
setwd("E:/SDM01")

# 加载必需的R包
required_packages <- c("tidyverse", "sf", "terra", "ggplot2")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 日志函数
log_message <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
}

################################################################################
# 可调参数配置
################################################################################

# 空间稀疏化参数（单位：度，1度≈111km）
# 作用：避免空间自相关，减少采样偏差，提高模型泛化能力
# 方法：将区域划分为网格，每个网格只保留1个观测点
SPATIAL_THINNING_GRID <- 0.09  # 约10km网格，适合全国尺度研究

################################################################################
# 创建输出目录
################################################################################

if (!dir.exists("output")) dir.create("output")
if (!dir.exists("output/01_data_preparation")) dir.create("output/01_data_preparation", recursive = TRUE)
if (!dir.exists("figures/01_data_preparation")) dir.create("figures/01_data_preparation", recursive = TRUE)

log_message("======================================")
log_message("开始处理Carassius auratus物种出现数据")
log_message("======================================")

################################################################################
# 第一步：读取GBIF数据
################################################################################

log_message("\n步骤 1/8: 读取数据...")

occ_file <- "Carassius_auratusOCC/Carassius_auratus.csv"

if (!file.exists(occ_file)) {
  stop(paste0("错误：找不到文件 ", occ_file))
}

# 读取数据
occ_raw <- read.csv(occ_file, stringsAsFactors = FALSE)
log_message(paste0("  - 原始记录数: ", nrow(occ_raw)))
log_message(paste0("  - 字段数: ", ncol(occ_raw)))

# 查看数据结构
log_message("  - 主要字段:")
log_message(paste0("    ", paste(head(names(occ_raw), 10), collapse = ", ")))

# 统一列名到 species/lon/lat，并进行基础清洗（仅适配 binomial/long/lat 简化格式）
occ <- occ_raw %>%
  rename(
    species = binomial,
    lon = long,
    lat = lat
  ) %>%
  mutate(
    lon = as.numeric(lon),
    lat = as.numeric(lat),
    source = "OCC_CSV"
  ) %>%
  # 移除无效/缺失坐标与越界坐标
  filter(!is.na(lon) & !is.na(lat) &
         lon >= -180 & lon <= 180 &
         lat >= -90 & lat <= 90)

log_message(paste0("  - 移除缺失/无效坐标后: ", nrow(occ), " 条记录"))

# 数据源已在上一步设置为 OCC_CSV

################################################################################
# 第三步：空间筛选 - 限制到中国范围
################################################################################

log_message("\n步骤 3/8: 空间筛选（限制到中国范围）...")

# 读取中国边界
china_boundary <- vect("earthenvstreams_china/china_boundary.shp")

# 转换出现点为空间对象
occ_sp <- vect(occ, geom = c("lon", "lat"), crs = "EPSG:4326")

# 使用is.related判断点是否在多边形内
log_message("  - 执行空间叠加分析...")
is_in_china <- relate(occ_sp, china_boundary, "intersects")

# is.related返回矩阵，取第一列（因为只有一个多边形）
if (is.matrix(is_in_china)) {
  is_in_china <- is_in_china[, 1]
}

occ_china <- occ[is_in_china, ]
log_message(paste0("  - 中国境内记录: ", nrow(occ_china), " 条"))
log_message(paste0("  - 移除境外记录: ", nrow(occ) - nrow(occ_china), " 条"))

# 更新occ为中国境内的数据
occ <- occ_china

################################################################################
# 第四步：基本坐标检查（简化清洗）
################################################################################

log_message("\n步骤 4/8: 基本坐标检查...")

# 由于已经限制在中国境内，只需进行基本检查
occ_clean <- occ %>%
  # 移除经纬度完全相同的记录（可能是数据错误）
  filter(lon != lat) %>%
  # 移除整数坐标（可能是粗略估计）
  filter(!(lon == round(lon) & lat == round(lat)))

n_removed <- nrow(occ) - nrow(occ_clean)
log_message(paste0("  - 移除异常坐标: ", n_removed, " 条"))
log_message(paste0("  - 剩余记录: ", nrow(occ_clean), " 条"))

################################################################################
# 第五步：移除重复坐标（SDM必要步骤）
################################################################################

log_message("\n步骤 5/8: 移除重复坐标...")

# 统计重复情况
coord_stats <- occ_clean %>%
  group_by(lon, lat) %>%
  summarise(n_records = n(), .groups = "drop") %>%
  arrange(desc(n_records))

log_message(paste0("  - 唯一坐标点: ", nrow(coord_stats), " 个"))
log_message(paste0("  - 平均每坐标记录数: ", round(mean(coord_stats$n_records), 2)))
log_message(paste0("  - 最多重复坐标: ", max(coord_stats$n_records), " 条记录"))

# 对于SDM，每个坐标只保留一条记录（简单去重）
occ_unique <- occ_clean %>%
  arrange(lon, lat) %>%
  distinct(lon, lat, .keep_all = TRUE)

log_message(paste0("  - 去重后: ", nrow(occ_unique), " 条记录"))
log_message(paste0("  - 移除重复: ", nrow(occ_clean) - nrow(occ_unique), " 条"))
log_message("  - 说明：对于SDM，同一坐标的多个记录不提供额外空间信息")

################################################################################
# 第六步：空间稀疏化
################################################################################

log_message("\n步骤 6/8: 空间稀疏化...")

# 创建网格并稀疏化
log_message(paste0("  - 使用网格大小: ", SPATIAL_THINNING_GRID, "° (约", 
                  round(SPATIAL_THINNING_GRID * 111, 1), " km)"))

# 使用简单的网格方法稀疏化
occ_thin <- occ_unique %>%
  mutate(
    grid_lon = floor(lon / SPATIAL_THINNING_GRID),
    grid_lat = floor(lat / SPATIAL_THINNING_GRID)
  ) %>%
  group_by(grid_lon, grid_lat) %>%
  slice_sample(n = 1) %>%  # 每个网格随机保留1个点
  ungroup() %>%
  select(-grid_lon, -grid_lat)

log_message(paste0("  - 稀疏化后: ", nrow(occ_thin), " 条记录"))
log_message(paste0("  - 移除冗余: ", nrow(occ_unique) - nrow(occ_thin), " 条"))

################################################################################
# 第七步：最终检查和统计
################################################################################

log_message("\n步骤 7/8: 最终检查和统计...")

# 坐标范围
lon_range <- range(occ_thin$lon)
lat_range <- range(occ_thin$lat)

log_message("  - 最终数据统计:")
log_message(paste0("    总记录数: ", nrow(occ_thin)))
log_message(paste0("    经度范围: ", round(lon_range[1], 2), " ~ ", round(lon_range[2], 2)))
log_message(paste0("    纬度范围: ", round(lat_range[1], 2), " ~ ", round(lat_range[2], 2)))

# 物种统计
species_counts <- occ_thin %>%
  count(species, sort = TRUE)

log_message(paste0("    涉及物种数: ", nrow(species_counts)))
log_message("    前5个物种:")
for (i in seq_len(min(5, nrow(species_counts)))) {
  log_message(paste0("      ", i, ". ", species_counts$species[i], " (", 
                    species_counts$n[i], " 条)"))
}

# 年份统计
if ("year" %in% names(occ_thin)) {
  year_counts <- occ_thin %>%
    filter(!is.na(year)) %>%
    count(year, sort = TRUE)
  
  if (nrow(year_counts) > 0) {
    log_message(paste0("    时间跨度: ", min(year_counts$year, na.rm = TRUE), 
                      " - ", max(year_counts$year, na.rm = TRUE)))
  }
}

################################################################################
# 第八步：保存结果
################################################################################

log_message("\n步骤 8/8: 保存结果...")

# 保存CSV文件
output_csv <- "output/01_data_preparation/species_occurrence_cleaned.csv"
write.csv(occ_thin, output_csv, row.names = FALSE)
log_message(paste0("  ✓ 已保存CSV: ", output_csv))

# 保存shapefile
output_shp <- "output/01_data_preparation/species_occurrence_cleaned.shp"
occ_sf <- st_as_sf(occ_thin, coords = c("lon", "lat"), crs = 4326)
st_write(occ_sf, output_shp, delete_dsn = TRUE, quiet = TRUE)
log_message(paste0("  ✓ 已保存Shapefile: ", output_shp))

# 保存处理日志
log_file <- "output/01_data_preparation/processing_log.txt"
sink(log_file)
cat("Carassius auratus物种出现数据处理日志\n")
cat("处理时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("数据源: Carassius_auratusOCC/Carassius_auratus.csv (GBIF)\n\n")
cat("处理流程:\n")
cat("  1. 原始记录数:", nrow(occ_raw), "\n")
cat("  2. 移除缺失/无效坐标后:", length(is_in_china), "\n")
cat("  3. 限制到中国境内后:", nrow(occ_china), "\n")
cat("  4. 基本坐标检查后:", nrow(occ_clean), "\n")
cat("  5. 去除重复记录后:", nrow(occ_unique), "\n")
cat("  6. 空间稀疏化后:", nrow(occ_thin), "\n\n")
cat("最终数据:\n")
cat("  记录数:", nrow(occ_thin), "\n")
cat("  物种数:", nrow(species_counts), "\n")
cat("  经度范围:", round(lon_range[1], 2), "~", round(lon_range[2], 2), "\n")
cat("  纬度范围:", round(lat_range[1], 2), "~", round(lat_range[2], 2), "\n\n")
cat("处理参数:\n")
cat("  空间稀疏化网格:", SPATIAL_THINNING_GRID, "° (约", 
    round(SPATIAL_THINNING_GRID * 111, 1), " km)\n")
cat("  研究区域: 中国境内\n")
sink()
log_message(paste0("  ✓ 已保存处理日志: ", log_file))

# 保存物种统计
species_stats <- "output/01_data_preparation/species_statistics.csv"
write.csv(species_counts, species_stats, row.names = FALSE)
log_message(paste0("  ✓ 已保存物种统计: ", species_stats))

################################################################################
# 第九步：生成简单的可视化（可选）
################################################################################

log_message("\n生成简单可视化...")

tryCatch({
  # 读取中国边界并转换为sf对象（用于绘图）
  china_sf <- st_as_sf(china_boundary)
  
  # 动态物种名（用于标题），默认取出现频次最高的第一个物种
  species_label <- tryCatch({
    if (exists("species_counts") && nrow(species_counts) > 0) species_counts$species[1] else "Target species"
  }, error = function(e) "Target species")
  
  # 绘制分布点图 - Nature期刊风格
  p <- ggplot() +
    geom_sf(data = china_sf, fill = "gray95", color = "black", size = 0.3) +
    geom_point(data = occ_thin, aes(x = lon, y = lat), 
               color = "#D62728", size = 1.5, alpha = 0.7) +
    theme_minimal(base_family = "Arial") +
    theme(
      panel.grid = element_line(color = "gray90", size = 0.2),
      panel.border = element_rect(color = "black", fill = NA, size = 0.5),
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    labs(
      title = paste0(species_label, " Species Occurrence in China"),
      x = "Longitude (°E)",
      y = "Latitude (°N)",
      caption = paste0("n = ", nrow(occ_thin), " occurrence records")
    ) +
    # 根据中国边界自动确定范围（参考viz_00脚本）
    coord_sf(expand = FALSE)
  
  # 保存PNG和PDF格式
  ggsave(
    filename = "figures/01_data_preparation/species_occurrence_map.png",
    plot = p,
    width = 10,
    height = 8,
    dpi = 1200,  # Nature期刊要求 ≥1200 dpi
    bg = "white"
  )
  
  log_message("  ✓ 已保存分布图: figures/01_data_preparation/species_occurrence_map.png")
  
  
}, error = function(e) {
  log_message(paste0("  - 可视化失败: ", e$message))
})

################################################################################
# 完成
################################################################################

log_message("\n======================================")
log_message("数据准备完成！")
log_message("======================================")
log_message(paste0("最终输出: ", nrow(occ_thin), " 条清洗后的物种出现记录"))
log_message(paste0("输出文件夹: output/01_data_preparation/"))
log_message("======================================\n")

log_message("脚本执行完毕！")

