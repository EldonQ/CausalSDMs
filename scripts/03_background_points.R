#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 03_background_points.R
# 功能说明: 在中国河网上生成随机均匀分布的背景点
# 策略: 河网随机采样 (Poisson-disk on river mask)，使用白名单/初筛后的合格变量
# 输入文件:
#   - output/02_env_extraction/occurrence_with_env_complete.csv
#   - output/01b_variable_prescreening/qualified_variables.csv
# 输出文件:
#   - output/03_background_points/background_points.csv
#   - output/03_background_points/combined_presence_absence.csv
# 作者: Nature级别科研项目
# 日期: 2025-10-20
# ==============================================================================

# 清空环境
rm(list = ls())
gc()

# 设置工作目录
setwd("E:/CausalSDMs")

# 加载必要的包
packages <- c("tidyverse", "raster", "sf", "sp", "dismo")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 统一可视化工具（Nature风格）
source("scripts/visualization/viz_utils.R")

# 创建输出目录
if (!dir.exists("output/03_background_points")) {
  dir.create("output/03_background_points", recursive = TRUE)
}

cat("======================================\n")
cat("生成背景点（中国河网）\n")
cat("策略: 全国水网均匀分布 (Poisson-disk on river mask)\n")
cat("======================================\n\n")

# 采样参数
TARGET_RATIO_BG_TO_PRES <- 5 # 背景点与出现点比例
LAEA_CRS <- "+proj=laea +lat_0=35 +lon_0=105 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"

# ------------------------------------------------------------------------------
# 1. 读取物种出现数据
# ------------------------------------------------------------------------------
cat("步骤 1/6: 读取物种出现数据...\n")

presence_data <- read.csv("output/02_env_extraction/occurrence_with_env_complete.csv")
cat("  - 出现点数量: ", nrow(presence_data), "\n", sep = "")
cat("  - 物种数量: ", length(unique(presence_data$species)), "\n", sep = "")

# ------------------------------------------------------------------------------
# 2. 读取初筛后的变量列表
# ------------------------------------------------------------------------------
cat("\n步骤 2/6: 读取初筛后的变量列表...\n")

qualified_vars <- read.csv("output/01b_variable_prescreening/qualified_variables.csv")
cat("  - 合格变量数: ", nrow(qualified_vars), "\n", sep = "")

# 按文件组织变量
vars_by_file <- qualified_vars %>%
  arrange(file, band) %>%
  group_by(file) %>%
  summarise(
    bands = list(band),
    var_names = list(variable),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 3. 构建河网采样掩膜
# ------------------------------------------------------------------------------
cat("\n步骤 3/6: 构建河网采样掩膜...\n")

# 读取flow_acc作为河网掩膜
flow_acc_brick <- brick("earthenvstreams_china/flow_acc.tif")
flow_acc_layer <- flow_acc_brick[[2]] # 第2波段是flow accumulation

# 创建河网掩膜（flow_acc > 0 的区域）
sampling_mask <- flow_acc_layer
mask_values <- getValues(sampling_mask)

# 转换NoData为NA
mask_values[mask_values == -127] <- NA
mask_values[mask_values == -999] <- NA
mask_values[mask_values == -9999] <- NA
mask_values[mask_values < -100] <- NA

# 只保留flow_acc > 0的区域（河网）
mask_values[!is.na(mask_values) & mask_values > 0] <- 1
mask_values[!is.na(mask_values) & mask_values <= 0] <- NA

sampling_mask <- setValues(sampling_mask, mask_values)

cat("  - 采样域: 中国河网（flow_acc > 0）\n")
cat("  - 河网有效像元数: ", sum(!is.na(mask_values)), "\n", sep = "")

# 可视化采样掩膜
try(
  {
    dir.create("figures/03_background_points", showWarnings = FALSE, recursive = TRUE)
    sampling_mask_spat <- terra::rast(sampling_mask)
    viz_save_raster_map(
      r = sampling_mask_spat,
      out_base = "figures/03_background_points/river_sampling_mask",
      title = "Sampling Mask (flow_acc > 0)",
      palette = "viridis",
      q_limits = c(0.00, 1.00),
      china_path = "earthenvstreams_china/china_boundary.shp",
      width_in = 6, height_in = 4
    )
    cat("  ✓ 采样掩膜预览: figures/03_background_points/river_sampling_mask.png/svg\n")
  },
  silent = TRUE
)


# ------------------------------------------------------------------------------
# 4. 生成背景点
#    策略: LAEA等面积域上的 Poisson-disk (蓝噪声) 采样，仅在河网掩膜内
# ------------------------------------------------------------------------------
cat("\n步骤 4/6: 在全国河网上分层均匀生成背景点 (Poisson-disk)...\n")

n_presence <- nrow(presence_data)
n_background <- n_presence * TARGET_RATIO_BG_TO_PRES

cat("  - 出现点数量: ", n_presence, "\n", sep = "")
cat("  - 目标背景点数量: ", n_background, " (5倍出现点)\n", sep = "")

# 中国边界用于估计尺度
china_sf <- sf::st_read("earthenvstreams_china/china_boundary.shp", quiet = TRUE)
china_laea <- sf::st_transform(china_sf, crs = LAEA_CRS)
bbox_laea <- sf::st_as_sfc(sf::st_bbox(china_laea))
area_bbox <- as.numeric(sf::st_area(bbox_laea))

# 预采样候选（仅从河网非NA像元中抽样）
cand_n <- min(2e6, max(50 * n_background, 50000))
cand <- raster::sampleRandom(sampling_mask, size = cand_n, xy = TRUE, na.rm = TRUE, sp = TRUE)
if (is.null(cand) || length(cand) == 0) stop("河网掩膜候选采样失败")

cand_sf <- sf::st_as_sf(cand)
cand_laea <- sf::st_transform(cand_sf, crs = LAEA_CRS)
coords <- as.matrix(sf::st_coordinates(cand_laea))

if (nrow(coords) < n_background) {
  stop("候选点数量不足，请增大 cand_n 或检查掩膜")
}

# 估算起始最小间距（米）；随后自适应放宽
d_min <- sqrt(area_bbox / (n_background * 10)) # 经验系数10
d_min <- max(d_min, 5000) # 不低于5km

select_with_dmin <- function(coords_laea, target_n, d_m) {
  idx_order <- sample(seq_len(nrow(coords_laea)))
  sel <- integer(0)
  for (ii in idx_order) {
    if (length(sel) == 0) {
      sel <- ii
    } else {
      # 计算到已选点的最小距离
      dx <- coords_laea[sel, 1] - coords_laea[ii, 1]
      dy <- coords_laea[sel, 2] - coords_laea[ii, 2]
      if (min(dx * dx + dy * dy) >= d_m * d_m) sel <- c(sel, ii)
    }
    if (length(sel) >= target_n) break
  }
  return(sel)
}

# 自适应寻找能达到目标数量的最小间距
sel_idx <- integer(0)
for (fac in c(1.0, 0.9, 0.8, 0.7, 0.6, 0.5)) {
  d_try <- d_min * fac
  cat("    尝试最小间距: ", round(d_try / 1000, 1), " km ... ")
  sel_idx <- select_with_dmin(coords, n_background, d_try)
  cat("选中 ", length(sel_idx), " 个\n", sep = "")
  if (length(sel_idx) >= n_background) break
}

if (length(sel_idx) == 0) stop("Poisson-disk 选择失败，请检查掩膜或参数")
if (length(sel_idx) > n_background) sel_idx <- sel_idx[seq_len(n_background)]

# 回到经纬度
chosen_laea <- cand_laea[sel_idx, , drop = FALSE]
chosen_ll <- sf::st_transform(chosen_laea, crs = 4326)
background_coords <- as.matrix(sf::st_coordinates(chosen_ll))

if (is.null(background_coords)) stop("未能生成任何背景点，请检查河网掩膜或参数设置。")
cat("  - 成功生成背景点: ", nrow(background_coords), "\n", sep = "")

# ------------------------------------------------------------------------------
# 5. 提取背景点的环境变量
# ------------------------------------------------------------------------------
cat("\n步骤 5/6: 提取背景点的环境变量...\n")

# 创建背景点数据框
background_df <- data.frame(
  id = (n_presence + 1):(n_presence + nrow(background_coords)),
  species = "background",
  lon = background_coords[, 1],
  lat = background_coords[, 2],
  source = "background"
)

background_sp <- SpatialPoints(
  background_coords,
  proj4string = CRS(proj4string(sampling_mask))
)

all_var_names <- c()
env_dir <- "earthenvstreams_china"

for (i in seq_len(nrow(vars_by_file))) {
  file_name <- vars_by_file$file[i]
  file_path <- file.path(env_dir, file_name)
  bands <- vars_by_file$bands[[i]]
  var_names <- vars_by_file$var_names[[i]]

  cat("  [", i, "/", nrow(vars_by_file), "] 处理: ", file_name, " (", length(bands), " 个变量)...\n", sep = "")

  tryCatch(
    {
      r <- brick(file_path)
      r_selected <- r[[bands]]
      values <- raster::extract(r_selected, background_sp)
      if (is.null(dim(values))) values <- matrix(values, ncol = 1)

      values[values == -127] <- NA
      values[values == -999] <- NA
      values[values == -9999] <- NA
      values[values < -1000] <- NA

      # 温度 / 10
      temp_vars <- grepl("^(tmin_|tmax_|hydro_wavg_0[1-9]|hydro_wavg_1[01])", var_names)
      if (any(temp_vars)) values[, which(temp_vars)] <- values[, which(temp_vars)] / 10

      # 坡度 / 100
      slope_vars <- grepl("^slope_", var_names)
      if (any(slope_vars)) values[, which(slope_vars)] <- values[, which(slope_vars)] / 100

      # pH / 10
      ph_vars <- var_names %in% c("soil_wavg_02")
      if (any(ph_vars)) values[, which(ph_vars)] <- values[, which(ph_vars)] / 10

      colnames(values) <- var_names
      background_df <- cbind(background_df, values)
      all_var_names <- c(all_var_names, var_names)
      
      cat("    ✓ 成功提取 ", length(var_names), " 个变量\n", sep = "")
      rm(r, r_selected, values); gc(verbose = FALSE)
    },
    error = function(e) cat("    ✗ 错误: ", e$message, "\n", sep = "")
  )
}

# 缺失处理：剔除缺失 > 10% 的点，其余中位数填充
n_generated <- nrow(background_df)
bg_missing_pct <- rowSums(is.na(background_df[, all_var_names])) / length(all_var_names) * 100
background_df <- background_df[bg_missing_pct < 10, ]

vars_to_impute <- all_var_names[colSums(is.na(background_df[, all_var_names])) > 0]
if (length(vars_to_impute) > 0) {
  for (vv in vars_to_impute) {
    background_df[[vv]][is.na(background_df[[vv]])] <- median(background_df[[vv]], na.rm = TRUE)
  }
}

n_final <- nrow(background_df)
cat("  - 最终保留背景点: ", n_final, " (保留率: ", round(n_final / n_generated * 100, 1), "%)\n", sep = "")

# ------------------------------------------------------------------------------
# 6. 合并与保存
# ------------------------------------------------------------------------------
cat("\n步骤 6/6: 合并与保存...\n")

presence_data$presence <- 1
background_df$presence <- 0

# 变量对齐
common_vars <- intersect(names(presence_data), names(background_df))
cols_keep <- c("id", "species", "lon", "lat", "source", setdiff(common_vars, c("id", "species", "lon", "lat", "source", "presence")), "presence")

combined_data <- rbind(presence_data[, cols_keep], background_df[, cols_keep])

write.csv(background_df, "output/03_background_points/background_points.csv", row.names = FALSE)
cat("  ✓ output/03_background_points/background_points.csv\n")

write.csv(combined_data, "output/03_background_points/combined_presence_absence.csv", row.names = FALSE)
cat("  ✓ output/03_background_points/combined_presence_absence.csv\n")

# 诊断图
try({
  bg_sf <- sf::st_as_sf(background_df[, c("lon", "lat")], coords = c("lon", "lat"), crs = 4326)
  bg_xy <- as.data.frame(sf::st_coordinates(sf::st_transform(bg_sf, crs = LAEA_CRS)))
  colnames(bg_xy) <- c("x", "y")
  
  p_hex <- ggplot(bg_xy, aes(x = x, y = y)) +
    stat_bin_2d(bins = 80) +
    scale_fill_viridis_c(option = "C", direction = -1) +
    coord_equal() +
    labs(title = "Background Points Density (LAEA)", x = "X (m)", y = "Y (m)", fill = "Count") +
    theme_minimal(base_family = "Arial") +
    theme(panel.grid = element_blank())
  
  ggsave("figures/03_background_points/bg_uniformity_laea_hex.png", p_hex, width = 4, height = 3.5, dpi = 600, bg = "white")
  cat("  ✓ 诊断图已保存\n")
}, silent = TRUE)

cat("\n脚本执行完成!\n\n")
