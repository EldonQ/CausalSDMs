#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 15_future_env_projection.R
# 功能说明: 处理和可视化未来气候情景下的生物气候变量
# 方法: 裁剪到中国区域，生成4个SSP情景的热图和趋势图
# 输入文件: E:/WorldClim/Future/*/wc2.1_30s_bioc_*.tif
# 输出文件: output/15_future_env/future_bioc_*_ssp*.tif
#          figures/15_future_env/bioc_*_comparison.png (空间分布)
#          figures/15_future_env/climate_change_trends.png (时间趋势)
# 作者: Nature级别科研项目
# 日期: 2025-10-20
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/SDM01")

# 加载必要的包
packages <- c("raster", "sf", "tidyverse", "viridis", "patchwork", "sysfonts", "showtext", "maxnet", "randomForest", "mgcv", "nnet", "terra", "svglite", "ggplot2")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

dir.create("output/15_future_env", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/15_future_env", showWarnings = FALSE, recursive = TRUE)

# 字体（确保后续制图使用Arial）
try({
  sysfonts::font_add(
    family = "Arial",
    regular = "C:/Windows/Fonts/arial.ttf",
    bold = "C:/Windows/Fonts/arialbd.ttf",
    italic = "C:/Windows/Fonts/ariali.ttf",
    bolditalic = "C:/Windows/Fonts/arialbi.ttf"
  )
  showtext::showtext_opts(dpi = 1200)
  showtext::showtext_auto(enable = TRUE)
}, silent = TRUE)

# 统一绘图工具
source("scripts/visualization/viz_utils.R")

cat("\n======================================\n")
cat("未来气候情景下的河网适生度预测（重训模型）\n")
cat("======================================\n\n")

# 可视化开关（附录性质）
# PLOT_CLIMATE_MAPS: 是否输出未来气候变量背景图（与预测无直接耦合）
#   - FALSE: 不输出
#   - TRUE: 输出，模式由 CLIMATE_MAP_MODE 决定
# CLIMATE_MAP_MODE: "all" 输出 bio01-19 全部；"auto_topK" 输出跨SSP变化最大的前K个
# CLIMATE_MAP_TOPK: auto_topK 模式下的 K 值
# PLOT_TRENDS: 是否输出趋势柱状图
PLOT_CLIMATE_MAPS <- FALSE
CLIMATE_MAP_MODE  <- "auto_topK"   # 或 "all"
CLIMATE_MAP_TOPK  <- 6
PLOT_TRENDS       <- FALSE

# 1. 读取中国边界
cat("步骤 1/5: 读取中国边界...\n")
china <- st_read("data-main/vector/china.shp", quiet = TRUE)
bbox <- st_bbox(china)
bbox_ext <- bbox + c(-2, -2, 2, 2)
cat("  ✓ 中国边界已加载\n")

# 2. 定义SSP情景
cat("\n步骤 2/5: 定义未来气候情景...\n")
ssp_scenarios <- c("SSP126", "SSP245", "SSP370", "SSP585")
cat("  - 4个SSP情景: SSP1-2.6 (低), SSP2-4.5 (中), SSP3-7.0 (高), SSP5-8.5 (极高)\n")
cat("  - GCM: BCC-CSM2-MR, 时间: 2041-2060\n")
cat("  - 生物气候变量: 19个 (bio01-bio19)\n")

# bioc变量信息
bioc_info <- data.frame(
  code = sprintf("bio%02d", 1:19),
  name = c("Ann Mean Temp", "Mean Diurnal Range", "Isothermality",
           "Temp Seasonality", "Max Temp Warmest", "Min Temp Coldest",
           "Temp Ann Range", "Mean Temp Wettest", "Mean Temp Driest",
           "Mean Temp Warmest", "Mean Temp Coldest", "Ann Precip",
           "Precip Wettest", "Precip Driest", "Precip Seasonality",
           "Precip Wettest Q", "Precip Driest Q", "Precip Warmest Q",
           "Precip Coldest Q"),
  stringsAsFactors = FALSE
)

# 3. 裁剪未来bioc数据
cat("\n步骤 3/5: 裁剪未来气候数据到中国区域...\n")
future_bioc_list <- list()
all_summaries <- list()  # 用于汇总四情景的河网均值趋势

for(ssp in ssp_scenarios) {
  cat("  处理情景: ", ssp, "\n", sep = "")
  
  bioc_file <- file.path("E:/WorldClim/Future", ssp,
                        paste0("wc2.1_30s_bioc_BCC-CSM2-MR_",
                               tolower(ssp), "_2041-2060.tif"))
  
  if(!file.exists(bioc_file)) {
    cat("    ✗ 文件不存在\n")
    next
  }
  
  # 读取并裁剪
  bioc_raster <- brick(bioc_file)
  bioc_china <- crop(bioc_raster, extent(bbox_ext))
  bioc_china <- mask(bioc_china, china)
  
  # 保存
  output_file <- paste0("output/15_future_env/future_bioc_china_", ssp, ".tif")
  writeRaster(bioc_china, output_file, format = "GTiff", overwrite = TRUE)
  
  future_bioc_list[[ssp]] <- bioc_china
  cat("    ✓ 完成 (", raster::nlayers(bioc_china), " 波段)\n", sep = "")
  
  rm(bioc_raster)
  gc(verbose = FALSE)
}

cat("  ✓ 共处理 ", length(future_bioc_list), " 个SSP情景\n\n", sep = "")

# 4. 计算统计摘要
cat("步骤 4/5: 计算统计摘要...\n")
stats_list <- list()
for(ssp in names(future_bioc_list)) {
  bioc_stack <- future_bioc_list[[ssp]]
  for(i in 1:nlayers(bioc_stack)) {
    bio_values <- getValues(bioc_stack[[i]])
    bio_values <- bio_values[!is.na(bio_values)]
    stats_list[[length(stats_list) + 1]] <- data.frame(
      scenario = ssp,
      bioc = sprintf("bio%02d", i),
      mean = mean(bio_values),
      sd = sd(bio_values),
      min = min(bio_values),
      max = max(bio_values)
    )
  }
}

stats_df <- do.call(rbind, stats_list)
write.csv(stats_df, "output/15_future_env/future_bioc_statistics.csv", row.names = FALSE)
cat("  ✓ 统计摘要已保存\n")

if (PLOT_CLIMATE_MAPS || PLOT_TRENDS) {
  cat("\n步骤 5/5: 生成气候背景图（可选）...\n")
}

if (PLOT_CLIMATE_MAPS) {
  # 自动选择或全部
  if (CLIMATE_MAP_MODE == "all") {
    key_biocs <- seq_len(19)
  } else {
    # 基于跨SSP的场均均值方差选择变化最大的前K个
    var_rank <- stats_df %>%
      dplyr::group_by(bioc) %>%
      dplyr::summarise(var_ssp = stats::var(mean), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(var_ssp))
    sel_codes <- head(var_rank$bioc, CLIMATE_MAP_TOPK)
    key_biocs <- as.integer(sub("bio", "", sel_codes))
  }

for(bio_idx in key_biocs) {
  bio_code <- sprintf("bio%02d", bio_idx)
  bio_name <- bioc_info$name[bio_idx]
    cat("  - ", bio_code, ": ", bio_name, "\n", sep = "")
  
  png(paste0("figures/15_future_env/", bio_code, "_comparison.png"),
        width = 4800, height = 3600, res = 1200, type = "cairo-png", family = "Arial")
  par(mfrow = c(2, 2), mar = c(2, 2, 3, 3), oma = c(0, 0, 2, 0))
  for(ssp in names(future_bioc_list)) {
    bio_layer <- future_bioc_list[[ssp]][[bio_idx]]
    plot(bio_layer, main = paste0(ssp, " (2041-2060)"),
           col = viridis::viridis(100, option = "D"), cex.main = 0.9,
         cex.axis = 0.6, legend = TRUE, axes = TRUE)
      plot(sf::st_geometry(china), add = TRUE, border = "grey30", lwd = 0.5)
  }
  mtext(paste0(bio_code, ": ", bio_name), outer = TRUE, cex = 1.1, font = 2)
  dev.off()
  
  }
}

if (PLOT_TRENDS) {
  cat("  (趋势图) 生成气候变化趋势图...\n")
climate_data <- stats_df %>%
    dplyr::filter(bioc %in% c("bio01", "bio05", "bio06", "bio12")) %>%
    dplyr::select(scenario, bioc, mean) %>%
    tidyr::pivot_wider(names_from = bioc, values_from = mean) %>%
    dplyr::mutate(
    scenario_name = factor(scenario, 
                          levels = c("SSP126", "SSP245", "SSP370", "SSP585"),
                          labels = c("SSP1-2.6\n(Low)", "SSP2-4.5\n(Medium)", 
                                   "SSP3-7.0\n(High)", "SSP5-8.5\n(Very High)")),
      temp_annual = bio01 / 10,
      temp_max = bio05 / 10,
      temp_min = bio06 / 10,
      precip_annual = bio12
    )

p_temp <- ggplot(climate_data, aes(x = scenario_name, y = temp_annual)) +
  geom_bar(stat = "identity", fill = "#D73027", alpha = 0.8) +
    labs(title = "Annual Mean Temperature (2041-2060)", x = "SSP Scenario", y = expression("Temperature ("*degree*"C)")) +
    theme_minimal(base_size = 8)

p_precip <- ggplot(climate_data, aes(x = scenario_name, y = precip_annual)) +
  geom_bar(stat = "identity", fill = "#2166AC", alpha = 0.8) +
    labs(title = "Annual Precipitation (2041-2060)", x = "SSP Scenario", y = "Precipitation (mm)") +
    theme_minimal(base_size = 8)

climate_range <- climate_data %>%
    dplyr::select(scenario_name, temp_max, temp_min) %>%
    tidyr::pivot_longer(cols = c(temp_max, temp_min), names_to = "type", values_to = "temp")
p_range <- ggplot(climate_range, aes(x = scenario_name, y = temp, fill = type)) +
  geom_bar(stat = "identity", position = "dodge", alpha = 0.8) +
    scale_fill_manual(values = c("temp_max" = "#D73027", "temp_min" = "#4575B4"), labels = c("Max (Warmest)", "Min (Coldest)"), name = "Temperature") +
    labs(title = "Temperature Range (2041-2060)", x = "SSP Scenario", y = expression("Temperature ("*degree*"C)")) +
    theme_minimal(base_size = 8)

  combined_plot <- (p_temp | p_precip) / p_range
  ggsave("figures/15_future_env/climate_change_trends.png", plot = combined_plot, width = 4, height = 3, dpi = 1200)
  cat("      ✓ 气候变化趋势图已保存\n")
}

# 6. 基于重训模型的未来河网适生度预测（每情景 × 4模型）
cat("\n步骤 6/6: 使用重训模型进行河网适生度预测...\n")

# 检查模型
models_path <- "output/15_future_env/models"
if(!dir.exists(models_path)) {
  stop("未找到重训模型，请先运行 scripts/15a_retrain_future_vars.R")
}
mdl_files <- c(Maxnet = file.path(models_path, "maxnet.rds"),
               NN = file.path(models_path, "nn.rds"),
               RF = file.path(models_path, "rf.rds"),
               GAM = file.path(models_path, "gam.rds"))
miss_m <- names(mdl_files)[!file.exists(unlist(mdl_files))]
if(length(miss_m) > 0) stop(paste0("缺少模型: ", paste(miss_m, collapse=", ")))

# 构建未来环境栅格（与重训变量一致）
build_future_env <- function(bioc_stack) {
  names(bioc_stack) <- sprintf("bio%02d", 1:raster::nlayers(bioc_stack))
  temp_ann_mean <- bioc_stack[["bio01"]] / 10; names(temp_ann_mean) <- "temp_ann_mean"
  prec_ann_sum  <- bioc_stack[["bio12"]];      names(prec_ann_sum)  <- "prec_ann_sum"
  static_list <- list()
  base_ref <- temp_ann_mean
  # dem_avg（band4）
  if(file.exists("earthenvstreams_china/elevation.tif")) {
    dem_br <- raster::brick("earthenvstreams_china/elevation.tif")[[4]]
    dem_br <- suppressWarnings(raster::projectRaster(dem_br, base_ref, method = "bilinear"))
    names(dem_br) <- "dem_avg"
    static_list[[length(static_list)+1]] <- dem_br
  }
  # slope_avg（band4，度×100 → 度）
  if(file.exists("earthenvstreams_china/slope.tif")) {
    slope_br <- raster::brick("earthenvstreams_china/slope.tif")[[4]] / 100
    slope_br <- suppressWarnings(raster::projectRaster(slope_br, base_ref, method = "bilinear"))
    names(slope_br) <- "slope_avg"
    static_list[[length(static_list)+1]] <- slope_br
  }
  # geology_total：所有band求和
  if(file.exists("earthenvstreams_china/geology_weighted_sum.tif")) {
    geo_total <- raster::calc(raster::brick("earthenvstreams_china/geology_weighted_sum.tif"), fun = function(x){
      if(all(is.na(x))) NA else sum(x, na.rm = TRUE)
    })
    geo_total <- suppressWarnings(raster::projectRaster(geo_total, base_ref, method = "bilinear"))
    names(geo_total) <- "geology_total"
    static_list[[length(static_list)+1]] <- geo_total
  }
  # soil_avg_01（SOC）
  if(file.exists("earthenvstreams_china/soil_average.tif")) {
    soil_soc <- raster::brick("earthenvstreams_china/soil_average.tif")[[1]]
    soil_soc <- suppressWarnings(raster::projectRaster(soil_soc, base_ref, method = "bilinear"))
    names(soil_soc) <- "soil_avg_01"
    static_list[[length(static_list)+1]] <- soil_soc
  }
  stk <- raster::stack(c(list(temp_ann_mean, prec_ann_sum), static_list))
  lon_r <- raster::init(stk[[1]], fun = 'x'); names(lon_r) <- 'lon'
  lat_r <- raster::init(stk[[1]], fun = 'y'); names(lat_r) <- 'lat'
  stk <- raster::addLayer(stk, lon_r, lat_r)
  return(stk)
}

# 河网掩膜
fa <- raster::brick("earthenvstreams_china/flow_acc.tif")[[2]]
fa_vals <- raster::getValues(fa); fa_vals[fa_vals <= 0] <- NA
river_mask <- raster::setValues(fa, fa_vals)
rm(fa_vals)

make_predict_fun <- function(model_name, model_obj) {
  if(model_name == "Maxnet") return(function(m, df) { as.numeric(predict(m, df, type = "logistic")) })
  if(model_name == "RF")     return(function(m, df) { as.numeric(predict(m, newdata = df, type = "prob")[, "1"]) })
  if(model_name == "GAM")    return(function(m, df) { as.numeric(predict(m, newdata = df, type = "response")) })
  if(model_name == "NN")     return(function(m, df) {
    mu <- m$mean; sdv <- m$sd; mod <- m$model; vars <- m$vars
    sdv[sdv == 0 | is.na(sdv)] <- 1
    x <- as.matrix(df[, vars, drop = FALSE])
    x <- sweep(x, 2, mu[vars], "-")
    x <- sweep(x, 2, sdv[vars], "/")
    as.numeric(nnet:::predict.nnet(mod, x, type = "raw"))
  })
}

for(ssp in names(future_bioc_list)) {
  bioc_china <- future_bioc_list[[ssp]]
  env_stk <- build_future_env(bioc_china)
  out_dir_ras <- file.path("output/15_future_env/rasters", ssp)
  out_dir_fig <- file.path("figures/15_future_env", ssp)
  dir.create(out_dir_ras, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)

  summary_rows <- list()
  for(mn in names(mdl_files)) {
    cat("  -> 情景 ", ssp, " | 模型 ", mn, " ...\n", sep = "")
    mdl <- readRDS(mdl_files[[mn]])
    pred_fun <- make_predict_fun(mn, mdl)
    tif_path <- file.path(out_dir_ras, paste0("pred_", tolower(mn), ".tif"))
    if(file.exists(tif_path)) { try({ file.remove(tif_path) }, silent = TRUE) }
    pred_r <- raster::predict(env_stk, model = mdl, fun = pred_fun, filename = tif_path, overwrite = TRUE, progress = "text")
    pred_r <- clamp(pred_r, lower = 0, upper = 1, useValues = TRUE)
    # 对齐河网掩膜到预测栅格，避免 extent/分辨率不一致
    river_mask_ref <- suppressWarnings(raster::projectRaster(river_mask, pred_r, method = "ngb"))
    pred_r_river <- raster::mask(pred_r, river_mask_ref)
    tif_mask_path <- file.path(out_dir_ras, paste0("pred_", tolower(mn), "_river.tif"))
    if(file.exists(tif_mask_path)) { try({ file.remove(tif_mask_path) }, silent = TRUE) }
    raster::writeRaster(pred_r_river, tif_mask_path, overwrite = TRUE)

    vals <- raster::getValues(pred_r_river); vals <- vals[!is.na(vals)]
    summary_rows[[length(summary_rows)+1]] <- data.frame(
      scenario = ssp, model = mn, n_pixels_river = length(vals),
      mean = mean(vals), sd = sd(vals), min = min(vals), max = max(vals),
      p10 = as.numeric(quantile(vals, 0.1)), p50 = as.numeric(quantile(vals, 0.5)), p90 = as.numeric(quantile(vals, 0.9))
    )

    # 图件
    out_base <- file.path(out_dir_fig, paste0("prediction_", tolower(mn)))
    viz_save_raster_map(r = pred_r_river, out_base = out_base,
                        title = paste0(ssp, " - ", mn),
                        palette = "magma", q_limits = c(0.01, 0.99),
                        china_path = "earthenvstreams_china/china_boundary.shp",
                        width_in = 8, height_in = 6)

    
  }
  if(length(summary_rows) > 0) {
    summary_df <- dplyr::bind_rows(summary_rows)
    # 保存本情景的统计摘要
    write.csv(summary_df, file.path(out_dir_ras, "prediction_summary.csv"), row.names = FALSE)
    # 记录到全局列表，供趋势图使用
    summary_df$scenario <- factor(summary_df$scenario, levels = c("SSP126","SSP245","SSP370","SSP585"))
    all_summaries[[ssp]] <- summary_df
  }
}

# 日志
sink("output/15_future_env/processing_log.txt")
cat("未来气候情景处理日志\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
cat("SSP情景: ", paste(ssp_scenarios, collapse = ", "), "\n", sep = "")
cat("GCM: BCC-CSM2-MR\n")
cat("时间段: 2041-2060\n")
cat("生物气候变量: 19个\n\n")
cat("生成图表:\n")
cat("  - 空间分布地图: 6张 (关键变量 × 4情景)\n")
cat("  - 气候变化趋势: 1张\n\n")
cat("气候统计 (bio01 - 年平均温度):\n")
print(stats_df[stats_df$bioc == "bio01", ])
sink()

# 追加输出：未来四情景生境变化趋势（各模型河网均值）
if(length(all_summaries) > 0) {
  trends <- dplyr::bind_rows(all_summaries)
  out_trend_csv <- "output/15_future_env/prediction_trends_all_models.csv"
  out_trend_png <- "figures/15_future_env/habitat_trends_all_models.png"
  write.csv(trends, out_trend_csv, row.names = FALSE)
  p_trend <- ggplot(trends, aes(x = scenario, y = mean, group = model, color = model)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.4) +
    labs(title = "Habitat Suitability Trend across SSP Scenarios",
         x = "Scenario", y = "Mean Predicted Suitability", color = "Model") +
    viz_theme_nature(base_size = 8, title_size = 9)
  ggsave(out_trend_png, p_trend, width = 4.8, height = 3.2, dpi = 1200, bg = "white")
}

cat("\n======================================\n")
cat("未来气候情景处理完成\n")
cat("======================================\n")
cat("处理情景: ", length(future_bioc_list), "\n", sep = "")
cat("输出文件:\n")
cat("  - 栅格数据: output/15_future_env/future_bioc_china_*.tif\n")
cat("  - 空间地图: figures/15_future_env/bio*_comparison.png (6张)\n")
cat("  - 趋势图: figures/15_future_env/climate_change_trends.png\n")
cat("  - 统计表: output/15_future_env/future_bioc_statistics.csv\n\n")

cat("✓ 脚本执行完成!\n\n")
