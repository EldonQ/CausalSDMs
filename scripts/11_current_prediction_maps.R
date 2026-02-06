#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 11_current_prediction_maps.R
# 功能说明: 生成当前气候条件下模型(Maxnet, RF)及其集成(Ensemble)的物种分布预测地图
# 升级: 所有输出结果均保存至 MaxnetRF 子文件夹，保持目录整洁
# 输入文件: output/05_model_maxnet/model.rds
#          output/06_model_rf/model.rds
# 输出文件: figures/11_prediction_maps/MaxnetRF/prediction_*.png
#          output/11_prediction_maps/MaxnetRF/rasters/pred_*.tif
#          output/11_prediction_maps/MaxnetRF/prediction_summary.csv
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/CausalSDMs")

# 加载必要的包
packages <- c("tidyverse", "sf", "ggplot2", "viridis", "rnaturalearth", "raster", "maxnet", "randomForest", "sysfonts", "showtext", "terra", "svglite")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 字体注册
try(
  {
    sysfonts::font_add(
      family = "Arial",
      regular = "C:/Windows/Fonts/arial.ttf",
      bold = "C:/Windows/Fonts/arialbd.ttf",
      italic = "C:/Windows/Fonts/ariali.ttf",
      bolditalic = "C:/Windows/Fonts/arialbi.ttf"
    )
    showtext::showtext_opts(dpi = 2400)
    showtext::showtext_auto(enable = TRUE)
  },
  silent = TRUE
)

# 统一绘图工具
source("scripts/visualization/viz_utils.R")

# 创建统一的 MaxnetRF 子目录
out_dir <- "output/11_prediction_maps/MaxnetRF"
fig_dir <- "figures/11_prediction_maps/MaxnetRF"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(out_dir, "rasters"), showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("当前气候预测地图绘制 (MaxnetRF 子目录版)\n")
cat("======================================\n\n")

# 1. 构建与训练变量一致的环境栅格栈
cat("步骤 1/6: 构建环境栅格栈 (与训练变量一致)...\n")
sel_vars <- read.csv("output/04_collinearity/selected_variables.csv", stringsAsFactors = FALSE)$variable
cat("  ✓ 入模变量数: ", length(sel_vars), "\n", sep = "")

var_map <- read.csv("output/02_env_extraction/extracted_variables.csv", stringsAsFactors = FALSE)
var_map <- var_map[var_map$variable %in% sel_vars, c("variable", "file", "band")]

build_env_stack <- function(var_map_df, base_dir = "earthenvstreams_china") {
  groups <- split(var_map_df, var_map_df$file)
  stk_list <- list()
  for (fn in names(groups)) {
    g <- groups[[fn]]
    r <- raster::brick(file.path(base_dir, fn))
    r_sel <- r[[g$band]]
    names(r_sel) <- g$variable
    stk_list[[length(stk_list) + 1]] <- r_sel
    rm(r, r_sel)
    gc(verbose = FALSE)
  }
  stk <- raster::stack(stk_list)
  stk <- stk[[sel_vars]]
  return(stk)
}

env_stack <- build_env_stack(var_map)

idx_temp <- which(grepl("^hydro_wavg_0[1-9]$|^hydro_wavg_1[01]$", names(env_stack)))
if (length(idx_temp) > 0) {
  for (i in idx_temp) env_stack[[i]] <- env_stack[[i]] / 10
  cat("  - 单位转换: 温度 ÷10\n")
}
idx_slope <- which(grepl("^slope_", names(env_stack)))
if (length(idx_slope) > 0) {
  for (i in idx_slope) env_stack[[i]] <- env_stack[[i]] / 100
  cat("  - 单位转换: 坡度 ÷100\n")
}
idx_ph <- which(names(env_stack) == "soil_wavg_02")
if (length(idx_ph) > 0) {
  env_stack[[idx_ph]] <- env_stack[[idx_ph]] / 10
  cat("  - 单位转换: 土壤pH ÷10\n")
}

# 2. 构建河网掩膜
cat("\n步骤 2/6: 生成河网掩膜...\n")
fa <- raster::brick("earthenvstreams_china/flow_acc.tif")[[2]]
fa_vals <- raster::getValues(fa)
fa_vals[fa_vals <= 0] <- NA
fa <- raster::setValues(fa, fa_vals)
river_mask <- fa
rm(fa_vals, fa)
gc()
cat("  ✓ 河网掩膜已生成\n")

# 3. 单模型预测
cat("\n步骤 3/6: 单模型预测...\n")
models <- c("Maxnet", "RF")
model_files <- c(
  Maxnet = "output/05_model_maxnet/model.rds",
  RF = "output/06_model_rf/model.rds"
)

make_predict_fun <- function(model_name, model_obj) {
  if (model_name == "Maxnet") {
    return(function(m, df) as.numeric(predict(m, df, type = "logistic")))
  }
  if (model_name == "RF") {
    return(function(m, df) as.numeric(predict(m, newdata = df, type = "prob")[, "1"]))
  }
}

summary_rows <- list()
pred_stack_list <- list()

for (m in models) {
  cat("  -> ", m, " ...\n", sep = "")
  if (!file.exists(model_files[[m]])) next
  
  mdl <- readRDS(model_files[[m]])
  pred_fun <- make_predict_fun(m, mdl)
  
  # 保存到 MaxnetRF/rasters
  tif_path <- file.path(out_dir, "rasters", paste0("pred_", tolower(m), ".tif"))
  
  if (file.exists(tif_path)) unlink(tif_path)
  pred_r <- raster::predict(env_stack, model = mdl, fun = pred_fun, filename = tif_path, overwrite = TRUE, progress = "text")
  
  pred_r <- clamp(pred_r, lower = 0, upper = 1, useValues = TRUE)
  pred_r_river <- mask(pred_r, river_mask)
  
  tif_mask_path <- file.path(out_dir, "rasters", paste0("pred_", tolower(m), "_river.tif"))
  writeRaster(pred_r_river, tif_mask_path, overwrite = TRUE)
  
  pred_stack_list[[m]] <- pred_r_river
  
  vals <- getValues(pred_r_river)
  vals <- vals[!is.na(vals)]
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    model = m, n_pixels = length(vals), mean = mean(vals), sd = sd(vals), min = min(vals), max = max(vals)
  )
  
  # 保存到 MaxnetRF 文件夹
  out_base <- file.path(fig_dir, paste0("prediction_", tolower(m)))
  viz_save_raster_map(
    r = pred_r_river, out_base = out_base,
    title = paste(m, "Prediction"),
    palette = "magma", q_limits = c(0.01, 0.99),
    china_path = "earthenvstreams_china/china_boundary.shp",
    width_in = 6.5, height_in = 5.0
  )
}

# 4. 集成预测
cat("\n步骤 4/6: 集成预测 (Ensemble)...\n")

if (length(pred_stack_list) == 2) {
  stk <- stack(pred_stack_list[["Maxnet"]], pred_stack_list[["RF"]])
  ensemble_mean <- calc(stk, fun = mean, na.rm = TRUE)
  
  tif_ens_path <- file.path(out_dir, "rasters", "pred_ensemble_river.tif")
  writeRaster(ensemble_mean, tif_ens_path, overwrite = TRUE)
  
  vals <- getValues(ensemble_mean)
  vals <- vals[!is.na(vals)]
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    model = "Ensemble", n_pixels = length(vals), mean = mean(vals), sd = sd(vals), min = min(vals), max = max(vals)
  )
  cat("  ✓ 集成栅格已保存\n")
  
  out_base <- file.path(fig_dir, "prediction_ensemble")
  viz_save_raster_map(
    r = ensemble_mean, out_base = out_base,
    title = "Ensemble Prediction",
    palette = "spectral",
    q_limits = c(0, 1),
    china_path = "earthenvstreams_china/china_boundary.shp",
    width_in = 6.5, height_in = 5.0
  )
  cat("  ✓ 集成地图已绘制: prediction_ensemble.png\n")
} else {
  cat("  ⚠ 预测模型不足2个，无法进行集成\n")
}

# 5. 保存统计表
cat("\n步骤 5/6: 保存统计表...\n")
summary_data <- bind_rows(summary_rows)
write.csv(summary_data, file.path(out_dir, "prediction_summary.csv"), row.names = FALSE)

# 6. 计算差值图
cat("\n步骤 6/6: 绘制模型差异图 (RF - Maxnet)...\n")
if (length(pred_stack_list) == 2) {
  diff_r <- pred_stack_list[["RF"]] - pred_stack_list[["Maxnet"]]
  out_base <- file.path(fig_dir, "prediction_difference")
  viz_save_raster_map(
    r = diff_r, out_base = out_base,
    title = "Model Difference (RF - Maxnet)",
    palette = "vik",
    q_limits = c(0.02, 0.98), # 使用分位数区间 (2%-98%) 自动适应数据范围
    china_path = "earthenvstreams_china/china_boundary.shp",
    width_in = 6.5, height_in = 5.0
  )
  cat("  ✓ 差异图已绘制\n")
}

cat("\n✓ 脚本执行完成! 结果保存在 MaxnetRF 子目录中。\n")
