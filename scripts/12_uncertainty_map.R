#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 12_uncertainty_map.R
# 功能说明: 分析4个模型预测的不确定性（标准差）
# 方法: 计算4个模型预测的标准差作为不确定性指标
# 输入文件: output/11_prediction_maps/all_predictions.csv
# 输出文件: figures/12_uncertainty/uncertainty_map.png
#          figures/12_uncertainty/model_agreement.png
#          output/12_uncertainty/uncertainty_summary.csv
# 作者: Nature级别科研项目
# 日期: 2025-10-20
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/SDM01")

# 加载必要的包（基于预测栅格计算并绘制河网热图）
packages <- c("tidyverse", "sf", "viridis", "rnaturalearth", "raster", "sysfonts", "showtext", "terra", "svglite", "ggplot2")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 字体设置：注册 Arial 并启用 showtext，确保 PNG/PDF 使用 Arial（期刊要求）
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

# 统一绘图工具（Nature风格/PNG+SVG/Arial）
source("scripts/visualization/viz_utils.R")

dir.create("output/12_uncertainty", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/12_uncertainty", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("模型不确定性分析（基于河网栅格）\n")
cat("======================================\n\n")

## 新流程：读取四个模型的河网概率栅格，逐像元计算均值/标准差/一致性指标
cat("步骤 1/3: 读取河网预测栅格...\n")
ras_paths <- list(
  Maxnet = "output/11_prediction_maps/rasters/pred_maxnet_river.tif",
  NN = "output/11_prediction_maps/rasters/pred_nn_river.tif",
  RF = "output/11_prediction_maps/rasters/pred_rf_river.tif",
  GAM = "output/11_prediction_maps/rasters/pred_gam_river.tif"
)

missing <- names(ras_paths)[!file.exists(unlist(ras_paths))]
if(length(missing) > 0) {
  stop(paste0("缺少预测栅格: ", paste(missing, collapse = ", "),
              "。请先运行 11_current_prediction_maps.R 生成河网预测。"))
}

stk <- raster::stack(unlist(ras_paths))
names(stk) <- names(ras_paths)

cat("  ✓ 输入层数: ", raster::nlayers(stk), "\n", sep = "")

cat("\n步骤 2/3: 计算不确定性指标 (逐像元)...\n")
mean_r <- raster::calc(stk, fun = function(x){ if(all(is.na(x))) NA else mean(x, na.rm = TRUE) })
sd_r   <- raster::calc(stk, fun = function(x){ if(all(is.na(x))) NA else sd(x, na.rm = TRUE) })
min_r  <- raster::calc(stk, fun = function(x){ if(all(is.na(x))) NA else min(x, na.rm = TRUE) })
max_r  <- raster::calc(stk, fun = function(x){ if(all(is.na(x))) NA else max(x, na.rm = TRUE) })
range_r <- max_r - min_r
agreement_r <- 1 - range_r  # 一致性：范围越小，一致性越高

# 保存GeoTIFF
raster::writeRaster(mean_r, "output/12_uncertainty/mean_prediction_river.tif", overwrite = TRUE)
raster::writeRaster(sd_r,   "output/12_uncertainty/sd_prediction_river.tif", overwrite = TRUE)
raster::writeRaster(agreement_r, "output/12_uncertainty/agreement_river.tif", overwrite = TRUE)

cat("  ✓ 指标栅格已保存\n")

cat("\n步骤 3/3: 绘制河网热图 (1200dpi, Arial)...\n")
china <- ne_countries(country = "China", scale = "medium", returnclass = "sf")

# 不确定性（标准差）图 —— 使用0~99%分位范围，提高对比度
sd_vals <- raster::getValues(sd_r)
sd_p99 <- as.numeric(stats::quantile(sd_vals, 0.99, na.rm = TRUE))
viz_save_raster_map(r = sd_r, out_base = "figures/12_uncertainty/uncertainty_map",
                    title = "Prediction Uncertainty (SD)",
                    palette = "magma", q_limits = c(0.01, 0.99),
                    china_path = "earthenvstreams_china/china_boundary.shp",
                    width_in = 8, height_in = 6)

 
cat("  ✓ 不确定性热图保存\n")

# 一致性图
viz_save_raster_map(r = agreement_r, out_base = "figures/12_uncertainty/model_agreement",
                    title = "Model Agreement",
                    palette = "viridis", q_limits = c(0.01, 0.99),
                    china_path = "earthenvstreams_china/china_boundary.shp",
                    width_in = 8, height_in = 6)

 
cat("  ✓ 一致性热图保存\n")

# 导出统计摘要CSV（仅河网像元）
vals_sd <- raster::getValues(sd_r); vals_sd <- vals_sd[!is.na(vals_sd)]
vals_ag <- raster::getValues(agreement_r); vals_ag <- vals_ag[!is.na(vals_ag)]
uncertainty_summary <- data.frame(
  metric = c("sd_prediction", "agreement"),
  mean = c(mean(vals_sd), mean(vals_ag)),
  sd = c(sd(vals_sd), sd(vals_ag)),
  min = c(min(vals_sd), min(vals_ag)),
  max = c(max(vals_sd), max(vals_ag)),
  p10 = c(as.numeric(quantile(vals_sd, 0.1)), as.numeric(quantile(vals_ag, 0.1))),
  p50 = c(as.numeric(quantile(vals_sd, 0.5)), as.numeric(quantile(vals_ag, 0.5))),
  p90 = c(as.numeric(quantile(vals_sd, 0.9)), as.numeric(quantile(vals_ag, 0.9)))
)

write.csv(uncertainty_summary, "output/12_uncertainty/uncertainty_summary.csv", row.names = FALSE)

# 日志
sink("output/12_uncertainty/processing_log.txt")
cat("不确定性分析日志\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
cat("输入预测栅格: 4 (Maxnet, NN, RF, GAM)\n\n")
cat("统计摘要:\n")
print(uncertainty_summary)
sink()

cat("\n======================================\n")
cat("不确定性分析完成\n")
cat("======================================\n\n")

cat("不确定性统计 (河网像元摘要) 已保存至 output/12_uncertainty/uncertainty_summary.csv\n")

cat("\n✓ 脚本执行完成!\n\n")
