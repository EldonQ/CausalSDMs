#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 11d_cate_maps.R
# 功能说明: 使用已训练的因果森林模型 (grf) 基于真实环境栅格生成 CATE 空间图
# 输入文件: output/14_causal/causal_forest_model.rds
#          earthenvstreams_china/*.tif, selected_variables.csv, extracted_variables.csv
# 输出文件: output/11_prediction_maps/rasters/cate_map.tif
#          figures/11_prediction_maps/cate_map.png
#          output/11_prediction_maps/cate_summary.csv
# 作者: Nature级别科研项目
# 日期: 2025-10-24
# ==============================================================================

# 初始化
rm(list = ls())
gc()
setwd("E:/SDM01")

packages <- c("tidyverse", "raster", "grf", "sf", "rnaturalearth", "viridis", "sysfonts", "showtext", "terra", "svglite", "ggplot2")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

dir.create("output/11_prediction_maps/rasters", showWarnings = FALSE, recursive = TRUE)

ttry <- try({
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

cat("\n======================================\n")
cat("CATE 空间图 (真实数据，河网掩膜)\n")
cat("======================================\n\n")

if(!file.exists("output/14_causal/causal_forest_model.rds")) {
  stop("缺少因果森林模型: output/14_causal/causal_forest_model.rds；请先运行 14b_causal_effects.R")
}

mobj <- readRDS("output/14_causal/causal_forest_model.rds")
cf <- mobj$model
feat <- mobj$features

# 环境栅格（确保仅使用训练特征）
sel_vars <- read.csv("output/04_collinearity/selected_variables.csv", stringsAsFactors = FALSE)$variable
var_map <- read.csv("output/02_env_extraction/extracted_variables.csv", stringsAsFactors = FALSE)
var_map <- var_map[var_map$variable %in% sel_vars & var_map$variable %in% feat, c("variable", "file", "band")]

build_env_stack <- function(var_map_df, base_dir = "earthenvstreams_china") {
  groups <- split(var_map_df, var_map_df$file)
  stk_list <- list()
  for(fn in names(groups)) {
    g <- groups[[fn]]
    r <- raster::brick(file.path(base_dir, fn))
    r_sel <- r[[g$band]]
    names(r_sel) <- g$variable
    stk_list[[length(stk_list) + 1]] <- r_sel
    rm(r, r_sel)
    gc(verbose = FALSE)
  }
  stk <- raster::stack(stk_list)
  # 中文注释：仅保留模型特征在栅格中存在的部分，并按训练顺序排列；避免无效名称触发警告
  present <- intersect(feat, names(stk))
  if(length(present) < length(feat)) {
    msg_miss <- setdiff(feat, present)
    if(length(msg_miss) > 0) {
      cat("  ! 下列特征在环境栅格中缺失，将被跳过: ", paste(msg_miss, collapse = ", "), "\n", sep = "")
    }
  }
  stk <- stk[[present]]
  return(stk)
}

env_stack <- build_env_stack(var_map)

# ===== 单位对齐（与 02/11 保持一致）=====
# 温度（hydro_wavg_01–11）：℃×10 → ℃
idx_temp <- which(grepl("^hydro_wavg_0[1-9]$|^hydro_wavg_1[01]$", names(env_stack)))
if(length(idx_temp) > 0) {
  for(i in idx_temp) env_stack[[i]] <- env_stack[[i]] / 10
  cat("  - 单位转换: 温度 hydrowavg_01–11 ÷10\n")
}
# 坡度：度×100 → 度
idx_slope <- which(grepl("^slope_", names(env_stack)))
if(length(idx_slope) > 0) {
  for(i in idx_slope) env_stack[[i]] <- env_stack[[i]] / 100
  cat("  - 单位转换: 坡度 ÷100\n")
}
# 土壤 pH：pH×10 → pH
idx_ph <- which(names(env_stack) == "soil_wavg_02")
if(length(idx_ph) > 0) {
  env_stack[[idx_ph]] <- env_stack[[idx_ph]] / 10
  cat("  - 单位转换: 土壤pH ÷10\n")
}

# 河网掩膜
fa <- raster::brick("earthenvstreams_china/flow_acc.tif")[[2]]
fa_vals <- raster::getValues(fa) ; fa_vals[fa_vals <= 0] <- NA
river_mask <- raster::setValues(fa, fa_vals)
rm(fa_vals)

# 中国边界
china <- rnaturalearth::ne_countries(country = "China", scale = "medium", returnclass = "sf")

# 分块预测
bs <- raster::blockSize(env_stack)
out_path <- "output/11_prediction_maps/rasters/cate_map.tif"
if(file.exists(out_path)) { try({ file.remove(out_path) }, silent = TRUE) }
out_r <- raster::raster(env_stack, layer = 1)
out_r <- raster::setValues(out_r, NA_real_)
wr <- raster::writeStart(out_r, filename = out_path, overwrite = TRUE)

for(i in seq_len(bs$n)) {
  X_block <- raster::getValues(env_stack, row = bs$row[i], nrows = bs$nrows[i])
  if(NROW(X_block) == 0) {
    wr <- raster::writeValues(wr, rep(NA_real_, 0), bs$row[i])
    next
  }
  X_df <- as.data.frame(X_block)
  # 预测阶段沿用训练阶段的缺失处理：NA → 0（见 14b 脚本）
  X_df[!is.finite(as.matrix(X_df))] <- 0
  pred <- try({ predict(cf, as.matrix(X_df))$predictions }, silent = TRUE)
  vec <- if(!inherits(pred, "try-error")) as.numeric(pred) else rep(NA_real_, nrow(X_df))
  wr <- raster::writeValues(wr, vec, bs$row[i])
}
res_r <- raster::writeStop(wr)

# 掩膜至河网
res_riv <- raster::mask(res_r, river_mask)
raster::writeRaster(res_riv, out_path, overwrite = TRUE)

# 统计与制图
vals <- raster::getValues(res_riv)
vals <- vals[is.finite(vals)]
if(length(vals) == 0) {
  summary_df <- data.frame(n_pixels_river = 0, mean = NA, sd = NA, min = NA, max = NA, p10 = NA, p50 = NA, p90 = NA)
} else {
  summary_df <- data.frame(
    n_pixels_river = length(vals),
    mean = mean(vals), sd = sd(vals),
    min = min(vals), max = max(vals),
    p10 = as.numeric(stats::quantile(vals, 0.1)),
    p50 = as.numeric(stats::quantile(vals, 0.5)),
    p90 = as.numeric(stats::quantile(vals, 0.9))
  )
}
write.csv(summary_df, "output/11_prediction_maps/cate_summary.csv", row.names = FALSE)

viz_save_raster_map(r = res_riv, out_base = "figures/11_prediction_maps/cate_map",
                    title = "CATE (Causal Forest)",
                    palette = "inferno", q_limits = c(0.01, 0.99),
                    china_path = "earthenvstreams_china/china_boundary.shp",
                    width_in = 8, height_in = 6)

 

cat("\n完成：CATE 空间图与统计输出\n\n")
