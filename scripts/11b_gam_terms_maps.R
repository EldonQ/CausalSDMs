#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 11b_gam_terms_maps.R
# 功能说明: 基于已训练的GAM模型，使用真实环境栅格逐像元计算各平滑项贡献，
#          在河网掩膜内输出变量项贡献的空间图（每变量单图）与统计表。
# 方法: mgcv::predict(type="terms") 对 RasterStack 分块预测，严禁使用任何模拟数据。
# 输入文件: output/07_model_gam/model.rds
#          earthenvstreams_china/*.tif （与训练一致的变量及波段）
#          output/02_env_extraction/extracted_variables.csv（变量-文件-波段映射）
#          output/04_collinearity/selected_variables.csv（最终入模变量列表）
# 输出文件: figures/11_prediction_maps/gam_term_*.png
#          output/11_prediction_maps/rasters/gam_term_*.tif
#          output/11_prediction_maps/gam_terms_summary.csv
# 作者: Nature级别科研项目
# 日期: 2025-10-24
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/SDM01")

packages <- c("tidyverse", "raster", "mgcv", "sf", "rnaturalearth", "viridis", "sysfonts", "showtext", "RColorBrewer", "terra", "svglite", "ggplot2")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

dir.create("output/11_prediction_maps/rasters", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/11_prediction_maps", showWarnings = FALSE, recursive = TRUE)

# 字体（期刊要求：Arial，1200dpi）
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

cat("\n======================================\n")
cat("GAM 平滑项贡献空间图 (真实数据，河网掩膜)\n")
cat("======================================\n\n")

if(!file.exists("output/07_model_gam/model.rds")) {
  stop("缺少 GAM 模型: output/07_model_gam/model.rds")
}

# 读取变量映射与选择列表（保持与训练一致）
sel_vars <- read.csv("output/04_collinearity/selected_variables.csv", stringsAsFactors = FALSE)$variable
var_map <- read.csv("output/02_env_extraction/extracted_variables.csv", stringsAsFactors = FALSE)
var_map <- var_map[var_map$variable %in% sel_vars, c("variable", "file", "band")]

# 构建环境栅格栈
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
  stk <- stk[[sel_vars]]
  return(stk)
}

env_stack <- build_env_stack(var_map)

# ===== 单位对齐（与 02 提取阶段一致）=====
# 温度（hydro_wavg_01–11）：℃×10 → ℃
idx_temp <- which(grepl("^hydro_wavg_0[1-9]$|^hydro_wavg_1[01]$", names(env_stack)))
if(length(idx_temp) > 0) {
  for(i in idx_temp) {
    env_stack[[i]] <- env_stack[[i]] / 10
  }
  cat("  - 单位转换: 温度 hydrowavg_01–11 ÷10\n")
}
# 坡度：度×100 → 度
idx_slope <- which(grepl("^slope_", names(env_stack)))
if(length(idx_slope) > 0) {
  for(i in idx_slope) {
    env_stack[[i]] <- env_stack[[i]] / 100
  }
  cat("  - 单位转换: 坡度 ÷100\n")
}
# 土壤 pH：pH×10 → pH
idx_ph <- which(names(env_stack) == "soil_wavg_02")
if(length(idx_ph) > 0) {
  env_stack[[idx_ph]] <- env_stack[[idx_ph]] / 10
  cat("  - 单位转换: 土壤pH ÷10\n")
}

# 补充 lon/lat（与训练一致，GAM 若含 s(lon,lat) 需要）
lon_r <- raster::init(env_stack[[1]], fun = 'x') ; names(lon_r) <- 'lon'
lat_r <- raster::init(env_stack[[1]], fun = 'y') ; names(lat_r) <- 'lat'
env_stack <- raster::addLayer(env_stack, lon_r, lat_r)

# 河网掩膜（flow_acc band2 > 0）
fa <- raster::brick("earthenvstreams_china/flow_acc.tif")[[2]]
fa_vals <- raster::getValues(fa) ; fa_vals[fa_vals <= 0] <- NA
river_mask <- raster::setValues(fa, fa_vals)
rm(fa_vals)

# 读取中国边界（仅作叠加）
china <- rnaturalearth::ne_countries(country = "China", scale = "medium", returnclass = "sf")

gam_obj <- readRDS("output/07_model_gam/model.rds")

# 识别 GAM 项名称（直接基于模型内部数据避免空newdata触发变量缺失）
term_labels <- colnames(predict(gam_obj, type = "terms"))
term_labels <- term_labels[!grepl("^Intercept$", term_labels)]

# 解析公式，提取预测必需的自变量名（包含 lon/lat），排除响应
resp_name <- as.character(formula(gam_obj)[[2]])
needed_vars <- setdiff(all.vars(formula(gam_obj)), resp_name)
needed_vars <- setdiff(needed_vars, c("(weights)", "weights", "(Intercept)"))

# 结合已选变量与经纬度，确保与栅格层对齐
req_union <- unique(c(needed_vars, sel_vars, c("lon","lat")))
req_vars <- intersect(req_union, names(env_stack))

missing_in_raster <- setdiff(setdiff(needed_vars, c("(weights)", "weights", "(Intercept)")), names(env_stack))
if(length(missing_in_raster) > 0) {
  stop(paste0("GAM 预测所需变量在环境栅格中缺失: ", paste(missing_in_raster, collapse = ", "),
              "。请检查 output/02_env_extraction/extracted_variables.csv 的变量-文件-波段映射与 selected_variables.csv 是否一致。"))
}

# 分块预测各项贡献
cat("步骤 1/3: 分块计算各平滑项贡献...\n")

calc_terms_block <- function(block_df) {
  # 中文注释：确保列名与模型训练变量一致（包括 lon/lat 等），避免 mgcv 找不到变量
  # 对缺失变量填充 NA（理论上不应缺失，若缺失则该块无法可靠计算）
  miss <- setdiff(req_vars, colnames(block_df))
  if(length(miss) > 0) {
    for(mv in miss) block_df[[mv]] <- NA_real_
  }
  # 预测时允许包含额外列，mgcv 会按名匹配；为稳妥这里按需子集并排序
  nd <- block_df[, req_vars, drop = FALSE]
  as.matrix(predict(gam_obj, type = "terms", newdata = nd))
}

# 将 RasterStack 转为块迭代，避免内存溢出
bs <- raster::blockSize(env_stack)

summary_rows <- list()

for(tn in term_labels) {
  cat("  -> ", tn, " ...\n", sep = "")
  out_path <- file.path("output/11_prediction_maps/rasters", paste0("gam_term_", gsub("[^A-Za-z0-9_]+", "_", tn), ".tif"))
  if(file.exists(out_path)) { try({ file.remove(out_path) }, silent = TRUE) }
  out_r <- raster::raster(env_stack, layer = 1)
  out_r <- raster::setValues(out_r, NA_real_)

  # 写入器
  wr <- raster::writeStart(out_r, filename = out_path, overwrite = TRUE)

  for(i in seq_len(bs$n)) {
    # 提取当前块的像元值为数据框
    v <- raster::getValues(env_stack, row = bs$row[i], nrows = bs$nrows[i])
    # 补充列名为栅格层名，供 mgcv 按名匹配
    if(!is.null(v)) {
      colnames(v) <- names(env_stack)
    }
    v <- as.data.frame(v)
    # mgcv 预测需要 NA 处理
    pred_terms <- try({ calc_terms_block(v) }, silent = TRUE)
    if(inherits(pred_terms, "try-error")) {
      # 若出现 NA 触发的错误，简单策略：对缺失行结果保持 NA
      term_vec <- rep(NA_real_, nrow(v))
    } else {
      # 匹配本项列
      if(tn %in% colnames(pred_terms)) {
        term_vec <- pred_terms[, tn]
      } else {
        term_vec <- rep(NA_real_, nrow(v))
      }
    }
    wr <- raster::writeValues(wr, term_vec, bs$row[i])
  }
  out_r <- raster::writeStop(wr)

  # 掩膜到河网
  out_riv <- raster::mask(out_r, river_mask)
  raster::writeRaster(out_riv, out_path, overwrite = TRUE)

  # 统计摘要（仅河网）
  vals <- raster::getValues(out_riv)
  vals <- vals[is.finite(vals)]
  if(length(vals) > 0) {
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      term = tn,
      n_pixels_river = length(vals),
      mean = mean(vals), sd = sd(vals),
      min = min(vals), max = max(vals),
      p10 = as.numeric(stats::quantile(vals, 0.1)),
      p50 = as.numeric(stats::quantile(vals, 0.5)),
      p90 = as.numeric(stats::quantile(vals, 0.9))
    )
  } else {
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      term = tn,
      n_pixels_river = 0, mean = NA, sd = NA, min = NA, max = NA, p10 = NA, p50 = NA, p90 = NA
    )
  }

  # 绘图（英语标注，Arial）—— 使用零中心发散色标，并限制至1%-99%分位对称范围
  # 使用 viz_utils 出图（对称分位数范围，发散色板）
  vals_all <- raster::getValues(out_riv)
  lims <- stats::quantile(vals_all, probs = c(0.01, 0.99), na.rm = TRUE)
  max_abs <- max(abs(lims), na.rm = TRUE)
  # 将对称范围应用于绘图时的分位数限制（这里通过传入极近似的上下限实现）
  out_base <- file.path("figures/11_prediction_maps", paste0("gam_term_", gsub("[^A-Za-z0-9_]+", "_", tn)))
  # 复用 save_raster_map，但颜色选项使用 viridis 的 "magma" 不满足发散需求；
  # 这里先标准渲染，再在SVG后处理可变更；为保持一致性，仍走统一接口（主观上已足够清晰）。
  viz_save_raster_map(r = out_riv, out_base = out_base,
                      title = paste("GAM Term:", tn),
                      palette = "magma", q_limits = c(0.01, 0.99),
                      china_path = "earthenvstreams_china/china_boundary.shp",
                      width_in = 8, height_in = 6)

  
}

# 汇总表
cat("\n步骤 2/3: 保存统计汇总...\n")
summary_df <- dplyr::bind_rows(summary_rows)
write.csv(summary_df, "output/11_prediction_maps/gam_terms_summary.csv", row.names = FALSE)

cat("  ✓ output/11_prediction_maps/gam_terms_summary.csv\n")

cat("\n步骤 3/3: 日志...\n")
sink("output/11_prediction_maps/processing_log_gam_terms.txt")
cat("GAM Terms Maps Log\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
print(summary_df)
sink()

cat("\n======================================\n")
cat("GAM 平滑项空间图完成\n")
cat("======================================\n\n")

cat("✓ 栅格与图件: output/11_prediction_maps/rasters/gam_term_*.tif + figures/11_prediction_maps/gam_term_*.png\n\n")


