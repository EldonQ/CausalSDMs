#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 11_current_prediction_maps.R
# 功能说明: 生成当前气候条件下所有模型的物种分布预测地图
# 方法: 使用4个模型（Maxnet, NN, RF, GAM）生成预测概率分布图
# 输入文件: output/05_model_maxnet/predictions.csv
#          output/05b_model_nn/predictions.csv
#          output/06_model_rf/predictions.csv
#          output/07_model_gam/predictions.csv
# 输出文件: figures/11_prediction_maps/prediction_*.png
#          output/11_prediction_maps/prediction_summary.csv
# 作者: Nature级别科研项目
# 日期: 2025-10-20
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/SDM01")

# 加载必要的包（本脚本改为基于河网的栅格热图绘制）
# 说明：为满足期刊级绘图要求，采用与训练变量一致的RasterStack进行整幅预测，
#       并用 flow_acc.tif 第2波段作为河网掩膜，仅在河网上显示概率热图。
packages <- c("tidyverse", "sf", "ggplot2", "viridis", "rnaturalearth", "raster", "maxnet", "nnet", "randomForest", "mgcv", "sysfonts", "showtext", "terra", "svglite")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 中文注释：注册并启用 Arial 字体，确保 PDF/PNG 输出嵌入或渲染为 Arial
# 在 Windows 上 Arial 通常位于 C:/Windows/Fonts；为稳妥提供常见字重路径
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

# 统一绘图工具（Nature风格、PNG+SVG导出、Arial）：用于所有地图输出
source("scripts/visualization/viz_utils.R")

dir.create("output/11_prediction_maps", showWarnings = FALSE, recursive = TRUE)
dir.create("output/11_prediction_maps/rasters", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/11_prediction_maps", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("当前气候预测地图绘制（河网栅格热图）\n")
cat("======================================\n\n")

### ========================= 新流程（核心重构） ========================= ###
# 1. 构建与训练变量一致的环境RasterStack（对齐 earthenvstreams_china）
cat("步骤 1/5: 构建环境栅格栈 (与训练变量一致)...\n")

# 读取最终用于建模的变量名（应为47个）
sel_vars <- read.csv("output/04_collinearity/selected_variables.csv", stringsAsFactors = FALSE)$variable
cat("  ✓ 入模变量数: ", length(sel_vars), "\n", sep = "")
if(length(sel_vars) != 47) {
  cat("  ⚠ 注意: 入模变量数 非 47 (当前=", length(sel_vars), ")，将按白名单继续，但请核对上游。\n", sep = "")
}

# 读取变量-文件-波段映射（由02脚本输出），用于从磁盘装载对应波段
var_map <- read.csv("output/02_env_extraction/extracted_variables.csv", stringsAsFactors = FALSE)
var_map <- var_map[var_map$variable %in% sel_vars, c("variable", "file", "band")]

# 辅助函数：按文件分组读取，避免重复IO
build_env_stack <- function(var_map_df, base_dir = "earthenvstreams_china") {
  # 中文注释：将同一文件内所需波段一次性读取，再按变量名重命名并合并为统一RasterStack
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
  # 仅保留需用变量并按sel_vars排序，保证与模型一致
  stk <- stk[[sel_vars]]
  return(stk)
}

env_stack <- build_env_stack(var_map)
cat("  ✓ 变量层数: ", raster::nlayers(env_stack), " (应与训练一致)\n", sep = "")

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

# 追加经纬度图层以兼容包含 s(lon, lat, ...) 的GAM模型
# 中文注释：部分模型（如GAM）在训练中显式使用了经纬度项，
#          因此需要在预测时为每个像元提供 lon/lat，
#          这里利用现有栅格的地理参考生成对应的经纬度栅格，并命名为 lon/lat。
lon_r <- raster::init(env_stack[[1]], fun = 'x')
lat_r <- raster::init(env_stack[[1]], fun = 'y')
names(lon_r) <- 'lon'
names(lat_r) <- 'lat'
env_stack <- raster::addLayer(env_stack, lon_r, lat_r)

# 2. 构建河网掩膜（flow_acc.tif 第2波段 > 0）
cat("\n步骤 2/5: 生成河网掩膜...\n")
fa <- raster::brick("earthenvstreams_china/flow_acc.tif")[[2]]
fa_vals <- raster::getValues(fa)
fa_vals[fa_vals <= 0] <- NA
fa <- raster::setValues(fa, fa_vals)  # NA为非河网
river_mask <- fa
rm(fa_vals)
cat("  ✓ 河网掩膜已生成\n")

# 3. 装载模型并进行整幅栅格预测（概率0-1），仅在河网上可见
cat("\n步骤 3/5: 模型栅格化预测（仅河网）...\n")
models <- c("Maxnet", "NN", "RF", "GAM")
model_files <- c(
  Maxnet = "output/05_model_maxnet/model.rds",
  NN = "output/05b_model_nn/model.rds",
  RF = "output/06_model_rf/model.rds",
  GAM = "output/07_model_gam/model.rds"
)

# 中文注释：为不同模型定义预测函数，使其可被 raster::predict 以块处理方式调用
make_predict_fun <- function(model_name, model_obj) {
  if(model_name == "Maxnet") {
    # maxnet 直接对数据框预测logistic概率；函数签名须为 (model, df)
    return(function(m, df) { as.numeric(predict(m, df, type = "logistic")) })
  }
  if(model_name == "RF") {
    return(function(m, df) { as.numeric(predict(m, newdata = df, type = "prob")[, "1"]) })
  }
  if(model_name == "GAM") {
    return(function(m, df) { as.numeric(predict(m, newdata = df, type = "response")) })
  }
  if(model_name == "NN") {
    # NN保存了标准化参数；对传入df逐列标准化后预测
    # 注意：nnet的predict不支持raster::predict的块处理，需要特殊处理
    return(function(m, df) {
      mu <- m$mean; sdv <- m$sd; mod <- m$model; vars <- m$vars
      sdv[sdv == 0 | is.na(sdv)] <- 1
      x <- as.matrix(df[, vars, drop = FALSE])
      x <- sweep(x, 2, mu[vars], "-")
      x <- sweep(x, 2, sdv[vars], "/")
      # 使用nnet包的predict，明确指定为矩阵输入
      as.numeric(nnet:::predict.nnet(mod, x, type = "raw"))
    })
  }
}

# 读取中国边界用于叠加绘图（仅用于美观，不影响数据）
china <- ne_countries(country = "China", scale = "medium", returnclass = "sf")

summary_rows <- list()

for(m in models) {
  cat("  -> ", m, " ...\n", sep = "")
  if(!file.exists(model_files[[m]])) {
    cat("     ✗ 模型文件不存在: ", model_files[[m]], "\n", sep = "")
    next
  }
  mdl <- readRDS(model_files[[m]])
  pred_fun <- make_predict_fun(m, mdl)

  # 栅格化预测（分块处理，避免内存溢出），并保存GeoTIFF
  tif_path <- paste0("output/11_prediction_maps/rasters/pred_", tolower(m), ".tif")
  # 中文注释：为避免Windows上文件锁导致无法覆盖，若目标文件已存在则先尝试删除
  if(file.exists(tif_path)) { try({ file.remove(tif_path) }, silent = TRUE) }
  pred_r <- raster::predict(env_stack, model = mdl, fun = pred_fun, filename = tif_path,
                            overwrite = TRUE, progress = "text")
  # 限定概率范围到[0,1]
  pred_r <- clamp(pred_r, lower = 0, upper = 1, useValues = TRUE)
  # 仅保留河网像元
  pred_r_river <- mask(pred_r, river_mask)

  # 保存遮罩后的GeoTIFF（用于不确定性分析等后续图件）
  tif_mask_path <- paste0("output/11_prediction_maps/rasters/pred_", tolower(m), "_river.tif")
  if(file.exists(tif_mask_path)) { try({ file.remove(tif_mask_path) }, silent = TRUE) }
  raster::writeRaster(pred_r_river, tif_mask_path, overwrite = TRUE)

  # 统计摘要（仅河网像元）
  vals <- raster::getValues(pred_r_river)
  vals <- vals[!is.na(vals)]
  summary_rows[[length(summary_rows) + 1]] <- data.frame(
    model = m,
    n_pixels_river = length(vals),
    mean = mean(vals), sd = sd(vals),
    min = min(vals), max = max(vals),
    p10 = as.numeric(quantile(vals, 0.1)),
    p50 = as.numeric(quantile(vals, 0.5)),
    p90 = as.numeric(quantile(vals, 0.9))
  )

  # 绘制河网热图（统一使用 ggplot2 + viz_utils，自动输出 PNG 与 SVG）
  out_base <- paste0("figures/11_prediction_maps/prediction_", tolower(m))
  viz_save_raster_map(r = pred_r_river, out_base = out_base,
                      title = paste(m, "Model Prediction"),
                      palette = "magma", q_limits = c(0.01, 0.99),
                      china_path = "earthenvstreams_china/china_boundary.shp",
                      width_in = 6.5, height_in = 5.0,
                      theme_base_size = 8, title_size = 9,
                      scale_trans = "sqrt")

  

  cat("     ✓ 栅格与图已保存\n")
}

# 4. 保存统计表
cat("\n步骤 4/5: 保存预测统计表...\n")
summary_data <- bind_rows(summary_rows)
write.csv(summary_data, "output/11_prediction_maps/prediction_summary.csv", row.names = FALSE)
cat("  ✓ output/11_prediction_maps/prediction_summary.csv\n")

# 5. 记录处理日志
cat("\n步骤 5/5: 记录处理日志...\n")
sink("output/11_prediction_maps/processing_log.txt")
cat("预测地图绘制日志\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
cat("说明: 本脚本输出为河网栅格热图 (flow_acc band2 > 0 掩膜)\n\n")
print(summary_data)
sink()

# 日志
sink("output/11_prediction_maps/processing_log.txt")
cat("预测地图绘制日志\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
cat("模型数量: 4\n")
cat("预测摘要:\n")
print(summary_data)
sink()

cat("\n======================================\n")
cat("预测地图绘制完成\n")
cat("======================================\n\n")

cat("预测概率统计(河网像元):\n")
for(i in seq_len(nrow(summary_data))) {
  cat(sprintf("  %-8s: n=%d, mean=%.3f ± %.3f (range: %.3f - %.3f)\n",
              summary_data$model[i],
              summary_data$n_pixels_river[i],
              summary_data$mean[i],
              summary_data$sd[i],
              summary_data$min[i],
              summary_data$max[i]))
}

cat("\n✓ 脚本执行完成!\n\n")
