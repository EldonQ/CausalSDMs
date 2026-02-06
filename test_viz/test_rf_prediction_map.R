#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: test_rf_prediction_map.R
# 功能说明: 快速测试河网可视化参数 - 仅使用RF模型
# 输出文件: test_viz/output/pred_rf.tif, test_viz/figures/prediction_rf.png
# 作者: 测试脚本
# 日期: 2025-12-05
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/CausalSDMs")

# 加载必要的包
packages <- c("tidyverse", "sf", "ggplot2", "viridis", "rnaturalearth", "raster", "randomForest", "sysfonts", "showtext", "terra", "svglite")
for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
}

# 字体设置
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

# 创建测试输出目录
dir.create("test_viz/output", showWarnings = FALSE, recursive = TRUE)
dir.create("test_viz/figures", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("快速测试：RF模型预测地图 + 河网增强\n")
cat("======================================\n\n")

# ===== 步骤 1: 构建环境栅格栈 =====
cat("步骤 1/4: 构建环境栅格栈...\n")

sel_vars <- read.csv("output/04_collinearity/selected_variables.csv", stringsAsFactors = FALSE)$variable
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
cat("  ✓ 变量层数: ", raster::nlayers(env_stack), "\n", sep = "")

# 单位对齐
idx_temp <- which(grepl("^hydro_wavg_0[1-9]$|^hydro_wavg_1[01]$", names(env_stack)))
if (length(idx_temp) > 0) {
    for (i in idx_temp) env_stack[[i]] <- env_stack[[i]] / 10
}
idx_slope <- which(grepl("^slope_", names(env_stack)))
if (length(idx_slope) > 0) {
    for (i in idx_slope) env_stack[[i]] <- env_stack[[i]] / 100
}
idx_ph <- which(names(env_stack) == "soil_wavg_02")
if (length(idx_ph) > 0) {
    env_stack[[idx_ph]] <- env_stack[[idx_ph]] / 10
}

# ===== 步骤 2: 河网掩膜 =====
cat("\n步骤 2/4: 生成河网掩膜...\n")
fa <- raster::brick("earthenvstreams_china/flow_acc.tif")[[2]]
fa_vals <- raster::getValues(fa)
fa_vals[fa_vals <= 0] <- NA
river_mask <- raster::setValues(fa, fa_vals)
rm(fa_vals)
cat("  ✓ 河网掩膜已生成\n")

# ===== 步骤 3: RF模型预测 =====
cat("\n步骤 3/4: RF模型预测...\n")

rf_model_path <- "output/06_model_rf/model.rds"
if (!file.exists(rf_model_path)) {
    stop("RF模型文件不存在: ", rf_model_path)
}

mdl <- readRDS(rf_model_path)

# RF预测函数
pred_fun <- function(m, df) {
    as.numeric(predict(m, newdata = df, type = "prob")[, "1"])
}

# 栅格化预测
tif_path <- "test_viz/output/pred_rf.tif"
if (file.exists(tif_path)) try(file.remove(tif_path), silent = TRUE)

pred_r <- raster::predict(env_stack,
    model = mdl, fun = pred_fun, filename = tif_path,
    overwrite = TRUE, progress = "text"
)
pred_r <- clamp(pred_r, lower = 0, upper = 1, useValues = TRUE)
pred_r_river <- mask(pred_r, river_mask)

# 保存掩膜后的栅格
tif_mask_path <- "test_viz/output/pred_rf_river.tif"
if (file.exists(tif_mask_path)) try(file.remove(tif_mask_path), silent = TRUE)
raster::writeRaster(pred_r_river, tif_mask_path, overwrite = TRUE)
cat("  ✓ 预测栅格已保存\n")

# ===== 步骤 4: 绘图 =====
cat("\n步骤 4/4: 绘制地图（使用河网增强）...\n")

# ==========================================
# ★★★ 在此调整河网参数快速测试 ★★★
# ==========================================
# 可以直接在这里覆盖全局参数进行测试：
# RIVER_ENHANCE_METHOD <- "vector"  # "vector" / "raster" / "both" / "none"
# RIVER_COLOR <- "dodgerblue3"
# RIVER_ALPHA <- 0.85
# RIVER_WIDTH_SCALE <- 1.0  # 增大让河流更粗

out_base <- "test_viz/figures/prediction_rf"
viz_save_raster_map(
    r = pred_r_river,
    out_base = out_base,
    title = "RF Model Prediction (Test)",
    palette = "magma",
    q_limits = c(0.01, 0.99),
    china_path = "earthenvstreams_china/china_boundary.shp",
    width_in = 6.5,
    height_in = 5.0,
    theme_base_size = 8,
    title_size = 9,
    scale_trans = "sqrt"
    # 可以在此覆盖河网参数：
    # river_enhance_method = "vector",
    # river_color = "steelblue",
    # river_alpha = 0.9,
    # river_width_scale = 1.5
)

cat("\n======================================\n")
cat("测试完成!\n")
cat("======================================\n\n")
cat("输出文件:\n")
cat("  - 栅格: test_viz/output/pred_rf_river.tif\n")
cat("  - 图像: test_viz/figures/prediction_rf.png\n")
cat("  - 矢量: test_viz/figures/prediction_rf.svg\n\n")
cat("提示: 修改 viz_utils.R 第 23-117 行的河网参数后重新运行此脚本测试效果\n\n")
