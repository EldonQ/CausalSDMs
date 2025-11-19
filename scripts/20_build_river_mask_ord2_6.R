#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 20_build_river_mask_ord2_6.R
# 功能说明: 采用 A 方案——直接使用 earthenvstreams_china/flow_acc.tif 的第2波段
#          （Flow Accumulation）构建“一像元宽”的河网骨架掩膜（>0即为河网）。
#          掩膜与参考波段严格对齐（CRS/范围/分辨率），并输出高质量预览图。
# 重要说明: 本掩膜仅用于“出图阶段”的可视化与裁剪；建模与预测时的主数据域
#          仍以 flow_acc 第2波段为准（>0）。本脚本不再依赖 HydroRIVERS。
# 输入文件: 
#   - 参考栅格: earthenvstreams_china/flow_acc.tif（第2波段）
#   - 线要素:   E:/HydroSHEDS/HydroRIVERS_v10_as_shp/*.shp（包含字段 ORD_FLOW）
# 输出文件:
#   - 掩膜:   earthenvstreams_china/river_mask_ord2_6.tif  (0/1)
#   - 预览图: figures/00_china_env_variables/river_mask_ord2_6_preview.png (1200 dpi)
# 作者: Nature级别科研项目
# 日期: 2025-11-01
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/SDM01")

# 加载必要包（优先使用 terra 以提速；引入 ggplot2/svglite 统一出图风格）
packages <- c("raster", "terra", "tidyverse", "ggplot2", "svglite")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

dir.create("figures/00_china_env_variables", showWarnings = FALSE, recursive = TRUE)

# 引入统一绘图工具（主题/导出/预览）
source("scripts/visualization/viz_utils.R")

cat("\n======================================\n")
cat("构建 2–6 级河网掩膜 (HydroRIVERS → 栅格)\n")
cat("======================================\n\n")

# 路径配置（如有变动，可在此修改）
REF_RASTER_PATH <- "earthenvstreams_china/flow_acc.tif"  # 参考：第2波段
OUT_MASK_PATH_A <- "earthenvstreams_china/river_mask_skeleton.tif"       # A方案命名
OUT_MASK_PATH_COMPAT <- "earthenvstreams_china/river_mask_ord2_6.tif"    # 兼容旧路径（覆盖写入）
OUT_PREVIEW_BASE <- "figures/00_china_env_variables/river_skeleton_preview"

# 读取参考栅格（第2波段作为对齐基准，使用 terra 提速）
ref_brick <- terra::rast(REF_RASTER_PATH)
ref_r     <- ref_brick[[2]]
ref_crs_terra <- terra::crs(ref_r)
ref_res   <- terra::res(ref_r)[1]  # 取经度方向像元宽度（度）

cat("参考栅格: ", REF_RASTER_PATH, " | 分辨率=", paste(terra::res(ref_r), collapse=","),
    " | 尺寸=", paste(ncol(ref_r), "x", nrow(ref_r)), "\n", sep = "")

# A 方案：直接基于 flow_acc band2 > 0 构建河网骨架掩膜
cat("构建河网骨架掩膜 (flow_acc band2 > 0)...\n")
fa2 <- ref_brick[[2]]
mask_terra <- fa2
mask_terra[mask_terra <= 0] <- NA
mask_terra[mask_terra > 0]  <- 1

# 保存掩膜（写入A方案命名与兼容旧路径两份文件）
cat("保存掩膜...\n")
terra::writeRaster(mask_terra, OUT_MASK_PATH_A, overwrite = TRUE, datatype = "INT1U", gdal = c("COMPRESS=LZW"))
terra::writeRaster(mask_terra, OUT_MASK_PATH_COMPAT, overwrite = TRUE, datatype = "INT1U", gdal = c("COMPRESS=LZW"))
cat("\n✓ 已生成掩膜: ", OUT_MASK_PATH_A, "\n", sep = "")
cat("  - 兼容写入: ", OUT_MASK_PATH_COMPAT, "\n", sep = "")
cat("  - 覆盖像元数: ", sum(!is.na(terra::values(mask_terra))), " / ", terra::ncell(mask_terra), "\n", sep = "")

# 预览图：用 log1p(flow_acc2) 的强度进行渲染（更贴近“热图”质感）
cat("生成预览图 (PNG+SVG, 1200dpi, Arial)...\n")
try({
  viz_preview_river_intensity(flow_acc_path = REF_RASTER_PATH,
                              out_base = OUT_PREVIEW_BASE,
                              china_path = "earthenvstreams_china/china_boundary.shp")
  cat("✓ 预览图: ", paste0(OUT_PREVIEW_BASE, ".png"), "\n", sep = "")
}, silent = TRUE)

cat("\n======================================\n")
cat("掩膜构建完成 (仅用于出图阶段)\n")
cat("======================================\n\n")


