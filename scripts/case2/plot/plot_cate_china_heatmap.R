################################################################################
# plot_cate_china_heatmap.R
# 中国范围 CATE 平滑热图：地理边界使用中国 SHP；CATE 点 → 栅格化 → 平滑渲染
#
# 为什么不用 gstat::idw 做全格网插值？
# - 0.02° 分辨率下格点数可达数百万，idw 单线程很慢，容易看起来“卡住”
# - 本脚本改用 terra 的栅格工作流：rasterize → approxNA → focal(gaussian)
#   速度更快、渲染更平滑，更接近你给的参考图（teal→cream→brown 连续渐变）
#
# 依赖（项目内已有同类用法）：
# - terra, sf, tidyterra, ggplot2, scales
# - 中国边界：earthenvstreams_china/china_boundary.shp（若无则用 rnaturalearth）
#
# 输入：
# - CATE CSV：必须包含 lon/lat/cate（中国边界为经纬度；若是 case2 的 x/y 投影坐标请先转换）
#
# 输出（高清）：
# - figures/case2/plot/cate_map_china_<var>.png (2400 dpi) + .svg
#
# 运行：
#   setwd("E:/CausalSDMs")
#   source("scripts/case2/plot/plot_cate_china_heatmap.R")
################################################################################

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

# -----------------------------------------------------------------------------
# 1. 包与路径
# -----------------------------------------------------------------------------
pkgs <- c("terra", "sf", "tidyterra", "ggplot2", "tidyverse", "scales", "rnaturalearth")
for (p in pkgs) {
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
        install.packages(p, repos = "https://cloud.r-project.org")
        library(p, character.only = TRUE)
    }
}

# 中国范围（与 16_revised_figure6 等一致）
CHINA_XLIM <- c(73.5, 135.0)
CHINA_YLIM <- c(18.0, 53.5)

# 可配置：CATE CSV 路径（默认 14_causal；可改为 case2/cate 下文件）
CATE_CSV <- "output/14_causal/d_cate/cate_result_flow_length.csv"
# 若 case2：CATE_CSV <- "output/case2/cate/cate_result_SWI_swi10_ann_mean_temp.csv"

FIG_DIR <- "figures/case2/plot"
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# 2. 中国边界（优先 SHP，与 scripts 内用法一致）
# -----------------------------------------------------------------------------
china_path <- "earthenvstreams_china/china_boundary.shp"
china_sf <- NULL
if (file.exists(china_path)) {
    china_sf <- sf::st_read(china_path, quiet = TRUE)
    if (!is.null(china_sf)) sf::st_crs(china_sf) <- 4326
}
if (is.null(china_sf) || nrow(china_sf) == 0) {
    if (requireNamespace("rnaturalearth", quietly = TRUE)) {
        china_sf <- rnaturalearth::ne_countries(country = "China", scale = "medium", returnclass = "sf")
    }
}
if (is.null(china_sf) || nrow(china_sf) == 0) {
    stop("Could not load China boundary (SHP or rnaturalearth).")
}

# -----------------------------------------------------------------------------
# 3. 读取 CATE 点数据
# -----------------------------------------------------------------------------
if (!file.exists(CATE_CSV)) {
    stop("CATE CSV not found: ", CATE_CSV)
}
dat <- read.csv(CATE_CSV, stringsAsFactors = FALSE)
if ("is_valid" %in% names(dat)) {
    dat <- dat[dat$is_valid == TRUE, ]
}
if (!all(c("lon", "lat", "cate") %in% names(dat))) {
    stop("中国边界使用经纬度坐标系：CATE CSV 必须包含 lon/lat/cate 列。")
}
pts <- dat[, c("lon", "lat", "cate")]
pts <- pts[complete.cases(pts), ]
if (nrow(pts) < 10) {
    stop("Too few valid CATE points.")
}

# 变量标签（用于标题与文件名）
treat_var <- if ("treatment_var" %in% names(dat)) dat$treatment_var[1] else "CATE"
treat_var <- gsub("^cate_result_", "", treat_var)
var_label <- treat_var
var_label <- gsub("flow_length", "Flow Length", var_label)
var_label <- gsub("lc_wavg_09|lc_09", "Urbanization", var_label)

# -----------------------------------------------------------------------------
# 4. 栅格化 + 缺失填充 + 高斯平滑（快、平滑、高清）
# -----------------------------------------------------------------------------
# 分辨率：越小越清晰，但越慢；建议 0.03~0.06 之间先出图，再逐步加密
RES_DEG <- 0.04
SMOOTH_WINDOW <- 9   # 高斯核大小（奇数）：9/11/15 更平滑但更慢
SMOOTH_SIGMA  <- 2.0 # 高斯核 sigma
FILL_ITER     <- 3  # 用邻域均值填洞的迭代次数（不依赖 approxNA）

cat(sprintf("Raster resolution: %.3f deg | smooth window=%d sigma=%.2f\n",
    RES_DEG, SMOOTH_WINDOW, SMOOTH_SIGMA))

# 中国边界转 terra 向量（用于 mask）
china_v <- terra::vect(china_sf)
terra::crs(china_v) <- "EPSG:4326"

# 构造模板栅格（经纬度）
r_tmpl <- terra::rast(
    xmin = CHINA_XLIM[1], xmax = CHINA_XLIM[2],
    ymin = CHINA_YLIM[1], ymax = CHINA_YLIM[2],
    resolution = RES_DEG, crs = "EPSG:4326"
)

# 点 → SpatVector
pts_v <- terra::vect(pts, geom = c("lon", "lat"), crs = "EPSG:4326")

# 栅格化（每像元取 mean CATE）
r0 <- terra::rasterize(pts_v, r_tmpl, field = "cate", fun = "mean", background = NA_real_)

# 掩膜到中国境内（先 mask 再补洞/平滑，减少不必要计算）
r0 <- terra::mask(r0, china_v)

# 补洞：用邻域均值填充 NA（不依赖 terra::approxNA，兼容旧版 terra）
w3 <- matrix(1, 3, 3) / 9
r1 <- r0
for (k in seq_len(FILL_ITER)) {
    r_f <- terra::focal(r1, w = w3, fun = "mean", na.rm = TRUE, fillvalue = NA_real_)
    r1 <- terra::ifel(is.na(r1), r_f, r1)
}

# 高斯平滑：得到参考图那种“尽量平滑”的面效果
w <- terra::focalMat(r1, d = SMOOTH_SIGMA * RES_DEG, type = "Gauss")
r_smooth <- terra::focal(r1, w = w, fun = "sum", na.rm = TRUE, fillvalue = NA_real_)

# 再次 mask，确保边界外为 NA
r_smooth <- terra::mask(r_smooth, china_v)

# 色阶范围：对称裁剪（95% 分位，避免极端值主导）
vlim <- as.numeric(stats::quantile(abs(terra::values(r_smooth)), 0.98, na.rm = TRUE))
if (!is.finite(vlim) || vlim <= 0) vlim <- max(abs(terra::values(r_smooth)), na.rm = TRUE)
if (!is.finite(vlim) || vlim <= 0) vlim <- 1

# -----------------------------------------------------------------------------
# 5. 绘图：参考图 2 风格（teal→cream→brown），高清平滑
# -----------------------------------------------------------------------------
# 配色：负值=teal，零=cream，正值=brown（连续渐变）
cols <- c("#004C5A", "#2A9D8F", "#BFE7DD", "#F4EAD5", "#E9D8A6", "#B08968", "#7F5539", "#5D3A1A")
vals <- c(0, 0.18, 0.33, 0.50, 0.58, 0.70, 0.85, 1.0)

p <- ggplot() +
    tidyterra::geom_spatraster(data = r_smooth, maxcell = 1e7) +
    scale_fill_gradientn(
        colours = cols,
        values = vals,
        limits = c(-vlim, vlim),
        oob = scales::squish,
        na.value = NA,
        name = "CATE",
        guide = guide_colorbar(
            title.position = "top",
            title.hjust = 0.5,
            barheight = unit(4, "cm"),
            barwidth = unit(0.55, "cm"),
            frame.colour = "black",
            frame.linewidth = 0.3,
            ticks.colour = "black",
            ticks.linewidth = 0.4
        )
    ) +
    tidyterra::geom_spatvector(data = china_v, fill = NA, color = "black", linewidth = 0.5) +
    coord_sf(xlim = CHINA_XLIM, ylim = CHINA_YLIM, expand = FALSE) +
    labs(
        title = paste0("CATE spatial map: ", var_label),
        x = NULL, y = NULL
    ) +
    theme_void(base_family = "Arial", base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 14, hjust = 0),
        legend.position = "right",
        legend.title = element_text(face = "bold", size = 11),
        legend.text  = element_text(size = 9)
    )

# -----------------------------------------------------------------------------
# 6. 保存（高清晰 2400 dpi PNG + SVG）
# -----------------------------------------------------------------------------
stem <- paste0("cate_map_china_", gsub("[^A-Za-z0-9_]", "_", treat_var))
out_png <- file.path(FIG_DIR, paste0(stem, ".png"))
out_svg <- file.path(FIG_DIR, paste0(stem, ".svg"))
ggsave(out_png, p, width = 10, height = 7.5, dpi = 2400, bg = "white")
ggsave(out_svg, p, width = 10, height = 7.5, bg = "white")
cat("Saved:", out_png, "\n")
cat("Saved:", out_svg, "\n")
