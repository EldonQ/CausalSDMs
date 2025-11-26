#!/usr/bin/env Rscript
# ==============================================================================
# 文件名称: viz_utils.R
# 功能说明: 提供Nature期刊级别的统一绘图工具函数（主题/导出/河网热图）
# 适用范围: 本工程所有需要出地图或热图的脚本均应引入本工具
# 使用方法: 在脚本开头添加 source("scripts/visualization/viz_utils.R")
# 重要规范: 图内文字一律英文标注；输出PNG(1200dpi)与SVG；字体为Arial
# 作者: Nature级别科研项目
# 日期: 2025-11-03
# ==============================================================================

# ------------------------------
# 依赖加载（按需安装）
# ------------------------------
required_pkgs <- c("terra", "sf", "ggplot2", "viridis", "RColorBrewer", "svglite", "sysfonts", "showtext", "rlang")
for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ------------------------------
# 字体与DPI设置：统一为Arial + 2400dpi
# ------------------------------
viz_ensure_arial <- function() {
  # 中文注释：在Windows环境注册Arial，启用showtext以保证嵌入与一致渲染
  try({
    sysfonts::font_add(
      family = "Arial",
      regular = "C:/Windows/Fonts/arial.ttf",
      bold = "C:/Windows/Fonts/arialbd.ttf",
      italic = "C:/Windows/Fonts/ariali.ttf",
      bolditalic = "C:/Windows/Fonts/arialbi.ttf"
    )
    showtext::showtext_opts(dpi = 2400)
    showtext::showtext_auto(enable = TRUE)
  }, silent = TRUE)
}

viz_ensure_arial()

# 避免R CMD check关于未绑定变量的提示（ggplot2美学映射）
utils::globalVariables(c("x", "y", "value", ".data"))

# ------------------------------
# Nature 期刊风格主题
# ------------------------------
viz_theme_nature <- function(base_size = 9, base_family = "Arial",
                             title_size = 9, subtitle_size = 8,
                             axis_title_size = 8, axis_text_size = 6,
                             legend_title_size = 7, legend_text_size = 7) {
  ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = title_size, face = "bold", hjust = 0.5),
      plot.subtitle = ggplot2::element_text(size = subtitle_size, hjust = 0.5),
      axis.text = ggplot2::element_text(size = axis_text_size),
      axis.title = ggplot2::element_text(size = axis_title_size, face = "bold"),
      legend.title = ggplot2::element_text(size = legend_title_size, face = "bold"),
      legend.text = ggplot2::element_text(size = legend_text_size),
      legend.position = "right",
      # 无网格线、透明背景（便于在期刊版面叠加使用）
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_blank(),
      plot.background = ggplot2::element_rect(fill = "transparent", color = NA),
      panel.background = ggplot2::element_rect(fill = "transparent", color = NA)
    )
}

# ------------------------------
# 读取中国边界（若可用），用于叠加边界线
# ------------------------------
viz_read_china <- function(path = "earthenvstreams_china/china_boundary.shp", crs_target = NULL) {
  china <- try(sf::st_read(path, quiet = TRUE), silent = TRUE)
  if (inherits(china, "try-error")) return(NULL)
  if (!is.null(crs_target)) {
    china <- sf::st_transform(china, crs = crs_target)
  }
  return(china)
}

# ------------------------------
# 将Raster/SpatRaster统一为SpatRaster
# ------------------------------
viz_as_spatraster <- function(r) {
  if (inherits(r, "SpatRaster")) return(r)
  return(terra::rast(r))
}

# ------------------------------
# 根据分位数裁剪色标范围（避免极端值主导）
# ------------------------------
viz_quantile_limits <- function(values, q = c(0.01, 0.99)) {
  v <- values[!is.na(values)]
  if (length(v) == 0) return(c(NA, NA))
  lo <- as.numeric(stats::quantile(v, q[1], na.rm = TRUE))
  hi <- as.numeric(stats::quantile(v, q[2], na.rm = TRUE))
  if (!is.finite(lo) || !is.finite(hi)) return(range(v, na.rm = TRUE))
  if (lo >= hi) return(range(v, na.rm = TRUE))
  c(lo, hi)
}

# ------------------------------
# 保存地图：PNG(2400dpi)+SVG，英文标注，Arial
# ------------------------------
viz_save_raster_map <- function(
  r,                       # RasterLayer/SpatRaster：待绘制的单层栅格（已按需要掩膜）
  out_base,                # 字符串：输出文件基名（不含扩展名）
  title = "",             # 字符串：图标题（英文）
  palette = "magma",      # viridis 选项："magma"/"viridis"/"plasma"/"inferno" 等
  q_limits = c(0.01,0.99), # 分位数裁剪区间
  china_path = "earthenvstreams_china/china_boundary.shp",
  width_in = 8, height_in = 6,
  theme_base_size = 9, title_size = 9,
  scale_trans = "identity",  # 连续色标的变换："identity"/"sqrt"/"log10" 等
  hydrorivers_sf = NULL
) {
  # 中文注释：该函数将输入栅格转换为数据框后使用ggplot2绘制；
  #           采用分位数裁剪色标范围，叠加中国国界，导出PNG与SVG两种格式。

  r_spat <- viz_as_spatraster(r)
  china <- viz_read_china(china_path, crs_target = terra::crs(r_spat))

  # 计算色标上下限
  vals <- terra::values(r_spat)
  lims <- viz_quantile_limits(vals, q_limits)
  # 转为数据框
  df <- as.data.frame(r_spat, xy = TRUE, na.rm = TRUE)
  colnames(df)[3] <- "value"
  if (nrow(df) == 0) {
    warning("待绘制栅格为空: ", out_base)
    return(invisible(FALSE))
  }

  # 颜色
  cols <- viridis::viridis(200, option = palette)

  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = df, ggplot2::aes_string(x = "x", y = "y", fill = "value")) +
    {
      if (!is.null(china)) ggplot2::geom_sf(data = china, fill = NA, color = "black", size = 0.3, alpha = 0.9) else NULL
    } +
    {
      if (!is.null(hydrorivers_sf) && inherits(hydrorivers_sf, "sf") &&
          nrow(hydrorivers_sf) > 0 && "ORD_FLOW" %in% names(hydrorivers_sf)) {
        hyd_sf <- hydrorivers_sf
        hyd_sf_ord2 <- hyd_sf[hyd_sf$ORD_FLOW == 2, ]
        hyd_sf_ord3 <- hyd_sf[hyd_sf$ORD_FLOW == 3, ]
        hyd_sf_ord4 <- hyd_sf[hyd_sf$ORD_FLOW == 4, ]
        hyd_sf_ord5 <- hyd_sf[hyd_sf$ORD_FLOW == 5, ]
        hyd_sf_ord6 <- hyd_sf[hyd_sf$ORD_FLOW == 6, ]
        list(  # 各阶河流的 linewidth（线宽）例如把全部乘以 1.5，让骨干线更粗
          if (nrow(hyd_sf_ord2) > 0) ggplot2::geom_sf(data = hyd_sf_ord2, color = "dodgerblue3",  # color = "dodgerblue3"：可改成你喜欢的蓝色
                                                      linewidth = 0.15, alpha = 0.9) else NULL,
          if (nrow(hyd_sf_ord3) > 0) ggplot2::geom_sf(data = hyd_sf_ord3, color = "dodgerblue3",
                                                      linewidth = 0.20, alpha = 0.9) else NULL,
          if (nrow(hyd_sf_ord4) > 0) ggplot2::geom_sf(data = hyd_sf_ord4, color = "dodgerblue3",
                                                      linewidth = 0.28, alpha = 0.9) else NULL,
          if (nrow(hyd_sf_ord5) > 0) ggplot2::geom_sf(data = hyd_sf_ord5, color = "dodgerblue3",
                                                      linewidth = 0.36, alpha = 0.9) else NULL,
          if (nrow(hyd_sf_ord6) > 0) ggplot2::geom_sf(data = hyd_sf_ord6, color = "dodgerblue3",
                                                      linewidth = 0.45, alpha = 0.9) else NULL
        )
      } else NULL
    } +
    ggplot2::scale_fill_gradientn(colors = cols, limits = lims, oob = scales::squish, name = "", trans = scale_trans) +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::labs(title = title, x = "Longitude (°E)", y = "Latitude (°N)") +
    viz_theme_nature(base_size = theme_base_size, title_size = title_size)

  # 导出PNG（透明背景 + 2400dpi）
  ggplot2::ggsave(filename = paste0(out_base, ".png"), plot = p,
                  width = width_in, height = height_in, dpi = 2400,
                  units = "in", bg = "transparent")
  # 导出SVG（透明背景）
  ggplot2::ggsave(filename = paste0(out_base, ".svg"), plot = p,
                  width = width_in, height = height_in, units = "in",
                  bg = "transparent", device = "svg")

  invisible(TRUE)
}

# ------------------------------
# 基于 flow_acc 第2波段构建河网骨架掩膜（A方案）
# ------------------------------
viz_build_river_skeleton <- function(flow_acc_path = "earthenvstreams_china/flow_acc.tif",
                                     out_mask_path = "earthenvstreams_china/river_mask_skeleton.tif") {
  # 中文注释：以 flow_acc.tif 的第2波段作为河网骨架；>0 视为河网，生成 0/1 掩膜（非河网为NA便于透明显示）
  r <- terra::rast(flow_acc_path)[[2]]
  mask <- r
  mask[mask <= 0] <- NA
  mask[mask > 0]  <- 1
  terra::writeRaster(mask, out_mask_path, overwrite = TRUE, datatype = "INT1U", gdal = c("COMPRESS=LZW"))
  return(mask)
}

# ------------------------------
# 预览：基于流量累积构造的“河网热度”图（对数增强 + 视觉加粗）
# ------------------------------
viz_preview_river_intensity <- function(flow_acc_path = "earthenvstreams_china/flow_acc.tif",
                                        out_base = "figures/00_china_env_variables/river_skeleton_preview",
                                        china_path = "earthenvstreams_china/china_boundary.shp") {
  # 中文注释：预览目的仅为观感评估。对 flow_acc 第2波段取 log1p 并归一化到 [0,1]，
  #           然后通过局部最大值卷积略微“加粗”河网，以便在全国范围查看时更显眼。
  r <- terra::rast(flow_acc_path)[[2]]
  # 仅限河网像元
  r[r <= 0] <- NA
  # 对数增强
  r_log <- log1p(r)
  # 归一化
  mx <- suppressWarnings(as.numeric(terra::global(r_log, "max", na.rm = TRUE)[1, 1]))
  if (is.finite(mx) && mx > 0) {
    r_log <- r_log / mx
  }

  # 通过局部最大值平滑略微加粗河网（仅用于预览，不影响掩膜本身的一像元宽结构）
  w <- matrix(1, nrow = 3, ncol = 3)  # 把 nrow = 3, ncol = 3 改成 5 / 7 → 预览骨架更粗
  r_thick <- terra::focal(r_log, w = w, fun = max, na.rm = TRUE)

  viz_save_raster_map(
    r_thick,
    out_base = out_base,
    title = "River Skeleton (from Flow Accumulation)",
    palette = "viridis",
    q_limits = c(0.05, 0.99),  # 左边值越大，背景越淡、亮的河段越突出
    china_path = china_path
  )
}

viz_thicken_river_raster <- function(r, w_size = 3) {
  r_spat <- viz_as_spatraster(r)
  w <- matrix(1, nrow = w_size, ncol = w_size)
  terra::focal(r_spat, w = w, fun = max, na.rm = TRUE)
}

viz_load_hydrorivers_ord2_6 <- function(
  shp_dir = "E:/HydroSHEDS/HydroRIVERS_v10_as_shp",
  ord_min = 2,
  ord_max = 6,
  crs_target = NULL
) {
  if (!dir.exists(shp_dir)) {
    warning("shp_dir not found: ", shp_dir)
    return(NULL)
  }
  shp_files <- list.files(shp_dir, pattern = "\\.shp$", full.names = TRUE)
  if (length(shp_files) == 0) {
    warning("no .shp files found in ", shp_dir)
    return(NULL)
  }
  rivers_list <- lapply(shp_files, function(f) {
    x <- try(sf::st_read(f, quiet = TRUE), silent = TRUE)
    if (inherits(x, "try-error")) return(NULL)
    x
  })
  rivers_list <- rivers_list[!vapply(rivers_list, is.null, logical(1))]
  if (length(rivers_list) == 0) return(NULL)
  rivers <- do.call(rbind, rivers_list)
  if (!"ORD_FLOW" %in% names(rivers)) {
    warning("ORD_FLOW field not found in HydroRIVERS data")
    return(NULL)
  }
  rivers <- rivers[!is.na(rivers$ORD_FLOW) & rivers$ORD_FLOW >= ord_min & rivers$ORD_FLOW <= ord_max, ]
  if (!is.null(crs_target)) {
    rivers <- sf::st_transform(rivers, crs = crs_target)
  }
  rivers
}

# ==============================================================================
# 结束
# ==============================================================================


