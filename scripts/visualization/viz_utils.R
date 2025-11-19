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
# 字体与DPI设置：统一为Arial + 1200dpi
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
    showtext::showtext_opts(dpi = 1200)
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
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(size = 0.3, color = "gray90"),
      plot.background = ggplot2::element_rect(fill = "white", color = NA),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
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
# 保存地图：PNG(1200dpi)+SVG，英文标注，Arial
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
  scale_trans = "identity"  # 连续色标的变换："identity"/"sqrt"/"log10" 等
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
    ggplot2::scale_fill_gradientn(colors = cols, limits = lims, oob = scales::squish, name = "", trans = scale_trans) +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::labs(title = title, x = "Longitude (°E)", y = "Latitude (°N)") +
    viz_theme_nature(base_size = theme_base_size, title_size = title_size)

  # 导出PNG
  ggplot2::ggsave(filename = paste0(out_base, ".png"), plot = p,
                  width = width_in, height = height_in, dpi = 1200, units = "in", bg = "white")
  # 导出SVG
  ggplot2::ggsave(filename = paste0(out_base, ".svg"), plot = p,
                  width = width_in, height = height_in, units = "in", bg = "white", device = "svg")

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
# 预览：基于流量累积构造的“河网热度”图（对数增强）
# ------------------------------
viz_preview_river_intensity <- function(flow_acc_path = "earthenvstreams_china/flow_acc.tif",
                                        out_base = "figures/00_china_env_variables/river_skeleton_preview",
                                        china_path = "earthenvstreams_china/china_boundary.shp") {
  # 中文注释：预览目的仅为观感评估。对flow_acc第2波段取 log1p 并线性归一化到[0,1] 后渲染。
  r <- terra::rast(flow_acc_path)[[2]]
  # 仅限河网像元
  r[r <= 0] <- NA
  # 对数增强
  r_log <- log1p(r)
  # 归一化
  mx <- suppressWarnings(as.numeric(terra::global(r_log, "max", na.rm = TRUE)[1,1]))
  if (is.finite(mx) && mx > 0) r_log <- r_log / mx
  viz_save_raster_map(r_log, out_base = out_base, title = "River Skeleton (from Flow Accumulation)",
                      palette = "viridis", q_limits = c(0.01, 0.99), china_path = china_path)
}

# ==============================================================================
# 结束
# ==============================================================================


