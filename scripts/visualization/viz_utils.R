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
required_pkgs <- c("terra", "sf", "ggplot2", "viridis", "RColorBrewer", "svglite", "sysfonts", "showtext", "rlang", "ggnewscale", "scales")
for (pkg in required_pkgs) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# ==============================================================================
# ★★★ 河网可视化增强参数配置 (River Network Visualization Enhancement) ★★★
# ==============================================================================
# 中文说明：以下参数用于控制地图输出时河网的渲染效果。
#           核心策略：通过形态学膨胀（Thickening）让河网线条变粗，从而在地图上清晰可见，
#           但**保留原始预测值的颜色**，不使用任何覆盖色。
# English: Parameters below control river network rendering.
#          Core strategy: Use morphological thickening to make river lines visible,
#          while PRESERVING the original prediction colors (no color overlay).
# ==============================================================================

# ------------------------------
# 河网增强开关 (River Enhancement Switch)
# ------------------------------
# 是否启用河网加粗增强？
# TRUE = 对输入栅格进行膨胀处理，使河网变粗可见
# FALSE = 保持原始栅格（1像元宽，可能看不清）
RIVER_ENHANCE_ENABLED <- TRUE

# ------------------------------
# 栅格加粗参数 (Raster Thickening Parameters)
# ------------------------------
# 加粗因子 (Thickening Factor)
# 控制河网线条的视觉宽度。
# 1 = 不加粗 (原始宽度)
# 3 = 3x3 窗口膨胀 (线条变粗)
# 5 = 5x5 窗口膨胀 (线条更粗)
# 推荐值: 3 或 5
RIVER_THICKEN_FACTOR <- 5

# 加粗阈值分位数 (Thickening Quantile Threshold)
# 控制仅对多少比例的“大河”进行加粗，避免小支流过粗导致画面杂乱。
# 0.95 = 仅对流量累积前 5% 的大河进行加粗，其余保持细线
# 0.00 = 对所有河流进行加粗
# 推荐值: 0.95 (仅突出干流) ~ 0.99 (仅突出特大河流)
RIVER_THICKEN_QUANTILE <- 0.95

# ------------------------------
# 矢量辅助参数 (Vector Helper Parameters)
# ------------------------------
# (已移除：不再使用矢量叠加)

# flow_acc.tif 路径 (Flow Accumulation Raster Path)
# 用于计算河流等级以进行选择性加粗
FLOW_ACC_PATH <- "earthenvstreams_china/flow_acc.tif"

# ==============================================================================
# ★★★ 河网参数配置结束 ★★★
# ==============================================================================

# ------------------------------
# 字体与DPI设置：统一为Arial + 2400dpi
# ------------------------------
viz_ensure_arial <- function() {
  # 中文注释：在Windows环境注册Arial，启用showtext以保证嵌入与一致渲染
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
  if (inherits(china, "try-error")) {
    return(NULL)
  }
  if (!is.null(crs_target)) {
    china <- sf::st_transform(china, crs = crs_target)
  }
  return(china)
}

# ------------------------------
# 将Raster/SpatRaster统一为SpatRaster
# ------------------------------
viz_as_spatraster <- function(r) {
  if (inherits(r, "SpatRaster")) {
    return(r)
  }
  return(terra::rast(r))
}

# ------------------------------
# 根据分位数裁剪色标范围（避免极端值主导）
# ------------------------------
viz_quantile_limits <- function(values, q = c(0.01, 0.99)) {
  v <- values[!is.na(values)]
  if (length(v) == 0) {
    return(c(NA, NA))
  }
  lo <- as.numeric(stats::quantile(v, q[1], na.rm = TRUE))
  hi <- as.numeric(stats::quantile(v, q[2], na.rm = TRUE))
  if (!is.finite(lo) || !is.finite(hi)) {
    return(range(v, na.rm = TRUE))
  }
  if (lo >= hi) {
    return(range(v, na.rm = TRUE))
  }
  c(lo, hi)
}

# ------------------------------
# 保存地图：PNG(2400dpi)+SVG，英文标注，Arial
# 支持河网视觉加粗（保留原始颜色）
# ------------------------------
viz_save_raster_map <- function(
    r, # RasterLayer/SpatRaster：待绘制的单层栅格（已按需要掩膜）
    out_base, # 字符串：输出文件基名（不含扩展名）
    title = "", # 字符串：图标题（英文）
    palette = "magma", # viridis 选项："magma"/"viridis"/"plasma"/"inferno" 等
    q_limits = c(0.01, 0.99), # 分位数裁剪区间
    china_path = "earthenvstreams_china/china_boundary.shp",
    width_in = 8, height_in = 6,
    theme_base_size = 9, title_size = 9,
    scale_trans = "identity", # 连续色标的变换："identity"/"sqrt"/"log10" 等
    # ========== 河网增强参数（可覆盖全局配置） ==========
    river_enhance_enabled = NULL, # TRUE/FALSE
    river_thicken_factor = NULL, # 加粗因子 (1, 3, 5...)
    river_thicken_quantile = NULL # 加粗分位数 (0.0~1.0)
    ) {
  # 中文注释：该函数将输入栅格转换为数据框后使用ggplot2绘制；
  #           若启用河网增强，则对输入栅格进行形态学膨胀（focal max），
  #           使河网线条变粗，但保留原始预测值的颜色。

  # ===== 参数回退到全局配置 =====
  if (is.null(river_enhance_enabled)) river_enhance_enabled <- RIVER_ENHANCE_ENABLED
  if (is.null(river_thicken_factor)) river_thicken_factor <- RIVER_THICKEN_FACTOR
  if (is.null(river_thicken_quantile)) river_thicken_quantile <- RIVER_THICKEN_QUANTILE

  r_spat <- viz_as_spatraster(r)
  china <- viz_read_china(china_path, crs_target = terra::crs(r_spat))

  # ===== 河网视觉加粗处理 (Selective Thickening) =====
  # 逻辑：如果启用增强且加粗因子 > 1
  if (river_enhance_enabled && river_thicken_factor > 1) {
    # 1. 全局膨胀（备用）
    w <- matrix(1, nrow = river_thicken_factor, ncol = river_thicken_factor)
    r_thick <- terra::focal(r_spat, w = w, fun = max, na.rm = TRUE)

    # 2. 选择性加粗：仅对大河应用膨胀结果
    if (river_thicken_quantile > 0 && file.exists(FLOW_ACC_PATH)) {
      # 加载流量累积栅格
      fa <- terra::rast(FLOW_ACC_PATH)[[2]]
      # 确保与预测栅格对齐（裁剪/重采样）
      # 注意：假设 flow_acc 与预测栅格原本就是对齐的（通常是），若不对齐需 project
      if (!terra::compareGeom(r_spat, fa, stopOnError = FALSE)) {
        fa <- terra::project(fa, r_spat)
        fa <- terra::crop(fa, r_spat, mask = TRUE)
      } else {
        fa <- terra::crop(fa, r_spat, mask = TRUE)
      }

      # 计算阈值
      fa_vals <- terra::values(fa, mat = FALSE)
      fa_vals <- fa_vals[fa_vals > 0 & !is.na(fa_vals)]
      if (length(fa_vals) > 0) {
        thresh <- quantile(fa_vals, probs = river_thicken_quantile, na.rm = TRUE)

        # 生成“大河掩膜”：流量 > 阈值
        # 同样需要对大河掩膜进行膨胀，以确定“加粗区域”
        mask_major <- (fa > thresh)
        mask_major[mask_major == 0] <- NA
        mask_major_thick <- terra::focal(mask_major, w = w, fun = max, na.rm = TRUE)

        # 组合：在 mask_major_thick 区域使用 r_thick，其他区域保持 r_spat
        # 找出需要加粗的像元索引
        idx_thick <- !is.na(terra::values(mask_major_thick))
        # 替换值
        r_combined <- r_spat
        r_combined[idx_thick] <- r_thick[idx_thick]
        r_spat <- r_combined
      } else {
        # 无法计算阈值（如无流量数据），回退到全量加粗
        r_spat <- r_thick
      }
    } else {
      # 不进行选择性加粗（全量加粗）
      r_spat <- r_thick
    }
  }

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

  # ===== 绘图 =====
  p <- ggplot2::ggplot() +
    ggplot2::geom_raster(data = df, ggplot2::aes_string(x = "x", y = "y", fill = "value")) +
    {
      if (!is.null(china)) ggplot2::geom_sf(data = china, fill = NA, color = "black", linewidth = 0.3, alpha = 0.9) else NULL
    } +
    ggplot2::scale_fill_gradientn(colors = cols, limits = lims, oob = scales::squish, name = "", trans = scale_trans) +
    ggplot2::coord_sf(expand = FALSE) +
    ggplot2::labs(title = title, x = "Longitude (°E)", y = "Latitude (°N)") +
    viz_theme_nature(base_size = theme_base_size, title_size = title_size)

  # 导出PNG（透明背景 + 2400dpi）
  ggplot2::ggsave(
    filename = paste0(out_base, ".png"), plot = p,
    width = width_in, height = height_in, dpi = 2400,
    units = "in", bg = "transparent"
  )
  # 导出SVG（透明背景）
  ggplot2::ggsave(
    filename = paste0(out_base, ".svg"), plot = p,
    width = width_in, height = height_in, units = "in",
    bg = "transparent", device = "svg"
  )

  invisible(TRUE)
}
