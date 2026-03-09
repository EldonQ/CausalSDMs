################################################################################
# Fig 6: 空间 CATE 热力图 — 批量生成，单图单文件
#
# 核心叙事: CAST 通过因果森林 (grf) 揭示每个环境变量对物种分布的
#           空间异质性因果效应 — 同一变量在不同地区影响方向和强度不同
#
# 渲染方式:
#   IDW 空间插值 → 栅格化 → 中国边界裁剪 → geom_tile + geom_sf
#   参考 scripts/CASTplot/fig5_spatial_cate_map.R 的渲染风格
#
# 输出: 每张 CATE 图独立保存为 fig6_cate_{species}_{variable}.png/svg
#       保存至 figures/case2_eco/plot/cate_maps/
#
# ═══ 可配置参数（脚本头部修改）═══
#   TARGET_SPECIES : 指定物种 (NULL = 全部)
#   TARGET_VARS    : 指定变量 (NULL = 自动选取每物种所有已计算 CATE 变量)
#   IDW_RES        : IDW 插值网格分辨率 (度)
#
# 数据来源:
#   output/case2_eco/all_spatial_cate_v3.csv  (真实因果森林 CATE 预测)
#   output/case2_eco/all_ate_results_v3.csv
#   chinashp/china.shp
#
# 运行: setwd("E:/CausalSDMs")
#       source("scripts/EcoISEA3H/plot/fig6_spatial_cate_maps.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

# ══════════════════════════════════════════════════════════════════════════════
# ★ 可配置参数 — 按需修改
# ══════════════════════════════════════════════════════════════════════════════
TARGET_SPECIES <- NULL        # NULL = 全部; 或 c("Alces_alces", "Ovis_ammon")
TARGET_VARS    <- NULL        # NULL = 每物种所有已有 CATE 变量; 或 c("elevation", "bio19")
IDW_RES        <- 0.15        # IDW 插值网格分辨率 (度, 越小越精细但越慢)
IDW_NMAX       <- 15          # IDW 最大邻域点数
IDW_IDP        <- 2.0         # IDW 反距离幂次
BUFFER_DEG     <- 0.5         # 插值范围超出数据边界的缓冲 (度)

fig_dir <- "figures/case2_eco/plot/cate_maps"
tbl_dir <- "figures/case2_eco/tables"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 加载依赖
# ══════════════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({
    library(tidyverse)
    library(sf)
    library(terra)
    library(gstat)
})

# ══════════════════════════════════════════════════════════════════════════════
# 地图主题 — 参考 CASTplot/fig5 的浮空渲染风格
# ══════════════════════════════════════════════════════════════════════════════
theme_cate <- function() {
    theme_void(base_family = "sans") +
        theme(
            plot.title        = element_text(face = "bold", hjust = 0.5, size = 14,
                                             margin = margin(b = 6)),
            plot.subtitle     = element_text(hjust = 0.5, color = "grey40", size = 10,
                                             margin = margin(b = 12)),
            legend.position   = "right",
            legend.title      = element_text(face = "bold", size = 10),
            legend.text       = element_text(size = 8),
            legend.key.height = unit(1.8, "cm"),
            legend.key.width  = unit(0.4, "cm"),
            plot.background   = element_rect(fill = "white", color = NA),
            panel.background  = element_rect(fill = "white", color = NA),
            plot.margin       = margin(10, 10, 10, 10)
        )
}

# ══════════════════════════════════════════════════════════════════════════════
# 变量英文显示名称
# ══════════════════════════════════════════════════════════════════════════════
var_labels <- c(
    "aridityindexthornthwaite" = "Aridity Index",
    "bio02"                    = "Diurnal Range (Bio02)",
    "bio15"                    = "Precip. Seasonality (Bio15)",
    "bio19"                    = "Precip. Coldest Qtr (Bio19)",
    "elevation"                = "Elevation",
    "etccdi_cwd"               = "Consecutive Wet Days",
    "landcover_igbp"           = "Land Cover (IGBP)",
    "maxtempcoldest"           = "Tmax Coldest Month",
    "nontree"                  = "Non-tree Vegetation",
    "topowet"                  = "Topographic Wetness",
    "tri"                      = "Terrain Ruggedness"
)

get_var_label <- function(x) {
    out <- var_labels[x]
    out[is.na(out)] <- x[is.na(out)]
    unname(out)
}

fmt_sp_display <- function(x) gsub("_", " ", x)

# ══════════════════════════════════════════════════════════════════════════════
# 读取数据
# ══════════════════════════════════════════════════════════════════════════════
cat("读取 CATE 数据...\n")
cate_all <- read.csv("output/case2_eco/all_spatial_cate_v3.csv",
                     stringsAsFactors = FALSE)

ate_all <- read.csv("output/case2_eco/all_ate_results_v3.csv",
                    stringsAsFactors = FALSE) %>%
    mutate(coef = as.numeric(coef),
           significant = as.logical(significant))

cat("读取中国边界...\n")
china_sf <- st_read("chinashp/china.shp", quiet = TRUE)

cat(sprintf("CATE 数据: %s 条, %d 物种, %d 变量\n",
            format(nrow(cate_all), big.mark = ","),
            n_distinct(cate_all$species),
            n_distinct(cate_all$variable)))

# ══════════════════════════════════════════════════════════════════════════════
# 确定待绘制的 (物种, 变量) 组合
# ══════════════════════════════════════════════════════════════════════════════
available_combos <- cate_all %>%
    distinct(species, variable)

if (!is.null(TARGET_SPECIES)) {
    available_combos <- available_combos %>%
        filter(species %in% TARGET_SPECIES)
}
if (!is.null(TARGET_VARS)) {
    available_combos <- available_combos %>%
        filter(variable %in% TARGET_VARS)
}

plot_tasks <- available_combos %>%
    left_join(
        ate_all %>% select(species, variable, coef, significant),
        by = c("species", "variable")
    ) %>%
    replace_na(list(coef = 0, significant = FALSE)) %>%
    arrange(species, desc(significant), desc(abs(coef)))

cat(sprintf("待绘制: %d 张 CATE 热力图 (%d 物种)\n",
            nrow(plot_tasks), n_distinct(plot_tasks$species)))

# 保存任务清单
write.csv(plot_tasks, file.path(tbl_dir, "fig6_cate_plot_tasks.csv"),
          row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 核心渲染函数: 单张 CATE 热力图
#   1. 提取 CATE 网格点
#   2. IDW 空间插值 → 平滑热力场
#   3. 栅格化 → 中国边界裁剪
#   4. geom_tile 渲染 + geom_sf 叠加边界
# ══════════════════════════════════════════════════════════════════════════════
render_cate_map <- function(cate_df, sp_name, var_name, ate_coef, china_boundary) {

    # ── 提取该物种×变量的 CATE 网格点 ────────────────────────────────────────
    pts <- cate_df %>%
        filter(species == sp_name, variable == var_name) %>%
        select(lon, lat, cate) %>%
        drop_na()

    if (nrow(pts) < 30) {
        cat(sprintf("  ⚠ 数据不足 (%d 点), 跳过\n", nrow(pts)))
        return(NULL)
    }

    # ── IDW 空间插值 ─────────────────────────────────────────────────────────
    pts_sf <- st_as_sf(pts, coords = c("lon", "lat"), crs = 4326)

    # 构建插值目标网格 (覆盖数据范围 + 缓冲)
    lon_rng <- range(pts$lon) + c(-BUFFER_DEG, BUFFER_DEG)
    lat_rng <- range(pts$lat) + c(-BUFFER_DEG, BUFFER_DEG)
    grid_lons <- seq(lon_rng[1], lon_rng[2], by = IDW_RES)
    grid_lats <- seq(lat_rng[1], lat_rng[2], by = IDW_RES)
    grid_df   <- expand.grid(lon = grid_lons, lat = grid_lats)
    grid_sf   <- st_as_sf(grid_df, coords = c("lon", "lat"), crs = 4326)

    # IDW 插值
    idw_out <- gstat::idw(
        cate ~ 1,
        locations  = pts_sf,
        newdata    = grid_sf,
        idp        = IDW_IDP,
        nmax       = IDW_NMAX,
        debug.level = 0
    )
    grid_df$CATE <- idw_out$var1.pred

    # ── 栅格化 + 中国边界裁剪 ────────────────────────────────────────────────
    r_interp <- rast(grid_df, type = "xyz", crs = "EPSG:4326")
    r_masked <- tryCatch(
        mask(r_interp, vect(china_boundary)),
        error = function(e) {
            cat(sprintf("  ⚠ mask 失败 (%s), 使用未裁剪数据\n", e$message))
            r_interp
        }
    )

    masked_df <- as.data.frame(r_masked, xy = TRUE, na.rm = TRUE)
    if (ncol(masked_df) >= 3) {
        names(masked_df)[3] <- "CATE"
    } else {
        cat("  ⚠ 裁剪后无数据, 跳过\n")
        return(NULL)
    }

    if (nrow(masked_df) < 10) {
        cat("  ⚠ 裁剪后数据不足, 跳过\n")
        return(NULL)
    }

    # ── 对称发散色标 (以 0 为中心) ───────────────────────────────────────────
    cate_abs_lim <- quantile(abs(masked_df$CATE), 0.98, na.rm = TRUE)
    if (cate_abs_lim < 1e-8) cate_abs_lim <- max(abs(masked_df$CATE)) + 1e-6

    # ── 渲染: geom_tile + geom_sf (参考 CASTplot/fig5 风格) ──────────────────
    p <- ggplot() +
        geom_tile(data = masked_df, aes(x = x, y = y, fill = CATE)) +
        geom_sf(data = china_boundary, fill = NA, color = "grey30",
                linewidth = 0.3) +
        scale_fill_gradient2(
            low      = "#2166AC",
            mid      = "#F7F7F7",
            high     = "#B2182B",
            midpoint = 0,
            limits   = c(-cate_abs_lim, cate_abs_lim),
            oob      = scales::squish,
            name     = "CATE"
        ) +
        coord_sf(expand = FALSE) +
        labs(
            title    = sprintf("Spatial CATE: %s", get_var_label(var_name)),
            subtitle = sprintf(
                "%s  |  ATE = %.4f  |  n = %s grid cells",
                fmt_sp_display(sp_name), ate_coef,
                format(nrow(pts), big.mark = ",")
            )
        ) +
        theme_cate()

    return(p)
}

# ══════════════════════════════════════════════════════════════════════════════
# 批量绘制
# ══════════════════════════════════════════════════════════════════════════════
cat(sprintf("\n%s\n", paste0(rep("=", 60), collapse = "")))
cat("开始批量绘制 CATE 热力图\n")
cat(sprintf("%s\n\n", paste0(rep("=", 60), collapse = "")))

success_count <- 0
t_start <- Sys.time()

for (i in seq_len(nrow(plot_tasks))) {
    sp  <- plot_tasks$species[i]
    var <- plot_tasks$variable[i]
    ate_val <- plot_tasks$coef[i]

    cat(sprintf("[%d/%d] %s x %s (ATE = %.4f)\n",
                i, nrow(plot_tasks), sp, var, ate_val))

    p <- tryCatch(
        render_cate_map(cate_all, sp, var, ate_val, china_sf),
        error = function(e) {
            cat(sprintf("  ✗ 错误: %s\n", e$message))
            return(NULL)
        }
    )

    if (!is.null(p)) {
        fname <- sprintf("fig6_cate_%s_%s", sp, var)

        ggsave(file.path(fig_dir, paste0(fname, ".png")),
               p, width = 10, height = 8, dpi = 1200, bg = "white")
        tryCatch(
            ggsave(file.path(fig_dir, paste0(fname, ".svg")),
                   p, width = 10, height = 8, bg = "white"),
            error = function(e) NULL
        )

        success_count <- success_count + 1
        cat(sprintf("  ✓ 已保存\n"))
    }
}

elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "mins"))

cat(sprintf("\n%s\n", paste0(rep("=", 60), collapse = "")))
cat(sprintf("  Fig 6 完成: %d / %d 张 CATE 热力图已保存\n",
            success_count, nrow(plot_tasks)))
cat(sprintf("  输出目录: %s\n", fig_dir))
cat(sprintf("  耗时: %.1f 分钟\n", elapsed))
cat(sprintf("%s\n", paste0(rep("=", 60), collapse = "")))
