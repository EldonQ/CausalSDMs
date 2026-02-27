# =============================================================================
# Fig S1: DAG密度分布（跨区域直方图 + 核密度估计）
# 数据来源: results/case2/all_dag_info_v3.csv
# 输出: figures/case2/plot/figS1_dag_density_distribution.{png,svg}
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(scales)
})

# ── 路径配置 ──────────────────────────────────────────────────────────────────
script_dir <- tryCatch(
    dirname(rstudioapi::getSourceEditorContext()$path),
    error = function(e) getwd()
)
proj_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
data_dir  <- file.path(proj_root, "results", "case2")
fig_dir   <- file.path(proj_root, "figures", "case2", "plot")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ── 调色板（Nature风格，6区域）──────────────────────────────────────────────
region_colors <- c(
    "AWT" = "#E64B35", "CAN" = "#4DBBD5", "NSW" = "#00A087",
    "NZ"  = "#3C5488", "SA"  = "#F39B7F", "SWI" = "#8491B4"
)

# ── 读取数据 ─────────────────────────────────────────────────────────────────
dag_info_file <- file.path(data_dir, "all_dag_info_v3.csv")
if (!file.exists(dag_info_file)) {
    stop("数据文件不存在，请先运行 02_multi_species_experiment.R: ", dag_info_file)
}
dag_info <- read.csv(dag_info_file, stringsAsFactors = FALSE)

# 计算各区域统计量
region_stats <- dag_info %>%
    group_by(region) %>%
    summarise(
        n_species    = n(),
        median_dens  = median(dag_density, na.rm = TRUE),
        mean_dens    = mean(dag_density, na.rm = TRUE),
        q25          = quantile(dag_density, 0.25, na.rm = TRUE),
        q75          = quantile(dag_density, 0.75, na.rm = TRUE),
        .groups = "drop"
    )

# 全局中位数
global_median <- median(dag_info$dag_density, na.rm = TRUE)

# ── 绘图 ─────────────────────────────────────────────────────────────────────
p <- ggplot(dag_info, aes(x = dag_density, fill = region, color = region)) +
    # 各区域直方图（半透明叠加）
    geom_histogram(
        aes(y = after_stat(density)),
        bins = 25, alpha = 0.45, position = "identity", linewidth = 0.3
    ) +
    # 核密度曲线
    geom_density(linewidth = 0.9, alpha = 0) +
    # 全局中位数参考线
    geom_vline(
        xintercept = global_median,
        linetype = "dashed", color = "grey40", linewidth = 0.7
    ) +
    annotate(
        "text", x = global_median + 0.01, y = Inf,
        label = sprintf("Global median = %.2f", global_median),
        hjust = 0, vjust = 1.5, size = 3.2, color = "grey35", family = "Arial"
    ) +
    # 各区域中位数点标注
    geom_point(
        data = region_stats,
        aes(x = median_dens, y = 0, color = region),
        shape = 25, size = 3, fill = NA, stroke = 1.2,
        inherit.aes = FALSE
    ) +
    scale_fill_manual(values = region_colors)  +
    scale_color_manual(values = region_colors) +
    facet_wrap(~region, ncol = 3, scales = "free_y") +
    labs(
        title    = "DAG Density Distribution across Regions",
        subtitle = sprintf(
            "Bootstrap Hill-Climbing DAG (B = 100 resamples) | N = %d species × 6 regions",
            nrow(dag_info)
        ),
        x = "DAG Density (edges / max possible edges)",
        y = "Density",
        caption = "Triangles mark regional medians. Dashed line = global median."
    ) +
    theme_classic(base_family = "Arial", base_size = 11) +
    theme(
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 11),
        legend.position  = "none",
        plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(size = 9, color = "grey40"),
        plot.caption     = element_text(size = 8, color = "grey50"),
        axis.title       = element_text(size = 10),
        panel.grid.major = element_line(color = "grey93", linewidth = 0.4)
    )

# ── 保存 ─────────────────────────────────────────────────────────────────────
out_prefix <- file.path(fig_dir, "figS1_dag_density_distribution")
ggsave(paste0(out_prefix, ".png"), p, width = 10, height = 7,
       dpi = 1200, bg = "white")
ggsave(paste0(out_prefix, ".svg"), p, width = 10, height = 7, bg = "white")
cat("Fig S1 saved:", fig_dir, "\n")
