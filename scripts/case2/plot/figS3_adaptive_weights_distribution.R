# =============================================================================
# Fig S3: 自适应权重分布（箱线图 + 抖点）
# 展示三个维度权重 w_dag / w_ate / w_imp 在物种间的分布与区域差异
# 数据来源: results/case2/all_screening_v3.csv
# 输出: figures/case2/plot/figS3_adaptive_weights_distribution.{png,svg}
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

region_colors <- c(
    "AWT" = "#E64B35", "CAN" = "#4DBBD5", "NSW" = "#00A087",
    "NZ"  = "#3C5488", "SA"  = "#F39B7F", "SWI" = "#8491B4"
)
weight_colors <- c("w_dag" = "#E64B35", "w_ate" = "#4DBBD5", "w_imp" = "#00A087")
weight_labels <- c("w_dag" = "DAG Structure\n(w_dag)", "w_ate" = "ATE Effect\n(w_ate)",
                   "w_imp" = "RF Importance\n(w_imp)")

# ── 读取数据 ─────────────────────────────────────────────────────────────────
scr_file <- file.path(data_dir, "all_screening_v3.csv")
if (!file.exists(scr_file)) {
    stop("数据文件不存在: ", scr_file)
}
scr <- read.csv(scr_file, stringsAsFactors = FALSE)

# 每物种取唯一权重行（每物种只有一组权重，取第一行）
weights_per_sp <- scr %>%
    select(region, species, w_dag, w_ate, w_imp) %>%
    distinct() %>%
    filter(!is.na(w_dag))

# 转为长格式
weights_long <- weights_per_sp %>%
    pivot_longer(
        cols      = c(w_dag, w_ate, w_imp),
        names_to  = "weight_type",
        values_to = "weight_value"
    ) %>%
    mutate(
        weight_type = factor(weight_type, levels = c("w_dag", "w_ate", "w_imp")),
        weight_label = weight_labels[weight_type]
    )

# 总体统计
global_stats <- weights_long %>%
    group_by(weight_type) %>%
    summarise(
        median_w = median(weight_value, na.rm = TRUE),
        mean_w   = mean(weight_value, na.rm = TRUE),
        .groups  = "drop"
    )

# ── Panel A: 全局权重分布（箱线 + jitter）──────────────────────────────────
pA <- ggplot(weights_long,
             aes(x = weight_label, y = weight_value, fill = weight_type)) +
    geom_jitter(
        aes(color = weight_type),
        width = 0.2, alpha = 0.25, size = 0.8
    ) +
    geom_boxplot(
        outlier.shape = NA, alpha = 0.75, linewidth = 0.6,
        width = 0.5
    ) +
    geom_hline(
        yintercept = 1 / 3, linetype = "dashed",
        color = "grey50", linewidth = 0.6
    ) +
    annotate("text", x = 3.4, y = 1 / 3 + 0.01, label = "Equal weight (1/3)",
             hjust = 1, size = 2.9, color = "grey45", family = "Arial") +
    scale_fill_manual(values  = weight_colors) +
    scale_color_manual(values = weight_colors) +
    scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
    labs(
        title    = "A  Global Adaptive Weight Distribution",
        x        = NULL,
        y        = "Adaptive Weight",
        subtitle = sprintf("N = %d species, all regions combined", nrow(weights_per_sp))
    ) +
    theme_classic(base_family = "Arial", base_size = 11) +
    theme(
        legend.position  = "none",
        plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(size = 8.5, color = "grey40"),
        axis.text.x      = element_text(size = 10),
        axis.title.y     = element_text(size = 10),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35)
    )

# ── Panel B: 逐区域权重分布（facet） ─────────────────────────────────────────
pB <- ggplot(weights_long,
             aes(x = weight_type, y = weight_value, fill = region, color = region)) +
    geom_violin(scale = "width", alpha = 0.35, linewidth = 0.4, trim = TRUE) +
    geom_boxplot(
        outlier.shape = NA, alpha = 0.8, linewidth = 0.5,
        width = 0.35
    ) +
    geom_hline(yintercept = 1 / 3, linetype = "dashed",
               color = "grey55", linewidth = 0.55) +
    scale_fill_manual(values  = region_colors) +
    scale_color_manual(values = region_colors) +
    scale_y_continuous(limits = c(0, 1), labels = percent_format(accuracy = 1)) +
    scale_x_discrete(labels = c("w_dag" = "w_dag", "w_ate" = "w_ate", "w_imp" = "w_imp")) +
    facet_wrap(~region, ncol = 3) +
    labs(
        title    = "B  Adaptive Weight Distribution by Region",
        x        = "Weight Dimension",
        y        = "Adaptive Weight"
    ) +
    theme_classic(base_family = "Arial", base_size = 10) +
    theme(
        strip.background   = element_blank(),
        strip.text         = element_text(face = "bold", size = 10),
        legend.position    = "none",
        plot.title         = element_text(face = "bold", size = 12),
        axis.title         = element_text(size = 10),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35)
    )

# ── 合并两图（patchwork）─────────────────────────────────────────────────────
suppressPackageStartupMessages(library(patchwork))
p_combined <- pA / pB + plot_layout(heights = c(1, 1.4)) +
    plot_annotation(
        title   = "Adaptive Weight Distribution in CAST Screening",
        subtitle = "CAST automatically balances DAG structure, ATE effect size, and RF importance;\ndashed line marks uniform weight (1/3) for reference",
        theme = theme(
            plot.title    = element_text(face = "bold", size = 14, family = "Arial"),
            plot.subtitle = element_text(size = 9, color = "grey40", family = "Arial")
        )
    )

# ── 保存 ─────────────────────────────────────────────────────────────────────
out_prefix <- file.path(fig_dir, "figS3_adaptive_weights_distribution")
ggsave(paste0(out_prefix, ".png"), p_combined, width = 10, height = 12,
       dpi = 1200, bg = "white")
ggsave(paste0(out_prefix, ".svg"), p_combined, width = 10, height = 12, bg = "white")
cat("Fig S3 saved:", fig_dir, "\n")
