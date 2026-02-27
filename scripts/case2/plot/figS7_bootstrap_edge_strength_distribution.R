# =============================================================================
# Fig S7: Bootstrap边强度分布（跨区域直方图 + 密度 + 阈值标注）
# 展示Bootstrap重采样得到的DAG边出现频率（即edge strength）的分布特征
# 数据来源: results/case2/all_dag_edges_v3.csv
# 输出: figures/case2/plot/figS7_bootstrap_edge_strength_distribution.{png,svg}
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(scales)
    library(patchwork)
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

# ── 读取数据 ─────────────────────────────────────────────────────────────────
edges_file <- file.path(data_dir, "all_dag_edges_v3.csv")
if (!file.exists(edges_file)) stop("数据文件不存在: ", edges_file)
edges <- read.csv(edges_file, stringsAsFactors = FALSE) %>%
    filter(!is.na(strength))

# CAST中使用的strong_edge阈值（bootstrap强度 >= 0.5 被视为强边）
strong_threshold <- 0.5

# 标注强弱边
edges <- edges %>%
    mutate(edge_type = ifelse(strength >= strong_threshold, "Strong edge", "Weak edge"))

# ── Panel A: 全局边强度分布 ──────────────────────────────────────────────────
global_strong_pct <- mean(edges$strength >= strong_threshold) * 100

pA <- ggplot(edges, aes(x = strength)) +
    geom_histogram(
        aes(fill = edge_type), bins = 40,
        position = "identity", alpha = 0.75, color = "white", linewidth = 0.2
    ) +
    geom_vline(xintercept = strong_threshold, linetype = "dashed",
               color = "#E64B35", linewidth = 0.8) +
    annotate(
        "text", x = strong_threshold + 0.02, y = Inf,
        label = sprintf("Threshold = %.1f\n(%.1f%% strong)", strong_threshold, global_strong_pct),
        hjust = 0, vjust = 1.3, size = 3.2, color = "#E64B35", family = "Arial"
    ) +
    scale_fill_manual(
        values = c("Strong edge" = "#E64B35", "Weak edge" = "grey70"),
        name   = NULL
    ) +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    labs(
        title    = "A  Global Bootstrap Edge Strength Distribution",
        subtitle = sprintf(
            "All DAG edges across %d species × 6 regions | B = 100 bootstrap resamples",
            n_distinct(paste(edges$region, edges$species))
        ),
        x = "Bootstrap Strength (proportion of resamples)",
        y = "Count"
    ) +
    theme_classic(base_family = "Arial", base_size = 11) +
    theme(
        legend.position  = c(0.8, 0.85),
        legend.background = element_blank(),
        plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(size = 8.5, color = "grey40"),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35)
    )

# ── Panel B: 逐区域密度曲线 ──────────────────────────────────────────────────
# 各区域强边比例
region_stats <- edges %>%
    group_by(region) %>%
    summarise(
        pct_strong   = mean(strength >= strong_threshold) * 100,
        median_str   = median(strength),
        n_edges      = n(),
        n_species    = n_distinct(species),
        .groups      = "drop"
    ) %>%
    mutate(region_label = sprintf("%s\n(%d sp, %.0f%% strong)", region, n_species, pct_strong))

region_label_map <- setNames(region_stats$region_label, region_stats$region)
edges$region_label <- region_label_map[edges$region]

pB <- ggplot(edges, aes(x = strength, fill = region, color = region)) +
    geom_histogram(
        aes(y = after_stat(density)), bins = 30,
        alpha = 0.35, position = "identity", linewidth = 0.2
    ) +
    geom_density(linewidth = 0.8, alpha = 0) +
    geom_vline(xintercept = strong_threshold, linetype = "dashed",
               color = "grey40", linewidth = 0.55) +
    scale_fill_manual(values  = region_colors) +
    scale_color_manual(values = region_colors) +
    scale_x_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1)) +
    facet_wrap(~region_label, ncol = 3, scales = "free_y") +
    labs(
        title = "B  Bootstrap Edge Strength by Region",
        x     = "Bootstrap Strength",
        y     = "Density"
    ) +
    theme_classic(base_family = "Arial", base_size = 10) +
    theme(
        strip.background   = element_blank(),
        strip.text         = element_text(face = "bold", size = 9),
        legend.position    = "none",
        plot.title         = element_text(face = "bold", size = 12),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.35)
    )

# ── Panel C: 强边比例 vs DAG密度散点 ──────────────────────────────────────────
sp_edge_stats <- edges %>%
    group_by(region, species) %>%
    summarise(
        n_strong = sum(strength >= strong_threshold),
        n_total  = n(),
        pct_strong = n_strong / n_total,
        .groups  = "drop"
    )

# 尝试加入DAG密度信息
dag_info_file <- file.path(data_dir, "all_dag_info_v3.csv")
if (file.exists(dag_info_file)) {
    dag_info <- read.csv(dag_info_file, stringsAsFactors = FALSE)
    sp_edge_stats <- sp_edge_stats %>%
        left_join(dag_info %>% select(region, species, dag_density),
                  by = c("region", "species"))

    pC <- ggplot(sp_edge_stats,
                 aes(x = pct_strong, y = dag_density, color = region)) +
        geom_point(alpha = 0.6, size = 1.5) +
        geom_smooth(aes(group = 1), method = "loess", se = TRUE,
                    color = "grey30", linewidth = 0.8, fill = "grey85") +
        scale_color_manual(values = region_colors, name = "Region") +
        scale_x_continuous(labels = percent_format(accuracy = 1)) +
        labs(
            title    = "C  Strong Edge Rate vs DAG Density (per Species)",
            x        = "Proportion of Strong Edges (strength ≥ 0.5)",
            y        = "DAG Density",
            subtitle = "Each point = one species; line = LOESS trend"
        ) +
        theme_classic(base_family = "Arial", base_size = 10) +
        theme(
            legend.position    = "right",
            plot.title         = element_text(face = "bold", size = 11),
            plot.subtitle      = element_text(size = 8, color = "grey40"),
            panel.grid.major   = element_line(color = "grey92", linewidth = 0.35),
            legend.key.size    = unit(0.4, "cm"),
            legend.title       = element_text(face = "bold", size = 9)
        )
} else {
    pC <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = "all_dag_info_v3.csv not found",
                 family = "Arial", size = 5, color = "grey50") +
        theme_void()
}

# ── 合并 ─────────────────────────────────────────────────────────────────────
p_combined <- (pA | pC) / pB + plot_layout(heights = c(1, 1.5)) +
    plot_annotation(
        title    = "Bootstrap DAG Edge Strength Analysis",
        subtitle = "Edges with strength ≥ 0.5 are retained as 'strong edges' for CAST feature engineering",
        theme    = theme(
            plot.title    = element_text(face = "bold", size = 14, family = "Arial"),
            plot.subtitle = element_text(size = 9, color = "grey40", family = "Arial")
        )
    )

# ── 保存 ─────────────────────────────────────────────────────────────────────
out_prefix <- file.path(fig_dir, "figS7_bootstrap_edge_strength_distribution")
ggsave(paste0(out_prefix, ".png"), p_combined, width = 12, height = 12,
       dpi = 1200, bg = "white")
ggsave(paste0(out_prefix, ".svg"), p_combined, width = 12, height = 12, bg = "white")
cat("Fig S7 saved:", fig_dir, "\n")
