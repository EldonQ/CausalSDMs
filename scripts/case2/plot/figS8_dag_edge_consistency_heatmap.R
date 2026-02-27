# =============================================================================
# Fig S8: 所有区域DAG边一致性完整热力图
# 变量×变量矩阵，色深 = 该有向边在N个物种中出现的比例，逐区域 facet
# 数据来源: results/case2/all_dag_edges_v3.csv
# 输出: figures/case2/plot/figS8_dag_edge_consistency_heatmap.{png,svg}
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

# ── 读取数据 ─────────────────────────────────────────────────────────────────
edges_file <- file.path(data_dir, "all_dag_edges_v3.csv")
if (!file.exists(edges_file)) stop("数据文件不存在: ", edges_file)
edges <- read.csv(edges_file, stringsAsFactors = FALSE) %>%
    filter(!is.na(strength))

# 只保留强边（避免噪声边污染一致性统计）
strong_threshold <- 0.5
strong_edges <- edges %>% filter(strength >= strong_threshold)

# 每区域物种数
sp_count <- strong_edges %>%
    group_by(region) %>%
    summarise(n_sp = n_distinct(species), .groups = "drop")

# 每区域每个有向边出现比例
edge_freq <- strong_edges %>%
    group_by(region, from, to) %>%
    summarise(n_appear = n_distinct(species), .groups = "drop") %>%
    left_join(sp_count, by = "region") %>%
    mutate(freq = n_appear / n_sp)

# 获取所有出现过的变量（合并 from 和 to）
all_vars <- sort(unique(c(edge_freq$from, edge_freq$to)))
n_vars   <- length(all_vars)

# 全区域补全完整格子（未出现 = 0）
all_regions <- unique(edge_freq$region)
full_grid <- expand.grid(
    region = all_regions,
    from   = all_vars,
    to     = all_vars,
    stringsAsFactors = FALSE
) %>%
    filter(from != to)

edge_freq_full <- full_grid %>%
    left_join(edge_freq, by = c("region", "from", "to")) %>%
    replace_na(list(freq = 0, n_appear = 0)) %>%
    # 添加区域物种数标签
    left_join(sp_count, by = "region") %>%
    replace_na(list(n_sp = 0)) %>%
    mutate(
        region_label = paste0(region, "\n(n=", n_sp, " sp)"),
        from = factor(from, levels = all_vars),
        to   = factor(to,   levels = rev(all_vars))  # Y轴反转使对角线从左上到右下
    )

# ── 绘图 ─────────────────────────────────────────────────────────────────────
# 若变量数超过20，自适应缩小字体
axis_font_size <- max(5, min(9, 120 / n_vars))

p <- ggplot(edge_freq_full, aes(x = from, y = to, fill = freq)) +
    geom_tile(color = "white", linewidth = 0.25) +
    # 仅在高频边上打印数值（频率>=0.3）
    geom_text(
        data = edge_freq_full %>% filter(freq >= 0.3),
        aes(label = sprintf("%.0f%%", freq * 100)),
        size = axis_font_size * 0.32, color = "white",
        family = "Arial", fontface = "bold"
    ) +
    facet_wrap(~region_label, ncol = 3) +
    scale_fill_gradientn(
        colors = c("white", "#FEF0D9", "#FDD49E", "#FC8D59", "#E34A33", "#B30000"),
        values = scales::rescale(c(0, 0.1, 0.2, 0.4, 0.6, 1)),
        limits = c(0, 1),
        labels = percent_format(accuracy = 1),
        name   = "Edge Frequency\n(prop. of species)"
    ) +
    labs(
        title    = "DAG Edge Consistency Across Species within Each Region",
        subtitle = sprintf(
            "Strong edges only (bootstrap strength ≥ %.1f) | %d variables",
            strong_threshold, n_vars
        ),
        x = "From (parent variable)",
        y = "To (child variable)"
    ) +
    theme_classic(base_family = "Arial", base_size = 10) +
    theme(
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 10),
        axis.text.x      = element_text(angle = 45, hjust = 1,
                                        size = axis_font_size),
        axis.text.y      = element_text(size = axis_font_size),
        axis.title       = element_text(size = 10),
        legend.title     = element_text(size = 9),
        legend.text      = element_text(size = 8),
        plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(size = 9, color = "grey40"),
        panel.spacing    = unit(0.8, "lines")
    )

# ── 附图：全区域合并（汇总一致性热力图）────────────────────────────────────
global_edge_freq <- strong_edges %>%
    group_by(from, to) %>%
    summarise(
        n_appear = n_distinct(paste(region, species)),
        .groups  = "drop"
    ) %>%
    mutate(freq = n_appear / n_distinct(paste(strong_edges$region, strong_edges$species))) %>%
    right_join(
        expand.grid(from = all_vars, to = all_vars,
                    stringsAsFactors = FALSE) %>% filter(from != to),
        by = c("from", "to")
    ) %>%
    replace_na(list(freq = 0)) %>%
    mutate(
        from = factor(from, levels = all_vars),
        to   = factor(to,   levels = rev(all_vars))
    )

p_global <- ggplot(global_edge_freq, aes(x = from, y = to, fill = freq)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(
        data = global_edge_freq %>% filter(freq >= 0.2),
        aes(label = sprintf("%.0f%%", freq * 100)),
        size = axis_font_size * 0.32, color = "white",
        family = "Arial", fontface = "bold"
    ) +
    scale_fill_gradientn(
        colors = c("white", "#FEF0D9", "#FDD49E", "#FC8D59", "#E34A33", "#B30000"),
        values = scales::rescale(c(0, 0.05, 0.1, 0.2, 0.4, 1)),
        limits = c(0, 1),
        labels = percent_format(accuracy = 1),
        name   = "Global Edge\nFrequency"
    ) +
    labs(
        title = "Global DAG Edge Consistency (All Regions Combined)",
        x = "From", y = "To"
    ) +
    theme_classic(base_family = "Arial", base_size = 10) +
    theme(
        axis.text.x  = element_text(angle = 45, hjust = 1, size = axis_font_size),
        axis.text.y  = element_text(size = axis_font_size),
        axis.title   = element_text(size = 10),
        plot.title   = element_text(face = "bold", size = 12),
        legend.title = element_text(size = 9)
    )

# ── 组合输出 ──────────────────────────────────────────────────────────────────
suppressPackageStartupMessages(library(patchwork))
p_combined <- p / p_global + plot_layout(heights = c(2.5, 1)) +
    plot_annotation(
        title    = "Complete DAG Edge Consistency Heatmaps",
        subtitle = "Top: by region | Bottom: all regions combined",
        theme    = theme(
            plot.title    = element_text(face = "bold", size = 14, family = "Arial"),
            plot.subtitle = element_text(size = 9, color = "grey40", family = "Arial")
        )
    )

# ── 保存 ─────────────────────────────────────────────────────────────────────
out_prefix <- file.path(fig_dir, "figS8_dag_edge_consistency_heatmap")

# 根据变量数自适应图幅
fig_w <- max(10, n_vars * 0.55 + 4)
fig_h <- max(16, n_vars * 0.45 * 2.5 + 5)

ggsave(paste0(out_prefix, ".png"), p_combined, width = fig_w, height = fig_h,
       dpi = 1200, bg = "white", limitsize = FALSE)
ggsave(paste0(out_prefix, ".svg"), p_combined, width = fig_w, height = fig_h,
       bg = "white", limitsize = FALSE)
cat("Fig S8 saved:", fig_dir, "\n")
