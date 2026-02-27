# =============================================================================
# Fig S4: 逐区域 CAST vs 每个基线详细对比（配对散点矩阵）
# 每个子图 X = baseline AUC，Y = CAST AUC，对角线 = 等效线
# 数据来源: results/case2/all_results_v3.csv
# 输出: figures/case2/plot/figS4_cast_vs_baselines_by_region.{png,svg}
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

# ── 读取并整理数据 ─────────────────────────────────────────────────────────────
res_file <- file.path(data_dir, "all_results_v3.csv")
if (!file.exists(res_file)) stop("数据文件不存在: ", res_file)
results <- read.csv(res_file, stringsAsFactors = FALSE)

# 提取 CAST 行
cast_auc <- results %>%
    filter(model == "CAST") %>%
    select(region, species, cast_auc = auc_mean)

# 基线模型列表
baselines <- c("RF", "MaxEnt", "BRT", "GAM", "MLP", "MLP_ATE")

# 合并CAST与各基线
paired <- results %>%
    filter(model %in% baselines) %>%
    select(region, species, model, base_auc = auc_mean) %>%
    inner_join(cast_auc, by = c("region", "species")) %>%
    filter(!is.na(base_auc), !is.na(cast_auc)) %>%
    mutate(delta = cast_auc - base_auc,
           cast_wins = delta > 0)

# 每 model×region 的胜率与mean delta
summary_stats <- paired %>%
    group_by(model, region) %>%
    summarise(
        n         = n(),
        win_rate  = mean(cast_wins, na.rm = TRUE),
        mean_delta = mean(delta, na.rm = TRUE),
        .groups   = "drop"
    ) %>%
    mutate(
        label = sprintf("Win: %d%%\nΔ=%.3f", round(win_rate * 100), mean_delta)
    )

# ── 绘图函数：单一基线的 6 区域散点 ──────────────────────────────────────────
make_baseline_plot <- function(df_base, model_name, stats_base) {
    # AUC 轴范围
    auc_range <- range(c(df_base$base_auc, df_base$cast_auc), na.rm = TRUE)
    lims <- c(max(0.4, auc_range[1] - 0.02), min(1, auc_range[2] + 0.02))

    ggplot(df_base, aes(x = base_auc, y = cast_auc, color = region)) +
        geom_abline(slope = 1, intercept = 0, color = "grey50",
                    linetype = "dashed", linewidth = 0.6) +
        geom_point(alpha = 0.65, size = 1.4) +
        # 每区域胜率文字标注
        geom_text(
            data = stats_base,
            aes(x = lims[1] + 0.01, y = lims[2] - 0.01, label = label,
                color = region),
            hjust = 0, vjust = 1, size = 2.3, family = "Arial",
            inherit.aes = FALSE
        ) +
        scale_color_manual(values = region_colors) +
        scale_x_continuous(limits = lims, labels = number_format(accuracy = 0.1)) +
        scale_y_continuous(limits = lims, labels = number_format(accuracy = 0.1)) +
        facet_wrap(~region, ncol = 3, scales = "fixed") +
        coord_fixed() +
        labs(
            title    = sprintf("CAST vs %s", model_name),
            x        = sprintf("%s AUC", model_name),
            y        = "CAST AUC",
            subtitle = "Dashed = parity line; points above = CAST wins"
        ) +
        theme_classic(base_family = "Arial", base_size = 9) +
        theme(
            strip.background = element_blank(),
            strip.text       = element_text(face = "bold", size = 9),
            legend.position  = "none",
            plot.title       = element_text(face = "bold", size = 11),
            plot.subtitle    = element_text(size = 7.5, color = "grey40"),
            axis.title       = element_text(size = 9)
        )
}

# 为每个基线生成图
plots <- lapply(baselines, function(m) {
    df_m    <- paired %>% filter(model == m)
    stats_m <- summary_stats %>% filter(model == m)
    make_baseline_plot(df_m, m, stats_m)
})
names(plots) <- baselines

# ── 拼接为 2×3 布局 ──────────────────────────────────────────────────────────
p_combined <- wrap_plots(plots, ncol = 2) +
    plot_annotation(
        title    = "CAST vs Baseline Models: Paired AUC Comparison by Region",
        subtitle = "Each panel: one baseline × six regions; points above diagonal = CAST outperforms",
        theme    = theme(
            plot.title    = element_text(face = "bold", size = 14, family = "Arial"),
            plot.subtitle = element_text(size = 9, color = "grey40", family = "Arial")
        )
    )

# ── 保存 ─────────────────────────────────────────────────────────────────────
out_prefix <- file.path(fig_dir, "figS4_cast_vs_baselines_by_region")
ggsave(paste0(out_prefix, ".png"), p_combined, width = 14, height = 18,
       dpi = 1200, bg = "white")
ggsave(paste0(out_prefix, ".svg"), p_combined, width = 14, height = 18, bg = "white")
cat("Fig S4 saved:", fig_dir, "\n")
