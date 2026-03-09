################################################################################
# Fig 4: CAST 性能鲁棒性 — 在所有物种上保持竞争力
#
# 核心叙事: CAST 不追求 AUC 的大幅超越，而是在提供因果可解释性的同时，
#           稳定地保持与 SOTA 模型可比的预测性能
#
# 三个独立子图（每图单独保存 1200 dpi PNG + SVG）:
#   (a) 全物种多模型 AUC 点阵图 (Cleveland dot plot)
#   (b) 模型排名分布 (stacked bar)
#   (c) 因果结构消融实验 (MLP → MLP_ATE → CAST)
#
# 数据来源:
#   output/case2_eco/all_results_v3.csv
#
# 运行: setwd("E:/CausalSDMs")
#       source("scripts/EcoISEA3H/plot/fig4_performance_robustness.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

fig_dir <- "figures/case2_eco/plot"
tbl_dir <- "figures/case2_eco/tables"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
    library(tidyverse)
})

# ══════════════════════════════════════════════════════════════════════════════
# 全局主题 · Nature 风格
# ══════════════════════════════════════════════════════════════════════════════
theme_pub <- function(base_size = 10) {
    theme_minimal(base_size = base_size, base_family = "sans") +
        theme(
            panel.grid.minor = element_blank(),
            axis.title       = element_text(face = "bold", size = 10),
            axis.text        = element_text(size = 8, color = "grey20"),
            plot.title       = element_text(face = "bold", size = 11, hjust = 0,
                                            margin = margin(b = 3)),
            plot.subtitle    = element_text(size = 8.5, color = "grey40", hjust = 0,
                                            margin = margin(b = 8)),
            legend.title     = element_text(face = "bold", size = 9),
            legend.text      = element_text(size = 8),
            plot.margin      = margin(12, 12, 8, 8)
        )
}

# ══════════════════════════════════════════════════════════════════════════════
# Nature 级模型配色
# ══════════════════════════════════════════════════════════════════════════════
model_colors <- c(
    "CAST"    = "#E64B35",
    "MLP_ATE" = "#F39B7F",
    "MLP"     = "#91D1C2",
    "RF"      = "#4DBBD5",
    "BRT"     = "#3C5488",
    "Maxent"  = "#B09C85"
)
models_ordered <- c("CAST", "MLP_ATE", "MLP", "RF", "BRT", "Maxent")

# 物种名格式化
fmt_sp <- function(x) {
    sapply(x, function(s) {
        p <- strsplit(s, "_")[[1]]
        if (length(p) >= 2) paste0(substr(p[1], 1, 1), ". ", p[2]) else s
    }, USE.NAMES = FALSE)
}

save_fig <- function(plt, name, w, h) {
    ggsave(file.path(fig_dir, paste0(name, ".png")),
           plt, width = w, height = h, dpi = 1200, bg = "white")
    tryCatch(
        ggsave(file.path(fig_dir, paste0(name, ".svg")),
               plt, width = w, height = h, bg = "white"),
        error = function(e) cat(sprintf("  [SVG 跳过: %s]\n", e$message))
    )
    cat(sprintf("  ✓ %s 已保存\n", name))
}

# ══════════════════════════════════════════════════════════════════════════════
# 读取数据
# ══════════════════════════════════════════════════════════════════════════════
d <- read.csv("output/case2_eco/all_results_v3.csv", stringsAsFactors = FALSE) %>%
    filter(!is.na(auc_mean), model %in% models_ordered) %>%
    mutate(model = factor(model, levels = models_ordered))

n_sp <- n_distinct(d$species)
cat(sprintf("数据加载完成: %d 物种, %d 模型\n", n_sp, length(models_ordered)))

# ══════════════════════════════════════════════════════════════════════════════
# (a) 全物种多模型 AUC 点阵图 (Cleveland dot plot)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── 绘制 (a) AUC 点阵图 ──\n")

# 物种按 CAST 的 AUC 排序
cast_order <- d %>%
    filter(model == "CAST") %>%
    arrange(auc_mean) %>%
    pull(species)

d_dot <- d %>%
    mutate(species = factor(species, levels = cast_order))

# 保存完整性能表
perf_wide <- d %>%
    select(species, model, auc_mean) %>%
    pivot_wider(names_from = model, values_from = auc_mean)
write.csv(perf_wide, file.path(tbl_dir, "fig4a_auc_all_species.csv"),
          row.names = FALSE)

# 汇总统计
model_summary <- d %>%
    group_by(model) %>%
    summarise(
        mean_auc = mean(auc_mean, na.rm = TRUE),
        sd_auc   = sd(auc_mean, na.rm = TRUE),
        min_auc  = min(auc_mean, na.rm = TRUE),
        .groups  = "drop"
    )
cat("模型汇总:\n")
print(model_summary)

pa <- ggplot(d_dot, aes(x = auc_mean, y = species, color = model, shape = model)) +
    geom_point(size = 2.2, alpha = 0.82) +
    scale_color_manual(values = model_colors, name = "Model") +
    scale_shape_manual(values = c(16, 17, 15, 18, 8, 4), name = "Model") +
    scale_y_discrete(labels = fmt_sp) +
    scale_x_continuous(
        limits = c(min(d$auc_mean, na.rm = TRUE) - 0.002, 1.001),
        breaks = seq(0.96, 1.00, by = 0.01)
    ) +
    labs(
        title = "(a) AUC across all species and models",
        subtitle = sprintf(
            "%d species — CAST consistently achieves competitive performance (mean AUC = %.4f)",
            n_sp, model_summary$mean_auc[model_summary$model == "CAST"]
        ),
        x = "AUC", y = ""
    ) +
    theme_pub() +
    theme(
        axis.text.y  = element_text(size = 5.5, face = "italic"),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
        legend.position = "right"
    )

save_fig(pa, "fig4a_auc_dotplot", w = 9, h = 8)

# ══════════════════════════════════════════════════════════════════════════════
# (b) 模型排名分布
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── 绘制 (b) 模型排名分布 ──\n")

rank_data <- d %>%
    group_by(species) %>%
    mutate(rank = rank(-auc_mean, ties.method = "min")) %>%
    ungroup()

rank_summary <- rank_data %>%
    count(model, rank) %>%
    group_by(model) %>%
    mutate(pct = n / sum(n) * 100) %>%
    ungroup()

# 保存排名表
write.csv(rank_summary, file.path(tbl_dir, "fig4b_rank_distribution.csv"),
          row.names = FALSE)

# 排名配色: 绿色(1st) → 红色(6th)
n_models <- length(models_ordered)
rank_colors <- c(
    "1" = "#1A9850", "2" = "#91CF60", "3" = "#D9EF8B",
    "4" = "#FEE08B", "5" = "#FC8D59", "6" = "#D73027"
)

# 模型排序: 按平均排名（好→差）
model_mean_rank <- rank_data %>%
    group_by(model) %>%
    summarise(mr = mean(rank), .groups = "drop") %>%
    arrange(mr)

rank_plot <- rank_summary %>%
    mutate(
        rank  = factor(rank),
        model = factor(model, levels = rev(model_mean_rank$model))
    )

pb <- ggplot(rank_plot, aes(x = model, y = pct, fill = rank)) +
    geom_col(position = "stack", width = 0.7) +
    geom_text(
        aes(label = ifelse(pct >= 8, sprintf("%.0f%%", pct), "")),
        position = position_stack(vjust = 0.5),
        size = 2.5, color = "grey20"
    ) +
    scale_fill_manual(values = rank_colors, name = "Rank") +
    coord_flip() +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(
        title = "(b) Model rank distribution across species",
        subtitle = sprintf(
            "Rank 1 = best AUC per species (n = %d species) — sorted by mean rank",
            n_sp
        ),
        x = "", y = "Proportion (%)"
    ) +
    theme_pub() +
    theme(panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3))

save_fig(pb, "fig4b_rank_distribution", w = 8, h = 4.5)

# ══════════════════════════════════════════════════════════════════════════════
# (c) 因果结构消融: MLP → MLP_ATE → CAST
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── 绘制 (c) 因果结构消融 ──\n")

ablation_colors <- c(
    "MLP"     = "#91D1C2",
    "MLP_ATE" = "#F39B7F",
    "CAST"    = "#E64B35"
)

d_abl <- d %>%
    filter(model %in% c("MLP", "MLP_ATE", "CAST")) %>%
    mutate(model = factor(model, levels = c("MLP", "MLP_ATE", "CAST")))

abl_stats <- d_abl %>%
    group_by(model) %>%
    summarise(mean_auc = mean(auc_mean), .groups = "drop")

# Wilcoxon 配对检验
wide_abl <- d %>%
    filter(model %in% c("MLP", "MLP_ATE", "CAST")) %>%
    select(species, model, auc_mean) %>%
    pivot_wider(names_from = model, values_from = auc_mean, values_fn = max) %>%
    filter(!is.na(CAST), !is.na(MLP), !is.na(MLP_ATE))

p_val1 <- tryCatch(
    wilcox.test(wide_abl$MLP_ATE, wide_abl$MLP, paired = TRUE)$p.value,
    error = function(e) NA_real_
)
p_val2 <- tryCatch(
    wilcox.test(wide_abl$CAST, wide_abl$MLP_ATE, paired = TRUE)$p.value,
    error = function(e) NA_real_
)

fmt_p <- function(p) {
    if (is.na(p)) return("")
    if (p < 0.001) return("***")
    if (p < 0.01)  return("**")
    if (p < 0.05)  return("*")
    return("ns")
}

cat(sprintf("  Wilcoxon: MLP→MLP_ATE p=%.4f (%s), MLP_ATE→CAST p=%.4f (%s)\n",
            p_val1, fmt_p(p_val1), p_val2, fmt_p(p_val2)))

# 保存消融结果
write.csv(abl_stats, file.path(tbl_dir, "fig4c_ablation_stats.csv"),
          row.names = FALSE)

y_max <- max(d_abl$auc_mean, na.rm = TRUE)

pc <- ggplot(d_abl, aes(x = model, y = auc_mean, fill = model)) +
    geom_violin(
        scale = "width", alpha = 0.5, trim = TRUE,
        linewidth = 0.3, color = "grey40"
    ) +
    geom_boxplot(
        width = 0.12, fill = "white", alpha = 0.85,
        outlier.size = 0.7, outlier.alpha = 0.5
    ) +
    stat_summary(
        fun = mean, geom = "point", size = 3.5,
        shape = 18, color = "black"
    ) +
    geom_text(
        data = abl_stats,
        aes(x = model, y = mean_auc,
            label = sprintf("%.4f", mean_auc)),
        vjust = -2.5, size = 3, fontface = "bold",
        color = "black", inherit.aes = FALSE
    ) +
    # 显著性标注
    annotate("segment", x = 1, xend = 2,
             y = y_max + 0.004, yend = y_max + 0.004,
             linewidth = 0.35, color = "grey30") +
    annotate("text", x = 1.5, y = y_max + 0.006,
             label = fmt_p(p_val1), size = 3.2, color = "grey30") +
    annotate("segment", x = 2, xend = 3,
             y = y_max + 0.009, yend = y_max + 0.009,
             linewidth = 0.35, color = "grey30") +
    annotate("text", x = 2.5, y = y_max + 0.011,
             label = fmt_p(p_val2), size = 3.2, color = "grey30") +
    scale_fill_manual(values = ablation_colors, guide = "none") +
    scale_x_discrete(labels = c(
        "MLP"     = "MLP\n(base)",
        "MLP_ATE" = "MLP_ATE\n(+ATE weighting)",
        "CAST"    = "CAST\n(+DAG interactions)"
    )) +
    labs(
        title = "(c) Causal structure ablation",
        subtitle = "Stepwise addition of causal information: base → +ATE weighting → +DAG interaction features",
        x = "", y = "AUC"
    ) +
    theme_pub() +
    theme(panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3))

save_fig(pc, "fig4c_causal_ablation", w = 7, h = 5)

cat("\n════════════════════════════════════════\n")
cat("  Fig 4 完成: 3 个子图 + 3 个数据表已保存\n")
cat("════════════════════════════════════════\n")
