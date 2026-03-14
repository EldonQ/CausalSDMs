################################################################################
# Fig 4: CAST 性能鲁棒性 — 在所有物种上保持竞争力
#
# 核心叙事: CAST 不追求 AUC 的大幅超越，而是在提供因果可解释性的同时，
#           稳定地保持与 SOTA 模型可比的预测性能
#
# 两个独立子图（每图单独保存 1200 dpi PNG + SVG）:
#   (a) 全物种多模型 AUC 点阵图 (Cleveland dot plot)
#   (b) 模型排名分布 (stacked bar)
#
# 数据来源 (Plant 案例, 03_run_Plant_multi_species.R 输出):
#   output/case4_plant/all_results_plant.csv
#
# 运行: setwd("E:/CausalSDMs")
#       source("scripts/Plant/plot/fig4_performance_robustness.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

fig_dir <- "figures/case4_plant/plot"
tbl_dir <- "figures/case4_plant/tables"
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
            axis.title = element_text(face = "bold", size = 10),
            axis.text = element_text(size = 8, color = "grey20"),
            plot.title = element_text(
                face = "bold", size = 11, hjust = 0,
                margin = margin(b = 3)
            ),
            plot.subtitle = element_text(
                size = 8.5, color = "grey40", hjust = 0,
                margin = margin(b = 8)
            ),
            legend.title = element_text(face = "bold", size = 9),
            legend.text = element_text(size = 8),
            plot.margin = margin(12, 12, 8, 8)
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
        plt,
        width = w, height = h, dpi = 1200, bg = "white"
    )
    tryCatch(
        ggsave(file.path(fig_dir, paste0(name, ".svg")),
            plt,
            width = w, height = h, bg = "white"
        ),
        error = function(e) cat(sprintf("  [SVG 跳过: %s]\n", e$message))
    )
    cat(sprintf("  ✓ %s 已保存\n", name))
}

# ══════════════════════════════════════════════════════════════════════════════
# 读取数据
# ══════════════════════════════════════════════════════════════════════════════
d <- read.csv("output/case4_plant/all_results_plant.csv", stringsAsFactors = FALSE) %>%
    mutate(auc_mean = as.numeric(auc_mean)) %>%
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
    row.names = FALSE
)

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
        axis.text.y = element_text(size = 5.5, face = "italic"),
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
    row.names = FALSE
)

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

cat("\n════════════════════════════════════════\n")
cat("  Fig 4 完成: 2 个子图 + 2 个数据表已保存\n")
cat("════════════════════════════════════════\n")
