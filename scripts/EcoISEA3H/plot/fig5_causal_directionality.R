################################################################################
# Fig 5: CAST 提供方向性因果洞察 — 超越相关性重要性
#
# 核心叙事: 传统 RF importance 只能给出"变量有多重要"(正值标量),
#           CAST 的 ATE 进一步揭示"因果效应是正还是负",
#           且同一变量在不同物种中方向可以完全翻转
#
# 三个独立子图（每图单独保存 1200 dpi PNG + SVG）:
#   (a) ATE 方向异质性展示: 选取翻转率最高的变量，展示其在所有物种中的正/负 ATE
#   (b) RF importance vs |ATE| 散点: 颜色区分因果方向, RF 无法给出的信息
#   (c) 物种间 Spearman ρ 分布: 预测重要性排名 ≠ 因果效应排名
#
# 数据来源:
#   output/case2_eco/all_ate_results_v3.csv
#   output/case2_eco/all_screening_v3.csv
#
# 运行: setwd("E:/CausalSDMs")
#       source("scripts/EcoISEA3H/plot/fig5_causal_directionality.R")
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
# 全局主题
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
# 变量显示名称
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
ate <- read.csv("output/case2_eco/all_ate_results_v3.csv", stringsAsFactors = FALSE)
scr <- read.csv("output/case2_eco/all_screening_v3.csv", stringsAsFactors = FALSE)

n_sp <- n_distinct(ate$species)
cat(sprintf("数据加载完成: %d 物种\n", n_sp))

# ══════════════════════════════════════════════════════════════════════════════
# (a) ATE 方向异质性展示 — 选取翻转率最高的 4 个变量
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── 绘制 (a) ATE 方向异质性 ──\n")

flip_stats <- ate %>%
    group_by(variable) %>%
    summarise(
        n_pos     = sum(coef > 0, na.rm = TRUE),
        n_neg     = sum(coef < 0, na.rm = TRUE),
        n_total   = n_pos + n_neg,
        flip_rate = min(n_pos, n_neg) / pmax(n_total, 1),
        .groups   = "drop"
    ) %>%
    arrange(desc(flip_rate))

write.csv(flip_stats, file.path(tbl_dir, "fig5a_ate_flip_rates.csv"),
          row.names = FALSE)

cat("变量方向翻转率:\n")
print(flip_stats)

# 选前 4 个
top4_vars <- flip_stats$variable[1:min(4, nrow(flip_stats))]

ate_top4 <- ate %>%
    filter(variable %in% top4_vars) %>%
    mutate(
        direction = ifelse(coef > 0, "Positive", "Negative"),
        var_label = get_var_label(variable),
        sp_label  = fmt_sp(species)
    )

# 每个 facet 内按 ATE 值排序物种
ate_top4 <- ate_top4 %>%
    group_by(variable) %>%
    mutate(sp_rank = rank(coef)) %>%
    ungroup() %>%
    arrange(variable, sp_rank)

# 使用 tidytext::reorder_within 的等效方法: 添加 facet 标签到排序 key
ate_top4 <- ate_top4 %>%
    mutate(sp_id = paste0(sp_label, "___", variable))

pa <- ggplot(ate_top4,
             aes(x = reorder(sp_id, coef), y = coef, fill = direction)) +
    geom_col(width = 0.78, show.legend = TRUE) +
    geom_hline(yintercept = 0, linewidth = 0.4, color = "grey30") +
    facet_wrap(~ var_label, scales = "free_x", ncol = 2) +
    scale_fill_manual(
        values = c("Positive" = "#B2182B", "Negative" = "#2166AC"),
        name   = "Causal\ndirection"
    ) +
    scale_x_discrete(labels = function(x) gsub("___.*", "", x)) +
    labs(
        title = "(a) ATE direction heterogeneity across species",
        subtitle = "Same environmental variable exerts opposite causal effects on different species",
        x = "", y = "ATE coefficient"
    ) +
    theme_pub() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 4.5, face = "italic"),
        strip.text  = element_text(face = "bold", size = 9),
        panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3)
    )

save_fig(pa, "fig5a_ate_direction_heterogeneity", w = 12, h = 8)

# ══════════════════════════════════════════════════════════════════════════════
# (b) RF importance vs |ATE| — 颜色编码因果方向
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── 绘制 (b) RF importance vs |ATE| ──\n")

combined <- scr %>%
    left_join(
        ate %>% select(species, variable, ate_coef = coef),
        by = c("species", "variable")
    ) %>%
    mutate(
        direction = ifelse(ate_coef > 0, "Positive", "Negative"),
        abs_ate_coef = abs(ate_coef),
        var_label = get_var_label(variable)
    ) %>%
    filter(!is.na(importance), !is.na(ate_coef))

pb <- ggplot(combined, aes(x = importance, y = abs_ate_coef, color = direction)) +
    geom_point(alpha = 0.55, size = 1.8) +
    scale_color_manual(
        values = c("Positive" = "#B2182B", "Negative" = "#2166AC"),
        name   = "ATE\ndirection"
    ) +
    labs(
        title = "(b) RF importance captures magnitude; ATE reveals direction",
        subtitle = paste0(
            "RF importance (x-axis) is always positive — ",
            "ATE direction (color) provides additional causal insight"
        ),
        x = "RF permutation importance (always positive)",
        y = "|ATE coefficient|"
    ) +
    theme_pub() +
    theme(
        panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
        legend.position  = "right"
    )

save_fig(pb, "fig5b_rf_vs_ate_direction", w = 7, h = 5.5)

# ══════════════════════════════════════════════════════════════════════════════
# (c) Spearman ρ 分布: 预测重要性 ≠ 因果效应
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── 绘制 (c) Spearman ρ 分布 ──\n")

rank_rho <- combined %>%
    group_by(species) %>%
    summarise(
        rho = tryCatch(
            cor(rank(-importance), rank(-abs_ate_coef),
                method = "spearman", use = "complete.obs"),
            error = function(e) NA_real_
        ),
        n_vars = n(),
        .groups = "drop"
    ) %>%
    filter(!is.na(rho))

write.csv(rank_rho, file.path(tbl_dir, "fig5c_spearman_rho.csv"),
          row.names = FALSE)

mean_rho <- mean(rank_rho$rho, na.rm = TRUE)
pct_low  <- 100 * mean(rank_rho$rho < 0.5, na.rm = TRUE)

cat(sprintf("  Mean Spearman ρ = %.3f, %.0f%% species with ρ < 0.5\n",
            mean_rho, pct_low))

pc <- ggplot(rank_rho, aes(x = rho)) +
    geom_histogram(
        binwidth = 0.08, fill = "#4DBBD5", alpha = 0.75,
        color = "white", linewidth = 0.3
    ) +
    geom_vline(
        xintercept = mean_rho, linetype = "dashed",
        color = "#E64B35", linewidth = 0.7
    ) +
    annotate("text",
             x = mean_rho + 0.03, y = Inf, vjust = 2,
             label = sprintf("Mean \u03c1 = %.3f", mean_rho),
             size = 3.5, fontface = "bold", color = "#E64B35", hjust = 0) +
    annotate("text",
             x = -0.1, y = Inf, vjust = 3.5,
             label = sprintf("%.0f%% species with \u03c1 < 0.5", pct_low),
             size = 3, fontface = "italic", color = "grey40", hjust = 0) +
    scale_x_continuous(
        limits = c(-0.5, 1), breaks = seq(-0.4, 1, by = 0.2)
    ) +
    labs(
        title = "(c) Rank concordance: RF importance vs |ATE|",
        subtitle = paste0(
            "Low Spearman \u03c1 = predictive importance ranking \u2260 ",
            "causal effect ranking"
        ),
        x = "Spearman \u03c1  (RF importance rank vs |ATE| rank per species)",
        y = "Number of species"
    ) +
    theme_pub() +
    theme(panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3))

save_fig(pc, "fig5c_rank_concordance", w = 7, h = 4.5)

cat("\n════════════════════════════════════════\n")
cat("  Fig 5 完成: 3 个子图 + 2 个数据表已保存\n")
cat("════════════════════════════════════════\n")
