################################################################################
# Fig 5: Cross-Region Robustness (3-panel)
# (a) Per-model AUC variability (IQR/SD across regions)
# (b) Region × Model AUC heatmap
# (c) Mean rank ± SD stability
#
# Requires: all_results_v3.csv
# Run: setwd("E:/CausalSDMs"); source("scripts/case2/plot/fig5_cross_region_robustness.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

fig_dir <- "figures/case2/plot"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(patchwork)

theme_pub <- function(base_size = 11) {
    theme_minimal(base_size = base_size, base_family = "sans") +
        theme(
            panel.grid.minor = element_blank(),
            axis.title = element_text(face = "bold"),
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, color = "grey40")
        )
}

# ---- Load data ----
d <- read.csv("output/case2/all_results_v3.csv", stringsAsFactors = FALSE) %>%
    filter(!is.na(auc_mean))

model_order <- c("CAST", "MLP_ATE", "MLP", "RF", "BRT", "Maxent")
model_colors <- c(
    "CAST" = "#2980B9", "MLP_ATE" = "#5DADE2", "MLP" = "#85C1E9",
    "RF" = "#27AE60", "BRT" = "#E67E22", "Maxent" = "#9B59B6"
)

# ══════════════════════════════════════════════════════════════
# (a) Per-model AUC variability across regions (IQR)
# ══════════════════════════════════════════════════════════════
# For each model, compute across-region AUC summary
region_model_auc <- d %>%
    group_by(region, model) %>%
    summarise(region_mean_auc = mean(auc_mean, na.rm = TRUE), .groups = "drop")

model_variability <- region_model_auc %>%
    group_by(model) %>%
    summarise(
        mean_auc = mean(region_mean_auc),
        sd_auc = sd(region_mean_auc),
        iqr_auc = IQR(region_mean_auc),
        min_auc = min(region_mean_auc),
        max_auc = max(region_mean_auc),
        range_auc = max_auc - min_auc,
        .groups = "drop"
    ) %>%
    mutate(model = factor(model, levels = model_order))

pa <- ggplot(model_variability, aes(x = reorder(model, sd_auc), y = sd_auc, fill = model)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_text(aes(label = sprintf("%.4f", sd_auc)), hjust = -0.15, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = model_colors, guide = "none") +
    coord_flip() +
    labs(
        title = "(a) Cross-region AUC variability",
        subtitle = "SD of mean AUC across 6 regions | Lower = more stable",
        x = "", y = "SD of region-mean AUC"
    ) +
    theme_pub()

# ══════════════════════════════════════════════════════════════
# (b) Region × Model AUC heatmap
# ══════════════════════════════════════════════════════════════
heatmap_data <- region_model_auc %>%
    mutate(model = factor(model, levels = model_order))

pb <- ggplot(heatmap_data, aes(x = model, y = region, fill = region_mean_auc)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", region_mean_auc)), size = 3.2, fontface = "bold") +
    scale_fill_gradient2(
        low = "#E74C3C", mid = "#F9E79F", high = "#27AE60",
        midpoint = median(heatmap_data$region_mean_auc),
        name = "Mean AUC"
    ) +
    labs(
        title = "(b) Region × Model performance heatmap",
        subtitle = "Mean AUC per (region, model) cell",
        x = "", y = ""
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")

# ══════════════════════════════════════════════════════════════
# (c) Mean rank ± SD (stability of ranking)
# ══════════════════════════════════════════════════════════════
rank_per_sp <- d %>%
    group_by(region, species) %>%
    mutate(rank = dense_rank(desc(auc_mean))) %>%
    ungroup()

rank_stability <- rank_per_sp %>%
    group_by(model) %>%
    summarise(
        mean_rank = mean(rank),
        sd_rank = sd(rank),
        median_rank = median(rank),
        .groups = "drop"
    ) %>%
    mutate(model = factor(model, levels = model_order))

pc <- ggplot(rank_stability, aes(x = reorder(model, mean_rank), y = mean_rank, fill = model)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_errorbar(aes(ymin = mean_rank - sd_rank, ymax = mean_rank + sd_rank),
        width = 0.2, linewidth = 0.6
    ) +
    geom_text(aes(label = sprintf("%.1f±%.1f", mean_rank, sd_rank)),
        vjust = -0.8, size = 3.2, fontface = "bold"
    ) +
    scale_fill_manual(values = model_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
        title = "(c) Mean rank ± SD across all species",
        subtitle = "Lower rank = better | Smaller SD = more consistent",
        x = "", y = "Mean rank (1 = best)"
    ) +
    theme_pub()

# ══════════════════════════════════════════════════════════════
# Combine
# ══════════════════════════════════════════════════════════════
fig5 <- (pa | pb) / pc +
    plot_annotation(
        title = "Fig 5  Cross-Region Robustness",
        theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    )

ggsave(file.path(fig_dir, "fig5_cross_region_robustness.png"),
    fig5,
    width = 13, height = 10, dpi = 300, bg = "white"
)
ggsave(file.path(fig_dir, "fig5_cross_region_robustness.svg"),
    fig5,
    width = 13, height = 10, bg = "white"
)
cat("Saved fig5_cross_region_robustness.png/.svg\n")
