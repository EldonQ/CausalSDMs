################################################################################
# Fig 5 (Eco): Cross-Family Robustness (3-panel)
# (a) Per-model AUC variability (IQR/SD across taxonomic families)
# (b) Family × Model AUC heatmap
# (c) Mean rank ± SD stability
#
# Requires: all_results_v3.csv, CAST_Species_Summary.csv
# Run: setwd("E:/CausalSDMs"); source("scripts/EcoISEA3H/plot/fig5_cross_family_robustness_eco.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

fig_dir <- "figures/case2_eco/plot"
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
d <- read.csv("output/case2_eco/all_results_v3.csv", stringsAsFactors = FALSE) %>%
    filter(!is.na(auc_mean))

sp_meta <- read.csv("outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv", stringsAsFactors = FALSE) %>%
    mutate(species = gsub(" ", "_", species))

d <- d %>%
    left_join(sp_meta %>% select(species, family, category), by = "species") %>%
    mutate(family = ifelse(is.na(family), "Unknown", family))

model_order <- c("CAST", "MLP_ATE", "MLP", "RF", "BRT", "Maxent")
model_colors <- c(
    "CAST" = "#2980B9", "MLP_ATE" = "#5DADE2", "MLP" = "#85C1E9",
    "RF" = "#27AE60", "BRT" = "#E67E22", "Maxent" = "#9B59B6"
)

# ══════════════════════════════════════════════════════════════
# (a) Per-model AUC variability across families (SD)
# ══════════════════════════════════════════════════════════════
fam_model_auc <- d %>%
    group_by(family, model) %>%
    summarise(fam_mean_auc = mean(auc_mean, na.rm = TRUE), .groups = "drop")

model_variability <- fam_model_auc %>%
    group_by(model) %>%
    summarise(
        mean_auc = mean(fam_mean_auc),
        sd_auc = sd(fam_mean_auc),
        iqr_auc = IQR(fam_mean_auc),
        min_auc = min(fam_mean_auc),
        max_auc = max(fam_mean_auc),
        range_auc = max_auc - min_auc,
        .groups = "drop"
    ) %>%
    mutate(model = factor(model, levels = model_order))

pa <- ggplot(model_variability, aes(x = reorder(model, sd_auc), y = sd_auc, fill = model)) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_text(aes(label = sprintf("%.4f", sd_auc)), hjust = -0.15, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = model_colors, guide = "none") +
    scale_y_continuous(limits = c(0, max(model_variability$sd_auc, na.rm = TRUE) * 1.2)) +
    coord_flip() +
    labs(
        title = "(a) Cross-family AUC variability",
        subtitle = "SD of mean AUC across Taxonomic Families | Lower = more stable",
        x = "", y = "SD of family-mean AUC"
    ) +
    theme_pub()

# ══════════════════════════════════════════════════════════════
# (b) Family × Model AUC heatmap
# ══════════════════════════════════════════════════════════════
heatmap_data <- fam_model_auc %>%
    mutate(model = factor(model, levels = model_order))

pb <- ggplot(heatmap_data, aes(x = model, y = family, fill = fam_mean_auc)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", fam_mean_auc)), size = 3.2, fontface = "bold", color = "black") +
    scale_fill_gradient2(
        low = "#E74C3C", mid = "#F9E79F", high = "#27AE60",
        midpoint = median(heatmap_data$fam_mean_auc, na.rm = TRUE),
        name = "Mean AUC"
    ) +
    labs(
        title = "(b) Taxonomic Family × Model performance heatmap",
        subtitle = "Mean AUC per (family, model) cell",
        x = "", y = ""
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1), legend.position = "right")

# ══════════════════════════════════════════════════════════════
# (c) Mean rank ± SD (stability of ranking)
# ══════════════════════════════════════════════════════════════
rank_per_sp <- d %>%
    group_by(species) %>%
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
        title = "Fig 5 (Eco) Cross-Family Robustness",
        theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    )

ggsave(file.path(fig_dir, "fig5_cross_family_robustness_eco.png"),
    fig5,
    width = 13, height = 10, dpi = 300, bg = "white"
)
cat("Saved fig5_cross_family_robustness_eco.png\n")
