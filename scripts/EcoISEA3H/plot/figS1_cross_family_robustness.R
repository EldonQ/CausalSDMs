################################################################################
# Fig S1 (Eco): Cross-Family Robustness — Supplementary
#   [Consolidates former fig5_cross_family_robustness_eco.R]
#
# Narrative: CAST predictive gains are consistent across taxonomic families,
#   demonstrating that the framework is not tuned to any specific clade.
#   This supplementary figure provides the detailed family-stratified evidence
#   that supports the cross-taxon generalisation claim in the main text.
#
# Panel (a): Boxplot per model, faceted by Taxonomic Family
#            (CAST vs all baselines within each family)
# Panel (b): Family × Model AUC heatmap
#            (mean AUC per cell; colour scale centred at dataset median)
# Panel (c): Model rank stability across families
#            (mean rank ± SD per model; lower rank = better)
#
# Data required:
#   output/case2_eco/all_results_v3.csv
#   outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv
#
# Run: setwd("E:/CausalSDMs")
#      source("scripts/EcoISEA3H/plot/figS1_cross_family_robustness.R")
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
            panel.grid.minor  = element_blank(),
            axis.title        = element_text(face = "bold"),
            plot.title        = element_text(face = "bold", hjust = 0.5),
            plot.subtitle     = element_text(hjust = 0.5, color = "grey40"),
            legend.background = element_rect(fill = "white", color = NA)
        )
}

# ── Palette ──────────────────────────────────────────────────────────────────
model_order <- c("CAST", "MLP_ATE", "MLP", "RF", "BRT", "Maxent")
model_colors <- c(
    "CAST"    = "#2980B9",
    "MLP_ATE" = "#5DADE2",
    "MLP"     = "#85C1E9",
    "RF"      = "#27AE60",
    "BRT"     = "#E67E22",
    "Maxent"  = "#9B59B6"
)

# ── Load data ─────────────────────────────────────────────────────────────────
d <- read.csv("output/case2_eco/all_results_v3.csv",
    stringsAsFactors = FALSE
) %>%
    filter(!is.na(auc_mean))

sp_meta <- read.csv(
    "outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv",
    stringsAsFactors = FALSE
) %>%
    mutate(species = gsub(" ", "_", species))

d <- d %>%
    left_join(sp_meta %>% select(species, family, category), by = "species") %>%
    mutate(family = ifelse(is.na(family), "Unknown", family))

# Exclude families with <2 species (would give degenerate boxplots)
valid_fams <- d %>%
    group_by(family) %>%
    summarise(n = n_distinct(species), .groups = "drop") %>%
    filter(n >= 2) %>%
    pull(family)

d <- d %>%
    filter(family %in% valid_fams) %>%
    mutate(model = factor(model, levels = model_order))

# ══════════════════════════════════════════════════════════════════════════════
# Panel (a): Faceted boxplot — model × family
# ══════════════════════════════════════════════════════════════════════════════
pa <- ggplot(d, aes(x = model, y = auc_mean, fill = model)) +
    geom_boxplot(
        alpha = 0.70, outlier.size = 0.8, outlier.alpha = 0.5,
        width = 0.6
    ) +
    stat_summary(
        fun = mean, geom = "point",
        shape = 18, size = 2.5, color = "black"
    ) +
    facet_wrap(~family, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = model_colors, guide = "none") +
    labs(
        title = "(a) CAST vs baselines across Taxonomic Families",
        subtitle = "Boxplot per model–family combination | Diamond = mean AUC",
        x = "", y = "AUC"
    ) +
    theme_pub() +
    theme(
        axis.text.x = element_text(angle = 40, hjust = 1, size = 8),
        strip.text = element_text(face = "bold"),
        panel.spacing = unit(0.8, "lines")
    )

# ══════════════════════════════════════════════════════════════════════════════
# Panel (b): Family × Model AUC heatmap
# ══════════════════════════════════════════════════════════════════════════════
fam_model_auc <- d %>%
    group_by(family, model) %>%
    summarise(fam_mean_auc = mean(auc_mean, na.rm = TRUE), .groups = "drop") %>%
    mutate(model = factor(model, levels = model_order))

median_auc <- median(fam_model_auc$fam_mean_auc, na.rm = TRUE)

pb <- ggplot(fam_model_auc, aes(x = model, y = family, fill = fam_mean_auc)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", fam_mean_auc)),
        size = 3.0, fontface = "bold", color = "grey10"
    ) +
    scale_fill_gradient2(
        low      = "#E74C3C",
        mid      = "#FDEBD0",
        high     = "#27AE60",
        midpoint = median_auc,
        name     = "Mean\nAUC"
    ) +
    labs(
        title = "(b) Family × Model performance heatmap",
        subtitle = sprintf("Colour centred at dataset median AUC = %.3f", median_auc),
        x = "", y = ""
    ) +
    theme_pub() +
    theme(
        axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "right"
    )

# ══════════════════════════════════════════════════════════════════════════════
# Panel (c): Model rank stability — mean rank ± SD across all species
# ══════════════════════════════════════════════════════════════════════════════
rank_per_sp <- d %>%
    group_by(species) %>%
    mutate(rank = dense_rank(desc(auc_mean))) %>%
    ungroup()

rank_stability <- rank_per_sp %>%
    group_by(model) %>%
    summarise(
        mean_rank   = mean(rank),
        sd_rank     = sd(rank),
        median_rank = median(rank),
        .groups     = "drop"
    ) %>%
    mutate(model = factor(model, levels = model_order))

pc <- ggplot(
    rank_stability,
    aes(x = reorder(model, mean_rank), y = mean_rank, fill = model)
) +
    geom_col(width = 0.6, alpha = 0.85) +
    geom_errorbar(
        aes(ymin = mean_rank - sd_rank, ymax = mean_rank + sd_rank),
        width = 0.22, linewidth = 0.7
    ) +
    geom_text(
        aes(label = sprintf("%.1f ± %.1f", mean_rank, sd_rank)),
        vjust = -0.9, size = 3.2, fontface = "bold"
    ) +
    scale_fill_manual(values = model_colors, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
        title = "(c) Model rank stability across all species",
        subtitle = "Mean rank ± SD | Lower rank = better | Smaller SD = more consistent",
        x = "", y = "Mean rank (1 = best)"
    ) +
    theme_pub()

# ══════════════════════════════════════════════════════════════════════════════
# Combine
# ══════════════════════════════════════════════════════════════════════════════
figS1 <- pa / (pb | pc) +
    plot_layout(heights = c(2, 1)) +
    plot_annotation(
        title = "Fig S1 (Eco)  Cross-Family Robustness of CAST",
        subtitle = paste0(
            "Performance comparisons stratified by Taxonomic Family demonstrate ",
            "that CAST's predictive gains are not confined to any particular clade ",
            "or ecological guild."
        ),
        theme = theme(
            plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
            plot.subtitle = element_text(
                face = "italic", size = 10, hjust = 0.5,
                color = "grey40"
            )
        )
    )

ggsave(file.path(fig_dir, "figS1_cross_family_robustness_eco.png"),
    figS1,
    width = 14, height = 14, dpi = 300, bg = "white"
)
cat("Saved figS1_cross_family_robustness_eco.png\n")
