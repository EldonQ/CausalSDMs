################################################################################
# Fig 5 (Eco): Causal Inference Decouples Predictive Importance from
#              Ecological Mechanism  [INTERPRETABILITY PANEL]
#
# Narrative: Traditional variable importance (RF permutation) mixes causal
#   effects with spurious collinearity. DML-based ATE isolates the direct
#   causal contribution of each environmental driver. The low / discordant
#   Spearman ρ between RF-rank and |ATE|-rank demonstrates that machine
#   learning importance is a poor proxy for ecological mechanism.
#
# Panel (a): Slope chart — mean RF-rank vs mean |ATE|-rank per variable
#            (aggregated across all species; crossing lines = discordance)
# Panel (b): Histogram of per-species Spearman ρ (RF rank vs |ATE| rank)
#            — distribution centred well below 1 confirms systematic divergence
# Panel (c): Lollipop — top-10 variables ranked by |ATE| vs ranked by RF
#            importance (side-by-side bars, difference highlighted)
#
# Data required:
#   output/case2_eco/all_screening_v3.csv
#   output/case2_eco/all_ate_results_v3.csv
#
# Run: setwd("E:/CausalSDMs")
#      source("scripts/EcoISEA3H/plot/fig5_importance_vs_causal_eco.R")
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

# ── Load data ─────────────────────────────────────────────────────────────────
screen <- read.csv("output/case2_eco/all_screening_v3.csv",
    stringsAsFactors = FALSE
)
ate <- read.csv("output/case2_eco/all_ate_results_v3.csv",
    stringsAsFactors = FALSE
) %>%
    mutate(
        coef        = as.numeric(coef),
        significant = as.logical(significant)
    )

# Join ATE coefficients into screening table
rank_data <- screen %>%
    left_join(
        ate %>% select(region, species, variable, ate_coef = coef),
        by = c("region", "species", "variable")
    ) %>%
    mutate(ate_coef = as.numeric(ate_coef)) %>%
    group_by(species) %>%
    mutate(
        rank_rf  = rank(-importance, ties.method = "average"),
        rank_ate = rank(-abs(ate_coef), ties.method = "average", na.last = TRUE)
    ) %>%
    ungroup()

# ══════════════════════════════════════════════════════════════════════════════
# Panel (a): Slope chart — mean RF-rank vs mean |ATE|-rank across ALL species
# ══════════════════════════════════════════════════════════════════════════════
slope_dat <- rank_data %>%
    group_by(variable) %>%
    summarise(
        mean_rank_rf  = mean(rank_rf, na.rm = TRUE),
        mean_rank_ate = mean(rank_ate, na.rm = TRUE),
        n_sp          = n(),
        .groups       = "drop"
    ) %>%
    filter(n_sp >= 5) %>% # require enough species for stability
    mutate(
        y_left = rank(mean_rank_rf),
        y_right = rank(mean_rank_ate),
        dir_change = sign(y_right - y_left)
    )

# Colour by direction of rank change
pa <- ggplot(slope_dat) +
    geom_segment(
        aes(
            x = 1, xend = 2, y = y_left, yend = y_right,
            color = factor(dir_change)
        ),
        linewidth = 1.1, alpha = 0.8
    ) +
    geom_point(aes(x = 1, y = y_left), size = 2.8, color = "#27AE60") +
    geom_point(aes(x = 2, y = y_right), size = 2.8, color = "#C0392B") +
    geom_text(aes(x = 0.90, y = y_left, label = variable),
        hjust = 1, size = 2.8, fontface = "bold"
    ) +
    geom_text(aes(x = 2.10, y = y_right, label = variable),
        hjust = 0, size = 2.8, fontface = "bold"
    ) +
    scale_x_continuous(
        limits = c(0.35, 2.65),
        breaks = c(1, 2),
        labels = c("RF importance\nrank", "|ATE| rank")
    ) +
    scale_color_manual(
        values = c("-1" = "#E74C3C", "0" = "grey60", "1" = "#2980B9"),
        labels = c("-1" = "Rank drops", "0" = "Unchanged", "1" = "Rank rises"),
        name   = "Rank change\n(RF → ATE)"
    ) +
    scale_y_reverse() +
    labs(
        title = "(a) RF importance rank vs |ATE| causal rank",
        subtitle = "Mean rank across all species | Line crossings indicate systematic discordance",
        x = "", y = "Mean rank across species (1 = highest)"
    ) +
    theme_pub() +
    theme(legend.position = "right")

# ══════════════════════════════════════════════════════════════════════════════
# Panel (b): Per-species Spearman ρ distribution
# ══════════════════════════════════════════════════════════════════════════════
spearman_per_sp <- rank_data %>%
    filter(!is.na(rank_rf), !is.na(rank_ate)) %>%
    group_by(species) %>%
    summarise(
        rho = cor(rank_rf, rank_ate, method = "spearman", use = "pairwise"),
        n_vars = n(),
        .groups = "drop"
    ) %>%
    filter(!is.na(rho))

mean_rho <- mean(spearman_per_sp$rho, na.rm = TRUE)
pct_low <- 100 * mean(spearman_per_sp$rho < 0.5, na.rm = TRUE)
n_sp_rho <- nrow(spearman_per_sp)

pb <- ggplot(spearman_per_sp, aes(x = rho)) +
    geom_histogram(
        binwidth = 0.08, fill = "#2980B9", alpha = 0.70,
        color = "white"
    ) +
    geom_vline(
        xintercept = mean_rho, linetype = "dashed",
        color = "#C0392B", linewidth = 0.9
    ) +
    geom_vline(
        xintercept = 1, linetype = "dotted",
        color = "grey50", linewidth = 0.6
    ) +
    geom_vline(
        xintercept = 0, linetype = "dotted",
        color = "grey50", linewidth = 0.6
    ) +
    annotate("text",
        x = mean_rho + 0.06, y = Inf, vjust = 2,
        label = sprintf("Mean ρ = %.2f", mean_rho),
        fontface = "bold", size = 3.4, color = "#C0392B"
    ) +
    annotate("text",
        x = -0.95, y = Inf, vjust = 2,
        label = sprintf("%.0f%% of\nspecies\nρ < 0.5", pct_low),
        fontface = "italic", size = 3.0, color = "#2980B9", hjust = 0
    ) +
    scale_x_continuous(
        limits = c(-1.05, 1.05),
        breaks = seq(-1, 1, 0.5)
    ) +
    labs(
        title = "(b) Per-species Spearman ρ: RF rank vs |ATE| rank",
        subtitle = sprintf(
            "n = %d species | Low ρ indicates causal rank ≠ predictive rank",
            n_sp_rho
        ),
        x = "Spearman ρ  (RF importance rank vs |ATE| rank)",
        y = "Number of species"
    ) +
    theme_pub()

# ══════════════════════════════════════════════════════════════════════════════
# Panel (c): Top-N variables — RF rank vs ATE rank side-by-side
# ══════════════════════════════════════════════════════════════════════════════
# Variables with highest absolute discrepancy between mean RF rank and mean |ATE| rank
top_vars <- slope_dat %>%
    mutate(rank_gap = abs(y_right - y_left)) %>%
    arrange(desc(rank_gap)) %>%
    slice_head(n = 12) %>%
    pull(variable)

top_dat <- slope_dat %>%
    filter(variable %in% top_vars) %>%
    pivot_longer(
        cols = c(y_left, y_right),
        names_to = "metric",
        values_to = "rank_val"
    ) %>%
    mutate(
        metric = factor(metric,
            levels = c("y_left", "y_right"),
            labels = c("RF importance rank", "|ATE| causal rank")
        )
    )

pc <- ggplot(top_dat, aes(
    x = rank_val,
    y = reorder(variable, -mean_rank_rf),
    color = metric, group = variable
)) +
    geom_line(color = "grey70", linewidth = 0.8) +
    geom_point(size = 3.5, alpha = 0.9) +
    scale_color_manual(
        values = c("RF importance rank" = "#27AE60", "|ATE| causal rank" = "#C0392B"),
        name   = ""
    ) +
    scale_x_continuous(
        breaks = scales::pretty_breaks(n = 6),
        trans  = "reverse" # lower rank = better = right side
    ) +
    labs(
        title = "(c) Rank discordance: top-discrepant variables",
        subtitle = "Variables with largest gap between RF rank and |ATE| rank",
        x = "Mean rank (lower rank = higher position)",
        y = ""
    ) +
    theme_pub() +
    theme(legend.position = "bottom")

# ══════════════════════════════════════════════════════════════════════════════
# Combine
# ══════════════════════════════════════════════════════════════════════════════
fig5 <- (pa | pb) / pc +
    plot_layout(heights = c(1, 1.1)) +
    plot_annotation(
        title = "Fig 5 (Eco)  Causal Inference Decouples Predictive Importance from Ecological Mechanism",
        subtitle = paste0(
            "RF importance ranks reflect collinearity structure; DML-ATE ranks measure ",
            "direct causal effects. The systematic discordance (low Spearman ρ) ",
            "reveals that correlational importance is a poor proxy for ecological mechanism."
        ),
        theme = theme(
            plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
            plot.subtitle = element_text(
                face = "italic", size = 9, hjust = 0.5,
                color = "grey40"
            )
        )
    )

ggsave(file.path(fig_dir, "fig5_importance_vs_causal_eco.png"),
    fig5,
    width = 14, height = 12, dpi = 300, bg = "white"
)
cat("Saved fig5_importance_vs_causal_eco.png\n")
