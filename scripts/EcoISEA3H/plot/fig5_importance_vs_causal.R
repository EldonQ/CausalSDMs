################################################################################
# Fig 5: Causal Inference Decouples Predictive Importance from Ecological
#          Mechanism
#
# Shows that RF permutation importance ranks and |ATE| causal effect ranks
# are frequently discordant, demonstrating that correlation-based importance
# ≠ causal mechanism.
#
# Panel (a): Slope plot — variable-level mean rank (RF vs |ATE|) across species
# Panel (b): Histogram of per-species Spearman ρ (RF rank vs |ATE| rank)
# Panel (c): Dumbbell plot — top variables comparing RF vs |ATE| rank
#
# Data required:
#   output/case2_eco/all_screening_v3.csv
#   output/case2_eco/all_ate_results_v3.csv
#
# Run: setwd("E:/CausalSDMs")
#      source("scripts/EcoISEA3H/plot/fig5_importance_vs_causal.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

fig_dir <- "figures/case2_eco/plot"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(patchwork)

# ── Theme ────────────────────────────────────────────────────────────────────
theme_pub <- function(base_size = 11) {
    theme_minimal(base_size = base_size, base_family = "sans") +
        theme(
            panel.grid.minor  = element_blank(),
            axis.title        = element_text(face = "bold"),
            plot.title        = element_text(face = "bold", hjust = 0),
            plot.subtitle     = element_text(hjust = 0, color = "grey40", size = 9),
            legend.background = element_rect(fill = "white", color = NA)
        )
}

# ── Load data ────────────────────────────────────────────────────────────────
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

# ── Join and compute ranks ──────────────────────────────────────────────────
rank_data <- screen %>%
    left_join(
        ate %>% select(region, species, variable, ate_coef = coef),
        by = c("region", "species", "variable")
    ) %>%
    mutate(ate_coef = as.numeric(ate_coef)) %>%
    group_by(species) %>%
    mutate(
        rank_rf  = rank(-importance, ties.method = "average"),
        rank_ate = rank(-abs(ate_coef), ties.method = "average"),
        n_vars   = n()
    ) %>%
    ungroup()

# ══════════════════════════════════════════════════════════════════════════════
# Panel (a): Slope plot — variable-level mean rank
# ══════════════════════════════════════════════════════════════════════════════
slope_dat <- rank_data %>%
    group_by(variable) %>%
    summarise(
        mean_rank_rf  = mean(rank_rf, na.rm = TRUE),
        mean_rank_ate = mean(rank_ate, na.rm = TRUE),
        n_sp          = n(),
        .groups       = "drop"
    ) %>%
    filter(n_sp >= 5) %>%
    mutate(
        y_left = rank(mean_rank_rf),
        y_right = rank(mean_rank_ate),
        dir_change = sign(y_right - y_left)
    )

pa <- ggplot(slope_dat) +
    geom_segment(
        aes(
            x = 1, xend = 2, y = y_left, yend = y_right,
            color = factor(dir_change)
        ),
        linewidth = 1.1, alpha = 0.8
    ) +
    geom_point(aes(x = 1, y = y_left), size = 2, color = "#27AE60") +
    geom_point(aes(x = 2, y = y_right), size = 2, color = "#C0392B") +
    geom_text(aes(x = 0.85, y = y_left, label = variable),
        size = 2.5, hjust = 1
    ) +
    geom_text(aes(x = 2.15, y = y_right, label = variable),
        size = 2.5, hjust = 0
    ) +
    scale_color_manual(
        values = c("-1" = "#E74C3C", "0" = "grey60", "1" = "#2980B9"),
        guide = "none"
    ) +
    scale_x_continuous(
        breaks = c(1, 2),
        labels = c("RF importance\nrank", "|ATE| causal\nrank"),
        limits = c(0.3, 2.7)
    ) +
    labs(
        title = "(a) Variable-level rank discordance",
        subtitle = "Mean rank across all species | Line crossings = systematic discordance",
        x = "", y = "Mean rank across species (1 = highest)"
    ) +
    theme_pub() +
    theme(legend.position = "right")

# ══════════════════════════════════════════════════════════════════════════════
# Panel (b): Per-species Spearman ρ histogram
# ══════════════════════════════════════════════════════════════════════════════
spearman_per_sp <- rank_data %>%
    group_by(species) %>%
    summarise(
        rho = tryCatch(
            cor(rank_rf, rank_ate, method = "spearman", use = "complete.obs"),
            error = function(e) NA_real_
        ),
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
        color = "#C0392B", linewidth = 0.7
    ) +
    annotate("text",
        x = mean_rho, y = Inf,
        label = sprintf("Mean ρ = %.3f", mean_rho),
        vjust = 1.5, hjust = -0.1, size = 3.5,
        fontface = "bold", color = "#C0392B"
    ) +
    annotate("text",
        x = 0.25, y = Inf,
        label = sprintf("%.0f%% species\nhave ρ < 0.5", pct_low),
        vjust = 2.5, size = 3.2, fontface = "italic", color = "grey40"
    ) +
    labs(
        title = "(b) Per-species rank concordance",
        subtitle = sprintf(
            "n = %d species | Low ρ = causal rank ≠ predictive rank",
            n_sp_rho
        ),
        x = "Spearman ρ  (RF importance rank vs |ATE| rank)",
        y = "Number of species"
    ) +
    theme_pub()

# ══════════════════════════════════════════════════════════════════════════════
# Panel (c): Dumbbell for top variables
# ══════════════════════════════════════════════════════════════════════════════
top_vars <- slope_dat %>%
    arrange(mean_rank_rf) %>%
    slice_head(n = 10) %>%
    pull(variable)

top_dat <- slope_dat %>%
    filter(variable %in% top_vars) %>%
    pivot_longer(
        cols = c(y_left, y_right),
        names_to = "metric", values_to = "rank_val"
    ) %>%
    mutate(
        metric = factor(
            metric,
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
        values = c(
            "RF importance rank" = "#27AE60",
            "|ATE| causal rank" = "#C0392B"
        ),
        name = ""
    ) +
    labs(
        title = "(c) Top-10 variables: RF vs causal rank",
        subtitle = "Large horizontal gaps = importance-mechanism discordance",
        x = "Rank position", y = ""
    ) +
    theme_pub() +
    theme(legend.position = "bottom")

# ══════════════════════════════════════════════════════════════════════════════
# Combine
# ══════════════════════════════════════════════════════════════════════════════
fig5 <- (pa | pb) / pc +
    plot_layout(heights = c(1, 1.1)) +
    plot_annotation(
        title = "Fig 5  Causal inference decouples predictive importance from ecological mechanism",
        subtitle = paste0(
            "RF permutation importance ranks and |ATE| causal effect ranks show ",
            "low concordance, demonstrating that predictive relevance ≠ causal mechanism"
        ),
        theme = theme(
            plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
            plot.subtitle = element_text(
                face = "italic", size = 9.5, hjust = 0.5,
                color = "grey40"
            )
        )
    )

ggsave(file.path(fig_dir, "fig5_importance_vs_causal.png"),
    fig5,
    width = 14, height = 12, dpi = 300, bg = "white"
)
cat("✓ Saved fig5_importance_vs_causal.png\n")
