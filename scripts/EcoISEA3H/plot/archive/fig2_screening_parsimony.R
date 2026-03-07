################################################################################
# Fig 2: Causal Screening Achieves Parsimony Without Accuracy Trade-offs
#
# Scientific Question Q1:
#   Can causal screening substantially reduce variable dimensionality
#   while maintaining or improving predictive accuracy?
#
# Panel (a): All-model dot-and-whisker — Mean AUC ± SE across models
#            Shows CAST-screened variables perform on par with full-variable sets
# Panel (b): Per-species ΔAUC distribution (histogram/density)
#            Shows distribution centred near zero ⟹ no systematic loss
# Panel (c): Variable reduction (%) vs ΔAUC scatter
#            Shows parsimony is not antagonistic to accuracy
#
# Data required:
#   output/case2_eco/all_results_v3.csv
#   output/case2_eco/all_screening_v3.csv
#   outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv
#
# Run: setwd("E:/CausalSDMs")
#      source("scripts/EcoISEA3H/plot/fig2_screening_parsimony.R")
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

# ── Palette ──────────────────────────────────────────────────────────────────
model_colors <- c(
    "CAST"    = "#2980B9",
    "MLP_ATE" = "#5DADE2",
    "MLP"     = "#AED6F1",
    "RF"      = "#27AE60",
    "BRT"     = "#E67E22",
    "Maxent"  = "#9B59B6"
)
model_order <- c("CAST", "MLP_ATE", "MLP", "RF", "BRT", "Maxent")

# ── Load data ────────────────────────────────────────────────────────────────
d <- read.csv("output/case2_eco/all_results_v3.csv", stringsAsFactors = FALSE) %>%
    filter(!is.na(auc_mean))

screen <- read.csv("output/case2_eco/all_screening_v3.csv",
    stringsAsFactors = FALSE
)

sp_meta <- read.csv(
    "outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv",
    stringsAsFactors = FALSE
) %>%
    mutate(species = gsub(" ", "_", species))

d <- d %>% left_join(sp_meta %>% select(species, family), by = "species")

n_sp_total <- n_distinct(d$species)

# ══════════════════════════════════════════════════════════════════════════════
# Panel (a): All-model Mean AUC ± SE — dot-and-whisker
# ══════════════════════════════════════════════════════════════════════════════
full_stats <- d %>%
    mutate(model = factor(model, levels = model_order)) %>%
    group_by(model) %>%
    summarise(
        mean_auc = mean(auc_mean, na.rm = TRUE),
        se_auc   = sd(auc_mean, na.rm = TRUE) / sqrt(n()),
        n        = n(),
        .groups  = "drop"
    )

pa <- ggplot(full_stats, aes(x = model, y = mean_auc, color = model)) +
    # Shade CAST family
    annotate("rect",
        xmin = 0.5, xmax = 3.5,
        ymin = -Inf, ymax = Inf, alpha = 0.05, fill = "#2980B9"
    ) +
    annotate("text",
        x = 2, y = min(full_stats$mean_auc - full_stats$se_auc) - 0.005,
        label = "CAST family", fontface = "italic",
        size = 2.8, color = "#2980B9"
    ) +
    geom_pointrange(
        aes(ymin = mean_auc - se_auc, ymax = mean_auc + se_auc),
        size = 0.9, linewidth = 0.9, fatten = 3.5
    ) +
    geom_text(aes(label = sprintf("%.4f", mean_auc)),
        vjust = -1.3, size = 3, fontface = "bold"
    ) +
    scale_color_manual(values = model_colors, guide = "none") +
    scale_x_discrete(limits = model_order) +
    labs(
        title = "(a) All models: Mean AUC ± SE",
        subtitle = sprintf(
            "n = %d species | CAST-screened variables maintain competitive accuracy",
            n_sp_total
        ),
        x = "", y = "Mean AUC"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 22, hjust = 1))

# ══════════════════════════════════════════════════════════════════════════════
# Panel (b): Distribution of per-species ΔAUC (CAST − MLP baseline)
# ══════════════════════════════════════════════════════════════════════════════
paired_wide <- d %>%
    filter(model %in% c("MLP", "CAST")) %>%
    select(species, family, model, auc_mean) %>%
    pivot_wider(names_from = model, values_from = auc_mean, values_fn = max) %>%
    filter(!is.na(CAST), !is.na(MLP)) %>%
    mutate(delta = CAST - MLP)

mean_delta <- mean(paired_wide$delta, na.rm = TRUE)
median_delta <- median(paired_wide$delta, na.rm = TRUE)
pct_positive <- 100 * mean(paired_wide$delta >= 0)

pb <- ggplot(paired_wide, aes(x = delta)) +
    geom_histogram(
        binwidth = 0.005, fill = "#2980B9", alpha = 0.65,
        color = "white", linewidth = 0.3
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.8) +
    geom_vline(
        xintercept = mean_delta, linetype = "solid",
        color = "#C0392B", linewidth = 0.7
    ) +
    annotate("text",
        x = mean_delta, y = Inf,
        label = sprintf("Mean = %+.4f", mean_delta),
        vjust = 1.8, hjust = -0.1, size = 3.2, fontface = "bold", color = "#C0392B"
    ) +
    labs(
        title = "(b) Per-species performance change",
        subtitle = sprintf(
            "ΔAUC (CAST − FlatNN) | %.0f%% of species ≥ 0 | Median = %+.4f",
            pct_positive, median_delta
        ),
        x = expression(Delta * "AUC (CAST − FlatNN)"),
        y = "Number of species"
    ) +
    theme_pub()

# ══════════════════════════════════════════════════════════════════════════════
# Panel (c): Variable reduction vs ΔAUC — parsimony scatter
# ══════════════════════════════════════════════════════════════════════════════
n_full_per_sp <- screen %>%
    group_by(species) %>%
    summarise(n_vars_full = n(), .groups = "drop")

n_cast_per_sp <- d %>%
    filter(model == "CAST") %>%
    group_by(species) %>%
    summarise(n_vars_cast = first(n_vars), .groups = "drop")

var_red <- n_full_per_sp %>%
    left_join(n_cast_per_sp, by = "species") %>%
    mutate(pct_removed = 100 * (1 - n_vars_cast / n_vars_full)) %>%
    filter(!is.na(pct_removed))

delta_df <- paired_wide %>% select(species, delta, family)

pd_data <- var_red %>%
    left_join(delta_df, by = "species") %>%
    filter(!is.na(delta), pct_removed >= 0)

n_total_d <- nrow(pd_data)
mean_pct <- mean(pd_data$pct_removed, na.rm = TRUE)

pc <- ggplot(pd_data, aes(x = pct_removed, y = delta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_point(aes(color = family), size = 2.2, alpha = 0.72) +
    geom_smooth(
        method = "loess", se = TRUE, color = "black",
        linetype = "solid", linewidth = 0.7, fill = "grey85"
    ) +
    scale_color_brewer(palette = "Set1", name = "Family") +
    labs(
        title = "(c) Parsimony vs predictive change",
        subtitle = sprintf(
            "n = %d species | Mean removal = %.1f%% | Above zero = accuracy gain from parsimony",
            n_total_d, mean_pct
        ),
        x = "Variables removed by CAST (%)",
        y = expression(Delta * "AUC (CAST − FlatNN)")
    ) +
    theme_pub() +
    theme(legend.position = "right")

# ══════════════════════════════════════════════════════════════════════════════
# Combine
# ══════════════════════════════════════════════════════════════════════════════
fig2 <- pa / (pb | pc) +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(
        title = "Fig 2  Causal screening achieves high parsimony without accuracy trade-offs",
        subtitle = paste0(
            "CAST's causally-informed variable screening substantially reduces feature dimensionality ",
            "while maintaining or improving SDM predictive accuracy across all model families"
        ),
        theme = theme(
            plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
            plot.subtitle = element_text(
                face = "italic", size = 9.5, hjust = 0.5,
                color = "grey40"
            )
        )
    )

ggsave(file.path(fig_dir, "fig2_screening_parsimony.png"),
    fig2,
    width = 13, height = 10, dpi = 300, bg = "white"
)
cat("✓ Saved fig2_screening_parsimony.png\n")
