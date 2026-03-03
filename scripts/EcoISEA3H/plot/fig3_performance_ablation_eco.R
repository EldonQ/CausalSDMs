################################################################################
# Fig 3 (Eco): Causal Enhancement — Ablation-Focused (4-panel)
#
# Adapted for EcoISEA3H dataset.
# Narrative: "因果信息注入不损害预测，同时提供精简变量集"
#   (a) Ablation violin: CAST vs MLP_ATE vs MLP ONLY
#   (b) Paired species lines: MLP → CAST per species (individual gain/loss)
#   (c) Variable reduction vs CAST/MLP AUC ratio (colored by Family)
#   (d) Compact full-model context dot plot (secondary, for transparency)
#
# Requires: all_results_v3.csv, all_screening_v3.csv, CAST_Species_Summary.csv
# Run: setwd("E:/CausalSDMs"); source("scripts/case2_eco/plot/fig3_performance_ablation_eco.R")
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
            axis.title       = element_text(face = "bold"),
            plot.title       = element_text(face = "bold", hjust = 0.5),
            plot.subtitle    = element_text(hjust = 0.5, color = "grey40")
        )
}

# ---- Load data ----
d <- read.csv("output/case2_eco/all_results_v3.csv", stringsAsFactors = FALSE) %>%
    filter(!is.na(auc_mean))
screen <- read.csv("output/case2_eco/all_screening_v3.csv", stringsAsFactors = FALSE)
sp_meta <- read.csv("outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv", stringsAsFactors = FALSE) %>%
    mutate(species = gsub(" ", "_", species))

d <- d %>% left_join(sp_meta %>% select(species, family), by = "species")
screen <- screen %>% left_join(sp_meta %>% select(species, family), by = "species")

ablation_colors <- c(
    "CAST" = "#2980B9", "MLP_ATE" = "#5DADE2", "MLP" = "#AED6F1"
)
all_model_colors <- c(
    "CAST" = "#2980B9", "MLP_ATE" = "#5DADE2", "MLP" = "#AED6F1",
    "RF" = "#27AE60", "BRT" = "#E67E22", "Maxent" = "#9B59B6"
)

# ══════════════════════════════════════════════════════════════
# (a) Ablation violin: ONLY MLP → MLP_ATE → CAST
# ══════════════════════════════════════════════════════════════
d_abl <- d %>%
    filter(model %in% c("MLP", "MLP_ATE", "CAST")) %>%
    mutate(model = factor(model, levels = c("MLP", "MLP_ATE", "CAST")))

abl_stats <- d_abl %>%
    group_by(model) %>%
    summarise(mean_auc = mean(auc_mean), .groups = "drop")

pa <- ggplot(d_abl, aes(x = model, y = auc_mean, fill = model)) +
    geom_violin(scale = "width", alpha = 0.6, trim = TRUE) +
    geom_boxplot(width = 0.15, fill = "white", alpha = 0.8, outlier.size = 0.8) +
    stat_summary(fun = mean, geom = "point", size = 3, color = "black", shape = 18) +
    geom_text(
        data = abl_stats, aes(x = model, y = mean_auc, label = sprintf("%.4f", mean_auc)),
        vjust = -1.5, size = 3.5, fontface = "bold", color = "black"
    ) +
    scale_fill_manual(values = ablation_colors, guide = "none") +
    labs(
        title = "(a) Causal ablation: predictive parity",
        subtitle = "MLP → + ATE weighting → + DAG structure (CAST)",
        x = "", y = "AUC"
    ) +
    theme_pub()

# ══════════════════════════════════════════════════════════════
# (b) Paired species lines: MLP → CAST
# ══════════════════════════════════════════════════════════════
paired_wide <- d %>%
    filter(model %in% c("MLP", "CAST")) %>%
    select(species, model, auc_mean) %>%
    pivot_wider(names_from = model, values_from = auc_mean, values_fn = max) %>%
    filter(!is.na(CAST), !is.na(MLP)) %>%
    mutate(
        delta     = CAST - MLP,
        direction = ifelse(delta >= 0, "CAST wins", "MLP wins")
    )

n_cast_wins <- sum(paired_wide$delta >= 0)
n_total <- nrow(paired_wide)

pb <- ggplot(paired_wide) +
    geom_segment(
        aes(x = 1, xend = 2, y = MLP, yend = CAST, color = direction),
        alpha = 0.5, linewidth = 0.5
    ) +
    geom_point(aes(x = 1, y = MLP), color = "#AED6F1", size = 1.8) +
    geom_point(aes(x = 2, y = CAST), color = "#2980B9", size = 1.8) +
    geom_hline(yintercept = mean(paired_wide$MLP), linetype = "dotted", color = "#AED6F1", linewidth = 0.6) +
    geom_hline(yintercept = mean(paired_wide$CAST), linetype = "dotted", color = "#2980B9", linewidth = 0.6) +
    scale_color_manual(
        values = c("CAST wins" = "#2980B9", "MLP wins" = "#E74C3C"),
        name = ""
    ) +
    scale_x_continuous(breaks = c(1, 2), labels = c("MLP", "CAST"), limits = c(0.7, 2.3)) +
    labs(
        title = "(b) Per-species causal enhancement",
        subtitle = sprintf("CAST ≥ MLP in %d/%d species (%.0f%%)", n_cast_wins, n_total, 100 * n_cast_wins / n_total),
        x = "", y = "AUC"
    ) +
    theme_pub() +
    theme(legend.position = "bottom")

# ══════════════════════════════════════════════════════════════
# (c) Interaction Complexity vs CAST/MLP AUC ratio
# ══════════════════════════════════════════════════════════════
cast_mlp <- d %>%
    filter(model %in% c("CAST", "MLP")) %>%
    select(species, family, model, auc_mean, n_interactions) %>%
    pivot_wider(names_from = model, values_from = c(auc_mean, n_interactions), names_glue = "{model}_{.value}", values_fn = max) %>%
    filter(!is.na(CAST_auc_mean), !is.na(MLP_auc_mean)) %>%
    mutate(auc_ratio = CAST_auc_mean / MLP_auc_mean)

mean_ratio <- mean(cast_mlp$auc_ratio, na.rm = TRUE)
mean_interactions <- mean(cast_mlp$CAST_n_interactions, na.rm = TRUE)

pc <- ggplot(cast_mlp, aes(x = CAST_n_interactions, y = auc_ratio)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40", linewidth = 0.6) +
    geom_point(aes(color = family), size = 2.5, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 0.5) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 5)) +
    scale_y_continuous(limits = c(min(cast_mlp$auc_ratio, na.rm = TRUE) * 0.99, max(max(cast_mlp$auc_ratio, na.rm = TRUE) * 1.01, 1.01))) +
    scale_color_brewer(palette = "Set1", name = "Family") +
    labs(
        title = "(c) Sparse interactions, competitive gain",
        subtitle = sprintf(
            "Mean: %.1f DAG interactions, mean AUC ratio = %.3f",
            mean_interactions, mean_ratio
        ),
        x = "Number of DAG-guided interactive features",
        y = "CAST AUC / MLP AUC"
    ) +
    theme_pub() +
    theme(legend.position = "right")

# ══════════════════════════════════════════════════════════════
# (d) Full-model context: compact mean ± SE dot plot
# ══════════════════════════════════════════════════════════════
model_order <- c("CAST", "MLP_ATE", "MLP", "RF", "BRT", "Maxent")
full_stats <- d %>%
    mutate(model = factor(model, levels = model_order)) %>%
    group_by(model) %>%
    summarise(
        mean_auc = mean(auc_mean, na.rm = TRUE),
        se_auc   = sd(auc_mean, na.rm = TRUE) / sqrt(n()),
        n        = n(),
        .groups  = "drop"
    )

pd <- ggplot(full_stats, aes(x = model, y = mean_auc, color = model)) +
    geom_pointrange(
        aes(ymin = mean_auc - se_auc, ymax = mean_auc + se_auc),
        size = 0.8, linewidth = 0.8, fatten = 3
    ) +
    geom_text(aes(label = sprintf("%.3f", mean_auc)), vjust = -1.2, size = 3, fontface = "bold") +
    scale_color_manual(values = all_model_colors, guide = "none") +
    annotate("rect",
        xmin = 0.5, xmax = 3.5, ymin = -Inf, ymax = Inf,
        alpha = 0.04, fill = "#2980B9"
    ) +
    annotate("text",
        x = 2, y = min(full_stats$mean_auc - full_stats$se_auc) - 0.005,
        label = "CAST family\n(causal components)",
        fontface = "italic", size = 2.8, color = "#2980B9"
    ) +
    labs(
        title = "(d) All models: mean AUC ± SE",
        subtitle = sprintf("n = %d species | CAST competitive with baselines", nrow(paired_wide)),
        x = "", y = "Mean AUC"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

# ══════════════════════════════════════════════════════════════
# Combine
# ══════════════════════════════════════════════════════════════
fig3 <- (pa | pb) / (pc | pd) +
    plot_annotation(
        title = "Fig 3 (Eco)  Causal Enhancement: Predictive Parity with Parsimony",
        theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    )

ggsave(file.path(fig_dir, "fig3_performance_ablation_eco.png"),
    fig3,
    width = 14, height = 10, dpi = 300, bg = "white"
)
cat("Saved fig3_performance_ablation_eco.png\n")
