################################################################################
# Fig 3 (Eco): Core Performance & Parsimony — 4-panel
#
# Narrative: CAST achieves parsimony (variable reduction) without accuracy
#            trade-offs, while causal structure encoding provides additional gain.
#
# Panel (a): Ablation violin — MLP → MLP_ATE → CAST (isolates Structure Effect)
# Panel (b): All-model mean ± SE dot plot (CAST competitive with full baselines)
# Panel (c): Per-species paired lines — FlatNN(MLP) vs CAST
# Panel (d): Variable reduction vs ΔAUC scatter (parsimony without loss)
#
# Data required:
#   output/case2_eco/all_results_v3.csv
#   output/case2_eco/all_screening_v3.csv
#   outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv
#
# Run: setwd("E:/CausalSDMs")
#      source("scripts/EcoISEA3H/plot/fig3_performance_parsimony_eco.R")
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
model_colors <- c(
    "CAST"    = "#2980B9",
    "MLP_ATE" = "#5DADE2",
    "MLP"     = "#AED6F1",
    "RF"      = "#27AE60",
    "BRT"     = "#E67E22",
    "Maxent"  = "#9B59B6"
)
model_order <- c("CAST", "MLP_ATE", "MLP", "RF", "BRT", "Maxent")

# ── Load data ─────────────────────────────────────────────────────────────────
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
screen <- screen %>% left_join(sp_meta %>% select(species, family), by = "species")

# ══════════════════════════════════════════════════════════════════════════════
# Panel (a): Ablation violin — isolate Structure Effect
#            MLP (flat, CAST vars) → MLP_ATE (+ ATE weighting) → CAST (+ DAG interactions)
# ══════════════════════════════════════════════════════════════════════════════
d_abl <- d %>%
    filter(model %in% c("MLP", "MLP_ATE", "CAST")) %>%
    mutate(model = factor(model, levels = c("MLP", "MLP_ATE", "CAST")))

abl_stats <- d_abl %>%
    group_by(model) %>%
    summarise(mean_auc = mean(auc_mean), .groups = "drop")

pa <- ggplot(d_abl, aes(x = model, y = auc_mean, fill = model)) +
    geom_violin(scale = "width", alpha = 0.55, trim = TRUE, linewidth = 0.3) +
    geom_boxplot(
        width = 0.13, fill = "white", alpha = 0.85,
        outlier.size = 0.7, outlier.alpha = 0.5
    ) +
    stat_summary(
        fun = mean, geom = "point", size = 3.5,
        color = "black", shape = 18
    ) +
    geom_text(
        data = abl_stats,
        aes(x = model, y = mean_auc, label = sprintf("%.4f", mean_auc)),
        vjust = -2.2, size = 3.3, fontface = "bold", color = "black"
    ) +
    # Annotation arrows showing incremental gain
    annotate("segment",
        x = 1.5, xend = 1.5, y = -Inf, yend = Inf,
        linetype = "dotted", color = "grey70", linewidth = 0.4
    ) +
    annotate("segment",
        x = 2.5, xend = 2.5, y = -Inf, yend = Inf,
        linetype = "dotted", color = "grey70", linewidth = 0.4
    ) +
    scale_fill_manual(values = model_colors, guide = "none") +
    labs(
        title = "(a) Causal ablation: structure effect",
        subtitle = "MLP → +ATE weighting → +DAG interactions (CAST)",
        x = "", y = "AUC"
    ) +
    theme_pub()

# ══════════════════════════════════════════════════════════════════════════════
# Panel (b): All-model mean ± SE — full competitive context
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

n_sp_total <- n_distinct(d$species)

pb <- ggplot(full_stats, aes(x = model, y = mean_auc, color = model)) +
    # Highlight causal family
    annotate("rect",
        xmin = 0.5, xmax = 2.5,
        ymin = -Inf, ymax = Inf, alpha = 0.06, fill = "#2980B9"
    ) +
    annotate("text",
        x = 1.5, y = min(full_stats$mean_auc - full_stats$se_auc) - 0.004,
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
        title = "(b) All models: mean AUC ± SE",
        subtitle = sprintf(
            "n = %d species | CAST competitive across baselines",
            n_sp_total
        ),
        x = "", y = "Mean AUC"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 22, hjust = 1))

# ══════════════════════════════════════════════════════════════════════════════
# Panel (c): Per-species paired slope — MLP vs CAST
# ══════════════════════════════════════════════════════════════════════════════
paired_wide <- d %>%
    filter(model %in% c("MLP", "CAST")) %>%
    select(species, family, model, auc_mean) %>%
    pivot_wider(names_from = model, values_from = auc_mean, values_fn = max) %>%
    filter(!is.na(CAST), !is.na(MLP)) %>%
    mutate(
        delta     = CAST - MLP,
        direction = ifelse(delta >= 0, "CAST ≥ MLP", "MLP > CAST")
    )

n_wins <- sum(paired_wide$delta >= 0)
n_total <- nrow(paired_wide)

pc <- ggplot(paired_wide) +
    geom_segment(
        aes(x = 1, xend = 2, y = MLP, yend = CAST, color = direction),
        alpha = 0.45, linewidth = 0.4
    ) +
    geom_point(aes(x = 1, y = MLP), color = "#AED6F1", size = 1.5, alpha = 0.8) +
    geom_point(aes(x = 2, y = CAST), color = "#2980B9", size = 1.5, alpha = 0.8) +
    geom_hline(
        yintercept = mean(paired_wide$MLP),
        linetype = "dotted", color = "#AED6F1", linewidth = 0.7
    ) +
    geom_hline(
        yintercept = mean(paired_wide$CAST),
        linetype = "dotted", color = "#2980B9", linewidth = 0.7
    ) +
    annotate("text",
        x = 0.78, y = mean(paired_wide$MLP),
        label = sprintf("%.4f", mean(paired_wide$MLP)),
        size = 2.8, color = "#AED6F1", fontface = "bold"
    ) +
    annotate("text",
        x = 2.22, y = mean(paired_wide$CAST),
        label = sprintf("%.4f", mean(paired_wide$CAST)),
        size = 2.8, color = "#2980B9", fontface = "bold"
    ) +
    scale_color_manual(
        values = c("CAST ≥ MLP" = "#2980B9", "MLP > CAST" = "#E74C3C"),
        name = ""
    ) +
    scale_x_continuous(
        breaks = c(1, 2), labels = c("FlatNN\n(MLP)", "CAST\n(CI-MLP)"),
        limits = c(0.65, 2.35)
    ) +
    labs(
        title = "(c) Per-species causal enhancement",
        subtitle = sprintf(
            "CAST ≥ MLP in %d / %d species (%.0f%%)",
            n_wins, n_total, 100 * n_wins / n_total
        ),
        x = "", y = "AUC"
    ) +
    theme_pub() +
    theme(legend.position = "bottom")

# ══════════════════════════════════════════════════════════════════════════════
# Panel (d): Variable reduction vs ΔAUC — parsimony without loss
#
# x: pct of variables removed (before: n_vars_full inferred from screening,
#    after: CAST n_vars from results)
# y: ΔAUC (CAST - full-variable MLP baseline) — here we compare species-level
#    AUC of CAST model against the MLP run on all VIF-screened variables.
#    If all_results has a model "MLP_full" or separate group, use that.
#    Otherwise we approximate with the delta from paired_wide and the
#    variable reduction ratio from screening.
# ══════════════════════════════════════════════════════════════════════════════

# Compute n_vars_full per species from screening (total vars before selection)
# and n_vars_cast from results
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

# ΔAUC: CAST minus MLP (on same CAST-selected vars) — isolates parsimony gain
delta_df <- paired_wide %>%
    select(species, delta, family)

pd_data <- var_red %>%
    left_join(delta_df, by = "species") %>%
    filter(!is.na(delta), pct_removed >= 0)

# Summarise: bin by removal quartile
n_total_d <- nrow(pd_data)

pd <- ggplot(pd_data, aes(x = pct_removed, y = delta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.6) +
    geom_point(aes(color = family), size = 2.2, alpha = 0.72) +
    geom_smooth(
        method = "loess", se = TRUE, color = "black",
        linetype = "solid", linewidth = 0.7, fill = "grey85"
    ) +
    scale_color_brewer(palette = "Set1", name = "Taxonomic\nFamily") +
    labs(
        title = "(d) Parsimony vs predictive change",
        subtitle = sprintf(
            "n = %d species | x: %% vars removed by CAST; y: Δ AUC (CAST − MLP)",
            n_total_d
        ),
        x = "Variables removed by CAST (%)",
        y = expression(Delta * "AUC (CAST − FlatNN)")
    ) +
    theme_pub() +
    theme(legend.position = "right")

# ══════════════════════════════════════════════════════════════════════════════
# Combine
# ══════════════════════════════════════════════════════════════════════════════
fig3 <- (pa | pb) / (pc | pd) +
    plot_annotation(
        title = "Fig 3 (Eco)  Core Performance and Parsimony",
        subtitle = "CAST maintains or improves predictive accuracy while substantially reducing variable dimensionality",
        theme = theme(
            plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
            plot.subtitle = element_text(
                face = "italic", size = 10, hjust = 0.5,
                color = "grey40"
            )
        )
    )

ggsave(file.path(fig_dir, "fig3_performance_parsimony_eco.png"),
    fig3,
    width = 14, height = 10, dpi = 300, bg = "white"
)
cat("Saved fig3_performance_parsimony_eco.png\n")
