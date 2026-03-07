################################################################################
# Fig 3: Causal Topology Encoding Provides Systematic Predictive Gains
#
# Scientific Question Q2:
#   Under identical variable sets and identical network architecture,
#   does encoding causal topology into feature space improve prediction?
#
# Panel (a): Ablation violin — MLP → MLP_ATE → CAST (Structure Effect)
# Panel (b): CI-MLP vs FlatNN per-species AUC scatter (45° diagonal)
# Panel (c): Paired slope plot — MLP vs CAST per species
#
# Data required:
#   output/case2_eco/all_results_v3.csv
#   outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv
#
# Run: setwd("E:/CausalSDMs")
#      source("scripts/EcoISEA3H/plot/fig3_topology_encoding_gain.R")
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
    "MLP"     = "#AED6F1"
)

# ── Load data ────────────────────────────────────────────────────────────────
d <- read.csv("output/case2_eco/all_results_v3.csv", stringsAsFactors = FALSE) %>%
    filter(!is.na(auc_mean))

sp_meta <- read.csv(
    "outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv",
    stringsAsFactors = FALSE
) %>%
    mutate(species = gsub(" ", "_", species))

d <- d %>% left_join(sp_meta %>% select(species, family), by = "species")

# ══════════════════════════════════════════════════════════════════════════════
# Panel (a): Ablation violin — isolate Structure Effect step by step
#   MLP (bare features) → MLP_ATE (+ ATE weighting) → CAST (+ DAG interactions)
# ══════════════════════════════════════════════════════════════════════════════
d_abl <- d %>%
    filter(model %in% c("MLP", "MLP_ATE", "CAST")) %>%
    mutate(model = factor(model, levels = c("MLP", "MLP_ATE", "CAST")))

abl_stats <- d_abl %>%
    group_by(model) %>%
    summarise(mean_auc = mean(auc_mean), .groups = "drop")

# Wilcoxon paired tests for statistical rigour
wide_abl <- d %>%
    filter(model %in% c("MLP", "MLP_ATE", "CAST")) %>%
    select(species, model, auc_mean) %>%
    pivot_wider(names_from = model, values_from = auc_mean, values_fn = max) %>%
    filter(!is.na(CAST), !is.na(MLP), !is.na(MLP_ATE))

p_mlp_ate <- tryCatch(
    wilcox.test(wide_abl$MLP_ATE, wide_abl$MLP, paired = TRUE)$p.value,
    error = function(e) NA_real_
)
p_ate_cast <- tryCatch(
    wilcox.test(wide_abl$CAST, wide_abl$MLP_ATE, paired = TRUE)$p.value,
    error = function(e) NA_real_
)

# Format p-value labels
fmt_p <- function(p) {
    if (is.na(p)) {
        return("")
    }
    if (p < 0.001) {
        return("***")
    }
    if (p < 0.01) {
        return("**")
    }
    if (p < 0.05) {
        return("*")
    }
    return("ns")
}

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
    # Show step-wise significance
    annotate("segment",
        x = 1, xend = 2, y = max(d_abl$auc_mean, na.rm = TRUE) + 0.005,
        yend = max(d_abl$auc_mean, na.rm = TRUE) + 0.005,
        linewidth = 0.4, color = "grey30"
    ) +
    annotate("text",
        x = 1.5, y = max(d_abl$auc_mean, na.rm = TRUE) + 0.007,
        label = fmt_p(p_mlp_ate), size = 3.5, color = "grey30"
    ) +
    annotate("segment",
        x = 2, xend = 3, y = max(d_abl$auc_mean, na.rm = TRUE) + 0.012,
        yend = max(d_abl$auc_mean, na.rm = TRUE) + 0.012,
        linewidth = 0.4, color = "grey30"
    ) +
    annotate("text",
        x = 2.5, y = max(d_abl$auc_mean, na.rm = TRUE) + 0.014,
        label = fmt_p(p_ate_cast), size = 3.5, color = "grey30"
    ) +
    scale_fill_manual(values = model_colors, guide = "none") +
    labs(
        title = "(a) Causal ablation: stepwise structure effect",
        subtitle = "FlatNN → +ATE weighting → +DAG interactions (CI-MLP)",
        x = "", y = "AUC"
    ) +
    theme_pub()

# ══════════════════════════════════════════════════════════════════════════════
# Panel (b): CI-MLP vs FlatNN — per-species scatter (45° line)
# ══════════════════════════════════════════════════════════════════════════════
scatter_wide <- d %>%
    filter(model %in% c("MLP", "CAST")) %>%
    select(species, family, model, auc_mean) %>%
    pivot_wider(names_from = model, values_from = auc_mean, values_fn = max) %>%
    filter(!is.na(CAST), !is.na(MLP))

# Extract DAG density for point sizing
dag_ref <- d %>%
    filter(model == "CAST", !is.na(dag_density)) %>%
    select(species, dag_density) %>%
    distinct(species, .keep_all = TRUE) %>%
    mutate(dag_density = as.numeric(dag_density))

scatter_wide <- scatter_wide %>%
    left_join(dag_ref, by = "species")

n_above <- sum(scatter_wide$CAST > scatter_wide$MLP)
n_scatter <- nrow(scatter_wide)

pb <- ggplot(scatter_wide, aes(x = MLP, y = CAST)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = family, size = dag_density), alpha = 0.7) +
    scale_color_brewer(palette = "Set1", name = "Family") +
    scale_size_continuous(
        name = "DAG\ndensity", range = c(1.5, 5),
        guide = guide_legend(override.aes = list(alpha = 0.8))
    ) +
    labs(
        title = "(b) CI-MLP vs FlatNN: per-species AUC",
        subtitle = sprintf(
            "%d / %d species above diagonal (%.0f%%) — systematic CI-MLP advantage",
            n_above, n_scatter, 100 * n_above / n_scatter
        ),
        x = "FlatNN AUC", y = "CI-MLP (CAST) AUC"
    ) +
    theme_pub() +
    theme(legend.position = "right")

# ══════════════════════════════════════════════════════════════════════════════
# Panel (c): Paired slope — MLP → CAST per species
# ══════════════════════════════════════════════════════════════════════════════
paired_wide <- scatter_wide %>%
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
        breaks = c(1, 2), labels = c("FlatNN\n(MLP)", "CI-MLP\n(CAST)"),
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
# Combine
# ══════════════════════════════════════════════════════════════════════════════
fig3 <- pa / (pb | pc) +
    plot_layout(heights = c(1, 1.1)) +
    plot_annotation(
        title = "Fig 3  Encoding causal topology provides systematic predictive gains",
        subtitle = paste0(
            "Under identical CAST-screened variables and identical 4-layer MLP architecture, ",
            "CI-MLP with ATE-weighted inputs and DAG-guided interaction features ",
            "systematically outperforms the structure-unaware FlatNN"
        ),
        theme = theme(
            plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
            plot.subtitle = element_text(
                face = "italic", size = 9.5, hjust = 0.5,
                color = "grey40"
            )
        )
    )

ggsave(file.path(fig_dir, "fig3_topology_encoding_gain.png"),
    fig3,
    width = 13, height = 10, dpi = 300, bg = "white"
)
cat("✓ Saved fig3_topology_encoding_gain.png\n")
