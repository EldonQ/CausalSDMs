################################################################################
# Fig 4 (Eco): DAG Density vs CAST Gain — Mechanistic Insight  [MAIN FIGURE]
#
# Narrative: The predictive advantage of CAST's causal structure encoding
#   (CI-MLP over FlatNN) is mediated by the informativeness (sparsity) of the
#   learned causal graph. Dense, near-saturated DAGs yield little additional
#   structure signal; sparse, hierarchically clear DAGs enable targeted
#   feature engineering that FlatNN cannot replicate.
#
# Panel (a): Species-level DAG density × ΔAUC scatter + regression
#            (coloured by taxonomic family)
# Panel (b): Family-level summary — mean DAG density vs mean ΔAUC
#            (bubble size = n species; labelled)
# Panel (c): Distribution of ΔAUC by DAG density tercile (violin/boxplot)
#
# Data required:
#   output/case2_eco/all_results_v3.csv
#   outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv
#
# Run: setwd("E:/CausalSDMs")
#      source("scripts/EcoISEA3H/plot/fig4_dag_density_mechanism_eco.R")
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
d <- read.csv("output/case2_eco/all_results_v3.csv", stringsAsFactors = FALSE) %>%
    filter(!is.na(auc_mean))

sp_meta <- read.csv(
    "outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv",
    stringsAsFactors = FALSE
) %>%
    mutate(species = gsub(" ", "_", species))

# ── Build species-level ΔAUC table ───────────────────────────────────────────
# NOTE: dag_density is a species-level property stored identically across model
# rows. Extract from CAST rows BEFORE pivot_wider to avoid list-column collapse.
dag_density_ref <- d %>%
    filter(model == "CAST", !is.na(dag_density)) %>%
    select(species, dag_density) %>%
    distinct(species, .keep_all = TRUE) %>%
    mutate(dag_density = as.numeric(dag_density))

abl_wide <- d %>%
    filter(model %in% c("MLP", "CAST")) %>%
    select(species, model, auc_mean) %>% # dag_density NOT in pivot
    pivot_wider(names_from = model, values_from = auc_mean, values_fn = max) %>%
    filter(!is.na(CAST), !is.na(MLP)) %>%
    mutate(delta_auc = CAST - MLP) %>%
    left_join(dag_density_ref, by = "species") %>%
    left_join(sp_meta %>% select(species, family, category), by = "species") %>%
    mutate(family = ifelse(is.na(family), "Unknown", family)) %>%
    filter(!is.na(dag_density))

n_sp <- nrow(abl_wide)

# ── Diagnostics ───────────────────────────────────────────────────────────────
cat(sprintf(
    "[Fig4 diag] n_species=%d  dag_density: min=%.4f  max=%.4f  sd=%.6f  unique=%d\n",
    n_sp,
    min(abl_wide$dag_density, na.rm = TRUE),
    max(abl_wide$dag_density, na.rm = TRUE),
    sd(abl_wide$dag_density, na.rm = TRUE),
    length(unique(abl_wide$dag_density))
))

# ── Linear regression (guarded against zero-variance predictor) ───────────────
dag_var <- var(abl_wide$dag_density, na.rm = TRUE)

if (!is.na(dag_var) && dag_var > 1e-10) {
    lm_fit <- lm(delta_auc ~ dag_density, data = abl_wide)
    lm_r2 <- summary(lm_fit)$r.squared
    lm_slope <- coef(lm_fit)["dag_density"]
    lm_pval <- summary(lm_fit)$coefficients["dag_density", "Pr(>|t|)"]
    pval_lab <- ifelse(lm_pval < 0.001, "p < 0.001",
        sprintf("p = %.3f", lm_pval)
    )
    stats_lab <- sprintf("R\u00b2 = %.3f\nslope = %.4f\n%s", lm_r2, lm_slope, pval_lab)
    has_lm <- TRUE
} else {
    warning(paste(
        "[Fig4] dag_density has near-zero variance across all species.",
        "The scatter plot will render without a regression line.",
        "This likely means all species share the same dag_density value in",
        "all_results_v3.csv — check that per-species DAG density was saved correctly."
    ))
    lm_r2 <- NA_real_
    lm_slope <- NA_real_
    stats_lab <- "(dag_density has\nno cross-species variation)"
    has_lm <- FALSE
}

# ══════════════════════════════════════════════════════════════════════════════
# Panel (a): Species-level scatter — the main evidence
# ══════════════════════════════════════════════════════════════════════════════
pa <- ggplot(abl_wide, aes(x = dag_density, y = delta_auc)) +
    geom_hline(
        yintercept = 0, linetype = "dotted",
        color = "grey50", linewidth = 0.6
    ) +
    geom_point(aes(color = family), size = 2.4, alpha = 0.75) +
    {
        if (has_lm) {
            list(
                geom_smooth(
                    method = "lm", se = TRUE,
                    color = "black", linetype = "dashed",
                    linewidth = 0.8, fill = "grey80", alpha = 0.35
                )
            )
        }
    } +
    annotate("text",
        x = max(abl_wide$dag_density, na.rm = TRUE) * 0.98,
        y = max(abl_wide$delta_auc, na.rm = TRUE) * 0.95,
        label = stats_lab,
        hjust = 1, vjust = 1, size = 3.2,
        fontface = "italic", color = "black"
    ) +
    scale_color_brewer(palette = "Set1", name = "Taxonomic\nFamily") +
    labs(
        title = "(a) DAG density moderates CAST structural gain",
        subtitle = sprintf(
            "n = %d species | Negative slope: sparser DAG \u2192 larger CAST advantage",
            n_sp
        ),
        x = "DAG density (proportion of possible edges retained)",
        y = expression(Delta * "AUC (CAST \u2212 FlatNN)")
    ) +
    theme_pub() +
    theme(legend.position = "right")

# ══════════════════════════════════════════════════════════════════════════════
# Panel (b): Family-level summary — mean DAG density vs mean ΔAUC
# ══════════════════════════════════════════════════════════════════════════════
fam_summary <- abl_wide %>%
    group_by(family) %>%
    summarise(
        mean_dag   = mean(dag_density, na.rm = TRUE),
        mean_delta = mean(delta_auc, na.rm = TRUE),
        n          = n(),
        .groups    = "drop"
    )

pb <- ggplot(fam_summary, aes(x = mean_dag, y = mean_delta, size = n)) +
    geom_hline(
        yintercept = 0, linetype = "dotted",
        color = "grey50", linewidth = 0.6
    ) +
    geom_point(aes(color = family), alpha = 0.85) +
    geom_text(aes(label = family), vjust = -1.0, size = 3.0, fontface = "italic") +
    scale_size_continuous(name = "# species", range = c(3, 11)) +
    scale_color_brewer(palette = "Set1", guide = "none") +
    labs(
        title = "(b) Family-level: sparser graphs yield larger gains",
        subtitle = "Mean DAG density vs mean \u0394AUC; bubble size = species count",
        x = "Mean DAG density",
        y = expression("Mean " * Delta * "AUC (CAST \u2212 FlatNN)")
    ) +
    theme_pub() +
    theme(legend.position = "right")

# ══════════════════════════════════════════════════════════════════════════════
# Panel (c): ΔAUC distribution by DAG-density tercile
# ══════════════════════════════════════════════════════════════════════════════
q_lo <- quantile(abl_wide$dag_density, 1 / 3, na.rm = TRUE)
q_hi <- quantile(abl_wide$dag_density, 2 / 3, na.rm = TRUE)

abl_wide <- abl_wide %>%
    mutate(
        density_tercile = case_when(
            dag_density <= q_lo ~ "Low density\n(sparse, informative)",
            dag_density <= q_hi ~ "Mid density",
            TRUE ~ "High density\n(dense, redundant)"
        ),
        density_tercile = factor(
            density_tercile,
            levels = c(
                "Low density\n(sparse, informative)",
                "Mid density",
                "High density\n(dense, redundant)"
            )
        )
    )

tercile_means <- abl_wide %>%
    group_by(density_tercile) %>%
    summarise(m = mean(delta_auc), .groups = "drop")

pc <- ggplot(abl_wide, aes(
    x = density_tercile, y = delta_auc,
    fill = density_tercile
)) +
    geom_hline(
        yintercept = 0, linetype = "dotted",
        color = "grey50", linewidth = 0.6
    ) +
    geom_violin(scale = "width", alpha = 0.5, trim = TRUE, linewidth = 0.3) +
    geom_boxplot(
        width = 0.14, fill = "white", alpha = 0.85,
        outlier.size = 0.8, outlier.alpha = 0.5
    ) +
    stat_summary(
        fun = mean, geom = "point", shape = 18,
        size = 3.5, color = "black"
    ) +
    geom_text(
        data = tercile_means,
        aes(
            x = density_tercile, y = m,
            label = sprintf("%.4f", m)
        ),
        vjust = -2.0, size = 3.3, fontface = "bold",
        color = "black", inherit.aes = FALSE
    ) +
    scale_fill_manual(
        values = c(
            "Low density\n(sparse, informative)" = "#2980B9",
            "Mid density" = "#85C1E9",
            "High density\n(dense, redundant)" = "#D5E8F3"
        ),
        guide = "none"
    ) +
    labs(
        title = "(c) \u0394AUC by DAG density tercile",
        subtitle = "Species grouped by causal graph sparsity | Mean shown as diamond",
        x = "", y = expression(Delta * "AUC (CAST \u2212 FlatNN)")
    ) +
    theme_pub()

# ══════════════════════════════════════════════════════════════════════════════
# Combine: wide top panel + two bottom panels
# ══════════════════════════════════════════════════════════════════════════════
fig4 <- pa / (pb | pc) +
    plot_layout(heights = c(1.3, 1)) +
    plot_annotation(
        title = "Fig 4 (Eco)  DAG Informativeness Drives Structural Encoding Gain",
        subtitle = paste0(
            "CAST's causal topology encoding provides the largest predictive edge ",
            "when the environmental causal graph is sparse and hierarchically structured"
        ),
        theme = theme(
            plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
            plot.subtitle = element_text(
                face = "italic", size = 10, hjust = 0.5,
                color = "grey40"
            )
        )
    )

ggsave(file.path(fig_dir, "fig4_dag_density_mechanism_eco.png"),
    fig4,
    width = 13, height = 12, dpi = 300, bg = "white"
)
cat("Saved fig4_dag_density_mechanism_eco.png\n")
