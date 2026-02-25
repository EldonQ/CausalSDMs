################################################################################
# 03_publication_figures.R
# CAST v3 Publication Figures — Multi-region CI-MLP benchmark
#
# Reads: output/case2/all_results_v3.csv (+ dag_info, screening, etc.)
# Saves: figures/case2/
#
# Figures:
#   1. Three-group AUC/TSS comparison (across all regions)
#   2. CI-MLP vs FlatNN_cast scatter (per species)
#   3. Screening effect (A→B) per algorithm
#   4. DAG density vs structure advantage (NEW — key diagnostic)
#   5. Per-region AUC comparison
#   6. Variable reduction & interaction feature count
#   7. Summary table
################################################################################

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

library(tidyverse)
library(patchwork)
library(scales)

dir.create("figures/case2", recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Load data
# ==============================================================================
results_file <- "output/case2/all_results_v3.csv"
if (!file.exists(results_file)) {
    stop("Results file not found. Run 02_multi_species_experiment.R first.")
}

all_results <- read.csv(results_file, stringsAsFactors = FALSE)
dag_info <- tryCatch(read.csv("output/case2/all_dag_info_v3.csv", stringsAsFactors = FALSE),
    error = function(e) data.frame()
)

cat(sprintf(
    "Loaded %d result rows for %d species across %d regions\n",
    nrow(all_results), length(unique(all_results$species)),
    length(unique(all_results$region))
))

# ==============================================================================
# Theming & Colors
# ==============================================================================
theme_pub <- function(base_size = 11) {
    theme_minimal(base_size = base_size) +
        theme(
            plot.title = element_text(face = "bold", size = base_size + 1),
            plot.subtitle = element_text(color = "grey40", size = base_size - 1),
            panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold"),
            legend.position = "bottom"
        )
}

# Group colors
pal <- list(
    full = "#95A5A6", # A: Full variables (grey)
    cast = "#3498DB", # B: CAST screened (blue)
    ci   = "#E74C3C" # C: CI-MLP (red)
)

group_colors <- c(
    "A: Full variables" = pal$full,
    "B: CAST screened" = pal$cast,
    "C: CI-MLP" = pal$ci
)

region_colors <- c(
    "AWT" = "#E74C3C", "CAN" = "#3498DB", "NSW" = "#2ECC71",
    "NZ"  = "#9B59B6", "SA"  = "#F39C12", "SWI" = "#1ABC9C"
)

# ==============================================================================
# Classify models into groups
# ==============================================================================
all_results <- all_results %>%
    mutate(
        group = case_when(
            model == "CI_MLP" ~ "C: CI-MLP",
            var_set == "full" ~ "A: Full variables",
            TRUE ~ "B: CAST screened"
        ),
        model_label = case_when(
            model == "CI_MLP" ~ "CI-MLP",
            grepl("_cast$", model) | grepl("_full$", model) ~
                gsub("_(cast|full)$", "", model),
            TRUE ~ model
        )
    )

all_results$group <- factor(all_results$group,
    levels = c("A: Full variables", "B: CAST screened", "C: CI-MLP")
)

# ==============================================================================
# Figure 1: Overall three-group AUC comparison
# ==============================================================================
cat("\nFig 1: Three-group AUC/TSS...\n")

summary_by_model <- all_results %>%
    group_by(model, group, model_label) %>%
    summarise(
        AUC_mean = mean(auc_mean, na.rm = TRUE),
        AUC_se = sd(auc_mean, na.rm = TRUE) / sqrt(n()),
        TSS_mean = mean(tss_mean, na.rm = TRUE),
        TSS_se = sd(tss_mean, na.rm = TRUE) / sqrt(n()),
        n = n(), .groups = "drop"
    ) %>%
    arrange(desc(AUC_mean))

model_order <- c("Maxent", "GAM", "RF", "BRT", "FlatNN", "CI-MLP")
summary_by_model$model_label <- factor(summary_by_model$model_label, levels = model_order)

fig1a <- ggplot(summary_by_model, aes(x = model_label, y = AUC_mean, fill = group)) +
    geom_col(
        position = position_dodge(width = 0.75), width = 0.65,
        alpha = 0.9, color = "white", linewidth = 0.3
    ) +
    geom_errorbar(aes(ymin = AUC_mean - AUC_se, ymax = AUC_mean + AUC_se),
        position = position_dodge(width = 0.75), width = 0.2
    ) +
    geom_text(aes(label = sprintf("%.3f", AUC_mean)),
        position = position_dodge(width = 0.75), vjust = -0.8, size = 2.5
    ) +
    scale_fill_manual(values = group_colors, name = "Group") +
    labs(title = "(a) AUC", x = "", y = "Mean AUC") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

fig1b <- ggplot(summary_by_model, aes(x = model_label, y = TSS_mean, fill = group)) +
    geom_col(
        position = position_dodge(width = 0.75), width = 0.65,
        alpha = 0.9, color = "white", linewidth = 0.3
    ) +
    geom_errorbar(aes(ymin = TSS_mean - TSS_se, ymax = TSS_mean + TSS_se),
        position = position_dodge(width = 0.75), width = 0.2
    ) +
    scale_fill_manual(values = group_colors, name = "Group") +
    labs(title = "(b) TSS", x = "", y = "Mean TSS") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

fig1 <- fig1a + fig1b + plot_layout(guides = "collect") &
    theme(legend.position = "bottom")
ggsave("figures/case2/fig1_three_group_performance.pdf", fig1,
    width = 12, height = 5.5, dpi = 300
)
cat("  → figures/case2/fig1_three_group_performance.pdf\n")

# ==============================================================================
# Figure 2: CI-MLP vs FlatNN_cast scatter (structure effect)
# ==============================================================================
cat("Fig 2: CI-MLP vs FlatNN scatter...\n")

wide_ci <- all_results %>%
    filter(model %in% c("CI_MLP", "FlatNN_cast")) %>%
    select(region, species, model, auc_mean, dag_density) %>%
    pivot_wider(names_from = model, values_from = auc_mean) %>%
    filter(!is.na(CI_MLP), !is.na(FlatNN_cast)) %>%
    mutate(delta = CI_MLP - FlatNN_cast)

fig2 <- ggplot(wide_ci, aes(x = FlatNN_cast, y = CI_MLP, color = region)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(size = dag_density), alpha = 0.7) +
    scale_color_manual(values = region_colors) +
    scale_size_continuous(name = "DAG density", range = c(1.5, 5)) +
    labs(
        title = "CI-MLP vs FlatNN (same CAST variables)",
        subtitle = sprintf(
            "Points above line = CI-MLP wins (%d/%d species, %.0f%%)",
            sum(wide_ci$delta > 0), nrow(wide_ci), mean(wide_ci$delta > 0) * 100
        ),
        x = "FlatNN_cast AUC", y = "CI-MLP AUC"
    ) +
    coord_equal() +
    theme_pub()
ggsave("figures/case2/fig2_cimlp_vs_flatnn.pdf", fig2, width = 7, height = 6, dpi = 300)
cat("  → figures/case2/fig2_cimlp_vs_flatnn.pdf\n")

# ==============================================================================
# Figure 3: Screening effect (A→B)
# ==============================================================================
cat("Fig 3: Screening effect...\n")

base_models_ab <- all_results %>%
    filter(
        model_label %in% c("Maxent", "GAM", "RF", "BRT", "FlatNN"),
        group != "C: CI-MLP"
    ) %>%
    select(region, species, model_label, group, auc_mean) %>%
    pivot_wider(names_from = group, values_from = auc_mean) %>%
    filter(!is.na(`A: Full variables`), !is.na(`B: CAST screened`)) %>%
    mutate(delta_screening = `B: CAST screened` - `A: Full variables`)

fig3 <- ggplot(base_models_ab, aes(x = model_label, y = delta_screening * 100)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_boxplot(fill = pal$cast, alpha = 0.3, outlier.alpha = 0) +
    geom_jitter(aes(color = region), width = 0.2, size = 2, alpha = 0.6) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
    scale_color_manual(values = region_colors, name = "Region") +
    labs(
        title = "CAST Screening Effect (A→B)",
        subtitle = "ΔAUC = CAST-screened − Full-variable (per species)",
        x = "Algorithm", y = "ΔAUC (percentage points)"
    ) +
    theme_pub()
ggsave("figures/case2/fig3_screening_effect.pdf", fig3, width = 8, height = 5.5, dpi = 300)
cat("  → figures/case2/fig3_screening_effect.pdf\n")

# ==============================================================================
# Figure 4: DAG density vs structure advantage (★ KEY DIAGNOSTIC)
# ==============================================================================
cat("Fig 4: DAG density vs structure advantage...\n")

if (nrow(wide_ci) > 5) {
    fig4 <- ggplot(wide_ci, aes(x = dag_density, y = delta * 100)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
        geom_point(aes(color = region), size = 3, alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
        scale_color_manual(values = region_colors, name = "Region") +
        labs(
            title = "DAG Density vs CI-MLP Advantage",
            subtitle = "Hypothesis: sparser DAG → stronger CI-MLP benefit",
            x = "DAG Density (fraction of possible edges present)",
            y = "ΔAUC (CI-MLP − FlatNN_cast, percentage points)"
        ) +
        annotate("text",
            x = min(wide_ci$dag_density) + 0.02,
            y = max(wide_ci$delta * 100) * 0.9,
            label = sprintf(
                "r = %.3f, p = %.3f",
                cor(wide_ci$dag_density, wide_ci$delta),
                cor.test(wide_ci$dag_density, wide_ci$delta)$p.value
            ),
            hjust = 0, size = 3.5
        ) +
        theme_pub()
    ggsave("figures/case2/fig4_dag_density_vs_advantage.pdf", fig4,
        width = 8, height = 6, dpi = 300
    )
    cat("  → figures/case2/fig4_dag_density_vs_advantage.pdf\n")
} else {
    cat("  Skipped (insufficient data)\n")
}

# ==============================================================================
# Figure 5: Per-region comparison
# ==============================================================================
cat("Fig 5: Per-region AUC...\n")

region_summary <- all_results %>%
    group_by(region, group) %>%
    summarise(
        AUC_mean = mean(auc_mean, na.rm = TRUE),
        AUC_se = sd(auc_mean, na.rm = TRUE) / sqrt(n()),
        n = n(), .groups = "drop"
    )

fig5 <- ggplot(region_summary, aes(x = region, y = AUC_mean, fill = group)) +
    geom_col(position = position_dodge(0.8), width = 0.7, alpha = 0.9) +
    geom_errorbar(aes(ymin = AUC_mean - AUC_se, ymax = AUC_mean + AUC_se),
        position = position_dodge(0.8), width = 0.2
    ) +
    scale_fill_manual(values = group_colors) +
    labs(
        title = "CAST Performance Across 6 Regions",
        subtitle = "disdat benchmark (Elith et al. 2020): birds, bats, reptiles, plants",
        x = "Region", y = "Mean AUC"
    ) +
    theme_pub()
ggsave("figures/case2/fig5_per_region.pdf", fig5, width = 9, height = 5.5, dpi = 300)
cat("  → figures/case2/fig5_per_region.pdf\n")

# ==============================================================================
# Figure 6: Variable reduction + interaction features
# ==============================================================================
cat("Fig 6: Variable reduction & interactions...\n")

var_info <- all_results %>%
    filter(model %in% c("CI_MLP", "FlatNN_full")) %>%
    select(region, species, model, n_vars, n_interactions, n_features_total) %>%
    pivot_wider(names_from = model, values_from = c(n_vars, n_interactions, n_features_total))

if (nrow(var_info) > 0 && "n_vars_FlatNN_full" %in% names(var_info) &&
    "n_vars_CI_MLP" %in% names(var_info)) {
    var_info <- var_info %>%
        mutate(reduction_pct = (1 - n_vars_CI_MLP / n_vars_FlatNN_full) * 100)

    fig6a <- ggplot(var_info, aes(x = region, y = reduction_pct)) +
        geom_boxplot(fill = pal$cast, alpha = 0.4) +
        geom_jitter(width = 0.15, size = 2, alpha = 0.5, color = pal$cast) +
        labs(
            title = "(a) Variable Reduction (%)",
            subtitle = "Post-VIF → CAST selected",
            x = "Region", y = "Variable reduction (%)"
        ) +
        theme_pub()

    fig6b <- ggplot(var_info, aes(x = region, y = n_interactions_CI_MLP)) +
        geom_boxplot(fill = pal$ci, alpha = 0.4) +
        geom_jitter(width = 0.15, size = 2, alpha = 0.5, color = pal$ci) +
        labs(
            title = "(b) DAG Interaction Features",
            subtitle = "Number of pairwise interactions from DAG edges",
            x = "Region", y = "# interaction features"
        ) +
        theme_pub()

    fig6 <- fig6a + fig6b
    ggsave("figures/case2/fig6_variables_interactions.pdf", fig6,
        width = 12, height = 5, dpi = 300
    )
    cat("  → figures/case2/fig6_variables_interactions.pdf\n")
}

# ==============================================================================
# Figure 7: CAST decomposition — screening + structure effects
# ==============================================================================
cat("Fig 7: CAST advantage decomposition...\n")

# Screening effect: A→B for FlatNN
delta_screening <- all_results %>%
    filter(model %in% c("FlatNN_cast", "FlatNN_full")) %>%
    select(region, species, model, auc_mean) %>%
    pivot_wider(names_from = model, values_from = auc_mean) %>%
    filter(!is.na(FlatNN_cast), !is.na(FlatNN_full)) %>%
    mutate(delta_screening = FlatNN_cast - FlatNN_full)

# Structure effect: B→C for CI-MLP vs FlatNN_cast
delta_structure <- wide_ci

if (nrow(delta_screening) > 0 && nrow(delta_structure) > 0) {
    fig7a <- ggplot(delta_screening, aes(x = region, y = delta_screening * 100)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
        geom_boxplot(fill = pal$cast, alpha = 0.4, outlier.alpha = 0) +
        geom_jitter(width = 0.15, size = 2, alpha = 0.5, color = pal$cast) +
        stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
        labs(
            title = "(a) Screening effect (A→B): FlatNN",
            subtitle = "ΔAUC = FlatNN_cast − FlatNN_full",
            x = "Region", y = "ΔAUC (pp)"
        ) +
        theme_pub()

    fig7b <- ggplot(delta_structure, aes(x = region, y = delta * 100)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
        geom_boxplot(fill = pal$ci, alpha = 0.4, outlier.alpha = 0) +
        geom_jitter(width = 0.15, size = 2, alpha = 0.5, color = pal$ci) +
        stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
        labs(
            title = "(b) Structure effect (B→C): CI-MLP",
            subtitle = "ΔAUC = CI-MLP − FlatNN_cast",
            x = "Region", y = "ΔAUC (pp)"
        ) +
        theme_pub()

    fig7 <- fig7a + fig7b
    ggsave("figures/case2/fig7_cast_decomposition.pdf", fig7,
        width = 12, height = 5.5, dpi = 300
    )
    cat("  → figures/case2/fig7_cast_decomposition.pdf\n")
}

# ==============================================================================
# Summary Table
# ==============================================================================
cat("\nSummary Table...\n")

summary_table <- all_results %>%
    group_by(model, group) %>%
    summarise(
        n_species = n(),
        AUC_mean = sprintf("%.3f", mean(auc_mean, na.rm = TRUE)),
        AUC_sd = sprintf("%.3f", sd(auc_mean, na.rm = TRUE)),
        TSS_mean = sprintf("%.3f", mean(tss_mean, na.rm = TRUE)),
        TSS_sd = sprintf("%.3f", sd(tss_mean, na.rm = TRUE)),
        .groups = "drop"
    ) %>%
    arrange(group, desc(as.numeric(AUC_mean)))

write.csv(summary_table, "figures/case2/summary_table.csv", row.names = FALSE)

cat("\n╔══════════════════════════════════════════════════════════════════════════╗\n")
cat("║  CAST v3 Summary Table (across all regions)                             ║\n")
cat("╠══════════════════════════════════════════════════════════════════════════╣\n")
for (i in 1:nrow(summary_table)) {
    r <- summary_table[i, ]
    mark <- if (grepl("C:", r$group)) " ★" else "  "
    cat(sprintf(
        "║%s %-14s | %-20s | AUC=%s±%s | TSS=%s±%s | n=%3d ║\n",
        mark, r$model, r$group, r$AUC_mean, r$AUC_sd,
        r$TSS_mean, r$TSS_sd, r$n_species
    ))
}
cat("╚══════════════════════════════════════════════════════════════════════════╝\n")

# Key result
if (nrow(wide_ci) > 0) {
    cat(sprintf(
        "\n★ CI-MLP wins in %d/%d species (%.0f%%), mean ΔAUC = %+.4f\n",
        sum(wide_ci$delta > 0), nrow(wide_ci),
        mean(wide_ci$delta > 0) * 100, mean(wide_ci$delta)
    ))
    if (nrow(wide_ci) > 5) {
        ct <- cor.test(wide_ci$dag_density, wide_ci$delta)
        cat(sprintf(
            "★ DAG density vs CI-MLP advantage: r=%.3f, p=%.4f\n",
            ct$estimate, ct$p.value
        ))
    }
}

cat("\nAll figures saved to figures/case2/\n")
