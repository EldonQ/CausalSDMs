################################################################################
# 01b_single_species_figures.R
# Publication-quality figures for a SINGLE species CAST pipeline result
#
# Prerequisite: Run 01_cast_pipeline.R first for the target species
#
# Figures:
#   Fig A: Performance comparison bar chart (AUC + TSS)
#   Fig B: CAST variable screening scores (3-component stacked)
#   Fig C: Causal role grouping diagram
#   Fig D: ATE forest plot (effect sizes with CI)
#   Fig E: DAG network edge plot
#   Table: Single-species performance summary
################################################################################

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

target_species <- "swi06" # <--- Must match 01_cast_pipeline.R
target_region <- "SWI" # <--- Must match 01_cast_pipeline.R

library(tidyverse)
library(ggplot2)
library(patchwork)
library(scales)

dir.create(sprintf("figures/case2/%s/single_species", target_region),
    recursive = TRUE, showWarnings = FALSE
)

# ---- Publication color scheme ----
pal <- list(
    ci_mlp  = "#E74C3C", # CI-MLP = red
    flat    = "#3498DB", # FlatNN = blue
    rf      = "#7F8C8D",
    maxent  = "#95A5A6",
    brt     = "#BDC3C7",
    gam     = "#D5DBDB",
    root    = "#E74C3C",
    med     = "#F39C12",
    term    = "#2980B9",
    pos     = "#27AE60",
    neg     = "#E74C3C"
)

theme_pub <- function(base_size = 11) {
    theme_minimal(base_size = base_size) +
        theme(
            text = element_text(family = "sans"),
            plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
            plot.subtitle = element_text(size = base_size, color = "grey40"),
            axis.title = element_text(face = "bold"),
            axis.text = element_text(color = "grey20"),
            legend.position = "bottom",
            legend.title = element_text(face = "bold", size = base_size - 1),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            strip.text = element_text(face = "bold"),
            plot.margin = margin(10, 15, 10, 15)
        )
}

# ==============================================================================
# Load Data
# ==============================================================================
cat(sprintf("Loading results for %s [%s]...\n", target_species, target_region))

# Try region-specific path first, then legacy SWI path
result_base <- sprintf("output/case2/%s/single_species", target_region)
if (!dir.exists(result_base)) {
    result_base <- "output/case2/single_species"
}

results <- read.csv(sprintf("%s/results_%s.csv", result_base, target_species),
    stringsAsFactors = FALSE
)
screening <- read.csv(sprintf("%s/screening_%s.csv", result_base, target_species),
    stringsAsFactors = FALSE
)
roles_file <- sprintf("%s/roles_%s.csv", result_base, target_species)
has_roles <- file.exists(roles_file)
if (has_roles) {
    roles <- read.csv(roles_file, stringsAsFactors = FALSE)
}

# Load ATE from saved CSV (proper DML SE from cross-fitting)
ate_file <- sprintf("%s/ate_%s.csv", result_base, target_species)
if (file.exists(ate_file)) {
    ate_data <- read.csv(ate_file, stringsAsFactors = FALSE) %>%
        mutate(ci_lower = coef - 1.96 * se, ci_upper = coef + 1.96 * se)
} else {
    ate_data <- screening %>%
        filter(!is.na(coef)) %>%
        select(variable, coef, p_value, significant) %>%
        mutate(
            se = abs(coef) / pmax(qnorm(1 - p_value / 2), 0.1),
            ci_lower = coef - 1.96 * se,
            ci_upper = coef + 1.96 * se
        )
}

# Standardize model labels for the three-group design
results <- results %>%
    mutate(
        # Extract base algorithm name
        base_model = gsub("_(cast|full)$", "", model),
        # Human-readable labels
        model_label = case_when(
            model == "CI_MLP" ~ "CI-MLP",
            model == "FlatNN_cast" ~ "FlatNN",
            model == "FlatNN_full" ~ "FlatNN",
            base_model == "RF" ~ "RF",
            base_model == "Maxent" ~ "Maxent",
            base_model == "BRT" ~ "BRT",
            base_model == "GAM" ~ "GAM",
            TRUE ~ model
        ),
        group = case_when(
            model == "CI_MLP" ~ "C: CI-MLP",
            var_set == "full" | grepl("_full$", model) ~ "A: Full variables",
            TRUE ~ "B: CAST screened"
        )
    )

model_order <- c("Maxent", "GAM", "RF", "BRT", "FlatNN", "CI-MLP")
results$model_label <- factor(results$model_label, levels = model_order)
group_order <- c("A: Full variables", "B: CAST screened", "C: CI-MLP")
results$group <- factor(results$group, levels = group_order)

cat(sprintf("  %d model-group combinations loaded\n", nrow(results)))

# ==============================================================================
# Figure A: Performance comparison (three-group design)
# ==============================================================================
cat("\nFig A: Performance comparison (3 groups)...\n")

color_map <- c(
    "A: Full variables" = "#95A5A6",
    "B: CAST screened" = "#3498DB",
    "C: CI-MLP" = "#E74C3C"
)

# Prepare long format
perf_long <- results %>%
    select(model_label, group, auc_mean, auc_sd, tss_mean, tss_sd) %>%
    pivot_longer(
        cols = c(auc_mean, tss_mean),
        names_to = "metric_raw", values_to = "value"
    ) %>%
    mutate(
        metric = ifelse(grepl("auc", metric_raw), "AUC", "TSS"),
        sd_val = ifelse(grepl("auc", metric_raw),
            results$auc_sd[match(paste(model_label, group), paste(results$model_label, results$group))],
            results$tss_sd[match(paste(model_label, group), paste(results$model_label, results$group))]
        )
    )

figA <- ggplot(perf_long, aes(x = model_label, y = value, fill = group)) +
    geom_col(
        position = position_dodge(width = 0.75), width = 0.65,
        alpha = 0.9, color = "white", linewidth = 0.3
    ) +
    geom_errorbar(aes(ymin = value - sd_val, ymax = value + sd_val),
        position = position_dodge(width = 0.75),
        width = 0.2, linewidth = 0.5
    ) +
    geom_text(aes(label = sprintf("%.3f", value)),
        position = position_dodge(width = 0.75),
        vjust = -0.8, size = 2.5, fontface = "bold"
    ) +
    facet_wrap(~metric, scales = "free_y") +
    scale_fill_manual(values = color_map, name = "Experimental Group") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
        title = sprintf("Model performance comparison: %s [%s]", target_species, target_region),
        subtitle = sprintf(
            "A: Full vars (%d) | B: CAST screened (%d) | C: CI-MLP (CAST + DAG interactions)",
            results$n_vars[results$group == "A: Full variables"][1],
            results$n_vars[results$group != "A: Full variables"][1]
        ),
        x = "", y = "Score"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))

fig_dir <- sprintf("figures/case2/%s/single_species", target_region)
ggsave(sprintf("%s/%s_figA_performance.png", fig_dir, target_species),
    figA,
    width = 14, height = 7, dpi = 300, bg = "white"
)
ggsave(sprintf("%s/%s_figA_performance.pdf", fig_dir, target_species),
    figA,
    width = 14, height = 7, bg = "white"
)
cat("  ✓ figA_performance\n")

# ==============================================================================
# Figure B: CAST Variable Screening Scores (3-component stacked bar)
# ==============================================================================
cat("\nFig B: Variable screening scores...\n")

# Compute weighted contributions for stacking visualization
# Weights: DAG=0.2, ATE=0.4, Importance=0.4
screen_long <- screening %>%
    select(variable, score_dag, score_ate, score_imp, score_total) %>%
    mutate(
        w_dag = 0.2 * score_dag,
        w_ate = 0.4 * score_ate,
        w_imp = 0.4 * score_imp
    ) %>%
    pivot_longer(
        cols = c(w_dag, w_ate, w_imp),
        names_to = "component", values_to = "weighted_score"
    ) %>%
    mutate(component = case_when(
        component == "w_dag" ~ "DAG Out-degree (×0.2)",
        component == "w_ate" ~ "ATE Effect (×0.4)",
        component == "w_imp" ~ "RF Importance (×0.4)"
    ))

# Determine threshold (elbow with median fallback, same logic as 01)
sorted_scores <- sort(screening$score_total, decreasing = TRUE)
if (length(sorted_scores) > 4) {
    score_drops <- -diff(sorted_scores)
    elbow_idx <- which.max(score_drops)
    elbow_threshold <- sorted_scores[elbow_idx]
    n_by_elbow <- sum(screening$score_total >= elbow_threshold)
    if (n_by_elbow < 3 || n_by_elbow > 0.8 * nrow(screening)) {
        sel_threshold <- median(sorted_scores)
        thr_label <- "median"
    } else {
        sel_threshold <- elbow_threshold
        thr_label <- "elbow"
    }
} else {
    sel_threshold <- median(sorted_scores)
    thr_label <- "median"
}
cast_vars <- screening$variable[screening$score_total >= sel_threshold]
if (length(cast_vars) < 3) cast_vars <- screening$variable[1:min(3, nrow(screening))]

figB <- ggplot(screen_long, aes(
    x = reorder(variable, -weighted_score),
    y = weighted_score, fill = component
)) +
    geom_col(position = "stack", width = 0.7, alpha = 0.85) +
    geom_hline(
        yintercept = sel_threshold, linetype = "dashed",
        color = pal$ci_mlp, linewidth = 0.8
    ) +
    annotate("text",
        x = nrow(screening) - 0.5, y = sel_threshold + 0.02,
        label = sprintf("Threshold = %.3f (%s)", sel_threshold, thr_label),
        color = pal$ci_mlp, fontface = "italic", hjust = 1, size = 3
    ) +
    scale_fill_manual(values = c(
        "DAG Out-degree (×0.2)" = "#2C3E50",
        "ATE Effect (×0.4)" = "#E67E22",
        "RF Importance (×0.4)" = "#27AE60"
    ), name = "Weighted Component") +
    labs(
        title = sprintf("CAST variable screening: %s", target_species),
        subtitle = sprintf(
            "Weighted score (DAG:0.2, ATE:0.4, Imp:0.4) → %d/%d selected (%s threshold)",
            length(cast_vars), nrow(screening), thr_label
        ),
        x = "Environmental variable", y = "Weighted screening score"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(sprintf("%s/%s_figB_screening.png", fig_dir, target_species),
    figB,
    width = 10, height = 6, dpi = 300, bg = "white"
)
ggsave(sprintf("%s/%s_figB_screening.pdf", fig_dir, target_species),
    figB,
    width = 10, height = 6, bg = "white"
)
cat("  ✓ figB_screening\n")

# ==============================================================================
# Figure C: Causal Role Grouping Diagram
# ==============================================================================
cat("\nFig C: Causal role grouping...\n")

group_colors <- c("Root" = pal$root, "Mediator" = pal$med, "Terminal" = pal$term)

roles_plot <- roles %>%
    mutate(group = factor(group, levels = c("Root", "Mediator", "Terminal")))

figC <- ggplot(roles_plot, aes(
    x = reorder(variable, -role_score),
    y = role_score, fill = group
)) +
    geom_col(width = 0.7, alpha = 0.85, color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("out=%d\nin=%d", out, inp)),
        vjust = -0.3, size = 2.8, lineheight = 0.8
    ) +
    scale_fill_manual(values = group_colors, name = "Causal Role") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(
        title = sprintf("Causal role assignment: %s", target_species),
        subtitle = "Role score = out_degree / (in_degree + 1); Root → Mediator → Terminal",
        x = "Variable", y = "Role score (higher = more upstream)"
    ) +
    theme_pub()

# Add flow diagram below
flow_text <- roles_plot %>%
    group_by(group) %>%
    summarise(vars = paste(variable, collapse = ", "), .groups = "drop") %>%
    arrange(match(group, c("Root", "Mediator", "Terminal")))

flow_label <- paste(
    sprintf("%s [%s]", flow_text$group, flow_text$vars),
    collapse = "  →  "
)

figC <- figC +
    labs(caption = sprintf("Information flow: %s", flow_label)) +
    theme(plot.caption = element_text(face = "bold", size = 10, color = "#2C3E50", hjust = 0.5))

ggsave(sprintf("%s/%s_figC_causal_roles.png", fig_dir, target_species),
    figC,
    width = 10, height = 6, dpi = 300, bg = "white"
)
ggsave(sprintf("%s/%s_figC_causal_roles.pdf", fig_dir, target_species),
    figC,
    width = 10, height = 6, bg = "white"
)
cat("  ✓ figC_causal_roles\n")

# ==============================================================================
# Figure D: ATE Forest Plot (effect sizes with 95% CI)
# ==============================================================================
cat("\nFig D: ATE forest plot...\n")

ate_plot <- ate_data %>%
    mutate(
        sig_marker = ifelse(significant, "Significant (p<0.05)", "Not significant"),
        var_label = sprintf("%s %s", variable, ifelse(significant, "*", ""))
    )

figD <- ggplot(ate_plot, aes(
    x = coef, y = reorder(variable, abs(coef)),
    color = sig_marker
)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.25, linewidth = 0.6) +
    scale_color_manual(values = c(
        "Significant (p<0.05)" = pal$ci_mlp,
        "Not significant" = "grey60"
    ), name = "") +
    labs(
        title = sprintf("Average Treatment Effect (ATE) estimates: %s", target_species),
        subtitle = "DML with 2-fold cross-fitting (Chernozhukov et al., 2018)",
        x = "ATE coefficient (effect on presence probability)", y = ""
    ) +
    theme_pub() +
    theme(panel.grid.major.y = element_line(color = "grey92"))

ggsave(sprintf("%s/%s_figD_ate_forest.png", fig_dir, target_species),
    figD,
    width = 10, height = 6, dpi = 300, bg = "white"
)
ggsave(sprintf("%s/%s_figD_ate_forest.pdf", fig_dir, target_species),
    figD,
    width = 10, height = 6, bg = "white"
)
cat("  ✓ figD_ate_forest\n")

# ==============================================================================
# Figure E: Combined summary panel (screening + roles + performance)
# ==============================================================================
cat("\nFig E: Combined summary panel...\n")

figE <- (figB | figC) / (figD | figA) +
    plot_annotation(
        title = sprintf("CAST v3 Pipeline Summary: %s [%s]", target_species, target_region),
        subtitle = "Causal screening → Role grouping → ATE estimation → Model comparison",
        theme = theme(
            plot.title = element_text(face = "bold", size = 16),
            plot.subtitle = element_text(size = 12, color = "grey40")
        )
    ) +
    plot_layout(heights = c(1, 1))

ggsave(sprintf("%s/%s_figE_summary_panel.png", fig_dir, target_species),
    figE,
    width = 18, height = 14, dpi = 300, bg = "white"
)
ggsave(sprintf("%s/%s_figE_summary_panel.pdf", fig_dir, target_species),
    figE,
    width = 18, height = 14, bg = "white"
)
cat("  ✓ figE_summary_panel\n")

# ==============================================================================
# TABLE: Single-species performance summary
# ==============================================================================
cat("\nGenerating performance table...\n")

# Find best model
best_model <- results$model_label[which.max(results$auc_mean)]

pub_table <- results %>%
    arrange(desc(auc_mean)) %>%
    mutate(
        AUC = sprintf("%.4f ± %.4f", auc_mean, auc_sd),
        TSS = sprintf("%.4f ± %.4f", tss_mean, tss_sd),
        Rank = row_number()
    ) %>%
    select(Rank, model, model_label, group, AUC, TSS)

write.csv(pub_table,
    sprintf("%s/table_%s.csv", result_base, target_species),
    row.names = FALSE
)

# Print formatted table
cat(sprintf("\n  ╔════════════════════════════════════════════════════════════════════════════════╗\n"))
cat(sprintf("  ║  Performance Summary: %-10s [%s]                                          ║\n", target_species, target_region))
cat("  ╠════════════════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("  ║  %-4s  %-16s  %-24s  %-14s  %-14s  ║\n", "Rank", "Model", "Group", "AUC", "TSS"))
cat("  ╠════════════════════════════════════════════════════════════════════════════════╣\n")
for (i in 1:nrow(pub_table)) {
    r <- pub_table[i, ]
    mk <- if (as.character(r$model) == "CI_MLP") " ★" else "  "
    cat(sprintf(
        "  ║%s%-4d  %-16s  %-24s  %-14s  %-14s  ║\n",
        mk, r$Rank, as.character(r$model), as.character(r$group), r$AUC, r$TSS
    ))
}
cat("  ╚════════════════════════════════════════════════════════════════════════════════╝\n")

# Pipeline summary table
n_cast_vars <- results$n_vars[results$model == "CI_MLP"][1]
if (is.na(n_cast_vars)) n_cast_vars <- results$n_vars[results$group == "B: CAST screened"][1]
n_total_vars <- nrow(screening)
cat("\n  ╔══════════════════════════════════════════════════════════════════════╗\n")
cat("  ║  CAST v3 Pipeline Summary                                          ║\n")
cat("  ╠══════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf(
    "  ║  VIF:        → %d variables                                        ║\n",
    n_total_vars
))
cat(sprintf(
    "  ║  ATE:        %d/%d significant (p<0.05)                             ║\n",
    sum(ate_data$significant), nrow(ate_data)
))
cat(sprintf(
    "  ║  Screening:  %d → %d CAST variables (%.0f%% reduction)               ║\n",
    n_total_vars, n_cast_vars,
    (1 - n_cast_vars / n_total_vars) * 100
))
if (has_roles) {
    flow_parts <- roles %>%
        group_by(group) %>%
        summarise(vars = paste(variable, collapse = ","), .groups = "drop")
    cat(sprintf(
        "  ║  Grouping:   %s                  ║\n",
        paste(sprintf("%s[%s]", flow_parts$group, flow_parts$vars), collapse = " → ")
    ))
}
cat("  ╚══════════════════════════════════════════════════════════════════════╝\n")

# ==============================================================================
# Summary
# ==============================================================================
cat(sprintf("\n======================================================================\n"))
cat(sprintf("  Single-Species Figures Complete: %s [%s]\n", target_species, target_region))
cat("======================================================================\n")
cat(sprintf("  Figures saved to %s/:\n", fig_dir))
cat(sprintf("    %s_figA_performance      (.png + .pdf)\n", target_species))
cat(sprintf("    %s_figB_screening        (.png + .pdf)\n", target_species))
cat(sprintf("    %s_figC_causal_roles     (.png + .pdf)\n", target_species))
cat(sprintf("    %s_figD_ate_forest       (.png + .pdf)\n", target_species))
cat(sprintf("    %s_figE_summary_panel    (.png + .pdf)\n", target_species))
cat(sprintf("  Tables saved to %s/:\n", result_base))
cat(sprintf("    table_%s.csv\n", target_species))
cat("======================================================================\n")
