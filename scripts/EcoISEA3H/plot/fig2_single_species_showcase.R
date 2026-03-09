################################################################################
# Fig 2: Single-Species CAST Workflow Showcase
#
# Demonstrates the complete CAST causal analysis pipeline on one species.
# This figure provides methodological context for the main text by showing
# each step of the CAST workflow in detail.
#
# Panel (a): Consensus DAG — directed acyclic graph learned via bnlearn
#            bootstrap (HC, R=100). Nodes colored by causal role.
# Panel (b): ATE forest plot — Average Treatment Effects with 95% CI
# Panel (c): CAST screening score decomposition (DAG + ATE + RF components)
#            with decision boundary line separating retained/dropped variables
# Panel (d): Model performance comparison — CAST vs ablations and baselines
#
# Species auto-selection: picks the species with most significant ATE
# variables AND DAG edge diversity for an informative showcase.
#
# Data required:
#   output/case2_eco/all_results_v3.csv
#   output/case2_eco/all_screening_v3.csv
#   output/case2_eco/all_ate_results_v3.csv
#   output/case2_eco/all_dag_edges_v3.csv
#   output/case2_eco/all_role_info_v3.csv
#   output/case2_eco/all_dag_info_v3.csv
#
# Run: setwd("E:/CausalSDMs")
#      source("scripts/EcoISEA3H/plot/fig2_single_species_showcase.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

fig_dir <- "figures/case2_eco/plot"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(patchwork)
library(igraph)
library(ggraph)

# ── Theme ────────────────────────────────────────────────────────────────────
theme_pub <- function(base_size = 11) {
    theme_minimal(base_size = base_size, base_family = "sans") +
        theme(
            panel.grid.minor  = element_blank(),
            axis.title        = element_text(face = "bold"),
            plot.title        = element_text(face = "bold", hjust = 0, size = 11),
            plot.subtitle     = element_text(hjust = 0, color = "grey40", size = 8.5),
            legend.background = element_rect(fill = "white", color = NA)
        )
}

pal <- list(
    root = "#2C3E50", med = "#E67E22", term = "#27AE60",
    sig = "#C0392B", nsig = "grey60", cast = "#2980B9"
)

# ── Load ALL data ────────────────────────────────────────────────────────────
results <- read.csv("output/case2_eco/all_results_v3.csv",
    stringsAsFactors = FALSE
)

screening <- read.csv("output/case2_eco/all_screening_v3.csv",
    stringsAsFactors = FALSE
)

ate_all <- read.csv("output/case2_eco/all_ate_results_v3.csv",
    stringsAsFactors = FALSE
) %>%
    mutate(
        coef        = as.numeric(coef),
        se          = as.numeric(se),
        significant = as.logical(significant),
        ci_lower    = coef - 1.96 * se,
        ci_upper    = coef + 1.96 * se
    )

dag_edges <- read.csv("output/case2_eco/all_dag_edges_v3.csv",
    stringsAsFactors = FALSE
)

role_info <- read.csv("output/case2_eco/all_role_info_v3.csv",
    stringsAsFactors = FALSE
)

dag_info <- read.csv("output/case2_eco/all_dag_info_v3.csv",
    stringsAsFactors = FALSE
)

sp_meta <- read.csv(
    "outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv",
    stringsAsFactors = FALSE
) %>%
    mutate(species = gsub(" ", "_", species))

# ── AUTO-SELECT best showcase species ────────────────────────────────────────
# Criteria: most significant ATEs + diverse DAG structure
species_quality <- ate_all %>%
    filter(significant == TRUE) %>%
    group_by(species) %>%
    summarise(n_sig = n(), .groups = "drop") %>%
    left_join(
        dag_info %>% select(species, n_edges, dag_density),
        by = "species"
    ) %>%
    mutate(
        showcase_score = n_sig * 2 + n_edges * 0.1
    ) %>%
    arrange(desc(showcase_score))

target_species <- species_quality$species[1]
target_region <- "China_Res9"

sp_display <- gsub("_", " ", target_species)
sp_family <- sp_meta %>%
    filter(species == target_species) %>%
    pull(family) %>%
    first()
if (is.null(sp_family) || is.na(sp_family)) sp_family <- ""
sp_info <- dag_info %>% filter(species == target_species)

cat(sprintf(
    "[Fig 2] Auto-selected species: %s (%s) | %d sig ATEs | %d DAG edges | density=%.3f\n",
    sp_display, sp_family,
    species_quality$n_sig[1],
    sp_info$n_edges[1],
    sp_info$dag_density[1]
))

# ── Filter species-specific data ─────────────────────────────────────────────
sp_results <- results %>% filter(species == target_species)
sp_screening <- screening %>% filter(species == target_species)
sp_ate <- ate_all %>% filter(species == target_species)
sp_edges <- dag_edges %>% filter(species == target_species)
sp_roles <- role_info %>% filter(species == target_species)

cast_res <- sp_results %>% filter(model == "CAST")
cast_n_vars <- ifelse(nrow(cast_res) > 0, cast_res$n_vars[1], 6)
cast_vars <- sp_screening %>%
    arrange(desc(score_total)) %>%
    slice_head(n = cast_n_vars) %>%
    pull(variable)

# ══════════════════════════════════════════════════════════════════════════════
# Panel (a): Real Consensus DAG network (ggraph)
# ══════════════════════════════════════════════════════════════════════════════
if (nrow(sp_edges) > 0) {
    # Build igraph object from real DAG edges
    g <- graph_from_data_frame(
        sp_edges %>% select(from, to),
        directed = TRUE,
        vertices = data.frame(name = unique(c(sp_edges$from, sp_edges$to)))
    )

    # Add causal role info to nodes
    node_names <- V(g)$name
    node_roles <- sp_roles %>%
        select(variable, group) %>%
        distinct()
    role_map <- setNames(node_roles$group, node_roles$variable)
    V(g)$role <- ifelse(node_names %in% names(role_map),
        role_map[node_names], "Unscreened"
    )

    # Add edge strength
    E(g)$edge_strength <- sp_edges$strength

    # Mark CAST-selected nodes
    V(g)$selected <- ifelse(node_names %in% cast_vars, "Selected", "Dropped")

    # Node size by total degree
    V(g)$deg <- degree(g, mode = "all")

    role_colors <- c(
        "Root" = pal$root, "Mediator" = pal$med,
        "Terminal" = pal$term, "Unscreened" = "grey75"
    )

    set.seed(42)
    pa <- ggraph(g, layout = "fr") +
        geom_edge_arc(
            aes(alpha = edge_strength),
            arrow = arrow(length = unit(2.5, "mm"), type = "closed"),
            end_cap = circle(4, "mm"),
            strength = 0.15,
            color = "grey40"
        ) +
        geom_node_point(
            aes(color = role, size = deg, shape = selected),
            alpha = 0.9
        ) +
        geom_node_text(
            aes(label = name),
            repel = TRUE, size = 2.8, fontface = "bold",
            max.overlaps = 20
        ) +
        scale_color_manual(values = role_colors, name = "Causal\nRole") +
        scale_size_continuous(range = c(3, 9), name = "Degree") +
        scale_shape_manual(
            values = c("Selected" = 16, "Dropped" = 1),
            name = "CAST"
        ) +
        scale_edge_alpha_continuous(range = c(0.3, 0.9), name = "Bootstrap\nstrength") +
        labs(
            title = sprintf("(a) Causal DAG — %s", sp_display),
            subtitle = sprintf(
                "%d edges | density = %.3f | HC bootstrap R=100",
                sp_info$n_edges[1], sp_info$dag_density[1]
            )
        ) +
        theme_void(base_size = 10) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0, size = 11),
            plot.subtitle = element_text(hjust = 0, color = "grey40", size = 8.5),
            legend.text = element_text(size = 7),
            legend.title = element_text(size = 8, face = "bold"),
            legend.position = "right"
        )
} else {
    pa <- ggplot() +
        theme_void() +
        annotate("text",
            x = 0.5, y = 0.5,
            label = "No DAG edge data available",
            size = 4, fontface = "italic", color = "grey50"
        ) +
        labs(title = "(a) Consensus DAG")
}

# ══════════════════════════════════════════════════════════════════════════════
# Panel (b): ATE forest plot with significance & effect direction
# ══════════════════════════════════════════════════════════════════════════════
ate_plot <- sp_ate %>%
    mutate(
        sig_marker = ifelse(significant,
            "Significant (p < 0.05)",
            "Not significant"
        ),
        direction = ifelse(coef > 0, "Positive", "Negative"),
        variable = factor(variable, levels = rev(variable[order(abs(coef))]))
    )

n_sig <- sum(sp_ate$significant)
n_total <- nrow(sp_ate)

pb <- ggplot(ate_plot, aes(
    x = coef, y = variable,
    color = sig_marker
)) +
    geom_vline(
        xintercept = 0, linetype = "dashed", color = "grey60",
        linewidth = 0.6
    ) +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper),
        width = 0.3, linewidth = 0.5, alpha = 0.7
    ) +
    geom_point(size = 3, alpha = 0.9) +
    scale_color_manual(
        values = c(
            "Significant (p < 0.05)" = pal$sig,
            "Not significant" = pal$nsig
        ),
        name = ""
    ) +
    labs(
        title = sprintf("(b) Causal effects (ATE)"),
        subtitle = sprintf(
            "DML 2-fold cross-fitting | %d/%d significant at α=0.05",
            n_sig, n_total
        ),
        x = "Average Treatment Effect", y = ""
    ) +
    theme_pub() +
    theme(
        legend.position = "bottom",
        legend.margin = margin(t = -5),
        panel.grid.major.y = element_line(color = "grey93", linewidth = 0.3)
    )

# ══════════════════════════════════════════════════════════════════════════════
# Panel (c): Screening score decomposition with decision boundary
# ══════════════════════════════════════════════════════════════════════════════
if (nrow(sp_screening) > 0 && "score_dag" %in% names(sp_screening)) {
    # Use stored weights from screening data
    w_dag <- sp_screening$w_dag[1]
    w_ate <- sp_screening$w_ate[1]
    w_imp <- sp_screening$w_imp[1]

    screen_long <- sp_screening %>%
        mutate(
            component_dag = score_dag * w_dag,
            component_ate = score_ate * w_ate,
            component_imp = score_imp * w_imp,
            is_selected = variable %in% cast_vars
        ) %>%
        arrange(desc(score_total)) %>%
        mutate(var_rank = row_number())

    # Decision boundary: threshold between selected and dropped
    n_selected <- length(cast_vars)
    if (n_selected < nrow(screen_long)) {
        threshold <- mean(c(
            screen_long$score_total[n_selected],
            screen_long$score_total[n_selected + 1]
        ))
    } else {
        threshold <- 0
    }

    screen_stacked <- screen_long %>%
        select(
            variable, score_total, var_rank, is_selected,
            component_dag, component_ate, component_imp
        ) %>%
        pivot_longer(
            cols = c(component_dag, component_ate, component_imp),
            names_to = "component", values_to = "score"
        ) %>%
        mutate(
            component = factor(component,
                levels = c("component_dag", "component_ate", "component_imp"),
                labels = c("DAG topology", "ATE causal", "RF importance")
            ),
            variable = factor(variable, levels = screen_long$variable)
        )

    pc <- ggplot(screen_stacked, aes(
        x = variable,
        y = score, fill = component
    )) +
        geom_col(position = "stack", width = 0.7, alpha = 0.85) +
        geom_hline(
            yintercept = threshold, linetype = "longdash",
            color = "#E74C3C", linewidth = 0.7
        ) +
        annotate("text",
            x = nrow(screen_long) * 0.85,
            y = threshold + 0.02,
            label = sprintf("Decision boundary (n=%d retained)", n_selected),
            size = 2.8, fontface = "bold.italic", color = "#E74C3C",
            hjust = 1
        ) +
        scale_fill_manual(values = c(
            "DAG topology" = "#3498DB",
            "ATE causal" = "#E74C3C",
            "RF importance" = "#2ECC71"
        ), name = "Score\ncomponent") +
        labs(
            title = "(c) CAST screening scores",
            subtitle = sprintf(
                "Adaptive weights: w_DAG=%.2f  w_ATE=%.2f  w_RF=%.2f",
                w_dag, w_ate, w_imp
            ),
            x = "", y = "Composite screening score"
        ) +
        theme_pub() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
            legend.position = "right"
        )
} else {
    pc <- ggplot() +
        theme_void() +
        labs(title = "(c) No screening data")
}

# ══════════════════════════════════════════════════════════════════════════════
# Panel (d): Model performance comparison for this species
# ══════════════════════════════════════════════════════════════════════════════
model_order <- c("CAST", "MLP_ATE", "MLP", "RF", "BRT", "Maxent")
model_labels <- c(
    "CAST" = "CAST\n(CI-MLP)",
    "MLP_ATE" = "MLP +\nATE wt",
    "MLP" = "Flat\nMLP",
    "RF" = "Random\nForest",
    "BRT" = "BRT",
    "Maxent" = "Maxent"
)
model_colors <- c(
    "CAST" = "#E74C3C", "MLP_ATE" = "#F39C12", "MLP" = "#3498DB",
    "RF" = "#27AE60", "BRT" = "#8E44AD", "Maxent" = "#16A085"
)

perf_data <- sp_results %>%
    filter(model %in% model_order, !is.na(auc_mean)) %>%
    mutate(
        model = factor(model, levels = model_order),
        is_cast = model == "CAST"
    )

cast_auc <- perf_data %>%
    filter(model == "CAST") %>%
    pull(auc_mean)

pd <- ggplot(perf_data, aes(x = model, y = auc_mean, fill = model)) +
    geom_col(alpha = 0.85, width = 0.65) +
    geom_errorbar(
        aes(ymin = auc_mean - auc_sd, ymax = auc_mean + auc_sd),
        width = 0.2, linewidth = 0.4, alpha = 0.7
    ) +
    geom_text(
        aes(label = sprintf("%.4f", auc_mean)),
        vjust = -0.8, size = 2.8, fontface = "bold"
    ) +
    geom_hline(
        yintercept = cast_auc, linetype = "dotted",
        color = "#E74C3C", linewidth = 0.5, alpha = 0.6
    ) +
    scale_fill_manual(values = model_colors, guide = "none") +
    scale_x_discrete(labels = model_labels) +
    coord_cartesian(ylim = c(
        min(perf_data$auc_mean, na.rm = TRUE) * 0.995,
        max(perf_data$auc_mean + perf_data$auc_sd, na.rm = TRUE) * 1.003
    )) +
    labs(
        title = "(d) Model performance comparison",
        subtitle = sprintf(
            "Test AUC (mean ± SD, %d runs) | CAST vs ablations & baselines",
            max(perf_data$n_success)
        ),
        x = "", y = "Test AUC"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(size = 8))

# ══════════════════════════════════════════════════════════════════════════════
# Combine  2×2  layout
# ══════════════════════════════════════════════════════════════════════════════
fig2 <- (pa | pb) / (pc | pd) +
    plot_layout(heights = c(1.1, 1)) +
    plot_annotation(
        title = sprintf(
            "Fig 2  CAST workflow showcase: %s (%s)",
            sp_display, sp_family
        ),
        subtitle = paste0(
            "Complete single-species pipeline: ",
            "DAG learning → ATE estimation → adaptive screening → model comparison"
        ),
        theme = theme(
            plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
            plot.subtitle = element_text(
                face = "italic", size = 10, hjust = 0.5, color = "grey40"
            )
        )
    )

ggsave(file.path(fig_dir, "fig2_single_species_showcase.png"),
    fig2,
    width = 16, height = 12, dpi = 1200, bg = "white"
)
cat(sprintf("✓ Saved fig2_single_species_showcase.png (species: %s)\n", sp_display))
