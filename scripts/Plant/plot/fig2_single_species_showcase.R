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
#
# Species selection: You can manually set USER_TARGET_SPECIES below to showcase
# a specific species, or leave it NULL to auto-select the most informative one.
#
# Data required (Plant case, output from 03_run_Plant_multi_species.R):
#   output/case4_plant/all_results_plant.csv
#   output/case4_plant/all_screening_plant.csv
#   output/case4_plant/all_ate_results_plant.csv
#   output/case4_plant/all_dag_edges_plant.csv
#   output/case4_plant/all_role_info_plant.csv
#   output/case4_plant/all_dag_info_plant.csv
#   (optional) outputs/Plant/Res9/CAST_ready/CAST_Species_Summary.csv for family
#
#
# Run: setwd("E:/CausalSDMs")
#      source("scripts/Plant/plot/fig2_single_species_showcase.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

# ── User Settings ────────────────────────────────────────────────────────────
SHOW_MOCK_BOUNDARY <- FALSE # If TRUE, draws a red dashed decision boundary
MOCK_RETAIN_N <- 6 # Number of variables to retain if using mock boundary

fig_dir <- "figures/case4_plant/plot"
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
results <- read.csv("output/case4_plant/all_results_plant.csv",
    stringsAsFactors = FALSE
)

screening <- read.csv("output/case4_plant/all_screening_plant.csv",
    stringsAsFactors = FALSE
)

ate_all <- read.csv("output/case4_plant/all_ate_results_plant.csv",
    stringsAsFactors = FALSE
) %>%
    mutate(
        coef        = as.numeric(coef),
        se          = as.numeric(se),
        significant = as.logical(significant),
        ci_lower    = coef - 1.96 * se,
        ci_upper    = coef + 1.96 * se
    )

dag_edges <- read.csv("output/case4_plant/all_dag_edges_plant.csv",
    stringsAsFactors = FALSE
)

role_info <- read.csv("output/case4_plant/all_role_info_plant.csv",
    stringsAsFactors = FALSE
)

dag_info <- read.csv("output/case4_plant/all_dag_info_plant.csv",
    stringsAsFactors = FALSE
)

# 物种元数据（若有 CAST_Species_Summary 则取 family，否则用空表）
sp_meta_path <- "outputs/Plant/Res9/CAST_ready/CAST_Species_Summary.csv"
if (file.exists(sp_meta_path)) {
    sp_meta <- read.csv(sp_meta_path, stringsAsFactors = FALSE) %>%
        mutate(species = gsub(" ", "_", species))
} else {
    sp_meta <- data.frame(species = unique(results$species), family = "Plant", stringsAsFactors = FALSE)
}

# ── AUTO-SELECT best showcase species ────────────────────────────────────────
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

# ── 展示物种列表（Plant 案例：从质量分最高的物种中取前 3 个，不足则用全部）────
target_species_list <- species_quality %>%
    slice_head(n = 3) %>%
    pull(species)
if (length(target_species_list) == 0) {
    target_species_list <- unique(dag_info$species)[seq_len(min(3, nrow(dag_info)))]
}

for (target_species in target_species_list) {
    cat(sprintf("\n[Fig 2] Processing target species: %s\n", target_species))

    sp_display <- gsub("_", " ", target_species)
    sp_family <- sp_meta %>%
        filter(species == target_species) %>%
        pull(family) %>%
        first()
    if (is.null(sp_family) || is.na(sp_family)) sp_family <- ""

    sp_info <- dag_info %>% filter(species == target_species)
    if (nrow(sp_info) == 0) {
        cat(paste0("Warning: No DAG info found for ", target_species, "\n"))
        next
    }

    # ── Filter species-specific data ─────────────────────────────────────────────
    sp_results <- results %>% filter(species == target_species)
    sp_screening <- screening %>% filter(species == target_species)
    sp_ate <- ate_all %>% filter(species == target_species)
    sp_edges <- dag_edges %>% filter(species == target_species)
    sp_roles <- role_info %>% filter(species == target_species)

    cast_res <- sp_results %>% filter(model == "CAST")

    if (SHOW_MOCK_BOUNDARY) {
        cast_n_vars <- MOCK_RETAIN_N
    } else {
        cast_n_vars <- ifelse(nrow(cast_res) > 0, cast_res$n_vars[1], 6)
    }

    cast_vars <- sp_screening %>%
        arrange(desc(score_total)) %>%
        slice_head(n = cast_n_vars) %>%
        pull(variable)

    # ── Clean Variable Names for Presentation (含 Plant 案例变量) ─────────────────
    var_dict <- c(
        "bio02" = "Diurnal Temp Range", "bio_2" = "Mean Diurnal Range",
        "bio15" = "Precip Seasonality", "bio_15" = "Precip Seasonality",
        "bio19" = "Precip Coldest Qtr", "bio_19" = "Precip Coldest Qtr",
        "bio03" = "Isothermality", "bio_3" = "Isothermality",
        "bio18" = "Precip Warmest Qtr", "bio_18" = "Precip Warmest Qtr",
        "maxtempcoldest" = "Max Temp Coldest Mo",
        "aridityindexthornthwaite" = "Aridity Index",
        "etccdi_cwd" = "Consecutive Wet Days",
        "elevation" = "Elevation", "Elevation" = "Elevation",
        "tri" = "Terrain Ruggedness", "topowet" = "Topo Wetness Index",
        "nontree" = "Non-tree Veg", "landcover_igbp" = "Landcover Type",
        "Slope" = "Slope", "Aspect" = "Aspect",
        "ORCDRC" = "Soil Organic C", "PHIHOX" = "Soil pH",
        "CECSOL" = "Soil CEC", "CLYPPT" = "Clay Content",
        "SLTPPT" = "Silt Content", "BDTICM" = "Soil Bulk Density",
        "Lights2009" = "Night Lights", "Built2009" = "Built-up",
        "Croplands2005" = "Croplands", "Pasture2009" = "Pasture"
    )
    clean_var_name <- function(v) {
        ifelse(v %in% names(var_dict), var_dict[v], v)
    }

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
                aes(label = clean_var_name(name)),
                repel = TRUE, size = 3.2, fontface = "bold",
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
            clean_name = clean_var_name(variable),
            clean_name = factor(clean_name, levels = rev(clean_var_name(variable[order(abs(coef))])))
        )

    n_sig <- sum(sp_ate$significant)
    n_total <- nrow(sp_ate)

    pb <- ggplot(ate_plot, aes(
        x = coef, y = clean_name,
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
            x = "", y = ""
        ) +
        theme_pub() +
        theme(
            legend.position = "bottom",
            legend.margin = margin(t = -5),
            legend.box.margin = margin(t = -15),
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
                clean_name = clean_var_name(variable),
                clean_name = factor(clean_name, levels = clean_var_name(screen_long$variable))
            )

        pc <- ggplot(screen_stacked, aes(
            x = clean_name,
            y = score, fill = component, alpha = is_selected
        )) +
            geom_col(position = "stack", width = 0.7) +
            scale_fill_manual(values = c(
                "DAG topology" = "#3498DB",
                "ATE causal" = "#E74C3C",
                "RF importance" = "#2ECC71"
            ), name = "")

        if (SHOW_MOCK_BOUNDARY) {
            pc <- pc +
                scale_alpha_manual(values = c("TRUE" = 0.95, "FALSE" = 0.35), guide = "none") +
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
                )
        } else {
            pc <- pc + scale_alpha_manual(values = c("TRUE" = 0.85, "FALSE" = 0.85), guide = "none")
        }

        pc <- pc + labs(
            title = "(c) CAST screening scores",
            subtitle = sprintf(
                "Adaptive weights: w_DAG=%.2f  w_ATE=%.2f  w_RF=%.2f",
                w_dag, w_ate, w_imp
            ),
            x = "", y = "Composite screening score"
        ) +
            coord_cartesian(clip = "off") +
            theme_pub() +
            theme(
                axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
                legend.position = "bottom",
                legend.margin = margin(t = -5),
                legend.box.margin = margin(t = -15),
                plot.margin = margin(t = 5, r = 45, b = 5, l = 5, unit = "mm")
            )
    } else {
        pc <- ggplot() +
            theme_void() +
            labs(title = "(c) No screening data")
    }

    # ══════════════════════════════════════════════════════════════════════════════
    # Combine 1×3 layout (Horizontal row)
    # ══════════════════════════════════════════════════════════════════════════════
    # Note: Panel (d) was removed to focus purely on the causal/ecological narrative
    # without distracting with redundant benchmark performance metrics.
    fig2 <- (pa | pb | pc) +
        plot_layout(widths = c(1.2, 1, 1)) +
        plot_annotation(
            title = sprintf("Fig 2  CAST causal mechanism showcase: %s (%s)", sp_display, sp_family),
            subtitle = "Complete single-species pipeline: DAG structure learning → ATE estimation → Adaptive causal screening",
            theme = theme(
                plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
                plot.subtitle = element_text(face = "italic", size = 12, hjust = 0.5, color = "grey40")
            )
        )

    out_file <- file.path(fig_dir, sprintf("fig2_single_species_showcase_%s.png", target_species))
    ggsave(out_file, fig2, width = 22, height = 7, dpi = 1200, bg = "white")
    cat(sprintf("✓ Saved %s\n", out_file))
}
