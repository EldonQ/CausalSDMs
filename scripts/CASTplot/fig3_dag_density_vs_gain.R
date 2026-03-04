# ==============================================================================
# Fig 3: Mechanistic Insight of Causal AI (Results 3.3)
# Narrative: The predictive gain of CAST stems from pruning the spurious dense
#            interaction space common in deep learning, proportional to the
#            sparseness (informativeness) of the true causal DAG.
# ==============================================================================

rm(list = ls())
setwd("E:/CausalSDMs")

out_dir <- "figures/CASTplot"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(patchwork)

# Theme definition
theme_cast <- function(base_size = 12) {
    theme_minimal(base_size = base_size, base_family = "sans") +
        theme(
            panel.grid.minor = element_blank(),
            axis.title       = element_text(face = "bold"),
            plot.title       = element_text(face = "bold", hjust = 0.5),
            plot.subtitle    = element_text(hjust = 0.5, color = "grey40")
        )
}

# ------------------------------------------------------------------------------
# Panel A: Visualizing "Dense Neural Network" vs "Causal Pruning"
# Mathematical layout for a circular conceptual network illustration
# ------------------------------------------------------------------------------
n_nodes <- 12
theta <- seq(0, 2 * pi, length.out = n_nodes + 1)[1:n_nodes]
nodes <- data.frame(
    id = 1:n_nodes,
    x = cos(theta),
    y = sin(theta),
    label = paste0("V", 1:n_nodes)
)

# FlatNN edges: Fully connected (representing uncontrolled dense interaction space)
edges_flat <- expand.grid(from = 1:n_nodes, to = 1:n_nodes) %>%
    filter(from < to) %>%
    left_join(nodes %>% select(from = id, x1 = x, y1 = y), by = "from") %>%
    left_join(nodes %>% select(to = id, x2 = x, y2 = y), by = "to")

# CI-MLP edges: Highly sparse DAG (simulating the CAST constraints)
edges_dag <- data.frame(
    from = c(1, 1, 1, 2, 3, 3, 4, 6, 6, 8, 8),
    to   = c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12)
) %>%
    left_join(nodes %>% select(from = id, x1 = x, y1 = y), by = "from") %>%
    left_join(nodes %>% select(to = id, x2 = x, y2 = y), by = "to")

p_hairball <- ggplot() +
    geom_segment(data = edges_flat, aes(x = x1, y = y1, xend = x2, yend = y2), color = "gray70", alpha = 0.6) +
    geom_point(data = nodes, aes(x = x, y = y), size = 6, color = "#34495E") +
    theme_void() +
    labs(title = "Traditional Flat NN: Dense Interaction", subtitle = "All factors interact (High Overfitting Risk)") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

p_dag <- ggplot() +
    geom_segment(data = edges_flat, aes(x = x1, y = y1, xend = x2, yend = y2), color = "gray90", alpha = 0.3) +
    geom_segment(data = edges_dag, aes(x = x1, y = y1, xend = x2, yend = y2), color = "#E74C3C", linewidth = 1.2, arrow = arrow(length = unit(0.3, "cm"), type = "closed")) +
    geom_point(data = nodes, aes(x = x, y = y), size = 6, color = "#E74C3C") +
    theme_void() +
    labs(title = "CAST CI-MLP: Causal Pruning", subtitle = "Only retain strict mechanistic graph structures") +
    theme(plot.title = element_text(face = "bold", hjust = 0.5), plot.subtitle = element_text(hjust = 0.5, color = "grey40"))

pA <- (p_hairball | p_dag) + plot_annotation(title = "(A) Mechanistic Concept: Pruning Spurious Interaction Spaces from the Ground Up", theme = theme(plot.title = element_text(size = 14, face = "bold")))

# ------------------------------------------------------------------------------
# Panel B: Pruning Intensity vs Gain (Empirical China EcoISEA3H evidence)
# ------------------------------------------------------------------------------
d_eco <- read.csv("output/case2_eco/all_results_v3.csv", stringsAsFactors = FALSE)
sp_meta_eco <- read.csv("outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv", stringsAsFactors = FALSE) %>%
    mutate(species = gsub(" ", "_", species))

abl_eco <- d_eco %>%
    filter(model %in% c("MLP", "CAST")) %>%
    select(species, model, auc_mean, n_interactions, n_vars) %>%
    pivot_wider(names_from = model, values_from = c(auc_mean, n_interactions, n_vars), values_fn = max) %>%
    filter(!is.na(auc_mean_CAST), !is.na(auc_mean_MLP)) %>%
    mutate(
        delta = auc_mean_CAST - auc_mean_MLP,
        possible_edges = (n_vars_CAST * (n_vars_CAST - 1)) / 2, # 理论全连接线数量
        pruning_intensity = 1 - (n_interactions_CAST / possible_edges) # 因果剪枝强度百分比
    ) %>%
    left_join(sp_meta_eco %>% select(species, family, category), by = "species") %>%
    mutate(family = ifelse(is.na(family), "Unknown", family))

pB <- ggplot(abl_eco %>% filter(!is.na(pruning_intensity)), aes(x = pruning_intensity, y = delta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    geom_point(aes(color = family), size = 3, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid", linewidth = 1) +
    scale_x_continuous(labels = scales::percent_format()) +
    scale_color_brewer(palette = "Set1", name = "Taxonomic Library") +
    labs(
        title = "(B) Empirical Proof: Predictive gains are bounded by structural pruning intensity",
        subtitle = "China EcoISEA3H (n = 27 species). X-axis: Spurious interaction hypothesis space eliminated by causal bounds.",
        x = "Causal Pruning Intensity (% of Dense Path Interactions Prevented)",
        y = expression(bold(Delta * " AUC (CI-MLP vs FlatNN)"))
    ) +
    theme_cast()

# === Combine and Save ===
final_plot <- pA / pB + plot_layout(heights = c(1, 1.2))
ggsave(file.path(out_dir, "fig3_mechanistic_pruning_vs_gain.png"), final_plot, width = 11, height = 10, dpi = 300, bg = "white")
cat("Saved plot to:", file.path(out_dir, "fig3_mechanistic_pruning_vs_gain.png"), "\n")
