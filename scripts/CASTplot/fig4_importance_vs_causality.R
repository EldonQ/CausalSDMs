# ==============================================================================
# Fig 4: 因果vs相关重要性的解耦 (Results 3.4)
# Narrative: 纯预测相关性 (RF Permutation Importance) 与 因果推断 (Absolute ATE) 的彻底解耦
# Format: Rank Slope Chart per species / or Spearman Correlation histogram across species.
# ==============================================================================

rm(list = ls())
setwd("E:/CausalSDMs")

out_dir <- "figures/CASTplot"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(patchwork)

theme_cast <- function(base_size = 12) {
    theme_minimal(base_size = base_size, base_family = "sans") +
        theme(
            panel.grid.minor = element_blank(),
            axis.title       = element_text(face = "bold"),
            plot.title       = element_text(face = "bold", hjust = 0.5),
            plot.subtitle    = element_text(hjust = 0.5, color = "grey40")
        )
}

# ---- Load Real EcoISEA3H Data ----
screen <- read.csv("output/case2_eco/all_screening_v3.csv", stringsAsFactors = FALSE)
ate <- read.csv("output/case2_eco/all_ate_results_v3.csv", stringsAsFactors = FALSE)

showcase_region <- "China_Res9"

rank_data <- screen %>%
    left_join(ate %>% select(region, species, variable, ate_coef = coef), by = c("region", "species", "variable")) %>%
    mutate(ate_coef = as.numeric(ate_coef)) %>%
    group_by(species, region) %>%
    mutate(
        rank_rf  = rank(-importance, ties.method = "average"),
        rank_ate = rank(-abs(ate_coef), ties.method = "average", na.last = TRUE)
    ) %>%
    ungroup()

slope_dat <- rank_data %>%
    group_by(variable) %>%
    summarise(
        mean_rank_rf = mean(rank_rf, na.rm = TRUE),
        mean_rank_ate = mean(rank_ate, na.rm = TRUE),
        n_sp = n(), .groups = "drop"
    ) %>%
    filter(n_sp >= 5) %>%
    mutate(y_left = rank(mean_rank_rf), y_right = rank(mean_rank_ate))

library(ggrepel)

# Set custom colors for lines (more than 12)
my_colors <- colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))(nrow(slope_dat))

pa <- ggplot(slope_dat) +
    geom_segment(aes(x = 1, xend = 2, y = y_left, yend = y_right, color = variable), linewidth = 1.1, alpha = 0.8) +
    geom_point(aes(x = 1, y = y_left), size = 2.5, color = "#27AE60") +
    geom_point(aes(x = 2, y = y_right), size = 2.5, color = "#C0392B") +
    geom_text_repel(aes(x = 1, y = y_left, label = variable), direction = "y", nudge_x = -0.1, hjust = 1, size = 3.5, fontface = "bold", segment.color = NA) +
    geom_text_repel(aes(x = 2, y = y_right, label = variable), direction = "y", nudge_x = 0.1, hjust = 0, size = 3.5, fontface = "bold", segment.color = NA) +
    scale_x_continuous(limits = c(0.4, 2.6), breaks = c(1, 2), labels = c("Predictive Rank\n(RF Permutation)", "Causal Rank\n(|ATE| from DML)")) +
    scale_color_manual(values = my_colors, guide = "none") +
    scale_y_reverse() +
    labs(title = "(A) Rank Discrepancy", subtitle = "Line crossings = RF importance vs Causal impact", x = "", y = "Mean rank across species") +
    theme_cast()

spearman_per_sp <- rank_data %>%
    filter(!is.na(rank_rf), !is.na(rank_ate)) %>%
    group_by(species) %>%
    summarise(rho = cor(rank_rf, rank_ate, method = "spearman", use = "pairwise"), n_vars = n(), .groups = "drop") %>%
    filter(!is.na(rho))

n_species_actual <- nrow(spearman_per_sp)

pb <- ggplot(spearman_per_sp, aes(x = rho)) +
    geom_histogram(binwidth = 0.1, fill = "#2980B9", alpha = 0.7, color = "white") +
    geom_vline(xintercept = mean(spearman_per_sp$rho, na.rm = TRUE), linetype = "dashed", color = "#C0392B", linewidth = 0.8) +
    annotate("text", x = mean(spearman_per_sp$rho, na.rm = TRUE) + 0.05, y = Inf, vjust = 2, label = sprintf("Mean ρ = %.2f", mean(spearman_per_sp$rho, na.rm = TRUE)), fontface = "bold", size = 4, color = "#C0392B") +
    scale_x_continuous(limits = c(-1, 1)) +
    labs(title = sprintf("(B) Global Rank Correlations (n = %d species)", n_species_actual), subtitle = "The decoupling of prediction and causality is systematic", x = "Spearman Rank Correlation (ρ)", y = "Number of species") +
    theme_cast()

final_plot <- (pa | pb) + plot_annotation(title = "Fig 4. Decoupling Predictive Importance from Ecological Mechanism", theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
ggsave(file.path(out_dir, "fig4_importance_vs_causality.png"), final_plot, width = 12, height = 6, dpi = 300, bg = "white")
cat("Saved plot to:", file.path(out_dir, "fig4_importance_vs_causality.png"), "\n")
