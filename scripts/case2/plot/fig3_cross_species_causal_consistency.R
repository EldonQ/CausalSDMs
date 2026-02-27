################################################################################
# Fig 3: Cross-species causal consistency (4-panel)
# (a) DAG edge frequency heatmap  (b) ATE sign consistency  (c) Variable selection frequency  (d) Adaptive weights ternary
# Requires: all_dag_edges_v3.csv, all_ate_results_v3.csv, all_screening_v3.csv
# Run from project root: setwd("E:/CausalSDMs"); source("scripts/case2/plot/fig3_cross_species_causal_consistency.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

fig_dir <- "figures/case2/plot"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_dag   <- "output/case2/all_dag_edges_v3.csv"
path_ate   <- "output/case2/all_ate_results_v3.csv"
path_screen <- "output/case2/all_screening_v3.csv"

if (!file.exists(path_dag))   stop("Missing ", path_dag)
if (!file.exists(path_ate))   stop("Missing ", path_ate)
if (!file.exists(path_screen)) stop("Missing ", path_screen)

library(tidyverse)
library(ggplot2)
library(patchwork)

theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40")
    )
}

dag   <- read.csv(path_dag,   stringsAsFactors = FALSE)
ate   <- read.csv(path_ate,   stringsAsFactors = FALSE)
screen <- read.csv(path_screen, stringsAsFactors = FALSE)

n_species <- length(unique(paste(dag$region, dag$species)))

# ---- (a) Edge frequency heatmap (variable x variable, facet by region) ----
region_n <- dag %>% distinct(region, species) %>% count(region, name = "n_sp")
edge_counts <- dag %>%
  group_by(region, from, to) %>%
  summarise(n_species = n(), .groups = "drop") %>%
  left_join(region_n, by = "region") %>%
  mutate(frac = n_species / n_sp)
edge_mat <- edge_counts
pa <- ggplot(edge_mat, aes(x = to, y = from, fill = frac)) +
  geom_tile(color = "white", linewidth = 0.3) +
  facet_wrap(~region, scales = "free", ncol = 3) +
  scale_fill_viridis_c(option = "plasma", name = "Fraction of species", limits = c(0, 1)) +
  labs(title = "(a) DAG edge frequency by region", subtitle = "Cell = fraction of species with that edge",
    x = "", y = "") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")

# ---- (b) ATE sign consistency (diverging bar: % negative vs % positive per variable) ----
ate_summ <- ate %>%
  filter(!is.na(coef)) %>%
  mutate(sign = case_when(coef > 0 ~ "positive", coef < 0 ~ "negative", TRUE ~ "zero")) %>%
  count(variable, region, sign, .drop = FALSE) %>%
  group_by(variable, region) %>%
  mutate(total = sum(n), pct = n / total) %>%
  ungroup()
ate_wide <- ate_summ %>%
  filter(sign %in% c("positive", "negative")) %>%
  select(variable, region, sign, pct) %>%
  pivot_wider(names_from = sign, values_from = pct, values_fill = 0)
pb <- ggplot(ate_wide, aes(x = reorder(variable, positive - negative))) +
  geom_col(aes(y = positive), fill = "#C0392B", alpha = 0.85, width = 0.6) +
  geom_col(aes(y = -negative), fill = "#2980B9", alpha = 0.85, width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  facet_wrap(~region, ncol = 3) +
  coord_flip() +
  scale_y_continuous(labels = function(x) paste0(abs(x) * 100, "%")) +
  labs(title = "(b) ATE sign consistency across species", subtitle = "Red = % species with positive ATE; Blue = negative",
    x = "", y = "Fraction of species") +
  theme_pub()

# ---- (c) Variable selection frequency (variable x region tile) ----
# Infer selected: top 6 or top 50% by score_total per (region, species)
sel_freq <- screen %>%
  group_by(region, species) %>%
  mutate(rank = row_number(desc(score_total)), n_v = n(),
    selected = rank <= pmin(6, max(1, ceiling(n_v * 0.5)))) %>%
  ungroup() %>%
  group_by(region, variable) %>%
  summarise(frac = mean(selected, na.rm = TRUE), .groups = "drop")
pc <- ggplot(sel_freq, aes(x = region, y = reorder(variable, frac), fill = frac)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_viridis_c(option = "viridis", name = "Selection frequency", limits = c(0, 1)) +
  labs(title = "(c) Variable selection frequency", subtitle = "Fraction of species in which variable was CAST-selected",
    x = "Region", y = "") +
  theme_pub() + theme(legend.position = "right")

# ---- (d) Adaptive weights ternary (w_dag, w_ate, w_imp) ----
if (all(c("w_dag", "w_ate", "w_imp") %in% names(screen))) {
  w_agg <- screen %>%
    distinct(region, species, w_dag, w_ate, w_imp) %>%
    group_by(region) %>%
    summarise(w_dag = mean(w_dag, na.rm = TRUE), w_ate = mean(w_ate, na.rm = TRUE), w_imp = mean(w_imp, na.rm = TRUE), .groups = "drop")
  # Ternary: x = w_imp, y = w_ate, z = w_dag (or use ggtern)
  # Simple 2D: x = w_ate, y = w_imp, color = region (w_dag = 1 - x - y)
  pd <- ggplot(w_agg, aes(x = w_ate, y = w_imp, color = region)) +
    geom_point(size = 4, alpha = 0.9) +
    geom_text(aes(label = region), hjust = -0.4, size = 3.5, fontface = "bold") +
    scale_x_continuous(limits = c(0, 1), name = "w_ate") +
    scale_y_continuous(limits = c(0, 1), name = "w_imp") +
    labs(title = "(d) Adaptive weight distribution by region",
      subtitle = "Mean (w_dag, w_ate, w_imp) per region; w_dag = 1 - w_ate - w_imp") +
    theme_pub() + theme(legend.position = "none")
} else {
  pd <- ggplot() + theme_void() + labs(title = "(d) Adaptive weights (no w_* in screening)")
}

# Combine
fig3 <- (pa | pb) / (pc | pd) +
  plot_annotation(
    title = "Fig 3  Cross-species causal consistency",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

ggsave(file.path(fig_dir, "fig3_cross_species_causal_consistency.png"),
  fig3, width = 14, height = 10, dpi = 1200, bg = "white")
ggsave(file.path(fig_dir, "fig3_cross_species_causal_consistency.svg"),
  fig3, width = 14, height = 10, bg = "white")
cat("Saved fig3_cross_species_causal_consistency.png and .svg to", fig_dir, "\n")
