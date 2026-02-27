################################################################################
# Fig 6: Ecological interpretability (3-panel)
# (a) Correlation vs causal rank: slope chart  (b) CATE placeholder or note  (c) Variable flow: Sankey/Alluvial
# Requires: all_ate_results_v3.csv, all_screening_v3.csv, all_dag_edges_v3.csv
# Run from project root: setwd("E:/CausalSDMs"); source("scripts/case2/plot/fig6_ecological_interpretability.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

fig_dir <- "figures/case2/plot"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_ate    <- "output/case2/all_ate_results_v3.csv"
path_screen <- "output/case2/all_screening_v3.csv"
path_dag    <- "output/case2/all_dag_edges_v3.csv"

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

# ---- (a) Correlation (RF importance) vs Causal (|ATE|) rank ----
# Use screening: importance = RF importance; ate = |coef| from ate. Merge by variable and species.
if (file.exists(path_ate) && file.exists(path_screen)) {
  ate    <- read.csv(path_ate, stringsAsFactors = FALSE)
  screen <- read.csv(path_screen, stringsAsFactors = FALSE)
  # Per-species rank: importance rank and |ATE| rank
  screen <- screen %>%
    group_by(region, species) %>%
    mutate(rank_imp = row_number(desc(importance)), rank_ate = row_number(desc(abs(coef)))) %>%
    ungroup()
  # Aggregate: mean rank across species per variable (or show one region)
  rank_agg <- screen %>%
    group_by(region, variable) %>%
    summarise(mean_rank_imp = mean(rank_imp, na.rm = TRUE), mean_rank_ate = mean(rank_ate, na.rm = TRUE), .groups = "drop")
  slope_dat <- rank_agg %>%
    filter(region == rank_agg$region[1]) %>%
    mutate(
      rank_imp_n = rank(mean_rank_imp),
      rank_ate_n = rank(mean_rank_ate),
      y_left  = max(rank_imp_n) - rank_imp_n + 1,
      y_right = max(rank_ate_n) - rank_ate_n + 1
    )
  # Slope chart: x=1 (RF rank) to x=2 (ATE rank), y = rank position
  pa <- ggplot(slope_dat) +
    geom_segment(aes(x = 1, xend = 2, y = y_left, yend = y_right, color = variable), linewidth = 1.2) +
    geom_text(aes(x = 0.98, y = y_left, label = variable), hjust = 1, size = 3) +
    geom_text(aes(x = 2.02, y = y_right, label = variable), hjust = 0, size = 3) +
    scale_x_continuous(limits = c(0.5, 2.5), breaks = c(1, 2), labels = c("RF importance rank", "|ATE| rank")) +
    scale_color_brewer(palette = "Set3", guide = "none") +
    labs(title = "(a) Correlation vs causal rank (one region)", subtitle = "Left = RF importance rank; Right = |ATE| rank", x = "", y = "Rank position") +
    theme_pub()
} else {
  pa <- ggplot() + theme_void() + labs(title = "(a) Rank comparison (data missing)")
}

# ---- (b) CATE spatial note (no coordinates in disdat by default) ----
pb <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "CATE maps require spatial coordinates.\nDisdat benchmark: use region-specific\ncoordinates if available.", size = 4, hjust = 0.5) +
  theme_void() +
  labs(title = "(b) CATE spatial heterogeneity", subtitle = "Placeholder: add maps when coordinates are available")

# ---- (c) Variable flow: variable -> role -> selected (bar or alluvial) ----
if (file.exists(path_screen)) {
  screen <- read.csv(path_screen, stringsAsFactors = FALSE)
  if ("w_dag" %in% names(screen)) {
    flow <- screen %>%
      group_by(region, variable) %>%
      summarise(
        mean_score = mean(score_total), mean_w_dag = mean(w_dag), mean_w_ate = mean(w_ate), mean_w_imp = mean(w_imp),
        selected_pct = mean(score_total >= quantile(score_total, 1 - 0.5)), .groups = "drop"
      ) %>%
      mutate(selected = selected_pct > 0.5)
    pc <- ggplot(flow %>% filter(region == flow$region[1]), aes(x = reorder(variable, -mean_score), y = mean_score, fill = selected)) +
      geom_col(width = 0.7, alpha = 0.85) +
      scale_fill_manual(values = c("TRUE" = "#2980B9", "FALSE" = "grey75"), name = "CAST-selected") +
      labs(title = "(c) Variable screening score (one region)", subtitle = "Mean composite score; blue = typically selected", x = "", y = "Mean score") +
      theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
  } else {
    pc <- ggplot() + theme_void() + labs(title = "(c) Variable flow (no w_* in screening)")
  }
} else {
  pc <- ggplot() + theme_void() + labs(title = "(c) Variable flow")
}

fig6 <- (pa | pb) / pc +
  plot_annotation(
    title = "Fig 6  Ecological interpretability",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

ggsave(file.path(fig_dir, "fig6_ecological_interpretability.png"), fig6, width = 12, height = 8, dpi = 1200, bg = "white")
ggsave(file.path(fig_dir, "fig6_ecological_interpretability.svg"), fig6, width = 12, height = 8, bg = "white")
cat("Saved fig6_ecological_interpretability.png and .svg to", fig_dir, "\n")
