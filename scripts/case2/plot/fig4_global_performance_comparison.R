################################################################################
# Fig 4: Global performance comparison (4-panel)
# (a) Violin + dot by model  (b) Bump chart rank by region  (c) CAST vs baselines scatter  (d) Delta AUC CDF
# Requires: all_results_v3.csv
# Run from project root: setwd("E:/CausalSDMs"); source("scripts/case2/plot/fig4_global_performance_comparison.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

fig_dir <- "figures/case2/plot"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

path_res <- "output/case2/all_results_v3.csv"
if (!file.exists(path_res)) stop("Missing ", path_res)

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

d <- read.csv(path_res, stringsAsFactors = FALSE)
d <- d %>% filter(!is.na(auc_mean))
model_order <- c("CAST", "MLP_ATE", "MLP", "RF", "BRT", "Maxent", "GAM")
d$model <- factor(d$model, levels = rev(model_order))

# ---- (a) Violin + mean dot by model ----
pa <- ggplot(d, aes(x = model, y = auc_mean, fill = model)) +
  geom_violin(scale = "width", alpha = 0.6, trim = TRUE) +
  stat_summary(fun = mean, geom = "point", size = 3, color = "black", shape = 18) +
  scale_fill_manual(values = setNames(c("#2980B9", "grey75", "grey85", "#27AE60", "#E67E22", "#9B59B6", "grey70"), model_order), guide = "none") +
  coord_flip() +
  labs(title = "(a) AUC distribution by model", subtitle = "Black diamond = mean; all species", x = "", y = "AUC") +
  theme_pub()

# ---- (b) Bump chart: rank by region ----
rank_by_region <- d %>%
  group_by(region, species) %>%
  mutate(rank = dense_rank(desc(auc_mean))) %>%
  ungroup() %>%
  group_by(region, model) %>%
  summarise(mean_rank = mean(rank), .groups = "drop")
pb <- ggplot(rank_by_region, aes(x = region, y = mean_rank, color = model, group = model)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  scale_y_reverse(breaks = 1:7, name = "Mean rank (1 = best)") +
  scale_color_manual(values = setNames(c("#2980B9", "grey50", "grey60", "#27AE60", "#E67E22", "#9B59B6", "grey50"), model_order), name = "Model") +
  labs(title = "(b) Mean rank by region", subtitle = "Bump chart: lower = better", x = "Region") +
  theme_pub() + theme(legend.position = "right")

# ---- (c) CAST vs each baseline scatter (one panel: CAST vs best baseline) ----
cast_wide <- d %>%
  select(region, species, model, auc_mean) %>%
  pivot_wider(names_from = model, values_from = auc_mean)
if ("CAST" %in% names(cast_wide)) {
  other_cols <- setdiff(names(cast_wide), c("region", "species", "CAST"))
  cast_wide$best_other <- apply(cast_wide[, other_cols, drop = FALSE], 1, function(x) max(x, na.rm = TRUE))
  pc <- ggplot(cast_wide, aes(x = best_other, y = CAST)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    geom_point(aes(color = region), alpha = 0.7, size = 2.5) +
    scale_color_brewer(palette = "Set2", name = "Region") +
    labs(title = "(c) CAST vs best baseline (per species)", subtitle = "Above line = CAST wins", x = "Best baseline AUC", y = "CAST AUC") +
    theme_pub()
} else {
  pc <- ggplot() + theme_void() + labs(title = "(c) CAST vs baseline (no CAST in results)")
}

# ---- (d) CDF of Delta AUC (CAST - best baseline) ----
if (exists("cast_wide") && "CAST" %in% names(cast_wide)) {
  cast_wide$delta <- cast_wide$CAST - cast_wide$best_other
  pd <- ggplot(cast_wide, aes(x = delta)) +
    stat_ecdf(geom = "step", linewidth = 1.2, color = "#2980B9") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    annotate("text", x = mean(cast_wide$delta), y = 0.5, label = sprintf("Mean %+.4f", mean(cast_wide$delta)), fontface = "bold") +
    labs(title = "(d) CDF of CAST minus best baseline", subtitle = "Fraction of species with AUC gain above x", x = expression(Delta * "AUC"), y = "Cumulative fraction") +
    theme_pub()
} else {
  pd <- ggplot() + theme_void() + labs(title = "(d) Delta AUC CDF")
}

fig4 <- (pa | pb) / (pc | pd) +
  plot_annotation(
    title = "Fig 4  Global performance comparison",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

ggsave(file.path(fig_dir, "fig4_global_performance_comparison.png"), fig4, width = 14, height = 10, dpi = 1200, bg = "white")
ggsave(file.path(fig_dir, "fig4_global_performance_comparison.svg"), fig4, width = 14, height = 10, bg = "white")
cat("Saved fig4_global_performance_comparison.png and .svg to", fig_dir, "\n")
