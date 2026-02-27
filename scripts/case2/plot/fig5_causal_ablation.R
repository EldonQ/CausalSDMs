################################################################################
# Fig 5: Causal enhancement ablation (3-panel)
# (a) Waterfall: MLP -> MLP_ATE -> CAST  (b) Step-wise AUC distribution  (c) DAG density vs CAST gain
# Requires: all_results_v3.csv (with MLP, MLP_ATE, CAST)
# Run from project root: setwd("E:/CausalSDMs"); source("scripts/case2/plot/fig5_causal_ablation.R")
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
d <- d %>% filter(model %in% c("MLP", "MLP_ATE", "CAST"), !is.na(auc_mean))

# Wide: one row per (region, species), columns MLP, MLP_ATE, CAST, dag_density
ablation_wide <- d %>%
  select(region, species, model, auc_mean, dag_density) %>%
  pivot_wider(names_from = model, values_from = auc_mean, values_fn = mean)
ablation_wide <- ablation_wide %>%
  filter(complete.cases(ablation_wide %>% select(any_of(c("MLP", "MLP_ATE", "CAST")))))

# ---- (a) Waterfall: mean AUC at each step ----
means <- c(
  MLP     = mean(ablation_wide$MLP, na.rm = TRUE),
  MLP_ATE = mean(ablation_wide$MLP_ATE, na.rm = TRUE),
  CAST    = mean(ablation_wide$CAST, na.rm = TRUE)
)
delta_ate  <- means["MLP_ATE"] - means["MLP"]
delta_dag  <- means["CAST"] - means["MLP_ATE"]
wf <- data.frame(
  step = factor(c("MLP (baseline)", "+ ATE weighting", "CAST (+ DAG interactions)"), levels = c("MLP (baseline)", "+ ATE weighting", "CAST (+ DAG interactions)")),
  auc  = c(means["MLP"], means["MLP_ATE"], means["CAST"]),
  delta = c(NA, delta_ate, delta_dag)
)
pa <- ggplot(wf, aes(x = step, y = auc, fill = step)) +
  geom_col(width = 0.6, alpha = 0.85) +
  geom_text(aes(label = sprintf("%.4f", auc)), vjust = -0.5, size = 3.5, fontface = "bold") +
  geom_segment(data = wf %>% filter(!is.na(delta)), aes(x = as.numeric(step) - 0.3, xend = as.numeric(step) + 0.3, y = auc, yend = auc), linewidth = 0.5) +
  scale_fill_manual(values = c("MLP (baseline)" = "grey60", "+ ATE weighting" = "#3498DB", "CAST (+ DAG interactions)" = "#2980B9"), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "(a) Ablation: mean AUC by step", subtitle = sprintf("ATE step: %+.4f  |  DAG step: %+.4f", delta_ate, delta_dag), x = "", y = "Mean AUC") +
  theme_pub() + theme(axis.text.x = element_text(angle = 25, hjust = 1))

# ---- (b) Step-wise AUC distribution (violin/box) ----
long <- d %>% filter(model %in% c("MLP", "MLP_ATE", "CAST")) %>%
  mutate(model = factor(model, levels = c("MLP", "MLP_ATE", "CAST")))
pb <- ggplot(long, aes(x = model, y = auc_mean, fill = model)) +
  geom_violin(alpha = 0.6, trim = TRUE) +
  geom_boxplot(width = 0.15, outlier.fill = NA, fill = "white", alpha = 0.7) +
  scale_fill_manual(values = c("MLP" = "grey60", "MLP_ATE" = "#3498DB", "CAST" = "#2980B9"), name = "") +
  labs(title = "(b) AUC distribution per step", subtitle = "All species", x = "", y = "AUC") +
  theme_pub()

# ---- (c) DAG density vs CAST - MLP gain ----
ablation_wide <- ablation_wide %>% mutate(delta_cast_mlp = CAST - MLP)
if ("dag_density" %in% names(ablation_wide) && all(!is.na(ablation_wide$dag_density))) {
  pc <- ggplot(ablation_wide, aes(x = dag_density, y = delta_cast_mlp)) +
    geom_point(aes(color = region), alpha = 0.7, size = 2.5) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    scale_color_brewer(palette = "Set2", name = "Region") +
    labs(title = "(c) DAG density vs CAST gain over MLP", subtitle = "Each point = one species; trend line", x = "DAG density", y = expression(Delta * "AUC (CAST - MLP)")) +
    theme_pub()
} else {
  pc <- ggplot() + theme_void() + labs(title = "(c) DAG density vs gain (dag_density missing)")
}

fig5 <- (pa | pb) / pc +
  plot_annotation(
    title = "Fig 5  Causal enhancement ablation",
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

ggsave(file.path(fig_dir, "fig5_causal_ablation.png"), fig5, width = 12, height = 8, dpi = 1200, bg = "white")
ggsave(file.path(fig_dir, "fig5_causal_ablation.svg"), fig5, width = 12, height = 8, bg = "white")
cat("Saved fig5_causal_ablation.png and .svg to", fig_dir, "\n")
