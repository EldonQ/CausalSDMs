################################################################################
# Fig 2: CAST single-species causal discovery showcase (4-panel)
# (a) DAG network  (b) ATE forest plot  (c) Screening scores  (d) SubDAG focus
# Run from project root: setwd("E:/CausalSDMs"); source("scripts/case2/plot/fig2_cast_single_species_showcase.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

target_species <- "swi10"
target_region  <- "SWI"

# Paths
result_base <- sprintf("output/case2/%s/single_species", target_region)
fig_dir     <- "figures/case2/plot"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Dependencies
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork")
if (!requireNamespace("igraph", quietly = TRUE)) install.packages("igraph")
if (!requireNamespace("ggraph", quietly = TRUE)) install.packages("ggraph")
library(tidyverse)
library(ggplot2)
library(patchwork)
library(igraph)
library(ggraph)

# Theme: Arial, publication-ready, no Chinese
theme_pub <- function(base_size = 11) {
  theme_minimal(base_size = base_size, base_family = "sans") +
    theme(
      panel.grid.minor = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold", hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40"),
      legend.position = "right"
    )
}

# Load data
results   <- read.csv(sprintf("%s/results_%s.csv", result_base, target_species), stringsAsFactors = FALSE)
screening <- read.csv(sprintf("%s/screening_%s.csv", result_base, target_species), stringsAsFactors = FALSE)
ate_raw   <- read.csv(sprintf("%s/ate_%s.csv", result_base, target_species), stringsAsFactors = FALSE)
roles_file <- sprintf("%s/roles_%s.csv", result_base, target_species)
dag_file   <- sprintf("%s/dag_edges_%s.csv", result_base, target_species)

has_roles <- file.exists(roles_file)
roles     <- if (has_roles) read.csv(roles_file, stringsAsFactors = FALSE) else NULL
dag_edges <- if (file.exists(dag_file)) read.csv(dag_file, stringsAsFactors = FALSE) else data.frame()

# ATE 95% CI
ate_data <- ate_raw %>%
  mutate(
    ci_lower = coef - 1.96 * se,
    ci_upper = coef + 1.96 * se,
    se = replace_na(se, 0)
  )

# Cast-selected variables (top by score_total or from roles)
cast_vars <- if (has_roles) roles$variable else screening$variable[1:min(6, nrow(screening))]

# Colour palette
pal <- list(
  root = "#2C3E50", med = "#E67E22", term = "#27AE60",
  sig = "#C0392B", nsig = "grey60", cast = "#2980B9"
)

# ---- Panel (a): DAG network ----
pa <- ggplot() + theme_void() + labs(title = "(a) Consensus DAG")
if (nrow(dag_edges) > 0 && requireNamespace("ggraph", quietly = TRUE)) {
  g <- graph_from_data_frame(dag_edges %>% select(from, to), directed = TRUE,
    vertices = data.frame(name = unique(c(dag_edges$from, dag_edges$to))))
  V(g)$role <- if (has_roles) roles$group[match(V(g)$name, roles$variable)] else "Other"
  V(g)$role[is.na(V(g)$role)] <- "Other"
  V(g)$is_cast <- V(g)$name %in% cast_vars
  E(g)$strength <- dag_edges$strength
  role_colors <- c("Root" = pal$root, "Mediator" = pal$med, "Terminal" = pal$term, "Other" = "grey70")
  pa <- ggraph(g, layout = "fr") +
    geom_edge_link(aes(edge_width = strength, edge_alpha = strength),
      arrow = arrow(length = unit(2.5, "mm"), type = "closed"), color = "grey35") +
    geom_node_point(aes(color = role, size = is_cast), alpha = 0.9) +
    geom_node_text(aes(label = name), repel = TRUE, size = 3, fontface = "bold") +
    scale_color_manual(values = role_colors, name = "Role") +
    scale_size_manual(values = c("TRUE" = 8, "FALSE" = 5), guide = "none") +
    scale_edge_width(range = c(0.4, 2)) +
    scale_edge_alpha(range = c(0.4, 0.95), guide = "none") +
    labs(subtitle = sprintf("%d edges | density %.2f", nrow(dag_edges), results$dag_density[1])) +
    theme_pub() + theme(legend.position = "right", axis.text = element_blank(), panel.grid = element_blank())
}

# ---- Panel (b): ATE forest plot ----
ate_plot <- ate_data %>%
  mutate(sig_marker = ifelse(significant, "Significant (p<0.05)", "Not significant"))
pb <- ggplot(ate_plot, aes(x = coef, y = reorder(variable, abs(coef)), color = sig_marker)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 2.5) +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), width = 0.2, linewidth = 0.5) +
  scale_color_manual(values = c("Significant (p<0.05)" = pal$sig, "Not significant" = pal$nsig), name = "") +
  labs(title = "(b) Average Treatment Effect", subtitle = "DML 2-fold | 95% CI", x = "ATE", y = "") +
  theme_pub() + theme(panel.grid.major.y = element_line(color = "grey92"))

# ---- Panel (c): Screening scores (stacked bar) ----
w_dag <- 0.15 + 0.15 * (1 - results$dag_density[1])
w_ate <- 0.25 + 0.25 * (sum(ate_data$significant) / nrow(ate_data))
w_imp <- 1 - w_dag - w_ate
screen_long <- screening %>%
  select(variable, score_dag, score_ate, score_imp, score_total) %>%
  mutate(
    w_dag = w_dag * score_dag,
    w_ate = w_ate * score_ate,
    w_imp = w_imp * score_imp
  ) %>%
  pivot_longer(cols = c(w_dag, w_ate, w_imp), names_to = "component", values_to = "score") %>%
  mutate(component = factor(component, levels = c("w_dag", "w_ate", "w_imp"),
    labels = c("DAG", "ATE", "RF")))
pc <- ggplot(screen_long, aes(x = reorder(variable, -score_total), y = score, fill = component)) +
  geom_col(position = "stack", width = 0.7, alpha = 0.85) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = pal$cast, linewidth = 0.6) +
  scale_fill_manual(values = c("DAG" = "#2C3E50", "ATE" = "#E67E22", "RF" = "#27AE60"), name = "Component") +
  labs(title = "(c) Adaptive screening score", subtitle = sprintf("Selected: %s", paste(cast_vars, collapse = ", ")),
    x = "", y = "Weighted score") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))

# ---- Panel (d): SubDAG (screened vars only) ----
pd <- ggplot() + theme_void() + labs(title = "(d) DAG among screened variables")
if (nrow(dag_edges) > 0 && length(cast_vars) > 0) {
  sub_edges <- dag_edges %>%
    filter(from %in% cast_vars, to %in% cast_vars)
  if (nrow(sub_edges) > 0) {
    gs <- graph_from_data_frame(sub_edges %>% select(from, to), directed = TRUE,
      vertices = data.frame(name = cast_vars))
    E(gs)$strength <- sub_edges$strength
    pd <- ggraph(gs, layout = "fr") +
      geom_edge_link(aes(edge_width = strength), arrow = arrow(length = unit(2.5, "mm"), type = "closed"),
        color = pal$cast, alpha = 0.8) +
      geom_node_point(size = 9, color = pal$cast, fill = "white", shape = 21, stroke = 1.5) +
      geom_node_text(aes(label = name), repel = TRUE, size = 3.2, fontface = "bold") +
      scale_edge_width(range = c(0.6, 2.2)) +
      labs(subtitle = sprintf("%d edges among %d screened vars", nrow(sub_edges), length(cast_vars))) +
      theme_pub() + theme(legend.position = "none", axis.text = element_blank(), panel.grid = element_blank())
  }
}

# Combine
fig2 <- (pa | pb) / (pc | pd) +
  plot_annotation(
    title = sprintf("Fig 2  CAST causal discovery: %s [%s]", target_species, target_region),
    theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
  )

ggsave(file.path(fig_dir, "fig2_cast_single_species_showcase.png"),
  fig2, width = 14, height = 10, dpi = 1200, bg = "white")
ggsave(file.path(fig_dir, "fig2_cast_single_species_showcase.svg"),
  fig2, width = 14, height = 10, bg = "white")
cat("Saved fig2_cast_single_species_showcase.png and .svg to", fig_dir, "\n")
