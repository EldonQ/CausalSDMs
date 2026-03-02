################################################################################
# Fig S5: Supplementary — DAG Density Distribution & ATE Effects by Region
# (a) DAG density violin by region  (b) ATE direction ratio per region
#
# This replaces old figS1 + figS2 into one combined supplementary figure.
# Requires: all_results_v3.csv, all_ate_results_v3.csv
# Run: setwd("E:/CausalSDMs"); source("scripts/case2/plot/figS5_dag_and_ate_by_region.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

fig_dir <- "figures/case2/plot"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
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

d <- read.csv("output/case2/all_results_v3.csv", stringsAsFactors = FALSE)
ate <- read.csv("output/case2/all_ate_results_v3.csv", stringsAsFactors = FALSE)

# ── (a) DAG density distribution by region ────────────────
dag_density <- d %>%
    filter(!is.na(dag_density)) %>%
    distinct(region, species, dag_density)

pa <- ggplot(dag_density, aes(x = region, y = dag_density, fill = region)) +
    geom_violin(alpha = 0.6, scale = "width") +
    geom_boxplot(width = 0.15, fill = "white", alpha = 0.8, outlier.size = 0.8) +
    geom_jitter(width = 0.1, alpha = 0.4, size = 1.2) +
    scale_fill_brewer(palette = "Set2", guide = "none") +
    labs(
        title = "(a) DAG density by region",
        subtitle = "Each point = one species | Higher density = more edges",
        x = "", y = "DAG density"
    ) +
    theme_pub()

# ── (b) ATE effect direction ratio per region ─────────────
ate_dir <- ate %>%
    filter(!is.na(coef)) %>%
    mutate(direction = ifelse(coef > 0, "Positive", "Negative")) %>%
    group_by(region, direction) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(region) %>%
    mutate(pct = n / sum(n))

pb <- ggplot(ate_dir, aes(x = region, y = pct, fill = direction)) +
    geom_col(position = "stack", width = 0.65, alpha = 0.85) +
    geom_text(aes(label = n), position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold") +
    scale_fill_manual(values = c("Positive" = "#27AE60", "Negative" = "#C0392B"), name = "ATE direction") +
    scale_y_continuous(labels = scales::percent) +
    labs(
        title = "(b) ATE effect direction by region",
        subtitle = "Proportion of positive vs negative effects",
        x = "", y = "Proportion"
    ) +
    theme_pub() +
    theme(legend.position = "right")

# Combine
figS5 <- pa / pb +
    plot_annotation(
        title = "Fig S5  DAG Structure & ATE Effects by Region",
        theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    )

ggsave(file.path(fig_dir, "figS5_dag_and_ate_by_region.png"),
    figS5,
    width = 10, height = 10, dpi = 300, bg = "white"
)
ggsave(file.path(fig_dir, "figS5_dag_and_ate_by_region.svg"),
    figS5,
    width = 10, height = 10, bg = "white"
)
cat("Saved figS5\n")
