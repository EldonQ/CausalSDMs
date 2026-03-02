################################################################################
# Fig S4: DAG Density vs CAST Gain (diagnostic)
# Requires: all_results_v3.csv (with dag_density)
# Run: setwd("E:/CausalSDMs"); source("scripts/case2/plot/figS4_dag_density_vs_gain.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

fig_dir <- "figures/case2/plot"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)

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

abl_wide <- d %>%
    filter(model %in% c("MLP", "CAST")) %>%
    select(region, species, model, auc_mean, dag_density) %>%
    pivot_wider(names_from = model, values_from = auc_mean, values_fn = max) %>%
    filter(!is.na(CAST), !is.na(MLP)) %>%
    mutate(delta = CAST - MLP)

# Fill dag_density where missing
abl_wide <- abl_wide %>%
    left_join(
        d %>% filter(!is.na(dag_density)) %>%
            distinct(region, species, dag_density),
        by = c("region", "species")
    ) %>%
    mutate(dag_density = coalesce(dag_density.x, dag_density.y)) %>%
    select(-dag_density.x, -dag_density.y)

p <- ggplot(
    abl_wide %>% filter(!is.na(dag_density)),
    aes(x = dag_density, y = delta)
) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    geom_point(aes(color = region), size = 2.5, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 0.6) +
    scale_color_brewer(palette = "Set2", name = "Region") +
    labs(
        title = "Fig S4  DAG Density vs CAST Gain over MLP",
        subtitle = "Each point = one species | Positive = CAST wins",
        x = "DAG density (proportion of possible edges)",
        y = expression(Delta * "AUC (CAST - MLP)")
    ) +
    theme_pub()

ggsave(file.path(fig_dir, "figS4_dag_density_vs_gain.png"),
    p,
    width = 8, height = 6, dpi = 300, bg = "white"
)
ggsave(file.path(fig_dir, "figS4_dag_density_vs_gain.svg"),
    p,
    width = 8, height = 6, bg = "white"
)
cat("Saved figS4\n")
