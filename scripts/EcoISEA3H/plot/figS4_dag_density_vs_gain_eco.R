################################################################################
# Fig S4 (Eco): DAG Density vs CAST Gain (diagnostic)
# Requires: all_results_v3.csv, CAST_Species_Summary.csv
# Run: setwd("E:/CausalSDMs"); source("scripts/EcoISEA3H/plot/figS4_dag_density_vs_gain_eco.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

fig_dir <- "figures/case2_eco/plot"
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

d <- read.csv("output/case2_eco/all_results_v3.csv", stringsAsFactors = FALSE)
sp_meta <- read.csv("outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv", stringsAsFactors = FALSE) %>%
    mutate(species = gsub(" ", "_", species))

abl_wide <- d %>%
    filter(model %in% c("MLP", "CAST")) %>%
    select(species, model, auc_mean, dag_density) %>%
    pivot_wider(names_from = model, values_from = auc_mean, values_fn = max) %>%
    filter(!is.na(CAST), !is.na(MLP)) %>%
    mutate(delta = CAST - MLP)

# Fill dag_density where missing
abl_wide <- abl_wide %>%
    left_join(
        d %>% filter(!is.na(dag_density)) %>%
            distinct(species, dag_density),
        by = "species"
    ) %>%
    mutate(dag_density = coalesce(dag_density.x, dag_density.y)) %>%
    select(-dag_density.x, -dag_density.y) %>%
    left_join(sp_meta %>% select(species, family, category), by = "species") %>%
    mutate(family = ifelse(is.na(family), "Unknown", family))

p <- ggplot(
    abl_wide %>% filter(!is.na(dag_density)),
    aes(x = dag_density, y = delta)
) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey50") +
    geom_point(aes(color = family), size = 2.5, alpha = 0.7) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed", linewidth = 0.6) +
    scale_color_brewer(palette = "Set1", name = "Taxonomic\nFamily") +
    labs(
        title = "Fig S4 (Eco) DAG Density vs CAST Gain over MLP",
        subtitle = "Each point = one species | Positive = CAST wins",
        x = "DAG density (proportion of possible edges)",
        y = expression(Delta * "AUC (CAST - MLP)")
    ) +
    theme_pub()

ggsave(file.path(fig_dir, "figS4_dag_density_vs_gain_eco.png"),
    p,
    width = 8, height = 6, dpi = 300, bg = "white"
)
cat("Saved figS4_dag_density_vs_gain_eco.png\n")
