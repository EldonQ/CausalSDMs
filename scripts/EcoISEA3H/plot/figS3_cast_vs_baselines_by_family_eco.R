################################################################################
# Fig S3 (Eco): CAST vs Baselines — Per-family detailed comparison
# Requires: all_results_v3.csv, CAST_Species_Summary.csv
# Run: setwd("E:/CausalSDMs"); source("scripts/EcoISEA3H/plot/figS3_cast_vs_baselines_by_family_eco.R")
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

d <- read.csv("output/case2_eco/all_results_v3.csv", stringsAsFactors = FALSE) %>%
    filter(!is.na(auc_mean))

sp_meta <- read.csv("outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv", stringsAsFactors = FALSE) %>%
    mutate(species = gsub(" ", "_", species))

d <- d %>%
    left_join(sp_meta %>% select(species, family), by = "species") %>%
    mutate(family = ifelse(is.na(family), "Unknown", family))

# Remove groups with less than 2 species to avoid empty boxplots or singular lines
fam_counts <- d %>%
    group_by(family) %>%
    summarise(n = n_distinct(species))
valid_fams <- fam_counts %>%
    filter(n >= 2) %>%
    pull(family)
d <- d %>% filter(family %in% valid_fams)

model_order <- c("CAST", "MLP_ATE", "MLP", "RF", "BRT", "Maxent")
model_colors <- c(
    "CAST" = "#2980B9", "MLP_ATE" = "#5DADE2", "MLP" = "#85C1E9",
    "RF" = "#27AE60", "BRT" = "#E67E22", "Maxent" = "#9B59B6"
)
d$model <- factor(d$model, levels = model_order)

p <- ggplot(d, aes(x = model, y = auc_mean, fill = model)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 2.5, color = "black") +
    facet_wrap(~family, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = model_colors, guide = "none") +
    labs(
        title = "Fig S3 (Eco) CAST vs Baselines by Taxonomic Family",
        subtitle = "Boxplot per model-family | Diamond = mean | Excludes singletons",
        x = "", y = "AUC"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "figS3_cast_vs_baselines_by_family_eco.png"),
    p,
    width = 12, height = 8, dpi = 300, bg = "white"
)
cat("Saved figS3_cast_vs_baselines_by_family_eco.png\n")
