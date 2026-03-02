################################################################################
# Fig S3: CAST vs Baselines — Per-region detailed comparison (6-panel faceted)
# Requires: all_results_v3.csv
# Run: setwd("E:/CausalSDMs"); source("scripts/case2/plot/figS3_cast_vs_baselines_by_region.R")
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

d <- read.csv("output/case2/all_results_v3.csv", stringsAsFactors = FALSE) %>%
    filter(!is.na(auc_mean))

model_order <- c("CAST", "MLP_ATE", "MLP", "RF", "BRT", "Maxent")
model_colors <- c(
    "CAST" = "#2980B9", "MLP_ATE" = "#5DADE2", "MLP" = "#85C1E9",
    "RF" = "#27AE60", "BRT" = "#E67E22", "Maxent" = "#9B59B6"
)
d$model <- factor(d$model, levels = model_order)

p <- ggplot(d, aes(x = model, y = auc_mean, fill = model)) +
    geom_boxplot(alpha = 0.7, outlier.size = 1) +
    stat_summary(fun = mean, geom = "point", shape = 18, size = 2.5, color = "black") +
    facet_wrap(~region, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = model_colors, guide = "none") +
    labs(
        title = "Fig S3  CAST vs Baselines by Region",
        subtitle = "Boxplot per model-region | Diamond = mean",
        x = "", y = "AUC"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(file.path(fig_dir, "figS3_cast_vs_baselines_by_region.png"),
    p,
    width = 12, height = 8, dpi = 300, bg = "white"
)
ggsave(file.path(fig_dir, "figS3_cast_vs_baselines_by_region.svg"),
    p,
    width = 12, height = 8, bg = "white"
)
cat("Saved figS3\n")
