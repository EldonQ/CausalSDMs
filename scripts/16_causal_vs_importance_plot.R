# ==============================================================================
# Script Name: 16_causal_vs_importance_plot.R
# Function: Generate "Lie Detector" plot (Importance vs ATE)
# Analysis: Identify variables that are statistically important but causally insignificant
# ==============================================================================

# Initialize
rm(list = ls())
setwd("E:/CausalSDMs")
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(sysfonts)
library(showtext)

# Font setup for Nature-style figures
try(
    {
        sysfonts::font_add(
            family = "Arial",
            regular = "C:/Windows/Fonts/arial.ttf",
            bold = "C:/Windows/Fonts/arialbd.ttf",
            italic = "C:/Windows/Fonts/ariali.ttf",
            bolditalic = "C:/Windows/Fonts/arialbi.ttf"
        )
        showtext::showtext_opts(dpi = 600)
        showtext::showtext_auto(enable = TRUE)
    },
    silent = TRUE
)

# Create output directory
dir.create("figures/16_causal_comparison", showWarnings = FALSE, recursive = TRUE)

# 1. Load Data
# ------------------------------------------------------------------------------
# ATE Results (Causal)
ate_df <- read.csv("output/14_causal/ate_all_variables.csv", stringsAsFactors = FALSE) %>%
    select(variable, ate_coef = coef, p_value, significant)

# Importance Results (Correlation-based)
imp_df <- read.csv("output/09_variable_importance/importance_summary.csv", stringsAsFactors = FALSE) %>%
    group_by(variable) %>%
    summarise(importance = mean(importance_normalized, na.rm = TRUE)) # Use RF importance

# 2. Merge and Process
# ------------------------------------------------------------------------------
# Normalize function
min_max_norm <- function(x) {
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
}

merged_df <- imp_df %>%
    inner_join(ate_df, by = "variable") %>%
    mutate(
        # Absolute ATE for magnitude comparison
        ate_abs = abs(ate_coef),
        # Normalize both to 0-1 scale for relative comparison
        imp_norm = min_max_norm(importance),
        ate_norm = min_max_norm(ate_abs),
        # Categorize variables
        category = case_when(
            ate_norm > 0.5 & imp_norm > 0.5 ~ "True Drivers (High Causal + High Imp)",
            ate_norm < 0.2 & imp_norm > 0.5 ~ "Spurious/Confounders (Low Causal + High Imp)",
            ate_norm > 0.5 & imp_norm < 0.2 ~ "Masked Drivers (High Causal + Low Imp)",
            TRUE ~ "Other"
        ),
        # Significance label
        sig_label = ifelse(significant, "Significant", "Non-significant")
    )

# Remove huge outliers if any (based on previous view)
merged_df <- merged_df %>% filter(ate_abs < 100) # outlier filter for that one huge value

# 3. Visualization
# ------------------------------------------------------------------------------
p <- ggplot(merged_df, aes(x = imp_norm, y = ate_norm)) +
    # Background zones
    annotate("rect",
        xmin = 0.5, xmax = 1, ymin = 0, ymax = 0.3,
        fill = "#FFDB6D", alpha = 0.2
    ) + # The "Lie" Zone
    annotate("text",
        x = 0.85, y = 0.15, label = "Spurious Correlations\n(High Imp, Low ATE)",
        color = "#D18800", fontface = "italic", size = 3
    ) +

    # Scatter points
    geom_point(aes(color = significant, size = imp_norm), alpha = 0.8) +

    # Diagonal reference
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey80") +

    # Labels for interesting points
    geom_text_repel(
        data = subset(merged_df, (imp_norm > 0.6 | ate_norm > 0.6) | (imp_norm > 0.5 & ate_norm < 0.3)),
        aes(label = variable),
        size = 3, fontface = "bold", box.padding = 0.5, max.overlaps = 20
    ) +

    # Styles
    scale_color_manual(values = c("Non-significant" = "grey60", "Significant" = "#E41A1C")) +
    scale_size(range = c(2, 6)) +
    labs(
        title = "The Lie Detector: Statistical Importance vs. Causal Effect",
        subtitle = "Variables in the yellow zone are 'Spurious Correlations' driven by confounders",
        x = "Random Forest Importance (Normalized)",
        y = "Causal Effect Magnitude (Normalized ATE)",
        caption = paste0("Data: ", nrow(merged_df), " variables. Significance threshold: P < 0.05 (FDR adj)")
    ) +
    theme_minimal() +
    theme(
        text = element_text(family = "Arial"),
        legend.position = "bottom",
        panel.grid.minor = element_blank()
    )

# Save
ggsave("figures/16_causal_comparison/lie_detector_plot.png", p, width = 8, height = 6, dpi = 600, bg = "white")
ggsave("figures/16_causal_comparison/lie_detector_plot.pdf", p, width = 8, height = 6)

# 4. Output Summary Table for User
# ------------------------------------------------------------------------------
liars <- merged_df %>%
    filter(imp_norm > 0.4 & ate_norm < 0.2) %>%
    arrange(desc(imp_norm)) %>%
    select(variable, RF_Importance = importance, Causal_ATE = ate_coef, P_value = p_value)

write.csv(liars, "output/15b_causal_retraining/identified_spurious_vars.csv", row.names = FALSE)

cat("==============================================================\n")
cat("Found", nrow(liars), "potential spurious variables (High Imp, Low ATE):\n")
print(head(liars))
cat("==============================================================\n")
