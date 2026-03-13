################################################################################
# Fig S2: DAG Complexity & Ecological Effect Heterogeneity (Supplementary)
#
# Renamed from former figS5. Shows ecological patterns in DAG structure
# and causal effects across taxonomic families.
#
# Panel (a): DAG density distribution by taxonomic family (boxplot + jitter)
# Panel (b): Diverging bar chart of significant ATE directions per family
# Panel (c): Ridge plot of ATE distributions for top variables
#
# Data required:
#   output/case4_plant/all_results_v3.csv
#   output/case4_plant/all_ate_results_v3.csv
#   outputs/Plant/Res9/CAST_ready/CAST_Species_Summary.csv
#
# Run: setwd("E:/CausalSDMs")
#      source("scripts/Plant/plot/figS2_dag_and_ate_ecology.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

fig_dir <- "figures/case4_plant/plot"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

if (!requireNamespace("ggridges", quietly = TRUE)) {
    install.packages("ggridges")
}

library(tidyverse)
library(patchwork)
library(ggridges)

# ── Theme ────────────────────────────────────────────────────────────────────
theme_pub <- function(base_size = 11) {
    theme_minimal(base_size = base_size, base_family = "sans") +
        theme(
            panel.grid.minor  = element_blank(),
            axis.title        = element_text(face = "bold"),
            plot.title        = element_text(face = "bold", hjust = 0),
            plot.subtitle     = element_text(hjust = 0, color = "grey40", size = 9)
        )
}

# ── Load data ────────────────────────────────────────────────────────────────
d <- read.csv("output/case4_plant/all_results_v3.csv", stringsAsFactors = FALSE)
ate <- read.csv("output/case4_plant/all_ate_results_v3.csv",
    stringsAsFactors = FALSE
) %>%
    mutate(coef = as.numeric(coef), significant = as.logical(significant))

sp_meta <- read.csv(
    "outputs/Plant/Res9/CAST_ready/CAST_Species_Summary.csv",
    stringsAsFactors = FALSE
) %>%
    mutate(species = gsub(" ", "_", species))

# ── Prepare family data ────────────────────────────────────────────────────
dag_data <- d %>%
    filter(!is.na(dag_density)) %>%
    distinct(species, dag_density) %>%
    left_join(sp_meta %>% select(species, family), by = "species") %>%
    mutate(family = ifelse(is.na(family), "Unknown", family))

family_order <- dag_data %>%
    group_by(family) %>%
    summarise(med = median(dag_density)) %>%
    arrange(med) %>%
    pull(family)

dag_data$family <- factor(dag_data$family, levels = family_order)

# ══════════════════════════════════════════════════════════════════════════════
# Panel (a): DAG density by family
# ══════════════════════════════════════════════════════════════════════════════
pa <- ggplot(dag_data, aes(x = family, y = dag_density, fill = family)) +
    geom_boxplot(alpha = 0.6, width = 0.5, outlier.shape = NA) +
    geom_jitter(
        width = 0.15, height = 0, alpha = 0.7, size = 1.8,
        color = "grey20"
    ) +
    scale_fill_brewer(palette = "Set3", guide = "none") +
    coord_flip() +
    labs(
        title = "(a) Causal network complexity by taxonomic family",
        subtitle = "DAG density distribution across species",
        x = "", y = "DAG Density"
    ) +
    theme_pub()

# ══════════════════════════════════════════════════════════════════════════════
# Panel (b): Diverging bar — ATE directions per family
# ══════════════════════════════════════════════════════════════════════════════
ate_fam <- ate %>%
    filter(significant == TRUE, !is.na(coef)) %>%
    left_join(sp_meta %>% select(species, family), by = "species") %>%
    mutate(
        family    = ifelse(is.na(family), "Unknown", family),
        direction = ifelse(coef > 0, "Positive", "Negative"),
        count_val = ifelse(coef > 0, 1, -1)
    )

fam_counts <- ate_fam %>%
    group_by(family, direction) %>%
    summarise(n = n(), count_val = sum(count_val), .groups = "drop")

all_combinations <- expand.grid(
    family = family_order,
    direction = c("Positive", "Negative"),
    stringsAsFactors = FALSE
)
fam_counts <- all_combinations %>%
    left_join(fam_counts, by = c("family", "direction")) %>%
    mutate(n = replace_na(n, 0), count_val = replace_na(count_val, 0))

fam_counts$family <- factor(fam_counts$family, levels = family_order)

pb <- ggplot(fam_counts, aes(x = family, y = count_val, fill = direction)) +
    geom_col(alpha = 0.9, width = 0.7) +
    scale_y_continuous(labels = abs) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.8) +
    scale_fill_manual(
        values = c("Positive" = "#27AE60", "Negative" = "#E74C3C"),
        name = "Causal\nEffect"
    ) +
    coord_flip() +
    labs(
        title = "(b) Direction of significant environmental effects",
        subtitle = "Count of positive vs negative ATEs per family",
        x = "", y = "Count of significant effects"
    ) +
    theme_pub() +
    theme(
        axis.text.y = element_blank(),
        legend.position = "bottom"
    )

# ══════════════════════════════════════════════════════════════════════════════
# Panel (c): Ridge plot — top variables ATE distributions
# ══════════════════════════════════════════════════════════════════════════════
top_vars <- ate %>%
    filter(significant == TRUE) %>%
    group_by(variable) %>%
    summarise(n = n()) %>%
    arrange(desc(n)) %>%
    slice_head(n = 7) %>%
    pull(variable)

ate_ridge <- ate %>%
    filter(significant == TRUE, variable %in% top_vars) %>%
    mutate(variable = factor(variable, levels = rev(top_vars)))

pc <- ggplot(ate_ridge, aes(x = coef, y = variable, fill = after_stat(x))) +
    geom_density_ridges_gradient(
        scale = 1.3, rel_min_height = 0.01,
        jittered_points = TRUE, point_shape = 21, point_size = 1.2,
        point_alpha = 0.7, alpha = 0.85
    ) +
    geom_vline(
        xintercept = 0, linetype = "dashed", color = "black",
        linewidth = 0.6
    ) +
    scale_fill_gradient2(
        low = "#E74C3C", mid = "white", high = "#27AE60",
        midpoint = 0, guide = "none"
    ) +
    labs(
        title = "(c) Causal effect size distributions for core variables",
        subtitle = sprintf(
            "Top %d most frequently significant variables across %d species",
            length(top_vars), length(unique(d$species))
        ),
        x = "Average Treatment Effect (ATE)", y = ""
    ) +
    theme_pub()

# ══════════════════════════════════════════════════════════════════════════════
# Combine
# ══════════════════════════════════════════════════════════════════════════════
design <- "
AABB
CCCC
"

figS2 <- pa + pb + pc +
    plot_layout(design = design, heights = c(1, 1.2)) +
    plot_annotation(
        title = "Fig S2  Causal graph complexity & ecological effect heterogeneity",
        theme = theme(
            plot.title = element_text(face = "bold", size = 13, hjust = 0.5)
        )
    )

ggsave(file.path(fig_dir, "figS2_dag_and_ate_ecology.png"),
    figS2,
    width = 13, height = 11, dpi = 300, bg = "white"
)
cat("✓ Saved figS2_dag_and_ate_ecology.png\n")
