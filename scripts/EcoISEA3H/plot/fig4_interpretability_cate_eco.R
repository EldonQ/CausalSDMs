################################################################################
# Fig 4 (Eco): Ecological Interpretability + CATE Spatial Map (4-panel)
# (a) RF importance rank vs |ATE| rank slope chart (China_Res9)
# (b) Spearman rho distribution across 27 species
# (c)(d) CATE spatial heatmap (Ovis_ammon, 2 key variables) over Chinese grid
#
# Requires: all_screening_v3.csv, all_ate_results_v3.csv, China_EnvData_Res9_Screened.csv
# Run: setwd("E:/CausalSDMs"); source("scripts/case2_eco/plot/fig4_interpretability_cate_eco.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

fig_dir <- "figures/case2_eco/plot"
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

# ---- Load data ----
screen <- read.csv("output/case2_eco/all_screening_v3.csv", stringsAsFactors = FALSE)
ate <- read.csv("output/case2_eco/all_ate_results_v3.csv", stringsAsFactors = FALSE)

# ══════════════════════════════════════════════════════════════
# (a) Slope chart: RF importance rank vs |ATE| rank (China_Res9)
# ══════════════════════════════════════════════════════════════
showcase_region <- "China_Res9"

rank_data <- screen %>%
    left_join(
        ate %>% select(region, species, variable, ate_coef = coef),
        by = c("region", "species", "variable")
    ) %>%
    mutate(ate_coef = as.numeric(ate_coef)) %>%
    group_by(species, region) %>%
    mutate(
        rank_rf  = rank(-importance, ties.method = "average"),
        rank_ate = rank(-abs(ate_coef), ties.method = "average", na.last = TRUE)
    ) %>%
    ungroup()

# Aggregate across species: mean rank per variable
slope_dat <- rank_data %>%
    group_by(variable) %>%
    summarise(
        mean_rank_rf = mean(rank_rf, na.rm = TRUE),
        mean_rank_ate = mean(rank_ate, na.rm = TRUE),
        n_sp = n(),
        .groups = "drop"
    ) %>%
    filter(n_sp >= 5) %>% # only variables appearing in enough species
    mutate(
        y_left  = rank(mean_rank_rf),
        y_right = rank(mean_rank_ate)
    )

pa <- ggplot(slope_dat) +
    geom_segment(aes(x = 1, xend = 2, y = y_left, yend = y_right, color = variable),
        linewidth = 1.1, alpha = 0.8
    ) +
    geom_point(aes(x = 1, y = y_left), size = 2.5, color = "#27AE60") +
    geom_point(aes(x = 2, y = y_right), size = 2.5, color = "#C0392B") +
    geom_text(aes(x = 0.92, y = y_left, label = variable), hjust = 1, size = 3, fontface = "bold") +
    geom_text(aes(x = 2.08, y = y_right, label = variable), hjust = 0, size = 3, fontface = "bold") +
    scale_x_continuous(
        limits = c(0.4, 2.6), breaks = c(1, 2),
        labels = c("RF importance\nrank", "|ATE|\nrank")
    ) +
    scale_color_brewer(palette = "Set3", guide = "none") +
    scale_y_reverse() +
    labs(
        title = sprintf("(a) Correlation vs causal rank [%s]", showcase_region),
        subtitle = "Line crossings = ranking discrepancy",
        x = "", y = "Mean rank across species"
    ) +
    theme_pub()

# ══════════════════════════════════════════════════════════════
# (b) Spearman rho distribution (RF rank vs ATE rank, per species)
# ══════════════════════════════════════════════════════════════
spearman_per_sp <- rank_data %>%
    filter(!is.na(rank_rf), !is.na(rank_ate)) %>%
    group_by(species) %>%
    summarise(
        rho = cor(rank_rf, rank_ate, method = "spearman", use = "pairwise"),
        n_vars = n(),
        .groups = "drop"
    ) %>%
    filter(!is.na(rho))

mean_rho <- mean(spearman_per_sp$rho, na.rm = TRUE)

pb <- ggplot(spearman_per_sp, aes(x = rho)) +
    geom_histogram(binwidth = 0.1, fill = "#2980B9", alpha = 0.7, color = "white") +
    geom_vline(xintercept = mean_rho, linetype = "dashed", color = "#C0392B", linewidth = 0.8) +
    annotate("text",
        x = mean_rho + 0.05, y = Inf, vjust = 2,
        label = sprintf("Mean ρ = %.2f", mean_rho),
        fontface = "bold", size = 3.5, color = "#C0392B"
    ) +
    scale_x_continuous(limits = c(-1, 1)) +
    labs(
        title = sprintf("(b) RF vs ATE rank agreement [%s]", showcase_region),
        subtitle = "Spearman ρ per species | Low ρ = causal ≠ correlational",
        x = "Spearman ρ (RF rank vs |ATE| rank)", y = "Number of species"
    ) +
    theme_pub()

# ══════════════════════════════════════════════════════════════
# (c)(d) CATE spatial heatmap: Ovis_ammon (Argali Sheep)
# ══════════════════════════════════════════════════════════════
target_sp <- "Ovis_ammon" # Classic widespread ungulate

# Load Environmental Grid for China (Screened variables)
env_grid_path <- "outputs/EcoISEA3H/Res9/CAST_ready/China_EnvData_Res9_Screened.csv"
if (file.exists(env_grid_path)) {
    env_grid <- read.csv(env_grid_path, stringsAsFactors = FALSE)

    # Load ATE results for this species to identify top causal variables
    sp_ate <- ate %>%
        mutate(coef = as.numeric(coef), significant = as.logical(significant)) %>%
        filter(species == target_sp, significant == TRUE) %>%
        arrange(desc(abs(coef)))

    if (nrow(sp_ate) >= 2) {
        # Pick top 2 significant causal variables
        cate_var1 <- sp_ate$variable[1]
        cate_var2 <- sp_ate$variable[2]

        make_cate_panel <- function(env_grid, var_name, ate_coef, panel_label) {
            if (!var_name %in% names(env_grid)) {
                return(ggplot() +
                    theme_void() +
                    labs(title = paste(panel_label, "Variable not in grid")))
            }
            plot_df <- env_grid %>%
                select(lon, lat, all_of(var_name)) %>%
                drop_na() %>%
                mutate(
                    var_scaled = as.numeric(scale(!!sym(var_name))[, 1]),
                    cate_proxy = var_scaled * ate_coef
                )

            ggplot(plot_df, aes(x = lon, y = lat, color = cate_proxy)) +
                geom_point(size = 0.5, alpha = 0.8) +
                scale_color_gradient2(
                    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                    midpoint = 0, name = "CATE\nproxy"
                ) +
                labs(
                    title = sprintf("%s ATE-weighted effect of %s", panel_label, var_name),
                    subtitle = sprintf("ATE = %.3f | %s (China)", ate_coef, gsub("_", " ", target_sp)),
                    x = "Longitude", y = "Latitude"
                ) +
                theme_pub() +
                theme(legend.position = "right", panel.grid = element_blank())
        }

        pc <- make_cate_panel(env_grid, cate_var1, sp_ate$coef[1], "(c)")
        pd <- make_cate_panel(env_grid, cate_var2, sp_ate$coef[2], "(d)")
    } else {
        pc <- ggplot() +
            theme_void() +
            labs(title = "(c) Insufficient significant ATE vars")
        pd <- ggplot() +
            theme_void() +
            labs(title = "(d) Insufficient significant ATE vars")
    }
} else {
    pc <- ggplot() +
        theme_void() +
        labs(title = "(c) Extracted Grid data missing")
    pd <- ggplot() +
        theme_void() +
        labs(title = "(d) Extracted Grid data missing")
}

# ══════════════════════════════════════════════════════════════
# Combine
# ══════════════════════════════════════════════════════════════
fig4 <- (pa | pb) / (pc | pd) +
    plot_annotation(
        title = "Fig 4 (Eco)  Ecological Interpretability & Causal Effect Heterogeneity",
        theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    )

ggsave(file.path(fig_dir, "fig4_interpretability_cate_eco.png"),
    fig4,
    width = 14, height = 11, dpi = 300, bg = "white"
)
cat("Saved fig4_interpretability_cate_eco.png\n")
