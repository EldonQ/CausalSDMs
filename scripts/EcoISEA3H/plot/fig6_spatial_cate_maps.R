################################################################################
# Fig 6: Spatial CATE Maps — Real Causal Forest Heterogeneous Effects
#
# Scientific Question Q4:
#   Can CAST provide ecological interpretability that traditional SDMs cannot?
#   Specifically, can it reveal WHERE each environmental driver exerts its
#   strongest causal influence?
#
# CRITICAL UPDATE: Now uses REAL grf::causal_forest CATE predictions saved
#   in all_spatial_cate_v3.csv, NOT the linear proxy (scale(var) * ATE).
#
# Panel (a)(b): Spatial CATE heatmaps for top-2 significant causal drivers
# Panel (c):    ATE forest plot — global causal effects with 95% CI
# Panel (d):    Causal role bar chart — DAG out-degree
#
# Species selection: Automatically selects the species with the most
#   significant ATE variables from the EcoISEA3H dataset
#
# Data required:
#   output/case2_eco/all_spatial_cate_v3.csv  (real Causal Forest CATE!)
#   output/case2_eco/all_ate_results_v3.csv
#   output/case2_eco/all_results_v3.csv
#   output/case2_eco/all_screening_v3.csv
#   outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv
#
# Run: setwd("E:/CausalSDMs")
#      source("scripts/EcoISEA3H/plot/fig6_spatial_cate_maps.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

fig_dir <- "figures/case2_eco/plot"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(patchwork)

# ── Themes ───────────────────────────────────────────────────────────────────
theme_map <- function(base_size = 11) {
    theme_minimal(base_size = base_size, base_family = "sans") +
        theme(
            panel.grid        = element_blank(),
            axis.text         = element_text(size = 8, color = "grey50"),
            axis.title        = element_text(face = "bold", size = 9),
            plot.title        = element_text(face = "bold", hjust = 0, size = 11),
            plot.subtitle     = element_text(hjust = 0, color = "grey40", size = 8.5),
            legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
            legend.title      = element_text(size = 8, face = "bold"),
            legend.text       = element_text(size = 7),
            legend.key.height = unit(0.8, "cm")
        )
}

theme_pub <- function(base_size = 11) {
    theme_minimal(base_size = base_size, base_family = "sans") +
        theme(
            panel.grid.minor  = element_blank(),
            axis.title        = element_text(face = "bold"),
            plot.title        = element_text(face = "bold", hjust = 0, size = 11),
            plot.subtitle     = element_text(hjust = 0, color = "grey40", size = 8.5),
            legend.background = element_rect(fill = "white", color = NA)
        )
}

# ── Load data ────────────────────────────────────────────────────────────────
cate_all <- read.csv("output/case2_eco/all_spatial_cate_v3.csv",
    stringsAsFactors = FALSE
)

ate <- read.csv("output/case2_eco/all_ate_results_v3.csv",
    stringsAsFactors = FALSE
) %>%
    mutate(
        coef        = as.numeric(coef),
        se          = as.numeric(se),
        significant = as.logical(significant),
        ci_lower    = coef - 1.96 * se,
        ci_upper    = coef + 1.96 * se
    )

results <- read.csv("output/case2_eco/all_results_v3.csv",
    stringsAsFactors = FALSE
)

screen <- read.csv("output/case2_eco/all_screening_v3.csv",
    stringsAsFactors = FALSE
)

sp_meta <- read.csv(
    "outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv",
    stringsAsFactors = FALSE
) %>%
    mutate(species = gsub(" ", "_", species))

# ── AUTO-SELECT best species for CATE showcase ──────────────────────────────
# Choose the species with the most CATE variables computed (visual richness)
# AND the most significant ATE variables
cate_species_vars <- cate_all %>%
    group_by(species) %>%
    summarise(
        n_cate_vars = n_distinct(variable),
        n_points = n(),
        .groups = "drop"
    )

species_sig_count <- ate %>%
    filter(significant == TRUE) %>%
    group_by(species) %>%
    summarise(n_sig = n(), .groups = "drop")

best_sp <- cate_species_vars %>%
    left_join(species_sig_count, by = "species") %>%
    mutate(n_sig = replace_na(n_sig, 0)) %>%
    arrange(desc(n_sig), desc(n_cate_vars)) %>%
    slice_head(n = 1) %>%
    pull(species)

target_sp <- best_sp
target_region <- "China_Res9"

sp_display <- gsub("_", " ", target_sp)
sp_family <- sp_meta %>%
    filter(species == target_sp) %>%
    pull(family) %>%
    first()
if (is.null(sp_family) || is.na(sp_family)) sp_family <- ""

n_sig_ate <- species_sig_count %>%
    filter(species == target_sp) %>%
    pull(n_sig)
if (length(n_sig_ate) == 0) n_sig_ate <- 0

cat(sprintf(
    "[Fig 6] Auto-selected species: %s (%s) with %d significant ATE variables\n",
    sp_display, sp_family, n_sig_ate
))

# ── Species CATE data (real Causal Forest predictions!) ─────────────────────
sp_cate <- cate_all %>% filter(species == target_sp)
sp_cate_vars <- unique(sp_cate$variable)

cat(sprintf(
    "[Fig 6] Real CATE data: %d grid points × %d variables (%s)\n",
    nrow(sp_cate) / length(sp_cate_vars),
    length(sp_cate_vars),
    paste(sp_cate_vars, collapse = ", ")
))

# ── Species ATE data ─────────────────────────────────────────────────────────
sp_ate <- ate %>%
    filter(species == target_sp) %>%
    arrange(desc(abs(coef)))

sp_ate_sig <- sp_ate %>% filter(significant == TRUE)

# ══════════════════════════════════════════════════════════════════════════════
# Panel (c): ATE forest plot — global causal effect context
# ══════════════════════════════════════════════════════════════════════════════
ate_plot <- sp_ate %>%
    mutate(
        sig_marker = ifelse(significant,
            "Significant (p < 0.05)",
            "Not significant"
        ),
        variable = factor(variable, levels = rev(variable)) # keep sorted by |coef|
    )

pc <- ggplot(
    ate_plot,
    aes(
        x = coef,
        y = variable,
        color = sig_marker
    )
) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(size = 2.5) +
    geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper),
        width = 0.25, linewidth = 0.5
    ) +
    scale_color_manual(
        values = c(
            "Significant (p < 0.05)" = "#C0392B",
            "Not significant" = "grey55"
        ),
        name = ""
    ) +
    labs(
        title = sprintf("(c) ATE forest plot"),
        subtitle = "DML cross-fitting | 95% CI | Red = significant causal driver",
        x = "Average Treatment Effect (ATE)", y = ""
    ) +
    theme_pub() +
    theme(
        legend.position = "bottom",
        panel.grid.major.y = element_line(color = "grey93")
    )

# ══════════════════════════════════════════════════════════════════════════════
# Panel (d): Causal roles bar chart
# ══════════════════════════════════════════════════════════════════════════════
sp_screen <- screen %>%
    filter(species == target_sp)

cast_res <- results %>%
    filter(species == target_sp, model == "CAST")
cast_vars <- sp_screen %>%
    arrange(desc(score_total)) %>%
    slice_head(n = ifelse(nrow(cast_res) > 0, cast_res$n_vars[1], 6)) %>%
    pull(variable)

role_data <- sp_screen %>%
    filter(variable %in% cast_vars) %>%
    mutate(
        role = case_when(
            out_degree >= 2 ~ "Root",
            out_degree == 1 ~ "Mediator",
            TRUE ~ "Terminal"
        )
    ) %>%
    select(variable, out_degree, role, score_total)

role_colors <- c("Root" = "#2C3E50", "Mediator" = "#E67E22", "Terminal" = "#27AE60")

pd <- ggplot(role_data, aes(
    x = reorder(variable, -out_degree), y = out_degree,
    fill = role
)) +
    geom_col(alpha = 0.85, width = 0.65) +
    scale_fill_manual(values = role_colors, name = "Causal Role") +
    labs(
        title = "(d) Causal roles",
        subtitle = "Out-degree = number of downstream variables influenced",
        x = "", y = "DAG out-degree"
    ) +
    theme_pub() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        legend.position = "bottom"
    )

# ══════════════════════════════════════════════════════════════════════════════
# Panels (a)(b): Spatial CATE heatmaps — REAL Causal Forest predictions
# ══════════════════════════════════════════════════════════════════════════════
# Select top-2 variables for CATE maps (prioritize significant ATEs that
# also have CATE data)
cate_candidates <- sp_ate %>%
    filter(variable %in% sp_cate_vars) %>%
    arrange(desc(significant), desc(abs(coef)))

if (nrow(cate_candidates) >= 2) {
    cate_top2 <- cate_candidates %>% slice_head(n = 2)
} else {
    cate_top2 <- sp_ate %>%
        filter(variable %in% sp_cate_vars) %>%
        slice_head(n = min(2, n()))
}

make_real_cate_panel <- function(sp_cate_df, var_name, ate_coef, panel_label) {
    plot_df <- sp_cate_df %>%
        filter(variable == var_name) %>%
        drop_na()

    if (nrow(plot_df) < 10) {
        return(ggplot() +
            theme_void() +
            labs(title = paste(panel_label, "— Insufficient CATE data")))
    }

    # Symmetric color scale for diverging palette
    cate_abs_max <- quantile(abs(plot_df$cate), 0.98, na.rm = TRUE)

    p <- ggplot(plot_df, aes(x = lon, y = lat, color = cate)) +
        geom_point(size = 0.5, alpha = 0.80, shape = 15) +
        coord_fixed() +
        scale_color_gradient2(
            low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
            midpoint = 0,
            limits = c(-cate_abs_max, cate_abs_max),
            oob = scales::squish,
            name = "CATE"
        ) +
        labs(
            title = sprintf("%s  Causal effect — %s", panel_label, var_name),
            subtitle = sprintf(
                "grf::causal_forest  |  ATE = %.4f  |  n = %s grid cells",
                ate_coef, format(nrow(plot_df), big.mark = ",")
            ),
            x = "Longitude", y = "Latitude"
        ) +
        theme_map()
    p
}

if (nrow(cate_top2) >= 2) {
    pa <- make_real_cate_panel(sp_cate, cate_top2$variable[1], cate_top2$coef[1], "(a)")
    pb <- make_real_cate_panel(sp_cate, cate_top2$variable[2], cate_top2$coef[2], "(b)")
} else if (nrow(cate_top2) == 1) {
    pa <- make_real_cate_panel(sp_cate, cate_top2$variable[1], cate_top2$coef[1], "(a)")
    pb <- ggplot() +
        theme_void() +
        labs(title = "(b) Only one CATE variable available")
} else {
    pa <- ggplot() +
        theme_void() +
        labs(title = "(a) No CATE data")
    pb <- ggplot() +
        theme_void() +
        labs(title = "(b) No CATE data")
}

# ══════════════════════════════════════════════════════════════════════════════
# Combine  2×2 layout
#   (a) CATE var1   |  (b) CATE var2
#   (c) ATE forest  |  (d) Causal roles
# ══════════════════════════════════════════════════════════════════════════════
fig6 <- (pa | pb) / (pc | pd) +
    plot_layout(heights = c(1.1, 1)) +
    plot_annotation(
        title = sprintf(
            "Fig 6  Spatially explicit causal effects: %s (%s)",
            sp_display, sp_family
        ),
        subtitle = paste0(
            "Real Causal Forest (grf) CATE predictions — heterogeneous ",
            "causal effects estimated at each grid cell. ",
            "Warm = strong positive effect; cool = negative/weak."
        ),
        theme = theme(
            plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
            plot.subtitle = element_text(
                face = "italic", size = 9.5, hjust = 0.5,
                color = "grey40"
            )
        )
    )

ggsave(file.path(fig_dir, "fig6_spatial_cate_maps.png"),
    fig6,
    width = 14, height = 11, dpi = 300, bg = "white"
)
cat(sprintf("✓ Saved fig6_spatial_cate_maps.png (species: %s)\n", sp_display))
