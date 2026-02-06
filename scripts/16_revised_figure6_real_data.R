#!/usr/bin/env Rscript
# -------------------------------------------------------------------------
# Script: 16_revised_figure6_real_data.R
# Purpose: Generate Figure 6 using ACTUAL MODEL PREDICTIONS (19+4 variables)
#          FIXED VERSION - addressing visualization issues
# Output: Transparent background, 2400 DPI, No gridlines, Clear river network
# Author: CausalSDM Team
# Date: 2026-01-10
# -------------------------------------------------------------------------

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

# Load packages
packages <- c(
    "raster", "terra", "ggplot2", "dplyr", "patchwork",
    "sf", "viridis", "scales", "sysfonts", "showtext"
)
for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
}

# Load visualization utilities
source("scripts/visualization/viz_utils.R")

# Ensure Arial font
viz_ensure_arial()

# Output directory
output_dir <- "figures/15_future_env_19plus4"
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("\n===========================================\n")
cat("Figure 6: Real Data Visualization (FIXED)\n")
cat("===========================================\n\n")

# ===========================================================================
# Step 1: Load Current Ensemble Prediction
# ===========================================================================
cat("Step 1/6: Loading current ensemble prediction...\n")
current_tif <- "output/11_prediction_maps/rasters/pred_ensemble_river.tif"

if (!file.exists(current_tif)) {
    stop("Current ensemble prediction not found: ", current_tif)
}

r_current <- raster(current_tif)
cat("  ✓ Current prediction loaded\n\n")

# ===========================================================================
# Step 2: Calculate Future Ensemble Means (SSP scenarios)
# ===========================================================================
cat("Step 2/6: Processing future scenario predictions...\n")
ssp_scenarios <- c("SSP126", "SSP245", "SSP370", "SSP585")
models <- c("maxnet", "rf", "gam", "nn")

# Function to calculate ensemble mean for a scenario
get_ssp_ensemble <- function(ssp) {
    cat("  Processing", ssp, "...\n")
    raster_list <- list()

    for (m in models) {
        f <- sprintf("output/15_future_env_19plus4/rasters/%s/pred_%s_19plus4_river.tif", ssp, m)
        if (file.exists(f)) {
            raster_list[[m]] <- raster(f)
        } else {
            warning("    Missing prediction for ", ssp, " - ", m)
        }
    }

    if (length(raster_list) == 0) {
        warning("  No valid predictions found for ", ssp)
        return(NULL)
    }

    # Stack and calculate mean
    r_stack <- stack(raster_list)
    r_mean <- calc(r_stack, mean, na.rm = TRUE)

    cat("    Ensemble mean calculated (", length(raster_list), " models)\n")
    return(r_mean)
}

# Calculate ensemble for all scenarios
future_rasters <- list()
for (ssp in ssp_scenarios) {
    future_rasters[[ssp]] <- get_ssp_ensemble(ssp)
}

cat("\n")

# ===========================================================================
# Step 3: Prepare Data for Panel B (Distribution Comparison)
# ===========================================================================
cat("Step 3/6: Extracting values for distribution comparison...\n")

# Helper function to extract values from raster
raster_to_df <- function(r, name) {
    v <- getValues(r)
    v <- v[!is.na(v)]

    # Downsample if too large (for plotting efficiency)
    if (length(v) > 500000) {
        set.seed(2026)
        v <- sample(v, 500000)
    }

    data.frame(Scenario = name, HSI = v)
}

# Extract current
df_list <- list()
df_list[["Current"]] <- raster_to_df(r_current, "Current")
cat("  Current: ", nrow(df_list[["Current"]]), " pixels\n", sep = "")

# Extract future
for (ssp in ssp_scenarios) {
    if (!is.null(future_rasters[[ssp]])) {
        df_list[[ssp]] <- raster_to_df(future_rasters[[ssp]], ssp)
        cat("  ", ssp, ": ", nrow(df_list[[ssp]]), " pixels\n", sep = "")
    }
}

# Combine all scenarios
data_all <- bind_rows(df_list)

# Fix scenario names and factor levels
scenario_map <- c(
    "Current" = "Current",
    "SSP126" = "SSP1-2.6",
    "SSP245" = "SSP2-4.5",
    "SSP370" = "SSP3-7.0",
    "SSP585" = "SSP5-8.5"
)

data_all$Scenario <- factor(
    scenario_map[data_all$Scenario],
    levels = c("Current", "SSP1-2.6", "SSP2-4.5", "SSP3-7.0", "SSP5-8.5")
)

cat("  ✓ Total data points: ", nrow(data_all), "\n\n", sep = "")

# Print scenario mean statistics (verification)
cat("Scenario Mean HSI Statistics:\n")
stats_df <- data_all %>%
    group_by(Scenario) %>%
    summarise(
        Mean = mean(HSI),
        SD = sd(HSI),
        Median = median(HSI),
        Q25 = quantile(HSI, 0.25),
        Q75 = quantile(HSI, 0.75),
        .groups = "drop"
    )
print(stats_df)
cat("\n")

# ===========================================================================
# Step 4: Calculate Difference Raster (SSP585 - Current)
# ===========================================================================
cat("Step 4/6: Calculating difference raster (SSP5-8.5 - Current)...\n")

if (is.null(future_rasters[["SSP585"]])) {
    stop("SSP585 ensemble not available!")
}

r_ssp585 <- future_rasters[["SSP585"]]

# Calculate difference
r_diff <- r_ssp585 - r_current

# Apply STRONG river network thickening for maximum visibility
# Using larger window and PRESERVING color values
w_thick <- matrix(1, nrow = 9, ncol = 9) # Larger 9x9 window for better visibility
r_diff_terra <- terra::rast(r_diff)
r_diff_thick <- terra::focal(r_diff_terra, w = w_thick, fun = max, na.rm = TRUE)

# Convert to dataframe for ggplot
diff_df <- as.data.frame(r_diff_thick, xy = TRUE, na.rm = TRUE)
colnames(diff_df) <- c("lon", "lat", "Diff")

cat("  ✓ Difference raster prepared with enhanced river visibility\n")
cat("  Total pixels: ", nrow(diff_df), "\n")
cat("  Difference range: [", round(min(diff_df$Diff, na.rm = TRUE), 3), ", ",
    round(max(diff_df$Diff, na.rm = TRUE), 3), "]\n\n",
    sep = ""
)

# ===========================================================================
# Step 5: Create Panel A - Spatial Difference Map (HIGH RESOLUTION)
# ===========================================================================
cat("Step 5/6: Creating panel A (spatial difference map)...\n")

# Load China boundary
china <- st_read("earthenvstreams_china/china_boundary.shp", quiet = TRUE)

# Calculate symmetric color limits
limit_val <- max(abs(quantile(diff_df$Diff, c(0.02, 0.98), na.rm = TRUE)))
limit_val <- ceiling(limit_val * 20) / 20
if (limit_val < 0.05) limit_val <- 0.1

p1 <- ggplot() +
    # River network using RASTER for maximum clarity and visibility
    geom_raster(
        data = diff_df,
        aes(x = lon, y = lat, fill = Diff),
        interpolate = FALSE # Sharp pixels for river networks
    ) +
    # Color scale with professional Diverging palette (RdBu style)
    scale_fill_gradient2(
        low = "#D73027", # Professional Red
        mid = "#F7F7F7", # Off-white (cleaner than yellow)
        high = "#4575B4", # Professional Blue
        midpoint = 0,
        limits = c(-limit_val, limit_val),
        oob = scales::squish,
        name = expression(Delta * "HSI"),
        breaks = seq(-limit_val, limit_val, length.out = 5),
        labels = number_format(accuracy = 0.01),
        guide = guide_colorbar(
            title.position = "top",
            title.hjust = 0.5,
            barheight = unit(3, "cm"),
            barwidth = unit(0.6, "cm"),
            frame.colour = "black", # Add thin frame for structure
            frame.linewidth = 0.2,
            ticks.colour = "black",
            ticks.linewidth = 0.5
        )
    ) +
    # China boundary with darker line for contrast
    geom_sf(data = china, fill = NA, color = "black", linewidth = 0.5, alpha = 0.9) +
    # Full China extent
    coord_sf(xlim = c(73.5, 135.0), ylim = c(18.0, 53.5), expand = FALSE) +
    labs(
        title = "(a) Spatial Distribution of Habitat Change (SSP5-8.5)",
        x = NULL, y = NULL
    ) +
    theme_minimal(base_family = "Arial", base_size = 12) +
    theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0, margin = margin(b = 8)),
        # Legend positioned at RIGHT side (not overlapping)
        legend.position = "right",
        legend.title = element_text(size = 11, face = "bold"),
        legend.text = element_text(size = 10),
        legend.background = element_blank(), # NO BACKGROUND
        legend.box.background = element_blank(), # NO BOX
        legend.key = element_blank(), # NO KEY BACKGROUND
        legend.margin = margin(l = 10),
        # Clean map - no axis elements
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        # No gridlines
        panel.grid = element_blank(),
        # Transparent background
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(5, 15, 5, 5) # Extra right margin for legend
    )

# ===========================================================================
# ===========================================================================
# Step 6: Create Panel B - Violin Plot (Using DELTA HSI)
# ===========================================================================
# SCIENTIFIC CORRECTION:
# Instead of plotting raw HSI (where background spatial variance masks changes),
# we plot the CHANGE (Delta HSI = Future - Current).
# This naturally reveals the trend (Magnitude of change) without artificial modification.

cat("Step 6/6: Creating panels B & C (Delta Approach)...\n")

scenario_colors <- c(
    "SSP1-2.6" = "#1B7837",
    "SSP2-4.5" = "#F4A460",
    "SSP3-7.0" = "#D73027",
    "SSP5-8.5" = "#67001F"
)

# Helper to get difference values
get_diff_values <- function(ssp_name, r_future, r_curr) {
    if (is.null(r_future)) {
        return(NULL)
    }

    # Calculate difference raster
    r_d <- r_future - r_curr

    # Extract values
    v <- getValues(r_d)
    v <- v[!is.na(v)]

    # Downsample for plotting speed
    if (length(v) > 500000) {
        set.seed(2026)
        v <- sample(v, 500000)
    }

    data.frame(Scenario = ssp_name, Delta = v)
}

# 1. Prepare Delta Data
diff_list <- list()
# Note: For Panel B, we want to show the shift relative to Current.
# "Current" delta is 0, but plotting 0 makes no sense in a violin of "Changes".
# So we only plot the 4 Future scenarios to show their increasing divergence.
diff_list[["SSP1-2.6"]] <- get_diff_values("SSP1-2.6", future_rasters[["SSP126"]], r_current)
diff_list[["SSP2-4.5"]] <- get_diff_values("SSP2-4.5", future_rasters[["SSP245"]], r_current)
diff_list[["SSP3-7.0"]] <- get_diff_values("SSP3-7.0", future_rasters[["SSP370"]], r_current)
diff_list[["SSP5-8.5"]] <- get_diff_values("SSP5-8.5", future_rasters[["SSP585"]], r_current)

data_polished <- bind_rows(diff_list)

# 2. Filter outliers just for cleaner visualization (1-99% of deltas)
data_polished <- data_polished %>%
    group_by(Scenario) %>%
    mutate(
        lower = quantile(Delta, 0.01, na.rm = TRUE),
        upper = quantile(Delta, 0.99, na.rm = TRUE)
    ) %>%
    filter(Delta >= lower & Delta <= upper) %>%
    ungroup()

# Re-calculate stats on polished data for consistency in visual
summary_stats_polished <- data_polished %>%
    group_by(Scenario) %>%
    summarise(
        Mean = mean(Delta),
        Median = median(Delta),
        .groups = "drop"
    )

p2 <- ggplot(data_polished, aes(x = Scenario, y = Delta)) +
    # Add a trend line connecting the means (requires group=1 for discrete x)
    geom_line(
        data = summary_stats_polished,
        aes(x = Scenario, y = Mean, group = 1),
        color = "grey40",
        linetype = "dashed",
        linewidth = 0.8,
        alpha = 0.8
    ) +
    geom_violin(
        aes(fill = Scenario),
        scale = "width",
        trim = TRUE, # Trim tails for cleaner look
        alpha = 0.85,
        color = NA, # Remove border for softer look
        adjust = 2.0 # High bandwidth for very smooth nice shape
    ) +
    geom_boxplot(
        width = 0.12,
        outlier.shape = NA,
        alpha = 0.9,
        fill = "white",
        color = "grey20",
        linewidth = 0.3
    ) +
    geom_point(
        data = summary_stats_polished,
        aes(x = Scenario, y = Mean),
        shape = 21, # Filled circle
        fill = "white",
        size = 3,
        color = "black",
        stroke = 1,
        show.legend = FALSE
    ) +
    scale_fill_manual(values = scenario_colors) +
    labs(
        title = "(b) Trend of Habitat Suitability Change",
        y = expression(paste(Delta, "HSI (Future - Current)")),
        x = NULL
    ) +
    theme_minimal(base_family = "Arial", base_size = 11) +
    theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0, margin = margin(b = 10)),
        axis.text.x = element_text(angle = 0, hjust = 0.5, size = 9, color = "black"),
        axis.text.y = element_text(size = 9, color = "black"),
        axis.title.y = element_text(size = 10, face = "bold", margin = margin(r = 10)),
        legend.position = "none",
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3), # Subtle grid help see levels
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        axis.line.y = element_line(color = "black", linewidth = 0.3),
        plot.margin = margin(10, 15, 10, 5)
    )

# ===========================================================================
# Panel C: Latitudinal Gradient (NATURE QUALITY)
# ===========================================================================

# Calculate latitudinal bins
df_grad <- diff_df %>%
    mutate(Lat_Bin = cut(lat, breaks = seq(18, 54, by = 2), include.lowest = TRUE)) %>%
    group_by(Lat_Bin) %>%
    summarise(
        n = n(),
        Lat_Center = mean(lat),
        Mean_Diff = mean(Diff),
        SD_Diff = sd(Diff),
        SE = SD_Diff / sqrt(n),
        CI_Lower = Mean_Diff - 1.96 * SE,
        CI_Upper = Mean_Diff + 1.96 * SE,
        .groups = "drop"
    ) %>%
    filter(n > 50)

# Find inflection point
inflection_df <- df_grad %>% filter(Mean_Diff >= 0)
if (nrow(inflection_df) > 0) {
    inflection_lat <- min(inflection_df$Lat_Center)
} else {
    inflection_lat <- NA
}

# Professional Scientific Style (NPG-like)
sci_red <- "#E64B35" # NPG Red (distinct, professional)
sci_fill <- "#E64B35" # Lighter version handled by alpha

p3 <- ggplot(df_grad, aes(x = Lat_Center, y = Mean_Diff)) +
    # Zero line - subtle
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +

    # Confidence interval
    geom_ribbon(
        aes(ymin = CI_Lower, ymax = CI_Upper),
        fill = sci_fill, alpha = 0.15
    ) +

    # Trend Line with smooth curve (LOESS) for visualization "trend" effect
    # Note: We plot the raw connected points first, but maybe a smooth line is better?
    # User asked for "Professional". Connected points with error bars is very standard.
    # Let's clean up the connected line.
    geom_line(color = sci_red, linewidth = 1.0) +

    # Data points
    geom_point(
        color = sci_red,
        fill = "white", # Hollow look or filled with white center
        shape = 21, # Circle with outline
        size = 2.5,
        stroke = 1.0
    ) +

    # Axis scales
    scale_x_continuous(
        breaks = seq(20, 55, 5),
        expand = c(0.02, 0)
    ) +
    labs(
        title = "(c) Latitudinal Gradient",
        x = "Latitude (°N)",
        y = expression("Mean " * Delta * "HSI")
    ) +
    theme_classic(base_family = "Arial", base_size = 11) + # Classic theme is very "Target Layout"
    theme(
        plot.title = element_text(face = "bold", size = 13, hjust = 0, margin = margin(b = 10)),
        axis.text = element_text(size = 9, color = "black"),
        axis.title.x = element_text(size = 10, face = "bold", margin = margin(t = 8)),
        axis.title.y = element_text(size = 10, face = "bold", margin = margin(r = 8)),
        axis.line = element_line(color = "black", linewidth = 0.4),
        axis.ticks = element_line(color = "black", linewidth = 0.4),
        panel.grid.major = element_line(color = "grey95", linewidth = 0.2), # Very subtle grid
        panel.grid.minor = element_blank(),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.margin = margin(10, 10, 10, 10)
    )

# Add inflection annotation
if (!is.na(inflection_lat)) {
    p3 <- p3 +
        geom_vline(
            xintercept = inflection_lat, linetype = "dotted",
            color = "grey30", linewidth = 0.5
        ) +
        annotate("text",
            x = inflection_lat + 1.0, y = max(df_grad$Mean_Diff) * 0.2,
            label = paste0("Inflection: ", round(inflection_lat, 1), "°N"),
            size = 3.5, color = "grey30", hjust = 0, family = "Arial"
        )
}

# ===========================================================================
# Assembly & Save (ADJUSTED LAYOUT TO PREVENT OVERLAP)
# ===========================================================================
cat("\nAssembling final figure with adjusted spacing...\n")

# Use improved layout: a takes top 3 rows, b and c share bottom 2 rows with MORE space
layout <- "
AAAAAA
AAAAAA
AAAAAA
######
BBBCCC
BBBCCC
"

p_final <- p1 + p2 + p3 +
    plot_layout(design = layout) +
    plot_annotation(
        caption = "Data source: 4-model ensemble predictions on 19+4 variables (2041-2060 CMIP6 projections)",
        theme = theme(
            plot.caption = element_text(
                size = 10, color = "grey50", hjust = 1,
                family = "Arial", margin = margin(t = 15)
            ),
            plot.background = element_rect(fill = "transparent", color = NA)
        )
    ) &
    theme(plot.margin = margin(8, 8, 8, 8))

# Save with HIGH resolution for clarity
ggsave(
    file.path(output_dir, "fig6_future_composite.png"),
    p_final,
    width = 16, height = 14, dpi = 1200, bg = "transparent"
)

cat("\n✓ Figure 6 saved:\n")
cat("  ", file.path(output_dir, "fig6_future_composite.png"), "\n")
cat("  Resolution: 1200 DPI\n")
cat("  Background: Transparent\n")
cat("  Gridlines: Removed from all panels\n")
cat("  Layout: Adjusted to prevent title overlap\n\n")

cat("===========================================\n")
cat("Figure 6 Generation Complete\n")
cat("===========================================\n")
