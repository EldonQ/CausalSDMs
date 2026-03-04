# ==============================================================================
# Fig 5: 真实因果森林空间异质性映射 (Results 3.5)
# Narrative: The actual Spatial Heterogeneity of Causal Effects computed
#            via Causal Forest (grf), showing localized, non-linear CATE.
# ==============================================================================

rm(list = ls())
setwd("E:/CausalSDMs")

out_dir <- "figures/CASTplot"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(grf)
library(patchwork)
library(sf)
library(terra)
library(tidyterra)
library(gstat)

# Theme for high-end transparent floating map
theme_map <- function() {
    theme_void(base_family = "sans") +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5, size = 15, margin = margin(b = 10)),
            plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 11, margin = margin(b = 20)),
            legend.position = "right",
            legend.title = element_text(face = "bold", size = 11),
            plot.background = element_rect(fill = "transparent", color = NA),
            panel.background = element_rect(fill = "transparent", color = NA),
            legend.background = element_rect(fill = "transparent", color = NA),
            legend.box.background = element_rect(fill = "transparent", color = NA)
        )
}

# ------------------------------------------------------------------------------
# 1. Configuration & Data Load
# ------------------------------------------------------------------------------
target_sp <- "Ovis_ammon" # Focus species with highly heterogeneous habitats

# Model outputs to extract top causal variables
ate_results <- read.csv("output/case2_eco/all_ate_results_v3.csv", stringsAsFactors = FALSE)

# Species actual training data for Causal Forest
sp_data_file <- paste0("outputs/EcoISEA3H/Res9/CAST_ready/species_data_screened/CAST_", target_sp, "_Res9_screened.csv")
# Prediction background grid
env_grid_path <- "outputs/EcoISEA3H/Res9/CAST_ready/China_EnvData_Res9_Screened.csv"

# Load China shapefile
china_shp_path <- "chinashp/china.shp"
if (file.exists(china_shp_path)) {
    china_sf <- st_read(china_shp_path, quiet = TRUE)
} else {
    stop("Cannot find China shapefile at ", china_shp_path)
}

if (file.exists(sp_data_file) && file.exists(env_grid_path)) {
    # --------------------------------------------------------------------------
    # 2. Get the Top 2 Causal Variables from ATE
    # --------------------------------------------------------------------------
    sp_ate <- ate_results %>%
        filter(species == target_sp, as.logical(significant) == TRUE) %>%
        arrange(desc(abs(as.numeric(coef))))

    if (nrow(sp_ate) >= 2) {
        var1 <- sp_ate$variable[1]
        var2 <- sp_ate$variable[2]

        # Load empirical species data
        sp_df <- read.csv(sp_data_file, stringsAsFactors = FALSE)
        meta_cols <- c("HID", "lon", "lat", "species", "sid", "family", "category", "presence", "fraction")
        env_cols <- setdiff(names(sp_df), meta_cols)

        # Prepare inputs for GRF
        X_train <- sp_df[, env_cols, drop = FALSE]
        Y_train <- sp_df$presence

        # Load map grid
        env_grid <- read.csv(env_grid_path, stringsAsFactors = FALSE)

        # Match grid feature names
        common_cols <- intersect(names(env_grid), env_cols)
        X_map <- env_grid[, common_cols, drop = FALSE]

        # We need to impute or drop NAs safely for GRF predictability
        env_grid$is_complete <- complete.cases(X_map)
        X_map_valid <- X_map[env_grid$is_complete, ]

        make_true_cate_panel <- function(treatment_var, panel_label) {
            cat("Training Real Causal Forest for:", treatment_var, "...\n")

            # W is our treatment. X are the confounders.
            W_train <- X_train[[treatment_var]]
            X_cf_train <- X_train[, setdiff(env_cols, treatment_var), drop = FALSE]

            # Train causal forest
            set.seed(42)
            c_forest <- causal_forest(
                X = as.matrix(X_cf_train),
                Y = Y_train,
                W = W_train,
                num.trees = 500
            )

            # Predict CATE on the whole spatial grid
            cat("Predicting spatial CATE for:", treatment_var, "...\n")
            X_cf_map <- as.matrix(X_map_valid[, setdiff(env_cols, treatment_var), drop = FALSE])
            cate_preds <- predict(c_forest, X_cf_map)$predictions

            # ------------------------------------------------------------------
            # Inverse Distance Weighting (IDW) Interpolation for Smooth Heatmap
            # ------------------------------------------------------------------
            plot_df <- env_grid[env_grid$is_complete, c("lon", "lat")]
            plot_df$true_cate <- as.numeric(cate_preds)

            # Convert to spatial points
            v_sf <- st_as_sf(plot_df, coords = c("lon", "lat"), crs = 4326)

            # Create a fine blank grid dataframe for IDW
            grid_xs <- seq(min(plot_df$lon) - 0.5, max(plot_df$lon) + 0.5, by = 0.1)
            grid_ys <- seq(min(plot_df$lat) - 0.5, max(plot_df$lat) + 0.5, by = 0.1)
            grid_df <- expand.grid(lon = grid_xs, lat = grid_ys)
            grid_sf <- st_as_sf(grid_df, coords = c("lon", "lat"), crs = 4326)

            # Interpolate globally using IDW
            cat("Interpolating spatial heatmap via IDW...\n")
            idw_res <- gstat::idw(true_cate ~ 1, locations = v_sf, newdata = grid_sf, idp = 2, nmax = 12, debug.level = 0)

            # Put back into dataframe
            grid_df$CATE <- idw_res$var1.pred

            # Convert to raster to apply exact China mask
            r_interp <- rast(grid_df, type = "xyz", crs = "EPSG:4326")
            r_masked <- mask(r_interp, china_sf)

            # Extract valid masked values back into plotting dataframe
            masked_df <- as.data.frame(r_masked, xy = TRUE, na.rm = TRUE)
            names(masked_df)[3] <- "CATE"

            # Build plot
            ggplot() +
                geom_tile(data = masked_df, aes(x = x, y = y, fill = CATE)) +
                geom_sf(data = china_sf, fill = NA, color = "black", linewidth = 0.3) +
                scale_fill_gradient2(
                    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                    midpoint = mean(plot_df$true_cate, na.rm = T), name = "CATE\n(grf)",
                    na.value = "transparent"
                ) +
                coord_sf(expand = FALSE) + # Prevent spatial squishing
                labs(
                    title = sprintf("%s Spatial CATE Map: %s", panel_label, treatment_var),
                    subtitle = sprintf("Target Species: %s (Geographic deviation from mean effect)", target_sp)
                ) +
                theme_map()
        }

        pA <- make_true_cate_panel(var1, "(A)")
        pB <- make_true_cate_panel(var2, "(B)")

        final_plot <- (pA | pB) +
            plot_annotation(
                title = "Fig 5. Authentic Spatial Heterogeneity via Causal Forests",
                theme = theme(
                    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
                    plot.background = element_rect(fill = "transparent", color = NA),
                    panel.background = element_rect(fill = "transparent", color = NA)
                )
            )

        ggsave(file.path(out_dir, "fig5_true_spatial_cate_map.png"), final_plot, width = 14, height = 7, dpi = 300, bg = "transparent")
        cat("Saved plot to:", file.path(out_dir, "fig5_true_spatial_cate_map.png"), "\n")
    } else {
        cat("Could not find enough significant ATE variables.\n")
    }
} else {
    cat("Missing required CSV files. Ensure CAST species data and grid are present.\n")
}
