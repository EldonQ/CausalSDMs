################################################################################
# Fig S3: Model Training Dynamics — Learning Curves & Convergence Diagnostics
#
# Shows that CAST, MLP_ATE, and MLP models converge stably across species,
# and that the CI-MLP architecture does not overfit despite additional
# causal interaction features.
#
# Panel (a): Aggregated training loss curves (median ± IQR across species)
# Panel (b): Aggregated validation AUC curves
# Panel (c): Best-epoch distribution per model (convergence speed comparison)
# Panel (d): Species-level final val AUC heatmap (model × species)
#
# Data required:
#   output/case2_eco/all_learning_curves_v3.csv
#   output/case2_eco/all_results_v3.csv
#
# Run: setwd("E:/CausalSDMs")
#      source("scripts/EcoISEA3H/plot/figS3_training_curves.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

fig_dir <- "figures/case2_eco/plot"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(patchwork)

# ── Theme ────────────────────────────────────────────────────────────────────
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

model_colors <- c(
    "CAST"    = "#E74C3C",
    "MLP_ATE" = "#F39C12",
    "MLP"     = "#3498DB"
)
model_labels <- c(
    "CAST"    = "CAST (CI-MLP)",
    "MLP_ATE" = "MLP + ATE weights",
    "MLP"     = "Flat MLP (baseline)"
)

# ── Load data ────────────────────────────────────────────────────────────────
lc <- read.csv("output/case2_eco/all_learning_curves_v3.csv",
    stringsAsFactors = FALSE
) %>%
    mutate(
        model = factor(model, levels = c("CAST", "MLP_ATE", "MLP")),
        train_loss = as.numeric(train_loss),
        val_auc = as.numeric(val_auc)
    ) %>%
    filter(!is.na(model))

results <- read.csv("output/case2_eco/all_results_v3.csv",
    stringsAsFactors = FALSE
)

n_species <- length(unique(lc$species))
n_runs <- length(unique(lc$run))
cat(sprintf(
    "[FigS3] Learning curves: %d species, %d runs/model, %d total rows\n",
    n_species, n_runs, nrow(lc)
))

# ══════════════════════════════════════════════════════════════════════════════
# Panel (a): Aggregated training loss curves (median ± IQR)
# ══════════════════════════════════════════════════════════════════════════════
loss_agg <- lc %>%
    filter(!is.na(train_loss)) %>%
    group_by(model, epoch) %>%
    summarise(
        med = median(train_loss, na.rm = TRUE),
        q25 = quantile(train_loss, 0.25, na.rm = TRUE),
        q75 = quantile(train_loss, 0.75, na.rm = TRUE),
        .groups = "drop"
    )

pa <- ggplot(loss_agg, aes(x = epoch, y = med, color = model, fill = model)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8, alpha = 0.9) +
    scale_color_manual(values = model_colors, labels = model_labels, name = "") +
    scale_fill_manual(values = model_colors, labels = model_labels, guide = "none") +
    scale_y_log10() +
    labs(
        title = "(a) Training loss convergence",
        subtitle = sprintf(
            "Median ± IQR across %d species × %d runs | log₁₀ scale",
            n_species, n_runs
        ),
        x = "Epoch", y = "Focal Loss (log₁₀)"
    ) +
    theme_pub() +
    theme(legend.position = c(0.75, 0.85))

# ══════════════════════════════════════════════════════════════════════════════
# Panel (b): Aggregated validation AUC curves
# ══════════════════════════════════════════════════════════════════════════════
auc_agg <- lc %>%
    filter(!is.na(val_auc), val_auc > 0) %>%
    group_by(model, epoch) %>%
    summarise(
        med = median(val_auc, na.rm = TRUE),
        q25 = quantile(val_auc, 0.25, na.rm = TRUE),
        q75 = quantile(val_auc, 0.75, na.rm = TRUE),
        .groups = "drop"
    )

pb <- ggplot(auc_agg, aes(x = epoch, y = med, color = model, fill = model)) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.15, color = NA) +
    geom_line(linewidth = 0.8, alpha = 0.9) +
    scale_color_manual(values = model_colors, labels = model_labels, name = "") +
    scale_fill_manual(values = model_colors, labels = model_labels, guide = "none") +
    labs(
        title = "(b) Validation AUC convergence",
        subtitle = "Median ± IQR | Higher = better predictive discrimination",
        x = "Epoch", y = "Validation AUC"
    ) +
    theme_pub() +
    theme(legend.position = c(0.75, 0.20))

# ══════════════════════════════════════════════════════════════════════════════
# Panel (c): Best-epoch distribution (convergence speed)
# ══════════════════════════════════════════════════════════════════════════════
best_epochs <- lc %>%
    filter(!is.na(val_auc)) %>%
    group_by(model, species, run) %>%
    slice_max(val_auc, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(model, species, run, epoch, val_auc)

model_medians <- best_epochs %>%
    group_by(model) %>%
    summarise(
        med_epoch = median(epoch),
        med_auc = sprintf("%.4f", median(val_auc)),
        .groups = "drop"
    )

pc <- ggplot(best_epochs, aes(x = model, y = epoch, fill = model)) +
    geom_violin(scale = "width", alpha = 0.5, linewidth = 0.3, trim = TRUE) +
    geom_boxplot(
        width = 0.12, fill = "white", alpha = 0.85,
        outlier.size = 0.6, outlier.alpha = 0.4
    ) +
    geom_text(
        data = model_medians,
        aes(x = model, y = med_epoch, label = sprintf("med=%d", as.integer(med_epoch))),
        vjust = -2.0, size = 3.0, fontface = "bold", color = "black",
        inherit.aes = FALSE
    ) +
    scale_fill_manual(values = model_colors, guide = "none") +
    scale_x_discrete(labels = model_labels) +
    labs(
        title = "(c) Convergence speed",
        subtitle = "Epoch at best validation AUC | Lower = faster convergence",
        x = "", y = "Best epoch"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 15, hjust = 1, size = 9))

# ══════════════════════════════════════════════════════════════════════════════
# Panel (d): Species × model AUC heatmap
# ══════════════════════════════════════════════════════════════════════════════
nn_results <- results %>%
    filter(model %in% c("CAST", "MLP_ATE", "MLP"), !is.na(auc_mean)) %>%
    select(species, model, auc_mean) %>%
    mutate(model = factor(model, levels = c("CAST", "MLP_ATE", "MLP")))

species_order <- nn_results %>%
    filter(model == "CAST") %>%
    arrange(desc(auc_mean)) %>%
    pull(species)

nn_results$species <- factor(nn_results$species, levels = rev(species_order))

pd <- ggplot(nn_results, aes(x = model, y = species, fill = auc_mean)) +
    geom_tile(color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("%.3f", auc_mean)),
        size = 2.2, color = "black"
    ) +
    scale_fill_gradient2(
        low = "#FDEBD0", mid = "#F5B041", high = "#C0392B",
        midpoint = median(nn_results$auc_mean, na.rm = TRUE),
        name = "AUC"
    ) +
    scale_x_discrete(labels = model_labels) +
    labs(
        title = "(d) Final test AUC by species",
        subtitle = "Heatmap: species (rows) × model (columns)",
        x = "", y = ""
    ) +
    theme_pub() +
    theme(
        axis.text.x = element_text(angle = 30, hjust = 1, size = 9),
        axis.text.y = element_text(size = 7),
        legend.position = "right"
    )

# ══════════════════════════════════════════════════════════════════════════════
# Combine
# ══════════════════════════════════════════════════════════════════════════════
figS3 <- (pa | pb) / (pc | pd) +
    plot_layout(heights = c(1, 1.2)) +
    plot_annotation(
        title = "Fig S3  Model training dynamics & convergence diagnostics",
        subtitle = paste0(
            "CI-MLP (CAST) converges stably despite additional causal interaction features. ",
            "Validation AUC plateaus consistently across all ", n_species, " species."
        ),
        theme = theme(
            plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
            plot.subtitle = element_text(
                face = "italic", size = 9.5, hjust = 0.5,
                color = "grey40"
            )
        )
    )

ggsave(file.path(fig_dir, "figS3_training_curves.png"),
    figS3,
    width = 15, height = 13, dpi = 300, bg = "white"
)
cat("✓ Saved figS3_training_curves.png\n")
