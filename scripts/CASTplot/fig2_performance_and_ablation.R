# ==============================================================================
# Fig 2: 因果筛选的简约性与结构拓扑特征编码增益 (Results 3.1 & 3.2)
# ==============================================================================

rm(list = ls())
setwd("E:/CausalSDMs")

out_dir <- "figures/CASTplot"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

library(tidyverse)
library(patchwork)

theme_cast <- function(base_size = 12) {
    theme_minimal(base_size = base_size, base_family = "sans") +
        theme(
            panel.grid.minor = element_blank(),
            axis.title       = element_text(face = "bold"),
            plot.title       = element_text(face = "bold", hjust = 0.5),
            plot.subtitle    = element_text(hjust = 0.5, color = "grey40")
        )
}

# ---- Load Real EcoISEA3H Data ----
d <- read.csv("output/case2_eco/all_results_v3.csv", stringsAsFactors = FALSE) %>% filter(!is.na(auc_mean))
sp_meta <- read.csv("outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv", stringsAsFactors = FALSE) %>%
    mutate(species = gsub(" ", "_", species))

d <- d %>%
    left_join(sp_meta %>% select(species, family), by = "species") %>%
    mutate(family = ifelse(is.na(family), "Unknown", family))

ablation_colors <- c("CAST" = "#E74C3C", "MLP_ATE" = "#5DADE2", "MLP" = "#AED6F1")
all_model_colors <- c("CAST" = "#E74C3C", "MLP_ATE" = "#5DADE2", "MLP" = "#AED6F1", "RF" = "#27AE60", "BRT" = "#E67E22", "Maxent" = "#9B59B6")

# ------------------------------------------------------------------------------
# (a) Ablation violin: MLP → MLP_ATE → CAST
# ------------------------------------------------------------------------------
d_abl <- d %>%
    filter(model %in% c("MLP", "MLP_ATE", "CAST")) %>%
    mutate(model = factor(model, levels = c("MLP", "MLP_ATE", "CAST")))

abl_stats <- d_abl %>%
    group_by(model) %>%
    summarise(mean_auc = mean(auc_mean), .groups = "drop")

pa <- ggplot(d_abl, aes(x = model, y = auc_mean, fill = model)) +
    geom_violin(scale = "width", alpha = 0.6, trim = TRUE) +
    geom_boxplot(width = 0.15, fill = "white", alpha = 0.8, outlier.size = 0.8) +
    stat_summary(fun = mean, geom = "point", size = 3, color = "black", shape = 18) +
    geom_text(data = abl_stats, aes(x = model, y = mean_auc, label = sprintf("%.4f", mean_auc)), vjust = -1.5, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = ablation_colors, guide = "none") +
    labs(title = "(A) Causal ablation: Structural Gain", subtitle = "Network constraint: MLP → + ATE w → + DAG info", x = "Neural Network Variant", y = "Test AUC") +
    theme_cast()

# ------------------------------------------------------------------------------
# (B) Paired species lines: MLP (FlatNN) → CAST (CI-MLP)
# ------------------------------------------------------------------------------
paired_wide <- d %>%
    filter(model %in% c("MLP", "CAST")) %>%
    select(species, model, auc_mean) %>%
    pivot_wider(names_from = model, values_from = auc_mean, values_fn = max) %>%
    drop_na() %>%
    mutate(delta = CAST - MLP, direction = ifelse(delta >= 0, "CAST better", "MLP better"))

n_wins <- sum(paired_wide$delta >= 0)

pb <- ggplot(paired_wide) +
    geom_segment(aes(x = 1, xend = 2, y = MLP, yend = CAST, color = direction), alpha = 0.5, linewidth = 0.5) +
    geom_point(aes(x = 1, y = MLP), color = "#AED6F1", size = 1.8) +
    geom_point(aes(x = 2, y = CAST), color = "#E74C3C", size = 1.8) +
    scale_color_manual(values = c("CAST better" = "#E74C3C", "MLP better" = "#2C3E50"), name = "") +
    scale_x_continuous(breaks = c(1, 2), labels = c("MLP (No Topology)", "CI-MLP (Feature Encoded)"), limits = c(0.7, 2.3)) +
    labs(title = "(B) Per-species structural effect", subtitle = sprintf("CI-MLP higher in %d/%d (%.0f%%) species", n_wins, nrow(paired_wide), 100 * n_wins / nrow(paired_wide)), x = "", y = "AUC") +
    theme_cast() +
    theme(legend.position = "bottom")

# ------------------------------------------------------------------------------
# (C) Interaction Complexity vs AUC Ratio
# ------------------------------------------------------------------------------
cast_mlp <- d %>%
    filter(model %in% c("CAST", "MLP")) %>%
    select(species, family, model, auc_mean, n_interactions) %>%
    pivot_wider(names_from = model, values_from = c(auc_mean, n_interactions), names_glue = "{model}_{.value}") %>%
    drop_na() %>%
    mutate(auc_ratio = CAST_auc_mean / MLP_auc_mean)

pc <- ggplot(cast_mlp, aes(x = CAST_n_interactions, y = auc_ratio)) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey40") +
    geom_point(aes(color = family), size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "dashed") +
    scale_color_brewer(palette = "Set1", name = "Taxonomic Family") +
    labs(title = "(C) Interactions & Performance", subtitle = "More explicit structural features drive ratio > 1", x = "Number of DAG-guided interactive features", y = "CI-MLP AUC / FlatNN AUC") +
    theme_cast()

# ------------------------------------------------------------------------------
# (D) Full-model context: compact mean ± SE dot plot
# ------------------------------------------------------------------------------
full_stats <- d %>%
    mutate(model = factor(model, levels = c("CAST", "MLP_ATE", "MLP", "RF", "BRT", "Maxent"))) %>%
    group_by(model) %>%
    summarise(mean_auc = mean(auc_mean), se_auc = sd(auc_mean) / sqrt(n()), .groups = "drop")

pd <- ggplot(full_stats, aes(x = model, y = mean_auc, color = model)) +
    geom_pointrange(aes(ymin = mean_auc - se_auc, ymax = mean_auc + se_auc), size = 0.8, fatten = 3) +
    geom_text(aes(label = sprintf("%.4f", mean_auc)), vjust = -1.2, size = 3, fontface = "bold") +
    scale_color_manual(values = all_model_colors, guide = "none") +
    annotate("rect", xmin = 0.5, xmax = 3.5, ymin = -Inf, ymax = Inf, alpha = 0.05, fill = "#E74C3C") +
    labs(title = "(D) All evaluated baselines", subtitle = "CAST family vs Traditional", x = "", y = "Mean AUC") +
    theme_cast()

# === Combine and Save ===
final_plot <- (pa | pb) / (pc | pd) +
    plot_annotation(title = "Fig 2. Causal Enhancement: Predictive Parity with Parsimony and Structural Gain", theme = theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)))
ggsave(file.path(out_dir, "fig2_performance_and_ablation.png"), final_plot, width = 14, height = 10, dpi = 300, bg = "white")
cat("Saved plot to:", file.path(out_dir, "fig2_performance_and_ablation.png"), "\n")
