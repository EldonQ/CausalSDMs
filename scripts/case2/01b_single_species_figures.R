################################################################################
# 01b_single_species_figures.R
# Publication-quality figures for a SINGLE species CAST pipeline result
#
# Prerequisite: Run 01_cast_pipeline.R first for the target species
#
# Figures produced (aligned with manuscript figure numbering):
#   Fig 2a: DAG network topology (consensus DAG with bootstrap edge strengths)
#   Fig 2b: Causal role grouping diagram (Root / Mediator / Terminal)
#   Fig 2c: ATE forest plot (effect sizes with 95% CI)
#   Fig 2d: CAST variable screening scores (3-component stacked)
#   Fig 3s: Performance comparison bar chart (AUC + TSS, single-species)
#   Fig 8s: CATE spatial heterogeneity heatmaps (causal forests)
#   Panel:  Combined summary panel
#   Table:  Single-species performance summary
################################################################################

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  ► CHANGE THESE TWO PARAMETERS ◄                                       ║
# ╚══════════════════════════════════════════════════════════════════════════╝
target_species <- "swi10"
target_region  <- "SWI"

# ---- Dependencies ----
pkgs <- c(
    "tidyverse", "ggplot2", "patchwork", "scales",
    "igraph", "ggraph",  # DAG network visualization
    "grf",               # Causal forests for CATE
    "ranger",            # RF for CATE nuisance
    "viridis"            # Color scales for CATE maps
)
for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
}

dir.create(sprintf("figures/case2/%s/single_species", target_region),
    recursive = TRUE, showWarnings = FALSE
)

# ==============================================================================
# Publication color scheme & theme (unified)
# ==============================================================================
pal <- list(
    cast    = "#E74C3C", # CAST = red (full pipeline)
    mlp     = "#3498DB", # MLP = blue
    full    = "#95A5A6", # Full variables = grey
    rf      = "#7F8C8D",
    maxent  = "#95A5A6",
    brt     = "#BDC3C7",
    gam     = "#D5DBDB",
    root    = "#E74C3C",
    med     = "#F39C12",
    term    = "#2980B9",
    pos     = "#27AE60",
    neg     = "#E74C3C",
    sig     = "#E74C3C",
    nsig    = "grey60"
)

# 大表配色: CAST突出，其他按变量集区分
model_colors <- c(
    "CAST"   = "#E74C3C",
    "MLP"    = "#3498DB",
    "RF"     = "#2ECC71",
    "Maxent" = "#9B59B6",
    "BRT"    = "#F39C12",
    "GAM"    = "#1ABC9C"
)
varset_shapes <- c("Full + Causal" = 16, "Full (VIF)" = 1)

theme_pub <- function(base_size = 11) {
    theme_minimal(base_size = base_size) +
        theme(
            text = element_text(family = "sans"),
            plot.title = element_text(face = "bold", size = base_size + 2, hjust = 0),
            plot.subtitle = element_text(size = base_size, color = "grey40"),
            axis.title = element_text(face = "bold"),
            axis.text = element_text(color = "grey20"),
            legend.position = "bottom",
            legend.title = element_text(face = "bold", size = base_size - 1),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            strip.text = element_text(face = "bold"),
            plot.margin = margin(10, 15, 10, 15)
        )
}

# ==============================================================================
# Load Data
# ==============================================================================
cat(sprintf("Loading results for %s [%s]...\n", target_species, target_region))

result_base <- sprintf("output/case2/%s/single_species", target_region)
if (!dir.exists(result_base)) {
    result_base <- "output/case2/single_species"
}

# Core data files
results <- read.csv(sprintf("%s/results_%s.csv", result_base, target_species),
    stringsAsFactors = FALSE
)
screening <- read.csv(sprintf("%s/screening_%s.csv", result_base, target_species),
    stringsAsFactors = FALSE
)
roles_file <- sprintf("%s/roles_%s.csv", result_base, target_species)
has_roles <- file.exists(roles_file)
if (has_roles) {
    roles <- read.csv(roles_file, stringsAsFactors = FALSE)
}

# ATE data
ate_file <- sprintf("%s/ate_%s.csv", result_base, target_species)
if (file.exists(ate_file)) {
    ate_data <- read.csv(ate_file, stringsAsFactors = FALSE) %>%
        mutate(ci_lower = coef - 1.96 * se, ci_upper = coef + 1.96 * se)
} else {
    ate_data <- screening %>%
        filter(!is.na(coef)) %>%
        select(variable, coef, p_value, significant) %>%
        mutate(
            se = abs(coef) / pmax(qnorm(1 - p_value / 2), 0.1),
            ci_lower = coef - 1.96 * se,
            ci_upper = coef + 1.96 * se
        )
}

# DAG edges (from pipeline export, or reconstruct from screening)
dag_edges_file <- sprintf("%s/dag_edges_%s.csv", result_base, target_species)
has_dag_edges <- file.exists(dag_edges_file)
if (has_dag_edges) {
    dag_edges <- read.csv(dag_edges_file, stringsAsFactors = FALSE)
    cat(sprintf("  Loaded %d DAG edges from file\n", nrow(dag_edges)))
} else {
    # Fallback: reconstruct from roles out/in degree (limited info)
    cat("  WARNING: dag_edges file not found. DAG network plot will be skipped.\n")
    cat("  Re-run 01_cast_pipeline.R to generate dag_edges CSV.\n")
    dag_edges <- NULL
}

# Training/test data for CATE estimation
train_file <- sprintf("output/case2/%s/train_data_%s.csv", target_region, target_species)
test_file <- sprintf("output/case2/%s/test_data_%s.csv", target_region, target_species)
if (!file.exists(train_file)) {
    train_file <- sprintf("output/case2/train_data_%s.csv", target_species)
    test_file <- sprintf("output/case2/test_data_%s.csv", target_species)
}
has_train_data <- file.exists(train_file)
if (has_train_data) {
    train_raw <- read.csv(train_file, stringsAsFactors = FALSE)
    test_raw <- read.csv(test_file, stringsAsFactors = FALSE)
    cat(sprintf("  Training data: %d rows\n", nrow(train_raw)))
}

# ==============================================================================
# 标准化模型标签 — 统一大表命名
# ==============================================================================
results <- results %>%
    mutate(
        model_label = case_when(
            model == "CAST"    ~ "CAST",
            model == "CI_MLP"  ~ "CAST",     # 兼容旧版
            model == "CGNet"   ~ "CAST",     # 兼容旧版
            model == "MLP"     ~ "MLP",
            model == "FlatNN_cast" ~ "MLP",  # 兼容旧版
            model == "FlatNN_full" ~ "MLP",  # 兼容旧版
            TRUE ~ model
        ),
        varset_label = case_when(
            model_label == "CAST" ~ "Full + Causal",
            var_set == "full+causal" ~ "Full + Causal",
            TRUE ~ "Full (VIF)"
        )
    )

model_order <- c("CAST", "MLP", "RF", "BRT", "Maxent", "GAM")
results$model_label <- factor(results$model_label, levels = model_order)

cat(sprintf("  %d model-group combinations loaded\n", nrow(results)))

fig_dir <- sprintf("figures/case2/%s/single_species", target_region)

# ==============================================================================
# Figure 2a: DAG Network Topology (consensus DAG with bootstrap edge strengths)
# ==============================================================================
cat("\nFig 2a: DAG network topology...\n")

if (!is.null(dag_edges) && nrow(dag_edges) > 0) {
    # Build igraph from edge list
    g <- graph_from_data_frame(
        dag_edges %>% select(from, to),
        directed = TRUE,
        vertices = data.frame(name = unique(c(dag_edges$from, dag_edges$to)))
    )

    # Assign node attributes from roles (if available)
    if (has_roles) {
        V(g)$role <- roles$group[match(V(g)$name, roles$variable)]
        V(g)$role[is.na(V(g)$role)] <- "Other"
        role_node_colors <- c(
            "Root"     = pal$root,
            "Mediator" = pal$med,
            "Terminal"  = pal$term,
            "Other"    = "grey70"
        )
    } else {
        V(g)$role <- "Other"
        role_node_colors <- c("Other" = "grey70")
    }

    # Edge attributes
    E(g)$strength <- dag_edges$strength
    E(g)$direction <- dag_edges$direction

    # Identify CAST-selected variables
    cast_vars <- if (has_roles) roles$variable else screening$variable[1:5]
    V(g)$is_cast <- V(g)$name %in% cast_vars

    fig2a <- ggraph(g, layout = "sugiyama") +
        geom_edge_link(
            aes(
                edge_width = strength,
                edge_alpha = strength
            ),
            arrow = arrow(length = unit(3, "mm"), type = "closed"),
            end_cap = circle(4, "mm"),
            color = "grey30"
        ) +
        geom_node_point(
            aes(color = role, shape = is_cast),
            size = 7
        ) +
        geom_node_text(
            aes(label = name),
            size = 3, repel = TRUE, fontface = "bold"
        ) +
        scale_color_manual(
            values = role_node_colors,
            name = "Causal Role"
        ) +
        scale_shape_manual(
            values = c("TRUE" = 16, "FALSE" = 1),
            name = "CAST Selected",
            labels = c("TRUE" = "Yes", "FALSE" = "No")
        ) +
        scale_edge_width(range = c(0.3, 2.5), name = "Edge Strength") +
        scale_edge_alpha(range = c(0.3, 0.9), guide = "none") +
        labs(
            title = sprintf("Consensus DAG: %s [%s]", target_species, target_region),
            subtitle = sprintf(
                "Bootstrap HC (R=200) | %d strong edges (str≥0.7, dir≥0.6) | Density=%.2f",
                nrow(dag_edges),
                results$dag_density[1]
            )
        ) +
        theme_pub() +
        theme(
            axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank()
        )

    ggsave(sprintf("%s/%s_fig2a_dag_network.png", fig_dir, target_species),
        fig2a, width = 12, height = 9, dpi = 1200, bg = "white"
    )
    ggsave(sprintf("%s/%s_fig2a_dag_network.svg", fig_dir, target_species),
        fig2a, width = 12, height = 9, bg = "white"
    )
    cat("  ✓ fig2a_dag_network\n")
} else {
    cat("  ✗ Skipped (no DAG edges data)\n")
    fig2a <- NULL
}

# ==============================================================================
# Figure 2b: Causal Role Grouping Diagram
# ==============================================================================
cat("\nFig 2b: Causal role grouping...\n")

role_colors <- c("Root" = pal$root, "Mediator" = pal$med, "Terminal" = pal$term)

roles_plot <- roles %>%
    mutate(group = factor(group, levels = c("Root", "Mediator", "Terminal")))

fig2b <- ggplot(roles_plot, aes(
    x = reorder(variable, -role_score),
    y = role_score, fill = group
)) +
    geom_col(width = 0.7, alpha = 0.85, color = "white", linewidth = 0.3) +
    geom_text(aes(label = sprintf("out=%d\nin=%d", out, inp)),
        vjust = -0.3, size = 2.8, lineheight = 0.8
    ) +
    scale_fill_manual(values = role_colors, name = "Causal Role") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(
        title = sprintf("Causal role assignment: %s", target_species),
        subtitle = "Role score = out_degree / (in_degree + 1); Root → Mediator → Terminal",
        x = "Variable", y = "Role score (higher = more upstream)"
    ) +
    theme_pub()

# Add flow label
flow_text <- roles_plot %>%
    group_by(group) %>%
    summarise(vars = paste(variable, collapse = ", "), .groups = "drop") %>%
    arrange(match(group, c("Root", "Mediator", "Terminal")))
flow_label <- paste(
    sprintf("%s [%s]", flow_text$group, flow_text$vars),
    collapse = "  →  "
)
fig2b <- fig2b +
    labs(caption = sprintf("Information flow: %s", flow_label)) +
    theme(plot.caption = element_text(face = "bold", size = 10, color = "#2C3E50", hjust = 0.5))

ggsave(sprintf("%s/%s_fig2b_causal_roles.png", fig_dir, target_species),
    fig2b, width = 10, height = 6, dpi = 1200, bg = "white"
)
ggsave(sprintf("%s/%s_fig2b_causal_roles.svg", fig_dir, target_species),
    fig2b, width = 10, height = 6, bg = "white"
)
cat("  ✓ fig2b_causal_roles\n")

# ==============================================================================
# Figure 2c: ATE Forest Plot (effect sizes with 95% CI)
# ==============================================================================
cat("\nFig 2c: ATE forest plot...\n")

ate_plot <- ate_data %>%
    mutate(
        sig_marker = ifelse(significant, "Significant (p<0.05)", "Not significant"),
        var_label = sprintf("%s %s", variable, ifelse(significant, "*", ""))
    )

fig2c <- ggplot(ate_plot, aes(
    x = coef, y = reorder(variable, abs(coef)),
    color = sig_marker
)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_point(size = 3) +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper),
        height = 0.25, linewidth = 0.6
    ) +
    scale_color_manual(values = c(
        "Significant (p<0.05)" = pal$sig,
        "Not significant" = pal$nsig
    ), name = "") +
    labs(
        title = sprintf("Average Treatment Effect (ATE): %s", target_species),
        subtitle = "DML with 2-fold cross-fitting | Median binarization | 95% CI",
        x = "ATE coefficient (effect on presence probability)", y = ""
    ) +
    theme_pub() +
    theme(panel.grid.major.y = element_line(color = "grey92"))

ggsave(sprintf("%s/%s_fig2c_ate_forest.png", fig_dir, target_species),
    fig2c, width = 10, height = 6, dpi = 1200, bg = "white"
)
ggsave(sprintf("%s/%s_fig2c_ate_forest.svg", fig_dir, target_species),
    fig2c, width = 10, height = 6, bg = "white"
)
cat("  ✓ fig2c_ate_forest\n")

# ==============================================================================
# Figure 2d: CAST Variable Screening Scores (3-component stacked bar)
# ==============================================================================
cat("\nFig 2d: Variable screening scores...\n")

# Retrieve adaptive weights from screening scores (reverse-engineer or use fixed)
# The pipeline uses adaptive weights; approximate from data
w_dag_approx <- 0.15 + 0.15 * (1 - results$dag_density[1])
r_sig <- sum(ate_data$significant) / nrow(ate_data)
w_ate_approx <- 0.25 + 0.25 * r_sig
w_imp_approx <- 1 - w_dag_approx - w_ate_approx

screen_long <- screening %>%
    select(variable, score_dag, score_ate, score_imp, score_total) %>%
    mutate(
        w_dag = w_dag_approx * score_dag,
        w_ate = w_ate_approx * score_ate,
        w_imp = w_imp_approx * score_imp
    ) %>%
    pivot_longer(
        cols = c(w_dag, w_ate, w_imp),
        names_to = "component", values_to = "weighted_score"
    ) %>%
    mutate(component = case_when(
        component == "w_dag" ~ sprintf("DAG Out-degree (×%.2f)", w_dag_approx),
        component == "w_ate" ~ sprintf("ATE Effect (×%.2f)", w_ate_approx),
        component == "w_imp" ~ sprintf("RF Importance (×%.2f)", w_imp_approx)
    ))

# Determine threshold
sorted_scores <- sort(screening$score_total, decreasing = TRUE)
if (length(sorted_scores) > 4) {
    score_drops <- -diff(sorted_scores)
    elbow_idx <- which.max(score_drops)
    elbow_threshold <- sorted_scores[elbow_idx]
    n_by_elbow <- sum(screening$score_total >= elbow_threshold)
    if (n_by_elbow < 3 || n_by_elbow > 0.8 * nrow(screening)) {
        sel_threshold <- median(sorted_scores)
        thr_label <- "median"
    } else {
        sel_threshold <- elbow_threshold
        thr_label <- "elbow"
    }
} else {
    sel_threshold <- median(sorted_scores)
    thr_label <- "median"
}
cast_vars_sel <- screening$variable[screening$score_total >= sel_threshold]
if (length(cast_vars_sel) < 3) {
    cast_vars_sel <- screening$variable[1:min(3, nrow(screening))]
}

comp_labels <- unique(screen_long$component)
fig2d <- ggplot(screen_long, aes(
    x = reorder(variable, -weighted_score),
    y = weighted_score, fill = component
)) +
    geom_col(position = "stack", width = 0.7, alpha = 0.85) +
    geom_hline(
        yintercept = sel_threshold, linetype = "dashed",
        color = pal$cast, linewidth = 0.8
    ) +
    annotate("text",
        x = nrow(screening) - 0.5, y = sel_threshold + 0.02,
        label = sprintf("Threshold = %.3f (%s)", sel_threshold, thr_label),
        color = pal$cast, fontface = "italic", hjust = 1, size = 3
    ) +
    scale_fill_manual(values = setNames(
        c("#2C3E50", "#E67E22", "#27AE60"),
        comp_labels
    ), name = "Weighted Component") +
    labs(
        title = sprintf("CAST variable screening: %s", target_species),
        subtitle = sprintf(
            "%d/%d variables selected (%s threshold)",
            length(cast_vars_sel), nrow(screening), thr_label
        ),
        x = "Environmental variable", y = "Weighted screening score"
    ) +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(sprintf("%s/%s_fig2d_screening.png", fig_dir, target_species),
    fig2d, width = 10, height = 6, dpi = 1200, bg = "white"
)
ggsave(sprintf("%s/%s_fig2d_screening.svg", fig_dir, target_species),
    fig2d, width = 10, height = 6, bg = "white"
)
cat("  ✓ fig2d_screening\n")

# ==============================================================================
# Figure 3s: 统一大表风格性能对比（按AUC排序，CAST居首）
# ==============================================================================
cat("\nFig 3s: Performance comparison (unified ranking)...\n")

# 构建排序标签: "Model (VarSet)"
results <- results %>%
    mutate(
        display_label = ifelse(model_label == "CAST",
            "CAST",
            sprintf("%s (%s)", model_label, varset_label)),
        is_cast = (model_label == "CAST")
    )

# 按AUC降序排列；因子水平去重（避免多变量集时重复 level 报错）
results_ranked <- results %>% arrange(desc(auc_mean))
levels_order <- rev(unique(as.character(results_ranked$display_label)))
results_ranked$display_label <- factor(results_ranked$display_label, levels = levels_order)

fig3s <- ggplot(results_ranked, aes(
    x = display_label, y = auc_mean,
    fill = as.character(model_label)
)) +
    geom_col(width = 0.7, alpha = 0.9, color = "white", linewidth = 0.3) +
    geom_errorbar(aes(ymin = auc_mean - auc_sd, ymax = auc_mean + auc_sd),
        width = 0.25, linewidth = 0.5, color = "grey30"
    ) +
    geom_text(aes(label = sprintf("%.4f", auc_mean)),
        hjust = -0.15, size = 3, fontface = "bold"
    ) +
    coord_flip() +
    scale_fill_manual(values = model_colors, name = "Model") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
        title = sprintf("Performance ranking: %s [%s]", target_species, target_region),
        subtitle = "CAST integrates causal screening and causally-informed feature engineering;\nall other models use standard feature spaces.",
        x = "", y = "AUC"
    ) +
    theme_pub() +
    theme(legend.position = "right")

ggsave(sprintf("%s/%s_fig3s_performance.png", fig_dir, target_species),
    fig3s, width = 12, height = 7, dpi = 1200, bg = "white"
)
ggsave(sprintf("%s/%s_fig3s_performance.svg", fig_dir, target_species),
    fig3s, width = 12, height = 7, bg = "white"
)
cat("  ✓ fig3s_performance\n")

# ==============================================================================
# Figure 8s: CATE Spatial Heterogeneity Maps (Causal Forests)
# ==============================================================================
cat("\nFig 8s: CATE spatial heterogeneity maps (causal forests)...\n")

# Identify significant causal driver variables
sig_vars <- ate_data %>% filter(significant) %>% pull(variable)

if (has_train_data && length(sig_vars) > 0) {
    # Prepare data
    env_cols <- setdiff(names(train_raw), "presence")
    for (col in env_cols) {
        train_raw[[col]] <- as.numeric(train_raw[[col]])
        test_raw[[col]] <- as.numeric(test_raw[[col]])
    }

    # Use all post-VIF variables for confounders
    all_vars <- screening$variable
    Y <- train_raw$presence
    X_all <- train_raw[, all_vars, drop = FALSE]
    X_all[is.na(X_all)] <- 0

    # Predict CATE on test data (which has independent spatial locations)
    X_test_all <- test_raw[, all_vars, drop = FALSE]
    X_test_all[is.na(X_test_all)] <- 0

    cate_results <- list()
    cate_plots <- list()

    # Limit to top 4 significant variables by |ATE|
    sig_ordered <- ate_data %>%
        filter(significant) %>%
        arrange(desc(abs(coef)))
    sig_vars_top <- head(sig_ordered$variable, 4)

    for (v in sig_vars_top) {
        cat(sprintf("  Training causal forest for: %s\n", v))
        tryCatch({
            # Binarize treatment at median (same as DML stage)
            W_bin <- as.integer(X_all[[v]] > median(X_all[[v]], na.rm = TRUE))
            X_confounders <- as.matrix(X_all[, setdiff(all_vars, v), drop = FALSE])

            # Train causal forest
            set.seed(42)
            cf <- grf::causal_forest(
                X = X_confounders,
                Y = Y,
                W = W_bin,
                num.trees = 2000,
                honesty = TRUE,
                seed = 42
            )

            # Predict CATE on test locations
            X_test_conf <- as.matrix(X_test_all[, setdiff(all_vars, v), drop = FALSE])
            cate_pred <- predict(cf, X_test_conf, estimate.variance = TRUE)
            tau_hat <- cate_pred$predictions
            tau_var <- cate_pred$variance.estimates

            # Store results
            cate_df <- data.frame(
                variable = v,
                tau_hat = tau_hat,
                tau_se = sqrt(tau_var),
                row_idx = 1:length(tau_hat)
            )

            # Add spatial coordinates if available (test data index as proxy)
            # Use first two principal env variables as spatial proxy
            pc <- prcomp(X_test_all, center = TRUE, scale. = TRUE)
            cate_df$pc1 <- pc$x[, 1]
            cate_df$pc2 <- pc$x[, 2]

            cate_results[[v]] <- cate_df

            # Global ATE from DML for reference
            global_ate <- ate_data$coef[ate_data$variable == v]

            p <- ggplot(cate_df, aes(x = pc1, y = pc2, color = tau_hat)) +
                geom_point(size = 1.5, alpha = 0.8) +
                scale_color_gradient2(
                    low = "#2166AC", mid = "white", high = "#B2182B",
                    midpoint = 0,
                    name = expression(hat(tau)(x[i]))
                ) +
                geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
                geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
                labs(
                    title = sprintf("CATE: %s", v),
                    subtitle = sprintf(
                        "Global ATE=%.4f | CATE range: [%.4f, %.4f]",
                        global_ate, min(tau_hat), max(tau_hat)
                    ),
                    x = "Environmental PC1", y = "Environmental PC2"
                ) +
                theme_pub() +
                theme(legend.position = "right")

            cate_plots[[v]] <- p
            cat(sprintf("    CATE range: [%.4f, %.4f], mean=%.4f\n",
                min(tau_hat), max(tau_hat), mean(tau_hat)))
        }, error = function(e) {
            cat(sprintf("    FAILED for %s: %s\n", v, e$message))
        })
    }

    # Combine CATE panels
    if (length(cate_plots) > 0) {
        n_plots <- length(cate_plots)
        ncols <- min(2, n_plots)
        nrows <- ceiling(n_plots / ncols)

        fig8s <- wrap_plots(cate_plots, ncol = ncols) +
            plot_annotation(
                title = sprintf(
                    "Spatially Heterogeneous Causal Effects (CATE): %s [%s]",
                    target_species, target_region
                ),
                subtitle = sprintf(
                    "Causal forests (2000 trees) | %d significant variables | Median binarization",
                    length(cate_plots)
                ),
                theme = theme(
                    plot.title = element_text(face = "bold", size = 14),
                    plot.subtitle = element_text(size = 11, color = "grey40")
                )
            )

        ggsave(sprintf("%s/%s_fig8s_cate_maps.png", fig_dir, target_species),
            fig8s,
            width = 7 * ncols, height = 6 * nrows, dpi = 1200, bg = "white"
        )
        ggsave(sprintf("%s/%s_fig8s_cate_maps.svg", fig_dir, target_species),
            fig8s,
            width = 7 * ncols, height = 6 * nrows, bg = "white"
        )

        # Save CATE data for manuscript reporting
        cate_all <- do.call(rbind, cate_results)
        write.csv(cate_all, sprintf("%s/cate_%s.csv", result_base, target_species),
            row.names = FALSE
        )

        cat(sprintf("  ✓ fig8s_cate_maps (%d panels)\n", n_plots))
    } else {
        cat("  ✗ No CATE plots generated (all causal forests failed)\n")
        fig8s <- NULL
    }
} else {
    if (!has_train_data) {
        cat("  ✗ Skipped (training data not found)\n")
    } else {
        cat("  ✗ Skipped (no significant ATE variables)\n")
    }
    fig8s <- NULL
    cate_plots <- list()
}

# ==============================================================================
# Combined summary panel (Fig 2 composite: DAG + Roles + ATE + Screening)
# ==============================================================================
cat("\nCombined panel: Fig 2 composite...\n")

# Build Fig 2 as a 2x2 panel
panel_list <- list()
panel_labels <- c()

if (!is.null(fig2a)) {
    # For the combined panel, use a simplified DAG
    panel_list[["a"]] <- fig2a + labs(title = NULL, subtitle = NULL) +
        ggtitle("(a) Consensus DAG")
    panel_labels <- c(panel_labels, "a")
}
panel_list[["b"]] <- fig2b + labs(title = NULL, subtitle = NULL) +
    ggtitle("(b) Causal Roles")
panel_list[["c"]] <- fig2c + labs(title = NULL, subtitle = NULL) +
    ggtitle("(c) ATE Forest Plot")
panel_list[["d"]] <- fig2d + labs(title = NULL, subtitle = NULL) +
    ggtitle("(d) CAST Screening Scores")

if (length(panel_list) == 4) {
    fig2_combined <- (panel_list[["a"]] | panel_list[["b"]]) /
        (panel_list[["c"]] | panel_list[["d"]]) +
        plot_annotation(
            title = sprintf(
                "CAST Causal Structure Analysis: %s [%s]",
                target_species, target_region
            ),
            theme = theme(
                plot.title = element_text(face = "bold", size = 16)
            )
        )
} else {
    fig2_combined <- (panel_list[["b"]] | panel_list[["d"]]) /
        (panel_list[["c"]] | fig3s) +
        plot_annotation(
            title = sprintf(
                "CAST Pipeline Summary: %s [%s]",
                target_species, target_region
            ),
            theme = theme(
                plot.title = element_text(face = "bold", size = 16)
            )
        )
}

ggsave(sprintf("%s/%s_fig2_combined.png", fig_dir, target_species),
    fig2_combined,
    width = 18, height = 14, dpi = 1200, bg = "white"
)
ggsave(sprintf("%s/%s_fig2_combined.svg", fig_dir, target_species),
    fig2_combined,
    width = 18, height = 14, bg = "white"
)
cat("  ✓ fig2_combined\n")

# ==============================================================================
# Full summary panel (all figures)
# ==============================================================================
cat("\nFull summary panel...\n")

fig_summary <- (fig2d | fig2b) / (fig2c | fig3s) +
    plot_annotation(
        title = sprintf("CAST Pipeline Summary: %s [%s]", target_species, target_region),
        subtitle = "Causal screening → Role grouping → ATE estimation → Model comparison",
        theme = theme(
            plot.title = element_text(face = "bold", size = 16),
            plot.subtitle = element_text(size = 12, color = "grey40")
        )
    ) +
    plot_layout(heights = c(1, 1))

ggsave(sprintf("%s/%s_summary_panel.png", fig_dir, target_species),
    fig_summary,
    width = 18, height = 14, dpi = 1200, bg = "white"
)
ggsave(sprintf("%s/%s_summary_panel.svg", fig_dir, target_species),
    fig_summary,
    width = 18, height = 14, bg = "white"
)
cat("  ✓ summary_panel\n")

# ==============================================================================
# TABLE: Single-species performance summary
# ==============================================================================
cat("\nGenerating performance table...\n")

pub_table <- results %>%
    arrange(desc(auc_mean)) %>%
    mutate(
        AUC = sprintf("%.4f ± %.4f", auc_mean, auc_sd),
        TSS = sprintf("%.4f ± %.4f", tss_mean, tss_sd),
        Rank = row_number()
    ) %>%
    select(Rank, model_label, varset_label, n_vars, AUC, TSS)

write.csv(pub_table,
    sprintf("%s/table_%s.csv", result_base, target_species),
    row.names = FALSE
)

# 打印统一大表
cat(sprintf("\n  ╔══════════════════════════════════════════════════════════════════════════════════╗\n"))
cat(sprintf("  ║  Performance Ranking: %-10s [%s]                                            ║\n",
    target_species, target_region))
cat("  ╠══════════════════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("  ║  %-4s  %-10s  %-16s  %-4s  %-16s  %-16s  ║\n",
    "Rank", "Model", "Variable Set", "Vars", "AUC", "TSS"))
cat("  ╠══════════════════════════════════════════════════════════════════════════════════╣\n")
for (i in 1:nrow(pub_table)) {
    r <- pub_table[i, ]
    mk <- if (as.character(r$model_label) == "CAST") " ★" else "  "
    cat(sprintf(
        "  ║%s%-4d  %-10s  %-16s  %4d  %-16s  %-16s  ║\n",
        mk, r$Rank, as.character(r$model_label),
        as.character(r$varset_label), r$n_vars, r$AUC, r$TSS
    ))
}
cat("  ╚══════════════════════════════════════════════════════════════════════════════════╝\n")

n_total_vars <- nrow(screening)
n_cast_features <- results$n_features_total[results$model_label == "CAST"][1]
n_interactions <- results$n_interactions[results$model_label == "CAST"][1]

cat("\n  ╔══════════════════════════════════════════════════════════════════════╗\n")
cat("  ║  CAST Pipeline Summary                                             ║\n")
cat("  ╠══════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf(
    "  ║  VIF filtering: → %d variables                                    ║\n",
    n_total_vars
))
cat(sprintf(
    "  ║  ATE:           %d/%d significant (p<0.05)                        ║\n",
    sum(ate_data$significant), nrow(ate_data)
))
cat(sprintf(
    "  ║  CAST features: %d base + %d DAG interactions = %d total          ║\n",
    n_total_vars, ifelse(is.na(n_interactions), 0, n_interactions),
    ifelse(is.na(n_cast_features), n_total_vars, n_cast_features)
))
if (has_roles) {
    flow_parts <- roles %>%
        group_by(group) %>%
        summarise(vars = paste(variable, collapse = ","), .groups = "drop")
    cat(sprintf(
        "  ║  Roles:         %s  ║\n",
        paste(sprintf("%s[%s]", flow_parts$group, flow_parts$vars), collapse = " → ")
    ))
}
if (length(cate_plots) > 0) {
    cat(sprintf(
        "  ║  CATE:          %d variables with spatial heterogeneity maps     ║\n",
        length(cate_plots)
    ))
}
if (!is.null(dag_edges)) {
    cat(sprintf(
        "  ║  DAG:           %d strong edges, density=%.3f                   ║\n",
        nrow(dag_edges), results$dag_density[1]
    ))
}
cat("  ╚══════════════════════════════════════════════════════════════════════╝\n")

# ==============================================================================
# Summary of all generated files
# ==============================================================================
cat(sprintf("\n======================================================================\n"))
cat(sprintf("  Single-Species Figures Complete: %s [%s]\n", target_species, target_region))
cat("======================================================================\n")
cat(sprintf("  Figures saved to %s/:\n", fig_dir))
if (!is.null(fig2a)) {
    cat(sprintf("    %s_fig2a_dag_network      (.png + .pdf)  [NEW: DAG topology]\n", target_species))
}
cat(sprintf("    %s_fig2b_causal_roles     (.png + .pdf)\n", target_species))
cat(sprintf("    %s_fig2c_ate_forest       (.png + .pdf)\n", target_species))
cat(sprintf("    %s_fig2d_screening        (.png + .pdf)\n", target_species))
cat(sprintf("    %s_fig3s_performance      (.png + .pdf)\n", target_species))
if (!is.null(fig8s)) {
    cat(sprintf("    %s_fig8s_cate_maps        (.png + .pdf)  [NEW: CATE heatmaps]\n", target_species))
}
cat(sprintf("    %s_fig2_combined          (.png + .pdf)\n", target_species))
cat(sprintf("    %s_summary_panel          (.png + .pdf)\n", target_species))
cat(sprintf("  Data saved to %s/:\n", result_base))
cat(sprintf("    table_%s.csv\n", target_species))
if (length(cate_plots) > 0) {
    cat(sprintf("    cate_%s.csv               [NEW: CATE estimates]\n", target_species))
}
cat("======================================================================\n")
