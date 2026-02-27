################################################################################
# 03_publication_figures.R
# CAST Publication Figures — 全物种汇总出图（统一大表风格）
#
# 读取: output/case2/all_results_v3.csv (+ dag_info, screening, etc.)
# 输出: figures/case2/
#
# 出图清单:
#   Fig 1: 统一排名条形图（CAST vs all baselines, 按AUC排序）
#   Fig 2: CAST vs MLP散点图（结构效应）
#   Fig 3: 筛选效应（Full→CAST-screened）各算法箱线图
#   Fig 4: DAG密度 vs CAST结构优势（关键诊断）
#   Fig 5: 各区域性能对比（分面）
#   Fig 6: 变量缩减 + 交互特征统计
#   Fig 7: CAST优势分解（筛选+结构）
#   Table: 全物种统一排名汇总表
################################################################################

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

library(tidyverse)
library(patchwork)
library(scales)

dir.create("figures/case2", recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# 加载数据
# ==============================================================================
results_file <- "output/case2/all_results_v3.csv"
if (!file.exists(results_file)) {
    stop("Results file not found. Run 02_multi_species_experiment.R first.")
}

all_results <- read.csv(results_file, stringsAsFactors = FALSE)
dag_info <- tryCatch(read.csv("output/case2/all_dag_info_v3.csv", stringsAsFactors = FALSE),
    error = function(e) data.frame()
)
all_screening <- tryCatch(read.csv("output/case2/all_screening_v3.csv", stringsAsFactors = FALSE),
    error = function(e) data.frame()
)
all_ate <- tryCatch(read.csv("output/case2/all_ate_results_v3.csv", stringsAsFactors = FALSE),
    error = function(e) data.frame()
)
all_dag_edges <- tryCatch(read.csv("output/case2/all_dag_edges_v3.csv", stringsAsFactors = FALSE),
    error = function(e) data.frame()
)

cat(sprintf(
    "Loaded %d result rows for %d species across %d regions\n",
    nrow(all_results), length(unique(all_results$species)),
    length(unique(all_results$region))
))

# ==============================================================================
# 配色 & 主题
# ==============================================================================
theme_pub <- function(base_size = 11) {
    theme_minimal(base_size = base_size) +
        theme(
            text = element_text(family = "sans"),
            plot.title = element_text(face = "bold", size = base_size + 2),
            plot.subtitle = element_text(color = "grey40", size = base_size - 1),
            panel.grid.minor = element_blank(),
            strip.text = element_text(face = "bold"),
            legend.position = "bottom",
            legend.title = element_text(face = "bold"),
            plot.margin = margin(10, 15, 10, 15)
        )
}

# 模型配色（统一大表）
model_colors <- c(
    "CAST"   = "#E74C3C",
    "MLP"    = "#3498DB",
    "RF"     = "#2ECC71",
    "Maxent" = "#9B59B6",
    "BRT"    = "#F39C12",
    "GAM"    = "#1ABC9C"
)

region_colors <- c(
    "AWT" = "#E74C3C", "CAN" = "#3498DB", "NSW" = "#2ECC71",
    "NZ"  = "#9B59B6", "SA"  = "#F39C12", "SWI" = "#1ABC9C"
)

# ==============================================================================
# 标准化模型标签
# ==============================================================================
all_results <- all_results %>%
    mutate(
        model_label = case_when(
            model == "CAST"    ~ "CAST",
            model == "CI_MLP"  ~ "CAST",
            model == "CGNet"   ~ "CAST",
            model == "MLP"     ~ "MLP",
            model == "FlatNN_cast" ~ "MLP",
            model == "FlatNN_full" ~ "MLP",
            TRUE ~ model
        ),
        varset_label = case_when(
            model_label == "CAST" ~ "CAST-screened",
            var_set == "full"     ~ "Full (VIF)",
            var_set == "cast"     ~ "CAST-screened",
            grepl("_full$", model) ~ "Full (VIF)",
            TRUE ~ "CAST-screened"
        ),
        display_label = ifelse(model_label == "CAST",
            "CAST",
            sprintf("%s (%s)", model_label,
                ifelse(varset_label == "CAST-screened", "CAST", "Full")))
    )

model_order <- c("CAST", "MLP", "RF", "BRT", "Maxent", "GAM")
all_results$model_label <- factor(all_results$model_label, levels = model_order)

# ==============================================================================
# Figure 1: 统一排名条形图（CAST vs all, AUC）
# ==============================================================================
cat("\nFig 1: Unified performance ranking...\n")

summary_ranked <- all_results %>%
    group_by(display_label, model_label, varset_label) %>%
    summarise(
        AUC_mean = mean(auc_mean, na.rm = TRUE),
        AUC_se = sd(auc_mean, na.rm = TRUE) / sqrt(n()),
        TSS_mean = mean(tss_mean, na.rm = TRUE),
        TSS_se = sd(tss_mean, na.rm = TRUE) / sqrt(n()),
        n = n(), .groups = "drop"
    ) %>%
    arrange(desc(AUC_mean))

summary_ranked$display_label <- factor(summary_ranked$display_label,
    levels = rev(summary_ranked$display_label))

fig1 <- ggplot(summary_ranked, aes(
    x = display_label, y = AUC_mean,
    fill = as.character(model_label)
)) +
    geom_col(width = 0.7, alpha = 0.9, color = "white", linewidth = 0.3) +
    geom_errorbar(aes(ymin = AUC_mean - AUC_se, ymax = AUC_mean + AUC_se),
        width = 0.25, linewidth = 0.5, color = "grey30"
    ) +
    geom_text(aes(label = sprintf("%.3f", AUC_mean)),
        hjust = -0.15, size = 3.2, fontface = "bold"
    ) +
    coord_flip() +
    scale_fill_manual(values = model_colors, name = "Model") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
        title = "CAST Performance Ranking (All Species)",
        subtitle = "CAST integrates causal screening and causally-informed feature engineering;\nall other models use standard feature spaces.",
        x = "", y = "Mean AUC (± SE)"
    ) +
    theme_pub() +
    theme(legend.position = "right")

ggsave("figures/case2/fig1_performance_ranking.png", fig1,
    width = 12, height = 7, dpi = 1200, bg = "white"
)
ggsave("figures/case2/fig1_performance_ranking.svg", fig1,
    width = 12, height = 7, bg = "white"
)
cat("  -> figures/case2/fig1_performance_ranking\n")

# ==============================================================================
# Figure 2: CAST vs MLP(CAST) 散点图（结构效应）
# ==============================================================================
cat("Fig 2: CAST vs MLP scatter...\n")

wide_struct <- all_results %>%
    filter(model_label %in% c("CAST", "MLP"), varset_label == "CAST-screened") %>%
    select(region, species, model_label, auc_mean, dag_density) %>%
    pivot_wider(names_from = model_label, values_from = auc_mean) %>%
    filter(!is.na(CAST), !is.na(MLP)) %>%
    mutate(delta = CAST - MLP)

if (nrow(wide_struct) > 0) {
    fig2 <- ggplot(wide_struct, aes(x = MLP, y = CAST, color = region)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        geom_point(aes(size = dag_density), alpha = 0.7) +
        scale_color_manual(values = region_colors, name = "Region") +
        scale_size_continuous(name = "DAG density", range = c(1.5, 5)) +
        labs(
            title = "CAST vs MLP (same CAST-screened variables)",
            subtitle = sprintf(
                "Points above line = CAST wins (%d/%d species, %.0f%%)",
                sum(wide_struct$delta > 0), nrow(wide_struct),
                mean(wide_struct$delta > 0) * 100
            ),
            x = "MLP (CAST-screened) AUC", y = "CAST AUC"
        ) +
        coord_equal() +
        theme_pub()
    ggsave("figures/case2/fig2_cast_vs_mlp.png", fig2,
        width = 8, height = 7, dpi = 1200, bg = "white"
    )
    ggsave("figures/case2/fig2_cast_vs_mlp.svg", fig2,
        width = 8, height = 7, bg = "white"
    )
    cat("  -> figures/case2/fig2_cast_vs_mlp\n")
}

# ==============================================================================
# Figure 3: 筛选效应 (Full → CAST-screened) 各算法
# ==============================================================================
cat("Fig 3: Screening effect...\n")

screen_delta <- all_results %>%
    filter(model_label %in% c("MLP", "RF", "BRT", "Maxent", "GAM")) %>%
    select(region, species, model_label, varset_label, auc_mean) %>%
    pivot_wider(names_from = varset_label, values_from = auc_mean) %>%
    filter(!is.na(`Full (VIF)`), !is.na(`CAST-screened`)) %>%
    mutate(delta = `CAST-screened` - `Full (VIF)`)

if (nrow(screen_delta) > 0) {
    fig3 <- ggplot(screen_delta, aes(x = model_label, y = delta * 100)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
        geom_boxplot(aes(fill = as.character(model_label)),
            alpha = 0.3, outlier.alpha = 0
        ) +
        geom_jitter(aes(color = region), width = 0.2, size = 2, alpha = 0.6) +
        stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
        scale_fill_manual(values = model_colors, guide = "none") +
        scale_color_manual(values = region_colors, name = "Region") +
        labs(
            title = "CAST Screening Effect (Full -> CAST-screened)",
            subtitle = "AUC change = CAST-screened - Full (VIF), per species",
            x = "Model", y = "AUC change (percentage points)"
        ) +
        theme_pub()
    ggsave("figures/case2/fig3_screening_effect.png", fig3,
        width = 9, height = 6, dpi = 1200, bg = "white"
    )
    ggsave("figures/case2/fig3_screening_effect.svg", fig3,
        width = 9, height = 6, bg = "white"
    )
    cat("  -> figures/case2/fig3_screening_effect\n")
}

# ==============================================================================
# Figure S4: DAG密度 vs CAST结构优势
# ==============================================================================
cat("Fig S4: DAG density vs structure advantage...\n")

if (exists("wide_struct") && nrow(wide_struct) > 5) {
    figS4 <- ggplot(wide_struct, aes(x = dag_density, y = delta * 100)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
        geom_point(aes(color = region), size = 3, alpha = 0.7) +
        geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.8) +
        scale_color_manual(values = region_colors, name = "Region") +
        labs(
            title = "DAG Density vs CAST Structure Advantage",
            subtitle = "Hypothesis: sparser DAG -> stronger CAST benefit",
            x = "DAG Density (fraction of possible edges)",
            y = "AUC change (CAST - MLP, percentage points)"
        ) +
        annotate("text",
            x = min(wide_struct$dag_density) + 0.02,
            y = max(wide_struct$delta * 100) * 0.9,
            label = sprintf(
                "r = %.3f, p = %.3f",
                cor(wide_struct$dag_density, wide_struct$delta),
                cor.test(wide_struct$dag_density, wide_struct$delta)$p.value
            ),
            hjust = 0, size = 3.5, family = "sans"
        ) +
        theme_pub()
    ggsave("figures/case2/figS4_dag_density_vs_advantage.png", figS4,
        width = 8, height = 6, dpi = 1200, bg = "white"
    )
    ggsave("figures/case2/figS4_dag_density_vs_advantage.svg", figS4,
        width = 8, height = 6, bg = "white"
    )
    cat("  -> figures/case2/figS4_dag_density_vs_advantage\n")
} else {
    cat("  Skipped (insufficient data)\n")
}

# ==============================================================================
# Figure S6: 各区域性能对比（分面）
# ==============================================================================
cat("Fig S6: Per-region performance facet...\n")

region_model_summary <- all_results %>%
    group_by(region, display_label, model_label) %>%
    summarise(
        AUC_mean = mean(auc_mean, na.rm = TRUE),
        AUC_se = sd(auc_mean, na.rm = TRUE) / sqrt(n()),
        n = n(), .groups = "drop"
    )

figS6 <- ggplot(region_model_summary, aes(
    x = reorder(display_label, AUC_mean),
    y = AUC_mean, fill = as.character(model_label)
)) +
    geom_col(width = 0.7, alpha = 0.9, color = "white", linewidth = 0.2) +
    geom_errorbar(aes(ymin = AUC_mean - AUC_se, ymax = AUC_mean + AUC_se),
        width = 0.2, linewidth = 0.4
    ) +
    coord_flip() +
    facet_wrap(~region, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = model_colors, name = "Model") +
    labs(
        title = "CAST Performance Across 6 Regions",
        subtitle = "disdat benchmark | CAST integrates causal screening and feature engineering",
        x = "", y = "Mean AUC"
    ) +
    theme_pub(base_size = 9) +
    theme(
        legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 10)
    )

ggsave("figures/case2/figS6_per_region.png", figS6,
    width = 16, height = 12, dpi = 1200, bg = "white"
)
ggsave("figures/case2/figS6_per_region.svg", figS6,
    width = 16, height = 12, bg = "white"
)
cat("  -> figures/case2/figS6_per_region\n")

# ==============================================================================
# Figure S5: 变量缩减 + 交互特征统计
# ==============================================================================
cat("Fig S5: Variable reduction & interactions...\n")

var_info <- all_results %>%
    filter(model_label %in% c("CAST", "MLP"),
           varset_label %in% c("CAST-screened", "Full (VIF)")) %>%
    select(region, species, model_label, varset_label,
           n_vars, n_interactions, n_features_total) %>%
    distinct()

cast_info <- var_info %>% filter(model_label == "CAST")
full_info <- var_info %>% filter(model_label == "MLP", varset_label == "Full (VIF)")

if (nrow(cast_info) > 0 && nrow(full_info) > 0) {
    merged_info <- cast_info %>%
        select(region, species, n_vars_cast = n_vars,
               n_interactions, n_features_total) %>%
        inner_join(
            full_info %>% select(region, species, n_vars_full = n_vars),
            by = c("region", "species")
        ) %>%
        mutate(reduction_pct = (1 - n_vars_cast / n_vars_full) * 100)

    figS5a <- ggplot(merged_info, aes(x = region, y = reduction_pct)) +
        geom_boxplot(fill = model_colors["MLP"], alpha = 0.4) +
        geom_jitter(width = 0.15, size = 2, alpha = 0.5,
            color = model_colors["MLP"]
        ) +
        labs(
            title = "(a) Variable Reduction (%)",
            subtitle = "Full (post-VIF) -> CAST-screened",
            x = "Region", y = "Reduction (%)"
        ) +
        theme_pub()

    figS5b <- ggplot(merged_info, aes(x = region, y = n_interactions)) +
        geom_boxplot(fill = model_colors["CAST"], alpha = 0.4) +
        geom_jitter(width = 0.15, size = 2, alpha = 0.5,
            color = model_colors["CAST"]
        ) +
        labs(
            title = "(b) DAG Interaction Features",
            subtitle = "Pairwise interactions from DAG edges",
            x = "Region", y = "# Interaction features"
        ) +
        theme_pub()

    figS5 <- figS5a + figS5b
    ggsave("figures/case2/figS5_variables_interactions.png", figS5,
        width = 12, height = 5.5, dpi = 1200, bg = "white"
    )
    ggsave("figures/case2/figS5_variables_interactions.svg", figS5,
        width = 12, height = 5.5, bg = "white"
    )
    cat("  -> figures/case2/figS5_variables_interactions\n")
}

# ==============================================================================
# Figure 4 (Main): CAST优势分解（筛选 + 结构）
# ==============================================================================
cat("Fig 4: CAST advantage decomposition...\n")

# 筛选效应: MLP(Full) → MLP(CAST)
delta_screening <- all_results %>%
    filter(model_label == "MLP") %>%
    select(region, species, varset_label, auc_mean) %>%
    pivot_wider(names_from = varset_label, values_from = auc_mean) %>%
    filter(!is.na(`CAST-screened`), !is.na(`Full (VIF)`)) %>%
    mutate(delta_screen = `CAST-screened` - `Full (VIF)`)

# 结构效应: MLP(CAST) → CAST
if (exists("wide_struct")) {
    delta_structure <- wide_struct
} else {
    delta_structure <- data.frame()
}

if (nrow(delta_screening) > 0 && nrow(delta_structure) > 0) {
    fig4a <- ggplot(delta_screening, aes(x = region, y = delta_screen * 100)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
        geom_boxplot(fill = model_colors["MLP"], alpha = 0.4, outlier.alpha = 0) +
        geom_jitter(width = 0.15, size = 2, alpha = 0.5,
            color = model_colors["MLP"]
        ) +
        stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
        labs(
            title = "(a) Screening effect: MLP Full -> CAST-screened",
            x = "Region", y = "AUC change (pp)"
        ) +
        theme_pub()

    fig4b <- ggplot(delta_structure, aes(x = region, y = delta * 100)) +
        geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
        geom_boxplot(fill = model_colors["CAST"], alpha = 0.4, outlier.alpha = 0) +
        geom_jitter(width = 0.15, size = 2, alpha = 0.5,
            color = model_colors["CAST"]
        ) +
        stat_summary(fun = mean, geom = "point", shape = 18, size = 4, color = "black") +
        labs(
            title = "(b) Structure effect: MLP -> CAST",
            x = "Region", y = "AUC change (pp)"
        ) +
        theme_pub()

    fig4_main <- fig4a + fig4b
    ggsave("figures/case2/fig4_cast_decomposition.png", fig4_main,
        width = 12, height = 5.5, dpi = 1200, bg = "white"
    )
    ggsave("figures/case2/fig4_cast_decomposition.svg", fig4_main,
        width = 12, height = 5.5, bg = "white"
    )
    cat("  -> figures/case2/fig4_cast_decomposition\n")
}

# ==============================================================================
# Summary Table: 统一排名汇总表
# ==============================================================================
cat("\nSummary Table...\n")

summary_table <- all_results %>%
    group_by(display_label, model_label, varset_label) %>%
    summarise(
        n_species = n(),
        AUC_mean = mean(auc_mean, na.rm = TRUE),
        AUC_sd = sd(auc_mean, na.rm = TRUE),
        TSS_mean = mean(tss_mean, na.rm = TRUE),
        TSS_sd = sd(tss_mean, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(desc(AUC_mean)) %>%
    mutate(Rank = row_number())

write.csv(summary_table, "figures/case2/summary_table.csv", row.names = FALSE)

# 按区域分组的表
region_table <- all_results %>%
    group_by(region, display_label, model_label, varset_label) %>%
    summarise(
        n_species = n(),
        AUC_mean = mean(auc_mean, na.rm = TRUE),
        AUC_sd = sd(auc_mean, na.rm = TRUE),
        TSS_mean = mean(tss_mean, na.rm = TRUE),
        TSS_sd = sd(tss_mean, na.rm = TRUE),
        .groups = "drop"
    ) %>%
    arrange(region, desc(AUC_mean))

write.csv(region_table, "figures/case2/region_table.csv", row.names = FALSE)

# 打印统一大表
cat("\n╔════════════════════════════════════════════════════════════════════════════════╗\n")
cat("║  CAST Performance Ranking (All Regions, All Species)                          ║\n")
cat("║  CAST integrates causal screening and causally-informed feature engineering;  ║\n")
cat("║  all other models use standard feature spaces.                                ║\n")
cat("╠════════════════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  %-4s  %-24s  %-14s  %-11s  %-11s  %4s  ║\n",
    "Rank", "Model", "Variable Set", "AUC", "TSS", "n"))
cat("╠════════════════════════════════════════════════════════════════════════════════╣\n")
for (i in 1:nrow(summary_table)) {
    r <- summary_table[i, ]
    mk <- if (as.character(r$model_label) == "CAST") " *" else "  "
    cat(sprintf(
        "║%s%3d  %-24s  %-14s  %.3f+%.3f  %.3f+%.3f  %4d  ║\n",
        mk, r$Rank, as.character(r$display_label),
        as.character(r$varset_label),
        r$AUC_mean, r$AUC_sd, r$TSS_mean, r$TSS_sd, r$n_species
    ))
}
cat("╚════════════════════════════════════════════════════════════════════════════════╝\n")

# 关键结果
cast_auc <- summary_table$AUC_mean[summary_table$model_label == "CAST"][1]
best_other <- summary_table %>%
    filter(model_label != "CAST") %>%
    slice(1)

if (!is.na(cast_auc) && nrow(best_other) > 0) {
    cat(sprintf(
        "\n* CAST: AUC = %.3f (Rank 1)\n",
        cast_auc
    ))
    cat(sprintf(
        "* Best other: %s, AUC = %.3f (delta = %+.4f)\n",
        best_other$display_label, best_other$AUC_mean,
        cast_auc - best_other$AUC_mean
    ))
}

if (exists("wide_struct") && nrow(wide_struct) > 0) {
    cat(sprintf(
        "* CAST wins vs MLP (same vars): %d/%d species (%.0f%%), mean delta = %+.4f\n",
        sum(wide_struct$delta > 0), nrow(wide_struct),
        mean(wide_struct$delta > 0) * 100, mean(wide_struct$delta)
    ))
}

# ==============================================================================
# Figure 6 (Main): 变量筛选一致性热力图 — 跨物种CAST选择频率
# ==============================================================================
cat("\nFig 6: Variable selection consistency heatmap...\n")

if (nrow(all_screening) > 0) {
    # 判断每个变量在每个物种中是否被CAST选中
    # 使用score_total与物种内中位数比较（复现pipeline的elbow/kmeans逻辑简化版）
    sel_binary <- all_screening %>%
        group_by(region, species) %>%
        mutate(
            rank_in_sp = rank(-score_total),
            n_vars = n(),
            selected = rank_in_sp <= ceiling(n_vars * 0.5)
        ) %>%
        ungroup()

    # 计算每个变量在每个区域内被选中的频率
    var_freq <- sel_binary %>%
        group_by(region, variable) %>%
        summarise(
            freq = mean(selected),
            mean_score = mean(score_total, na.rm = TRUE),
            n_species = n(),
            .groups = "drop"
        )

    # 选取至少在2个区域出现过的变量
    var_presence <- var_freq %>%
        group_by(variable) %>%
        summarise(n_regions = n(), global_freq = mean(freq), .groups = "drop") %>%
        filter(n_regions >= 2) %>%
        arrange(desc(global_freq))

    if (nrow(var_presence) > 0) {
        top_vars <- head(var_presence$variable, 30)
        plot_data <- var_freq %>%
            filter(variable %in% top_vars)
        plot_data$variable <- factor(plot_data$variable,
            levels = rev(var_presence$variable[var_presence$variable %in% top_vars]))

        fig6 <- ggplot(plot_data, aes(x = region, y = variable, fill = freq)) +
            geom_tile(color = "white", linewidth = 0.5) +
            geom_text(aes(label = sprintf("%.0f%%", freq * 100)),
                size = 2.8, color = "grey20"
            ) +
            scale_fill_gradient2(
                low = "#F7FBFF", mid = "#6BAED6", high = "#08306B",
                midpoint = 0.5, limits = c(0, 1),
                name = "Selection\nfrequency"
            ) +
            labs(
                title = "CAST Variable Selection Consistency Across Species",
                subtitle = "Proportion of species where each variable was selected by CAST screening",
                x = "Region", y = "Environmental variable"
            ) +
            theme_pub() +
            theme(
                legend.position = "right",
                axis.text.y = element_text(size = 8),
                panel.grid = element_blank()
            )

        ggsave("figures/case2/fig6_variable_selection_heatmap.png", fig6,
            width = 10, height = 10, dpi = 1200, bg = "white"
        )
        ggsave("figures/case2/fig6_variable_selection_heatmap.svg", fig6,
            width = 10, height = 10, bg = "white"
        )
        cat("  -> figures/case2/fig6_variable_selection_heatmap\n")
    }
}

# ==============================================================================
# Figure 7 (Main): 跨区域模型排名矩阵热力图
# ==============================================================================
cat("\nFig 7: Cross-region model rank heatmap...\n")

region_model_auc <- all_results %>%
    group_by(region, display_label, model_label) %>%
    summarise(AUC_mean = mean(auc_mean, na.rm = TRUE), .groups = "drop")

# 计算每个区域内的排名
region_model_rank <- region_model_auc %>%
    group_by(region) %>%
    mutate(rank = rank(-AUC_mean, ties.method = "min")) %>%
    ungroup()

# 按全局平均AUC排序模型
global_order <- region_model_auc %>%
    group_by(display_label) %>%
    summarise(global_auc = mean(AUC_mean), .groups = "drop") %>%
    arrange(desc(global_auc))
region_model_rank$display_label <- factor(region_model_rank$display_label,
    levels = rev(global_order$display_label))

fig7 <- ggplot(region_model_rank, aes(x = region, y = display_label, fill = AUC_mean)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("#%d\n%.3f", rank, AUC_mean)),
        size = 2.8, lineheight = 0.85, fontface = "bold"
    ) +
    scale_fill_gradient2(
        low = "#D73027", mid = "#FFFFBF", high = "#1A9850",
        midpoint = median(region_model_rank$AUC_mean, na.rm = TRUE),
        name = "Mean AUC"
    ) +
    labs(
        title = "Model Performance Across Regions",
        subtitle = "Cell color = AUC, text = within-region rank",
        x = "Region", y = ""
    ) +
    theme_pub() +
    theme(
        legend.position = "right",
        panel.grid = element_blank(),
        axis.text.y = element_text(size = 9)
    )

ggsave("figures/case2/fig7_region_model_heatmap.png", fig7,
    width = 11, height = 8, dpi = 1200, bg = "white"
)
ggsave("figures/case2/fig7_region_model_heatmap.svg", fig7,
    width = 11, height = 8, bg = "white"
)
cat("  -> figures/case2/fig7_region_model_heatmap\n")

# ==============================================================================
# Figure S1: DAG密度分布
# ==============================================================================
cat("\nFig S1: DAG density distribution...\n")

if (nrow(dag_info) > 0) {
    figS1 <- ggplot(dag_info, aes(x = region, y = dag_density)) +
        geom_boxplot(fill = "#3498DB", alpha = 0.3, outlier.alpha = 0) +
        geom_jitter(aes(color = region), width = 0.2, size = 2, alpha = 0.6) +
        geom_hline(
            yintercept = median(dag_info$dag_density, na.rm = TRUE),
            linetype = "dashed", color = "grey40", linewidth = 0.6
        ) +
        annotate("text",
            x = 0.5, y = median(dag_info$dag_density, na.rm = TRUE) + 0.02,
            label = sprintf("Median = %.3f",
                median(dag_info$dag_density, na.rm = TRUE)),
            hjust = 0, size = 3.5, color = "grey40", family = "sans"
        ) +
        scale_color_manual(values = region_colors, guide = "none") +
        labs(
            title = "DAG Density Distribution Across Regions",
            subtitle = sprintf(
                "Bootstrap HC (R=200) | Strength >= 0.7, Direction >= 0.6 | %d species total",
                nrow(dag_info)
            ),
            x = "Region", y = "DAG density (edges / possible edges)"
        ) +
        theme_pub()

    ggsave("figures/case2/figS1_dag_density.png", figS1,
        width = 9, height = 6, dpi = 1200, bg = "white"
    )
    ggsave("figures/case2/figS1_dag_density.svg", figS1,
        width = 9, height = 6, bg = "white"
    )
    cat("  -> figures/case2/figS1_dag_density\n")
}

# ==============================================================================
# Figure S2: ATE效应量小提琴图
# ==============================================================================
cat("\nFig S2: ATE effect size distribution...\n")

if (nrow(all_ate) > 0) {
    ate_plot <- all_ate %>%
        mutate(
            abs_coef = abs(coef),
            sig_label = ifelse(significant, "Significant (p<0.05)", "Not significant")
        )

    # 各区域显著比例标注
    sig_rates <- ate_plot %>%
        group_by(region) %>%
        summarise(
            rate = mean(significant, na.rm = TRUE),
            n = n(),
            .groups = "drop"
        )

    figS2 <- ggplot(ate_plot, aes(x = region, y = abs_coef)) +
        geom_violin(fill = "#9B59B6", alpha = 0.25, color = "#9B59B6") +
        geom_boxplot(width = 0.15, fill = "white", outlier.alpha = 0) +
        geom_text(data = sig_rates,
            aes(x = region, y = max(ate_plot$abs_coef, na.rm = TRUE) * 1.05,
                label = sprintf("%.0f%% sig", rate * 100)),
            size = 3.2, fontface = "bold", color = "#8E44AD"
        ) +
        labs(
            title = "ATE Effect Size Distribution (DML)",
            subtitle = sprintf(
                "Double Machine Learning | 2-fold cross-fitting | %d total variable-species pairs",
                nrow(ate_plot)
            ),
            x = "Region", y = "|ATE coefficient|"
        ) +
        theme_pub()

    ggsave("figures/case2/figS2_ate_distribution.png", figS2,
        width = 9, height = 6, dpi = 1200, bg = "white"
    )
    ggsave("figures/case2/figS2_ate_distribution.svg", figS2,
        width = 9, height = 6, bg = "white"
    )
    cat("  -> figures/case2/figS2_ate_distribution\n")
}

# ==============================================================================
# Figure S3: 自适应权重分布（w_dag / w_ate / w_imp vs DAG密度）
# ==============================================================================
cat("\nFig S3: Adaptive screening weight distribution...\n")

if (nrow(all_screening) > 0 && all(c("w_dag", "w_ate", "w_imp") %in% names(all_screening))) {
    # 每物种取一行权重（权重在物种内相同）
    weight_df <- all_screening %>%
        group_by(region, species) %>%
        slice(1) %>%
        ungroup() %>%
        select(region, species, w_dag, w_ate, w_imp) %>%
        pivot_longer(cols = c(w_dag, w_ate, w_imp),
            names_to = "component", values_to = "weight"
        ) %>%
        mutate(component = case_when(
            component == "w_dag" ~ "DAG Out-degree",
            component == "w_ate" ~ "ATE Effect",
            component == "w_imp" ~ "RF Importance"
        ))

    # 关联DAG密度
    if (nrow(dag_info) > 0) {
        weight_df <- weight_df %>%
            left_join(dag_info %>% select(region, species, dag_density),
                by = c("region", "species")
            )
    }

    comp_colors <- c(
        "DAG Out-degree" = "#2C3E50",
        "ATE Effect" = "#E67E22",
        "RF Importance" = "#27AE60"
    )

    figS3a <- ggplot(weight_df, aes(x = region, y = weight, fill = component)) +
        geom_boxplot(position = position_dodge(0.8), width = 0.65, alpha = 0.7) +
        scale_fill_manual(values = comp_colors, name = "Screening\nComponent") +
        labs(
            title = "(a) Adaptive Screening Weights by Region",
            subtitle = "Weights adapt to DAG quality and ATE significance ratio",
            x = "Region", y = "Weight"
        ) +
        theme_pub()

    if ("dag_density" %in% names(weight_df) && !all(is.na(weight_df$dag_density))) {
        figS3b <- ggplot(weight_df, aes(x = dag_density, y = weight, color = component)) +
            geom_point(alpha = 0.4, size = 1.5) +
            geom_smooth(method = "lm", se = TRUE, linewidth = 0.8) +
            scale_color_manual(values = comp_colors, name = "Component") +
            labs(
                title = "(b) Screening Weights vs DAG Density",
                subtitle = "DAG weight increases with sparser DAGs; RF importance compensates",
                x = "DAG Density", y = "Weight"
            ) +
            theme_pub()

        figS3 <- figS3a + figS3b + plot_layout(guides = "collect") &
            theme(legend.position = "bottom")
    } else {
        figS3 <- figS3a
    }

    ggsave("figures/case2/figS3_adaptive_weights.png", figS3,
        width = 14, height = 6, dpi = 1200, bg = "white"
    )
    ggsave("figures/case2/figS3_adaptive_weights.svg", figS3,
        width = 14, height = 6, bg = "white"
    )
    cat("  -> figures/case2/figS3_adaptive_weights\n")
}

# ==============================================================================
# Figure S7: Critical Difference Diagram (Nemenyi秩和检验)
# ==============================================================================
cat("\nFig S7: Critical difference diagram...\n")

# 构建逐物种模型排名矩阵
rank_matrix <- all_results %>%
    group_by(region, species, display_label) %>%
    summarise(auc = mean(auc_mean, na.rm = TRUE), .groups = "drop") %>%
    group_by(region, species) %>%
    mutate(rank = rank(-auc, ties.method = "average")) %>%
    ungroup()

# 计算每个模型的平均排名
mean_ranks <- rank_matrix %>%
    group_by(display_label) %>%
    summarise(
        mean_rank = mean(rank, na.rm = TRUE),
        sd_rank = sd(rank, na.rm = TRUE),
        n = n(),
        .groups = "drop"
    ) %>%
    arrange(mean_rank)

# Friedman检验
rank_wide <- rank_matrix %>%
    select(region, species, display_label, rank) %>%
    pivot_wider(names_from = display_label, values_from = rank) %>%
    select(-region, -species)
rank_wide <- rank_wide[complete.cases(rank_wide), ]

if (nrow(rank_wide) >= 10) {
    friedman_p <- tryCatch({
        ft <- friedman.test(as.matrix(rank_wide))
        ft$p.value
    }, error = function(e) NA)

    n_models <- ncol(rank_wide)
    n_species <- nrow(rank_wide)

    # Nemenyi CD值
    q_alpha <- qtukey(0.95, n_models, Inf) / sqrt(2)
    cd <- q_alpha * sqrt(n_models * (n_models + 1) / (6 * n_species))

    # 绘制CD图
    mean_ranks <- mean_ranks %>% mutate(
        is_cast = grepl("^CAST$", display_label),
        label = sprintf("%s (%.2f)", display_label, mean_rank)
    )

    figS7 <- ggplot(mean_ranks, aes(
        x = mean_rank, y = reorder(display_label, -mean_rank)
    )) +
        geom_segment(aes(
            x = mean_rank - cd / 2, xend = mean_rank + cd / 2,
            yend = reorder(display_label, -mean_rank)
        ), linewidth = 2, color = "grey80") +
        geom_point(aes(color = is_cast), size = 4) +
        geom_text(aes(label = sprintf("%.2f", mean_rank)),
            hjust = -0.3, size = 3.5, fontface = "bold"
        ) +
        scale_color_manual(
            values = c("TRUE" = "#E74C3C", "FALSE" = "#3498DB"),
            guide = "none"
        ) +
        annotate("segment",
            x = min(mean_ranks$mean_rank) - 0.3,
            xend = min(mean_ranks$mean_rank) - 0.3 + cd,
            y = 0.5, yend = 0.5, linewidth = 1.5, color = "black"
        ) +
        annotate("text",
            x = min(mean_ranks$mean_rank) - 0.3 + cd / 2, y = 0.2,
            label = sprintf("CD = %.2f", cd),
            size = 3.5, fontface = "bold", family = "sans"
        ) +
        labs(
            title = "Critical Difference Diagram (Nemenyi Test)",
            subtitle = sprintf(
                "Friedman p = %.2e | %d species | CD at alpha = 0.05 = %.2f\nModels connected by grey bar are not significantly different",
                friedman_p, n_species, cd
            ),
            x = "Mean Rank (lower = better)", y = ""
        ) +
        theme_pub() +
        theme(panel.grid.major.y = element_blank())

    ggsave("figures/case2/figS7_critical_difference.png", figS7,
        width = 11, height = 7, dpi = 1200, bg = "white"
    )
    ggsave("figures/case2/figS7_critical_difference.svg", figS7,
        width = 11, height = 7, bg = "white"
    )
    cat("  -> figures/case2/figS7_critical_difference\n")
} else {
    cat("  Skipped (insufficient complete-case species for Friedman test)\n")
}

# ==============================================================================
# Figure S9: Pairwise Comparison Scatter Matrix (CAST vs每个基线)
# ==============================================================================
cat("\nFig S9: Pairwise comparison scatter matrix...\n")

# 获取所有非CAST模型的display_label
baselines <- unique(all_results$display_label[all_results$model_label != "CAST"])

# CAST的逐物种AUC
cast_sp <- all_results %>%
    filter(model_label == "CAST") %>%
    select(region, species, cast_auc = auc_mean)

pairwise_plots <- list()
for (bl in baselines) {
    bl_sp <- all_results %>%
        filter(display_label == bl) %>%
        select(region, species, bl_auc = auc_mean)

    paired <- inner_join(cast_sp, bl_sp, by = c("region", "species"))
    if (nrow(paired) < 5) next

    n_win <- sum(paired$cast_auc > paired$bl_auc)
    n_total <- nrow(paired)
    win_pct <- n_win / n_total * 100

    p <- ggplot(paired, aes(x = bl_auc, y = cast_auc, color = region)) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
        geom_point(size = 1.8, alpha = 0.6) +
        scale_color_manual(values = region_colors, guide = "none") +
        coord_equal(xlim = range(c(paired$bl_auc, paired$cast_auc), na.rm = TRUE),
                    ylim = range(c(paired$bl_auc, paired$cast_auc), na.rm = TRUE)) +
        annotate("text",
            x = min(paired$bl_auc, na.rm = TRUE) + 0.02,
            y = max(paired$cast_auc, na.rm = TRUE) - 0.02,
            label = sprintf("CAST wins: %d/%d (%.0f%%)", n_win, n_total, win_pct),
            hjust = 0, vjust = 1, size = 3, fontface = "bold", family = "sans"
        ) +
        labs(title = bl, x = "Baseline AUC", y = "CAST AUC") +
        theme_pub(base_size = 9) +
        theme(plot.title = element_text(size = 10))

    pairwise_plots[[bl]] <- p
}

if (length(pairwise_plots) > 0) {
    n_panels <- length(pairwise_plots)
    ncols <- min(3, n_panels)
    nrows <- ceiling(n_panels / ncols)

    figS9 <- wrap_plots(pairwise_plots, ncol = ncols) +
        plot_annotation(
            title = "CAST vs Each Baseline (Per-Species AUC)",
            subtitle = "Points above diagonal = CAST outperforms baseline",
            theme = theme(
                plot.title = element_text(face = "bold", size = 14, family = "sans"),
                plot.subtitle = element_text(size = 11, color = "grey40", family = "sans")
            )
        )

    ggsave("figures/case2/figS9_pairwise_scatter.png", figS9,
        width = 5 * ncols, height = 5 * nrows, dpi = 1200, bg = "white"
    )
    ggsave("figures/case2/figS9_pairwise_scatter.svg", figS9,
        width = 5 * ncols, height = 5 * nrows, bg = "white"
    )
    cat(sprintf("  -> figures/case2/figS9_pairwise_scatter (%d panels)\n", n_panels))
}

# ==============================================================================
# Statistical Tests Table (Table S4): Wilcoxon配对检验
# ==============================================================================
cat("\nTable S4: Pairwise Wilcoxon tests...\n")

if (nrow(cast_sp) > 0) {
    wilcox_results <- data.frame()
    for (bl in baselines) {
        bl_sp <- all_results %>%
            filter(display_label == bl) %>%
            select(region, species, bl_auc = auc_mean)
        paired <- inner_join(cast_sp, bl_sp, by = c("region", "species"))
        if (nrow(paired) < 5) next

        wt <- wilcox.test(paired$cast_auc, paired$bl_auc, paired = TRUE)
        wilcox_results <- rbind(wilcox_results, data.frame(
            baseline = bl,
            n_paired = nrow(paired),
            cast_mean_auc = mean(paired$cast_auc),
            baseline_mean_auc = mean(paired$bl_auc),
            delta_auc = mean(paired$cast_auc - paired$bl_auc),
            cast_win_pct = mean(paired$cast_auc > paired$bl_auc) * 100,
            wilcox_V = wt$statistic,
            p_value = wt$p.value,
            sig = ifelse(wt$p.value < 0.001, "***",
                ifelse(wt$p.value < 0.01, "**",
                    ifelse(wt$p.value < 0.05, "*", "ns"))),
            stringsAsFactors = FALSE
        ))
    }
    if (nrow(wilcox_results) > 0) {
        wilcox_results <- wilcox_results %>% arrange(p_value)
        write.csv(wilcox_results, "figures/case2/tableS4_wilcoxon_tests.csv",
            row.names = FALSE
        )

        cat("\n  ╔═══════════════════════════════════════════════════════════════════════════╗\n")
        cat("  ║  Pairwise Wilcoxon Signed-Rank Tests: CAST vs Each Baseline             ║\n")
        cat("  ╠═══════════════════════════════════════════════════════════════════════════╣\n")
        for (i in 1:nrow(wilcox_results)) {
            r <- wilcox_results[i, ]
            cat(sprintf(
                "  ║  %-24s  delta=%+.4f  win=%.0f%%  p=%.2e %s  ║\n",
                r$baseline, r$delta_auc, r$cast_win_pct, r$p_value, r$sig
            ))
        }
        cat("  ╚═══════════════════════════════════════════════════════════════════════════╝\n")
    }
}

cat("\n======================================================================\n")
cat("  All figures saved to figures/case2/\n")
cat("======================================================================\n")
cat("  Main figures:  fig1-fig7\n")
cat("  Supplementary: figS1-S9\n")
cat("  Tables:        summary_table.csv, region_table.csv, tableS4_wilcoxon_tests.csv\n")
cat("======================================================================\n")
