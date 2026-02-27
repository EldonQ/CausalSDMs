# =============================================================================
# Fig S6: 全部物种 AUC-TSS 双指标散点图（逐模型×区域）
# 数据来源: results/case2/all_results_v3.csv
# 输出: figures/case2/plot/figS6_auc_tss_scatter_all_species.{png,svg}
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(scales)
    library(patchwork)
})

# ── 路径配置 ──────────────────────────────────────────────────────────────────
script_dir <- tryCatch(
    dirname(rstudioapi::getSourceEditorContext()$path),
    error = function(e) getwd()
)
proj_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
data_dir  <- file.path(proj_root, "results", "case2")
fig_dir   <- file.path(proj_root, "figures", "case2", "plot")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

region_colors <- c(
    "AWT" = "#E64B35", "CAN" = "#4DBBD5", "NSW" = "#00A087",
    "NZ"  = "#3C5488", "SA"  = "#F39B7F", "SWI" = "#8491B4"
)
model_order <- c("CAST", "MLP_ATE", "MLP", "RF", "BRT", "GAM", "MaxEnt")

# ── 读取数据 ─────────────────────────────────────────────────────────────────
res_file <- file.path(data_dir, "all_results_v3.csv")
if (!file.exists(res_file)) stop("数据文件不存在: ", res_file)
results <- read.csv(res_file, stringsAsFactors = FALSE) %>%
    filter(model %in% model_order, !is.na(auc_mean), !is.na(tss_mean)) %>%
    mutate(model = factor(model, levels = model_order),
           is_cast = model == "CAST")

# 各模型相关系数
cor_df <- results %>%
    group_by(model) %>%
    summarise(r = cor(auc_mean, tss_mean, use = "complete.obs"),
              n = n(), .groups = "drop") %>%
    mutate(label = sprintf("r = %.2f (n=%d)", r, n))

# ── Panel A: 全数据散点（模型颜色，区域形状） ─────────────────────────────────
model_colors <- c(
    "CAST"    = "#E64B35", "MLP_ATE" = "#F39B7F", "MLP"     = "#FBBC05",
    "RF"      = "#4DBBD5", "BRT"     = "#3C5488",  "GAM"     = "#00A087",
    "MaxEnt"  = "#8491B4"
)

pA <- ggplot(results,
             aes(x = auc_mean, y = tss_mean, color = model, shape = region)) +
    geom_abline(slope = 1, intercept = -0.5, linetype = "dotted",
                color = "grey60", linewidth = 0.5) +
    geom_point(alpha = 0.55, size = 1.8) +
    # 参考性能阈值
    geom_hline(yintercept = 0.4, linetype = "dashed", color = "grey50",
               linewidth = 0.5) +
    geom_vline(xintercept = 0.7, linetype = "dashed", color = "grey50",
               linewidth = 0.5) +
    annotate("text", x = 0.705, y = 0.02, label = "AUC=0.7",
             hjust = 0, size = 2.8, color = "grey45", family = "Arial") +
    annotate("text", x = 0.42, y = 0.42, label = "TSS=0.4",
             hjust = 0, size = 2.8, color = "grey45", family = "Arial") +
    scale_color_manual(values = model_colors, name = "Model") +
    scale_shape_manual(values = c(16, 17, 15, 18, 8, 7),
                       name = "Region") +
    labs(
        title    = "A  AUC vs TSS for All Models and Species",
        x        = "AUC (mean over runs)",
        y        = "TSS (mean over runs)",
        subtitle = "Dashed lines mark moderate performance thresholds (AUC=0.7, TSS=0.4)"
    ) +
    coord_fixed(ratio = 0.6) +
    theme_classic(base_family = "Arial", base_size = 10) +
    theme(
        plot.title       = element_text(face = "bold", size = 11),
        plot.subtitle    = element_text(size = 8, color = "grey40"),
        legend.title     = element_text(face = "bold", size = 9),
        legend.text      = element_text(size = 8),
        legend.key.size  = unit(0.4, "cm"),
        panel.grid.major = element_line(color = "grey92", linewidth = 0.35)
    )

# ── Panel B: 逐模型 facet（区域着色）─────────────────────────────────────────
pB <- ggplot(results,
             aes(x = auc_mean, y = tss_mean, color = region)) +
    geom_abline(slope = 1, intercept = -0.5, linetype = "dotted",
                color = "grey70", linewidth = 0.4) +
    geom_point(alpha = 0.55, size = 1.3) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.7, alpha = 0.7) +
    geom_text(
        data  = cor_df,
        aes(x = 0.52, y = 0.92, label = label),
        hjust = 0, size = 2.5, color = "grey30", family = "Arial",
        inherit.aes = FALSE
    ) +
    scale_color_manual(values = region_colors, name = "Region") +
    facet_wrap(~model, ncol = 4, scales = "fixed") +
    labs(
        title    = "B  AUC vs TSS by Model (colored by Region)",
        x        = "AUC",
        y        = "TSS",
        subtitle = "Lines = OLS fit per region; r = Pearson correlation for all species"
    ) +
    theme_classic(base_family = "Arial", base_size = 9) +
    theme(
        strip.background   = element_blank(),
        strip.text         = element_text(face = "bold", size = 9),
        legend.position    = "right",
        legend.title       = element_text(face = "bold", size = 9),
        plot.title         = element_text(face = "bold", size = 11),
        plot.subtitle      = element_text(size = 8, color = "grey40"),
        panel.grid.major   = element_line(color = "grey92", linewidth = 0.3),
        legend.key.size    = unit(0.4, "cm")
    )

# ── 合并 ─────────────────────────────────────────────────────────────────────
p_combined <- pA / pB + plot_layout(heights = c(1, 1.6)) +
    plot_annotation(
        title    = "Dual-Metric Performance: AUC vs TSS",
        subtitle = "All species across 6 biogeographic regions",
        theme    = theme(
            plot.title    = element_text(face = "bold", size = 14, family = "Arial"),
            plot.subtitle = element_text(size = 9, color = "grey40", family = "Arial")
        )
    )

# ── 保存 ─────────────────────────────────────────────────────────────────────
out_prefix <- file.path(fig_dir, "figS6_auc_tss_scatter_all_species")
ggsave(paste0(out_prefix, ".png"), p_combined, width = 13, height = 14,
       dpi = 1200, bg = "white")
ggsave(paste0(out_prefix, ".svg"), p_combined, width = 13, height = 14, bg = "white")
cat("Fig S6 saved:", fig_dir, "\n")
