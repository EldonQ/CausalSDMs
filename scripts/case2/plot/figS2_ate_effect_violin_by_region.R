# =============================================================================
# Fig S2: ATE效应量逐区域Violin图
# 数据来源: results/case2/all_ate_results_v3.csv
# 输出: figures/case2/plot/figS2_ate_effect_violin_by_region.{png,svg}
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(ggforce)  # geom_sina
    library(scales)
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

# ── 读取数据 ─────────────────────────────────────────────────────────────────
ate_file <- file.path(data_dir, "all_ate_results_v3.csv")
if (!file.exists(ate_file)) {
    stop("数据文件不存在，请先运行 02: ", ate_file)
}
ate <- read.csv(ate_file, stringsAsFactors = FALSE)

# 取ATE绝对值，标注显著性
ate <- ate %>%
    mutate(
        abs_ate   = abs(coef),
        sig_label = ifelse(significant, "Significant", "Non-significant")
    ) %>%
    filter(!is.na(abs_ate))

# 每个区域-变量的中位绝对ATE，用于排序变量
var_order <- ate %>%
    group_by(variable) %>%
    summarise(med_abs = median(abs_ate, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(med_abs)) %>%
    pull(variable)

ate$variable <- factor(ate$variable, levels = var_order)

# 每个区域中显著比例（用于facet标题补充）
sig_rate <- ate %>%
    group_by(region) %>%
    summarise(
        pct_sig = round(100 * mean(significant, na.rm = TRUE), 1),
        n       = n(),
        .groups = "drop"
    )
ate <- ate %>%
    left_join(
        sig_rate %>% mutate(region_label = sprintf("%s\n(%s%% sig.)", region, pct_sig)),
        by = "region"
    )

# ── 绘图 ─────────────────────────────────────────────────────────────────────
p <- ggplot(ate, aes(x = abs_ate, y = variable, fill = region, color = region)) +
    geom_violin(
        scale = "width", alpha = 0.35, linewidth = 0.4,
        trim = TRUE, bw = "nrd0"
    ) +
    # 箱线摘要
    stat_summary(
        fun.data = function(x) {
            data.frame(
                x    = median(x, na.rm = TRUE),
                xmin = quantile(x, 0.25, na.rm = TRUE),
                xmax = quantile(x, 0.75, na.rm = TRUE)
            )
        },
        geom = "crossbar", width = 0.45, linewidth = 0.5,
        aes(color = region)
    ) +
    # 零效应参考线
    geom_vline(xintercept = 0, linetype = "dotted", color = "grey50", linewidth = 0.5) +
    scale_fill_manual(values = region_colors)  +
    scale_color_manual(values = region_colors) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
    facet_wrap(~region_label, ncol = 3, scales = "free_x") +
    labs(
        title    = "ATE Effect Size Distribution by Region and Variable",
        subtitle = "Absolute DML-estimated Average Treatment Effect | Crossbar = IQR, center = median",
        x        = "|ATE| (Absolute Average Treatment Effect)",
        y        = "Environmental Variable"
    ) +
    theme_classic(base_family = "Arial", base_size = 10) +
    theme(
        strip.background = element_blank(),
        strip.text       = element_text(face = "bold", size = 9),
        legend.position  = "none",
        plot.title       = element_text(face = "bold", size = 12),
        plot.subtitle    = element_text(size = 8.5, color = "grey40"),
        axis.text.y      = element_text(size = 8),
        axis.title       = element_text(size = 10),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.35)
    )

# ── 保存 ─────────────────────────────────────────────────────────────────────
out_prefix <- file.path(fig_dir, "figS2_ate_effect_violin_by_region")
ggsave(paste0(out_prefix, ".png"), p, width = 12, height = 10,
       dpi = 1200, bg = "white")
ggsave(paste0(out_prefix, ".svg"), p, width = 12, height = 10, bg = "white")
cat("Fig S2 saved:", fig_dir, "\n")
