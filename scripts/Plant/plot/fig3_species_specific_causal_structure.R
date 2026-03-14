################################################################################
# Fig 3: CAST 挖掘物种特异性因果结构 — 核心证据图
#
# 核心叙事: 同一套 CAST 流程在不同物种上自动挖掘出不同的因果结构、
#           因果效应方向、因果角色和变量筛选优先级
#
# 四个独立子图（每图单独保存 1200 dpi PNG + SVG）:
#   (a) DAG 边频谱: 少量通用边 vs 大量物种特有边
#   (b) ATE 系数热力图: 同一变量在不同物种中正/负方向截然不同
#   (c) 因果角色热力图: 同一变量在不同物种中扮演不同因果角色
#   (d) 自适应筛选得分热力图: 不同物种对变量重视程度不同
#
# 数据来源 (Plant 案例, 03_run_Plant_multi_species.R 输出):
#   output/case4_plant/all_dag_edges_plant.csv
#   output/case4_plant/all_ate_results_plant.csv
#   output/case4_plant/all_role_info_plant.csv
#   output/case4_plant/all_screening_plant.csv
#
# 运行: setwd("E:/CausalSDMs")
#       source("scripts/Plant/plot/fig3_species_specific_causal_structure.R")
################################################################################

rm(list = ls())
setwd("E:/CausalSDMs")

fig_dir <- "figures/case4_plant/plot"
tbl_dir <- "figures/case4_plant/tables"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tbl_dir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
    library(tidyverse)
})

# ══════════════════════════════════════════════════════════════════════════════
# 全局主题 · Nature 风格
# ══════════════════════════════════════════════════════════════════════════════
theme_pub <- function(base_size = 10) {
    theme_minimal(base_size = base_size, base_family = "sans") +
        theme(
            panel.grid       = element_blank(),
            axis.title       = element_text(face = "bold", size = 10),
            axis.text        = element_text(size = 8, color = "grey20"),
            plot.title       = element_text(face = "bold", size = 11, hjust = 0,
                                            margin = margin(b = 3)),
            plot.subtitle    = element_text(size = 8.5, color = "grey40", hjust = 0,
                                            margin = margin(b = 8)),
            legend.title     = element_text(face = "bold", size = 9),
            legend.text      = element_text(size = 8),
            plot.margin      = margin(12, 12, 8, 8)
        )
}

# ══════════════════════════════════════════════════════════════════════════════
# 变量 → 英文显示名称映射 (含 Plant 案例变量)
# ══════════════════════════════════════════════════════════════════════════════
var_labels <- c(
    "aridityindexthornthwaite" = "Aridity Index",
    "bio02" = "Diurnal Range", "bio_2" = "Mean Diurnal Range",
    "bio15" = "Precip. Seasonality", "bio_15" = "Precip. Seasonality",
    "bio19" = "Precip. Coldest Qtr", "bio_19" = "Precip. Coldest Qtr",
    "bio03" = "Isothermality", "bio_3" = "Isothermality",
    "bio18" = "Precip. Warmest Qtr", "bio_18" = "Precip. Warmest Qtr",
    "elevation" = "Elevation", "Elevation" = "Elevation",
    "etccdi_cwd" = "Consecutive Wet Days",
    "landcover_igbp" = "Land Cover (IGBP)",
    "maxtempcoldest" = "Tmax Coldest Month",
    "nontree" = "Non-tree Vegetation",
    "topowet" = "Topographic Wetness",
    "tri" = "Terrain Ruggedness",
    "Slope" = "Slope", "Aspect" = "Aspect",
    "ORCDRC" = "Soil Organic C", "PHIHOX" = "Soil pH",
    "CECSOL" = "Soil CEC", "CLYPPT" = "Clay Content",
    "SLTPPT" = "Silt Content", "BDTICM" = "Bulk Density",
    "Lights2009" = "Night Lights", "Built2009" = "Built-up",
    "Croplands2005" = "Croplands", "Pasture2009" = "Pasture"
)

# 获取显示名称（缺失则原样返回）
get_var_label <- function(x) {
    out <- var_labels[x]
    out[is.na(out)] <- x[is.na(out)]
    unname(out)
}

# 物种名格式化: Genus_species → G. species
fmt_sp <- function(x) {
    sapply(x, function(s) {
        p <- strsplit(s, "_")[[1]]
        if (length(p) >= 2) paste0(substr(p[1], 1, 1), ". ", p[2]) else s
    }, USE.NAMES = FALSE)
}

# 保存工具
save_fig <- function(plt, name, w, h) {
    ggsave(file.path(fig_dir, paste0(name, ".png")),
           plt, width = w, height = h, dpi = 1200, bg = "white")
    tryCatch(
        ggsave(file.path(fig_dir, paste0(name, ".svg")),
               plt, width = w, height = h, bg = "white"),
        error = function(e) cat(sprintf("  [SVG 跳过: %s]\n", e$message))
    )
    cat(sprintf("  ✓ %s 已保存\n", name))
}

# ══════════════════════════════════════════════════════════════════════════════
# 读取数据
# ══════════════════════════════════════════════════════════════════════════════
edges <- read.csv("output/case4_plant/all_dag_edges_plant.csv", stringsAsFactors = FALSE)
ate   <- read.csv("output/case4_plant/all_ate_results_plant.csv", stringsAsFactors = FALSE) %>%
    mutate(coef = as.numeric(coef), p_value = as.numeric(p_value))
roles <- read.csv("output/case4_plant/all_role_info_plant.csv", stringsAsFactors = FALSE)
scr   <- read.csv("output/case4_plant/all_screening_plant.csv", stringsAsFactors = FALSE)

all_species <- sort(unique(ate$species))
all_vars    <- sort(unique(ate$variable))
n_sp        <- length(all_species)
n_var       <- length(all_vars)

cat(sprintf("数据加载完成: %d 物种, %d 环境变量\n", n_sp, n_var))

# ══════════════════════════════════════════════════════════════════════════════
# 统一排序: 物种按 ATE 轮廓 Ward 层次聚类; 变量按跨物种 ATE 方差
# ══════════════════════════════════════════════════════════════════════════════
ate_mat <- ate %>%
    select(species, variable, coef) %>%
    pivot_wider(names_from = variable, values_from = coef, values_fill = 0) %>%
    column_to_rownames("species") %>%
    as.matrix()

sp_hc    <- hclust(dist(ate_mat), method = "ward.D2")
sp_order <- sp_hc$labels[sp_hc$order]

var_heterogeneity <- ate %>%
    group_by(variable) %>%
    summarise(v = var(coef, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(v))
var_order <- var_heterogeneity$variable

cat("物种聚类排序与变量异质性排序完成\n")

# ══════════════════════════════════════════════════════════════════════════════
# (a) DAG 边频谱 — 物种特有边 vs 通用边
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── 绘制 (a) DAG 边频谱 ──\n")

edge_freq <- edges %>%
    distinct(species, from, to) %>%
    group_by(from, to) %>%
    summarise(n_sp = n(), .groups = "drop")

# 保存边频率表
write.csv(edge_freq, file.path(tbl_dir, "fig3a_edge_frequency.csv"),
          row.names = FALSE)

freq_tab <- edge_freq %>%
    count(n_sp, name = "n_edges") %>%
    right_join(tibble(n_sp = 1:n_sp), by = "n_sp") %>%
    replace_na(list(n_edges = 0))

n_total_edges <- nrow(edge_freq)
n_univ  <- sum(edge_freq$n_sp == n_sp)
n_rare  <- sum(edge_freq$n_sp <= ceiling(n_sp * 0.25))

pa <- ggplot(freq_tab, aes(x = n_sp, y = n_edges, fill = n_sp)) +
    geom_col(width = 0.82, show.legend = FALSE) +
    scale_fill_gradient(low = "#AED6F1", high = "#154360") +
    scale_x_continuous(breaks = seq(0, n_sp, by = 4)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
    # 标注: 物种特有区域
    annotate("rect", xmin = 0.4, xmax = ceiling(n_sp * 0.25) + 0.5,
             ymin = -Inf, ymax = Inf, fill = "#AED6F1", alpha = 0.12) +
    annotate("text", x = ceiling(n_sp * 0.25) / 2 + 0.5,
             y = max(freq_tab$n_edges) * 1.08,
             label = sprintf("Species-specific\n(%d edges, %.0f%%)",
                             n_rare, 100 * n_rare / n_total_edges),
             size = 2.8, color = "#154360", fontface = "bold") +
    # 标注: 通用区域
    annotate("rect", xmin = n_sp - 3.5, xmax = n_sp + 0.5,
             ymin = -Inf, ymax = Inf, fill = "#154360", alpha = 0.06) +
    annotate("text", x = n_sp - 1.5,
             y = max(freq_tab$n_edges) * 1.08,
             label = sprintf("Universal\n(%d edges, %.0f%%)",
                             n_univ, 100 * n_univ / n_total_edges),
             size = 2.8, color = "#154360", fontface = "bold") +
    labs(
        title = "(a) DAG edge occurrence spectrum",
        subtitle = sprintf(
            "%d unique causal edges across %d species — most are species-specific",
            n_total_edges, n_sp),
        x = "Number of species sharing the edge",
        y = "Number of edges"
    ) +
    theme_pub() +
    theme(panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3))

save_fig(pa, "fig3a_dag_edge_spectrum", w = 7, h = 4.5)

# ══════════════════════════════════════════════════════════════════════════════
# (b) ATE 系数热力图 — 发散配色 (蓝负/红正)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── 绘制 (b) ATE 系数热力图 ──\n")

ate_plot <- ate %>%
    mutate(
        species  = factor(species, levels = sp_order),
        variable = factor(variable, levels = rev(var_order)),
        sig_mark = ifelse(p_value < 0.05, "*", "")
    )

# 保存 ATE 宽表
ate_wide <- ate %>%
    select(species, variable, coef) %>%
    pivot_wider(names_from = species, values_from = coef)
write.csv(ate_wide, file.path(tbl_dir, "fig3b_ate_matrix.csv"),
          row.names = FALSE)

ate_lim <- max(abs(ate$coef), na.rm = TRUE) * 1.02

pb <- ggplot(ate_plot, aes(x = species, y = variable, fill = coef)) +
    geom_tile(color = "white", linewidth = 0.25) +
    # 显著的 ATE 标记 *
    geom_text(aes(label = sig_mark), size = 2.5, color = "grey20", vjust = 0.8) +
    scale_fill_gradient2(
        low = "#2166AC", mid = "white", high = "#B2182B",
        midpoint = 0, limits = c(-ate_lim, ate_lim),
        name = "ATE\ncoefficient"
    ) +
    scale_y_discrete(labels = get_var_label) +
    scale_x_discrete(labels = fmt_sp) +
    labs(
        title = "(b) Species-specific ATE coefficients",
        subtitle = "Same variable exerts opposite causal effects on different species (* p < 0.05)",
        x = "", y = ""
    ) +
    theme_pub() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 5.5, face = "italic"),
        legend.key.height = unit(1.2, "cm")
    )

save_fig(pb, "fig3b_ate_heatmap", w = 12, h = 5)

# ══════════════════════════════════════════════════════════════════════════════
# (c) 因果角色热力图 — Root / Mediator / Terminal 三色分类
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── 绘制 (c) 因果角色热力图 ──\n")

role_colors <- c(
    "Root"     = "#4DAF4A",
    "Mediator" = "#FF7F00",
    "Terminal" = "#377EB8"
)

# 完整网格 (含缺失 = 不在 DAG 中)
role_full <- expand.grid(
    species  = all_species,
    variable = all_vars,
    stringsAsFactors = FALSE
) %>%
    left_join(roles %>% select(species, variable, group), by = c("species", "variable")) %>%
    mutate(
        group    = factor(group, levels = c("Root", "Mediator", "Terminal")),
        species  = factor(species, levels = sp_order),
        variable = factor(variable, levels = rev(var_order))
    )

# 保存角色宽表
role_wide <- role_full %>%
    select(species, variable, group) %>%
    mutate(group = as.character(group)) %>%
    pivot_wider(names_from = species, values_from = group)
write.csv(role_wide, file.path(tbl_dir, "fig3c_role_matrix.csv"),
          row.names = FALSE)

# 统计: 每个变量有多少种不同角色
n_distinct_roles <- role_full %>%
    filter(!is.na(group)) %>%
    group_by(variable) %>%
    summarise(n_roles = n_distinct(group), .groups = "drop")

pc <- ggplot(role_full, aes(x = species, y = variable, fill = group)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_manual(
        values = role_colors,
        name   = "Causal role",
        na.value = "#F0F0F0",
        drop   = FALSE
    ) +
    scale_y_discrete(labels = get_var_label) +
    scale_x_discrete(labels = fmt_sp) +
    labs(
        title = "(c) Species-specific causal role assignments",
        subtitle = "Same variable plays different causal roles (Root / Mediator / Terminal) across species; grey = not in DAG",
        x = "", y = ""
    ) +
    theme_pub() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 5.5, face = "italic"),
        legend.key.size = unit(0.45, "cm")
    )

save_fig(pc, "fig3c_causal_role_heatmap", w = 12, h = 5)

# ══════════════════════════════════════════════════════════════════════════════
# (d) 自适应筛选得分热力图 — 物种间归一化连续色
# ══════════════════════════════════════════════════════════════════════════════
cat("\n── 绘制 (d) 筛选得分热力图 ──\n")

scr_plot <- scr %>%
    group_by(species) %>%
    mutate(
        score_norm = (score_total - min(score_total, na.rm = TRUE)) /
            (max(score_total, na.rm = TRUE) - min(score_total, na.rm = TRUE) + 1e-9)
    ) %>%
    ungroup() %>%
    mutate(
        species  = factor(species, levels = sp_order),
        variable = factor(variable, levels = rev(var_order))
    )

# 保存筛选得分宽表
scr_wide <- scr %>%
    select(species, variable, score_total) %>%
    pivot_wider(names_from = species, values_from = score_total)
write.csv(scr_wide, file.path(tbl_dir, "fig3d_screening_score_matrix.csv"),
          row.names = FALSE)

pd <- ggplot(scr_plot, aes(x = species, y = variable, fill = score_norm)) +
    geom_tile(color = "white", linewidth = 0.25) +
    scale_fill_gradient(
        low  = "#F7FBFF",
        high = "#08306B",
        name = "Screening\nscore\n(normalized)"
    ) +
    scale_y_discrete(labels = get_var_label) +
    scale_x_discrete(labels = fmt_sp) +
    labs(
        title = "(d) CAST adaptive screening scores across species",
        subtitle = "Different species prioritize different variables — the pipeline adapts automatically",
        x = "", y = ""
    ) +
    theme_pub() +
    theme(
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 5.5, face = "italic"),
        legend.key.height = unit(1.2, "cm")
    )

save_fig(pd, "fig3d_screening_score_heatmap", w = 12, h = 5)

cat("\n════════════════════════════════════════\n")
cat("  Fig 3 完成: 4 个子图 + 4 个数据表已保存\n")
cat("════════════════════════════════════════\n")
