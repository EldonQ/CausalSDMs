#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 14b_dag_pc.R
# 版本: V2.0 (适配Bootstrap结果)
# 功能说明: 基于Bootstrap PC算法结果生成专业DAG可视化
# 特点:
#   - 支持有向边 (基于Bootstrap方向概率)
#   - 边的宽度映射连接强度 (Strength)
#   - 边的透明度/样式映射方向置信度 (Direction)
# 输入文件: output/14_causal/a_structure_pc/edges_pc.csv
# 输出文件: figures/14_causal/b_dag/dag_pc.png/svg
# 作者: CausalSDMs项目
# 日期: 2026-02-06
# ==============================================================================

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

options(repos = c(CRAN = "https://mirrors.sustech.edu.cn/CRAN/"))

ensure_package <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}

ensure_package("tidyverse")
ensure_package("igraph")
ensure_package("viridis")

suppressPackageStartupMessages({
    library(tidyverse)
    library(igraph)
    library(viridis)
})

dir.create("figures/14_causal/b_dag", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("生成专业级DAG可视化 (PC算法 Bootstrap版)\n")
cat("======================================\n\n")

# ============================================================================
# 可调参数
# ============================================================================
layout_type <- "graphopt"
vertex_size_range <- c(8, 25)
vertex_label_cex <- 0.9
edge_width_range <- c(1, 5) # 边宽映射 strength
edge_arrow_size <- 0.6
fig_width <- 28
fig_height <- 26
dpi_output <- 600

# 方向置信度阈值：低于此值的方向将被视为无向边
DIR_CONF_THRESHOLD <- 0.55

# ============================================================================
# 数据读取
# ============================================================================
edges_file <- "output/14_causal/a_structure_pc/edges_pc.csv"
if (!file.exists(edges_file)) {
    stop("错误: 未找到 edges_pc.csv，请先运行 14a_structure_learning_pc.R")
}

edges_df <- read.csv(edges_file, stringsAsFactors = FALSE)
var_groups_raw <- read.csv("scripts/variables_selected_47.csv", stringsAsFactors = FALSE)

cat(sprintf("  - 读取边数: %d\n", nrow(edges_df)))

# ============================================================================
# 变量名称映射与分组颜色
# ============================================================================
name_mapping <- c(
    "dem_avg" = "Elev", "dem_range" = "ElevRange", "slope_avg" = "Slope",
    "slope_range" = "SlopeRange", "flow_acc" = "FlowAcc", "flow_length" = "FlowLen",
    "hydro_wavg_01" = "BIO1", "hydro_wavg_02" = "BIO2", "hydro_wavg_03" = "BIO3",
    "hydro_wavg_04" = "BIO4", "hydro_wavg_05" = "BIO5", "hydro_wavg_06" = "BIO6",
    "hydro_wavg_07" = "BIO7", "hydro_wavg_08" = "BIO8", "hydro_wavg_09" = "BIO9",
    "hydro_wavg_10" = "BIO10", "hydro_wavg_11" = "BIO11", "hydro_wavg_12" = "BIO12",
    "hydro_wavg_13" = "BIO13", "hydro_wavg_14" = "BIO14", "hydro_wavg_15" = "BIO15",
    "hydro_wavg_16" = "BIO16", "hydro_wavg_17" = "BIO17", "hydro_wavg_18" = "BIO18",
    "hydro_wavg_19" = "BIO19",
    "lc_wavg_01" = "LC_Conif", "lc_wavg_02" = "LC_EBL", "lc_wavg_03" = "LC_DBL",
    "lc_wavg_04" = "LC_Mixed", "lc_wavg_05" = "LC_Shrub", "lc_wavg_06" = "LC_Herb",
    "lc_wavg_07" = "LC_Agri", "lc_wavg_08" = "LC_Flood", "lc_wavg_09" = "LC_Urban",
    "lc_wavg_10" = "LC_Snow", "lc_wavg_11" = "LC_Barren", "lc_wavg_12" = "LC_Water",
    "soil_wavg_01" = "SOC", "soil_wavg_02" = "pH", "soil_wavg_03" = "Sand",
    "soil_wavg_04" = "Silt", "soil_wavg_05" = "Clay", "soil_wavg_06" = "Coite",
    "soil_wavg_07" = "CEC", "soil_wavg_08" = "BulkDen", "soil_wavg_09" = "BedDepth",
    "soil_wavg_10" = "BedProb"
)

var_groups <- var_groups_raw %>%
    dplyr::select(variable, category) %>%
    dplyr::mutate(
        new_name = ifelse(variable %in% names(name_mapping), name_mapping[variable], variable),
        group_label = dplyr::case_when(
            category %in% c("地形", "水文") ~ "Topography",
            category == "水文气候" ~ "Climate",
            category == "土地覆盖" ~ "Land Cover",
            category == "土壤" ~ "Soil",
            TRUE ~ "Other"
        )
    )

group_colors <- c(
    "Topography" = "#E41A1C", "Climate" = "#377EB8",
    "Land Cover" = "#4DAF4A", "Soil" = "#984EA3", "Other" = "#666666"
)

# ============================================================================
# 构建 igraph 对象
# ============================================================================
cat("构建网络 ...\n")

# 应用名称映射
edges_plot <- edges_df %>%
    mutate(
        from = ifelse(from %in% names(name_mapping), name_mapping[from], from),
        to = ifelse(to %in% names(name_mapping), name_mapping[to], to),
        # 定义边类型：如果有向概率高，则有向；否则无向
        is_directed = direction > DIR_CONF_THRESHOLD | direction < (1 - DIR_CONF_THRESHOLD),
        # 统一方向：确保 direction > 0.5
        final_strength = strength
    )

# 重新整理边的方向，确保 from->to 是高概率方向
# 如果原始数据中 direction < 0.5，说明实际主要是 to->from
edges_plot_final <- edges_plot %>%
    rowwise() %>%
    mutate(
        temp_from = ifelse(direction < 0.5, to, from),
        temp_to = ifelse(direction < 0.5, from, to),
        final_direction = ifelse(direction < 0.5, 1 - direction, direction)
    ) %>%
    ungroup() %>%
    select(
        from = temp_from, to = temp_to, strength = final_strength,
        direction = final_direction, is_directed
    )

# 构建各有向属性的图
g <- graph_from_data_frame(edges_plot_final, directed = TRUE)

# 节点属性
node_names <- V(g)$name
node_data <- data.frame(name = node_names) %>%
    left_join(var_groups, by = c("name" = "new_name")) %>%
    mutate(group_label = ifelse(is.na(group_label), "Other", group_label))

V(g)$group <- node_data$group_label
V(g)$color <- group_colors[V(g)$group]

# 节点大小 (基于度) - 此时基于 Total Degree (In + Out)
V(g)$degree <- degree(g, mode = "all")
deg_norm <- (V(g)$degree - min(V(g)$degree)) / (max(V(g)$degree) - min(V(g)$degree) + 1e-6)
V(g)$size <- vertex_size_range[1] + deg_norm * (vertex_size_range[2] - vertex_size_range[1])

# Hub节点标记
top_degree <- quantile(V(g)$degree, 0.85)
V(g)$frame.color <- ifelse(V(g)$degree >= top_degree, "gold", "white")
V(g)$frame.width <- ifelse(V(g)$degree >= top_degree, 3, 1.5)

# 边属性设置
# 宽度映射 Strength
E(g)$width <- edge_width_range[1] +
    (edges_plot_final$strength - min(edges_plot_final$strength)) /
        (max(edges_plot_final$strength) - min(edges_plot_final$strength) + 1e-6) *
        (edge_width_range[2] - edge_width_range[1])

# 箭头模式
# 2 = 箭头 (有向), 0 = 无箭头 (无向)
E(g)$arrow.mode <- ifelse(edges_plot_final$is_directed, 2, 0)
E(g)$arrow.size <- edge_arrow_size

# 边颜色: 有向边深色，无向边浅色/虚线感
E(g)$color <- ifelse(edges_plot_final$is_directed,
    adjustcolor("#2C3E50", alpha.f = 0.8),
    adjustcolor("grey50", alpha.f = 0.6)
)

cat(sprintf("  - 节点数: %d\n", vcount(g)))
cat(sprintf("  - 边数: %d (其中有向边: %d)\n", ecount(g), sum(edges_plot_final$is_directed)))

# ============================================================================
# 计算布局
# ============================================================================
cat("计算布局 ...\n")
set.seed(42)
if (layout_type == "graphopt") {
    layout_coords <- layout_with_graphopt(g, niter = 5000, charge = 0.03, mass = 30)
} else {
    layout_coords <- layout_with_fr(g, niter = 2000)
}
layout_coords <- layout.norm(layout_coords, xmin = -1, xmax = 1, ymin = -1, ymax = 1) * 6

# ============================================================================
# 绘图输出
# ============================================================================
cat("生成图件 ...\n")
png_path <- "figures/14_causal/b_dag/dag_pc.png"
png(png_path, width = fig_width, height = fig_height, units = "in", res = dpi_output, bg = "transparent")
par(mar = c(2, 1, 4, 8))

plot(g,
    layout = layout_coords,
    vertex.label.font = 2,
    vertex.label.color = "black",
    edge.curved = 0.2, # 增加弯曲度，显示双向关系
    main = "Probabilistic PC Algorithm: Bootstrapped Causal Network"
)

# 标题信息
mtext(
    paste0("Bootstrap R=800 | Strength Threshold=0.55 | Direction Confidence > ", DIR_CONF_THRESHOLD),
    side = 3, line = 0, cex = 1.0, col = "grey40"
)

# 图例1: 变量类别
legend("bottomright",
    legend = names(group_colors), pch = 21, pt.bg = group_colors,
    pt.cex = 3.0, bty = "n", title = "Variable Category", cex = 1.5
)

# 图例2: 边含义
legend("topright",
    legend = c(
        "Directed (High Confidence)",
        "Undirected/Ambiguous"
    ),
    lwd = 3,
    col = c(adjustcolor("#2C3E50", 0.8), adjustcolor("grey50", 0.6)),
    lty = c(1, 1),
    bty = "n", title = "Edge Type (PC-Stable)", cex = 1.3
)

dev.off()
cat(sprintf("✓ 图件已生成: %s\n", png_path))

# SVG
svg_path <- "figures/14_causal/b_dag/dag_pc.svg"
svg(svg_path, width = fig_width, height = fig_height)
plot(g,
    layout = layout_coords, vertex.label.font = 2, edge.curved = 0.2,
    main = "Probabilistic PC Algorithm: Bootstrapped Causal Network"
)
dev.off()

cat("\n脚本执行完成!\n")
