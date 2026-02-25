#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 14b_dag_hc.R
# 功能说明: 基于HC算法结果生成专业DAG可视化 (复用14b_dag_publication_v2.R逻辑)
# 输入文件: output/14_causal/a_structure/edges_hc.csv, graph_hc_avg.rds
# 输出文件: figures/14_causal/b_dag/dag_hc.png/svg
# 作者: CausalSDMs项目
# 日期: 2026-02-06
# ==============================================================================

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

options(repos = c(CRAN = "https://mirrors.sustech.edu.cn/CRAN/"))

# ============================================================================
# 安装和加载必要的包
# ============================================================================
ensure_package <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("正在安装包:", pkg))
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
cat("生成专业级DAG可视化 (HC算法)\n")
cat("======================================\n\n")

# ============================================================================
# 可调参数
# ============================================================================
threshold <- 0.65
layout_type <- "graphopt"
vertex_size_range <- c(8, 25)
vertex_label_cex <- 0.9
vertex_label_dist <- 0
edge_width_range <- c(0.7, 7)
edge_arrow_size <- 0.6
edge_curved <- 0
fig_width <- 28
fig_height <- 26
dpi_output <- 600

# ============================================================================
# 数据读取
# ============================================================================
cat("读取HC算法输出 ...\n")

if (!file.exists("output/14_causal/a_structure/edges_hc.csv")) {
    stop("错误: 未找到 edges_hc.csv，请先运行 14a_structure_learning.R")
}

edges_strength <- read.csv("output/14_causal/a_structure/edges_hc.csv", stringsAsFactors = FALSE)
avg_hc <- readRDS("output/14_causal/a_structure/graph_hc_avg.rds")
var_groups_raw <- read.csv("scripts/variables_selected_47.csv", stringsAsFactors = FALSE)

# ============================================================================
# 变量名称映射
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
# 构建igraph对象
# ============================================================================
cat("构建网络 ...\n")

edges_filtered <- edges_strength %>%
    filter(strength >= threshold) %>%
    mutate(
        from = ifelse(from %in% names(name_mapping), name_mapping[from], from),
        to = ifelse(to %in% names(name_mapping), name_mapping[to], to)
    )

cat("  - 阈值:", threshold, "\n")
cat("  - 筛选后边数:", nrow(edges_filtered), "\n")

g <- graph_from_data_frame(d = edges_filtered[, c("from", "to")], directed = TRUE)
E(g)$weight <- edges_filtered$strength

V(g)$out_degree <- degree(g, mode = "out")
V(g)$in_degree <- degree(g, mode = "in")
V(g)$total_degree <- degree(g, mode = "all")

node_df <- data.frame(name = V(g)$name) %>%
    left_join(var_groups, by = c("name" = "new_name")) %>%
    mutate(group_label = ifelse(is.na(group_label), "Other", group_label))

V(g)$group <- node_df$group_label
V(g)$color <- group_colors[V(g)$group]

out_degrees <- V(g)$out_degree
if (max(out_degrees) > min(out_degrees)) {
    V(g)$size <- vertex_size_range[1] +
        (vertex_size_range[2] - vertex_size_range[1]) *
            (out_degrees - min(out_degrees)) / (max(out_degrees) - min(out_degrees))
} else {
    V(g)$size <- mean(vertex_size_range)
}

top_out_degree <- quantile(V(g)$out_degree, 0.85)
V(g)$frame.color <- ifelse(V(g)$out_degree >= top_out_degree, "gold", "white")
V(g)$frame.width <- ifelse(V(g)$out_degree >= top_out_degree, 3, 1.5)

weights <- E(g)$weight
if (max(weights) > min(weights)) {
    E(g)$width <- edge_width_range[1] +
        (edge_width_range[2] - edge_width_range[1]) *
            (weights - min(weights)) / (max(weights) - min(weights))
    E(g)$arrow.size <- E(g)$width * 0.15
} else {
    E(g)$width <- mean(edge_width_range)
    E(g)$arrow.size <- edge_arrow_size
}

edge_colors <- sapply(weights, function(w) {
    if (w >= 0.95) {
        adjustcolor("#B22222", alpha.f = 1)
    } else if (w >= 0.85) {
        adjustcolor("#FF6347", alpha.f = 1)
    } else if (w >= 0.75) {
        adjustcolor("#FF8C00", alpha.f = 1)
    } else {
        adjustcolor("grey50", alpha.f = 1)
    }
})
E(g)$color <- edge_colors

cat("  - 节点数:", vcount(g), "\n")
cat("  - 边数:", ecount(g), "\n")

# ============================================================================
# 计算布局
# ============================================================================
set.seed(42)
if (layout_type == "graphopt") {
    layout_coords <- layout_with_graphopt(g,
        niter = 5000, charge = 0.03, mass = 50,
        spring.length = 2, spring.constant = 0.5
    )
} else {
    layout_coords <- layout_with_fr(g, niter = 2000)
}
layout_coords <- layout.norm(layout_coords, xmin = -1, xmax = 1, ymin = -1, ymax = 1) * 6

# ============================================================================
# 绘图 - PNG
# ============================================================================
cat("生成DAG图 ...\n")

png_path <- "figures/14_causal/b_dag/dag_hc.png"
png(png_path, width = fig_width, height = fig_height, units = "in", res = dpi_output, bg = "transparent")
par(mar = c(2, 1, 4, 8))

plot(g,
    layout = layout_coords,
    vertex.color = V(g)$color, vertex.size = V(g)$size,
    vertex.frame.color = V(g)$frame.color, vertex.frame.width = V(g)$frame.width,
    vertex.label = V(g)$name, vertex.label.cex = vertex_label_cex,
    vertex.label.color = "black", vertex.label.font = 2, vertex.label.dist = vertex_label_dist,
    edge.width = E(g)$width, edge.color = E(g)$color,
    edge.arrow.size = E(g)$arrow.size, edge.curved = edge_curved, edge.arrow.mode = 2,
    main = "HC Algorithm Causal Network (DAG)"
)

mtext(paste0("Edges >= ", threshold, " bootstrap stability | N=", vcount(g), " nodes, E=", ecount(g), " edges"),
    side = 3, line = 0, cex = 1.0, col = "grey40"
)

legend("bottomright",
    legend = names(group_colors), pch = 21, pt.bg = group_colors,
    pt.cex = 3.0, col = "white", bty = "n", title = "Variable Category", title.font = 2, cex = 1.5
)
legend("right",
    legend = c("Low", "Medium", "High"), pch = 21, pt.bg = "grey60",
    pt.cex = c(2.0, 3.0, 4.0), col = "white", bty = "n",
    title = "Out-degree\n(Causal Influence)", title.font = 2, cex = 1.3, y.intersp = 1.5
)
legend("topright",
    legend = c("≥0.95 (Very High)", "0.85-0.95 (High)", "0.75-0.85 (Medium)", "<0.75 (Lower)"),
    lwd = c(6, 4.5, 3, 1.5),
    col = c(
        adjustcolor("#B22222", 0.85), adjustcolor("#FF6347", 0.75),
        adjustcolor("#FF8C00", 0.65), adjustcolor("grey50", 0.45)
    ),
    bty = "n", title = "Edge Color\n(Bootstrap Stability)", title.font = 2, cex = 1.2, seg.len = 2.5
)
legend("topleft",
    legend = c("Core Node (Top 15% Out-degree)"), pch = 21, pt.bg = "grey70",
    pt.cex = 3.0, col = "gold", lwd = 3, bty = "n", title = "Node Highlight", title.font = 2, cex = 1.2
)

dev.off()

# SVG
svg_path <- "figures/14_causal/b_dag/dag_hc.svg"
svg(svg_path, width = fig_width, height = fig_height, bg = "white")
par(mar = c(1, 1, 4, 1))
plot(g,
    layout = layout_coords,
    vertex.color = V(g)$color, vertex.size = V(g)$size,
    vertex.frame.color = "white", vertex.frame.width = 2,
    vertex.label = V(g)$name, vertex.label.cex = vertex_label_cex,
    vertex.label.color = "black", vertex.label.font = 2, vertex.label.dist = vertex_label_dist,
    edge.width = E(g)$width, edge.color = E(g)$color,
    edge.arrow.size = edge_arrow_size, edge.curved = edge_curved,
    main = "HC Algorithm Causal Network (DAG)"
)
legend("bottomright",
    legend = names(group_colors), pch = 21, pt.bg = group_colors,
    pt.cex = 2, col = "white", bty = "n", title = "Variable Category", title.font = 2, cex = 1.2
)
dev.off()

cat("\n✓ DAG图已生成:\n")
cat("  -", png_path, "\n")
cat("  -", svg_path, "\n")

# ============================================================================
# 输出关键节点信息
# ============================================================================
cat("\n★ 关键节点分析 ★\n")
cat("【因果源头节点】(出度 Top 5):\n")
source_nodes <- data.frame(name = V(g)$name, out_degree = V(g)$out_degree, group = V(g)$group) %>%
    arrange(desc(out_degree)) %>%
    head(5)
print(source_nodes)

cat("\n完成!\n")
