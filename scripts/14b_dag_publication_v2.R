#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 14b_dag_publication_v2.R
# 功能说明: 使用 ggdag 风格生成专业的因果DAG图
# 特点:
#   - 清晰的圆形布局或力导向布局
#   - 节点标签清晰可见
#   - 简洁的配色方案
#   - 适合学术发表
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

dir.create("figures/14_causal", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("生成专业级DAG可视化 (igraph原生版)\n")
cat("======================================\n\n")

# ============================================================================
# 可调参数
# ============================================================================
threshold <- 0.65 # 边稳定性阈值 (稍微增大以减少边数)
layout_type <- "graphopt" # 布局类型: "graphopt"(均匀分布), "fr", "kk", "circle"

# 节点参数
vertex_size_range <- c(8, 25) # 节点大小范围 (缩小以减少重叠)
vertex_label_cex <- 0.9 # 标签大小比例 (稍小)
vertex_label_dist <- 0 # 标签距离节点中心的距离 (0=在节点上)

# 边参数
edge_width_range <- c(0.7, 7) # 边宽度范围
edge_arrow_size <- 0.6 # 箭头大小
edge_curved <- 0 # 边曲率 (0=直线)

# 输出参数 - 增大画布
fig_width <- 28 # 图片宽度 (大幅增大)
fig_height <- 26 # 图片高度 (大幅增大)
dpi_output <- 600 # 分辨率

# ============================================================================
# 数据读取
# ============================================================================
cat("读取因果分析输出 ...\n")

edges_strength <- read.csv("output/14_causal/edges_summary.csv", stringsAsFactors = FALSE)
avg_hc <- readRDS("output/14_causal/graph_hc_avg.rds")

# 读取变量分组
var_groups_raw <- read.csv("scripts/variables_selected_47.csv", stringsAsFactors = FALSE)

# ============================================================================
# 定义变量名称映射: 旧名 -> 新名（符合学术规范的简短名称）
# ============================================================================
name_mapping <- c(
    # 地形与河网拓扑 (6变量)
    "dem_avg"       = "Elev",
    "dem_range"     = "ElevRange",
    "slope_avg"     = "Slope",
    "slope_range"   = "SlopeRange",
    "flow_acc"      = "FlowAcc",
    "flow_length"   = "FlowLen",
    # 上游加权水文气候 (19变量)
    "hydro_wavg_01" = "BIO1",
    "hydro_wavg_02" = "BIO2",
    "hydro_wavg_03" = "BIO3",
    "hydro_wavg_04" = "BIO4",
    "hydro_wavg_05" = "BIO5",
    "hydro_wavg_06" = "BIO6",
    "hydro_wavg_07" = "BIO7",
    "hydro_wavg_08" = "BIO8",
    "hydro_wavg_09" = "BIO9",
    "hydro_wavg_10" = "BIO10",
    "hydro_wavg_11" = "BIO11",
    "hydro_wavg_12" = "BIO12",
    "hydro_wavg_13" = "BIO13",
    "hydro_wavg_14" = "BIO14",
    "hydro_wavg_15" = "BIO15",
    "hydro_wavg_16" = "BIO16",
    "hydro_wavg_17" = "BIO17",
    "hydro_wavg_18" = "BIO18",
    "hydro_wavg_19" = "BIO19",
    # 上游加权土地覆盖 (12变量)
    "lc_wavg_01"    = "LC_Conif",
    "lc_wavg_02"    = "LC_EBL",
    "lc_wavg_03"    = "LC_DBL",
    "lc_wavg_04"    = "LC_Mixed",
    "lc_wavg_05"    = "LC_Shrub",
    "lc_wavg_06"    = "LC_Herb",
    "lc_wavg_07"    = "LC_Agri",
    "lc_wavg_08"    = "LC_Flood",
    "lc_wavg_09"    = "LC_Urban",
    "lc_wavg_10"    = "LC_Snow",
    "lc_wavg_11"    = "LC_Barren",
    "lc_wavg_12"    = "LC_Water",
    # 上游加权土壤属性 (10变量)
    "soil_wavg_01"  = "SOC",
    "soil_wavg_02"  = "pH",
    "soil_wavg_03"  = "Sand",
    "soil_wavg_04"  = "Silt",
    "soil_wavg_05"  = "Clay",
    "soil_wavg_06"  = "Coite",
    "soil_wavg_07"  = "CEC",
    "soil_wavg_08"  = "BulkDen",
    "soil_wavg_09"  = "BedDepth",
    "soil_wavg_10"  = "BedProb"
)

cat("  - 变量名称映射已定义 (共", length(name_mapping), "个变量)\n")

# 构建分组映射（使用新变量名）
var_groups <- var_groups_raw %>%
    dplyr::select(variable, category) %>%
    dplyr::mutate(
        # 将旧名替换为新名
        new_name = ifelse(variable %in% names(name_mapping),
            name_mapping[variable],
            variable
        ),
        group_label = dplyr::case_when(
            category %in% c("地形", "水文") ~ "Topography",
            category == "水文气候" ~ "Climate",
            category == "土地覆盖" ~ "Land Cover",
            category == "土壤" ~ "Soil",
            TRUE ~ "Other"
        )
    )

# 配色方案 - 更清晰的颜色
group_colors <- c(
    "Topography"  = "#E41A1C", # 红色
    "Climate"     = "#377EB8", # 蓝色
    "Land Cover"  = "#4DAF4A", # 绿色
    "Soil"        = "#984EA3", # 紫色
    "Other"       = "#666666" # 灰色
)

# ============================================================================
# 构建igraph对象
# ============================================================================
cat("构建网络 ...\n")

# 筛选高稳定性边
edges_filtered <- edges_strength %>%
    filter(strength >= threshold)

cat("  - 阈值:", threshold, "\n")
cat("  - 筛选后边数:", nrow(edges_filtered), "\n")

# ============================================================================
# 将边的变量名从旧名替换为新名
# ============================================================================
edges_filtered <- edges_filtered %>%
    mutate(
        from = ifelse(from %in% names(name_mapping), name_mapping[from], from),
        to   = ifelse(to %in% names(name_mapping), name_mapping[to], to)
    )

cat("  - 变量名称已转换为新命名规范\n")

# 创建有向图
g <- graph_from_data_frame(
    d = edges_filtered[, c("from", "to")],
    directed = TRUE
)

# 添加边权重
E(g)$weight <- edges_filtered$strength

# 计算节点度
V(g)$out_degree <- degree(g, mode = "out")
V(g)$in_degree <- degree(g, mode = "in")
V(g)$total_degree <- degree(g, mode = "all")

# 映射节点颜色 (使用新变量名)
node_df <- data.frame(name = V(g)$name)
node_df <- node_df %>%
    left_join(var_groups, by = c("name" = "new_name")) %>%
    mutate(group_label = ifelse(is.na(group_label), "Other", group_label))

V(g)$group <- node_df$group_label
V(g)$color <- group_colors[V(g)$group]

# 计算节点大小 (基于出度)
out_degrees <- V(g)$out_degree
if (max(out_degrees) > min(out_degrees)) {
    V(g)$size <- vertex_size_range[1] +
        (vertex_size_range[2] - vertex_size_range[1]) *
            (out_degrees - min(out_degrees)) / (max(out_degrees) - min(out_degrees))
} else {
    V(g)$size <- mean(vertex_size_range)
}

# ============================================================================
# 突出核心节点：高出度节点添加金色边框
# ============================================================================
top_out_degree <- quantile(V(g)$out_degree, 0.85) # Top 15% 出度节点
V(g)$frame.color <- ifelse(V(g)$out_degree >= top_out_degree, "gold", "white")
V(g)$frame.width <- ifelse(V(g)$out_degree >= top_out_degree, 3, 1.5)

# ============================================================================
# 边宽度和箭头大小 (基于稳定性，箭头与宽度成比例)
# ============================================================================
weights <- E(g)$weight

if (max(weights) > min(weights)) {
    # 边宽度
    E(g)$width <- edge_width_range[1] +
        (edge_width_range[2] - edge_width_range[1]) *
            (weights - min(weights)) / (max(weights) - min(weights))
    # 箭头大小与边宽度成比例
    E(g)$arrow.size <- E(g)$width * 0.15
} else {
    E(g)$width <- mean(edge_width_range)
    E(g)$arrow.size <- edge_arrow_size
}

# ============================================================================
# 边颜色：根据稳定性强度分级着色
# ============================================================================
edge_colors <- sapply(weights, function(w) {
    if (w >= 0.95) {
        # 极高强度: 深红色
        adjustcolor("#B22222", alpha.f = 1)
    } else if (w >= 0.85) {
        # 高强度: 橙红色
        adjustcolor("#FF6347", alpha.f = 1)
    } else if (w >= 0.75) {
        # 中等强度: 橙色
        adjustcolor("#FF8C00", alpha.f = 1)
    } else {
        # 较低强度: 灰色
        adjustcolor("grey50", alpha.f = 1)
    }
})
E(g)$color <- edge_colors

cat("  - 节点数:", vcount(g), "\n")
cat("  - 边数:", ecount(g), "\n")
cat("  - 高强度边 (≥0.85):", sum(weights >= 0.85), "\n")
cat("  - 核心节点 (Top出度):", sum(V(g)$out_degree >= top_out_degree), "\n")

# ============================================================================
# 计算布局
# ============================================================================
cat("计算布局 (", layout_type, ") ...\n")

set.seed(42) # 可重复性

if (layout_type == "graphopt") {
    # graphopt布局 - 产生更均匀的节点分布
    layout_coords <- layout_with_graphopt(
        g,
        niter = 5000, # 更多迭代次数
        charge = 0.03, # 节点排斥力 (增大使节点分散)
        mass = 50, # 节点质量
        spring.length = 2, # 弹簧长度 (边的理想长度)
        spring.constant = 0.5 # 弹簧常数
    )
} else if (layout_type == "fr") {
    # Fruchterman-Reingold布局
    layout_coords <- layout_with_fr(g, niter = 3000, area = vcount(g)^5)
} else if (layout_type == "kk") {
    # Kamada-Kawai布局 - 边长度更一致
    layout_coords <- layout_with_kk(g, maxiter = 5000)
} else if (layout_type == "circle") {
    layout_coords <- layout_in_circle(g)
} else {
    layout_coords <- layout_with_fr(g, niter = 2000)
}

# 标准化布局坐标到 [-1, 1] 范围
layout_coords <- layout.norm(layout_coords, xmin = -1, xmax = 1, ymin = -1, ymax = 1)

# 扩展布局以增加节点间距
layout_coords <- layout_coords * 6 # 大幅扩展确保节点分开

# ============================================================================
# 绘图 - 使用igraph原生绘图
# ============================================================================
cat("生成DAG图 ...\n")

# PNG输出
png_path <- "figures/14_causal/dag_publication_v2.png"
png(png_path, width = fig_width, height = fig_height, units = "in", res = dpi_output, bg = "transparent")

# 设置绘图参数 (增加右侧边距以放置图例)
par(mar = c(2, 1, 4, 8))

# 绘制网络
plot(g,
    layout = layout_coords,

    # 节点样式
    vertex.color = V(g)$color,
    vertex.size = V(g)$size,
    vertex.frame.color = V(g)$frame.color, # 使用自适应边框颜色
    vertex.frame.width = V(g)$frame.width, # 使用自适应边框宽度

    # 标签样式
    vertex.label = V(g)$name,
    vertex.label.cex = vertex_label_cex,
    vertex.label.color = "black",
    vertex.label.font = 2, # 粗体
    vertex.label.dist = vertex_label_dist,

    # 边样式 - 使用每边自适应的箭头大小
    edge.width = E(g)$width,
    edge.color = E(g)$color,
    edge.arrow.size = E(g)$arrow.size, # 箭头大小与边宽度成比例
    edge.curved = edge_curved,
    edge.arrow.mode = 2,

    # 其他
    main = "Consensus Causal Network (DAG) of Environmental Variables"
)

# 添加副标题
mtext(
    paste0("Edges >= ", threshold, " bootstrap stability | N=", vcount(g), " nodes, E=", ecount(g), " edges"),
    side = 3, line = 0, cex = 1.0, col = "grey40"
)

# 添加图例1: 变量类别 (右下)
legend("bottomright",
    legend = names(group_colors),
    pch = 21,
    pt.bg = group_colors,
    pt.cex = 3.0,
    col = "white",
    bty = "n",
    title = "Variable Category",
    title.font = 2,
    cex = 1.5
)

# 添加图例2: 节点大小 (右中)
legend("right",
    legend = c("Low", "Medium", "High"),
    pch = 21,
    pt.bg = "grey60",
    pt.cex = c(2.0, 3.0, 4.0),
    col = "white",
    bty = "n",
    title = "Out-degree\n(Causal Influence)",
    title.font = 2,
    cex = 1.3,
    y.intersp = 1.5
)

# 添加图例3: 边颜色分级 (右上)
legend("topright",
    legend = c("≥0.95 (Very High)", "0.85-0.95 (High)", "0.75-0.85 (Medium)", "<0.75 (Lower)"),
    lwd = c(6, 4.5, 3, 1.5),
    col = c(
        adjustcolor("#B22222", 0.85),
        adjustcolor("#FF6347", 0.75),
        adjustcolor("#FF8C00", 0.65),
        adjustcolor("grey50", 0.45)
    ),
    bty = "n",
    title = "Edge Color\n(Bootstrap Stability)",
    title.font = 2,
    cex = 1.2,
    seg.len = 2.5
)

# 添加图例4: 核心节点标识 (左上)
legend("topleft",
    legend = c("Core Node (Top 15% Out-degree)"),
    pch = 21,
    pt.bg = "grey70",
    pt.cex = 3.0,
    col = "gold",
    lwd = 3,
    bty = "n",
    title = "Node Highlight",
    title.font = 2,
    cex = 1.2
)

dev.off()

cat("\n✓ DAG图已生成:\n")
cat("  -", png_path, "\n")

# ============================================================================
# SVG输出 (矢量格式)
# ============================================================================
svg_path <- "figures/14_causal/dag_publication_v2.svg"
svg(svg_path, width = fig_width, height = fig_height, bg = "white")

par(mar = c(1, 1, 4, 1))

plot(g,
    layout = layout_coords,
    vertex.color = V(g)$color,
    vertex.size = V(g)$size,
    vertex.frame.color = "white",
    vertex.frame.width = 2,
    vertex.label = V(g)$name,
    vertex.label.cex = vertex_label_cex,
    vertex.label.color = "black",
    vertex.label.font = 2,
    vertex.label.dist = vertex_label_dist,
    edge.width = E(g)$width,
    edge.color = E(g)$color,
    edge.arrow.size = edge_arrow_size,
    edge.curved = edge_curved,
    main = "Consensus Causal Network (DAG) of Environmental Variables"
)

legend("bottomright",
    legend = names(group_colors),
    pch = 21,
    pt.bg = group_colors,
    pt.cex = 2,
    col = "white",
    bty = "n",
    title = "Variable Category",
    title.font = 2,
    cex = 1.2
)

mtext(
    paste0("Edges ≥ ", threshold, " stability | N=", vcount(g), " nodes, E=", ecount(g), " edges"),
    side = 3, line = 0, cex = 1.0, col = "grey40"
)

dev.off()

cat("  -", svg_path, "\n")

# ============================================================================
# 输出关键节点信息
# ============================================================================
cat("\n")
cat("======================================\n")
cat("★ 关键节点分析 ★\n")
cat("======================================\n")

# 源头节点（高出度）
cat("\n【因果源头节点】(出度 Top 10):\n")
source_nodes <- data.frame(
    name = V(g)$name,
    out_degree = V(g)$out_degree,
    group = V(g)$group
) %>%
    arrange(desc(out_degree)) %>%
    head(10)
print(source_nodes)

# 终点节点（高入度）
cat("\n【响应终点节点】(入度 Top 10):\n")
sink_nodes <- data.frame(
    name = V(g)$name,
    in_degree = V(g)$in_degree,
    group = V(g)$group
) %>%
    arrange(desc(in_degree)) %>%
    head(10)
print(sink_nodes)

cat("\n完成!\n")
