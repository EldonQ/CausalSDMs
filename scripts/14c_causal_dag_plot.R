#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 14c_causal_dag_plot.R
# 功能说明: 整合 14a 与 14b 的优点，基于真实因果DAG与bootstrap强度，
#           使用47变量真实分组，生成专业Nature级别的因果网络全套图件。
# 输出图件: figures/14_causal/
#   1) dag_hc_avg_full.png/svg            - 平均HC-DAG（去孤立点）
#   2) dag_hc_avg_lcc.png/svg             - 平均HC-DAG最大连通分量
#   3) dag_core_topN.png/svg              - 核心路径（Top N强边）
#   4) dag_cross_group.png/svg            - 跨组因果边
#   5) dag_strength_matrix.png/svg        - 边强度矩阵（≥min_strength）
#   6) dag_pc_layout.png/svg              - PC-DAG（与分组上色）
# 依赖输入: output/14_causal/graph_hc_avg.rds, edges_summary.csv, graph_pc.rds
#         scripts/variables_selected_47.csv（变量—组别）
# 备注: 图件英文标注、Arial、≥1200dpi、PNG+SVG，矢量优先
# ==============================================================================

rm(list = ls())
gc()
setwd("E:/SDM01")

options(repos = c(CRAN = "https://mirrors.sustech.edu.cn/CRAN/"))

# 加载包
pkgs <- c("tidyverse", "bnlearn", "igraph", "ggraph", "tidygraph",
          "ggrepel", "reshape2", "sysfonts", "showtext")
for (p in pkgs) {
  if (!require(p, character.only = TRUE)) {
    install.packages(p, dependencies = TRUE)
    library(p, character.only = TRUE)
  }
}

# 字体
try({
  sysfonts::font_add(family = "Arial",
                      regular = "C:/Windows/Fonts/arial.ttf",
                      bold    = "C:/Windows/Fonts/arialbd.ttf")
  showtext::showtext_opts(dpi = 1200)
  showtext::showtext_auto(enable = TRUE)
}, silent = TRUE)

dir.create("figures/14_causal", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("整合因果DAG绘图 (14c)\n")
cat("======================================\n\n")

# ------------------------ 参数 ------------------------
threshold <- 0.55     # 平均网络保留边阈值（与14_discovery一致）
topN      <- 50       # 核心路径Top N边
min_strength_matrix <- 0.30 # 矩阵热图的最小强度

# ------------------------ 数据读取 ------------------------
cat("读取 discovery 输出 ...\n")

edges_strength <- read.csv("output/14_causal/edges_summary.csv", stringsAsFactors = FALSE)

avg_hc <- NULL
if (file.exists("output/14_causal/graph_hc_avg.rds")) {
  avg_hc <- readRDS("output/14_causal/graph_hc_avg.rds")
} else {
  warning("找不到 graph_hc_avg.rds，将只基于 edges_summary 作图")
}

pc_fit <- NULL
if (file.exists("output/14_causal/graph_pc.rds")) {
  pc_fit <- readRDS("output/14_causal/graph_pc.rds")
} else {
  warning("找不到 graph_pc.rds，PC图将跳过")
}

# 真实变量分组（首选CSV；否则按前缀回退）
cat("读取变量分组 ...\n")
if (file.exists("scripts/variables_selected_47.csv")) {
  var_groups <- read.csv("scripts/variables_selected_47.csv", stringsAsFactors = FALSE) %>%
    dplyr::select(variable, category) %>%
    dplyr::rename(name = variable, group = category) %>%
    dplyr::mutate(
      group_label = dplyr::case_when(
        group == "G1_TopoSlopeFlow" ~ "Topography & Flow",
        group == "G2_Hydroclim_wavg" ~ "Hydroclimatic",
        group == "G3_Landcover_wavg" ~ "Land Cover",
        group == "G4_Soil_wavg" ~ "Soil Properties",
        TRUE ~ "Other")
    )
} else {
  warning("找不到 variables_selected_47.csv，按变量前缀推断分组")
  all_vars <- unique(c(edges_strength$from, edges_strength$to))
  var_groups <- tibble::tibble(name = all_vars) %>%
    dplyr::mutate(group_label = dplyr::case_when(
      grepl("^(dem_|slope_|flow_)", name) ~ "Topography & Flow",
      grepl("^hydro_", name) ~ "Hydroclimatic",
      grepl("^lc_", name) ~ "Land Cover",
      grepl("^soil_", name) ~ "Soil Properties",
      TRUE ~ "Other"))
}

# 配色（Nature风格）
group_colors <- c(
  "Topography & Flow" = "#E41A1C",
  "Hydroclimatic"     = "#377EB8",
  "Land Cover"        = "#4DAF4A",
  "Soil Properties"   = "#984EA3",
  "Other"             = "#999999"
)

# ------------------------ 构造平均网络边集 ------------------------
cat("构造平均网络边集 ...\n")

if (!is.null(avg_hc)) {
  arcs_avg <- bnlearn::arcs(avg_hc) %>% as.data.frame(stringsAsFactors = FALSE)
  colnames(arcs_avg) <- c("from","to")
} else {
  # 回退：使用强度≥threshold的边作为平均网络近似
  arcs_avg <- edges_strength %>% dplyr::filter(strength >= threshold) %>% dplyr::select(from,to)
}

# 合并强度（双向合并取最大）
strength_df <- edges_strength %>% dplyr::select(from,to,strength)
strength_rev <- strength_df %>% dplyr::transmute(from = to, to = from, strength_rev = strength)
edges_avg <- arcs_avg %>%
  dplyr::left_join(strength_df, by = c("from","to")) %>%
  dplyr::left_join(strength_rev, by = c("from","to")) %>%
  dplyr::mutate(strength = dplyr::coalesce(strength, strength_rev, threshold))

vertices_df <- var_groups %>% dplyr::select(name, group_label)
g_avg <- igraph::graph_from_data_frame(edges_avg, directed = TRUE, vertices = vertices_df)

# 去除孤立点及LCC
deg_all <- igraph::degree(g_avg, mode = "all")
g_full <- igraph::induced_subgraph(g_avg, vids = which(deg_all > 0))
comp <- igraph::components(igraph::as.undirected(g_full))
g_lcc <- if (comp$no > 1) igraph::induced_subgraph(g_full, which(comp$membership == which.max(comp$csize))) else g_full

to_tidy <- function(g) {
  tidygraph::as_tbl_graph(g) %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate(
      degree = tidygraph::centrality_degree(mode = "all"),
      group_label = factor(group_label, levels = names(group_colors))
    ) %>%
    tidygraph::activate(edges) %>%
    dplyr::mutate(strength = strength)
}

tg_full <- to_tidy(g_full)
tg_lcc  <- to_tidy(g_lcc)

# 通用节点/边绘制组件
plot_network <- function(tg, title_text, file_stub, layout = "kk") {
  set.seed(42)
  p <- ggraph::ggraph(tg, layout = layout) +
    ggraph::geom_edge_link(aes(width = strength, alpha = strength),
                           arrow = grid::arrow(length = grid::unit(3, "mm"), type = "closed"),
                           end_cap = ggraph::circle(3.5, "mm"), colour = "grey30", lineend = "round") +
    ggraph::scale_edge_width_continuous(range = c(0.5, 3.2), name = "Edge strength") +
    ggraph::scale_edge_alpha_continuous(range = c(0.35, 0.95), guide = "none") +
    ggraph::geom_node_point(aes(fill = group_label, size = degree), shape = 21, colour = "black", stroke = 0.6) +
    scale_fill_manual(values = group_colors, name = "Variable group") +
    ggplot2::scale_size_continuous(range = c(3, 9), name = "Node degree") +
    ggrepel::geom_text_repel(aes(x = x, y = y, label = name), size = 2.8, family = "Arial",
                             box.padding = grid::unit(0.35, "lines"), point.padding = grid::unit(0.3, "lines"),
                             segment.color = "grey60", segment.size = 0.25, max.overlaps = 100) +
    labs(title = title_text, subtitle = paste0("Edges ≥ ", threshold, "; layout=", layout)) +
    theme_void(base_family = "Arial")

  ggsave(paste0("figures/14_causal/", file_stub, ".png"), p, width = 14, height = 10, units = "in", dpi = 1200, bg = "white")
  ggsave(paste0("figures/14_causal/", file_stub, ".svg"), p, width = 14, height = 10, units = "in", bg = "white")
}

# 1) 完整平均网络（非孤立点）
cat("绘制: 平均网络(去孤立点) ...\n")
plot_network(tg_full, "Averaged HC-DAG (Non-isolated nodes)", "dag_hc_avg_full", layout = "kk")

# 2) 最大连通分量
cat("绘制: 平均网络(LCC) ...\n")
plot_network(tg_lcc,  "Averaged HC-DAG (Largest Connected Component)", "dag_hc_avg_lcc", layout = "kk")

# 3) 核心路径（TopN）
cat("绘制: 核心路径 TopN ...\n")
edges_top <- edges_strength %>% dplyr::arrange(dplyr::desc(strength)) %>% dplyr::slice(1:topN)
g_core <- igraph::graph_from_data_frame(edges_top, directed = TRUE, vertices = vertices_df)
tg_core <- to_tidy(g_core)
set.seed(123)
p_core <- ggraph::ggraph(tg_core, layout = "kk") +
  ggraph::geom_edge_arc(aes(width = strength, color = strength),
                        arrow = grid::arrow(length = grid::unit(3, "mm"), type = "closed"),
                        end_cap = ggraph::circle(3.5, "mm"), strength = 0.12, lineend = "round") +
  ggraph::scale_edge_width_continuous(range = c(0.6, 3.6), name = "Edge strength") +
  ggraph::scale_edge_color_gradient(low = "grey70", high = "grey10", name = "Edge strength") +
  ggraph::geom_node_point(aes(fill = group_label, size = degree), shape = 21, colour = "black", stroke = 0.6) +
  scale_fill_manual(values = group_colors, name = "Variable group") +
  ggplot2::scale_size_continuous(range = c(4, 11), name = "Node degree") +
  ggrepel::geom_text_repel(aes(x = x, y = y, label = name), size = 3.0, family = "Arial",
                           box.padding = grid::unit(0.45, "lines"), point.padding = grid::unit(0.4, "lines"),
                           segment.color = "grey60", segment.size = 0.3, max.overlaps = 120, fontface = "bold") +
  labs(title = paste0("Core Causal Pathways (Top ", topN, ")"),
       subtitle = "Edge width/color = bootstrap strength") +
  theme_void(base_family = "Arial")
ggsave("figures/14_causal/dag_core_topN.png", p_core, width = 16, height = 12, units = "in", dpi = 1200, bg = "white")
ggsave("figures/14_causal/dag_core_topN.svg", p_core, width = 16, height = 12, units = "in", bg = "white")

# 4) 跨组因果关系（取TopN中跨组边）
cat("绘制: 跨组因果边 ...\n")
edges_cross <- edges_top %>%
  dplyr::left_join(var_groups, by = c("from" = "name")) %>% dplyr::rename(from_group = group_label) %>%
  dplyr::left_join(var_groups, by = c("to"   = "name")) %>% dplyr::rename(to_group   = group_label) %>%
  dplyr::filter(from_group != to_group)
g_cross <- igraph::graph_from_data_frame(edges_cross, directed = TRUE, vertices = vertices_df)
tg_cross <- to_tidy(g_cross)
p_cross <- ggraph::ggraph(tg_cross, layout = "graphopt") +
  ggraph::geom_edge_link(aes(width = strength, alpha = strength),
                         arrow = grid::arrow(length = grid::unit(3, "mm"), type = "closed"),
                         end_cap = ggraph::circle(3.5, "mm"), colour = "grey25", lineend = "round") +
  ggraph::scale_edge_width_continuous(range = c(0.6, 3.5), name = "Edge strength") +
  ggraph::scale_edge_alpha_continuous(range = c(0.4, 1), guide = "none") +
  ggraph::geom_node_point(aes(fill = group_label, size = degree), shape = 21, colour = "black", stroke = 0.6) +
  scale_fill_manual(values = group_colors, name = "Variable group") +
  ggplot2::scale_size_continuous(range = c(4, 10), name = "Node degree") +
  ggrepel::geom_text_repel(aes(x = x, y = y, label = name), size = 3.0, family = "Arial",
                           box.padding = grid::unit(0.4, "lines"), point.padding = grid::unit(0.35, "lines"),
                           segment.color = "grey60", segment.size = 0.25, max.overlaps = 100) +
  labs(title = "Cross-Domain Causal Relationships",
       subtitle = "Edges among different variable groups (Top N set)") +
  theme_void(base_family = "Arial")
ggsave("figures/14_causal/dag_cross_group.png", p_cross, width = 16, height = 12, units = "in", dpi = 1200, bg = "white")
ggsave("figures/14_causal/dag_cross_group.svg", p_cross, width = 16, height = 12, units = "in", bg = "white")

# 5) 边强度矩阵
cat("绘制: 边强度矩阵 ...\n")
edges_mat <- edges_strength %>% dplyr::filter(strength >= min_strength_matrix) %>% dplyr::select(from,to,strength)
nodes_all <- unique(c(edges_mat$from, edges_mat$to))
adj <- matrix(0, nrow = length(nodes_all), ncol = length(nodes_all), dimnames = list(nodes_all, nodes_all))
if (nrow(edges_mat) > 0) {
  for (i in seq_len(nrow(edges_mat))) adj[edges_mat$from[i], edges_mat$to[i]] <- edges_mat$strength[i]
}
adj_long <- reshape2::melt(adj) %>% dplyr::filter(value > 0) %>% dplyr::rename(from = Var1, to = Var2, strength = value)
p_mat <- ggplot(adj_long, aes(x = to, y = from, fill = strength)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(low = "white", mid = "lightblue", high = "darkblue", midpoint = 0.65, limits = c(min_strength_matrix, 1), name = "Edge strength") +
  labs(title = "Causal Edge Strength Matrix", subtitle = paste0("Directed edges ≥ ", min_strength_matrix)) +
  theme_minimal(base_family = "Arial", base_size = 9) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
        axis.text.y = element_text(size = 7),
        axis.title  = element_text(size = 9, face = "bold"),
        plot.title  = element_text(face = "bold", hjust = 0.5),
        plot.subtitle= element_text(color = "grey30", hjust = 0.5),
        legend.position = "right")
ggsave("figures/14_causal/dag_strength_matrix.png", p_mat, width = 12, height = 11, units = "in", dpi = 1200, bg = "white")
ggsave("figures/14_causal/dag_strength_matrix.svg", p_mat, width = 12, height = 11, units = "in", bg = "white")

# 6) PC-DAG
if (!is.null(pc_fit)) {
  cat("绘制: PC-DAG ...\n")
  amat <- as(pc_fit@graph, "matrix")
  g_pc <- igraph::graph_from_adjacency_matrix(amat, mode = "directed", diag = FALSE)
  # 合并分组（按名称匹配）
  vert_pc <- var_groups %>% dplyr::select(name, group_label)
  g_pc <- igraph::set_vertex_attr(g_pc, "group_label", index = V(g_pc), value = vert_pc$group_label[match(V(g_pc)$name, vert_pc$name)])
  tg_pc <- tidygraph::as_tbl_graph(g_pc) %>%
    tidygraph::activate(nodes) %>%
    dplyr::mutate(group_label = factor(group_label, levels = names(group_colors)))
  p_pc <- ggraph::ggraph(tg_pc, layout = "sugiyama") +
    ggraph::geom_edge_link(arrow = grid::arrow(length = grid::unit(2, "mm"), type = "closed"),
                           end_cap = ggraph::circle(2, "mm"), edge_colour = "grey35", lineend = "round", width = 0.3) +
    ggraph::geom_node_point(aes(fill = group_label), shape = 21, size = 2.4, stroke = 0.4, colour = "black") +
    scale_fill_manual(values = group_colors, name = "Variable group") +
    ggraph::geom_node_text(aes(label = name), size = 2.6, family = "Arial", vjust = -0.8) +
    labs(title = "PC-DAG") +
    theme_void(base_family = "Arial")
  ggsave("figures/14_causal/dag_pc_layout.png", p_pc, width = 10, height = 7.5, units = "in", dpi = 1200, bg = "white")
  ggsave("figures/14_causal/dag_pc_layout.svg", p_pc, width = 10, height = 7.5, units = "in", bg = "white")
} else {
  cat("跳过PC-DAG（无graph_pc.rds）\n")
}

cat("\n✓ 全部因果DAG图件已生成：figures/14_causal/ 下查看。\n\n")


