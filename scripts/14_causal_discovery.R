#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 14_causal_discovery.R
# 功能说明: 基于约束与评分结合的方法进行因果结构学习（DAG），并输出图件与表格
# 方法: PC算法 (pcalg) + 评分驱动HC (bnlearn)，子样本稳定性评估
# 输入文件: output/04_collinearity/collinearity_removed.csv
# 输出文件: output/14_causal/graph_pc.rds, graph_ges.rds, edges_summary.csv
#          figures/14_causal/dag_pc.png, dag_ges.png
#          figures/14_causal/edge_stability.png
# 作者: Nature级别科研项目
# 日期: 2025-10-24
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/SDM01")

# 设定 CRAN 镜像，避免交互式选择
options(repos = c(CRAN = "https://mirrors.sustech.edu.cn/CRAN/"))

# 先安装/加载 Bioconductor 依赖（graph/RBGL/Rgraphviz）
try({
  if(!require("BiocManager", character.only = TRUE)) {
    install.packages("BiocManager")
    library(BiocManager)
  }
  for(bpkg in c("graph", "RBGL", "Rgraphviz")) {
    if(!require(bpkg, character.only = TRUE)) {
      BiocManager::install(bpkg, ask = FALSE, update = FALSE)
      library(bpkg, character.only = TRUE)
    }
  }
}, silent = TRUE)

# 加载必要的 CRAN 包（中文注释：全部使用英文标注出图，Arial字体，1200dpi；DAG绘图优先专业方案）
packages <- c("tidyverse", "bnlearn", "pcalg", "igraph", "ggraph", "tidygraph", "ggrepel", "sysfonts", "showtext")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

dir.create("output/14_causal", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/14_causal", showWarnings = FALSE, recursive = TRUE)

# 字体设置
try({
  sysfonts::font_add(
    family = "Arial",
    regular = "C:/Windows/Fonts/arial.ttf",
    bold = "C:/Windows/Fonts/arialbd.ttf",
    italic = "C:/Windows/Fonts/ariali.ttf",
    bolditalic = "C:/Windows/Fonts/arialbi.ttf"
  )
  showtext::showtext_opts(dpi = 1200)
  showtext::showtext_auto(enable = TRUE)
}, silent = TRUE)

cat("\n======================================\n")
cat("因果结构学习 (PC + HC)\n")
cat("======================================\n\n")

# 1. 读取数据并准备变量
cat("步骤 1/4: 读取数据...\n")
dat <- read.csv("output/04_collinearity/collinearity_removed.csv")
exclude_cols <- c("id", "species", "lon", "lat", "source", "presence", "presence.1")
env_vars <- setdiff(colnames(dat), exclude_cols)

X <- dat[, env_vars, drop = FALSE]

# 连续变量标准化（PC算法更稳健）
X_scaled <- as.data.frame(scale(X))
X_scaled[is.na(X_scaled)] <- 0

# 2. PC算法 (pcalg)
cat("步骤 2/4: PC 算法...\n")
## 估计协方差与样本量
S <- cor(X_scaled)
pc_fit <- pcalg::pc(suffStat = list(C = S, n = nrow(X_scaled)), indepTest = pcalg::gaussCItest,
                    alpha = 0.01, labels = colnames(X_scaled))
saveRDS(pc_fit, file = "output/14_causal/graph_pc.rds")

# 3. 评分驱动结构学习（HC, BIC-G 对连续数据）
cat("步骤 3/4: 评分驱动结构学习 (HC, BIC-G)...\n")
## 连续数据使用 BIC-G；若失败可回退 BGE
score_fit <- tryCatch({ bnlearn::hc(X_scaled, score = "bic-g") }, error = function(e) {
  bnlearn::hc(X_scaled, score = "bge")
})
saveRDS(score_fit, file = "output/14_causal/graph_hc.rds")

# 4. 边稳定性 (Bootstrap; bnlearn 专业实现)
cat("步骤 4/4: 子样本稳定性 (boot.strength)...\n")
set.seed(20251024)
R <- 300  # 中文注释：重复次数；如需更稳健可提升到500或1000，但会更耗时

# 使用 bnlearn::boot.strength 计算边强度（每次使用约80%样本）
boot_hc <- bnlearn::boot.strength(
  data = X_scaled,
  R = R,
  algorithm = "hc",
  algorithm.args = list(score = "bic-g"),
  m = floor(0.8 * nrow(X_scaled))
)

# 整理与保存边强度表（英文标注，便于后续制图）
edges_strength <- boot_hc %>%
  dplyr::select(from, to, strength) %>%
  dplyr::arrange(dplyr::desc(strength))
write.csv(edges_strength, "output/14_causal/edges_summary.csv", row.names = FALSE)

# 依据阈值构建平均网络（避免过度连线导致信息拥挤）
threshold <- 0.55  # 中文注释：可按需要调整；0.5-0.7 常用
avg_hc <- bnlearn::averaged.network(boot_hc, threshold = threshold)
saveRDS(avg_hc, file = "output/14_causal/graph_hc_avg.rds")

# ==========================================================
# 4a. 基于“真实DAG对象”的专业网络绘制（HC平均网络）
# - 使用 averaged.network 的弧集合（bnlearn::arcs）作为“真实DAG”
# - 将 boot.strength 的强度合并到弧上，用线宽表达强度
# - 读取 47 变量分组信息，按组着色
# - 同时输出完整网络与“最大连通分量(LCC)”两张图，确保网络结构清晰
# ==========================================================
cat("\n附加: 真实DAG(averaged HC) 的专业网络绘制...\n")

# 读取变量分组（真实来源：scripts/variables_selected_47.csv）
var_groups <- tryCatch({
  read.csv("scripts/variables_selected_47.csv", stringsAsFactors = FALSE) %>%
    dplyr::select(variable, category) %>%
    dplyr::rename(name = variable, group = category) %>%
    dplyr::mutate(
      group_label = dplyr::case_when(
        group == "G1_TopoSlopeFlow" ~ "Topography & Flow",
        group == "G2_Hydroclim_wavg" ~ "Hydroclimatic",
        group == "G3_Landcover_wavg" ~ "Land Cover",
        group == "G4_Soil_wavg" ~ "Soil Properties",
        TRUE ~ "Other"
      )
    )
}, error = function(e) {
  # 若读取失败，使用env_vars构造默认分组
  message("WARNING: 未能读取 variables_selected_47.csv, 将使用默认分组。")
  data.frame(name = env_vars, group_label = "Other", stringsAsFactors = FALSE)
})

# 构造平均网络的弧集合（真实DAG边）
arcs_avg <- bnlearn::arcs(avg_hc) %>% as.data.frame(stringsAsFactors = FALSE)
colnames(arcs_avg) <- c("from", "to")

# 将 boot.strength 的强度合并到平均网络的边上（双向合并，取最大值）
strength_df <- boot_hc %>% dplyr::select(from, to, strength)
strength_rev <- strength_df %>% dplyr::transmute(from = to, to = from, strength_rev = strength)
edges_avg <- arcs_avg %>%
  dplyr::left_join(strength_df, by = c("from", "to")) %>%
  dplyr::left_join(strength_rev, by = c("from", "to")) %>%
  dplyr::mutate(strength = dplyr::coalesce(strength, strength_rev, 0.55))

# 构建 igraph 对象，顶点包含分组信息
vertices_df <- var_groups %>% dplyr::select(name, group_label)
g_avg <- igraph::graph_from_data_frame(d = edges_avg, directed = TRUE, vertices = vertices_df)

# 去除孤立点以凸显网络结构（完整网络图仅保留有边的节点）
deg_all <- igraph::degree(g_avg, mode = "all")
g_avg_noniso <- igraph::induced_subgraph(g_avg, vids = which(deg_all > 0))

# 计算最大连通分量（LCC），便于呈现最“网络型”的子图
comp <- igraph::components(igraph::as.undirected(g_avg_noniso))
if (comp$no > 1) {
  lcc_vids <- which(comp$membership == which.max(comp$csize))
  g_avg_lcc <- igraph::induced_subgraph(g_avg_noniso, vids = lcc_vids)
} else {
  g_avg_lcc <- g_avg_noniso
}

# 将图转换为 tidygraph 以便 ggraph 作图
tg_full <- tidygraph::as_tbl_graph(g_avg_noniso) %>%
  tidygraph::activate(nodes) %>%
  dplyr::mutate(
    degree = tidygraph::centrality_degree(mode = "all"),
    group_label = factor(group_label, levels = c("Topography & Flow", "Hydroclimatic", "Land Cover", "Soil Properties", "Other"))
  ) %>%
  tidygraph::activate(edges) %>%
  dplyr::mutate(strength = strength)

tg_lcc <- tidygraph::as_tbl_graph(g_avg_lcc) %>%
  tidygraph::activate(nodes) %>%
  dplyr::mutate(
    degree = tidygraph::centrality_degree(mode = "all"),
    group_label = factor(group_label, levels = c("Topography & Flow", "Hydroclimatic", "Land Cover", "Soil Properties", "Other"))
  ) %>%
  tidygraph::activate(edges) %>%
  dplyr::mutate(strength = strength)

# 组颜色（Nature风格）
group_colors <- c(
  "Topography & Flow" = "#E41A1C",
  "Hydroclimatic" = "#377EB8",
  "Land Cover" = "#4DAF4A",
  "Soil Properties" = "#984EA3",
  "Other" = "#999999"
)

# 专业出图：完整平均网络（去除孤立点）
set.seed(42)
p_avg_full <- ggraph::ggraph(tg_full, layout = "kk") +
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
  labs(title = "Averaged HC-DAG (Non-isolated nodes)", subtitle = paste0("Edges ≥ threshold ", threshold, "; layout=KK")) +
  theme_void(base_family = "Arial")

ggsave("figures/14_causal/dag_hc_avg_network_full.png", plot = p_avg_full,
       width = 14, height = 10, units = "in", dpi = 1200, bg = "white")
ggsave("figures/14_causal/dag_hc_avg_network_full.svg", plot = p_avg_full,
       width = 14, height = 10, units = "in", bg = "white")

# 专业出图：最大连通分量（LCC）
set.seed(43)
p_avg_lcc <- ggraph::ggraph(tg_lcc, layout = "kk") +
  ggraph::geom_edge_link(aes(width = strength, alpha = strength),
                         arrow = grid::arrow(length = grid::unit(3, "mm"), type = "closed"),
                         end_cap = ggraph::circle(3.5, "mm"), colour = "grey30", lineend = "round") +
  ggraph::scale_edge_width_continuous(range = c(0.6, 3.4), name = "Edge strength") +
  ggraph::scale_edge_alpha_continuous(range = c(0.4, 0.95), guide = "none") +
  ggraph::geom_node_point(aes(fill = group_label, size = degree), shape = 21, colour = "black", stroke = 0.6) +
  scale_fill_manual(values = group_colors, name = "Variable group") +
  ggplot2::scale_size_continuous(range = c(3.5, 10), name = "Node degree") +
  ggrepel::geom_text_repel(aes(x = x, y = y, label = name), size = 3.0, family = "Arial",
                           box.padding = grid::unit(0.4, "lines"), point.padding = grid::unit(0.35, "lines"),
                           segment.color = "grey60", segment.size = 0.25, max.overlaps = 100) +
  labs(title = "Averaged HC-DAG (Largest Connected Component)", subtitle = paste0("Edges ≥ threshold ", threshold, "; layout=KK")) +
  theme_void(base_family = "Arial")

ggsave("figures/14_causal/dag_hc_avg_network_lcc.png", plot = p_avg_lcc,
       width = 14, height = 10, units = "in", dpi = 1200, bg = "white")
ggsave("figures/14_causal/dag_hc_avg_network_lcc.svg", plot = p_avg_lcc,
       width = 14, height = 10, units = "in", bg = "white")

# -----------------------------------------
# 专业制图：PC-DAG（ggraph，层次布局）
# -----------------------------------------
amat_pc <- as(pc_fit@graph, "matrix")
g_pc <- igraph::graph_from_adjacency_matrix(amat_pc, mode = "directed", diag = FALSE)
g_pc_tbl <- tidygraph::as_tbl_graph(g_pc)

p_pc <- ggraph(g_pc_tbl, layout = "sugiyama") +
  geom_edge_link(
    arrow = grid::arrow(length = grid::unit(2, "mm"), type = "closed"),
    end_cap = circle(2, "mm"),
    edge_colour = "grey35",
    lineend = "round",
    width = 0.3
  ) +
  geom_node_point(shape = 21, size = 2.2, stroke = 0.3, fill = "white", color = "black") +
  geom_node_text(aes(label = name), size = 2.6, family = "Arial", vjust = -0.8) +
  labs(title = "PC-DAG") +
  theme_void(base_family = "Arial")

ggsave("figures/14_causal/dag_pc.png", plot = p_pc, width = 8, height = 6, units = "in", dpi = 1200, bg = "white")
ggsave("figures/14_causal/dag_pc.svg", plot = p_pc, width = 8, height = 6, units = "in", bg = "white")

# -----------------------------------------
# 专业制图：HC-DAG（边强度可视化）
# 若可用 Rgraphviz，使用 bnlearn::strength.plot；否则回退到 ggraph 宽度映射
# -----------------------------------------
if (requireNamespace("Rgraphviz", quietly = TRUE)) {
  svg("figures/14_causal/dag_hc.svg", width = 8, height = 6, onefile = FALSE, family = "Arial")
  bnlearn::strength.plot(avg_hc, strength = boot_hc, threshold = threshold,
                         shape = "circle", layout = "dot")
  dev.off()
  png("figures/14_causal/dag_hc.png", width = 3000, height = 2400, res = 1200, type = "cairo-png", family = "Arial")
  bnlearn::strength.plot(avg_hc, strength = boot_hc, threshold = threshold,
                         shape = "circle", layout = "dot")
  dev.off()
} else {
  edges_plot <- edges_strength %>% dplyr::filter(strength >= threshold)
  if(nrow(edges_plot) == 0) edges_plot <- edges_strength %>% dplyr::top_n(30, strength)
  g_hc <- igraph::graph_from_data_frame(edges_plot, directed = TRUE, vertices = data.frame(name = env_vars))
  g_hc_tbl <- tidygraph::as_tbl_graph(g_hc)
  p_hc <- ggraph(g_hc_tbl, layout = "sugiyama") +
    geom_edge_link(aes(width = strength, alpha = strength),
                   arrow = grid::arrow(length = grid::unit(2, "mm"), type = "closed"),
                   end_cap = circle(2, "mm"), colour = "grey25", lineend = "round") +
    scale_edge_width(range = c(0.2, 2.6)) +
    scale_edge_alpha(range = c(0.4, 1)) +
    geom_node_point(shape = 21, size = 2.2, stroke = 0.3, fill = "white", color = "black") +
    geom_node_text(aes(label = name), size = 2.6, family = "Arial", vjust = -0.8) +
    labs(title = "HC-DAG (edge strength)") +
    theme_void(base_family = "Arial")
  ggsave("figures/14_causal/dag_hc.png", plot = p_hc, width = 8, height = 6, units = "in", dpi = 1200, bg = "white")
  ggsave("figures/14_causal/dag_hc.svg", plot = p_hc, width = 8, height = 6, units = "in", bg = "white")
}

# -----------------------------------------
# 边稳定性 Top 30（ggplot 专业风格；英文标注；Arial；1200dpi）
# -----------------------------------------
topN <- edges_strength %>% dplyr::top_n(30, strength) %>% dplyr::arrange(strength)
p_bar <- ggplot(topN, aes(x = reorder(paste(from, "->", to), strength), y = strength)) +
  geom_col(fill = "#377EB8", width = 0.7) +
  coord_flip() +
  labs(title = "Edge Stability (Top 30)", x = "", y = "Strength (0-1)") +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid = element_line(color = "grey90", linewidth = 0.2),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    plot.title = element_text(face = "bold")
  )
ggsave("figures/14_causal/edge_stability.png", plot = p_bar, width = 8, height = 6, units = "in", dpi = 1200, bg = "white")
ggsave("figures/14_causal/edge_stability.svg", plot = p_bar, width = 8, height = 6, units = "in", bg = "white")

 

cat("\n======================================\n")
cat("因果结构学习完成\n")
cat("======================================\n\n")

cat("✓ 结果表: output/14_causal/edges_summary.csv\n")
cat("✓ 图件: figures/14_causal/dag_pc.png / dag_hc.png / edge_stability.png\n\n")


