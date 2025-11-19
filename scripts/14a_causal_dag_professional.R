#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 14a_causal_dag_professional.R
# 功能说明: 专业绘制因果DAG网络图，优化布局与边权重可视化
# 改进点: 
#   1. 使用力导向布局避免节点重叠
#   2. 边粗细映射稳定性权重
#   3. 按变量类型分组着色
#   4. 使用ggrepel避免标签重叠
#   5. 创建高清晰度、可读性强的图件
# ==============================================================================

rm(list = ls())
gc()
setwd("E:/SDM01")

options(repos = c(CRAN = "https://mirrors.sustech.edu.cn/CRAN/"))

# 加载包
packages <- c("tidyverse", "igraph", "ggraph", "tidygraph", 
              "ggrepel", "RColorBrewer", "sysfonts", "showtext")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 字体设置
try({
  sysfonts::font_add(
    family = "Arial",
    regular = "C:/Windows/Fonts/arial.ttf",
    bold = "C:/Windows/Fonts/arialbd.ttf"
  )
  showtext::showtext_opts(dpi = 1200)
  showtext::showtext_auto(enable = TRUE)
}, silent = TRUE)

cat("\n======================================\n")
cat("专业因果DAG网络图绘制\n")
cat("======================================\n\n")

# 读取边强度数据
edges_data <- read.csv("output/14_causal/edges_summary.csv")

# 定义变量分组（基于变量名前缀）
variable_groups <- data.frame(
  var = c(),
  group = c(),
  stringsAsFactors = FALSE
)

all_vars <- unique(c(edges_data$from, edges_data$to))

for(v in all_vars) {
  if(grepl("^hydro_", v)) {
    variable_groups <- rbind(variable_groups, 
                            data.frame(var = v, group = "Hydroclimatic", stringsAsFactors = FALSE))
  } else if(grepl("^(dem_|slope_|tpi_|tri_)", v)) {
    variable_groups <- rbind(variable_groups, 
                            data.frame(var = v, group = "Topographic", stringsAsFactors = FALSE))
  } else if(grepl("^lc_", v)) {
    variable_groups <- rbind(variable_groups, 
                            data.frame(var = v, group = "Land cover", stringsAsFactors = FALSE))
  } else if(grepl("^soil_", v)) {
    variable_groups <- rbind(variable_groups, 
                            data.frame(var = v, group = "Soil", stringsAsFactors = FALSE))
  } else if(grepl("^flow_", v)) {
    variable_groups <- rbind(variable_groups, 
                            data.frame(var = v, group = "Hydrological", stringsAsFactors = FALSE))
  } else {
    variable_groups <- rbind(variable_groups, 
                            data.frame(var = v, group = "Other", stringsAsFactors = FALSE))
  }
}

# ========================================
# 图1: HC-DAG 完整网络（强边，threshold ≥ 0.55）
# ========================================
cat("绘制图1: HC-DAG (strength ≥ 0.55)...\n")

edges_strong <- edges_data %>% 
  filter(strength >= 0.55) %>%
  arrange(desc(strength))

if(nrow(edges_strong) > 0) {
  # 构建图对象
  g1 <- igraph::graph_from_data_frame(
    d = edges_strong, 
    directed = TRUE, 
    vertices = variable_groups
  )
  
  # 转换为tidygraph对象
  tg1 <- tidygraph::as_tbl_graph(g1) %>%
    activate(nodes) %>%
    mutate(
      degree = centrality_degree(mode = "all"),
      group = factor(group)
    ) %>%
    activate(edges) %>%
    mutate(strength = strength)
  
  # 配色方案
  group_colors <- c(
    "Hydroclimatic" = "#E41A1C",
    "Topographic" = "#377EB8",
    "Land cover" = "#4DAF4A",
    "Soil" = "#984EA3",
    "Hydrological" = "#FF7F00",
    "Other" = "#999999"
  )
  
  # 使用Fruchterman-Reingold力导向布局
  p1 <- ggraph(tg1, layout = "fr") +
    # 绘制边，宽度映射强度
    geom_edge_link(
      aes(width = strength, alpha = strength),
      arrow = arrow(length = unit(3, "mm"), type = "closed"),
      end_cap = circle(4, "mm"),
      color = "grey40",
      lineend = "round"
    ) +
    scale_edge_width_continuous(
      range = c(0.3, 2.5),
      name = "Edge\nstrength"
    ) +
    scale_edge_alpha_continuous(
      range = c(0.3, 0.9),
      guide = "none"
    ) +
    # 绘制节点
    geom_node_point(
      aes(fill = group, size = degree),
      shape = 21,
      color = "black",
      stroke = 0.5
    ) +
    scale_fill_manual(
      values = group_colors,
      name = "Variable\ntype"
    ) +
    scale_size_continuous(
      range = c(3, 8),
      name = "Node\ndegree"
    ) +
    # 绘制标签（使用ggrepel避免重叠）
    geom_node_text(
      aes(label = name),
      size = 2.5,
      family = "Arial",
      repel = TRUE,
      point.padding = unit(0.3, "lines"),
      box.padding = unit(0.3, "lines"),
      segment.color = "grey70",
      segment.size = 0.2,
      max.overlaps = 30
    ) +
    labs(
      title = "Causal Network: HC-DAG (Bootstrap Strength ≥ 0.55)",
      subtitle = paste0("300 bootstrap replicates, ", nrow(edges_strong), " stable edges")
    ) +
    theme_void(base_family = "Arial", base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # 保存高分辨率图件
  ggsave(
    "figures/14_causal/dag_hc_professional_v2.png",
    plot = p1,
    width = 14,
    height = 10,
    units = "in",
    dpi = 1200,
    bg = "white"
  )
  
  ggsave(
    "figures/14_causal/dag_hc_professional_v2.svg",
    plot = p1,
    width = 14,
    height = 10,
    units = "in",
    bg = "white"
  )
  
  cat("  ✓ 已保存: dag_hc_professional_v2.png/svg\n")
} else {
  cat("  警告: 无边强度 ≥ 0.55 的边\n")
}

# ========================================
# 图2: 核心因果路径（top 50强边，简化视图）
# ========================================
cat("绘制图2: 核心因果路径 (top 50 edges)...\n")

edges_core <- edges_data %>%
  arrange(desc(strength)) %>%
  head(50)

g2 <- igraph::graph_from_data_frame(
  d = edges_core,
  directed = TRUE,
  vertices = variable_groups
)

tg2 <- tidygraph::as_tbl_graph(g2) %>%
  activate(nodes) %>%
  mutate(
    degree = centrality_degree(mode = "all"),
    betweenness = centrality_betweenness(directed = TRUE),
    group = factor(group)
  ) %>%
  activate(edges) %>%
  mutate(strength = strength)

# 使用Kamada-Kawai布局（适合中等规模网络）
p2 <- ggraph(tg2, layout = "kk") +
  geom_edge_arc(
    aes(width = strength, alpha = strength),
    arrow = arrow(length = unit(3, "mm"), type = "closed"),
    end_cap = circle(4, "mm"),
    color = "grey30",
    strength = 0.1,  # 轻微弯曲避免重叠
    lineend = "round"
  ) +
  scale_edge_width_continuous(
    range = c(0.5, 3),
    name = "Edge\nstrength"
  ) +
  scale_edge_alpha_continuous(
    range = c(0.4, 1),
    guide = "none"
  ) +
  geom_node_point(
    aes(fill = group, size = betweenness),
    shape = 21,
    color = "black",
    stroke = 0.6
  ) +
  scale_fill_manual(
    values = group_colors,
    name = "Variable\ntype"
  ) +
  scale_size_continuous(
    range = c(4, 10),
    name = "Betweenness\ncentrality"
  ) +
  geom_node_text(
    aes(label = name),
    size = 3,
    family = "Arial",
    repel = TRUE,
    point.padding = unit(0.4, "lines"),
    box.padding = unit(0.4, "lines"),
    segment.color = "grey60",
    segment.size = 0.3,
    max.overlaps = 50,
    fontface = "bold"
  ) +
  labs(
    title = "Core Causal Pathways (Top 50 Strongest Edges)",
    subtitle = "Node size reflects betweenness centrality; edge width reflects bootstrap strength"
  ) +
  theme_void(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  "figures/14_causal/dag_core_pathways_v2.png",
  plot = p2,
  width = 14,
  height = 10,
  units = "in",
  dpi = 1200,
  bg = "white"
)

ggsave(
  "figures/14_causal/dag_core_pathways_v2.svg",
  plot = p2,
  width = 14,
  height = 10,
  units = "in",
  bg = "white"
)

cat("  ✓ 已保存: dag_core_pathways_v2.png/svg\n")

# ========================================
# 图3: 按变量组分面的子网络
# ========================================
cat("绘制图3: 分组子网络...\n")

# 筛选组间边（跨组因果关系）
edges_cross_group <- edges_core %>%
  left_join(variable_groups, by = c("from" = "var")) %>%
  rename(from_group = group) %>%
  left_join(variable_groups, by = c("to" = "var")) %>%
  rename(to_group = group) %>%
  filter(from_group != to_group) %>%
  arrange(desc(strength))

if(nrow(edges_cross_group) > 0) {
  # 创建简化的跨组关系图
  edges_cross_top <- edges_cross_group %>% head(30)
  
  g3 <- igraph::graph_from_data_frame(
    d = edges_cross_top,
    directed = TRUE,
    vertices = variable_groups
  )
  
  tg3 <- tidygraph::as_tbl_graph(g3) %>%
    activate(nodes) %>%
    mutate(
      degree = centrality_degree(mode = "all"),
      group = factor(group)
    ) %>%
    activate(edges) %>%
    mutate(strength = strength)
  
  p3 <- ggraph(tg3, layout = "graphopt") +
    geom_edge_link(
      aes(width = strength, color = strength),
      arrow = arrow(length = unit(3, "mm"), type = "closed"),
      end_cap = circle(4, "mm"),
      lineend = "round"
    ) +
    scale_edge_width_continuous(
      range = c(0.5, 3),
      name = "Edge\nstrength"
    ) +
    scale_edge_color_gradient(
      low = "grey70",
      high = "grey10",
      name = "Edge\nstrength"
    ) +
    geom_node_point(
      aes(fill = group, size = degree),
      shape = 21,
      color = "black",
      stroke = 0.6
    ) +
    scale_fill_manual(
      values = group_colors,
      name = "Variable\ntype"
    ) +
    scale_size_continuous(
      range = c(4, 9),
      name = "Node\ndegree"
    ) +
    geom_node_text(
      aes(label = name),
      size = 2.8,
      family = "Arial",
      repel = TRUE,
      point.padding = unit(0.3, "lines"),
      box.padding = unit(0.3, "lines"),
      segment.color = "grey70",
      segment.size = 0.2,
      max.overlaps = 40
    ) +
    labs(
      title = "Cross-Domain Causal Relationships",
      subtitle = "Top 30 inter-group edges showing coupling between environmental domains"
    ) +
    theme_void(base_family = "Arial", base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey30"),
      legend.position = "right",
      legend.title = element_text(size = 9, face = "bold"),
      legend.text = element_text(size = 8),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  ggsave(
    "figures/14_causal/dag_cross_domain_v2.png",
    plot = p3,
    width = 14,
    height = 10,
    units = "in",
    dpi = 1200,
    bg = "white"
  )
  
  ggsave(
    "figures/14_causal/dag_cross_domain_v2.svg",
    plot = p3,
    width = 14,
    height = 10,
    units = "in",
    bg = "white"
  )
  
  cat("  ✓ 已保存: dag_cross_domain_v2.png/svg\n")
}

# ========================================
# 图4: 边强度分布热图
# ========================================
cat("绘制图4: 边强度矩阵热图...\n")

# 创建邻接矩阵（仅显示强边）
edges_matrix <- edges_data %>%
  filter(strength >= 0.3) %>%
  select(from, to, strength)

# 获取所有涉及的变量
all_nodes <- unique(c(edges_matrix$from, edges_matrix$to))

# 创建完整的邻接矩阵
adj_matrix <- matrix(0, nrow = length(all_nodes), ncol = length(all_nodes))
rownames(adj_matrix) <- all_nodes
colnames(adj_matrix) <- all_nodes

for(i in 1:nrow(edges_matrix)) {
  adj_matrix[edges_matrix$from[i], edges_matrix$to[i]] <- edges_matrix$strength[i]
}

# 转换为长格式用于ggplot
adj_long <- reshape2::melt(adj_matrix) %>%
  filter(value > 0) %>%
  rename(from = Var1, to = Var2, strength = value)

# 添加变量组信息
adj_long <- adj_long %>%
  left_join(variable_groups, by = c("from" = "var")) %>%
  rename(from_group = group) %>%
  left_join(variable_groups, by = c("to" = "var")) %>%
  rename(to_group = group)

# 绘制热图
p4 <- ggplot(adj_long, aes(x = to, y = from, fill = strength)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_gradient2(
    low = "white",
    mid = "lightblue",
    high = "darkblue",
    midpoint = 0.65,
    limits = c(0.3, 1),
    name = "Edge\nstrength"
  ) +
  labs(
    title = "Causal Edge Strength Matrix",
    subtitle = "Showing directed edges with bootstrap strength ≥ 0.3",
    x = "Target variable (to)",
    y = "Source variable (from)"
  ) +
  theme_minimal(base_family = "Arial", base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title = element_text(size = 9, face = "bold"),
    plot.title = element_text(face = "bold", size = 12, hjust = 0.5),
    plot.subtitle = element_text(size = 9, hjust = 0.5, color = "grey30"),
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text = element_text(size = 8),
    panel.grid = element_blank(),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  "figures/14_causal/dag_strength_matrix_v2.png",
  plot = p4,
  width = 12,
  height = 11,
  units = "in",
  dpi = 1200,
  bg = "white"
)

ggsave(
  "figures/14_causal/dag_strength_matrix_v2.svg",
  plot = p4,
  width = 12,
  height = 11,
  units = "in",
  bg = "white"
)

cat("  ✓ 已保存: dag_strength_matrix_v2.png/svg\n")

cat("\n======================================\n")
cat("专业DAG图绘制完成！\n")
cat("======================================\n\n")
cat("输出文件:\n")
cat("1. dag_hc_professional_v2.png/svg - 完整网络(strength≥0.55)\n")
cat("2. dag_core_pathways_v2.png/svg - 核心路径(top 50)\n")
cat("3. dag_cross_domain_v2.png/svg - 跨域因果关系\n")
cat("4. dag_strength_matrix_v2.png/svg - 强度矩阵热图\n\n")

