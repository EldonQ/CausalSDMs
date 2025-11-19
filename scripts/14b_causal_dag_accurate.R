#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 14b_causal_dag_accurate.R
# 功能说明: 基于实际47变量分组绘制精确的因果DAG网络图
# 改进点: 
#   1. 从variables_selected_47.csv读取真实分组
#   2. 边粗细精确映射bootstrap强度
#   3. 四组变量准确着色
#   4. 力导向布局避免节点重叠
#   5. 使用ggrepel避免标签重叠
# 输出: 4张专业Nature级别因果网络图
# ==============================================================================

rm(list = ls())
gc()
setwd("E:/SDM01")

options(repos = c(CRAN = "https://mirrors.sustech.edu.cn/CRAN/"))

# 加载必需的包
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
cat("精确因果DAG网络图绘制\n")
cat("基于47变量真实分组\n")
cat("======================================\n\n")

# ========================================
# 第1步: 读取真实变量分组
# ========================================
cat("步骤 1/5: 读取变量分组信息...\n")

var_groups <- read.csv("scripts/variables_selected_47.csv", stringsAsFactors = FALSE)
var_groups <- var_groups %>%
  filter(!is.na(variable) & variable != "") %>%
  select(variable, category) %>%
  rename(var = variable, group = category)

# 创建更友好的组名
var_groups <- var_groups %>%
  mutate(
    group_label = case_when(
      group == "G1_TopoSlopeFlow" ~ "Topography & Flow",
      group == "G2_Hydroclim_wavg" ~ "Hydroclimatic",
      group == "G3_Landcover_wavg" ~ "Land Cover",
      group == "G4_Soil_wavg" ~ "Soil Properties",
      TRUE ~ "Other"
    )
  )

cat(sprintf("  ✓ 读取到 %d 个变量\n", nrow(var_groups)))
cat("  分组统计:\n")
group_summary <- var_groups %>% 
  count(group_label) %>% 
  arrange(desc(n))
for(i in 1:nrow(group_summary)) {
  cat(sprintf("    - %s: %d 变量\n", 
              group_summary$group_label[i], 
              group_summary$n[i]))
}

# ========================================
# 第2步: 读取边强度数据并校验
# ========================================
cat("\n步骤 2/5: 读取边强度数据...\n")

edges_data <- read.csv("output/14_causal/edges_summary.csv", stringsAsFactors = FALSE)

# 校验所有边的变量是否在变量清单中
all_edge_vars <- unique(c(edges_data$from, edges_data$to))
missing_vars <- setdiff(all_edge_vars, var_groups$var)
if(length(missing_vars) > 0) {
  cat("  警告: 以下变量在edges但不在变量清单中:\n")
  cat(paste("   ", missing_vars, collapse = "\n"), "\n")
}

cat(sprintf("  ✓ 边数据: %d 条边，强度范围 [%.3f, %.3f]\n", 
            nrow(edges_data), 
            min(edges_data$strength), 
            max(edges_data$strength)))

# ========================================
# 第3步: 配置专业配色方案（Nature风格）
# ========================================
cat("\n步骤 3/5: 配置配色方案...\n")

# Nature期刊风格配色（饱和度适中，对比明显）
group_colors <- c(
  "Topography & Flow" = "#E41A1C",      # 红色 - 地形流向
  "Hydroclimatic" = "#377EB8",          # 蓝色 - 水文气候
  "Land Cover" = "#4DAF4A",             # 绿色 - 土地覆盖
  "Soil Properties" = "#984EA3"         # 紫色 - 土壤属性
)

cat("  配色方案:\n")
for(gn in names(group_colors)) {
  cat(sprintf("    %s: %s\n", gn, group_colors[gn]))
}

# ========================================
# 图1: 完整HC-DAG网络（strength ≥ 0.55）
# ========================================
cat("\n步骤 4/5: 绘制图1 - 完整网络...\n")

edges_strong <- edges_data %>% 
  filter(strength >= 0.55) %>%
  arrange(desc(strength))

cat(sprintf("  筛选强边 (≥0.55): %d 条\n", nrow(edges_strong)))

if(nrow(edges_strong) > 0) {
  # 构建图对象
  g1 <- igraph::graph_from_data_frame(
    d = edges_strong, 
    directed = TRUE, 
    vertices = var_groups %>% select(var, group_label)
  )
  
  # 转换为tidygraph
  tg1 <- tidygraph::as_tbl_graph(g1) %>%
    activate(nodes) %>%
    mutate(
      degree = centrality_degree(mode = "all"),
      group_label = factor(group_label, levels = names(group_colors))
    ) %>%
    activate(edges) %>%
    mutate(strength = strength)
  
  # 计算每组节点数（用于图例说明）
  node_counts <- tg1 %>%
    activate(nodes) %>%
    as_tibble() %>%
    count(group_label)
  
  # 使用FR力导向布局
  set.seed(42)
  p1 <- ggraph(tg1, layout = "fr") +
    # 边：宽度和透明度映射强度
    geom_edge_link(
      aes(width = strength, alpha = strength),
      arrow = arrow(length = unit(3, "mm"), type = "closed"),
      end_cap = circle(3.5, "mm"),
      color = "grey30",
      lineend = "round"
    ) +
    scale_edge_width_continuous(
      range = c(0.4, 3.0),
      name = "Edge strength\n(bootstrap)",
      breaks = c(0.55, 0.7, 0.85, 1.0)
    ) +
    scale_edge_alpha_continuous(
      range = c(0.35, 0.95),
      guide = "none"
    ) +
    # 节点：按组着色，大小反映度中心性
    geom_node_point(
      aes(fill = group_label, size = degree),
      shape = 21,
      color = "black",
      stroke = 0.5
    ) +
    scale_fill_manual(
      values = group_colors,
      name = "Variable group",
      labels = function(x) {
        sapply(x, function(g) {
          n <- node_counts$n[node_counts$group_label == g]
          if(length(n) > 0) {
            sprintf("%s (n=%d)", g, n)
          } else {
            g
          }
        })
      }
    ) +
    scale_size_continuous(
      range = c(3, 9),
      name = "Node degree"
    ) +
    # 标签：使用ggrepel避免重叠
    geom_node_text(
      aes(label = name),
      size = 2.8,
      family = "Arial",
      repel = TRUE,
      point.padding = unit(0.4, "lines"),
      box.padding = unit(0.4, "lines"),
      segment.color = "grey60",
      segment.size = 0.25,
      max.overlaps = 50,
      force = 2,
      force_pull = 0.5
    ) +
    labs(
      title = "Causal Network Structure: HC-DAG with Bootstrap Validation",
      subtitle = sprintf("47 environmental variables, %d stable edges (strength ≥ 0.55), 300 bootstrap replicates", 
                        nrow(edges_strong))
    ) +
    theme_void(base_family = "Arial", base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey20", margin = margin(b = 15)),
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      legend.key.size = unit(0.8, "lines"),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  # 保存高分辨率图件
  ggsave(
    "figures/14_causal/dag_hc_complete_accurate.png",
    plot = p1,
    width = 16,
    height = 12,
    units = "in",
    dpi = 1200,
    bg = "white"
  )
  
  ggsave(
    "figures/14_causal/dag_hc_complete_accurate.svg",
    plot = p1,
    width = 16,
    height = 12,
    units = "in",
    bg = "white"
  )
  
  cat("  ✓ 已保存: dag_hc_complete_accurate.png/svg\n")
} else {
  cat("  警告: 无强度≥0.55的边\n")
}

# ========================================
# 图2: 核心路径网络（top 50强边）
# ========================================
cat("\n步骤 4/5: 绘制图2 - 核心路径...\n")

edges_top50 <- edges_data %>%
  arrange(desc(strength)) %>%
  head(50)

cat(sprintf("  筛选top 50边，强度范围 [%.3f, %.3f]\n", 
            min(edges_top50$strength), 
            max(edges_top50$strength)))

g2 <- igraph::graph_from_data_frame(
  d = edges_top50,
  directed = TRUE,
  vertices = var_groups %>% select(var, group_label)
)

tg2 <- tidygraph::as_tbl_graph(g2) %>%
  activate(nodes) %>%
  mutate(
    degree = centrality_degree(mode = "all"),
    betweenness = centrality_betweenness(directed = TRUE),
    group_label = factor(group_label, levels = names(group_colors))
  ) %>%
  activate(edges) %>%
  mutate(strength = strength)

# 使用KK布局
set.seed(123)
p2 <- ggraph(tg2, layout = "kk") +
  geom_edge_arc(
    aes(width = strength, color = strength),
    arrow = arrow(length = unit(3, "mm"), type = "closed"),
    end_cap = circle(3.5, "mm"),
    strength = 0.12,
    lineend = "round"
  ) +
  scale_edge_width_continuous(
    range = c(0.5, 3.5),
    name = "Edge strength"
  ) +
  scale_edge_color_gradient(
    low = "grey70",
    high = "grey10",
    name = "Edge strength"
  ) +
  geom_node_point(
    aes(fill = group_label, size = betweenness),
    shape = 21,
    color = "black",
    stroke = 0.6
  ) +
  scale_fill_manual(
    values = group_colors,
    name = "Variable group"
  ) +
  scale_size_continuous(
    range = c(4, 11),
    name = "Betweenness\ncentrality"
  ) +
  geom_node_text(
    aes(label = name),
    size = 3.2,
    family = "Arial",
    repel = TRUE,
    point.padding = unit(0.5, "lines"),
    box.padding = unit(0.5, "lines"),
    segment.color = "grey50",
    segment.size = 0.3,
    max.overlaps = 60,
    fontface = "bold",
    force = 3
  ) +
  labs(
    title = "Core Causal Pathways (Top 50 Strongest Edges)",
    subtitle = "Node size reflects betweenness centrality; edge width and color reflect bootstrap strength"
  ) +
  theme_void(base_family = "Arial", base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey20", margin = margin(b = 15)),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  "figures/14_causal/dag_core_pathways_accurate.png",
  plot = p2,
  width = 16,
  height = 12,
  units = "in",
  dpi = 1200,
  bg = "white"
)

ggsave(
  "figures/14_causal/dag_core_pathways_accurate.svg",
  plot = p2,
  width = 16,
  height = 12,
  units = "in",
  bg = "white"
)

cat("  ✓ 已保存: dag_core_pathways_accurate.png/svg\n")

# ========================================
# 图3: 跨组因果关系（组间耦合）
# ========================================
cat("\n步骤 5/5: 绘制图3 - 跨组因果关系...\n")

# 为边添加分组信息
edges_with_groups <- edges_top50 %>%
  left_join(var_groups %>% select(var, group_label), by = c("from" = "var")) %>%
  rename(from_group = group_label) %>%
  left_join(var_groups %>% select(var, group_label), by = c("to" = "var")) %>%
  rename(to_group = group_label)

# 筛选跨组边
edges_cross <- edges_with_groups %>%
  filter(from_group != to_group) %>%
  arrange(desc(strength))

cat(sprintf("  跨组边数量: %d\n", nrow(edges_cross)))

if(nrow(edges_cross) > 0) {
  # 取前30条跨组边
  edges_cross_top <- edges_cross %>% head(30)
  
  g3 <- igraph::graph_from_data_frame(
    d = edges_cross_top,
    directed = TRUE,
    vertices = var_groups %>% select(var, group_label)
  )
  
  tg3 <- tidygraph::as_tbl_graph(g3) %>%
    activate(nodes) %>%
    mutate(
      degree = centrality_degree(mode = "all"),
      group_label = factor(group_label, levels = names(group_colors))
    ) %>%
    activate(edges) %>%
    mutate(strength = strength)
  
  set.seed(456)
  p3 <- ggraph(tg3, layout = "graphopt") +
    geom_edge_link(
      aes(width = strength, alpha = strength),
      arrow = arrow(length = unit(3, "mm"), type = "closed"),
      end_cap = circle(3.5, "mm"),
      color = "grey25",
      lineend = "round"
    ) +
    scale_edge_width_continuous(
      range = c(0.6, 3.5),
      name = "Edge strength"
    ) +
    scale_edge_alpha_continuous(
      range = c(0.4, 1),
      guide = "none"
    ) +
    geom_node_point(
      aes(fill = group_label, size = degree),
      shape = 21,
      color = "black",
      stroke = 0.6
    ) +
    scale_fill_manual(
      values = group_colors,
      name = "Variable group"
    ) +
    scale_size_continuous(
      range = c(4, 10),
      name = "Node degree"
    ) +
    geom_node_text(
      aes(label = name),
      size = 3,
      family = "Arial",
      repel = TRUE,
      point.padding = unit(0.4, "lines"),
      box.padding = unit(0.4, "lines"),
      segment.color = "grey60",
      segment.size = 0.25,
      max.overlaps = 50,
      force = 2.5
    ) +
    labs(
      title = "Cross-Domain Causal Relationships",
      subtitle = "Top 30 inter-group edges showing coupling between environmental domains"
    ) +
    theme_void(base_family = "Arial", base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 5)),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey20", margin = margin(b = 15)),
      legend.position = "right",
      legend.box = "vertical",
      legend.title = element_text(size = 10, face = "bold"),
      legend.text = element_text(size = 9),
      plot.margin = margin(10, 10, 10, 10)
    )
  
  ggsave(
    "figures/14_causal/dag_cross_domain_accurate.png",
    plot = p3,
    width = 16,
    height = 12,
    units = "in",
    dpi = 1200,
    bg = "white"
  )
  
  ggsave(
    "figures/14_causal/dag_cross_domain_accurate.svg",
    plot = p3,
    width = 16,
    height = 12,
    units = "in",
    bg = "white"
  )
  
  cat("  ✓ 已保存: dag_cross_domain_accurate.png/svg\n")
}

# ========================================
# 图4: 分组统计与边强度分布
# ========================================
cat("\n生成补充图4: 边强度统计...\n")

# 按来源组和目标组统计边数和平均强度
edge_group_summary <- edges_with_groups %>%
  filter(strength >= 0.3) %>%
  group_by(from_group, to_group) %>%
  summarise(
    n_edges = n(),
    mean_strength = mean(strength),
    median_strength = median(strength),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_strength))

# 绘制组间边强度热图
p4 <- ggplot(edge_group_summary, 
             aes(x = to_group, y = from_group, fill = mean_strength, size = n_edges)) +
  geom_point(shape = 21, color = "black", stroke = 0.5) +
  scale_fill_gradient2(
    low = "white",
    mid = "lightblue",
    high = "darkblue",
    midpoint = 0.65,
    limits = c(0.3, 1),
    name = "Mean edge\nstrength"
  ) +
  scale_size_continuous(
    range = c(3, 15),
    name = "Number\nof edges"
  ) +
  geom_text(
    aes(label = n_edges),
    size = 3,
    family = "Arial",
    fontface = "bold",
    color = "white"
  ) +
  labs(
    title = "Inter-Group Causal Relationship Summary",
    subtitle = "Bubble size = number of edges; color = mean bootstrap strength",
    x = "Target variable group",
    y = "Source variable group"
  ) +
  theme_minimal(base_family = "Arial", base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    axis.title = element_text(size = 11, face = "bold"),
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey20", margin = margin(b = 15)),
    legend.position = "right",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    plot.margin = margin(10, 10, 10, 10)
  )

ggsave(
  "figures/14_causal/dag_group_summary_accurate.png",
  plot = p4,
  width = 10,
  height = 8,
  units = "in",
  dpi = 1200,
  bg = "white"
)

ggsave(
  "figures/14_causal/dag_group_summary_accurate.svg",
  plot = p4,
  width = 10,
  height = 8,
  units = "in",
  bg = "white"
)

cat("  ✓ 已保存: dag_group_summary_accurate.png/svg\n")

# ========================================
# 保存变量分组信息供参考
# ========================================
write.csv(var_groups, 
          "output/14_causal/variable_groups_validated.csv", 
          row.names = FALSE)

cat("\n======================================\n")
cat("精确DAG图绘制完成！\n")
cat("======================================\n\n")

cat("输出文件:\n")
cat("1. dag_hc_complete_accurate.png/svg\n")
cat("   - 完整网络 (strength≥0.55)\n")
cat("   - 4组变量准确着色\n")
cat("   - 边粗细映射bootstrap强度\n\n")

cat("2. dag_core_pathways_accurate.png/svg\n")
cat("   - 核心路径 (top 50边)\n")
cat("   - 节点大小映射介数中心性\n")
cat("   - 边颜色+粗细双重映射强度\n\n")

cat("3. dag_cross_domain_accurate.png/svg\n")
cat("   - 跨域因果关系 (top 30跨组边)\n")
cat("   - 突出组间耦合\n\n")

cat("4. dag_group_summary_accurate.png/svg\n")
cat("   - 组间边统计汇总\n")
cat("   - 气泡图展示边数与平均强度\n\n")

cat("变量分组验证:\n")
for(i in 1:nrow(group_summary)) {
  cat(sprintf("  %s: %d 变量\n", 
              group_summary$group_label[i], 
              group_summary$n[i]))
}
cat("\n")

