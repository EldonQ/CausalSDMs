#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 15c_venn_diagram.R
# 功能说明: 生成因果变量筛选的Venn图(三维度: DAG出度、模型重要性、ATE显著性)
# 输入文件: output/15b_causal_retraining/selection_criteria_summary.csv
# 输出文件: figures/15b_causal_retraining/venn_selection.png
# 作者: Nature级别科研项目
# 日期: 2025-12-15
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/CausalSDMs")

# 设定 CRAN 镜像
options(repos = c(CRAN = "https://mirrors.sustech.edu.cn/CRAN/"))

# 加载必要的包
packages <- c("tidyverse", "ggVennDiagram", "sysfonts", "showtext")
for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
}

# 字体设置
try(
    {
        sysfonts::font_add(
            family = "Arial",
            regular = "C:/Windows/Fonts/arial.ttf",
            bold = "C:/Windows/Fonts/arialbd.ttf",
            italic = "C:/Windows/Fonts/ariali.ttf",
            bolditalic = "C:/Windows/Fonts/arialbi.ttf"
        )
        showtext::showtext_opts(dpi = 2400)
        showtext::showtext_auto(enable = TRUE)
    },
    silent = TRUE
)

dir.create("figures/15b_causal_retraining", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("生成因果变量筛选Venn图\n")
cat("======================================\n\n")

# 读取筛选标准汇总数据
criteria <- read.csv("output/15b_causal_retraining/selection_criteria_summary.csv", stringsAsFactors = FALSE)

# 根据source列构建三个集合
dag_vars <- criteria %>%
    filter(source == "DAG") %>%
    pull(variable) %>%
    unique()
importance_vars <- criteria %>%
    filter(source == "Importance") %>%
    pull(variable) %>%
    unique()
ate_vars <- criteria %>%
    filter(source == "ATE") %>%
    pull(variable) %>%
    unique()

cat("DAG Top 15 变量数:", length(dag_vars), "\n")
cat("模型重要性 Top 15 变量数:", length(importance_vars), "\n")
cat("ATE 显著变量数:", length(ate_vars), "\n")

# 构建Venn图输入列表
venn_list <- list(
    "DAG Topology\n(Out-degree Top 15)" = dag_vars,
    "Model Importance\n(Top 15)" = importance_vars,
    "Causal Effect\n(ATE P<0.05)" = ate_vars
)

# 计算并集
all_drivers <- unique(c(dag_vars, importance_vars, ate_vars))
cat("核心驱动因子总数 (并集):", length(all_drivers), "\n")

# 绘制Venn图
p_venn <- ggVennDiagram(
    venn_list,
    label_alpha = 0,
    label = "count",
    edge_size = 0.8
) +
    scale_fill_gradient(low = "#F4F4F9", high = "#377EB8", name = "Count") +
    scale_color_manual(values = c("#E41A1C", "#4DAF4A", "#984EA3")) +
    labs(
        title = "Core Driver Selection: Multi-criteria Triangulation",
        subtitle = paste0(
            "Total: ", length(all_drivers), " variables from 47 (",
            round((1 - length(all_drivers) / 47) * 100, 1), "% reduction)"
        )
    ) +
    theme_void(base_family = "Arial") +
    theme(
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5, margin = margin(b = 10)),
        plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40", margin = margin(b = 20)),
        legend.position = "right",
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(30, 80, 30, 80) # 大幅增加左右边距 (top, right, bottom, left)
    ) +
    coord_cartesian(clip = "off") # 防止标签被裁切

# 保存图像
ggsave("figures/15b_causal_retraining/venn_selection.png",
    plot = p_venn,
    width = 12, height = 8, units = "in", dpi = 2400, bg = "white"
)
ggsave("figures/15b_causal_retraining/venn_selection.svg",
    plot = p_venn,
    width = 12, height = 8, units = "in", bg = "white"
)

cat("\n======================================\n")
cat("Venn图生成完成\n")
cat("======================================\n\n")
cat("✓ 图件: figures/15b_causal_retraining/venn_selection.png\n")
cat("✓ 图件: figures/15b_causal_retraining/venn_selection.svg\n\n")
