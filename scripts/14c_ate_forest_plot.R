#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 02_ate_forest_plot.R
# 功能说明: 绘制因果效应(ATE)森林图 - 使用 forestploter 包
# 风格: Nature/Lancet 顶刊风格 (左表右图)
# ==============================================================================

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

options(repos = c(CRAN = "https://cloud.r-project.org"))

# ============================================================================
# 0. 自动安装依赖包
# ============================================================================
ensure_package <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("正在安装包:", pkg))
        install.packages(pkg)
    }
}

ensure_package("tidyverse")
ensure_package("grid")

# forestploter 可能需要从 Github 或 CRAN 安装
if (!requireNamespace("forestploter", quietly = TRUE)) {
    message("正在安装 forestploter...")
    tryCatch(
        {
            install.packages("forestploter")
        },
        error = function(e) {
            message("CRAN安装失败，尝试从GitHub安装...")
            ensure_package("remotes")
            remotes::install_github("alandipert/forestploter")
        }
    )
}

suppressPackageStartupMessages({
    library(tidyverse)
    library(forestploter)
    library(grid)
})

# 创建输出目录
dir.create("scripts/drawing/figures", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/14_causal", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("生成专业级 ATE Forest Plot (forestploter)\n")
cat("======================================\n\n")

# ============================================================================
# 1. 数据准备
# ============================================================================
ate_data <- read.csv("output/14_causal/ate_all_variables.csv")

# 筛选显著变量并处理
dt <- ate_data %>%
    filter(significant == TRUE) %>%
    arrange(desc(coef)) %>% # 按ATE降序排列
    mutate(
        # 1. 变量名称映射 (与论文 Table 1 一致)
        Variable = case_when(
            variable == "lc_wavg_09" ~ "LC_Urban (Urban/built-up)",
            variable == "lc_wavg_12" ~ "LC_Water (Open water)",
            variable == "flow_acc" ~ "FlowAcc (Flow accumulation)",
            variable == "flow_length" ~ "FlowLen (Flow length)",
            variable == "hydro_wavg_18" ~ "BIO18 (Precip. warmest qtr)",
            variable == "dem_range" ~ "ElevRange (Elevation range)",
            variable == "slope_range" ~ "SlopeRange (Slope range)",
            variable == "hydro_wavg_16" ~ "BIO16 (Precip. wettest qtr)",
            variable == "hydro_wavg_07" ~ "BIO7 (Temp. annual range)",
            variable == "soil_wavg_04" ~ "Silt (Silt content)",
            variable == "hydro_wavg_19" ~ "BIO19 (Precip. coldest qtr)",
            variable == "dem_avg" ~ "Elev (Mean elevation)",
            variable == "slope_avg" ~ "Slope (Mean slope)",
            TRUE ~ variable
        ),

        # 2. 格式化数值列 (增加缩进以区分层级，如果有的话)
        # 这里为了美观，我们将 Variable 这一列作为主要的文本列
        Variable = ifelse(
            variable %in% c("flow_acc", "flow_length"),
            paste0("  ", Variable, "*"), # 关键变量加星号
            paste0("  ", Variable)
        ),
        ATE = sprintf("%.3f", coef),
        `95% CI` = sprintf("[%.3f, %.3f]", ci_lower, ci_upper),
        `P Value` = case_when(
            p_value < 0.001 ~ "< 0.001",
            p_value < 0.01 ~ sprintf("%.3f", p_value),
            TRUE ~ sprintf("%.3f", p_value)
        ),

        # 3. 定义绘图所需的空列 (用于放置森林图)
        ` ` = paste(rep(" ", 20), collapse = " "),

        # 4. 辅助列
        se = (ci_upper - coef) / 1.96,
        is_topology = variable %in% c("flow_acc", "flow_length"),
        direction = ifelse(coef > 0, "Positive", "Negative")
    ) %>%
    select(Variable, ATE, `95% CI`, `P Value`, ` `, coef, ci_lower, ci_upper, is_topology, direction)

# ============================================================================
# 2. 绘图参数设置
# ============================================================================

# 定义简洁主题 (使用统一颜色，避免 legend_value 错误)
tm <- forest_theme(
    base_size = 10,
    core = list(
        bg_params = list(fill = c("white")),
        fg_params = list(hjust = 0, x = 0.05) # 文本左对齐
    ),

    # 统一颜色设置 (使用深蓝色)
    ci_col = "#2C3E50",
    ci_fill = "#2C3E50",
    refline_col = "grey50",
    vertline_col = "grey90",

    # 点的形状和大小
    ci_pch = 16, # 统一使用圆形
    ci_lwd = 1.5,
    ci_Theight = 0.2,

    # 表格样式
    colhead = list(fg_params = list(hjust = 0.5, fontface = "bold")),
    summary_col = "black",
    footnote_col = "grey40",
    footnote_fontface = "italic"
)

# ============================================================================
# 3. 生成森林图
# ============================================================================

# 计算X轴范围
x_lim <- c(min(dt$ci_lower) - 0.05, max(dt$ci_upper) + 0.05)
ticks_at <- round(seq(x_lim[1], x_lim[2], length.out = 5), 2)

p <- forest(
    data = dt[, c(1:5)], # 显示前5列文本 (Variable, ATE, CI, P, Space)
    est = dt$coef,
    lower = dt$ci_lower,
    upper = dt$ci_upper,
    sizes = ifelse(dt$is_topology, 0.8, 0.6), # 关键变量点大一点
    ci_column = 5, # 森林图画在第5列(空格列)
    ref_line = 0,
    xlim = x_lim,
    ticks_at = ticks_at,
    theme = tm
)

# ============================================================================
# 4. 保存图片
# ============================================================================
# 设置输出路径
out_png <- "figures/14_causal/ate_all_variables_forest.png"
out_svg <- "figures/14_causal/ate_all_variables_forest.svg"
out_pdf <- "figures/14_causal/ate_all_variables_forest.pdf"

# 计算图片尺寸 (高度自适应)
w <- 10
h <- nrow(dt) * 0.35 + 1.5

# 保存 PNG
png(out_png, width = w, height = h, units = "in", res = 2400)
print(p)
grid.text("Causal Effects of Environmental Variables (DML)",
    x = 0.5, y = 0.96,
    gp = gpar(fontsize = 12, fontface = "bold")
)
dev.off()

# 保存 SVG
svg(out_svg, width = w, height = h)
print(p)
grid.text("Causal Effects of Environmental Variables (DML)",
    x = 0.5, y = 0.96,
    gp = gpar(fontsize = 12, fontface = "bold")
)
dev.off()

# 保存 PDF
pdf(out_pdf, width = w, height = h)
print(p)
grid.text("Causal Effects of Environmental Variables (DML)",
    x = 0.5, y = 0.96,
    gp = gpar(fontsize = 12, fontface = "bold")
)
dev.off()

cat("\n✓ 图件已生成:\n")
cat("  -", out_png, "\n")
cat("  -", out_svg, "\n")
cat("  -", out_pdf, "\n")
