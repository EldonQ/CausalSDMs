#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 14c_ate_hc.R
# 功能说明: 基于HC算法的DAG生成ATE森林图 (复用14c_ate_forest_plot.R逻辑)
# 输入文件: output/14_causal/c_ate/ate_hc.csv (由ATE估计脚本生成)
# 输出文件: figures/14_causal/c_ate/ate_forest_hc.png/svg/pdf
# 作者: CausalSDMs项目
# 日期: 2026-02-06
# 备注: 需要先运行ATE估计脚本生成ate_hc.csv
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

if (!requireNamespace("forestploter", quietly = TRUE)) {
    message("正在安装 forestploter...")
    tryCatch(
        install.packages("forestploter"),
        error = function(e) {
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

dir.create("figures/14_causal/c_ate", showWarnings = FALSE, recursive = TRUE)
dir.create("output/14_causal/c_ate", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("生成专业级 ATE Forest Plot (HC算法)\n")
cat("======================================\n\n")

# ============================================================================
# 1. 数据准备
# ============================================================================
# 优先读取HC专用的ATE结果，如果不存在则回退到通用ATE结果
ate_file <- "output/14_causal/c_ate/ate_hc.csv"
if (!file.exists(ate_file)) {
    ate_file <- "output/14_causal/ate_all_variables.csv"
    if (!file.exists(ate_file)) {
        stop("错误: 未找到ATE数据文件，请先运行ATE估计脚本")
    }
    cat("  使用通用ATE数据:", ate_file, "\n")
} else {
    cat("  使用HC专用ATE数据:", ate_file, "\n")
}

ate_data <- read.csv(ate_file)

# 筛选显著变量并处理
dt <- ate_data %>%
    filter(significant == TRUE) %>%
    arrange(desc(coef)) %>%
    mutate(
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
        Variable = ifelse(
            variable %in% c("flow_acc", "flow_length"),
            paste0("  ", Variable, "*"),
            paste0("  ", Variable)
        ),
        ATE = sprintf("%.3f", coef),
        `95% CI` = sprintf("[%.3f, %.3f]", ci_lower, ci_upper),
        `P Value` = case_when(
            p_value < 0.001 ~ "<0.001",
            p_value < 0.01 ~ sprintf("%.3f", p_value),
            TRUE ~ sprintf("%.3f", p_value)
        ),
        ` ` = paste(rep(" ", 20), collapse = " "),
        se = (ci_upper - coef) / 1.96,
        is_topology = variable %in% c("flow_acc", "flow_length"),
        direction = ifelse(coef > 0, "Positive", "Negative")
    ) %>%
    select(Variable, ATE, `95% CI`, `P Value`, ` `, coef, ci_lower, ci_upper, is_topology, direction)

if (nrow(dt) == 0) {
    stop("错误: 没有显著的ATE结果可用于绘图")
}

cat("  显著变量数:", nrow(dt), "\n")

# ============================================================================
# 2. 绘图参数设置
# ============================================================================
tm <- forest_theme(
    base_size = 10,
    core = list(
        bg_params = list(fill = c("white")),
        fg_params = list(hjust = 0, x = 0.05)
    ),
    ci_col = "#2C3E50",
    ci_fill = "#2C3E50",
    refline_col = "grey50",
    vertline_col = "grey90",
    ci_pch = 16,
    ci_lwd = 1.5,
    ci_Theight = 0.2,
    colhead = list(fg_params = list(hjust = 0.5, fontface = "bold")),
    summary_col = "black",
    footnote_col = "grey40",
    footnote_fontface = "italic"
)

# ============================================================================
# 3. 生成森林图
# ============================================================================
x_lim <- c(min(dt$ci_lower) - 0.05, max(dt$ci_upper) + 0.05)
ticks_at <- round(seq(x_lim[1], x_lim[2], length.out = 5), 2)

p <- forest(
    data = dt[, c(1:5)],
    est = dt$coef,
    lower = dt$ci_lower,
    upper = dt$ci_upper,
    sizes = ifelse(dt$is_topology, 0.8, 0.6),
    ci_column = 5,
    ref_line = 0,
    xlim = x_lim,
    ticks_at = ticks_at,
    theme = tm
)

# ============================================================================
# 4. 保存图片
# ============================================================================
out_png <- "figures/14_causal/c_ate/ate_forest_hc.png"
out_svg <- "figures/14_causal/c_ate/ate_forest_hc.svg"
out_pdf <- "figures/14_causal/c_ate/ate_forest_hc.pdf"

w <- 10
h <- nrow(dt) * 0.35 + 1.5

png(out_png, width = w, height = h, units = "in", res = 2400)
print(p)
grid.text("Causal Effects of Environmental Variables (HC-based DML)",
    x = 0.5, y = 0.96,
    gp = gpar(fontsize = 12, fontface = "bold")
)
dev.off()

svg(out_svg, width = w, height = h)
print(p)
grid.text("Causal Effects of Environmental Variables (HC-based DML)",
    x = 0.5, y = 0.96,
    gp = gpar(fontsize = 12, fontface = "bold")
)
dev.off()

pdf(out_pdf, width = w, height = h)
print(p)
grid.text("Causal Effects of Environmental Variables (HC-based DML)",
    x = 0.5, y = 0.96,
    gp = gpar(fontsize = 12, fontface = "bold")
)
dev.off()

cat("\n✓ 图件已生成:\n")
cat("  -", out_png, "\n")
cat("  -", out_svg, "\n")
cat("  -", out_pdf, "\n")
cat("\n完成!\n")
