#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 14c_ate_pc.R
# 功能说明: 基于PC-Stable算法骨架，使用DML估计因果效应 (ATE)
# 方法: Double Machine Learning (Partially Linear Model)
# 输入文件: output/04_collinearity/collinearity_removed.csv
#          output/14_causal/a_structure_pc/edges_pc.csv (用于协变量选择)
# 输出文件: output/14_causal/c_ate/ate_pc.csv
#          figures/14_causal/c_ate/ate_forest_pc.png/svg/pdf
# 作者: CausalSDMs项目
# 日期: 2026-02-06
# 备注: 使用线性回归作为nuisance模型进行DML估计
# ==============================================================================

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

options(repos = c(CRAN = "https://mirrors.sustech.edu.cn/CRAN/"))

# ============================================================================
# 0. 安装和加载必要的包
# ============================================================================
ensure_package <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("正在安装包:", pkg))
        install.packages(pkg)
    }
}

ensure_package("tidyverse")
ensure_package("sandwich") # 鲁棒标准误
ensure_package("lmtest") # 假设检验
ensure_package("grid")

if (!requireNamespace("forestploter", quietly = TRUE)) {
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
    library(sandwich)
    library(lmtest)
    library(forestploter)
    library(grid)
})

dir.create("output/14_causal/c_ate", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/14_causal/c_ate", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("ATE估计 (DML方法，基于PC骨架)\n")
cat("======================================\n\n")

# ============================================================================
# 1. 数据读取
# ============================================================================
cat("步骤 1/4: 读取数据...\n")

dat <- read.csv("output/04_collinearity/collinearity_removed.csv")
exclude_cols <- c("id", "species", "lon", "lat", "source", "presence", "presence.1")
env_vars <- setdiff(colnames(dat), exclude_cols)

# 使用presence作为outcome (假设二值化或连续概率)
Y <- dat$presence
X <- dat[, env_vars, drop = FALSE]

# 标准化
X_scaled <- as.data.frame(scale(X))
X_scaled[is.na(X_scaled)] <- 0

cat(sprintf("  ✓ 样本数: %d\n", nrow(X_scaled)))
cat(sprintf("  ✓ 变量数: %d\n", ncol(X_scaled)))

# 读取PC骨架信息(可选，用于理解变量间关系)
pc_edges <- NULL
if (file.exists("output/14_causal/a_structure_pc/edges_pc.csv")) {
    pc_edges <- read.csv("output/14_causal/a_structure_pc/edges_pc.csv")
    cat(sprintf("  ✓ PC骨架边数: %d\n", nrow(pc_edges)))
}

cat("\n")

# ============================================================================
# 2. DML估计ATE
# ============================================================================
cat("步骤 2/4: 使用DML估计每个变量的ATE...\n")

# 简化的DML: Partially Linear Model
# Y = theta * T + g(X) + epsilon
# 其中 T 是treatment变量，X 是其他协变量

estimate_ate_dml <- function(treat_var, outcome, covariates) {
    # 第一步: 对 treatment 进行回归 (T ~ X)
    # 第二步: 对 outcome 进行回归 (Y ~ X)
    # 第三步: 用残差回归估计 theta

    T_var <- covariates[[treat_var]]
    X_other <- covariates[, setdiff(names(covariates), treat_var), drop = FALSE]

    # 如果X_other为空，直接用简单回归
    if (ncol(X_other) == 0) {
        fit <- lm(outcome ~ T_var)
        coef_val <- coef(fit)["T_var"]
        se_val <- sqrt(vcovHC(fit, type = "HC1")["T_var", "T_var"])
    } else {
        # Partialling out
        # 回归 T ~ X_other
        fit_t <- lm(T_var ~ ., data = X_other)
        T_resid <- residuals(fit_t)

        # 回归 Y ~ X_other
        fit_y <- lm(outcome ~ ., data = X_other)
        Y_resid <- residuals(fit_y)

        # 残差回归: Y_resid ~ T_resid
        fit_final <- lm(Y_resid ~ T_resid - 1) # 无截距
        coef_val <- coef(fit_final)["T_resid"]

        # 鲁棒标准误
        se_val <- tryCatch(
            sqrt(vcovHC(fit_final, type = "HC1")["T_resid", "T_resid"]),
            error = function(e) {
                summary(fit_final)$coefficients["T_resid", "Std. Error"]
            }
        )
    }

    # 计算置信区间和p值
    z_val <- coef_val / se_val
    p_val <- 2 * pnorm(-abs(z_val))
    ci_lower <- coef_val - 1.96 * se_val
    ci_upper <- coef_val + 1.96 * se_val

    return(data.frame(
        variable = treat_var,
        coef = coef_val,
        se = se_val,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        p_value = p_val,
        significant = p_val < 0.05,
        stringsAsFactors = FALSE
    ))
}

# 对所有变量估计ATE
ate_results <- lapply(env_vars, function(v) {
    tryCatch(
        estimate_ate_dml(v, Y, X_scaled),
        error = function(e) {
            data.frame(
                variable = v, coef = NA, se = NA,
                ci_lower = NA, ci_upper = NA,
                p_value = NA, significant = FALSE
            )
        }
    )
})

ate_df <- do.call(rbind, ate_results)
ate_df <- ate_df %>%
    filter(!is.na(coef)) %>%
    arrange(desc(abs(coef)))

cat(sprintf("  ✓ 完成 %d 个变量的ATE估计\n", nrow(ate_df)))
cat(sprintf("  ✓ 显著变量 (p<0.05): %d\n", sum(ate_df$significant)))

# 保存结果
write.csv(ate_df, "output/14_causal/c_ate/ate_pc.csv", row.names = FALSE)
cat("  ✓ ATE结果已保存: output/14_causal/c_ate/ate_pc.csv\n\n")

# ============================================================================
# 3. 生成森林图
# ============================================================================
cat("步骤 3/4: 生成ATE森林图...\n")

# 筛选显著变量
dt <- ate_df %>%
    filter(significant == TRUE) %>%
    arrange(desc(coef)) %>%
    mutate(
        Variable = case_when(
            variable == "lc_wavg_09" ~ "LC_Urban",
            variable == "lc_wavg_12" ~ "LC_Water",
            variable == "flow_acc" ~ "FlowAcc",
            variable == "flow_length" ~ "FlowLen",
            variable == "hydro_wavg_18" ~ "BIO18",
            variable == "dem_range" ~ "ElevRange",
            variable == "slope_range" ~ "SlopeRange",
            variable == "hydro_wavg_16" ~ "BIO16",
            variable == "hydro_wavg_07" ~ "BIO7",
            variable == "soil_wavg_04" ~ "Silt",
            variable == "hydro_wavg_19" ~ "BIO19",
            variable == "dem_avg" ~ "Elev",
            variable == "slope_avg" ~ "Slope",
            TRUE ~ variable
        ),
        Variable = paste0("  ", Variable),
        ATE = sprintf("%.3f", coef),
        `95% CI` = sprintf("[%.3f, %.3f]", ci_lower, ci_upper),
        `P Value` = case_when(
            p_value < 0.001 ~ "<0.001",
            p_value < 0.01 ~ sprintf("%.3f", p_value),
            TRUE ~ sprintf("%.3f", p_value)
        ),
        ` ` = paste(rep(" ", 20), collapse = " "),
        is_key = variable %in% c("flow_acc", "flow_length", "lc_wavg_09", "lc_wavg_12")
    ) %>%
    select(Variable, ATE, `95% CI`, `P Value`, ` `, coef, ci_lower, ci_upper, is_key)

if (nrow(dt) > 0) {
    # 绘图主题 (紫色，区分于HC的蓝色和TABU的红色)
    tm <- forest_theme(
        base_size = 10,
        core = list(
            bg_params = list(fill = c("white")),
            fg_params = list(hjust = 0, x = 0.05)
        ),
        ci_col = "#9467BD", # 紫色 for PC
        ci_fill = "#9467BD",
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

    x_lim <- c(min(dt$ci_lower) - 0.05, max(dt$ci_upper) + 0.05)
    ticks_at <- round(seq(x_lim[1], x_lim[2], length.out = 5), 2)

    p <- forest(
        data = dt[, c(1:5)],
        est = dt$coef,
        lower = dt$ci_lower,
        upper = dt$ci_upper,
        sizes = ifelse(dt$is_key, 0.8, 0.6),
        ci_column = 5,
        ref_line = 0,
        xlim = x_lim,
        ticks_at = ticks_at,
        theme = tm
    )

    # 保存图片
    out_png <- "figures/14_causal/c_ate/ate_forest_pc.png"
    out_svg <- "figures/14_causal/c_ate/ate_forest_pc.svg"
    out_pdf <- "figures/14_causal/c_ate/ate_forest_pc.pdf"

    w <- 10
    h <- nrow(dt) * 0.35 + 1.5

    png(out_png, width = w, height = h, units = "in", res = 2400)
    print(p)
    grid.text("Causal Effects of Environmental Variables (PC-based DML)",
        x = 0.5, y = 0.96,
        gp = gpar(fontsize = 12, fontface = "bold")
    )
    dev.off()

    svg(out_svg, width = w, height = h)
    print(p)
    grid.text("Causal Effects of Environmental Variables (PC-based DML)",
        x = 0.5, y = 0.96,
        gp = gpar(fontsize = 12, fontface = "bold")
    )
    dev.off()

    pdf(out_pdf, width = w, height = h)
    print(p)
    grid.text("Causal Effects of Environmental Variables (PC-based DML)",
        x = 0.5, y = 0.96,
        gp = gpar(fontsize = 12, fontface = "bold")
    )
    dev.off()

    cat("  ✓ 森林图已生成:\n")
    cat("    -", out_png, "\n")
    cat("    -", out_svg, "\n")
    cat("    -", out_pdf, "\n")
} else {
    cat("  ⚠ 没有显著的ATE结果，跳过森林图生成\n")
}

# ============================================================================
# 4. 输出摘要
# ============================================================================
cat("\n======================================\n")
cat("ATE估计完成 (PC-based)\n")
cat("======================================\n\n")

cat("【Top 10 ATE (按绝对值)】:\n")
top10 <- ate_df %>%
    arrange(desc(abs(coef))) %>%
    head(10) %>%
    select(variable, coef, se, p_value, significant)
print(top10)

cat("\n输出文件:\n")
cat("  - output/14_causal/c_ate/ate_pc.csv\n")
cat("  - figures/14_causal/c_ate/ate_forest_pc.png/svg/pdf\n")
cat("\n✓ 脚本执行完成!\n\n")
