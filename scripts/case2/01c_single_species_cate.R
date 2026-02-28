################################################################################
# 01c_single_species_cate.R
# 单物种 CATE 估计 + 空间地图（case2：disdat 基准数据）
#
# 目标：
# - 选择一个区域 + 物种 + 若干环境变量作为 Treatment
# - 使用 Causal Forest 估计条件平均处理效应 (CATE)
# - 输出：
#   1) output/case2/cate/cate_result_<region>_<species>_<var>.csv
#   2) figures/case2/plot/cate_map_<region>_<species>_<var>.png/svg  (1200 dpi)
# - 地图：
#   - 使用原始 disdat 的 x,y 坐标（投影坐标），自动计算数据外框
#   - 尝试叠加世界底图 / 区域边界（若 rnaturalearth 可用）
#
# 注意：
# - 本脚本与 01_cast_pipeline.R 解耦，直接从 disdat 重新读取数据，
#   确保拥有完整的坐标信息 (x,y)。
# - 绘图中的文字全部为英文；代码注释全部为中文。
################################################################################

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

options(repos = c(CRAN = "https://cloud.r-project.org"))

## 安装 / 加载所需 R 包（若缺失则自动安装）
ensure_package <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        message(paste("正在安装包:", pkg))
        install.packages(pkg)
    }
}

ensure_package("disdat")
ensure_package("grf")
ensure_package("tidyverse")
ensure_package("ggplot2")
ensure_package("viridis")
ensure_package("patchwork")
ensure_package("sf")
ensure_package("rnaturalearth")
ensure_package("rnaturalearthdata")

suppressPackageStartupMessages({
    library(disdat)
    library(grf)
    library(tidyverse)
    library(ggplot2)
    library(viridis)
    library(patchwork)
    library(sf)
    library(rnaturalearth)
    library(rnaturalearthdata)
})

## 输出路径配置 ---------------------------------------------------------------
dir.create("output/case2/cate", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/case2/plot", showWarnings = FALSE, recursive = TRUE)

## 用户可配置参数 -------------------------------------------------------------
# 区域与物种（与 01_cast_pipeline 一致；示例使用 SWI 区域中的 swi10）
TARGET_REGION  <- "SWI"     # 例如 "SWI"
TARGET_SPECIES <- "swi10"   # 例如 "swi10"

# Treatment 变量列表（须在该区域的 env_cols 中存在；不同区域变量名不同）
# 若留空或与当前区域不匹配，将自动使用该区域前 2 个环境变量
TREATMENT_VARS <- c(
    "flow_length",
    "lc_wavg_09"
)
# 例如 SWI 区域常见变量: "ann_mean_temp", "ann_precip", "min_temp_coldest_month" 等

# Causal Forest 参数（与 14d_causal_forest.R 类似，但为单物种精简版）
CF_PARAMS <- list(
    num.trees      = 2000L,
    min.node.size  = 5L,
    sample.fraction = 0.5,
    honesty        = TRUE,
    tune.parameters = "all",
    seed           = 42L
)

## 1. 读取 disdat 原始数据（含坐标）-------------------------------------------
cat("======================================================================\n")
cat("  01c_single_species_cate.R — 单物种 CATE + 空间地图\n")
cat("======================================================================\n\n")

cat(sprintf("  Region:  %s\n", TARGET_REGION))
cat(sprintf("  Species: %s\n\n", TARGET_SPECIES))

cat("步骤 1/4: 读取 disdat 原始数据（包含坐标）...\n")

po  <- disdat::disPo(TARGET_REGION)
bg  <- disdat::disBg(TARGET_REGION)
pa  <- disdat::disPa(TARGET_REGION)
env <- disdat::disEnv(TARGET_REGION)

if (!"siteid" %in% names(pa)) {
    stop("disdat::disPa 返回的数据中缺少 siteid 列。")
}

## 环境变量列：与 00_data_preparation 一致，用 po/bg 定义，再取与 env 的交集
meta_cols <- c("siteid", "spid", "x", "y", "occ", "group")
env_cols  <- setdiff(names(po), meta_cols)
env_cols  <- intersect(env_cols, names(bg))
env_cols  <- intersect(env_cols, names(env))

# 确定实际用于 CATE 的 Treatment 变量；若用户指定的不在本区域则用前 2 个 env 变量
treat_vars_to_use <- intersect(TREATMENT_VARS, env_cols)
if (length(treat_vars_to_use) == 0) {
    treat_vars_to_use <- head(env_cols, 2)
    cat(sprintf("  环境变量数: %d\n", length(env_cols)))
    cat(sprintf("  本区域 env 变量: %s\n", paste(env_cols, collapse = ", ")))
    cat(sprintf("  未匹配到 TREATMENT_VARS，改用前 2 个: %s\n\n",
        paste(treat_vars_to_use, collapse = ", ")))
} else {
    cat(sprintf("  环境变量数: %d\n", length(env_cols)))
    cat(sprintf("  可用 Treatment 变量: %s\n\n", paste(treat_vars_to_use, collapse = ", ")))
}

## 2. 构建单物种样本表（含坐标 + 响应 + 环境）---------------------------------
cat("步骤 2/4: 构建单物种数据集...\n")

# 物种列名通常是类似 "swi10" 的列，包含 0/1 或 NA
if (!TARGET_SPECIES %in% names(pa)) {
    stop(sprintf("在 disPa(%s) 中未找到物种列: %s", TARGET_REGION, TARGET_SPECIES))
}

# 与 00 一致：用 disEnv 与 pa 按共有列合并（pa 中无 spid 时仍可合并）
merge_keys <- intersect(names(pa), names(env))
pa_full    <- merge(pa, env, by = merge_keys, all.x = TRUE)
env_in_pa  <- intersect(env_cols, names(pa_full))
if (length(env_in_pa) < length(env_cols)) env_cols <- env_in_pa

# 提取该物种的 presence 及坐标、环境变量（坐标优先从 pa，若无则从 env）
need_cols <- c("siteid", "x", "y", "group", TARGET_SPECIES, env_cols)
need_cols <- intersect(need_cols, names(pa_full))
dat_pa    <- pa_full %>%
    select(all_of(need_cols)) %>%
    rename(presence = !!TARGET_SPECIES) %>%
    drop_na(presence)
if (!"x" %in% names(dat_pa) || !"y" %in% names(dat_pa)) {
    stop("合并后的数据中缺少 x 或 y 列，无法绘制空间图。请检查 disPa/disEnv 返回的列名。")
}

# 简单过滤：至少要有一定数量的 presence 样本
if (sum(dat_pa$presence == 1, na.rm = TRUE) < 30) {
    stop("该物种 presence 样本数 < 30，不适合做 CATE 空间图。")
}

cat(sprintf("  总样本数: %d (presence=%d, background=%d)\n",
    nrow(dat_pa),
    sum(dat_pa$presence == 1), sum(dat_pa$presence == 0)))

## 3. 对每个 Treatment 变量运行 Causal Forest，得到 CATE ----------------------
cat("\n步骤 3/4: 估计 CATE（Causal Forest）...\n")

results_list <- list()

for (treat_var in treat_vars_to_use) {
    cat("\n------------------------------------------------------------------\n")
    cat(sprintf("  Treatment: %s\n", treat_var))

    # Outcome
    Y <- dat_pa$presence

    # Treatment 原始与标准化
    W_raw <- dat_pa[[treat_var]]
    if (all(is.na(W_raw))) {
        cat("  ⚠ 该变量全为 NA，跳过。\n")
        next
    }
    W <- as.vector(scale(W_raw))

    # 协变量：其余环境变量
    covariate_vars <- setdiff(env_cols, treat_var)
    X_all <- dat_pa[, covariate_vars, drop = FALSE]
    X_scaled <- as.matrix(scale(X_all))
    X_scaled[is.na(X_scaled)] <- 0

    set.seed(CF_PARAMS$seed)
    cf <- tryCatch(
        {
            do.call(grf::causal_forest, c(
                list(X = X_scaled, Y = Y, W = W),
                CF_PARAMS[names(CF_PARAMS) %in% formalArgs(grf::causal_forest)]
            ))
        },
        error = function(e) {
            message("  ❌ Causal Forest 训练失败: ", e$message)
            return(NULL)
        }
    )

    if (is.null(cf)) next

    cat("  ✓ 模型训练完成，开始预测 CATE...\n")

    pred <- predict(cf, estimate.variance = TRUE)
    cate    <- as.numeric(pred$predictions)
    cate_se <- sqrt(as.numeric(pred$variance.estimates))

    # 估计整体平均效应
    ate_res <- tryCatch(
        {
            grf::average_treatment_effect(cf, target.sample = "all")
        },
        error = function(e) {
            c(estimate = mean(cate, na.rm = TRUE),
              std.err  = sd(cate, na.rm = TRUE) / sqrt(sum(!is.na(cate))))
        }
    )

    cat(sprintf("  平均效应 (ATE): %.4f (SE=%.4f)\n",
        ate_res["estimate"], ate_res["std.err"]))

    # 整理结果表（含坐标）
    res_df <- dat_pa %>%
        mutate(
            treatment_var  = treat_var,
            W              = W_raw,
            W_std          = W,
            cate           = cate,
            cate_se        = cate_se,
            cate_ci_lower  = cate - 1.96 * cate_se,
            cate_ci_upper  = cate + 1.96 * cate_se
        )

    # 输出 CSV 至 output/case2/cate/
    out_csv <- sprintf(
        "output/case2/cate/cate_result_%s_%s_%s.csv",
        TARGET_REGION, TARGET_SPECIES, treat_var
    )
    write.csv(res_df, out_csv, row.names = FALSE)
    cat("  ✓ 已保存 CATE 结果表: ", out_csv, "\n")

    # 存入列表以便后续统一出图
    results_list[[treat_var]] <- list(
        data   = res_df,
        ate    = ate_res,
        forest = cf
    )
}

if (length(results_list) == 0) {
    stop("没有成功估计的 Treatment 变量，脚本结束。")
}

## 4. 空间地图绘制（自动检查地理范围与底图）-----------------------------------
cat("\n步骤 4/4: 绘制 CATE 空间地图...\n")

# 原始坐标是 disdat 的 x,y 投影坐标（单位 ~ km），无法直接当经纬度使用。
# 为了快速出图，这里采用：
#   - 使用 x,y 直接作为平面坐标（保持 disdat 原生坐标系）
#   - 自动计算数据的 bounding box，用于设置坐标轴范围
#   - 可选：叠加世界轮廓的粗略底图（只做参考，不用于精确制图）

for (treat_var in names(results_list)) {
    res_df <- results_list[[treat_var]]$data

    # 计算 CATE 颜色范围（对称）
    vmax <- quantile(abs(res_df$cate), probs = 0.95, na.rm = TRUE)
    if (!is.finite(vmax) || vmax <= 0) vmax <- max(abs(res_df$cate), na.rm = TRUE)
    if (!is.finite(vmax) || vmax <= 0) vmax <- 1

    # x,y 范围（加一点 padding）
    x_range <- range(res_df$x, na.rm = TRUE)
    y_range <- range(res_df$y, na.rm = TRUE)
    pad_x <- diff(x_range) * 0.02
    pad_y <- diff(y_range) * 0.02
    x_lims <- c(x_range[1] - pad_x, x_range[2] + pad_x)
    y_lims <- c(y_range[1] - pad_y, y_range[2] + pad_y)

    # 地理区域边界：用样本点的凸包作为研究区轮廓（与 disdat 投影坐标一致）
    hull_idx <- chull(res_df$x, res_df$y)
    hull_df <- res_df[c(hull_idx, hull_idx[1]), c("x", "y")]

    # 热图风格、高对比配色：蓝 — 黄 — 红（参考附图）
    p_base <- ggplot() +
        # 研究区边界（黑色轮廓，无填充）
        geom_polygon(
            data = hull_df, aes(x = x, y = y),
            fill = NA, color = "black", linewidth = 0.5
        ) +
        geom_point(
            data = res_df,
            aes(x = x, y = y, color = cate),
            size = 0.8, alpha = 0.92, shape = 16
        ) +
        scale_color_gradient2(
            low = "#2166AC", mid = "#FFFFBF", high = "#B2182B",
            midpoint = 0,
            limits = c(-vmax, vmax),
            oob = scales::squish,
            name = "CATE"
        ) +
        coord_equal(xlim = x_lims, ylim = y_lims, expand = FALSE) +
        labs(
            title = sprintf("CATE map — %s (%s, %s)", treat_var, TARGET_REGION, TARGET_SPECIES),
            subtitle = "Red = positive causal effect; Blue = negative causal effect"
        ) +
        theme_void(base_family = "Arial", base_size = 11) +
        theme(
            plot.title = element_text(face = "bold", hjust = 0.5),
            plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 9),
            legend.position = "right",
            # 不显示横轴、纵轴与网格线
            axis.line = element_blank(),
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank()
        )

    out_png <- sprintf(
        "figures/case2/plot/cate_map_%s_%s_%s.png",
        TARGET_REGION, TARGET_SPECIES, treat_var
    )
    out_svg <- sprintf(
        "figures/case2/plot/cate_map_%s_%s_%s.svg",
        TARGET_REGION, TARGET_SPECIES, treat_var
    )

    ggsave(out_png, p_base, width = 6, height = 4.5, dpi = 1200, bg = "white")
    ggsave(out_svg, p_base, width = 6, height = 4.5, dpi = 1200, bg = "white")

    cat("  ✓ 已保存 CATE 空间图: ", out_png, "\n")
}

cat("\n全部 CATE 估计与空间地图生成完成。\n")

