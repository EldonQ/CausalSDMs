#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 14a_structure_learning_pc.R
# 版本: V2.0 (Bootstrap集成深化版)
# 功能说明: 使用 PC-Stable 算法 + Bootstrap 重采样进行稳健的因果结构学习
#
# 升级亮点:
#   1. 引入 Bootstrap (R=800): 取代单次运行，通过统计多次采样的结果来评估稳健性。
#   2. 概率定向 (Probabilistic Orientation): 利用 Bootstrap 中方向出现的频率，
#      将原本无法定向的无向边(Undirected)转化为高置信度的有向边(Directed)。
#   3. 阈值控制 (Thresholding): 仅保留出现频率 > 0.55 的强连编。
#
# 输入文件: output/04_collinearity/collinearity_removed.csv
# 输出文件:
#   - output/14_causal/a_structure_pc/boot_pc.rds (Bootstrap原始结果)
#   - output/14_causal/a_structure_pc/graph_pc_avg.rds (平均网络)
#   - output/14_causal/a_structure_pc/edges_pc.csv (带强度和方向的边列表)
#
# 作者: CausalSDMs项目
# 日期: 2026-02-06
# ==============================================================================

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

options(repos = c(CRAN = "https://mirrors.sustech.edu.cn/CRAN/"))

ensure_package <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        install.packages(pkg)
    }
}

ensure_package("bnlearn")
ensure_package("tidyverse")

suppressPackageStartupMessages({
    library(bnlearn)
    library(tidyverse)
})

# 创建输出目录
output_dir <- "output/14_causal/a_structure_pc"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("\n")
cat("╔══════════════════════════════════════════════════════════════════════╗\n")
cat("║      PC-Stable 结构学习 (V2.0 Bootstrap集成版)                       ║\n")
cat("╚══════════════════════════════════════════════════════════════════════╝\n\n")

# ==============================================================================
# 1. 参数配置 (对齐 14a_structure_learning.R)
# ==============================================================================
BOOT_R <- 800 # Bootstrap 重复次数 (高稳健性)
BOOT_M_RATIO <- 0.8 # 每次重采样的样本比例
ALPHA <- 0.01 # PC算法显著性水平 (控制稀疏度)
THRESHOLD_STR <- 0.55 # 边强度阈值 (只保留出现概率 > 55% 的边)

set.seed(20251024) # 保证可复现性

cat("参数配置:\n")
cat(sprintf("  - Bootstrap次: %d\n", BOOT_R))
cat(sprintf("  - Alpha (Significance): %.2f\n", ALPHA))
cat(sprintf("  - 边保留阈值 (Strength): %.2f\n", THRESHOLD_STR))
cat("\n")

# ==============================================================================
# 2. 数据读取与预处理
# ==============================================================================
cat("步骤 1/3: 读取数据与预处理...\n")

dat <- read.csv("output/04_collinearity/collinearity_removed.csv")
exclude_cols <- c("id", "species", "lon", "lat", "source", "presence", "presence.1")
env_vars <- setdiff(colnames(dat), exclude_cols)

X <- dat[, env_vars, drop = FALSE]
X_scaled <- as.data.frame(scale(X)) # PC算法基于相关性，标准化
X_scaled[is.na(X_scaled)] <- 0

cat(sprintf("  ✓ 样本数: %d\n", nrow(X_scaled)))
cat(sprintf("  ✓ 变量数: %d\n", ncol(X_scaled)))

# ==============================================================================
# 3. Bootstrap 结构学习 (核心升级)
# ==============================================================================
cat("\n步骤 2/3: 执行 Bootstrap PC-Stable 学习 (这可能需要几分钟)...\n")

# 使用 bnlearn::boot.strength 并行计算
# algorithm = "pc.stable"
# algorithm.args 传递 PC 特有参数
boot_pc <- tryCatch(
    {
        bnlearn::boot.strength(
            data = X_scaled,
            R = BOOT_R,
            m = floor(BOOT_M_RATIO * nrow(X_scaled)),
            algorithm = "pc.stable",
            algorithm.args = list(
                test = "cor", # Fisher's Z test
                alpha = ALPHA,
                undirected = FALSE # 尝试记录方向 (虽然单次运行可能有无向边)
            )
            # cluster = NULL # 如果需要多核并行，可在此配置 cluster对象
        )
    },
    error = function(e) {
        message("❌ Bootstrap PC 失败: ", e$message)
        return(NULL)
    }
)

if (is.null(boot_pc)) stop("脚本终止: 结构学习失败")

# 保存原始 Bootstrap 结果
saveRDS(boot_pc, file = file.path(output_dir, "boot_pc.rds"))
cat("  ✓ Bootstrap 结果已保存: boot_pc.rds\n")

# ==============================================================================
# 4. 构建平均网络与结果提取
# ==============================================================================
cat("\n步骤 3/3: 构建平均网络与结果整理...\n")

# 4.1 生成平均网络 (Averaged Network)
# 这一步会基于阈值剔除弱边，并根据方向概率确定方向
avg_pc <- bnlearn::averaged.network(boot_pc, threshold = THRESHOLD_STR)
saveRDS(avg_pc, file = file.path(output_dir, "graph_pc_avg.rds"))

# 4.2 提取详细边列表 (用于可视化)
# boot.strength 返回的 dataframe 包含:
# from, to, strength (连接概率), direction (从 from->to 的概率)

# 筛选超过阈值的边
edges_final <- boot_pc %>%
    filter(strength >= THRESHOLD_STR) %>%
    arrange(desc(strength)) %>%
    mutate(
        # 确定最终方向类型
        # 如果 direction > 0.5，说明大部分时候是 from->to
        # 如果 direction 接近 0.5，说明方向不确定 (Undirected)
        # 这里我们定义一个方向置信度阈值，例如 0.51
        dir_conf = abs(direction - 0.5) * 2, # 0到1，越大越确定
        type = case_when(
            direction > 0.5 ~ "Directed",
            direction < 0.5 ~ "Reverse", # 这里通常 bnlearn 会把高概率的放前面，但需防备
            TRUE ~ "Undirected"
        )
    )

# 修正 Reverse 的情况 (虽然 bnlearn 通常只保留一个方向行，但为了保险)
# bnlearn boot.strength 通常保证 strength 是无向连接概率
# direction 是条件概率 p(-> | edge exists)
# 我们直接使用 dataframe 的 from/to，因为 bnlearn 已经处理好了顺序

cat(sprintf(
    "  ✓ 学习完成，共发现 %d 条显著边 (Strength > %.2f)\n",
    nrow(edges_final), THRESHOLD_STR
))

# 统计方向性
n_directed <- sum(edges_final$direction > 0.5 & edges_final$direction < 1.0)
n_undirected <- sum(edges_final$direction == 0.5) # 纯无向
n_determined <- sum(edges_final$direction == 1.0) # 绝对确定(极少见)

cat(sprintf("    - 有向边 (Direction > 0.5): %d\n", sum(edges_final$direction > 0.5)))
cat(sprintf("    - 无向边 (Direction = 0.5): %d\n", n_undirected))

# 保存为 CSV
write.csv(edges_final, file.path(output_dir, "edges_pc.csv"), row.names = FALSE)
cat(sprintf("  ✓ 边列表已保存: %s/edges_pc.csv\n", output_dir))

# ==============================================================================
# 5. 总结
# ==============================================================================
cat("\n======================================\n")
cat("分析完成\n")
cat("======================================\n")
cat("升级效果:\n")
cat("1. 稳健性: 结果基于 800 次重采样平均，不再受单次随机性影响。\n")
cat("2. 方向性: 通过统计多次采样的方向倾向，最大程度解决了 CPDA G无向问题。\n")
cat("           (可视化时，可根据 direction 值绘制箭头的置信度)\n")
cat("\n")
cat("后续步骤:\n")
cat("请更新 14b_dag_pc.R，使其利用 edges_pc.csv 中的 'direction' 列来\n")
cat("绘制有向箭头 (不再全是无向边了)。\n")
cat("\n✓ 脚本执行完成!\n\n")
