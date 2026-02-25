#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 14a_structure_learning.R
# 功能说明: 使用Hill-Climbing (HC)算法进行因果结构学习（DAG）
# 方法: 评分驱动HC (bnlearn)，Bootstrap子样本稳定性评估
# 输入文件: output/04_collinearity/collinearity_removed.csv
# 输出文件: output/14_causal/graph_hc.rds, graph_hc_avg.rds
#          output/14_causal/boot_hc.rds, edges_summary.csv
# 作者: CausalSDMs项目
# 日期: 2025-10-24
# 更新: 2026-02-06 (移除PC算法，仅保留HC；可视化由14b系列脚本负责)
# 备注: 14a系列仅负责因果算法计算，不含可视化
# ==============================================================================

# ==============================================================================
# 0. 环境初始化
# ==============================================================================
rm(list = ls())
gc()
setwd("E:/CausalSDMs")

options(repos = c(CRAN = "https://mirrors.sustech.edu.cn/CRAN/"))

# 加载必要的包
packages <- c("tidyverse", "bnlearn")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

dir.create("output/14_causal/a_structure", showWarnings = FALSE, recursive = TRUE)

# ==============================================================================
# 1. 全局参数配置
# ==============================================================================
cat("\n======================================\n")
cat("因果结构学习 (HC算法)\n")
cat("======================================\n\n")

# Bootstrap参数
BOOT_R <- 800 # Bootstrap重复次数 (推荐500-1000)
BOOT_M_RATIO <- 0.8 # 每次子采样比例

# 边强度阈值
THRESHOLD_SINGLE <- 0.55

# 随机种子
set.seed(20251024)

cat("参数配置:\n")
cat(sprintf("  - Bootstrap重复次数: %d\n", BOOT_R))
cat(sprintf("  - 子采样比例: %.0f%%\n", BOOT_M_RATIO * 100))
cat(sprintf("  - 边强度阈值: %.2f\n", THRESHOLD_SINGLE))
cat("\n")

# ==============================================================================
# 2. 数据读取与预处理
# ==============================================================================
cat("步骤 1/3: 读取数据与预处理...\n")

dat <- read.csv("output/04_collinearity/collinearity_removed.csv")
exclude_cols <- c("id", "species", "lon", "lat", "source", "presence", "presence.1")
env_vars <- setdiff(colnames(dat), exclude_cols)

X <- dat[, env_vars, drop = FALSE]

# 连续变量标准化
X_scaled <- as.data.frame(scale(X))
X_scaled[is.na(X_scaled)] <- 0

cat(sprintf("  ✓ 样本数: %d\n", nrow(X_scaled)))
cat(sprintf("  ✓ 变量数: %d\n", ncol(X_scaled)))
cat(sprintf("  ✓ Bootstrap子样本量: %d\n", floor(BOOT_M_RATIO * nrow(X_scaled))))
cat("\n")

# ==============================================================================
# 3. HC算法结构学习与稳定性评估
# ==============================================================================
cat("步骤 2/3: HC算法 (Hill-Climbing, BIC-G)...\n")

hc_fit <- tryCatch(
  {
    bnlearn::hc(X_scaled, score = "bic-g")
  },
  error = function(e) {
    message("  HC with bic-g failed, trying bge...")
    bnlearn::hc(X_scaled, score = "bge")
  }
)
saveRDS(hc_fit, file = "output/14_causal/a_structure/graph_hc.rds")

cat(sprintf("  ✓ HC学习完成，边数: %d\n", nrow(bnlearn::arcs(hc_fit))))

# Bootstrap稳定性评估
cat("步骤 3/3: 计算Bootstrap稳定性...\n")
boot_hc <- bnlearn::boot.strength(
  data = X_scaled,
  R = BOOT_R,
  algorithm = "hc",
  algorithm.args = list(score = "bic-g"),
  m = floor(BOOT_M_RATIO * nrow(X_scaled))
)
saveRDS(boot_hc, file = "output/14_causal/a_structure/boot_hc.rds")

# 保存边强度表
edges_summary <- boot_hc %>%
  dplyr::select(from, to, strength, direction) %>%
  dplyr::arrange(dplyr::desc(strength))
write.csv(edges_summary, "output/14_causal/a_structure/edges_hc.csv", row.names = FALSE)

# 平均网络
avg_hc <- bnlearn::averaged.network(boot_hc, threshold = THRESHOLD_SINGLE)
saveRDS(avg_hc, file = "output/14_causal/a_structure/graph_hc_avg.rds")

cat(sprintf(
  "  ✓ HC Bootstrap完成，强边(≥%.2f): %d\n",
  THRESHOLD_SINGLE, sum(boot_hc$strength >= THRESHOLD_SINGLE)
))

# ==============================================================================
# 4. 输出摘要
# ==============================================================================
cat("\n======================================\n")
cat("因果结构学习 (HC算法) 完成\n")
cat("======================================\n\n")

cat("算法结果摘要:\n")
cat(sprintf(
  "  - HC: %d 强边 (平均强度 %.3f)\n",
  sum(boot_hc$strength >= THRESHOLD_SINGLE),
  mean(boot_hc$strength[boot_hc$strength >= THRESHOLD_SINGLE])
))
cat("\n")

cat("输出文件:\n")
cat("  数据:\n")
cat("    - output/14_causal/a_structure/edges_hc.csv (HC算法边强度列表)\n")
cat("  模型:\n")
cat("    - output/14_causal/a_structure/graph_hc.rds (HC算法原始网络)\n")
cat("    - output/14_causal/a_structure/graph_hc_avg.rds (HC算法平均网络)\n")
cat("    - output/14_causal/a_structure/boot_hc.rds (HC算法Bootstrap结果)\n")
cat("\n")

cat("后续步骤:\n")
cat("  - 运行 14b_dag_hc.R 生成DAG可视化\n")
cat("  - 运行 14c_ate_hc.R 生成ATE可视化\n")
cat("\n")

cat("✓ 脚本执行完成!\n\n")
