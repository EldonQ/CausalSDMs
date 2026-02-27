# =============================================================================
# Fig S5: Critical Difference Diagram (Friedman + Nemenyi post-hoc test)
# 对所有物种的AUC进行非参数多模型检验，展示模型间的显著性差异
# 数据来源: results/case2/all_results_v3.csv
# 依赖包: PMCMRplus (非参数多重比较), ggplot2
# 输出: figures/case2/plot/figS5_critical_difference_nemenyi.{png,svg}
# =============================================================================

suppressPackageStartupMessages({
    library(ggplot2)
    library(dplyr)
    library(tidyr)
    library(scales)
})

# 检查并安装 PMCMRplus
if (!requireNamespace("PMCMRplus", quietly = TRUE)) {
    install.packages("PMCMRplus", repos = "https://cloud.r-project.org")
}
library(PMCMRplus)

# ── 路径配置 ──────────────────────────────────────────────────────────────────
script_dir <- tryCatch(
    dirname(rstudioapi::getSourceEditorContext()$path),
    error = function(e) getwd()
)
proj_root <- normalizePath(file.path(script_dir, "..", "..", ".."))
data_dir  <- file.path(proj_root, "results", "case2")
fig_dir   <- file.path(proj_root, "figures", "case2", "plot")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ── 读取数据 ─────────────────────────────────────────────────────────────────
res_file <- file.path(data_dir, "all_results_v3.csv")
if (!file.exists(res_file)) stop("数据文件不存在: ", res_file)
results <- read.csv(res_file, stringsAsFactors = FALSE)

# 构造 species_id = region + "_" + species（唯一标识）
results <- results %>%
    mutate(sp_id = paste(region, species, sep = "_")) %>%
    filter(!is.na(auc_mean))

# 筛选关注的模型
focus_models <- c("CAST", "MLP_ATE", "MLP", "RF", "BRT", "GAM", "MaxEnt")
results_focus <- results %>%
    filter(model %in% focus_models)

# 仅保留所有模型都有记录的物种
complete_sps <- results_focus %>%
    group_by(sp_id) %>%
    summarise(n_models = n_distinct(model), .groups = "drop") %>%
    filter(n_models == length(focus_models)) %>%
    pull(sp_id)

mat <- results_focus %>%
    filter(sp_id %in% complete_sps) %>%
    select(sp_id, model, auc_mean) %>%
    pivot_wider(names_from = model, values_from = auc_mean,
                values_fn = mean) %>%
    tibble::column_to_rownames("sp_id") %>%
    as.matrix()

n_sp <- nrow(mat)
cat(sprintf("用于统计检验的物种数: %d\n", n_sp))

# ── Friedman检验 ─────────────────────────────────────────────────────────────
friedman_result <- friedman.test(mat)
cat(sprintf("Friedman test: chi^2 = %.3f, df = %d, p = %.4e\n",
            friedman_result$statistic, friedman_result$parameter,
            friedman_result$p.value))

# ── Nemenyi多重比较 ──────────────────────────────────────────────────────────
nemenyi_result <- frdAllPairsNemenyiTest(mat)
p_matrix <- nemenyi_result$p.value
cat("Nemenyi p-value matrix:\n")
print(round(p_matrix, 4))

# ── 计算平均排名 ──────────────────────────────────────────────────────────────
# 每行（物种）排名（rank 1 = 最好）
rank_mat <- t(apply(-mat, 1, rank, ties.method = "average"))
mean_ranks <- sort(colMeans(rank_mat))

rank_df <- data.frame(
    model      = names(mean_ranks),
    mean_rank  = as.numeric(mean_ranks),
    stringsAsFactors = FALSE
)
rank_df$model <- factor(rank_df$model, levels = rank_df$model)

# ── Critical Difference 计算 (Nemenyi, alpha = 0.05) ─────────────────────────
# CD = q_alpha * sqrt(k*(k+1) / (6*n))
# q_alpha for k=7 at alpha=0.05 ≈ 2.949 (from Nemenyi table)
k <- length(focus_models)
q_alpha <- qtukey(0.95, nmeans = k, df = Inf) / sqrt(2)  # 近似
CD <- q_alpha * sqrt(k * (k + 1) / (6 * n_sp))
cat(sprintf("Critical Difference (CD) = %.4f at alpha=0.05\n", CD))

# ── 构建显著性分组（不可区分的模型连线） ────────────────────────────────────
# 对 p_matrix 整理为对称矩阵
all_models <- rank_df$model
p_sym <- matrix(1, nrow = k, ncol = k,
                dimnames = list(all_models, all_models))
for (i in 1:(k - 1)) {
    for (j in (i + 1):k) {
        mi <- as.character(all_models[i])
        mj <- as.character(all_models[j])
        pval <- tryCatch(p_matrix[mj, mi], error = function(e)
            tryCatch(p_matrix[mi, mj], error = function(e2) 1))
        p_sym[mi, mj] <- p_sym[mj, mi] <- ifelse(is.na(pval), 1, pval)
    }
}

# 不显著对（p >= 0.05）用横线连接
non_sig_pairs <- which(p_sym >= 0.05, arr.ind = TRUE)
non_sig_pairs <- non_sig_pairs[non_sig_pairs[, 1] < non_sig_pairs[, 2], , drop = FALSE]
cd_segments <- data.frame(
    x    = mean_ranks[non_sig_pairs[, 1]],
    xend = mean_ranks[non_sig_pairs[, 2]],
    y    = seq(0.3, 0.3 + 0.15 * max(1, nrow(non_sig_pairs) - 1),
               length.out = max(1, nrow(non_sig_pairs))),
    stringsAsFactors = FALSE
)

# ── 绘制 CD Diagram ───────────────────────────────────────────────────────────
# 高亮CAST
rank_df$highlight <- rank_df$model == "CAST"
rank_df$color_group <- ifelse(rank_df$model == "CAST", "CAST",
                       ifelse(rank_df$model %in% c("MLP_ATE", "MLP"), "Ablation",
                              "Baseline"))
group_colors <- c("CAST" = "#E64B35", "Ablation" = "#4DBBD5", "Baseline" = "grey50")

p_cd <- ggplot(rank_df, aes(x = mean_rank, y = 0)) +
    # CD区间（以CAST为中心）
    annotate("rect",
             xmin = mean_ranks["CAST"] - CD, xmax = mean_ranks["CAST"] + CD,
             ymin = -0.18, ymax = 0.18,
             fill = "#E64B351A", color = "#E64B3566", linewidth = 0.5) +
    annotate("text",
             x = mean_ranks["CAST"] + CD + 0.02, y = 0.05,
             label = sprintf("CD = %.2f", CD),
             hjust = 0, size = 3, color = "#E64B35", family = "Arial") +
    # 不显著连接线
    {
        if (nrow(cd_segments) > 0) {
            geom_segment(data = cd_segments,
                         aes(x = x, xend = xend, y = y, yend = y),
                         color = "grey60", linewidth = 1.2, inherit.aes = FALSE)
        }
    } +
    # 模型点
    geom_point(aes(color = color_group), size = 4, shape = 16) +
    # 模型名称标注
    geom_text(
        aes(label = sprintf("%s\n(%.2f)", model, mean_rank),
            color = color_group, vjust = ifelse(seq_len(k) %% 2 == 0, 2.2, -1.2)),
        size = 3.2, family = "Arial"
    ) +
    # 排名轴
    geom_hline(yintercept = 0, linewidth = 0.5, color = "grey30") +
    scale_color_manual(values = group_colors) +
    scale_x_continuous(
        name = "Average Rank (lower = better)",
        breaks = 1:k,
        expand = expansion(mult = 0.12)
    ) +
    scale_y_continuous(limits = c(-0.5, 0.6)) +
    labs(
        title    = "Critical Difference Diagram (Nemenyi Test, α = 0.05)",
        subtitle = sprintf(
            "Friedman: χ² = %.2f, p = %.2e | N = %d species | Horizontal bars = non-significant differences",
            friedman_result$statistic, friedman_result$p.value, n_sp
        ),
        y = NULL
    ) +
    theme_classic(base_family = "Arial", base_size = 11) +
    theme(
        axis.text.y      = element_blank(),
        axis.ticks.y     = element_blank(),
        axis.line.y      = element_blank(),
        legend.position  = "none",
        plot.title       = element_text(face = "bold", size = 13),
        plot.subtitle    = element_text(size = 9, color = "grey40"),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.4)
    )

# ── 附：p值热力图 ─────────────────────────────────────────────────────────────
p_long <- as.data.frame(as.table(p_sym)) %>%
    rename(model1 = Var1, model2 = Var2, pval = Freq) %>%
    filter(as.character(model1) != as.character(model2)) %>%
    mutate(
        sig_level = case_when(
            pval < 0.001 ~ "p<0.001",
            pval < 0.01  ~ "p<0.01",
            pval < 0.05  ~ "p<0.05",
            TRUE         ~ "n.s."
        ),
        sig_level = factor(sig_level, levels = c("p<0.001","p<0.01","p<0.05","n.s."))
    )

model_order <- as.character(rank_df$model)
p_long$model1 <- factor(p_long$model1, levels = model_order)
p_long$model2 <- factor(p_long$model2, levels = rev(model_order))

p_heatmap <- ggplot(p_long, aes(x = model1, y = model2, fill = sig_level)) +
    geom_tile(color = "white", linewidth = 0.5) +
    geom_text(aes(label = ifelse(pval < 0.001, "***",
                          ifelse(pval < 0.01,  "**",
                          ifelse(pval < 0.05,  "*", "")))),
              size = 4, color = "white", family = "Arial") +
    scale_fill_manual(
        values = c("p<0.001" = "#E64B35", "p<0.01" = "#F39B7F",
                   "p<0.05"  = "#FBDECF", "n.s."   = "grey88"),
        name = "Significance"
    ) +
    labs(
        title = "Pairwise Nemenyi Post-hoc p-values",
        x = NULL, y = NULL
    ) +
    theme_classic(base_family = "Arial", base_size = 10) +
    theme(
        axis.text.x  = element_text(angle = 45, hjust = 1),
        plot.title   = element_text(face = "bold", size = 11),
        legend.title = element_text(size = 9)
    )

# ── 合并 ─────────────────────────────────────────────────────────────────────
suppressPackageStartupMessages(library(patchwork))
p_combined <- p_cd / p_heatmap + plot_layout(heights = c(1.2, 1)) +
    plot_annotation(
        title = "Statistical Comparison of Model Performance (AUC)",
        theme = theme(
            plot.title = element_text(face = "bold", size = 14, family = "Arial")
        )
    )

# ── 保存 ─────────────────────────────────────────────────────────────────────
out_prefix <- file.path(fig_dir, "figS5_critical_difference_nemenyi")
ggsave(paste0(out_prefix, ".png"), p_combined, width = 10, height = 11,
       dpi = 1200, bg = "white")
ggsave(paste0(out_prefix, ".svg"), p_combined, width = 10, height = 11, bg = "white")
cat("Fig S5 saved:", fig_dir, "\n")
