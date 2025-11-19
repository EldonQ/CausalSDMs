#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 14c_batch_ate_estimation.R
# 功能说明: 批量计算多个候选变量的平均处理效应(ATE)，为因果驱动建模提供依据
# 方法: Double Machine Learning (DoubleML) - IRM模型
# 输入文件: output/04_collinearity/collinearity_removed.csv
#          output/09_variable_importance/importance_summary.csv (候选变量)
# 输出文件: output/14_causal/ate_all_variables.csv
#          figures/14_causal/ate_all_variables_forest.png
# 作者: Nature级别科研项目
# 日期: 2025-11-10
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/SDM01")

# 设定 CRAN 镜像
options(repos = c(CRAN = "https://mirrors.sustech.edu.cn/CRAN/"))

# 加载必要的包
packages <- c("tidyverse", "DoubleML", "mlr3", "mlr3learners", "ranger", 
              "ggplot2", "sysfonts", "showtext")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

dir.create("output/14_causal", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/14_causal", showWarnings = FALSE, recursive = TRUE)

# 字体设置
try({
  sysfonts::font_add(
    family = "Arial",
    regular = "C:/Windows/Fonts/arial.ttf",
    bold = "C:/Windows/Fonts/arialbd.ttf",
    italic = "C:/Windows/Fonts/ariali.ttf",
    bolditalic = "C:/Windows/Fonts/arialbi.ttf"
  )
  showtext::showtext_opts(dpi = 1200)
  showtext::showtext_auto(enable = TRUE)
}, silent = TRUE)

cat("\n======================================\n")
cat("批量ATE估计 (Top候选变量)\n")
cat("======================================\n\n")

# ============================================================================
# 步骤1: 选择候选变量（Top 20）
# ============================================================================
cat("步骤 1/3: 选择候选处理变量...\n")

# 读取数据
dat <- read.csv("output/04_collinearity/collinearity_removed.csv", stringsAsFactors = FALSE)
exclude_cols <- c("id", "species", "lon", "lat", "source", "presence", "presence.1")
env_vars <- setdiff(colnames(dat), exclude_cols)

# 读取变量重要性，选择Top 20
imp_path <- "output/09_variable_importance/importance_summary.csv"
if(file.exists(imp_path)) {
  imp_df <- read.csv(imp_path, stringsAsFactors = FALSE)
  candidate_vars <- imp_df %>% 
    group_by(variable) %>%
    summarise(mean_imp = mean(importance_normalized, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(mean_imp)) %>%
    head(20) %>%
    pull(variable)
  candidate_vars <- intersect(candidate_vars, env_vars)
} else {
  # 如果没有重要性文件，随机选择前20个
  candidate_vars <- head(env_vars, 20)
}

cat("  ✓ 候选变量数: ", length(candidate_vars), "\n")
cat("  变量列表:\n")
for(i in seq_along(candidate_vars)) {
  cat(sprintf("    %2d. %s\n", i, candidate_vars[i]))
}

# ============================================================================
# 步骤2: 批量估计ATE
# ============================================================================
cat("\n步骤 2/3: 批量估计ATE...\n")
cat("  (这可能需要几分钟...)\n\n")

y <- dat$presence
ate_results <- list()

for(i in seq_along(candidate_vars)) {
  treat_var <- candidate_vars[i]
  cat(sprintf("  [%2d/%2d] %s ... ", i, length(candidate_vars), treat_var))
  
  tryCatch({
    # 准备处理变量（二值化：中位数阈值）
    vx <- as.numeric(dat[[treat_var]])
    vx[!is.finite(vx)] <- NA
    thr <- median(vx, na.rm = TRUE)
    T_var <- as.numeric(vx > thr)
    
    # 协变量（所有其他环境变量）
    X <- dat[, setdiff(env_vars, treat_var), drop = FALSE]
    X[is.na(X)] <- 0
    
    # 构建DoubleML数据
    df <- data.frame(y = y, d = T_var, X)
    task <- DoubleML::DoubleMLData$new(df, y_col = "y", d_cols = "d")
    
    # 机器学习模型（回归树 + 分类树）
    ml_g <- mlr3::lrn("regr.ranger", num.trees = 300)
    ml_m <- mlr3::lrn("classif.ranger", num.trees = 300, predict_type = "prob")
    
    # 拟合IRM
    dml_irm <- DoubleML::DoubleMLIRM$new(task, ml_g = ml_g, ml_m = ml_m, n_folds = 3)
    dml_irm$fit()
    
    # 提取结果
    summ <- dml_irm$summary()
    summ_df <- as.data.frame(summ)
    
    # 标准化列名（不同版本的DoubleML可能列名不同）
    names(summ_df) <- tolower(gsub("[^a-z0-9]+", "_", names(summ_df)))
    
    # 识别系数和标准误列
    coef_col <- which(names(summ_df) %in% c("coef", "estimate", "theta"))[1]
    se_col <- which(names(summ_df) %in% c("std_error", "std_err", "se", "stderr"))[1]
    pval_col <- which(names(summ_df) %in% c("p_value", "pval", "p_val"))[1]
    
    if(is.na(coef_col)) coef_col <- 1
    if(is.na(se_col)) se_col <- 2
    if(is.na(pval_col)) pval_col <- ncol(summ_df)
    
    ate_results[[treat_var]] <- data.frame(
      variable = treat_var,
      coef = as.numeric(summ_df[1, coef_col]),
      std_error = as.numeric(summ_df[1, se_col]),
      p_value = as.numeric(summ_df[1, pval_col]),
      threshold = thr,
      n_treated = sum(T_var == 1, na.rm = TRUE),
      n_control = sum(T_var == 0, na.rm = TRUE)
    )
    
    cat("✓\n")
    
  }, error = function(e) {
    cat("✗ (", e$message, ")\n", sep = "")
    ate_results[[treat_var]] <<- data.frame(
      variable = treat_var,
      coef = NA,
      std_error = NA,
      p_value = NA,
      threshold = NA,
      n_treated = NA,
      n_control = NA
    )
  })
  
  # 清理内存
  gc(verbose = FALSE)
}

# 合并结果
ate_all <- bind_rows(ate_results)

# 计算置信区间
ate_all <- ate_all %>%
  mutate(
    ci_lower = coef - 1.96 * std_error,
    ci_upper = coef + 1.96 * std_error,
    significant = p_value < 0.05
  ) %>%
  arrange(desc(abs(coef)))

# 保存
write.csv(ate_all, "output/14_causal/ate_all_variables.csv", row.names = FALSE)

cat("\n  ✓ ATE估计完成\n")
cat("    显著变量数 (p<0.05): ", sum(ate_all$significant, na.rm = TRUE), "\n")

# ============================================================================
# 步骤3: 可视化
# ============================================================================
cat("\n步骤 3/3: 生成森林图...\n")

# 筛选有效结果（排除估计失败的）
ate_valid <- ate_all %>% filter(!is.na(coef))

if(nrow(ate_valid) > 0) {
  # 排序并限制显示数量（Top 15）
  ate_plot <- ate_valid %>%
    arrange(desc(abs(coef))) %>%
    head(15) %>%
    mutate(variable = factor(variable, levels = rev(variable)))
  
  p_forest <- ggplot(ate_plot, aes(x = coef, y = variable, color = significant)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    geom_errorbarh(aes(xmin = ci_lower, xmax = ci_upper), height = 0.3, linewidth = 0.5) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c("TRUE" = "#E41A1C", "FALSE" = "grey50"),
                       labels = c("TRUE" = "p < 0.05", "FALSE" = "p ≥ 0.05"),
                       name = "Significance") +
    labs(title = "Average Treatment Effects (ATE) - Top 15 Variables",
         subtitle = "Error bars: 95% confidence intervals",
         x = "ATE Estimate", 
         y = "Environmental Variable") +
    theme_minimal(base_family = "Arial", base_size = 8) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
      plot.title = element_text(face = "bold", size = 9),
      plot.subtitle = element_text(size = 7, color = "grey40"),
      legend.position = "top"
    )
  
  ggsave("figures/14_causal/ate_all_variables_forest.png", 
         plot = p_forest, width = 6, height = 5, dpi = 1200, bg = "white")
  ggsave("figures/14_causal/ate_all_variables_forest.svg", 
         plot = p_forest, width = 6, height = 5, bg = "white")
  
  cat("  ✓ 森林图已保存\n")
}

# ============================================================================
# 日志输出
# ============================================================================
sink("output/14_causal/ate_batch_log.txt")
cat("批量ATE估计日志\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("候选变量数: ", length(candidate_vars), "\n")
cat("成功估计数: ", sum(!is.na(ate_all$coef)), "\n")
cat("显著变量数 (p<0.05): ", sum(ate_all$significant, na.rm = TRUE), "\n\n")

cat("Top 10 最大ATE (绝对值):\n")
print(head(ate_all %>% select(variable, coef, std_error, p_value, significant), 10))
sink()

cat("\n======================================\n")
cat("批量ATE估计完成\n")
cat("======================================\n\n")
cat("结果摘要:\n")
cat("  - 候选变量: ", length(candidate_vars), "\n")
cat("  - 成功估计: ", sum(!is.na(ate_all$coef)), "\n")
cat("  - 显著变量 (p<0.05): ", sum(ate_all$significant, na.rm = TRUE), "\n\n")

cat("输出文件:\n")
cat("  - output/14_causal/ate_all_variables.csv\n")
cat("  - figures/14_causal/ate_all_variables_forest.png\n\n")

cat("✓ 脚本执行完成!\n\n")

