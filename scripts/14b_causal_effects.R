#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 14b_causal_effects.R
# 功能说明: 估计平均处理效应ATE与个体异质效应CATE，并输出图表与表格
# 方法: Double Machine Learning (DoubleML) 与 Generalized Random Forest (grf)
# 输入文件: output/04_collinearity/collinearity_removed.csv
# 输出文件: output/14_causal/ate_summary.csv, cate_summary.csv
#          figures/14_causal/ate_forest.png, cate_importance.png
# 作者: Nature级别科研项目
# 日期: 2025-10-24
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/SDM01")

# 加载必要的包（中文注释：全部使用英文标注出图，Arial字体，1200dpi）
packages <- c("tidyverse", "DoubleML", "grf", "ranger", "mlr3", "mlr3learners", "sysfonts", "showtext")
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
cat("因果效应估计 (ATE/CATE)\n")
cat("======================================\n\n")

# 1. 读取数据与定义处理变量
cat("步骤 1/4: 读取数据与定义处理...\n")
dat <- read.csv("output/04_collinearity/collinearity_removed.csv")
exclude_cols <- c("id", "species", "lon", "lat", "source", "presence", "presence.1")
env_vars <- setdiff(colnames(dat), exclude_cols)

# 目标Y：presence（出现概率/标签）
y <- dat$presence

# 处理T：严格使用原始变量构造，不引入任何外部噪声或模拟数据
# 支持命令行参数：--treat=变量名 --cutoff=median|q75|q25|value:数值
args <- commandArgs(trailingOnly = TRUE)
arg_treat <- NA_character_
arg_cutoff <- "median"
if(length(args) > 0) {
  for(a in args) {
    if(grepl("^--treat=", a)) arg_treat <- sub("^--treat=", "", a)
    if(grepl("^--cutoff=", a)) arg_cutoff <- sub("^--cutoff=", "", a)
  }
}

if(!is.na(arg_treat) && arg_treat %in% env_vars) {
  treat_var <- arg_treat
} else {
  # 若未指定，则默认选择在变量重要性中排名靠前者（若存在该文件），否则env_vars首个
  imp_path <- "output/09_variable_importance/importance_summary.csv"
  if(file.exists(imp_path)) {
    imp_df <- read.csv(imp_path)
    treat_candidates <- imp_df %>% dplyr::group_by(variable) %>%
      dplyr::summarise(mean_imp = mean(importance_normalized, na.rm = TRUE), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(mean_imp)) %>%
      dplyr::pull(variable)
    treat_var <- intersect(treat_candidates, env_vars)[1]
    if(is.na(treat_var)) treat_var <- env_vars[1]
  } else {
    treat_var <- env_vars[1]
  }
}

# 基于原始变量分位点/固定值阈值进行二分化（不生成任何模拟样本）
vx <- dat[[treat_var]]
vx <- as.numeric(vx)
vx[!is.finite(vx)] <- NA

if (arg_cutoff == "median") {
  thr <- median(vx, na.rm = TRUE)
} else if (arg_cutoff == "q75") {
  thr <- quantile(vx, 0.75, na.rm = TRUE)
} else if (arg_cutoff == "q25") {
  thr <- quantile(vx, 0.25, na.rm = TRUE)
} else if (grepl("^value:", arg_cutoff)) {
  thr <- as.numeric(sub("^value:", "", arg_cutoff))
} else {
  thr <- median(vx, na.rm = TRUE)
}

T_var <- as.numeric(vx > thr)

X <- dat[, setdiff(env_vars, c(treat_var)), drop = FALSE]
X[is.na(X)] <- 0

cat("  - Treatment variable: ", treat_var, "; cutoff=", arg_cutoff, " (thr=", round(thr, 4), ")\n", sep = "")

# 2. ATE: Double Machine Learning（IRM 适用于二值处理）
cat("步骤 2/4: ATE (DoubleML IRM) ...\n")
df <- data.frame(y = y, d = T_var, X)
task <- DoubleML::DoubleMLData$new(df, y_col = "y", d_cols = "d")
# 中文注释：结果回归使用回归树（regr.ranger），倾向得分使用分类树（classif.ranger，输出概率）
ml_g <- mlr3::lrn("regr.ranger", num.trees = 500)
ml_m <- mlr3::lrn("classif.ranger", num.trees = 500, predict_type = "prob")
dml_irm <- DoubleML::DoubleMLIRM$new(task, ml_g = ml_g, ml_m = ml_m, n_folds = 5)
dml_irm$fit()
ate_est <- dml_irm$summary()
ate_tbl <- as.data.frame(ate_est)
write.csv(ate_tbl, "output/14_causal/ate_summary.csv", row.names = FALSE)

# 3. CATE: grf 的因果森林 (causal_forest)
cat("步骤 3/4: CATE (causal forest) ...\n")
cf <- grf::causal_forest(X = as.matrix(X), Y = as.numeric(y), W = as.numeric(T_var), num.trees = 2000)
tau_hat <- predict(cf)$predictions

cate_df <- data.frame(id = dat$id, tau_hat = tau_hat)
write.csv(cate_df, "output/14_causal/cate_summary.csv", row.names = FALSE)

# 保存模型与配置，便于空间映射重现（严格基于真实数据训练）
saveRDS(list(
  model = cf,
  treat_var = treat_var,
  cutoff = arg_cutoff,
  threshold = thr,
  features = colnames(X)
), "output/14_causal/causal_forest_model.rds")

# 4. 图件：ATE森林图与CATE分布（ggplot，Nature风格）
cat("步骤 4/4: 图件...\n")

nm <- colnames(ate_tbl)
nm_norm <- tolower(gsub("[^a-z0-9]+", "", nm))
idx_coef <- which(nm_norm %in% c("coef", "estimate", "theta"))
idx_se   <- which(nm_norm %in% c("stderr", "se", "stderr", "stderror"))
if(length(idx_coef) == 0) idx_coef <- 1
if(length(idx_se) == 0) idx_se <- min(which(sapply(ate_tbl, is.numeric) & (seq_along(nm) != idx_coef)))
est <- as.numeric(ate_tbl[1, idx_coef])
se  <- as.numeric(ate_tbl[1, idx_se])
ci <- 1.96 * se

df_ate <- data.frame(term = paste0("T: ", treat_var), estimate = est, ymin = est - ci, ymax = est + ci)
p_ate <- ggplot(df_ate, aes(x = term, y = estimate, ymin = ymin, ymax = ymax)) +
  geom_pointrange(color = "#377EB8", linewidth = 0.5) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  labs(title = "ATE (DoubleML)", x = "", y = "ATE Estimate") +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid = element_line(color = "grey90", linewidth = 0.2),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    plot.title = element_text(face = "bold")
  )

ggsave("figures/14_causal/ate_forest.png", plot = p_ate, width = 3.6, height = 2.7, units = "in", dpi = 1200, bg = "white")

df_cate <- data.frame(tau_hat = tau_hat)
p_cate <- ggplot(df_cate, aes(x = tau_hat)) +
  geom_histogram(bins = 50, fill = "#4DAF4A", color = "white", linewidth = 0.2) +
  labs(title = "CATE Distribution (causal forest)", x = "Estimated CATE", y = "Count") +
  theme_minimal(base_family = "Arial") +
  theme(
    panel.grid = element_line(color = "grey90", linewidth = 0.2),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.3),
    plot.title = element_text(face = "bold")
  )

ggsave("figures/14_causal/cate_distribution.png", plot = p_cate, width = 3.6, height = 2.7, units = "in", dpi = 1200, bg = "white")

cat("\n======================================\n")
cat("因果效应估计完成\n")
cat("======================================\n\n")

cat("✓ ATE: output/14_causal/ate_summary.csv\n")
cat("✓ CATE: output/14_causal/cate_summary.csv\n")
cat("✓ 图件: figures/14_causal/ate_forest.png / cate_distribution.png\n\n")


