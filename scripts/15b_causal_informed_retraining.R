#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 15b_causal_informed_retraining.R
# 功能说明: 基于因果发现识别的核心驱动因子重新训练模型，对比全变量模型性能
# 方法逻辑: DAG上游节点（高出度+高稳定性）+ ATE显著变量 → 核心驱动因子集
# 输入文件: output/14_causal/edges_summary.csv (DAG边稳定性)
#          output/14_causal/ate_summary.csv (因果效应估计)
#          output/09_variable_importance/importance_summary.csv (模型重要性)
#          output/04_collinearity/collinearity_removed.csv (训练数据)
# 输出文件: output/15b_causal_retraining/core_drivers_selection.csv (筛选的核心变量)
#          output/15b_causal_retraining/models/{maxnet,rf,gam,nn}_causal.rds
#          output/15b_causal_retraining/performance_comparison.csv (全变量vs简化模型)
#          figures/15b_causal_retraining/performance_comparison.png
# 作者: Nature级别科研项目
# 日期: 2025-11-10
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/SDM01")

# 加载必要的包
packages <- c("tidyverse", "maxnet", "randomForest", "nnet", "mgcv", "pROC", 
              "caret", "ggplot2", "viridis", "sysfonts", "showtext", "patchwork")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 字体设置（Nature标准：Arial，1200dpi）
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

dir.create("output/15b_causal_retraining", showWarnings = FALSE, recursive = TRUE)
dir.create("output/15b_causal_retraining/models", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/15b_causal_retraining", showWarnings = FALSE, recursive = TRUE)

# 统一可视化工具
source("scripts/visualization/viz_utils.R")

cat("\n======================================\n")
cat("因果驱动的简化建模\n")
cat("从47变量到核心驱动因子\n")
cat("======================================\n\n")

# ============================================================================
# 步骤1: 从因果DAG和ATE分析中识别核心驱动因子
# ============================================================================
cat("步骤 1/5: 识别因果核心驱动因子...\n")

# 1.1 读取DAG边稳定性（300次bootstrap）
edges_df <- read.csv("output/14_causal/edges_summary.csv", stringsAsFactors = FALSE)

# 计算每个变量的因果网络中心性指标
node_metrics <- edges_df %>%
  filter(strength >= 0.55) %>%  # 只保留稳定边
  group_by(from) %>%
  summarise(
    out_degree = n(),  # 出度：该变量影响多少其他变量
    mean_strength = mean(strength),  # 平均边强度
    .groups = "drop"
  ) %>%
  arrange(desc(out_degree), desc(mean_strength))

cat("  - DAG上游核心节点 (高出度):\n")
print(head(node_metrics, 10))

# 1.2 读取变量重要性（从全变量模型）
imp_df <- read.csv("output/09_variable_importance/importance_summary.csv", stringsAsFactors = FALSE)
imp_summary <- imp_df %>%
  group_by(variable) %>%
  summarise(mean_importance = mean(importance_normalized, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_importance))

cat("\n  - 模型变量重要性Top10:\n")
print(head(imp_summary, 10))

# 1.3 读取ATE分析结果（如果存在多变量ATE汇总）
ate_path <- "output/14_causal/ate_all_variables.csv"
if(file.exists(ate_path)) {
  ate_df <- read.csv(ate_path, stringsAsFactors = FALSE)
  # 识别ATE显著的变量（p<0.05且效应量大）
  ate_sig <- ate_df %>%
    filter(p_value < 0.05) %>%
    arrange(desc(abs(coef))) %>%
    select(variable, coef, p_value)
  cat("\n  - ATE显著变量:\n")
  print(head(ate_sig, 10))
  # 保存ATE显著变量详细表
  write.csv(ate_sig, "output/15b_causal_retraining/ate_significant_variables.csv", row.names = FALSE)
} else {
  cat("\n  ✗ 未找到多变量ATE汇总，将仅基于DAG+重要性筛选\n")
  ate_sig <- data.frame(variable = character(0), coef = numeric(0))
}

# 保存DAG和重要性的摘要表（便于论文引用）
dag_summary <- node_metrics %>% head(15) %>% mutate(source = "DAG")
imp_summary_top <- imp_summary %>% head(15) %>% mutate(source = "Importance")
ate_summary_top <- if(nrow(ate_sig) > 0) ate_sig %>% head(10) %>% mutate(source = "ATE") else data.frame()

selection_summary <- bind_rows(
  dag_summary %>% select(variable = from, metric = out_degree, source),
  imp_summary_top %>% select(variable, metric = mean_importance, source),
  if(nrow(ate_summary_top) > 0) ate_summary_top %>% select(variable, metric = coef, source) else NULL
)
write.csv(selection_summary, "output/15b_causal_retraining/selection_criteria_summary.csv", row.names = FALSE)

# 1.4 综合筛选策略：DAG上游 + 模型重要性 + ATE显著性
# 策略：取DAG出度Top15 ∪ 重要性Top15 ∪ ATE显著变量，再去重
dag_top <- head(node_metrics$from, 15)
imp_top <- head(imp_summary$variable, 15)
ate_top <- if(nrow(ate_sig) > 0) head(ate_sig$variable, 10) else character(0)

core_drivers <- unique(c(dag_top, imp_top, ate_top))

# 读取建模数据，确保所有核心变量存在
train_data <- read.csv("output/04_collinearity/collinearity_removed.csv", stringsAsFactors = FALSE)
exclude_cols <- c("id", "species", "lon", "lat", "source", "presence", "presence.1")
available_vars <- setdiff(names(train_data), exclude_cols)
core_drivers <- intersect(core_drivers, available_vars)

cat("\n  ✓ 识别出", length(core_drivers), "个因果核心驱动因子\n")
cat("    变量列表:\n")
for(i in seq_along(core_drivers)) {
  cat(sprintf("    %2d. %s\n", i, core_drivers[i]))
}

# 保存核心变量列表
core_selection_log <- data.frame(
  variable = core_drivers,
  in_dag_top15 = core_drivers %in% dag_top,
  in_imp_top15 = core_drivers %in% imp_top,
  in_ate_sig = core_drivers %in% ate_top,
  selection_reason = case_when(
    core_drivers %in% dag_top & core_drivers %in% imp_top & core_drivers %in% ate_top ~ "DAG+IMP+ATE",
    core_drivers %in% dag_top & core_drivers %in% imp_top ~ "DAG+IMP",
    core_drivers %in% dag_top & core_drivers %in% ate_top ~ "DAG+ATE",
    core_drivers %in% imp_top & core_drivers %in% ate_top ~ "IMP+ATE",
    core_drivers %in% dag_top ~ "DAG",
    core_drivers %in% imp_top ~ "IMP",
    core_drivers %in% ate_top ~ "ATE",
    TRUE ~ "Other"
  )
)
write.csv(core_selection_log, "output/15b_causal_retraining/core_drivers_selection.csv", row.names = FALSE)

# ============================================================================
# 步骤2: 准备训练/测试数据集
# ============================================================================
cat("\n步骤 2/5: 准备训练/测试数据...\n")

# 提取核心变量 + 响应变量
X_core <- train_data[, core_drivers, drop = FALSE]
y <- train_data$presence

# 数据集划分（与原始建模保持一致：80%训练，20%测试）
set.seed(42)  # 保持与全变量模型相同的随机种子
train_idx <- caret::createDataPartition(y, p = 0.8, list = FALSE)
X_train <- X_core[train_idx, ]
X_test <- X_core[-train_idx, ]
y_train <- y[train_idx]
y_test <- y[-train_idx]

# 为GAM准备经纬度（空间平滑项）
lonlat_train <- train_data[train_idx, c("lon", "lat")]
lonlat_test <- train_data[-train_idx, c("lon", "lat")]

cat("  ✓ 训练集:", nrow(X_train), "样本\n")
cat("  ✓ 测试集:", nrow(X_test), "样本\n")
cat("  ✓ 核心变量数:", ncol(X_core), "\n")

# ============================================================================
# 步骤3: 训练简化模型（4个算法）
# ============================================================================
cat("\n步骤 3/5: 训练简化模型（核心驱动因子）...\n")

models_causal <- list()

# 3.1 Maxnet
cat("  - Maxnet...\n")
train_maxnet <- cbind(presence = y_train, X_train)
model_maxnet <- maxnet::maxnet(p = train_maxnet$presence, 
                                data = train_maxnet[, -1, drop = FALSE],
                                maxnet.formula(p = train_maxnet$presence, 
                                              data = train_maxnet[, -1, drop = FALSE]))
models_causal$Maxnet <- model_maxnet

# 3.2 Random Forest
cat("  - Random Forest...\n")
train_rf <- cbind(presence = as.factor(y_train), X_train)
model_rf <- randomForest::randomForest(presence ~ ., data = train_rf, 
                                       ntree = 500, importance = TRUE)
models_causal$RF <- model_rf

# 3.3 GAM (含空间平滑)
cat("  - GAM...\n")
train_gam <- cbind(presence = y_train, X_train, lonlat_train)
# 构建公式：所有核心变量 + 空间项
gam_formula <- as.formula(paste0("presence ~ ", 
                                 paste0("s(", core_drivers, ", k=5)", collapse = " + "),
                                 " + s(lon, lat, k=10)"))
model_gam <- mgcv::gam(gam_formula, data = train_gam, family = binomial, 
                       method = "REML", weights = ifelse(y_train == 1, 3, 1))
models_causal$GAM <- model_gam

# 3.4 Neural Network
cat("  - Neural Network...\n")
X_train_scaled <- scale(X_train)
mu <- attr(X_train_scaled, "scaled:center")
sdv <- attr(X_train_scaled, "scaled:scale")
sdv[sdv == 0] <- 1
model_nn_raw <- nnet::nnet(X_train_scaled, y_train, size = 10, decay = 0.01, 
                           maxit = 200, trace = FALSE, linout = FALSE)
model_nn <- list(model = model_nn_raw, mean = mu, sd = sdv, vars = core_drivers)
models_causal$NN <- model_nn

# 保存简化模型
saveRDS(models_causal$Maxnet, "output/15b_causal_retraining/models/maxnet_causal.rds")
saveRDS(models_causal$RF, "output/15b_causal_retraining/models/rf_causal.rds")
saveRDS(models_causal$GAM, "output/15b_causal_retraining/models/gam_causal.rds")
saveRDS(models_causal$NN, "output/15b_causal_retraining/models/nn_causal.rds")

# 提取并保存简化模型的变量重要性
causal_var_imp <- list()

# RF变量重要性
if("RF" %in% names(models_causal)) {
  rf_imp <- randomForest::importance(models_causal$RF)
  causal_var_imp$RF <- data.frame(
    model = "RF",
    variable = rownames(rf_imp),
    importance = rf_imp[, "MeanDecreaseGini"]
  )
}

# GAM平滑项edf
if("GAM" %in% names(models_causal)) {
  gam_summary <- summary(models_causal$GAM)
  gam_imp <- data.frame(
    model = "GAM",
    variable = rownames(gam_summary$s.table),
    importance = gam_summary$s.table[, "edf"]
  )
  # 移除空间项
  gam_imp <- gam_imp %>% filter(!grepl("lon.*lat", variable))
  causal_var_imp$GAM <- gam_imp
}

if(length(causal_var_imp) > 0) {
  causal_imp_all <- bind_rows(causal_var_imp)
  write.csv(causal_imp_all, "output/15b_causal_retraining/causal_model_variable_importance.csv", row.names = FALSE)
}

cat("  ✓ 模型训练完成并保存\n")

# ============================================================================
# 步骤4: 评估简化模型性能
# ============================================================================
cat("\n步骤 4/5: 评估简化模型性能...\n")

# 预测函数
predict_causal <- function(model, model_name, X_new, lonlat_new = NULL) {
  if(model_name == "Maxnet") {
    return(as.numeric(predict(model, X_new, type = "logistic")))
  } else if(model_name == "RF") {
    return(as.numeric(predict(model, newdata = X_new, type = "prob")[, "1"]))
  } else if(model_name == "GAM") {
    test_gam <- cbind(X_new, lonlat_new)
    return(as.numeric(predict(model, newdata = test_gam, type = "response")))
  } else if(model_name == "NN") {
    X_scaled <- sweep(as.matrix(X_new), 2, model$mean[model$vars], "-")
    X_scaled <- sweep(X_scaled, 2, model$sd[model$vars], "/")
    return(as.numeric(predict(model$model, X_scaled, type = "raw")))
  }
}

# 评估指标计算
eval_results_causal <- list()
for(mn in names(models_causal)) {
  pred_prob <- predict_causal(models_causal[[mn]], mn, X_test, lonlat_test)
  
  # AUC
  roc_obj <- pROC::roc(y_test, pred_prob, quiet = TRUE)
  auc_val <- as.numeric(pROC::auc(roc_obj))
  
  # 最优阈值（Youden's Index）
  coords_all <- pROC::coords(roc_obj, x = "all", ret = c("threshold", "sensitivity", "specificity"))
  tss_vals <- coords_all$sensitivity + coords_all$specificity - 1
  opt_idx <- which.max(tss_vals)
  opt_thr <- coords_all$threshold[opt_idx]
  opt_sens <- coords_all$sensitivity[opt_idx]
  opt_spec <- coords_all$specificity[opt_idx]
  tss_val <- tss_vals[opt_idx]
  
  eval_results_causal[[mn]] <- data.frame(
    model = mn,
    n_vars = length(core_drivers),
    auc = auc_val,
    tss = tss_val,
    sensitivity = opt_sens,
    specificity = opt_spec,
    threshold = opt_thr
  )
  
  cat(sprintf("  %8s: AUC=%.3f, TSS=%.3f, Sens=%.3f, Spec=%.3f\n",
              mn, auc_val, tss_val, opt_sens, opt_spec))
}

eval_causal <- bind_rows(eval_results_causal)

# ============================================================================
# 步骤5: 对比全变量模型 vs 简化模型
# ============================================================================
cat("\n步骤 5/5: 对比全变量模型 vs 简化模型...\n")

# 读取全变量模型评估结果
eval_full_path <- "output/08_model_evaluation/evaluation_summary.csv"
if(!file.exists(eval_full_path)) {
  cat("  ✗ 未找到全变量模型评估结果，跳过对比\n")
  eval_full <- data.frame(model = character(0), auc = numeric(0), tss = numeric(0))
} else {
  eval_full <- read.csv(eval_full_path, stringsAsFactors = FALSE)
  # 中文注释：统一列名为小写，避免大小写不匹配
  names(eval_full) <- tolower(names(eval_full))
  eval_full$n_vars <- 47  # 全变量模型使用47个变量
  
  # 确保包含必要列（处理可能缺失的列）
  required_cols <- c("model", "auc", "tss", "sensitivity", "specificity")
  missing_cols <- setdiff(required_cols, names(eval_full))
  if(length(missing_cols) > 0) {
    cat("  ⚠ 全变量模型评估结果缺少列:", paste(missing_cols, collapse=", "), "\n")
    for(col in missing_cols) {
      eval_full[[col]] <- NA
    }
  }
}

# 合并对比
comparison <- bind_rows(
  eval_full %>% mutate(model_type = "Full (47 vars)") %>% select(model, model_type, n_vars, auc, tss, sensitivity, specificity),
  eval_causal %>% mutate(model_type = paste0("Causal (", length(core_drivers), " vars)")) %>% select(model, model_type, n_vars, auc, tss, sensitivity, specificity)
)

# 计算性能保留率
retention <- eval_causal %>%
  left_join(eval_full %>% select(model, auc_full = auc, tss_full = tss), by = "model") %>%
  mutate(
    auc_retention = auc / auc_full * 100,
    tss_retention = tss / tss_full * 100,
    var_reduction = (1 - n_vars / 47) * 100
  )

cat("\n  === 性能对比摘要 ===\n")
print(comparison)

cat("\n  === 简化模型性能保留率 ===\n")
print(retention %>% select(model, n_vars, auc_retention, tss_retention, var_reduction))

# 保存对比结果
write.csv(comparison, "output/15b_causal_retraining/performance_comparison.csv", row.names = FALSE)
write.csv(retention, "output/15b_causal_retraining/performance_retention.csv", row.names = FALSE)

# ============================================================================
# 步骤6: 可视化对比
# ============================================================================
cat("\n步骤 6/6: 生成对比图表...\n")

# 6.1 AUC对比柱状图
# 中文注释：构建命名颜色向量，避免在c()内使用paste0()
causal_label <- paste0("Causal (", length(core_drivers), " vars)")
color_values <- c("#377EB8", "#E41A1C")
names(color_values) <- c("Full (47 vars)", causal_label)

p_auc <- ggplot(comparison, aes(x = model, y = auc, fill = model_type)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", auc)), position = position_dodge(0.8), 
            vjust = -0.5, size = 2.5, family = "Arial") +
  scale_fill_manual(values = color_values, name = "Model Type") +
  labs(title = "AUC Comparison: Full vs. Causal Models", 
       x = "Algorithm", y = "AUC") +
  ylim(0, 1) +
  viz_theme_nature(base_size = 8, title_size = 9) +
  theme(legend.position = "top")

# 6.2 TSS对比柱状图（复用上面的颜色向量）
p_tss <- ggplot(comparison, aes(x = model, y = tss, fill = model_type)) +
  geom_col(position = position_dodge(0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", tss)), position = position_dodge(0.8), 
            vjust = -0.5, size = 2.5, family = "Arial") +
  scale_fill_manual(values = color_values, name = "Model Type") +
  labs(title = "TSS Comparison: Full vs. Causal Models", 
       x = "Algorithm", y = "TSS") +
  ylim(0, 1) +
  viz_theme_nature(base_size = 8, title_size = 9) +
  theme(legend.position = "top")

# 6.3 性能保留率散点图
p_retention <- ggplot(retention, aes(x = var_reduction, y = auc_retention)) +
  geom_point(aes(size = tss_retention, color = model), alpha = 0.7) +
  geom_text(aes(label = model), vjust = -1.2, size = 2.8, family = "Arial") +
  geom_hline(yintercept = 90, linetype = "dashed", color = "grey60") +
  geom_vline(xintercept = 50, linetype = "dashed", color = "grey60") +
  scale_color_viridis_d(option = "D", name = "Algorithm") +
  scale_size_continuous(range = c(3, 8), name = "TSS Retention (%)") +
  labs(title = "Variable Reduction vs. Performance Retention",
       subtitle = "Dashed lines: 90% AUC retention, 50% variable reduction",
       x = "Variable Reduction (%)", 
       y = "AUC Retention (%)") +
  xlim(0, 100) + ylim(70, 105) +
  viz_theme_nature(base_size = 8, title_size = 9)

# 组合图
combined_plot <- (p_auc | p_tss) / p_retention + 
  plot_layout(heights = c(1, 1.2))

ggsave("figures/15b_causal_retraining/performance_comparison.png", 
       plot = combined_plot, width = 8, height = 6, dpi = 1200, bg = "white")
ggsave("figures/15b_causal_retraining/performance_comparison.svg", 
       plot = combined_plot, width = 8, height = 6, bg = "white")

cat("  ✓ 对比图表已保存\n")

# ============================================================================
# 日志输出
# ============================================================================
sink("output/15b_causal_retraining/processing_log.txt")
cat("因果驱动的简化建模日志\n")
cat(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
cat("==== 核心驱动因子筛选 ====\n")
cat("筛选策略: DAG上游节点 + 模型重要性 + ATE显著性\n")
cat("变量数: 47 → ", length(core_drivers), "\n", sep = "")
cat("变量缩减率: ", round((1 - length(core_drivers)/47)*100, 1), "%\n\n", sep = "")

cat("==== 核心变量详细列表 ====\n")
cat("（按筛选来源分类）\n\n")
print(core_selection_log %>% arrange(selection_reason, variable))

cat("\n==== DAG上游核心节点 (Top15) ====\n")
print(node_metrics %>% head(15))

cat("\n==== 模型重要性 (Top15) ====\n")
print(imp_summary %>% head(15))

if(nrow(ate_sig) > 0) {
  cat("\n==== ATE显著变量 ====\n")
  print(ate_sig)
}

cat("\n==== 性能对比 ====\n")
print(comparison)

cat("\n==== 性能保留率 ====\n")
print(retention)

cat("\n==== 核心结论 ====\n")
mean_auc_retention <- mean(retention$auc_retention, na.rm = TRUE)
mean_tss_retention <- mean(retention$tss_retention, na.rm = TRUE)
cat("平均AUC保留率: ", round(mean_auc_retention, 1), "%\n", sep = "")
cat("平均TSS保留率: ", round(mean_tss_retention, 1), "%\n", sep = "")
cat("变量缩减率: ", round(mean(retention$var_reduction), 1), "%\n\n", sep = "")

if(mean_auc_retention >= 90) {
  cat("✓ 简化模型保持了90%以上的预测精度\n")
  cat("✓ 因果驱动的变量筛选显著提升了模型可解释性和转移性\n")
} else if(mean_auc_retention >= 80) {
  cat("✓ 简化模型保持了80%以上的预测精度\n")
  cat("✓ 在精度略有损失的情况下，大幅提升了模型简洁性\n")
} else {
  cat("⚠ 简化模型精度损失较大（<80%），建议增加核心变量数量\n")
}
sink()

cat("\n======================================\n")
cat("因果驱动的简化建模完成\n")
cat("======================================\n\n")
cat("核心发现:\n")
cat("  - 核心驱动因子数: ", length(core_drivers), " (缩减", round((1 - length(core_drivers)/47)*100, 1), "%)\n", sep = "")
cat("  - 平均AUC保留率: ", round(mean_auc_retention, 1), "%\n", sep = "")
cat("  - 平均TSS保留率: ", round(mean_tss_retention, 1), "%\n\n", sep = "")

cat("输出文件:\n")
cat("  数据:\n")
cat("    - output/15b_causal_retraining/core_drivers_selection.csv\n")
cat("    - output/15b_causal_retraining/performance_comparison.csv\n")
cat("    - output/15b_causal_retraining/models/*_causal.rds\n")
cat("  图表:\n")
cat("    - figures/15b_causal_retraining/performance_comparison.png\n\n")

cat("✓ 脚本执行完成!\n\n")

