#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 16_screening_benchmark.R
# 功能说明: 变量筛选方法基准对照实验
#           比较 5 种变量选择策略的 SDM 预测性能
#
# 对比方案:
#   1. Full       — 全部 47 个变量 (基线)
#   2. Causal     — DAG拓扑 + ATE显著性 + 重要性联合筛选 (本文方法)
#   3. VIF        — 方差膨胀因子逐步剔除 (传统方法)
#   4. LASSO      — L1正则化自动选择 (机器学习方法)
#   5. Random     — 随机选取与Causal等量变量 (零假设基线)
#
# 评估维度:
#   - 标准测试集 AUC / TSS
#   - 5折空间交叉验证 (Spatial CV) — 检验可转移性
#
# 输入: output/04_collinearity/collinearity_removed.csv
#       output/15b_causal_retraining/core_drivers_selection.csv
# 输出: output/16_benchmark/benchmark_results.csv
#       figures/16_benchmark/benchmark_comparison.png
#
# 作者: CausalSDM 项目
# 日期: 2025
# ==============================================================================

# === 初始化 ===
rm(list = ls())
gc()
setwd("E:/CausalSDMs")

packages <- c(
    "tidyverse", "maxnet", "randomForest", "pROC", "caret",
    "glmnet", # LASSO
    "usdm", # VIF
    "ggplot2", "viridis", "patchwork",
    "sysfonts", "showtext"
)

for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
}

# 字体
try(
    {
        sysfonts::font_add(
            family = "Arial",
            regular = "C:/Windows/Fonts/arial.ttf",
            bold = "C:/Windows/Fonts/arialbd.ttf"
        )
        showtext::showtext_opts(dpi = 2400)
        showtext::showtext_auto(enable = TRUE)
    },
    silent = TRUE
)

dir.create("output/16_benchmark", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/16_benchmark", showWarnings = FALSE, recursive = TRUE)

source("scripts/visualization/viz_utils.R")

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("     变量筛选方法基准对照实验\n")
cat("     Full vs Causal vs VIF vs LASSO vs Random\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# ==============================================================================
# 1. 数据准备
# ==============================================================================
cat("步骤 1/6: 数据准备...\n")

dat <- read.csv("output/04_collinearity/collinearity_removed.csv",
    stringsAsFactors = FALSE
)
exclude_cols <- c("id", "species", "lon", "lat", "source", "presence", "presence.1")
env_vars <- setdiff(colnames(dat), exclude_cols)
Y <- dat$presence
coords <- dat[, c("lon", "lat")]
X_all <- dat[, env_vars, drop = FALSE]

cat("  样本数:", nrow(dat), "\n")
cat("  环境变量数:", length(env_vars), "\n")
cat("  出现点:", sum(Y == 1), " | 背景点:", sum(Y == 0), "\n")

# 统一训练-测试划分 (与全变量模型一致)
set.seed(42)
train_idx <- caret::createDataPartition(Y, p = 0.8, list = FALSE)[, 1]
X_train_all <- X_all[train_idx, ]
X_test_all <- X_all[-train_idx, ]
y_train <- Y[train_idx]
y_test <- Y[-train_idx]
coords_train <- coords[train_idx, ]
coords_test <- coords[-train_idx, ]

cat("  训练集:", length(y_train), " | 测试集:", length(y_test), "\n\n")

# ==============================================================================
# 2. 五种变量选择策略
# ==============================================================================
cat("步骤 2/6: 执行变量选择策略...\n\n")

strategies <- list()

# ---- 2.1 Full (47 vars) ----
strategies[["Full"]] <- list(
    vars = env_vars,
    label = paste0("Full (", length(env_vars), " vars)"),
    color = "#636363"
)
cat("  [Full] 全部", length(env_vars), "变量\n")

# ---- 2.2 Causal (本文方法) ----
causal_path <- "output/15b_causal_retraining/core_drivers_selection.csv"
if (file.exists(causal_path)) {
    causal_df <- read.csv(causal_path, stringsAsFactors = FALSE)
    causal_vars <- intersect(causal_df$variable, env_vars)
} else {
    # 如果文件不存在，从DAG+ATE+Importance三源重建
    cat("    ⚠ 未找到因果筛选结果，将从原始数据重建...\n")

    # DAG Top15
    edges_df <- read.csv("output/14_causal/edges_summary.csv", stringsAsFactors = FALSE)
    dag_top <- edges_df %>%
        filter(strength >= 0.55) %>%
        group_by(from) %>%
        summarise(out_degree = n(), .groups = "drop") %>%
        arrange(desc(out_degree)) %>%
        head(15) %>%
        pull(from)

    # Importance Top15
    imp_df <- read.csv("output/09_variable_importance/importance_summary.csv", stringsAsFactors = FALSE)
    imp_top <- imp_df %>%
        group_by(variable) %>%
        summarise(mean_imp = mean(importance_normalized, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(mean_imp)) %>%
        head(15) %>%
        pull(variable)

    # ATE significant
    ate_path <- "output/14_causal/ate_all_variables.csv"
    if (file.exists(ate_path)) {
        ate_df <- read.csv(ate_path, stringsAsFactors = FALSE)
        ate_top <- ate_df %>%
            filter(p_value < 0.05) %>%
            pull(variable)
    } else {
        ate_top <- character(0)
    }

    causal_vars <- intersect(unique(c(dag_top, imp_top, ate_top)), env_vars)
}

strategies[["Causal"]] <- list(
    vars = causal_vars,
    label = paste0("Causal (", length(causal_vars), " vars)"),
    color = "#E41A1C"
)
cat("  [Causal] 因果筛选:", length(causal_vars), "变量\n")

# ---- 2.3 VIF (传统方法) ----
cat("  [VIF] 执行VIF逐步剔除...\n")
tryCatch(
    {
        vif_result <- usdm::vifstep(X_train_all, th = 10)
        vif_vars <- vif_result@results$Variables
    },
    error = function(e) {
        cat("    ⚠ VIF计算失败, 使用correlation-based替代\n")
        # 后备方案: 用相关性阈值|r| < 0.8筛选
        cor_mat <- cor(X_train_all, use = "complete.obs")
        vif_vars <<- colnames(X_train_all)
        repeat {
            cor_upper <- cor_mat
            cor_upper[lower.tri(cor_upper, diag = TRUE)] <- 0
            max_cor <- max(abs(cor_upper))
            if (max_cor < 0.8 || length(vif_vars) < 5) break
            idx <- which(abs(cor_upper) == max_cor, arr.ind = TRUE)[1, ]
            remove_var <- colnames(cor_mat)[idx[2]]
            vif_vars <<- setdiff(vif_vars, remove_var)
            cor_mat <- cor_mat[-idx[2], -idx[2], drop = FALSE]
        }
    }
)

strategies[["VIF"]] <- list(
    vars = vif_vars,
    label = paste0("VIF (", length(vif_vars), " vars)"),
    color = "#377EB8"
)
cat("  [VIF] 保留:", length(vif_vars), "变量\n")

# ---- 2.4 LASSO (L1正则化) ----
cat("  [LASSO] 执行LASSO变量选择...\n")
set.seed(42)
X_lasso <- as.matrix(scale(X_train_all))
X_lasso[is.na(X_lasso)] <- 0

cv_lasso <- glmnet::cv.glmnet(
    x = X_lasso, y = y_train,
    family = "binomial", alpha = 1,
    nfolds = 10, type.measure = "auc"
)

# 使用 lambda.1se (更简约的模型)
lasso_coefs <- coef(cv_lasso, s = "lambda.1se")
lasso_vars <- rownames(lasso_coefs)[which(lasso_coefs != 0)]
lasso_vars <- setdiff(lasso_vars, "(Intercept)")
lasso_vars <- intersect(lasso_vars, env_vars)

# 如果LASSO选出变量太少(<5), 使用lambda.min
if (length(lasso_vars) < 5) {
    lasso_coefs <- coef(cv_lasso, s = "lambda.min")
    lasso_vars <- rownames(lasso_coefs)[which(lasso_coefs != 0)]
    lasso_vars <- setdiff(lasso_vars, "(Intercept)")
    lasso_vars <- intersect(lasso_vars, env_vars)
    cat("    (使用lambda.min, 选出更多变量)\n")
}

strategies[["LASSO"]] <- list(
    vars = lasso_vars,
    label = paste0("LASSO (", length(lasso_vars), " vars)"),
    color = "#4DAF4A"
)
cat("  [LASSO] 保留:", length(lasso_vars), "变量\n")

# ---- 2.5 Random (零假设基线) ----
# 随机选取与Causal等量的变量, 重复10次取平均
n_causal <- length(causal_vars)
set.seed(123)
random_vars <- sample(env_vars, n_causal) # 取一组用于模型训练

strategies[["Random"]] <- list(
    vars = random_vars,
    label = paste0("Random (", n_causal, " vars)"),
    color = "#984EA3"
)
cat("  [Random] 随机选取:", n_causal, "变量 (与Causal等量)\n\n")

# 打印变量选择摘要
cat("  === 变量选择摘要 ===\n")
for (s in names(strategies)) {
    cat(sprintf("  %-10s: %2d 变量\n", s, length(strategies[[s]]$vars)))
}
cat("\n")

# ==============================================================================
# 3. 训练模型与评估 (标准测试集)
# ==============================================================================
cat("步骤 3/6: 训练4种算法 × 5种策略 = 20个模型...\n\n")

# 通用训练-评估函数
train_and_eval <- function(X_tr, y_tr, X_te, y_te, strategy_name) {
    results <- list()

    # --- Maxent ---
    tryCatch(
        {
            train_df <- cbind(presence = y_tr, X_tr)
            model_mx <- maxnet::maxnet(
                p = train_df$presence,
                data = train_df[, -1, drop = FALSE],
                maxnet.formula(p = train_df$presence, data = train_df[, -1, drop = FALSE])
            )
            pred_mx <- as.numeric(predict(model_mx, X_te, type = "logistic"))
            roc_mx <- pROC::roc(y_te, pred_mx, quiet = TRUE)
            auc_mx <- as.numeric(pROC::auc(roc_mx))
            coords_mx <- pROC::coords(roc_mx, x = "all", ret = c("threshold", "sensitivity", "specificity"))
            tss_vals <- coords_mx$sensitivity + coords_mx$specificity - 1
            opt_idx <- which.max(tss_vals)

            results[["Maxent"]] <- data.frame(
                strategy = strategy_name, algorithm = "Maxent",
                n_vars = ncol(X_tr),
                auc = auc_mx, tss = tss_vals[opt_idx],
                sensitivity = coords_mx$sensitivity[opt_idx],
                specificity = coords_mx$specificity[opt_idx],
                stringsAsFactors = FALSE
            )
            cat(sprintf("    %-8s × %-10s: AUC=%.3f TSS=%.3f\n", strategy_name, "Maxent", auc_mx, tss_vals[opt_idx]))
        },
        error = function(e) {
            cat(sprintf("    %-8s × %-10s: FAILED (%s)\n", strategy_name, "Maxent", e$message))
        }
    )

    # --- Random Forest ---
    tryCatch(
        {
            train_rf <- cbind(presence = as.factor(y_tr), X_tr)
            model_rf <- randomForest::randomForest(presence ~ ., data = train_rf, ntree = 500, importance = TRUE)
            pred_rf <- as.numeric(predict(model_rf, newdata = X_te, type = "prob")[, "1"])
            roc_rf <- pROC::roc(y_te, pred_rf, quiet = TRUE)
            auc_rf <- as.numeric(pROC::auc(roc_rf))
            coords_rf <- pROC::coords(roc_rf, x = "all", ret = c("threshold", "sensitivity", "specificity"))
            tss_vals <- coords_rf$sensitivity + coords_rf$specificity - 1
            opt_idx <- which.max(tss_vals)

            results[["RF"]] <- data.frame(
                strategy = strategy_name, algorithm = "RF",
                n_vars = ncol(X_tr),
                auc = auc_rf, tss = tss_vals[opt_idx],
                sensitivity = coords_rf$sensitivity[opt_idx],
                specificity = coords_rf$specificity[opt_idx],
                stringsAsFactors = FALSE
            )
            cat(sprintf("    %-8s × %-10s: AUC=%.3f TSS=%.3f\n", strategy_name, "RF", auc_rf, tss_vals[opt_idx]))
        },
        error = function(e) {
            cat(sprintf("    %-8s × %-10s: FAILED (%s)\n", strategy_name, "RF", e$message))
        }
    )

    # --- GAM ---
    tryCatch(
        {
            train_gam <- data.frame(presence = y_tr, X_tr)
            # 自动构建GAM公式 (每个变量一个平滑项)
            var_names <- colnames(X_tr)
            smooth_terms <- paste0("s(", var_names, ", k=5)")
            gam_formula <- as.formula(paste("presence ~", paste(smooth_terms, collapse = " + ")))

            model_gam <- mgcv::bam(gam_formula,
                data = train_gam, family = binomial(), method = "fREML",
                select = TRUE, gamma = 1.2
            )
            pred_gam <- as.numeric(predict(model_gam, newdata = data.frame(X_te), type = "response"))
            roc_gam <- pROC::roc(y_te, pred_gam, quiet = TRUE)
            auc_gam <- as.numeric(pROC::auc(roc_gam))
            coords_gam <- pROC::coords(roc_gam, x = "all", ret = c("threshold", "sensitivity", "specificity"))
            tss_vals <- coords_gam$sensitivity + coords_gam$specificity - 1
            opt_idx <- which.max(tss_vals)

            results[["GAM"]] <- data.frame(
                strategy = strategy_name, algorithm = "GAM",
                n_vars = ncol(X_tr),
                auc = auc_gam, tss = tss_vals[opt_idx],
                sensitivity = coords_gam$sensitivity[opt_idx],
                specificity = coords_gam$specificity[opt_idx],
                stringsAsFactors = FALSE
            )
            cat(sprintf("    %-8s × %-10s: AUC=%.3f TSS=%.3f\n", strategy_name, "GAM", auc_gam, tss_vals[opt_idx]))
        },
        error = function(e) {
            cat(sprintf("    %-8s × %-10s: FAILED (%s)\n", strategy_name, "GAM", e$message))
        }
    )

    # --- Neural Network ---
    tryCatch(
        {
            X_tr_scaled <- as.data.frame(scale(X_tr))
            X_te_scaled <- as.data.frame(scale(X_te,
                center = attr(scale(X_tr), "scaled:center"),
                scale  = attr(scale(X_tr), "scaled:scale")
            ))
            X_tr_scaled[is.na(X_tr_scaled)] <- 0
            X_te_scaled[is.na(X_te_scaled)] <- 0

            train_nn <- data.frame(presence = y_tr, X_tr_scaled)
            model_nn <- nnet::nnet(presence ~ .,
                data = train_nn, size = max(3, ncol(X_tr) %/% 5),
                decay = 5e-4, maxit = 500, trace = FALSE
            )
            pred_nn <- as.numeric(predict(model_nn, newdata = data.frame(X_te_scaled), type = "raw"))
            roc_nn <- pROC::roc(y_te, pred_nn, quiet = TRUE)
            auc_nn <- as.numeric(pROC::auc(roc_nn))
            coords_nn <- pROC::coords(roc_nn, x = "all", ret = c("threshold", "sensitivity", "specificity"))
            tss_vals <- coords_nn$sensitivity + coords_nn$specificity - 1
            opt_idx <- which.max(tss_vals)

            results[["NN"]] <- data.frame(
                strategy = strategy_name, algorithm = "NN",
                n_vars = ncol(X_tr),
                auc = auc_nn, tss = tss_vals[opt_idx],
                sensitivity = coords_nn$sensitivity[opt_idx],
                specificity = coords_nn$specificity[opt_idx],
                stringsAsFactors = FALSE
            )
            cat(sprintf("    %-8s × %-10s: AUC=%.3f TSS=%.3f\n", strategy_name, "NN", auc_nn, tss_vals[opt_idx]))
        },
        error = function(e) {
            cat(sprintf("    %-8s × %-10s: FAILED (%s)\n", strategy_name, "NN", e$message))
        }
    )

    return(bind_rows(results))
}

# 执行所有策略
all_results <- list()
for (s_name in names(strategies)) {
    cat(sprintf("  [%s]\n", s_name))
    s_vars <- strategies[[s_name]]$vars

    X_tr <- X_train_all[, s_vars, drop = FALSE]
    X_te <- X_test_all[, s_vars, drop = FALSE]

    res <- train_and_eval(X_tr, y_train, X_te, y_test, s_name)
    all_results[[s_name]] <- res
    cat("\n")
}

benchmark_df <- bind_rows(all_results)
write.csv(benchmark_df, "output/16_benchmark/benchmark_results.csv", row.names = FALSE)
cat("  ✓ 标准测试集评估完成\n\n")

# ==============================================================================
# 4. 空间交叉验证 (5折)
# ==============================================================================
cat("步骤 4/6: 空间交叉验证 (5折, 经度分层)...\n\n")

# 将训练数据按经度分成5个空间块
lon_breaks <- quantile(coords_train$lon, probs = seq(0, 1, 0.2))
spatial_folds <- cut(coords_train$lon, breaks = lon_breaks, labels = FALSE, include.lowest = TRUE)
# 确保没有 NA
spatial_folds[is.na(spatial_folds)] <- 1

# 只对 Causal / VIF / LASSO 三种策略做空间CV (加上Full基线)
cv_strategies <- c("Full", "Causal", "VIF", "LASSO")
spatial_cv_results <- list()

for (s_name in cv_strategies) {
    s_vars <- strategies[[s_name]]$vars
    cat(sprintf("  [%s] 空间CV...", s_name))

    fold_aucs <- c()
    fold_tss <- c()

    for (fold in 1:5) {
        val_idx <- which(spatial_folds == fold)
        tr_idx <- which(spatial_folds != fold)

        if (length(val_idx) < 10) next

        X_tr <- X_train_all[tr_idx, s_vars, drop = FALSE]
        X_va <- X_train_all[val_idx, s_vars, drop = FALSE]
        y_tr <- y_train[tr_idx]
        y_va <- y_train[val_idx]

        # 使用RF作为代表算法
        tryCatch(
            {
                rf_data <- cbind(presence = as.factor(y_tr), X_tr)
                rf_mod <- randomForest::randomForest(presence ~ ., data = rf_data, ntree = 300)
                rf_pred <- as.numeric(predict(rf_mod, newdata = X_va, type = "prob")[, "1"])
                roc_cv <- pROC::roc(y_va, rf_pred, quiet = TRUE)
                auc_cv <- as.numeric(pROC::auc(roc_cv))

                coords_cv <- pROC::coords(roc_cv, x = "all", ret = c("sensitivity", "specificity"))
                tss_cv <- max(coords_cv$sensitivity + coords_cv$specificity - 1, na.rm = TRUE)

                fold_aucs <- c(fold_aucs, auc_cv)
                fold_tss <- c(fold_tss, tss_cv)
            },
            error = function(e) {}
        )
    }

    if (length(fold_aucs) > 0) {
        spatial_cv_results[[s_name]] <- data.frame(
            strategy = s_name, algorithm = "RF_SpatialCV",
            n_vars = length(s_vars),
            auc_mean = mean(fold_aucs), auc_sd = sd(fold_aucs),
            tss_mean = mean(fold_tss), tss_sd = sd(fold_tss),
            n_folds = length(fold_aucs),
            stringsAsFactors = FALSE
        )
        cat(sprintf(
            " AUC=%.3f±%.3f TSS=%.3f±%.3f\n",
            mean(fold_aucs), sd(fold_aucs), mean(fold_tss), sd(fold_tss)
        ))
    } else {
        cat(" FAILED\n")
    }
}

spatial_cv_df <- bind_rows(spatial_cv_results)
write.csv(spatial_cv_df, "output/16_benchmark/spatial_cv_results.csv", row.names = FALSE)
cat("\n  ✓ 空间交叉验证完成\n\n")

# ==============================================================================
# 5. Random Baseline: 多次重复实验
# ==============================================================================
cat("步骤 5/6: 随机基线重复实验 (10次)...\n")

random_repeat_results <- list()
for (rep_i in 1:10) {
    set.seed(rep_i * 100)
    rand_vars <- sample(env_vars, n_causal)
    X_tr <- X_train_all[, rand_vars, drop = FALSE]
    X_te <- X_test_all[, rand_vars, drop = FALSE]

    tryCatch(
        {
            rf_data <- cbind(presence = as.factor(y_train), X_tr)
            rf_mod <- randomForest::randomForest(presence ~ ., data = rf_data, ntree = 300)
            rf_pred <- as.numeric(predict(rf_mod, newdata = X_te, type = "prob")[, "1"])
            roc_r <- pROC::roc(y_test, rf_pred, quiet = TRUE)

            coords_r <- pROC::coords(roc_r, x = "all", ret = c("sensitivity", "specificity"))
            tss_r <- max(coords_r$sensitivity + coords_r$specificity - 1, na.rm = TRUE)

            random_repeat_results[[rep_i]] <- data.frame(
                rep = rep_i, n_vars = n_causal,
                auc = as.numeric(pROC::auc(roc_r)),
                tss = tss_r
            )
        },
        error = function(e) {}
    )
}

random_repeat_df <- bind_rows(random_repeat_results)
write.csv(random_repeat_df, "output/16_benchmark/random_baseline_repeats.csv", row.names = FALSE)
cat(sprintf(
    "  随机基线 (RF, %d vars, 10次): AUC=%.3f±%.3f\n\n",
    n_causal, mean(random_repeat_df$auc), sd(random_repeat_df$auc)
))

# ==============================================================================
# 6. 可视化
# ==============================================================================
cat("步骤 6/6: 生成对比图表...\n")

# 6.1 策略 × 算法 热力图 (AUC)
strategy_order <- c("Full", "Causal", "VIF", "LASSO", "Random")
algo_order <- c("Maxent", "RF", "GAM", "NN")

benchmark_df$strategy <- factor(benchmark_df$strategy, levels = strategy_order)
benchmark_df$algorithm <- factor(benchmark_df$algorithm, levels = algo_order)

p1 <- ggplot(benchmark_df, aes(x = algorithm, y = strategy, fill = auc)) +
    geom_tile(color = "white", linewidth = 0.8) +
    geom_text(aes(label = sprintf("%.3f", auc)), size = 3, fontface = "bold", family = "Arial") +
    scale_fill_gradient2(
        low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
        midpoint = median(benchmark_df$auc, na.rm = TRUE),
        name = "AUC"
    ) +
    labs(
        title = "AUC: Screening Strategy × Algorithm",
        x = "SDM Algorithm", y = "Variable Screening Strategy"
    ) +
    theme_minimal(base_size = 10) +
    theme(
        text = element_text(family = "Arial"),
        axis.text.x = element_text(angle = 0, hjust = 0.5),
        panel.grid = element_blank(),
        plot.title = element_text(face = "bold", size = 11)
    )

# 6.2 策略比较柱状图 (算法平均)
strategy_means <- benchmark_df %>%
    group_by(strategy) %>%
    summarise(
        auc_mean = mean(auc, na.rm = TRUE),
        auc_se = sd(auc, na.rm = TRUE) / sqrt(n()),
        tss_mean = mean(tss, na.rm = TRUE),
        tss_se = sd(tss, na.rm = TRUE) / sqrt(n()),
        n_vars = first(n_vars),
        .groups = "drop"
    )

color_map <- c(
    "Full" = "#636363", "Causal" = "#E41A1C",
    "VIF"  = "#377EB8", "LASSO"  = "#4DAF4A", "Random" = "#984EA3"
)

p2 <- ggplot(strategy_means, aes(x = strategy, y = auc_mean, fill = strategy)) +
    geom_col(width = 0.7, alpha = 0.9) +
    geom_errorbar(aes(ymin = auc_mean - auc_se, ymax = auc_mean + auc_se),
        width = 0.2, linewidth = 0.5
    ) +
    geom_text(aes(label = sprintf("%.3f\n(%d vars)", auc_mean, n_vars)),
        vjust = -0.3, size = 2.8, family = "Arial"
    ) +
    scale_fill_manual(values = color_map, guide = "none") +
    labs(
        title = "Mean AUC Across Algorithms",
        subtitle = "Error bars: ±1 SE across 4 algorithms",
        x = "Screening Strategy", y = "Mean AUC"
    ) +
    ylim(0.75, 1.0) +
    theme_minimal(base_size = 10) +
    theme(
        text = element_text(family = "Arial"),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, color = "grey40")
    )

# 6.3 空间CV对比
if (nrow(spatial_cv_df) > 0) {
    spatial_cv_df$strategy <- factor(spatial_cv_df$strategy, levels = cv_strategies)

    p3 <- ggplot(spatial_cv_df, aes(x = strategy, y = auc_mean, fill = strategy)) +
        geom_col(width = 0.7, alpha = 0.9) +
        geom_errorbar(aes(ymin = auc_mean - auc_sd, ymax = auc_mean + auc_sd),
            width = 0.2, linewidth = 0.5
        ) +
        geom_text(aes(label = sprintf("%.3f±%.3f", auc_mean, auc_sd)),
            vjust = -0.3, size = 2.8, family = "Arial"
        ) +
        scale_fill_manual(values = color_map, guide = "none") +
        labs(
            title = "Spatial Cross-validation (5-fold, Longitude-stratified)",
            subtitle = "RF algorithm | Error bars: ±1 SD across folds",
            x = "Screening Strategy", y = "Mean AUC"
        ) +
        ylim(0.5, 1.0) +
        theme_minimal(base_size = 10) +
        theme(
            text = element_text(family = "Arial"),
            plot.title = element_text(face = "bold", size = 11),
            plot.subtitle = element_text(size = 8, color = "grey40")
        )
} else {
    p3 <- ggplot() +
        theme_void() +
        ggtitle("Spatial CV: No data")
}

# 6.4 维度-性能权衡散点图
p4 <- ggplot(strategy_means, aes(x = n_vars, y = auc_mean, color = strategy)) +
    geom_point(size = 5, alpha = 0.9) +
    geom_text(aes(label = strategy), vjust = -1.2, size = 3, family = "Arial") +
    geom_hline(
        yintercept = strategy_means$auc_mean[strategy_means$strategy == "Full"],
        linetype = "dashed", color = "grey50", alpha = 0.6
    ) +
    scale_color_manual(values = color_map, guide = "none") +
    labs(
        title = "Parsimony-Performance Trade-off",
        subtitle = "Dashed line: Full model baseline",
        x = "Number of Variables", y = "Mean AUC"
    ) +
    theme_minimal(base_size = 10) +
    theme(
        text = element_text(family = "Arial"),
        plot.title = element_text(face = "bold", size = 11),
        plot.subtitle = element_text(size = 8, color = "grey40")
    )

# 组合
combined <- (p1 | p2) / (p3 | p4) +
    plot_annotation(
        title = "Variable Screening Benchmark: Causal vs Traditional Methods",
        theme = theme(
            plot.title = element_text(face = "bold", size = 14, family = "Arial", hjust = 0.5)
        )
    )

ggsave("figures/16_benchmark/benchmark_comparison.png",
    plot = combined, width = 14, height = 10, dpi = 300, bg = "white"
)
ggsave("figures/16_benchmark/benchmark_comparison.svg",
    plot = combined, width = 14, height = 10, bg = "white"
)

cat("  ✓ 图表已保存\n\n")

# ==============================================================================
# 最终汇总
# ==============================================================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("                 基准对照实验结果汇总\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

cat("  === 标准测试集 (各算法平均) ===\n")
print(strategy_means %>% dplyr::select(strategy, n_vars, auc_mean, tss_mean) %>% as.data.frame())

cat("\n  === 空间交叉验证 (RF算法) ===\n")
if (nrow(spatial_cv_df) > 0) {
    print(spatial_cv_df %>% dplyr::select(strategy, n_vars, auc_mean, auc_sd, tss_mean) %>% as.data.frame())
}

cat("\n  === 随机基线 (RF, 10次重复) ===\n")
cat(sprintf("  AUC: %.3f ± %.3f\n", mean(random_repeat_df$auc), sd(random_repeat_df$auc)))
cat(sprintf("  TSS: %.3f ± %.3f\n\n", mean(random_repeat_df$tss), sd(random_repeat_df$tss)))

# 保存变量选择详细列表
var_selection_detail <- data.frame(
    variable = env_vars,
    in_Full = TRUE,
    in_Causal = env_vars %in% strategies[["Causal"]]$vars,
    in_VIF = env_vars %in% strategies[["VIF"]]$vars,
    in_LASSO = env_vars %in% strategies[["LASSO"]]$vars,
    in_Random = env_vars %in% strategies[["Random"]]$vars
)
write.csv(var_selection_detail, "output/16_benchmark/variable_selection_detail.csv", row.names = FALSE)

cat("输出文件:\n")
cat("  output/16_benchmark/benchmark_results.csv\n")
cat("  output/16_benchmark/spatial_cv_results.csv\n")
cat("  output/16_benchmark/random_baseline_repeats.csv\n")
cat("  output/16_benchmark/variable_selection_detail.csv\n")
cat("  figures/16_benchmark/benchmark_comparison.png\n")
cat("  figures/16_benchmark/benchmark_comparison.svg\n\n")

cat("✓ 基准对照实验完成!\n")
