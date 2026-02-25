#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 17_dagnet_model.R
# 功能说明: DAG-Net — 因果结构化神经网络 SDM
#
# 核心架构 (3个原创设计):
#   1. 网络层 = 因果层: DAG拓扑排序决定层0/层1/层2
#   2. 残差连接 = 因果信息传递: h₁ = MLP(h₀ ⊕ x_clim) + h₀
#   3. ATE初始化输入权重: 因果效应作为先验注入模型
#
# 实验设计:
#   A. DAG-Net (完整): 因果结构 + ATE初始化
#   B. DAG-Net (无ATE): 因果结构, 无ATE初始化 (消融)
#   C. Flat NN: 同参数量, 无因果结构 (消融基线)
#   D. 传统SDM: RF, Maxent, BRT, GAM (对比基线)
#
# 依赖: torch, ranger, maxnet, gbm, mgcv, pROC, caret, tidyverse
# ==============================================================================

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

# ---- 加载依赖 ----
packages <- c(
    "tidyverse", "maxnet", "ranger", "gbm", "mgcv",
    "pROC", "caret", "ggplot2", "patchwork"
)
for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
}
if (!require(torch, quietly = TRUE)) {
    install.packages("torch")
    library(torch)
    if (!torch_is_installed()) torch::install_torch()
}

dir.create("output/17_dagnet", recursive = TRUE, showWarnings = FALSE)
dir.create("figures/17_dagnet", recursive = TRUE, showWarnings = FALSE)

cat("======================================================================\n")
cat("           DAG-Net: Causally-Structured Neural Network SDM\n")
cat("======================================================================\n\n")

# ==============================================================================
# 1. 加载数据 + 因果层级分组
# ==============================================================================
cat("步骤 1/7: 数据加载与因果层级分组...\n")

dat <- read.csv("output/04_collinearity/collinearity_removed.csv", stringsAsFactors = FALSE)
groups <- read.csv("output/14_causal/variable_groups_validated.csv", stringsAsFactors = FALSE)
ate_df <- read.csv("output/14_causal/ate_all_variables.csv", stringsAsFactors = FALSE)

# 按因果层级分组 (由DAG拓扑排序自动确定)
topo_vars <- intersect(groups$var[groups$group == "G1_TopoSlopeFlow"], colnames(dat))
clim_vars <- intersect(groups$var[groups$group == "G2_HydroClimate"], colnames(dat))
surf_vars <- intersect(groups$var[groups$group == "G3_LandSoil"], colnames(dat))
all_vars <- c(topo_vars, clim_vars, surf_vars)

n_topo <- length(topo_vars)
n_clim <- length(clim_vars)
n_surf <- length(surf_vars)
n_vars <- length(all_vars)

cat(sprintf(
    "  因果层级: Layer0(根)=%d | Layer1(中介)=%d | Layer2(终端)=%d | 总计=%d\n",
    n_topo, n_clim, n_surf, n_vars
))

# ATE权重提取
ate_weights <- setNames(rep(0.1, n_vars), all_vars)
for (i in seq_len(nrow(ate_df))) {
    v <- ate_df$variable[i]
    if (v %in% all_vars) ate_weights[v] <- abs(ate_df$coef[i])
}
ate_weights <- 0.1 + 0.9 * (ate_weights - min(ate_weights)) /
    (max(ate_weights) - min(ate_weights) + 1e-8)

cat("  ATE权重 (Top 5):\n")
top_ate <- sort(ate_weights, decreasing = TRUE)[1:5]
for (i in seq_along(top_ate)) cat(sprintf("    %s: %.4f\n", names(top_ate)[i], top_ate[i]))

# ==============================================================================
# 2. 数据准备
# ==============================================================================
cat("\n步骤 2/7: 数据准备...\n")

Y <- dat$presence
X <- dat[, all_vars, drop = FALSE]
X_means <- colMeans(X, na.rm = TRUE)
X_sds <- apply(X, 2, sd, na.rm = TRUE)
X_sds[X_sds < 1e-10] <- 1
X_scaled <- as.data.frame(scale(X, center = X_means, scale = X_sds))
X_scaled[is.na(X_scaled)] <- 0

set.seed(42)
train_idx <- caret::createDataPartition(Y, p = 0.8, list = FALSE)[, 1]
X_train <- X_scaled[train_idx, ]
y_train <- Y[train_idx]
X_test <- X_scaled[-train_idx, ]
y_test <- Y[-train_idx]

set.seed(123)
val_idx <- sample(seq_len(nrow(X_train)), size = round(0.2 * nrow(X_train)))
X_val <- X_train[val_idx, ]
y_val <- y_train[val_idx]
X_tr <- X_train[-val_idx, ]
y_tr <- y_train[-val_idx]

idx_topo <- match(topo_vars, all_vars)
idx_clim <- match(clim_vars, all_vars)
idx_surf <- match(surf_vars, all_vars)

pos_weight <- sum(y_tr == 0) / max(sum(y_tr == 1), 1)
cat(sprintf("  训练: %d | 验证: %d | 测试: %d\n", nrow(X_tr), nrow(X_val), nrow(X_test)))
cat(sprintf("  正类权重: %.2f\n", pos_weight))

# ==============================================================================
# 3. 模型架构定义
# ==============================================================================
cat("\n步骤 3/7: 定义模型架构...\n")

# ---- Swish 激活 ----
swish_module <- nn_module(
    "Swish",
    forward = function(x) x * torch_sigmoid(x)
)

# ---- DAG-Net: 因果结构化神经网络 ----
# 核心设计:
#   Layer 0: MLP(root_vars) → h₀
#   Layer 1: MLP(h₀ ⊕ mediator_vars) + h₀ → h₁   (残差 = 因果流)
#   Layer 2: MLP(h₁ ⊕ response_vars) + h₁ → h₂   (残差 = 因果流)
#   Output:  Head(h₂) → P(presence)

DAGNet <- nn_module(
    "DAGNet",
    initialize = function(n_topo, n_clim, n_surf, hidden = 128,
                          dropout = 0.2, ate_topo = NULL, ate_clim = NULL, ate_surf = NULL) {
        # ATE缩放因子 (因果先验)
        if (!is.null(ate_topo)) {
            self$ate_scale_topo <- nn_parameter(torch_tensor(ate_topo, dtype = torch_float()))
        } else {
            self$ate_scale_topo <- nn_parameter(torch_ones(n_topo))
        }
        if (!is.null(ate_clim)) {
            self$ate_scale_clim <- nn_parameter(torch_tensor(ate_clim, dtype = torch_float()))
        } else {
            self$ate_scale_clim <- nn_parameter(torch_ones(n_clim))
        }
        if (!is.null(ate_surf)) {
            self$ate_scale_surf <- nn_parameter(torch_tensor(ate_surf, dtype = torch_float()))
        } else {
            self$ate_scale_surf <- nn_parameter(torch_ones(n_surf))
        }

        # Layer 0: 因果根节点 (地形/水文拓扑)
        self$layer0 <- nn_sequential(
            nn_linear(n_topo, hidden),
            nn_batch_norm1d(hidden),
            swish_module(),
            nn_dropout(dropout)
        )

        # Layer 1: 因果中介节点 (水文气候, 接收h₀)
        self$layer1 <- nn_sequential(
            nn_linear(hidden + n_clim, hidden),
            nn_batch_norm1d(hidden),
            swish_module(),
            nn_dropout(dropout)
        )

        # Layer 2: 因果终端节点 (地表+土壤, 接收h₁)
        self$layer2 <- nn_sequential(
            nn_linear(hidden + n_surf, hidden),
            nn_batch_norm1d(hidden),
            swish_module(),
            nn_dropout(dropout)
        )

        # Output head
        self$head <- nn_sequential(
            nn_linear(hidden, 64),
            swish_module(),
            nn_dropout(dropout * 0.5),
            nn_linear(64, 32),
            swish_module(),
            nn_linear(32, 1)
        )
    },
    forward = function(x_topo, x_clim, x_surf) {
        # ATE缩放 (因果先验注入)
        x_topo <- x_topo * self$ate_scale_topo$unsqueeze(1)
        x_clim <- x_clim * self$ate_scale_clim$unsqueeze(1)
        x_surf <- x_surf * self$ate_scale_surf$unsqueeze(1)

        # Layer 0: 因果根节点
        h0 <- self$layer0(x_topo)

        # Layer 1: 因果中介 + 残差连接 (因果信息传递)
        h1 <- self$layer1(torch_cat(list(h0, x_clim), dim = 2)) + h0

        # Layer 2: 因果终端 + 残差连接 (因果信息传递)
        h2 <- self$layer2(torch_cat(list(h1, x_surf), dim = 2)) + h1

        # Output
        self$head(h2)
    }
)

# ---- Flat NN: 无因果结构基线 (参数量匹配) ----
FlatNN <- nn_module(
    "FlatNN",
    initialize = function(n_input, hidden = 128, dropout = 0.2) {
        self$net <- nn_sequential(
            nn_linear(n_input, hidden),
            nn_batch_norm1d(hidden),
            swish_module(),
            nn_dropout(dropout),
            nn_linear(hidden, hidden),
            nn_batch_norm1d(hidden),
            swish_module(),
            nn_dropout(dropout),
            nn_linear(hidden, hidden),
            nn_batch_norm1d(hidden),
            swish_module(),
            nn_dropout(dropout),
            nn_linear(hidden, 64),
            swish_module(),
            nn_dropout(dropout * 0.5),
            nn_linear(64, 32),
            swish_module(),
            nn_linear(32, 1)
        )
    },
    forward = function(x) self$net(x)
)

cat("  ✓ DAG-Net + Flat NN 架构定义完成\n")

# ==============================================================================
# 4. 训练引擎
# ==============================================================================
cat("\n步骤 4/7: 训练引擎...\n")

# ---- Datasets ----
sdm_dataset <- dataset(
    "SDMDataset",
    initialize = function(x_topo, x_clim, x_surf, y) {
        self$x_topo <- torch_tensor(as.matrix(x_topo), dtype = torch_float())
        self$x_clim <- torch_tensor(as.matrix(x_clim), dtype = torch_float())
        self$x_surf <- torch_tensor(as.matrix(x_surf), dtype = torch_float())
        self$y <- torch_tensor(y, dtype = torch_float())$unsqueeze(2)
    },
    .getitem = function(i) {
        list(
            x_topo = self$x_topo[i, ], x_clim = self$x_clim[i, ],
            x_surf = self$x_surf[i, ], y = self$y[i, ]
        )
    },
    .length = function() self$y$size(1)
)

flat_dataset <- dataset(
    "FlatDataset",
    initialize = function(x, y) {
        self$x <- torch_tensor(as.matrix(x), dtype = torch_float())
        self$y <- torch_tensor(y, dtype = torch_float())$unsqueeze(2)
    },
    .getitem = function(i) list(x = self$x[i, ], y = self$y[i, ]),
    .length = function() self$y$size(1)
)

# ---- DAG-Net 训练函数 ----
train_dagnet <- function(model, train_dl, X_val_list, y_val_vec,
                         epochs = 300, lr = 1e-3, weight_decay = 1e-4,
                         patience = 30, pw = 1, verbose = TRUE) {
    optimizer <- optim_adamw(model$parameters, lr = lr, weight_decay = weight_decay)
    scheduler <- lr_cosine_annealing(optimizer, T_max = epochs, eta_min = 1e-5)

    pw_tensor <- torch_tensor(pw, dtype = torch_float())
    loss_fn <- function(pred, target) {
        nn_bce_with_logits_loss(pos_weight = pw_tensor)(pred, target)
    }

    best_val_auc <- 0
    best_state <- NULL
    epochs_no_improve <- 0
    history <- data.frame()

    for (epoch in seq_len(epochs)) {
        model$train()
        train_loss <- 0
        n_batch <- 0
        coro::loop(for (batch in train_dl) {
            optimizer$zero_grad()
            logits <- model(batch$x_topo, batch$x_clim, batch$x_surf)
            loss <- loss_fn(logits, batch$y)
            loss$backward()
            nn_utils_clip_grad_norm_(model$parameters, max_norm = 1.0)
            optimizer$step()
            train_loss <- train_loss + loss$item()
            n_batch <- n_batch + 1
        })
        scheduler$step()

        model$eval()
        with_no_grad({
            val_logits <- model(X_val_list$topo, X_val_list$clim, X_val_list$surf)
            val_pred <- as.numeric(torch_sigmoid(val_logits)$squeeze()$cpu())
        })
        val_auc <- tryCatch(
            as.numeric(pROC::auc(pROC::roc(y_val_vec, val_pred, quiet = TRUE))),
            error = function(e) 0
        )
        history <- rbind(history, data.frame(
            epoch = epoch, train_loss = train_loss / n_batch,
            val_auc = val_auc, lr = optimizer$param_groups[[1]]$lr
        ))
        if (verbose && epoch %% 10 == 0) {
            cat(sprintf(
                "  Epoch %3d: loss=%.4f val_AUC=%.4f lr=%.6f\n",
                epoch, train_loss / n_batch, val_auc, optimizer$param_groups[[1]]$lr
            ))
        }
        if (val_auc > best_val_auc + 1e-4) {
            best_val_auc <- val_auc
            best_state <- lapply(model$state_dict(), function(p) p$clone())
            epochs_no_improve <- 0
        } else {
            epochs_no_improve <- epochs_no_improve + 1
        }
        if (epochs_no_improve >= patience) {
            if (verbose) {
                cat(sprintf(
                    "  Early stop @ epoch %d (best val_AUC=%.4f)\n",
                    epoch, best_val_auc
                ))
            }
            break
        }
    }
    if (!is.null(best_state)) model$load_state_dict(best_state)
    list(model = model, history = history, best_val_auc = best_val_auc)
}

# ---- Flat NN 训练函数 ----
train_flatnn <- function(model, train_dl, X_val_tensor, y_val_vec,
                         epochs = 300, lr = 1e-3, weight_decay = 1e-4,
                         patience = 30, pw = 1, verbose = TRUE) {
    optimizer <- optim_adamw(model$parameters, lr = lr, weight_decay = weight_decay)
    scheduler <- lr_cosine_annealing(optimizer, T_max = epochs, eta_min = 1e-5)

    pw_tensor <- torch_tensor(pw, dtype = torch_float())
    loss_fn <- function(pred, target) {
        nn_bce_with_logits_loss(pos_weight = pw_tensor)(pred, target)
    }

    best_val_auc <- 0
    best_state <- NULL
    epochs_no_improve <- 0

    for (epoch in seq_len(epochs)) {
        model$train()
        train_loss <- 0
        n_batch <- 0
        coro::loop(for (batch in train_dl) {
            optimizer$zero_grad()
            logits <- model(batch$x)
            loss <- loss_fn(logits, batch$y)
            loss$backward()
            nn_utils_clip_grad_norm_(model$parameters, max_norm = 1.0)
            optimizer$step()
            train_loss <- train_loss + loss$item()
            n_batch <- n_batch + 1
        })
        scheduler$step()

        model$eval()
        with_no_grad({
            val_pred <- as.numeric(torch_sigmoid(model(X_val_tensor))$squeeze()$cpu())
        })
        val_auc <- tryCatch(
            as.numeric(pROC::auc(pROC::roc(y_val_vec, val_pred, quiet = TRUE))),
            error = function(e) 0
        )
        if (verbose && epoch %% 10 == 0) {
            cat(sprintf(
                "  Epoch %3d: loss=%.4f val_AUC=%.4f\n",
                epoch, train_loss / n_batch, val_auc
            ))
        }
        if (val_auc > best_val_auc + 1e-4) {
            best_val_auc <- val_auc
            best_state <- lapply(model$state_dict(), function(p) p$clone())
            epochs_no_improve <- 0
        } else {
            epochs_no_improve <- epochs_no_improve + 1
        }
        if (epochs_no_improve >= patience) {
            if (verbose) {
                cat(sprintf(
                    "  Early stop @ epoch %d (best=%.4f)\n",
                    epoch, best_val_auc
                ))
            }
            break
        }
    }
    if (!is.null(best_state)) model$load_state_dict(best_state)
    list(model = model, best_val_auc = best_val_auc)
}

# ---- 预测辅助 ----
predict_dagnet <- function(model, x_topo, x_clim, x_surf) {
    model$eval()
    with_no_grad({
        t_topo <- torch_tensor(as.matrix(x_topo), dtype = torch_float())
        t_clim <- torch_tensor(as.matrix(x_clim), dtype = torch_float())
        t_surf <- torch_tensor(as.matrix(x_surf), dtype = torch_float())
        as.numeric(torch_sigmoid(model(t_topo, t_clim, t_surf))$squeeze()$cpu())
    })
}

predict_flatnn <- function(model, x) {
    model$eval()
    with_no_grad({
        as.numeric(torch_sigmoid(
            model(torch_tensor(as.matrix(x), dtype = torch_float()))
        )$squeeze()$cpu())
    })
}

evaluate_model <- function(pred, y_true, model_name) {
    roc_obj <- pROC::roc(y_true, pred, quiet = TRUE)
    auc_val <- as.numeric(pROC::auc(roc_obj))
    coords_all <- pROC::coords(roc_obj,
        x = "all",
        ret = c("threshold", "sensitivity", "specificity")
    )
    tss_vals <- coords_all$sensitivity + coords_all$specificity - 1
    opt_idx <- which.max(tss_vals)
    data.frame(
        model = model_name, auc = auc_val, tss = tss_vals[opt_idx],
        sensitivity = coords_all$sensitivity[opt_idx],
        specificity = coords_all$specificity[opt_idx],
        stringsAsFactors = FALSE
    )
}

cat("  ✓ 训练引擎定义完成\n")

# ==============================================================================
# 5. 训练 DAG-Net 多次运行 (公平集成)
# ==============================================================================
cat("\n步骤 5/7: 训练神经网络模型 (5-run ensemble)...\n\n")

n_runs <- 5
run_seeds <- c(42, 55, 71, 89, 103)

# --- 数据加载器 ---
train_ds <- sdm_dataset(X_tr[, idx_topo], X_tr[, idx_clim], X_tr[, idx_surf], y_tr)
train_dl <- dataloader(train_ds, batch_size = 128, shuffle = TRUE, drop_last = TRUE)

X_val_list <- list(
    topo = torch_tensor(as.matrix(X_val[, idx_topo]), dtype = torch_float()),
    clim = torch_tensor(as.matrix(X_val[, idx_clim]), dtype = torch_float()),
    surf = torch_tensor(as.matrix(X_val[, idx_surf]), dtype = torch_float())
)

flat_train_ds <- flat_dataset(X_tr, y_tr)
flat_train_dl <- dataloader(flat_train_ds, batch_size = 128, shuffle = TRUE, drop_last = TRUE)
X_val_flat <- torch_tensor(as.matrix(X_val), dtype = torch_float())

# ---- 5a. DAG-Net (完整: 因果结构 + ATE初始化) ----
cat("  === DAG-Net (Full: Structure + ATE-Init) ===\n")
dagnet_preds <- matrix(0, nrow = nrow(X_test), ncol = n_runs)
dagnet_val_aucs <- numeric(n_runs)
dagnet_histories <- list()

for (run_i in seq_len(n_runs)) {
    cat(sprintf("\n  --- Run %d/%d (seed=%d) ---\n", run_i, n_runs, run_seeds[run_i]))
    torch_manual_seed(run_seeds[run_i])
    set.seed(run_seeds[run_i])

    model_i <- DAGNet(
        n_topo = n_topo, n_clim = n_clim, n_surf = n_surf,
        hidden = 128, dropout = 0.2,
        ate_topo = ate_weights[topo_vars],
        ate_clim = ate_weights[clim_vars],
        ate_surf = ate_weights[surf_vars]
    )
    if (run_i == 1) {
        n_params <- sum(sapply(model_i$parameters, function(p) p$numel()))
        cat(sprintf("  DAG-Net 参数量: %s\n", format(n_params, big.mark = ",")))
    }

    result_i <- train_dagnet(
        model_i, train_dl, X_val_list, y_val,
        epochs = 300, lr = 1e-3, weight_decay = 1e-4,
        patience = 30, pw = pos_weight
    )

    pred_i <- predict_dagnet(
        result_i$model,
        X_test[, idx_topo], X_test[, idx_clim], X_test[, idx_surf]
    )
    dagnet_preds[, run_i] <- pred_i
    dagnet_val_aucs[run_i] <- result_i$best_val_auc
    dagnet_histories[[run_i]] <- result_i$history

    test_auc_i <- as.numeric(pROC::auc(pROC::roc(y_test, pred_i, quiet = TRUE)))
    cat(sprintf("  Run %d: test_AUC=%.4f val_AUC=%.4f\n", run_i, test_auc_i, result_i$best_val_auc))
}

# DAG-Net ensemble (val-AUC加权)
ens_w <- dagnet_val_aucs / sum(dagnet_val_aucs)
pred_dagnet_ens <- as.numeric(dagnet_preds %*% ens_w)
pred_dagnet_best <- dagnet_preds[, which.max(dagnet_val_aucs)]

# ---- 5b. DAG-Net (消融: 无ATE初始化) ----
cat("\n\n  === DAG-Net (Ablation: No ATE-Init) ===\n")
noate_preds <- matrix(0, nrow = nrow(X_test), ncol = n_runs)
noate_val_aucs <- numeric(n_runs)

for (run_i in seq_len(n_runs)) {
    cat(sprintf("\n  --- Run %d/%d (seed=%d) ---\n", run_i, n_runs, run_seeds[run_i]))
    torch_manual_seed(run_seeds[run_i])
    set.seed(run_seeds[run_i])

    model_i <- DAGNet(
        n_topo = n_topo, n_clim = n_clim, n_surf = n_surf,
        hidden = 128, dropout = 0.2,
        ate_topo = NULL, ate_clim = NULL, ate_surf = NULL # 无ATE
    )

    result_i <- train_dagnet(
        model_i, train_dl, X_val_list, y_val,
        epochs = 300, lr = 1e-3, weight_decay = 1e-4,
        patience = 30, pw = pos_weight
    )

    pred_i <- predict_dagnet(
        result_i$model,
        X_test[, idx_topo], X_test[, idx_clim], X_test[, idx_surf]
    )
    noate_preds[, run_i] <- pred_i
    noate_val_aucs[run_i] <- result_i$best_val_auc

    test_auc_i <- as.numeric(pROC::auc(pROC::roc(y_test, pred_i, quiet = TRUE)))
    cat(sprintf("  Run %d: test_AUC=%.4f val_AUC=%.4f\n", run_i, test_auc_i, result_i$best_val_auc))
}

ens_w_noate <- noate_val_aucs / sum(noate_val_aucs)
pred_noate_ens <- as.numeric(noate_preds %*% ens_w_noate)

# ---- 5c. Flat NN (消融: 无因果结构) ----
cat("\n\n  === Flat NN (Ablation: No Causal Structure) ===\n")
flat_preds <- matrix(0, nrow = nrow(X_test), ncol = n_runs)
flat_val_aucs <- numeric(n_runs)

for (run_i in seq_len(n_runs)) {
    cat(sprintf("\n  --- Run %d/%d (seed=%d) ---\n", run_i, n_runs, run_seeds[run_i]))
    torch_manual_seed(run_seeds[run_i])
    set.seed(run_seeds[run_i])

    model_i <- FlatNN(n_input = n_vars, hidden = 128, dropout = 0.2)
    if (run_i == 1) {
        n_params_flat <- sum(sapply(model_i$parameters, function(p) p$numel()))
        cat(sprintf("  Flat NN 参数量: %s\n", format(n_params_flat, big.mark = ",")))
    }

    result_i <- train_flatnn(
        model_i, flat_train_dl, X_val_flat, y_val,
        epochs = 300, lr = 1e-3, weight_decay = 1e-4,
        patience = 30, pw = pos_weight
    )

    pred_i <- predict_flatnn(result_i$model, X_test)
    flat_preds[, run_i] <- pred_i
    flat_val_aucs[run_i] <- result_i$best_val_auc

    test_auc_i <- as.numeric(pROC::auc(pROC::roc(y_test, pred_i, quiet = TRUE)))
    cat(sprintf("  Run %d: test_AUC=%.4f val_AUC=%.4f\n", run_i, test_auc_i, result_i$best_val_auc))
}

ens_w_flat <- flat_val_aucs / sum(flat_val_aucs)
pred_flat_ens <- as.numeric(flat_preds %*% ens_w_flat)
pred_flat_best <- flat_preds[, which.max(flat_val_aucs)]

# ==============================================================================
# 6. 传统 SDM 基线
# ==============================================================================
cat("\n\n步骤 6/7: 传统SDM基线...\n\n")

results <- list()

# ---- NN 结果 (集成 + 单最优) ----
results[["DAGNet_Ens"]] <- evaluate_model(pred_dagnet_ens, y_test, "DAG-Net (Ensemble)")
results[["DAGNet_Best"]] <- evaluate_model(pred_dagnet_best, y_test, "DAG-Net (Best)")
results[["NoATE_Ens"]] <- evaluate_model(pred_noate_ens, y_test, "DAG-Net (No ATE)")
results[["FlatNN_Ens"]] <- evaluate_model(pred_flat_ens, y_test, "Flat NN (Ensemble)")
results[["FlatNN_Best"]] <- evaluate_model(pred_flat_best, y_test, "Flat NN (Best)")

cat(sprintf(
    "  DAG-Net Ensemble:     AUC=%.4f TSS=%.4f\n",
    results[["DAGNet_Ens"]]$auc, results[["DAGNet_Ens"]]$tss
))
cat(sprintf(
    "  DAG-Net (No ATE):     AUC=%.4f TSS=%.4f\n",
    results[["NoATE_Ens"]]$auc, results[["NoATE_Ens"]]$tss
))
cat(sprintf(
    "  Flat NN Ensemble:     AUC=%.4f TSS=%.4f\n\n",
    results[["FlatNN_Ens"]]$auc, results[["FlatNN_Ens"]]$tss
))

# ---- RF ----
cat("  [RF] Random Forest...\n")
tryCatch(
    {
        set.seed(42)
        rf_model <- ranger::ranger(
            presence ~ .,
            data = cbind(presence = as.factor(y_train), X_train),
            num.trees = 1000, probability = TRUE,
            min.node.size = 5, seed = 42
        )
        pred_rf <- predict(rf_model, data = X_test)$predictions[, "1"]
        results[["RF"]] <- evaluate_model(pred_rf, y_test, "Random Forest")
        cat(sprintf("  RF: AUC=%.4f TSS=%.4f\n", results[["RF"]]$auc, results[["RF"]]$tss))
    },
    error = function(e) cat("  RF FAILED:", e$message, "\n")
)

# ---- Maxent ----
cat("  [Maxent]...\n")
tryCatch(
    {
        mx_model <- maxnet::maxnet(
            p = y_train, data = X_train,
            maxnet.formula(p = y_train, data = X_train)
        )
        pred_mx <- as.numeric(predict(mx_model, X_test, type = "logistic"))
        results[["Maxent"]] <- evaluate_model(pred_mx, y_test, "Maxent")
        cat(sprintf("  Maxent: AUC=%.4f TSS=%.4f\n", results[["Maxent"]]$auc, results[["Maxent"]]$tss))
    },
    error = function(e) cat("  Maxent FAILED:", e$message, "\n")
)

# ---- BRT (Boosted Regression Trees) ----
cat("  [BRT] Gradient Boosted Trees...\n")
tryCatch(
    {
        set.seed(42)
        brt_data <- cbind(presence = y_train, X_train)
        brt_model <- gbm::gbm(
            presence ~ .,
            data = brt_data,
            distribution = "bernoulli",
            n.trees = 1000,
            interaction.depth = 5,
            shrinkage = 0.01,
            n.minobsinnode = 10,
            cv.folds = 5,
            verbose = FALSE
        )
        best_trees <- gbm::gbm.perf(brt_model, method = "cv", plot.it = FALSE)
        pred_brt <- predict(brt_model, newdata = X_test, n.trees = best_trees, type = "response")
        results[["BRT"]] <- evaluate_model(pred_brt, y_test, "BRT")
        cat(sprintf(
            "  BRT: AUC=%.4f TSS=%.4f (best_trees=%d)\n",
            results[["BRT"]]$auc, results[["BRT"]]$tss, best_trees
        ))
    },
    error = function(e) cat("  BRT FAILED:", e$message, "\n")
)

# ---- GAM ----
cat("  [GAM] Generalized Additive Model...\n")
tryCatch(
    {
        # 选前15个变量避免GAM过慢
        imp_vars <- names(sort(ate_weights, decreasing = TRUE))[1:min(15, n_vars)]
        gam_formula <- as.formula(paste(
            "presence ~",
            paste(paste0("s(", imp_vars, ", k=5)"), collapse = " + ")
        ))
        gam_data <- cbind(presence = y_train, X_train[, imp_vars, drop = FALSE])
        gam_model <- mgcv::gam(gam_formula, data = gam_data, family = binomial(), method = "REML")
        pred_gam <- as.numeric(predict(gam_model,
            newdata = X_test[, imp_vars, drop = FALSE], type = "response"
        ))
        results[["GAM"]] <- evaluate_model(pred_gam, y_test, "GAM")
        cat(sprintf("  GAM: AUC=%.4f TSS=%.4f\n", results[["GAM"]]$auc, results[["GAM"]]$tss))
    },
    error = function(e) cat("  GAM FAILED:", e$message, "\n")
)

# ==============================================================================
# 7. 结果汇总 + 消融分析 + 可视化
# ==============================================================================
cat("\n步骤 7/7: 结果汇总...\n\n")

# --- 汇总表 ---
comparison_df <- bind_rows(results) %>%
    mutate(type = case_when(
        grepl("DAG-Net", model) ~ "DAG-Net (ours)",
        grepl("Flat", model) ~ "NN Baseline",
        TRUE ~ "Traditional SDM"
    )) %>%
    arrange(desc(auc))

cat("  ╔═══════════════════════════════════════════════════════════════╗\n")
cat("  ║          Model Performance Comparison (Test Set)            ║\n")
cat("  ╠═══════════════════════════════════════════════════════════════╣\n")
for (i in seq_len(nrow(comparison_df))) {
    r <- comparison_df[i, ]
    marker <- if (grepl("DAG-Net", r$model)) " ★" else "  "
    cat(sprintf(
        "  ║%s %-32s AUC=%.4f TSS=%.4f ║\n",
        marker, r$model, r$auc, r$tss
    ))
}
cat("  ╚═══════════════════════════════════════════════════════════════╝\n\n")

write.csv(comparison_df, "output/17_dagnet/comparison_table.csv", row.names = FALSE)

# --- 消融分析 ---
cat("  ══════════════════════════════════════\n")
cat("  消融分析 (Ablation Study)\n")
cat("  ══════════════════════════════════════\n")

auc_dagnet <- results[["DAGNet_Ens"]]$auc
auc_noate <- results[["NoATE_Ens"]]$auc
auc_flat <- results[["FlatNN_Ens"]]$auc

cat(sprintf("  DAG-Net (Full):   AUC = %.4f\n", auc_dagnet))
cat(sprintf("  DAG-Net (No ATE): AUC = %.4f\n", auc_noate))
cat(sprintf("  Flat NN:          AUC = %.4f\n", auc_flat))
cat("  ──────────────────────────────────────\n")
cat(sprintf(
    "  因果结构贡献:     ΔAUC = %+.4f (DAG-Net No-ATE vs Flat NN)\n",
    auc_noate - auc_flat
))
cat(sprintf(
    "  ATE初始化贡献:    ΔAUC = %+.4f (DAG-Net Full vs No-ATE)\n",
    auc_dagnet - auc_noate
))
cat(sprintf(
    "  总因果贡献:       ΔAUC = %+.4f (DAG-Net Full vs Flat NN)\n",
    auc_dagnet - auc_flat
))

# 保存消融表
ablation_df <- data.frame(
    model = c("DAG-Net (Full)", "DAG-Net (No ATE-Init)", "Flat NN"),
    auc = c(auc_dagnet, auc_noate, auc_flat),
    delta_vs_flat = c(auc_dagnet - auc_flat, auc_noate - auc_flat, 0),
    component = c("Structure + ATE", "Structure only", "None (baseline)"),
    stringsAsFactors = FALSE
)
write.csv(ablation_df, "output/17_dagnet/ablation_results.csv", row.names = FALSE)

# --- 可视化 ---
color_map <- c(
    "DAG-Net (ours)" = "#E41A1C", "NN Baseline" = "#377EB8",
    "Traditional SDM" = "#636363"
)

p1 <- ggplot(comparison_df, aes(x = reorder(model, auc), y = auc, fill = type)) +
    geom_col(width = 0.7, alpha = 0.9) +
    geom_text(aes(label = sprintf("%.4f", auc)), hjust = -0.1, size = 3) +
    scale_fill_manual(values = color_map, name = "Model Type") +
    coord_flip(ylim = c(min(comparison_df$auc) - 0.02, 1.0)) +
    labs(
        title = "Test Set AUC Comparison",
        subtitle = "DAG-Net (Causal Structure + ATE Init) vs Baselines",
        x = "", y = "AUC"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

# 训练曲线 (best run)
best_run_idx <- which.max(dagnet_val_aucs)
hist_df <- dagnet_histories[[best_run_idx]]
p2 <- ggplot(hist_df, aes(x = epoch)) +
    geom_line(aes(y = val_auc), color = "#E41A1C", linewidth = 0.8) +
    geom_line(aes(y = train_loss), color = "#377EB8", linewidth = 0.5, linetype = "dashed") +
    labs(
        title = sprintf("DAG-Net Training (Best Run %d)", best_run_idx),
        x = "Epoch", y = "Value (red=val_AUC, blue=train_loss)"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

# 消融图
ablation_df$model <- factor(ablation_df$model,
    levels = c("Flat NN", "DAG-Net (No ATE-Init)", "DAG-Net (Full)")
)
p3 <- ggplot(ablation_df, aes(x = model, y = auc, fill = component)) +
    geom_col(width = 0.6, alpha = 0.9) +
    geom_text(aes(label = sprintf("%.4f", auc)), vjust = -0.5, size = 3.5) +
    scale_fill_manual(
        values = c(
            "None (baseline)" = "#636363",
            "Structure only" = "#FF7F00",
            "Structure + ATE" = "#E41A1C"
        ),
        name = "Components"
    ) +
    labs(
        title = "Ablation Study: Causal Components",
        subtitle = "Quantifying the contribution of each causal component",
        x = "", y = "AUC"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

combined <- (p1 | p3) / p2 +
    plot_annotation(
        title = "DAG-Net: Causally-Structured Neural Network for SDM",
        theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    )

ggsave("figures/17_dagnet/performance_comparison.png",
    plot = combined, width = 14, height = 10, dpi = 300, bg = "white"
)

# 保存训练历史
all_hist_df <- do.call(rbind, lapply(seq_along(dagnet_histories), function(i) {
    h <- dagnet_histories[[i]]
    h$run <- i
    h
}))
write.csv(all_hist_df, "output/17_dagnet/dagnet_training_history.csv", row.names = FALSE)

cat("  ✓ 图表已保存: figures/17_dagnet/performance_comparison.png\n\n")

cat(paste(rep("=", 70), collapse = ""), "\n")
cat("               DAG-Net 实验完成\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

best <- comparison_df[1, ]
cat(sprintf("  最优模型: %s (AUC=%.4f, TSS=%.4f)\n", best$model, best$auc, best$tss))
cat("\n  输出文件:\n")
cat("    output/17_dagnet/comparison_table.csv\n")
cat("    output/17_dagnet/ablation_results.csv\n")
cat("    output/17_dagnet/dagnet_training_history.csv\n")
cat("    figures/17_dagnet/performance_comparison.png\n\n")
cat("✓ 完成!\n")
