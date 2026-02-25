#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 17_dagnet_model.R
# 功能说明: DAG-Net — 因果结构化神经网络 SDM (torch 专业实现)
#
# 核心创新: 网络架构由因果DAG决定。环境变量按因果层级分批进入不同网络层，
#           信息沿因果方向自上而下流动。ATE估计值用于缩放输入权重。
#
# 架构 (torch nn_module):
#   Layer 0 (Causal Roots):     地形/拓扑 (6 vars)  → h₀
#   Layer 1 (Causal Mediators): [h₀; 水文气候 19 vars] → h₁ (+残差)
#   Layer 2 (Causal Responses): [h₁; 土地覆盖+土壤 22 vars] → h₂ (+残差)
#   Output:  h₂ → ŷ ∈ [0,1]
#
# 对照: DAG-Net Ensemble vs Flat-NN vs Maxent vs RF vs CPF
#
# 依赖: torch, ranger, maxnet, pROC, caret, tidyverse, ggplot2, patchwork
# ==============================================================================

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

# ---- 安装并加载依赖 ----
packages <- c(
    "tidyverse", "maxnet", "randomForest", "ranger", "pROC", "caret",
    "ggplot2", "viridis", "patchwork", "sysfonts", "showtext"
)
for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
}

# 安装 torch (核心深度学习框架)
if (!require(torch, quietly = TRUE)) {
    cat("  正在安装 torch...\n")
    install.packages("torch")
    library(torch)
    if (!torch_is_installed()) {
        cat("  正在安装 libtorch 后端...\n")
        torch::install_torch()
    }
}
cat(sprintf("  torch version: %s\n", as.character(packageVersion("torch"))))

try(
    {
        sysfonts::font_add(
            family = "Arial",
            regular = "C:/Windows/Fonts/arial.ttf",
            bold = "C:/Windows/Fonts/arialbd.ttf"
        )
        showtext::showtext_opts(dpi = 300)
        showtext::showtext_auto(enable = TRUE)
    },
    silent = TRUE
)

dir.create("output/17_dagnet", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/17_dagnet", showWarnings = FALSE, recursive = TRUE)
try(source("scripts/visualization/viz_utils.R"), silent = TRUE)

cat("\n", paste(rep("=", 70), collapse = ""), "\n")
cat("     DAG-Net: Causal Architecture Neural Network for SDM\n")
cat("     Powered by R torch\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

# ==============================================================================
# 1. 数据加载与变量分层
# ==============================================================================
cat("步骤 1/8: 数据加载与因果层级分组...\n")

dat <- read.csv("output/04_collinearity/collinearity_removed.csv", stringsAsFactors = FALSE)
groups <- read.csv("output/14_causal/variable_groups_validated.csv", stringsAsFactors = FALSE)
ate_df <- read.csv("output/14_causal/ate_all_variables.csv", stringsAsFactors = FALSE)

topo_vars <- intersect(groups$var[groups$group == "G1_TopoSlopeFlow"], colnames(dat))
clim_vars <- intersect(groups$var[groups$group == "G2_Hydroclim_wavg"], colnames(dat))
lc_vars <- intersect(groups$var[groups$group == "G3_Landcover_wavg"], colnames(dat))
soil_vars <- intersect(groups$var[groups$group == "G4_Soil_wavg"], colnames(dat))
surf_vars <- c(lc_vars, soil_vars)
all_vars <- c(topo_vars, clim_vars, surf_vars)

cat(sprintf(
    "  Layer 0 (Causal Roots):     %d vars [%s]\n",
    length(topo_vars), paste(topo_vars, collapse = ", ")
))
cat(sprintf("  Layer 1 (Causal Mediators): %d vars\n", length(clim_vars)))
cat(sprintf(
    "  Layer 2 (Causal Responses): %d vars (%d LC + %d Soil)\n",
    length(surf_vars), length(lc_vars), length(soil_vars)
))
cat(sprintf("  Total: %d vars\n\n", length(all_vars)))

# ATE权重 → 归一化到 [0.1, 1]
ate_weights <- setNames(rep(0.1, length(all_vars)), all_vars)
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
cat("\n步骤 2/8: 数据准备...\n")

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

cat(sprintf("  训练: %d | 验证: %d | 测试: %d\n", nrow(X_tr), nrow(X_val), nrow(X_test)))

# ==============================================================================
# 3. DAG-Net 定义 (torch nn_module)
# ==============================================================================
cat("\n步骤 3/8: 定义DAG-Net (torch)...\n")

# Swish 激活函数模块
swish_module <- nn_module(
    "Swish",
    forward = function(x) x * torch_sigmoid(x)
)

# ---- DAG-Net 核心模型 ----
DAGNet <- nn_module(
    "DAGNet",
    initialize = function(n_topo, n_clim, n_surf, hidden = 128,
                          dropout = 0.2, ate_topo = NULL, ate_clim = NULL, ate_surf = NULL) {
        # 可学习的ATE缩放因子 (因果先验嵌入参数)
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

        # Layer 0: Causal Roots (地形/拓扑)
        self$layer0 <- nn_sequential(
            nn_linear(n_topo, hidden),
            nn_batch_norm1d(hidden),
            swish_module(),
            nn_dropout(dropout)
        )

        # Layer 1: Causal Mediators (水文气候, 接收h0)
        self$layer1 <- nn_sequential(
            nn_linear(hidden + n_clim, hidden),
            nn_batch_norm1d(hidden),
            swish_module(),
            nn_dropout(dropout)
        )

        # Layer 2: Causal Responses (土地覆盖+土壤, 接收h1)
        self$layer2 <- nn_sequential(
            nn_linear(hidden + n_surf, hidden),
            nn_batch_norm1d(hidden),
            swish_module(),
            nn_dropout(dropout)
        )

        # Output head (较深)
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
        # ATE缩放输入 (因果先验)
        x_topo <- x_topo * self$ate_scale_topo$unsqueeze(1)
        x_clim <- x_clim * self$ate_scale_clim$unsqueeze(1)
        x_surf <- x_surf * self$ate_scale_surf$unsqueeze(1)

        # Layer 0: 因果源头
        h0 <- self$layer0(x_topo)

        # Layer 1: 因果中介 + 残差连接
        h1 <- self$layer1(torch_cat(list(h0, x_clim), dim = 2)) + h0

        # Layer 2: 因果终端 + 残差连接
        h2 <- self$layer2(torch_cat(list(h1, x_surf), dim = 2)) + h1

        # Output
        self$head(h2)
    }
)

# ---- Flat NN (对照: 同参数量, 无因果结构) ----
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
    forward = function(x) {
        self$net(x)
    }
)

# ==============================================================================
# 4. 训练引擎
# ==============================================================================
cat("\n步骤 4/8: 训练引擎...\n")

# 创建torch数据集
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
    .length = function() {
        self$y$size(1)
    }
)

# 全变量数据集 (给Flat NN)
flat_dataset <- dataset(
    "FlatDataset",
    initialize = function(x, y) {
        self$x <- torch_tensor(as.matrix(x), dtype = torch_float())
        self$y <- torch_tensor(y, dtype = torch_float())$unsqueeze(2)
    },
    .getitem = function(i) {
        list(x = self$x[i, ], y = self$y[i, ])
    },
    .length = function() {
        self$y$size(1)
    }
)

# ---- 训练函数 (DAG-Net) ----
train_dagnet_torch <- function(model, train_dl, val_dl, X_val_list, y_val_vec,
                               epochs = 300, lr = 1e-3, weight_decay = 1e-4,
                               patience = 30, pos_weight_val = 1, verbose = TRUE) {
    optimizer <- optim_adamw(model$parameters, lr = lr, weight_decay = weight_decay)
    scheduler <- lr_cosine_annealing(optimizer, T_max = epochs, eta_min = 1e-5)

    pw_tensor <- torch_tensor(pos_weight_val, dtype = torch_float())
    loss_fn <- function(pred, target) {
        nn_bce_with_logits_loss(pos_weight = pw_tensor)(pred, target)
    }

    best_val_auc <- 0
    best_state <- NULL
    epochs_no_improve <- 0
    history <- data.frame()

    for (epoch in seq_len(epochs)) {
        # --- Training ---
        model$train()
        train_loss <- 0
        n_batch <- 0

        coro::loop(for (batch in train_dl) {
            optimizer$zero_grad()
            logits <- model(batch$x_topo, batch$x_clim, batch$x_surf)
            loss <- loss_fn(logits, batch$y)
            loss$backward()

            # 梯度裁剪
            nn_utils_clip_grad_norm_(model$parameters, max_norm = 1.0)

            optimizer$step()
            train_loss <- train_loss + loss$item()
            n_batch <- n_batch + 1
        })
        scheduler$step()

        # --- Validation ---
        model$eval()
        with_no_grad({
            val_logits <- model(X_val_list$topo, X_val_list$clim, X_val_list$surf)
            val_pred <- torch_sigmoid(val_logits)$squeeze()$to(device = "cpu")
            val_pred_r <- as.numeric(val_pred)
        })

        val_auc <- tryCatch(
            as.numeric(pROC::auc(pROC::roc(y_val_vec, val_pred_r, quiet = TRUE))),
            error = function(e) 0
        )

        history <- rbind(history, data.frame(
            epoch = epoch, train_loss = train_loss / n_batch, val_auc = val_auc,
            lr = optimizer$param_groups[[1]]$lr
        ))

        if (verbose && epoch %% 10 == 0) {
            cat(sprintf(
                "  Epoch %3d: loss=%.4f val_AUC=%.4f lr=%.6f\n",
                epoch, train_loss / n_batch, val_auc, optimizer$param_groups[[1]]$lr
            ))
        }

        # Early stopping (基于val AUC)
        if (val_auc > best_val_auc + 1e-4) {
            best_val_auc <- val_auc
            best_state <- lapply(model$state_dict(), function(p) p$clone())
            epochs_no_improve <- 0
        } else {
            epochs_no_improve <- epochs_no_improve + 1
        }
        if (epochs_no_improve >= patience) {
            if (verbose) cat(sprintf("  Early stop @ epoch %d (best val_AUC=%.4f)\n", epoch, best_val_auc))
            break
        }
    }

    # 恢复最优参数
    if (!is.null(best_state)) model$load_state_dict(best_state)
    list(model = model, history = history, best_val_auc = best_val_auc)
}

# ---- 训练函数 (Flat NN) ----
train_flatnn_torch <- function(model, train_dl, val_dl, X_val_tensor, y_val_vec,
                               epochs = 300, lr = 1e-3, weight_decay = 1e-4,
                               patience = 30, pos_weight_val = 1, verbose = TRUE) {
    optimizer <- optim_adamw(model$parameters, lr = lr, weight_decay = weight_decay)
    scheduler <- lr_cosine_annealing(optimizer, T_max = epochs, eta_min = 1e-5)

    pw_tensor <- torch_tensor(pos_weight_val, dtype = torch_float())
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
            val_logits <- model(X_val_tensor)
            val_pred_r <- as.numeric(torch_sigmoid(val_logits)$squeeze()$cpu())
        })
        val_auc <- tryCatch(
            as.numeric(pROC::auc(pROC::roc(y_val_vec, val_pred_r, quiet = TRUE))),
            error = function(e) 0
        )
        if (verbose && epoch %% 10 == 0) {
            cat(sprintf("  Epoch %3d: loss=%.4f val_AUC=%.4f\n", epoch, train_loss / n_batch, val_auc))
        }
        if (val_auc > best_val_auc + 1e-4) {
            best_val_auc <- val_auc
            best_state <- lapply(model$state_dict(), function(p) p$clone())
            epochs_no_improve <- 0
        } else {
            epochs_no_improve <- epochs_no_improve + 1
        }
        if (epochs_no_improve >= patience) {
            if (verbose) cat(sprintf("  Early stop @ epoch %d (best=%.4f)\n", epoch, best_val_auc))
            break
        }
    }
    if (!is.null(best_state)) model$load_state_dict(best_state)
    list(model = model, best_val_auc = best_val_auc)
}

# ---- 预测辅助函数 ----
predict_dagnet <- function(model, x_topo, x_clim, x_surf) {
    model$eval()
    with_no_grad({
        t_topo <- torch_tensor(as.matrix(x_topo), dtype = torch_float())
        t_clim <- torch_tensor(as.matrix(x_clim), dtype = torch_float())
        t_surf <- torch_tensor(as.matrix(x_surf), dtype = torch_float())
        logits <- model(t_topo, t_clim, t_surf)
        as.numeric(torch_sigmoid(logits)$squeeze()$cpu())
    })
}

predict_flatnn <- function(model, x) {
    model$eval()
    with_no_grad({
        t_x <- torch_tensor(as.matrix(x), dtype = torch_float())
        as.numeric(torch_sigmoid(model(t_x))$squeeze()$cpu())
    })
}

# ---- 评估函数 ----
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

# ==============================================================================
# 5. 训练 DAG-Net Ensemble (5次独立训练)
# ==============================================================================
cat("\n步骤 5/8: 训练 DAG-Net Ensemble...\n\n")

pos_weight <- sum(y_tr == 0) / max(sum(y_tr == 1), 1)

# 创建数据加载器
train_ds <- sdm_dataset(X_tr[, idx_topo], X_tr[, idx_clim], X_tr[, idx_surf], y_tr)
train_dl <- dataloader(train_ds, batch_size = 64, shuffle = TRUE, drop_last = TRUE)
val_ds <- sdm_dataset(X_val[, idx_topo], X_val[, idx_clim], X_val[, idx_surf], y_val)
val_dl <- dataloader(val_ds, batch_size = nrow(X_val), shuffle = FALSE)

# 验证集tensor (用于快速评估)
X_val_list <- list(
    topo = torch_tensor(as.matrix(X_val[, idx_topo]), dtype = torch_float()),
    clim = torch_tensor(as.matrix(X_val[, idx_clim]), dtype = torch_float()),
    surf = torch_tensor(as.matrix(X_val[, idx_surf]), dtype = torch_float())
)

n_ensemble <- 5
dagnet_preds <- matrix(0, nrow = nrow(X_test), ncol = n_ensemble)
dagnet_val_aucs <- numeric(n_ensemble)
all_histories <- list()

for (run_i in seq_len(n_ensemble)) {
    cat(sprintf("  === DAG-Net Run %d/%d ===\n", run_i, n_ensemble))
    torch_manual_seed(42 + run_i * 13)
    set.seed(42 + run_i * 13)

    model_i <- DAGNet(
        n_topo = length(idx_topo), n_clim = length(idx_clim), n_surf = length(idx_surf),
        hidden = 128, dropout = 0.2,
        ate_topo = ate_weights[topo_vars],
        ate_clim = ate_weights[clim_vars],
        ate_surf = ate_weights[surf_vars]
    )

    result_i <- train_dagnet_torch(
        model_i, train_dl, val_dl, X_val_list, y_val,
        epochs = 300, lr = 1e-3, weight_decay = 1e-4,
        patience = 35, pos_weight_val = pos_weight
    )

    pred_i <- predict_dagnet(
        result_i$model,
        X_test[, idx_topo], X_test[, idx_clim], X_test[, idx_surf]
    )
    dagnet_preds[, run_i] <- pred_i
    dagnet_val_aucs[run_i] <- result_i$best_val_auc
    all_histories[[run_i]] <- result_i$history

    test_auc_i <- as.numeric(pROC::auc(pROC::roc(y_test, pred_i, quiet = TRUE)))
    cat(sprintf("  Run %d: test_AUC=%.4f val_AUC=%.4f\n\n", run_i, test_auc_i, result_i$best_val_auc))
}

# Ensemble: val-AUC加权平均
ens_weights <- dagnet_val_aucs / sum(dagnet_val_aucs)
pred_dagnet_ens <- as.numeric(dagnet_preds %*% ens_weights)

# Best single
best_run <- which.max(dagnet_val_aucs)
pred_dagnet_best <- dagnet_preds[, best_run]

cat("  Ensemble weights:", sprintf("%.3f", ens_weights), "\n")

results <- list()
results[["DAGNet_Ensemble"]] <- evaluate_model(pred_dagnet_ens, y_test, "DAG-Net Ensemble")
results[["DAGNet_Best"]] <- evaluate_model(pred_dagnet_best, y_test, "DAG-Net (Best Single)")
cat(sprintf(
    "  DAG-Net Ensemble: AUC=%.4f TSS=%.4f\n",
    results[["DAGNet_Ensemble"]]$auc, results[["DAGNet_Ensemble"]]$tss
))
cat(sprintf(
    "  DAG-Net Best:     AUC=%.4f TSS=%.4f\n",
    results[["DAGNet_Best"]]$auc, results[["DAGNet_Best"]]$tss
))

# ==============================================================================
# 6. Flat NN 对照
# ==============================================================================
cat("\n步骤 6/8: Flat NN + Baselines...\n\n")

cat("  [Flat NN] Training...\n")
flat_train_ds <- flat_dataset(X_tr, y_tr)
flat_train_dl <- dataloader(flat_train_ds, batch_size = 64, shuffle = TRUE, drop_last = TRUE)
X_val_flat <- torch_tensor(as.matrix(X_val), dtype = torch_float())

torch_manual_seed(42)
flat_model <- FlatNN(n_input = length(all_vars), hidden = 128, dropout = 0.2)
flat_result <- train_flatnn_torch(
    flat_model, flat_train_dl, NULL, X_val_flat, y_val,
    epochs = 300, lr = 1e-3, weight_decay = 1e-4,
    patience = 35, pos_weight_val = pos_weight
)
pred_flat <- predict_flatnn(flat_result$model, X_test)
results[["FlatNN"]] <- evaluate_model(pred_flat, y_test, "Flat NN")
cat(sprintf("  Flat NN: AUC=%.4f TSS=%.4f\n\n", results[["FlatNN"]]$auc, results[["FlatNN"]]$tss))

# ---- Causal Prior Forest (CPF): 因果加权随机森林 ----
cat("  [CPF] Causal Prior Forest...\n")
edges_df <- tryCatch(
    read.csv("output/14_causal/edges_summary.csv", stringsAsFactors = FALSE),
    error = function(e) data.frame(from = character(), to = character(), strength = numeric())
)
node_outdeg <- edges_df %>%
    filter(strength >= 0.55) %>%
    group_by(from) %>%
    summarise(out_degree = n(), .groups = "drop")

cpf_weights <- sapply(all_vars, function(v) {
    od <- node_outdeg$out_degree[node_outdeg$from == v]
    if (length(od) == 0) od <- 0
    ate_row <- ate_df[ate_df$variable == v, ]
    ate_sig <- if (nrow(ate_row) > 0 && ate_row$p_value[1] < 0.05) abs(ate_row$coef[1]) * 10 else 0
    max(od + ate_sig, 0.5)
})
cpf_weights <- 0.1 + 0.9 * (cpf_weights - min(cpf_weights)) /
    (max(cpf_weights) - min(cpf_weights) + 1e-8)

set.seed(42)
cpf_model <- ranger::ranger(
    presence ~ .,
    data = cbind(presence = as.factor(y_train), X_train),
    num.trees = 1000, probability = TRUE, importance = "permutation",
    split.select.weights = cpf_weights, min.node.size = 5, seed = 42
)
pred_cpf <- predict(cpf_model, data = X_test)$predictions[, "1"]
results[["CPF"]] <- evaluate_model(pred_cpf, y_test, "Causal Prior Forest")
cat(sprintf("  CPF: AUC=%.4f TSS=%.4f\n\n", results[["CPF"]]$auc, results[["CPF"]]$tss))

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

# ---- Random Forest ----
cat("  [RF]...\n")
tryCatch(
    {
        rf_model <- ranger::ranger(presence ~ .,
            data = cbind(presence = as.factor(y_train), X_train),
            num.trees = 1000, probability = TRUE, min.node.size = 5, seed = 42
        )
        pred_rf <- predict(rf_model, data = X_test)$predictions[, "1"]
        results[["RF"]] <- evaluate_model(pred_rf, y_test, "Random Forest")
        cat(sprintf("  RF: AUC=%.4f TSS=%.4f\n", results[["RF"]]$auc, results[["RF"]]$tss))
    },
    error = function(e) cat("  RF FAILED:", e$message, "\n")
)

# ==============================================================================
# 7. CAST Super-Ensemble: DAG-Net + CPF
# ==============================================================================
cat("\n步骤 7/8: CAST Super-Ensemble...\n")

# 将DAG-Net Ensemble和CPF按0.5:0.5混合
pred_cast <- 0.5 * pred_dagnet_ens + 0.5 * pred_cpf
results[["CAST_Super"]] <- evaluate_model(pred_cast, y_test, "CAST (DAG-Net + CPF)")
cat(sprintf(
    "  CAST Super: AUC=%.4f TSS=%.4f\n",
    results[["CAST_Super"]]$auc, results[["CAST_Super"]]$tss
))

# ==============================================================================
# 8. 结果汇总 + 可视化
# ==============================================================================
cat("\n步骤 8/8: 结果汇总...\n\n")

comparison_df <- bind_rows(results) %>%
    mutate(type = case_when(
        grepl("CAST|DAG|CPF", model) ~ "CAST (ours)",
        grepl("Flat", model) ~ "NN Baseline",
        TRUE ~ "Traditional SDM"
    )) %>%
    arrange(desc(auc))

cat("  ╔══════════════════════════════════════════════════════════════╗\n")
cat("  ║          Model Performance Comparison (Test Set)           ║\n")
cat("  ╠══════════════════════════════════════════════════════════════╣\n")
for (i in seq_len(nrow(comparison_df))) {
    r <- comparison_df[i, ]
    marker <- if (grepl("CAST", r$model)) " ★" else "  "
    cat(sprintf(
        "  ║%s %-30s AUC=%.4f TSS=%.4f ║\n",
        marker, r$model, r$auc, r$tss
    ))
}
cat("  ╚══════════════════════════════════════════════════════════════╝\n\n")

write.csv(comparison_df, "output/17_dagnet/comparison_table.csv", row.names = FALSE)

# 消融: DAGNet vs Flat
cat("  消融分析:\n")
if (!is.null(results[["DAGNet_Ensemble"]]) && !is.null(results[["FlatNN"]])) {
    cat(sprintf(
        "    DAG结构贡献:  ΔAUC = %+.4f (DAG-Net Ensemble vs Flat NN)\n",
        results[["DAGNet_Ensemble"]]$auc - results[["FlatNN"]]$auc
    ))
}
if (!is.null(results[["CPF"]]) && !is.null(results[["RF"]])) {
    cat(sprintf(
        "    因果权重贡献: ΔAUC = %+.4f (CPF vs RF)\n",
        results[["CPF"]]$auc - results[["RF"]]$auc
    ))
}

# ---- 可视化 ----
color_map <- c(
    "CAST (ours)" = "#E41A1C", "NN Baseline" = "#377EB8",
    "Traditional SDM" = "#636363"
)

p1 <- ggplot(comparison_df, aes(x = reorder(model, auc), y = auc, fill = type)) +
    geom_col(width = 0.7, alpha = 0.9) +
    geom_text(aes(label = sprintf("%.4f", auc)), hjust = -0.1, size = 3) +
    scale_fill_manual(values = color_map, name = "Model Type") +
    coord_flip() +
    labs(
        title = "Test Set AUC Comparison",
        subtitle = "CAST DAG-Net vs Traditional SDMs",
        x = "", y = "AUC"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"), legend.position = "bottom")

# 训练曲线 (best run)
hist_df <- all_histories[[best_run]]
p2 <- ggplot(hist_df, aes(x = epoch)) +
    geom_line(aes(y = val_auc), color = "#E41A1C", linewidth = 0.8) +
    geom_line(aes(y = train_loss), color = "#377EB8", linewidth = 0.5, linetype = "dashed") +
    labs(
        title = sprintf("DAG-Net Training (Run %d)", best_run),
        x = "Epoch", y = "Value (red=val_AUC, blue=train_loss)"
    ) +
    theme_minimal(base_size = 11) +
    theme(plot.title = element_text(face = "bold"))

combined <- p1 / p2 +
    plot_annotation(
        title = "CAST DAG-Net: Performance Analysis",
        theme = theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5))
    )

ggsave("figures/17_dagnet/performance_comparison.png",
    plot = combined, width = 12, height = 10, dpi = 300, bg = "white"
)

cat("  ✓ 图表已保存: figures/17_dagnet/performance_comparison.png\n\n")

# ==============================================================================
# 最终日志
# ==============================================================================
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("               DAG-Net 实验完成 (torch版)\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

best <- comparison_df[1, ]
cat(sprintf("  最优模型: %s (AUC=%.4f, TSS=%.4f)\n", best$model, best$auc, best$tss))
cat("\n  输出文件:\n")
cat("    output/17_dagnet/comparison_table.csv\n")
cat("    figures/17_dagnet/performance_comparison.png\n\n")
cat("✓ 完成!\n")
