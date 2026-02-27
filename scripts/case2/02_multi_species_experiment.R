################################################################################
# 02_multi_species_experiment.R
# CAST v3 Multi-Region Multi-Species Experiment
# disdat benchmark — 6 regions, 226 species
#
# CAST: Causal Analysis for Species distribution modelling Toolkit
#
# Pipeline per species (与01_cast_pipeline.R完全一致):
#   Step 1: VIF filtering (threshold=10)
#   Step 2: DAG learning (HC bootstrap R=100, strength≥0.7, direction≥0.6)
#   Step 3: ATE estimation (DML 2-fold cross-fitting)
#   Step 4: Adaptive CAST Screening v2 (用于解释，不限制建模)
#   Step 5: Causal Role Grouping
#   Step 6: Unified comparison table:
#       CAST (全变量ATE加权 + DAG交互特征, CI-MLP架构)
#       MLP  (全变量, 标准MLP)
#       RF / Maxent / BRT / GAM (全变量)
#
# Prerequisite: Run 00_data_preparation.R first
################################################################################

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

# ---- Dependencies ----
pkgs <- c("tidyverse", "bnlearn", "pROC", "caret", "ranger", "maxnet", "gbm", "mgcv", "torch")
for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
}
if (!torch_is_installed()) torch::install_torch()

# ---- Configuration ----
# Which regions to run (default: all 6)
REGIONS <- c("AWT", "CAN", "NSW", "NZ", "SA", "SWI")
n_runs <- 3
seeds <- c(42, 71, 103)

# ==============================================================================
# Utility Functions
# ==============================================================================
calc_vif <- function(data) {
    sapply(names(data), function(v) {
        others <- setdiff(names(data), v)
        if (length(others) == 0) {
            return(1)
        }
        r2 <- tryCatch(summary(lm(reformulate(others, v), data = data))$r.squared,
            error = function(e) 0
        )
        1 / (1 - r2)
    })
}

normalize01 <- function(x) {
    r <- range(x, na.rm = TRUE)
    if (r[2] - r[1] < 1e-10) {
        return(rep(0.5, length(x)))
    }
    (x - r[1]) / (r[2] - r[1])
}

evaluate_model <- function(pred, obs) {
    pred <- pmin(pmax(pred, 1e-7), 1 - 1e-7)
    auc_val <- tryCatch(as.numeric(pROC::auc(pROC::roc(obs, pred, quiet = TRUE))),
        error = function(e) NA_real_
    )
    tss_val <- tryCatch(
        {
            roc_obj <- pROC::roc(obs, pred, quiet = TRUE)
            coords <- pROC::coords(roc_obj, "best", ret = c("sensitivity", "specificity"))
            as.numeric(coords$sensitivity + coords$specificity - 1)
        },
        error = function(e) NA_real_
    )
    c(auc = auc_val, tss = tss_val)
}

# ==============================================================================
# DML ATE (Chernozhukov et al., 2018)
# ==============================================================================
dml_ate <- function(Y, T_var, W, K = 2, num_trees = 300) {
    n <- length(Y)
    folds <- sample(rep(1:K, length.out = n))
    y_res <- numeric(n)
    t_res <- numeric(n)
    for (k in 1:K) {
        train_idx <- which(folds != k)
        test_idx <- which(folds == k)
        W_train <- W[train_idx, , drop = FALSE]
        W_test <- W[test_idx, , drop = FALSE]
        rf_y <- ranger::ranger(y ~ .,
            data = cbind(y = Y[train_idx], W_train),
            num.trees = num_trees, verbose = FALSE
        )
        y_res[test_idx] <- Y[test_idx] - predict(rf_y, data = W_test)$predictions
        rf_t <- ranger::ranger(y ~ .,
            data = cbind(y = T_var[train_idx], W_train),
            num.trees = num_trees, verbose = FALSE
        )
        t_res[test_idx] <- T_var[test_idx] - predict(rf_t, data = W_test)$predictions
    }
    ate <- sum(t_res * y_res) / sum(t_res^2)
    residuals <- y_res - ate * t_res
    se <- sqrt(mean(residuals^2 * t_res^2) / (mean(t_res^2)^2) / n)
    p_value <- 2 * pnorm(-abs(ate / se))
    list(ate = ate, se = se, p_value = p_value, significant = p_value < 0.05)
}

# ==============================================================================
# Causal Role Grouping
# ==============================================================================
assign_causal_roles <- function(selected_vars, strong_edges, n_groups = 3) {
    out_deg <- strong_edges %>%
        group_by(from) %>%
        summarise(out = n(), .groups = "drop") %>%
        rename(variable = from)
    in_deg <- strong_edges %>%
        group_by(to) %>%
        summarise(inp = n(), .groups = "drop") %>%
        rename(variable = to)
    role_df <- data.frame(variable = selected_vars, stringsAsFactors = FALSE) %>%
        left_join(out_deg, by = "variable") %>%
        left_join(in_deg, by = "variable") %>%
        mutate(
            out = replace_na(out, 0), inp = replace_na(inp, 0),
            role_score = ifelse(inp == 0, out + 1, out / (inp + 1))
        ) %>%
        arrange(desc(role_score))
    n_vars <- nrow(role_df)
    if (n_vars < n_groups) {
        role_df$group <- c("Root", "Mediator", "Terminal")[1:n_vars]
        role_df$group_idx <- 1:n_vars
    } else {
        sizes <- c(
            ceiling(n_vars / n_groups),
            ceiling((n_vars - ceiling(n_vars / n_groups)) / 2),
            0
        )
        sizes[3] <- n_vars - sizes[1] - sizes[2]
        sizes <- pmax(sizes, 0)
        if (sum(sizes) != n_vars) sizes[3] <- n_vars - sizes[1] - sizes[2]
        role_df$group <- c(
            rep("Root", sizes[1]), rep("Mediator", sizes[2]),
            rep("Terminal", sizes[3])
        )
        role_df$group_idx <- c(rep(1L, sizes[1]), rep(2L, sizes[2]), rep(3L, sizes[3]))
    }
    role_df
}

build_adj_matrix <- function(cast_vars, strong_edges) {
    p <- length(cast_vars)
    adj <- matrix(0, nrow = p, ncol = p)
    rownames(adj) <- colnames(adj) <- cast_vars
    for (i in 1:nrow(strong_edges)) {
        if (strong_edges$from[i] %in% cast_vars && strong_edges$to[i] %in% cast_vars) {
            adj[strong_edges$from[i], strong_edges$to[i]] <- 1
        }
    }
    adj
}

# ==============================================================================
# CAST特征工程：全变量ATE加权 + 仅筛选变量间的DAG交互（与01完全一致）
# ==============================================================================
# 基础特征：all_vars（全变量ATE加权）；交互：仅 cast_vars 之间的DAG边
build_cast_features <- function(X_full_sc, all_vars, cast_vars, strong_edges,
                                ate_results, boot_str = NULL) {
    X_base <- X_full_sc[, all_vars, drop = FALSE]
    p <- length(all_vars)

    # 1. ATE加权: 对所有变量施加因果权重（显著变量放大）
    ate_weights <- rep(1.0, p)
    names(ate_weights) <- all_vars
    for (v in all_vars) {
        idx <- which(ate_results$variable == v)
        if (length(idx) > 0 && ate_results$significant[idx[1]]) {
            ate_weights[v] <- 1.0 + abs(ate_results$coef[idx[1]])
        }
    }
    X_weighted <- X_base
    for (v in all_vars) {
        X_weighted[[v]] <- X_weighted[[v]] * ate_weights[v]
    }

    # 2. DAG导向交互: 仅筛选变量(cast_vars)之间的边，避免噪声交互过拟合
    interaction_cols <- list()
    edge_names <- c()
    if (nrow(strong_edges) > 0 && length(cast_vars) > 0) {
        for (k in 1:nrow(strong_edges)) {
            from_v <- strong_edges$from[k]
            to_v <- strong_edges$to[k]
            if (from_v %in% cast_vars && to_v %in% cast_vars &&
                from_v %in% all_vars && to_v %in% all_vars) {
                col_name <- paste0("int_", from_v, "_", to_v)
                edge_weight <- strong_edges$strength[k]
                interaction_cols[[col_name]] <- X_base[[from_v]] *
                    X_base[[to_v]] * edge_weight
                edge_names <- c(edge_names, col_name)
            }
        }
    }

    if (length(interaction_cols) > 0) {
        int_df <- as.data.frame(interaction_cols)
        X_out <- cbind(X_weighted, int_df)
    } else {
        X_out <- X_weighted
    }

    list(
        data = X_out, n_base = p, n_interactions = length(interaction_cols),
        n_total = ncol(X_out), ate_weights = ate_weights,
        interaction_names = edge_names
    )
}

# ==============================================================================
# CI-MLP Architecture (Causally-Informed MLP)
# ==============================================================================
# CAST的核心神经网络: 因果信息以特征形式注入，而非架构约束
# - 与MLP共享完全相同的骨架结构 → 公平对比
# - DAG导向的交互特征: 仅对DAG边 i→j 构建 x_i * x_j
# - ATE加权特征缩放: 因果显著变量被放大
#
CI_MLP <- nn_module("CI_MLP",
    initialize = function(n_input, hidden = 64, dropout = 0.2) {
        self$net <- nn_sequential(
            nn_linear(n_input, hidden), nn_layer_norm(hidden),
            nn_silu(), nn_dropout(dropout),
            nn_linear(hidden, hidden), nn_layer_norm(hidden),
            nn_silu(), nn_dropout(dropout),
            nn_linear(hidden, hidden), nn_layer_norm(hidden),
            nn_silu(), nn_dropout(dropout),
            nn_linear(hidden, as.integer(hidden %/% 2)),
            nn_layer_norm(as.integer(hidden %/% 2)),
            nn_silu(), nn_dropout(dropout * 0.5),
            nn_linear(as.integer(hidden %/% 2), 1L)
        )
    },
    forward = function(x) self$net(x)
)

# MLP: 标准多层感知机（无因果特征）
MLP <- nn_module("MLP",
    initialize = function(n_input, hidden = 64, dropout = 0.2) {
        self$net <- nn_sequential(
            nn_linear(n_input, hidden), nn_layer_norm(hidden),
            nn_silu(), nn_dropout(dropout),
            nn_linear(hidden, hidden), nn_layer_norm(hidden),
            nn_silu(), nn_dropout(dropout),
            nn_linear(hidden, hidden), nn_layer_norm(hidden),
            nn_silu(), nn_dropout(dropout),
            nn_linear(hidden, as.integer(hidden %/% 2)),
            nn_layer_norm(as.integer(hidden %/% 2)),
            nn_silu(), nn_dropout(dropout * 0.5),
            nn_linear(as.integer(hidden %/% 2), 1L)
        )
    },
    forward = function(x) self$net(x)
)

# ==============================================================================
# Focal Loss + Training
# ==============================================================================
focal_loss <- function(logits, targets, alpha = 0.25, gamma = 2.0) {
    bce <- nn_bce_with_logits_loss(reduction = "none")(logits, targets)
    probs <- torch_sigmoid(logits)
    pt <- targets * probs + (1 - targets) * (1 - probs)
    focal_weight <- alpha * (1 - pt)^gamma
    (focal_weight * bce)$mean()
}

train_nn <- function(model, train_dl, val_pred_fn, y_val_vec,
                     epochs = 200, lr = 1e-3, wd = 1e-4, patience = 30,
                     warmup_epochs = 10, focal_alpha = 0.25) {
    optimizer <- optim_adamw(model$parameters, lr = lr, weight_decay = wd)
    best_auc <- 0
    best_state <- NULL
    no_imp <- 0
    nan_count <- 0
    for (epoch in seq_len(epochs)) {
        if (epoch <= warmup_epochs) {
            current_lr <- lr * epoch / warmup_epochs
        } else {
            progress <- (epoch - warmup_epochs) / (epochs - warmup_epochs)
            current_lr <- 1e-5 + 0.5 * (lr - 1e-5) * (1 + cos(pi * progress))
        }
        for (pg in optimizer$param_groups) pg$lr <- current_lr
        model$train()
        tl <- 0
        nb <- 0
        coro::loop(for (batch in train_dl) {
            optimizer$zero_grad()
            logits <- model(batch$x)
            loss <- focal_loss(logits, batch$y, alpha = focal_alpha, gamma = 2.0)
            if (is.nan(loss$item())) {
                nan_count <- nan_count + 1
                if (nan_count > 5) {
                    return(list(model = model, best_val_auc = 0))
                }
                next
            }
            loss$backward()
            nn_utils_clip_grad_norm_(model$parameters, max_norm = 1.0)
            optimizer$step()
            tl <- tl + loss$item()
            nb <- nb + 1
        })
        model$eval()
        with_no_grad({
            vp <- val_pred_fn(model)
        })
        if (any(is.nan(vp))) {
            nan_count <- nan_count + 1
            if (nan_count > 5) {
                return(list(model = model, best_val_auc = 0))
            }
            next
        }
        va <- tryCatch(as.numeric(pROC::auc(pROC::roc(y_val_vec, vp, quiet = TRUE))),
            error = function(e) 0
        )
        if (va > best_auc + 1e-4) {
            best_auc <- va
            best_state <- lapply(model$state_dict(), function(p) p$clone())
            no_imp <- 0
        } else {
            no_imp <- no_imp + 1
        }
        if (no_imp >= patience) break
    }
    if (!is.null(best_state)) model$load_state_dict(best_state)
    list(model = model, best_val_auc = best_auc)
}

# Dataset helper
flat_dataset <- dataset("FlatDS",
    initialize = function(X, y) {
        self$x <- torch_tensor(as.matrix(X), dtype = torch_float())
        self$y <- torch_tensor(y, dtype = torch_float())$unsqueeze(2)
    },
    .getitem = function(i) list(x = self$x[i, ], y = self$y[i, ]),
    .length = function() self$y$size(1)
)

# ==============================================================================
# Traditional SDM helper
# ==============================================================================
train_sdm <- function(sdm_name, X_tr_raw, y_tr_raw, X_te_raw, gam_vars = NULL) {
    switch(sdm_name,
        "RF" = {
            set.seed(42)
            rf <- ranger::ranger(presence ~ .,
                data = cbind(presence = as.factor(y_tr_raw), X_tr_raw),
                num.trees = 1000, probability = TRUE, seed = 42
            )
            predict(rf, data = X_te_raw)$predictions[, "1"]
        },
        "Maxent" = {
            mx <- maxnet::maxnet(
                p = y_tr_raw, data = X_tr_raw,
                maxnet.formula(p = y_tr_raw, data = X_tr_raw)
            )
            as.numeric(predict(mx, X_te_raw, type = "logistic"))
        },
        "BRT" = {
            set.seed(42)
            brt <- gbm::gbm(presence ~ .,
                data = cbind(presence = y_tr_raw, X_tr_raw),
                distribution = "bernoulli", n.trees = 2000, interaction.depth = 5,
                shrinkage = 0.01, cv.folds = 5, verbose = FALSE
            )
            bt <- gbm::gbm.perf(brt, method = "cv", plot.it = FALSE)
            predict(brt, X_te_raw, n.trees = bt, type = "response")
        },
        "GAM" = {
            gv <- if (!is.null(gam_vars)) gam_vars else names(X_tr_raw)
            gf <- as.formula(paste(
                "presence ~",
                paste(paste0("s(", gv, ", k=5)"), collapse = " + ")
            ))
            gm <- mgcv::gam(gf,
                data = cbind(presence = y_tr_raw, X_tr_raw[, gv, drop = FALSE]),
                family = binomial(), method = "REML"
            )
            as.numeric(predict(gm, X_te_raw[, gv, drop = FALSE], type = "response"))
        }
    )
}

# ==============================================================================
# Main Loop — across all regions and species
# ==============================================================================
cat("======================================================================\n")
cat(sprintf("  CAST Multi-Region Experiment: %d regions\n", length(REGIONS)))
cat("  Models: CAST, MLP, RF, Maxent, BRT, GAM\n")
cat("======================================================================\n\n")

all_results <- data.frame()
all_ate_results <- data.frame()
all_dag_info <- data.frame()
all_dag_edges <- data.frame()
all_screening <- data.frame()
all_role_info <- data.frame()

global_sp_idx <- 0

for (region in REGIONS) {
    cat(sprintf("\n╔══════════════════════════════════════════════════════╗\n"))
    cat(sprintf("║  Region: %-4s                                        ║\n", region))
    cat(sprintf("╚══════════════════════════════════════════════════════╝\n"))

    summary_file <- sprintf("output/case2/%s/species_summary.csv", region)
    if (!file.exists(summary_file)) {
        cat(sprintf("  ⚠ Species summary not found for %s → skip\n", region))
        next
    }
    sp_summary <- read.csv(summary_file, stringsAsFactors = FALSE)
    region_species <- sp_summary$species
    cat(sprintf("  %d species to process\n", length(region_species)))

    dir.create(sprintf("output/case2/%s/multi_species", region),
        recursive = TRUE, showWarnings = FALSE
    )

    for (sp_idx in seq_along(region_species)) {
        sp <- region_species[sp_idx]
        global_sp_idx <- global_sp_idx + 1
        cat(sprintf(
            "\n  ─── [%s %d/%d] Species: %s ───\n",
            region, sp_idx, length(region_species), sp
        ))

        # Load data
        train_file <- sprintf("output/case2/%s/train_data_%s.csv", region, sp)
        test_file <- sprintf("output/case2/%s/test_data_%s.csv", region, sp)
        if (!file.exists(train_file) || !file.exists(test_file)) {
            cat("    ⚠ Data files missing → skip\n")
            next
        }
        train_raw <- read.csv(train_file, stringsAsFactors = FALSE)
        test_raw <- read.csv(test_file, stringsAsFactors = FALSE)
        env_cols <- setdiff(names(train_raw), "presence")
        for (col in env_cols) {
            train_raw[[col]] <- as.numeric(train_raw[[col]])
            test_raw[[col]] <- as.numeric(test_raw[[col]])
        }

        # --- Step 1: VIF ---
        X_env <- train_raw[, env_cols, drop = FALSE]
        sds <- apply(X_env, 2, sd, na.rm = TRUE)
        X_env <- X_env[, sds > 1e-10, drop = FALSE]
        repeat {
            vifs <- calc_vif(X_env)
            if (max(vifs, na.rm = TRUE) <= 10 || ncol(X_env) <= 3) break
            X_env <- X_env[, setdiff(names(X_env), names(which.max(vifs))), drop = FALSE]
        }
        selected_vars <- names(X_env)
        train_data <- cbind(
            presence = train_raw$presence,
            train_raw[, selected_vars, drop = FALSE]
        )
        test_data <- cbind(
            presence = test_raw$presence,
            test_raw[, selected_vars, drop = FALSE]
        )
        train_data <- train_data[complete.cases(train_data), ]
        test_data <- test_data[complete.cases(test_data), ]

        # --- Step 2: DAG ---
        env_for_dag <- train_data[, selected_vars, drop = FALSE]
        if (nrow(env_for_dag) > 8000) {
            set.seed(42)
            env_for_dag <- env_for_dag[sample(nrow(env_for_dag), 8000), ]
        }
        set.seed(42)
        boot_str <- bnlearn::boot.strength(env_for_dag,
            R = 100, algorithm = "hc",
            algorithm.args = list(score = "bic-g")
        )
        strong_edges <- boot_str %>% filter(strength >= 0.7, direction >= 0.6)
        n_possible <- length(selected_vars) * (length(selected_vars) - 1) / 2
        dag_density <- nrow(strong_edges) / max(n_possible, 1)
        node_outdeg <- strong_edges %>%
            group_by(from) %>%
            summarise(out_degree = n(), .groups = "drop")

        all_dag_info <- rbind(all_dag_info, data.frame(
            region = region, species = sp, n_edges = nrow(strong_edges),
            dag_density = dag_density, n_vars_after_vif = length(selected_vars),
            stringsAsFactors = FALSE
        ))

        # 保存逐物种DAG边信息（供03出图: Fig S8 跨物种DAG边一致性）
        if (nrow(strong_edges) > 0) {
            sp_edges <- strong_edges %>%
                select(from, to, strength, direction) %>%
                mutate(region = region, species = sp)
            all_dag_edges <- rbind(all_dag_edges, sp_edges)
        }

        # --- Step 3: ATE ---
        Y_full <- train_data$presence
        X_full <- train_data[, selected_vars, drop = FALSE]
        X_full[is.na(X_full)] <- 0
        ate_results <- data.frame()
        for (v in selected_vars) {
            T_bin <- as.integer(X_full[[v]] > median(X_full[[v]], na.rm = TRUE))
            W <- X_full[, setdiff(selected_vars, v), drop = FALSE]
            tryCatch(
                {
                    set.seed(42)
                    res <- dml_ate(Y = Y_full, T_var = T_bin, W = W, K = 2, num_trees = 300)
                    ate_results <- rbind(ate_results, data.frame(
                        variable = v, coef = res$ate, se = res$se,
                        p_value = res$p_value, significant = res$significant,
                        stringsAsFactors = FALSE
                    ))
                },
                error = function(e) {}
            )
        }
        n_sig <- sum(ate_results$significant)
        ate_results$species <- sp
        ate_results$region <- region
        all_ate_results <- rbind(all_ate_results, ate_results)

        # --- Step 4: Adaptive CAST Screening v2 ---
        set.seed(42)
        rf_imp <- ranger::ranger(presence ~ .,
            data = cbind(presence = as.factor(Y_full), X_full),
            num.trees = 500, importance = "permutation", verbose = FALSE
        )$variable.importance

        dag_quality <- 1 - dag_density
        ate_sig_ratio <- if (nrow(ate_results) > 0) {
            sum(ate_results$significant) / nrow(ate_results)
        } else {
            0
        }
        w_dag <- 0.15 + 0.15 * dag_quality
        w_ate <- 0.25 + 0.25 * ate_sig_ratio
        w_imp <- 1 - w_dag - w_ate

        screening_df <- data.frame(variable = selected_vars, stringsAsFactors = FALSE) %>%
            left_join(node_outdeg %>% rename(variable = from), by = "variable") %>%
            left_join(ate_results %>% select(variable, coef, p_value, significant),
                by = "variable"
            ) %>%
            mutate(
                out_degree = replace_na(out_degree, 0),
                abs_ate = abs(replace_na(coef, 0)),
                p_val = replace_na(p_value, 1),
                sig = replace_na(significant, FALSE),
                importance = rf_imp[variable],
                score_dag = normalize01(out_degree),
                score_ate_raw = normalize01(abs_ate),
                ate_penalty = pmin(1.0, -log10(pmax(p_val, 1e-10)) / 3),
                score_ate = score_ate_raw * ate_penalty,
                score_imp = normalize01(importance),
                score_total = w_dag * score_dag + w_ate * score_ate + w_imp * score_imp
            ) %>%
            arrange(desc(score_total))

        min_keep_n <- max(5L, ceiling(length(selected_vars) * 0.5))
        if (length(unique(screening_df$score_total)) >= 2) {
            km <- kmeans(screening_df$score_total, centers = 2, nstart = 10)
            high_cluster <- which.max(km$centers)
            cast_vars_km <- screening_df$variable[km$cluster == high_cluster]
        } else {
            cast_vars_km <- screening_df$variable
        }
        if (length(cast_vars_km) < min_keep_n) {
            cast_vars <- screening_df$variable[1:min(min_keep_n, nrow(screening_df))]
        } else {
            cast_vars <- cast_vars_km
        }

        screening_df$species <- sp
        screening_df$region <- region
        screening_df$w_dag <- w_dag
        screening_df$w_ate <- w_ate
        screening_df$w_imp <- w_imp
        all_screening <- rbind(all_screening, screening_df)

        # --- Step 5: Causal Roles ---
        role_df <- assign_causal_roles(cast_vars, strong_edges, n_groups = 3)
        role_df$species <- sp
        role_df$region <- region
        all_role_info <- rbind(all_role_info, role_df)

        cat(sprintf(
            "    VIF:%d→%d | DAG:%d(d=%.2f) | ATE:%d/%d | screened:%d\n",
            length(env_cols), length(selected_vars), nrow(strong_edges), dag_density,
            n_sig, nrow(ate_results), length(cast_vars)
        ))

        # ═══ Step 6: Modeling ═══
        y_train_all <- train_data$presence
        y_test_all <- test_data$presence

        # 标准化全变量集（所有模型共用）
        X_train_full <- train_data[, selected_vars, drop = FALSE]
        X_test_full <- test_data[, selected_vars, drop = FALSE]
        X_means_full <- colMeans(X_train_full, na.rm = TRUE)
        X_sds_full <- apply(X_train_full, 2, sd, na.rm = TRUE)
        X_sds_full[X_sds_full < 1e-10] <- 1
        X_train_full_sc <- as.data.frame(scale(X_train_full,
            center = X_means_full, scale = X_sds_full
        ))
        X_test_full_sc <- as.data.frame(scale(X_test_full,
            center = X_means_full, scale = X_sds_full
        ))
        X_train_full_sc[is.na(X_train_full_sc)] <- 0
        X_test_full_sc[is.na(X_test_full_sc)] <- 0

        # CAST特征工程：全变量ATE加权 + 仅筛选变量间的DAG交互
        cast_train_info <- build_cast_features(
            X_train_full_sc, selected_vars, cast_vars, strong_edges, ate_results, boot_str
        )
        cast_test_info <- build_cast_features(
            X_test_full_sc, selected_vars, cast_vars, strong_edges, ate_results, boot_str
        )
        X_train_cast <- cast_train_info$data
        X_test_cast <- cast_test_info$data
        n_cast_interactions <- cast_train_info$n_interactions

        # Stratified val split
        set.seed(123)
        pos_idx <- which(y_train_all == 1)
        neg_idx <- which(y_train_all == 0)
        val_pos <- sample(pos_idx, round(0.2 * length(pos_idx)))
        val_neg <- sample(neg_idx, round(0.2 * length(neg_idx)))
        val_idx <- c(val_pos, val_neg)

        y_val <- y_train_all[val_idx]
        y_tr <- y_train_all[-val_idx]
        focal_alpha <- 1 - mean(y_tr)
        batch_size <- min(128L, max(32L, as.integer(length(y_tr) / 100)))

        X_tr_cast <- X_train_cast[-val_idx, ]
        X_val_cast <- X_train_cast[val_idx, ]
        X_tr_full <- X_train_full_sc[-val_idx, ]
        X_val_full <- X_train_full_sc[val_idx, ]

        hidden_size_cast <- max(32L, min(128L, as.integer(cast_train_info$n_total * 4)))
        hidden_size_full <- max(32L, min(128L, as.integer(length(selected_vars) * 8)))
        sp_results <- data.frame()

        # ---- CAST（全变量ATE加权 + DAG交互特征） ----
        {
            aucs <- numeric(n_runs)
            tsss <- numeric(n_runs)
            for (ri in 1:n_runs) {
                torch_manual_seed(seeds[ri])
                set.seed(seeds[ri])
                tryCatch(
                    {
                        ds <- flat_dataset(X_tr_cast, y_tr)
                        dl <- dataloader(ds,
                            batch_size = batch_size,
                            shuffle = TRUE, drop_last = TRUE
                        )
                        vt <- torch_tensor(as.matrix(X_val_cast), dtype = torch_float())
                        m <- CI_MLP(ncol(X_tr_cast), hidden_size_cast, 0.2)
                        res <- train_nn(m, dl,
                            function(m) as.numeric(torch_sigmoid(m(vt))$squeeze()$cpu()),
                            y_val,
                            epochs = 200, patience = 20, focal_alpha = focal_alpha
                        )
                        res$model$eval()
                        with_no_grad({
                            tt <- torch_tensor(as.matrix(X_test_cast), dtype = torch_float())
                            pred <- as.numeric(torch_sigmoid(res$model(tt))$squeeze()$cpu())
                        })
                        ev <- evaluate_model(pred, y_test_all)
                        aucs[ri] <- ev["auc"]
                        tsss[ri] <- ev["tss"]
                    },
                    error = function(e) {
                        cat(sprintf("      CAST run %d FAILED: %s\n", ri, e$message))
                        aucs[ri] <<- NA
                        tsss[ri] <<- NA
                    }
                )
            }
            sp_results <- rbind(sp_results, data.frame(
                region = region, species = sp, model = "CAST", var_set = "full+causal",
                n_vars = length(selected_vars), n_interactions = n_cast_interactions,
                n_features_total = cast_train_info$n_total,
                n_dag_edges = nrow(strong_edges), dag_density = dag_density,
                auc_mean = mean(aucs, na.rm = TRUE), auc_sd = sd(aucs, na.rm = TRUE),
                tss_mean = mean(tsss, na.rm = TRUE), tss_sd = sd(tsss, na.rm = TRUE),
                n_success = sum(!is.na(aucs)), stringsAsFactors = FALSE
            ))
            cat(sprintf(
                "      CAST:  AUC=%.4f±%.4f (%dv+%d int)\n",
                mean(aucs, na.rm = TRUE), sd(aucs, na.rm = TRUE),
                length(selected_vars), n_cast_interactions
            ))
        }

        # ---- MLP_ATE（消融：仅ATE加权，无DAG交互） ----
        ate_only_train_info <- build_cast_features(
            X_train_full_sc, selected_vars, character(0), strong_edges, ate_results, boot_str
        )
        ate_only_test_info <- build_cast_features(
            X_test_full_sc, selected_vars, character(0), strong_edges, ate_results, boot_str
        )
        X_tr_ate <- ate_only_train_info$data[-val_idx, ]
        X_val_ate <- ate_only_train_info$data[val_idx, ]
        X_test_ate <- ate_only_test_info$data
        hidden_size_ate <- max(32L, min(128L, as.integer(ncol(X_tr_ate) * 8)))
        {
            aucs <- numeric(n_runs)
            tsss <- numeric(n_runs)
            for (ri in 1:n_runs) {
                torch_manual_seed(seeds[ri])
                set.seed(seeds[ri])
                tryCatch(
                    {
                        ds <- flat_dataset(X_tr_ate, y_tr)
                        dl <- dataloader(ds,
                            batch_size = batch_size,
                            shuffle = TRUE, drop_last = TRUE
                        )
                        vt <- torch_tensor(as.matrix(X_val_ate), dtype = torch_float())
                        m <- MLP(ncol(X_tr_ate), hidden_size_ate, 0.2)
                        res <- train_nn(m, dl,
                            function(m) as.numeric(torch_sigmoid(m(vt))$squeeze()$cpu()),
                            y_val,
                            epochs = 200, patience = 20, focal_alpha = focal_alpha
                        )
                        res$model$eval()
                        with_no_grad({
                            tt <- torch_tensor(as.matrix(X_test_ate), dtype = torch_float())
                            pred <- as.numeric(torch_sigmoid(res$model(tt))$squeeze()$cpu())
                        })
                        ev <- evaluate_model(pred, y_test_all)
                        aucs[ri] <- ev["auc"]
                        tsss[ri] <- ev["tss"]
                    },
                    error = function(e) {
                        aucs[ri] <<- NA
                        tsss[ri] <<- NA
                    }
                )
            }
            sp_results <- rbind(sp_results, data.frame(
                region = region, species = sp, model = "MLP_ATE", var_set = "full+ate",
                n_vars = length(selected_vars), n_interactions = 0L,
                n_features_total = ncol(X_tr_ate),
                n_dag_edges = NA, dag_density = dag_density,
                auc_mean = mean(aucs, na.rm = TRUE), auc_sd = sd(aucs, na.rm = TRUE),
                tss_mean = mean(tsss, na.rm = TRUE), tss_sd = sd(tsss, na.rm = TRUE),
                n_success = sum(!is.na(aucs)), stringsAsFactors = FALSE
            ))
            cat(sprintf("      MLP_ATE: AUC=%.4f±%.4f\n",
                mean(aucs, na.rm = TRUE), sd(aucs, na.rm = TRUE)))
        }

        # ---- MLP（全变量，无因果增强） ----
        {
            aucs <- numeric(n_runs)
            tsss <- numeric(n_runs)
            for (ri in 1:n_runs) {
                torch_manual_seed(seeds[ri])
                set.seed(seeds[ri])
                tryCatch(
                    {
                        ds <- flat_dataset(X_tr_full, y_tr)
                        dl <- dataloader(ds,
                            batch_size = batch_size,
                            shuffle = TRUE, drop_last = TRUE
                        )
                        vt <- torch_tensor(as.matrix(X_val_full), dtype = torch_float())
                        m <- MLP(length(selected_vars), hidden_size_full, 0.2)
                        res <- train_nn(m, dl,
                            function(m) as.numeric(torch_sigmoid(m(vt))$squeeze()$cpu()),
                            y_val,
                            epochs = 200, patience = 20, focal_alpha = focal_alpha
                        )
                        res$model$eval()
                        with_no_grad({
                            tt <- torch_tensor(as.matrix(X_test_full_sc), dtype = torch_float())
                            pred <- as.numeric(torch_sigmoid(res$model(tt))$squeeze()$cpu())
                        })
                        ev <- evaluate_model(pred, y_test_all)
                        aucs[ri] <- ev["auc"]
                        tsss[ri] <- ev["tss"]
                    },
                    error = function(e) {
                        aucs[ri] <<- NA
                        tsss[ri] <<- NA
                    }
                )
            }
            sp_results <- rbind(sp_results, data.frame(
                region = region, species = sp, model = "MLP", var_set = "full",
                n_vars = length(selected_vars), n_interactions = 0L,
                n_features_total = length(selected_vars),
                n_dag_edges = NA, dag_density = dag_density,
                auc_mean = mean(aucs, na.rm = TRUE), auc_sd = sd(aucs, na.rm = TRUE),
                tss_mean = mean(tsss, na.rm = TRUE), tss_sd = sd(tsss, na.rm = TRUE),
                n_success = sum(!is.na(aucs)), stringsAsFactors = FALSE
            ))
            cat(sprintf(
                "      MLP:   AUC=%.4f±%.4f\n",
                mean(aucs, na.rm = TRUE), sd(aucs, na.rm = TRUE)
            ))
        }

        # ---- 传统SDM (RF, Maxent, BRT, GAM) ——全变量 ----
        y_tr_raw <- train_data$presence
        for (sdm_name in c("RF", "Maxent", "BRT", "GAM")) {
            tryCatch(
                {
                    pred <- train_sdm(sdm_name,
                        train_data[, selected_vars, drop = FALSE],
                        y_tr_raw,
                        test_data[, selected_vars, drop = FALSE],
                        gam_vars = selected_vars
                    )
                    ev <- evaluate_model(pred, y_test_all)
                    sp_results <- rbind(sp_results, data.frame(
                        region = region, species = sp, model = sdm_name,
                        var_set = "full", n_vars = length(selected_vars),
                        n_interactions = 0L, n_features_total = length(selected_vars),
                        n_dag_edges = NA, dag_density = dag_density,
                        auc_mean = ev["auc"], auc_sd = 0,
                        tss_mean = ev["tss"], tss_sd = 0,
                        n_success = 1, stringsAsFactors = FALSE
                    ))
                    cat(sprintf("      %s:    AUC=%.4f\n", sdm_name, ev["auc"]))
                },
                error = function(e) {
                    cat(sprintf(
                        "      %s FAILED: %s\n", sdm_name,
                        substr(e$message, 1, 60)
                    ))
                }
            )
        }

        all_results <- rbind(all_results, sp_results)

        # Checkpoint
        if (global_sp_idx %% 5 == 0) {
            write.csv(all_results, "output/case2/all_results_v3.csv", row.names = FALSE)
            write.csv(all_ate_results, "output/case2/all_ate_results_v3.csv", row.names = FALSE)
            write.csv(all_dag_info, "output/case2/all_dag_info_v3.csv", row.names = FALSE)
            write.csv(all_dag_edges, "output/case2/all_dag_edges_v3.csv", row.names = FALSE)
            write.csv(all_screening, "output/case2/all_screening_v3.csv", row.names = FALSE)
            write.csv(all_role_info, "output/case2/all_role_info_v3.csv", row.names = FALSE)
            cat(sprintf("    [Checkpoint: %d species total saved]\n", global_sp_idx))
        }
    } # end species loop

    # Per-region checkpoint
    write.csv(all_results %>% filter(region == !!region),
        sprintf("output/case2/%s/multi_species/all_results.csv", region),
        row.names = FALSE
    )
} # end region loop

# ==============================================================================
# Final save
# ==============================================================================
write.csv(all_results, "output/case2/all_results_v3.csv", row.names = FALSE)
write.csv(all_ate_results, "output/case2/all_ate_results_v3.csv", row.names = FALSE)
write.csv(all_dag_info, "output/case2/all_dag_info_v3.csv", row.names = FALSE)
write.csv(all_dag_edges, "output/case2/all_dag_edges_v3.csv", row.names = FALSE)
write.csv(all_screening, "output/case2/all_screening_v3.csv", row.names = FALSE)
write.csv(all_role_info, "output/case2/all_role_info_v3.csv", row.names = FALSE)

# ==============================================================================
# Summary
# ==============================================================================
cat("\n======================================================================\n")
cat("  CAST v3 Multi-Region Experiment Complete!\n")
cat("======================================================================\n")

# 按区域汇总
for (r in REGIONS) {
    n_sp <- sum(all_results$region == r & all_results$model == "CAST")
    if (n_sp > 0) {
        cat(sprintf("\n  ── %s (%d species) ──\n", r, n_sp))
        region_summary <- all_results %>%
            filter(region == r) %>%
            group_by(model) %>%
            summarise(
                mean_auc = mean(auc_mean, na.rm = TRUE),
                sd_auc = sd(auc_mean, na.rm = TRUE),
                n = n(), .groups = "drop"
            ) %>%
            arrange(desc(mean_auc))
        for (i in 1:nrow(region_summary)) {
            row <- region_summary[i, ]
            mk <- if (row$model == "CAST") " ★" else "  "
            cat(sprintf(
                "  %s %-12s AUC=%.4f±%.4f (n=%d)\n",
                mk, row$model, row$mean_auc, row$sd_auc, row$n
            ))
        }
    }
}

# 关键诊断: CAST vs 最优基线
cat("\n  ── CAST vs Best Baseline ──\n")
if (nrow(all_results) == 0) {
    cat("  No results to analyze.\n")
} else {
    cast_vs_others <- all_results %>%
        mutate(model_label = ifelse(model == "CAST", "CAST", model)) %>%
        select(region, species, model_label, auc_mean) %>%
        pivot_wider(names_from = model_label, values_from = auc_mean,
                    values_fn = max)

    if ("CAST" %in% names(cast_vs_others)) {
        other_cols <- setdiff(names(cast_vs_others), c("region", "species", "CAST"))
        cast_vs_others$best_other <- apply(cast_vs_others[, other_cols, drop = FALSE], 1,
            function(x) max(x, na.rm = TRUE))
        cast_vs_others$delta <- cast_vs_others$CAST - cast_vs_others$best_other
        valid <- cast_vs_others %>% filter(!is.na(CAST), is.finite(best_other))

        cat(sprintf(
            "  CAST wins in %d/%d species (%.0f%%)\n",
            sum(valid$delta > 0), nrow(valid),
            mean(valid$delta > 0) * 100
        ))
        cat(sprintf("  Mean ΔAUC vs best baseline: %+.4f\n", mean(valid$delta)))
    }

    # 因果增强效应: CAST vs MLP (同为全变量，区别在于因果特征)
    cat("\n  ── Causal Enhancement Effect (CAST vs MLP) ──\n")
    causal_eff <- all_results %>%
        filter(model %in% c("CAST", "MLP")) %>%
        select(region, species, model, auc_mean) %>%
        pivot_wider(names_from = model, values_from = auc_mean)

    if (all(c("CAST", "MLP") %in% names(causal_eff))) {
        causal_eff <- causal_eff %>%
            filter(!is.na(CAST), !is.na(MLP)) %>%
            mutate(delta = CAST - MLP)
        cat(sprintf("  Mean causal enhancement: %+.4f AUC (wins %d/%d, %.0f%%)\n",
            mean(causal_eff$delta), sum(causal_eff$delta > 0), nrow(causal_eff),
            mean(causal_eff$delta > 0) * 100))
    }
} # end if nrow(all_results) > 0

cat("\n  Next: source('scripts/case2/03_publication_figures.R')\n")
