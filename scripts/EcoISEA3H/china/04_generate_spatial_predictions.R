# 04_generate_spatial_predictions.R
# ==============================================================================
# CAST v3 — Spatial Habitat Suitability Predictions for ALL Models (China Region)
#
# Purpose:
#   For each species, retrain ALL 6 models (CAST, MLP_ATE, MLP, RF, Maxent, BRT)
#   using the same pipeline as 03_run_Eco_multi_species.R, then predict habitat
#   suitability scores (HSS) on every hexagon in the China spatial grid.
#
# Output per species:
#   pred_{species}.csv with columns:
#     HID, lon, lat, presence (if known), HSS_CAST, HSS_MLP_ATE, HSS_MLP,
#     HSS_RF, HSS_Maxent, HSS_BRT
#
# Combined output:
#   all_spatial_predictions.csv — long format for all species × models
#
# Input:
#   - Species data: E:/CausalSDMs/outputs/EcoISEA3H/Res9/CAST_ready/species_data_screened/
#   - Full env grid: E:/CausalSDMs/outputs/EcoISEA3H/Res9/CAST_ready/China_EnvData_Res9_Screened.csv
#
# ==============================================================================

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

# ---- Dependencies ----
pkgs <- c("tidyverse", "data.table", "bnlearn", "pROC", "caret",
          "ranger", "maxnet", "gbm", "torch")
for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
}
if (!torch_is_installed()) torch::install_torch()

# ---- Configuration ----
REGION <- "China_Res9"
SEED   <- 42   # single seed for final prediction (not ensemble of 3 runs)

data_dir <- "E:/CausalSDMs/outputs/EcoISEA3H/Res9/CAST_ready/species_data_screened"
env_file <- "E:/CausalSDMs/outputs/EcoISEA3H/Res9/CAST_ready/China_EnvData_Res9_Screened.csv"
out_dir  <- "E:/CausalSDMs/output/case2_eco/spatial_predictions"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Utility Functions (identical to 03_run)
# ==============================================================================
normalize01 <- function(x) {
    r <- range(x, na.rm = TRUE)
    if (r[2] - r[1] < 1e-10) return(rep(0.5, length(x)))
    (x - r[1]) / (r[2] - r[1])
}

# ==============================================================================
# DML ATE (identical to 03_run)
# ==============================================================================
dml_ate <- function(Y, T_var, W, K = 2, num_trees = 300) {
    n <- length(Y)
    folds <- sample(rep(1:K, length.out = n))
    y_res <- numeric(n)
    t_res <- numeric(n)
    for (k in 1:K) {
        train_idx <- which(folds != k)
        test_idx  <- which(folds == k)
        W_train <- W[train_idx, , drop = FALSE]
        W_test  <- W[test_idx, , drop = FALSE]
        rf_y <- ranger::ranger(y ~ ., data = cbind(y = Y[train_idx], W_train),
                               num.trees = num_trees, verbose = FALSE)
        y_res[test_idx] <- Y[test_idx] - predict(rf_y, data = W_test)$predictions
        rf_t <- ranger::ranger(y ~ ., data = cbind(y = T_var[train_idx], W_train),
                               num.trees = num_trees, verbose = FALSE)
        t_res[test_idx] <- T_var[test_idx] - predict(rf_t, data = W_test)$predictions
    }
    ate <- sum(t_res * y_res) / sum(t_res^2)
    residuals <- y_res - ate * t_res
    se <- sqrt(mean(residuals^2 * t_res^2) / (mean(t_res^2)^2) / n)
    p_value <- 2 * pnorm(-abs(ate / se))
    list(ate = ate, se = se, p_value = p_value, significant = p_value < 0.05)
}

# ==============================================================================
# Causal Role Grouping (identical to 03_run)
# ==============================================================================
assign_causal_roles <- function(selected_vars, strong_edges, n_groups = 3) {
    out_deg <- strong_edges %>% group_by(from) %>% summarise(out = n(), .groups = "drop") %>% rename(variable = from)
    in_deg  <- strong_edges %>% group_by(to)   %>% summarise(inp = n(), .groups = "drop") %>% rename(variable = to)
    role_df <- data.frame(variable = selected_vars, stringsAsFactors = FALSE) %>%
        left_join(out_deg, by = "variable") %>% left_join(in_deg, by = "variable") %>%
        mutate(out = replace_na(out, 0), inp = replace_na(inp, 0),
               role_score = ifelse(inp == 0, out + 1, out / (inp + 1))) %>%
        arrange(desc(role_score))
    n_vars <- nrow(role_df)
    if (n_vars < n_groups) {
        role_df$group <- c("Root", "Mediator", "Terminal")[1:n_vars]
    } else {
        sizes <- c(ceiling(n_vars / n_groups), ceiling((n_vars - ceiling(n_vars / n_groups)) / 2), 0)
        sizes[3] <- max(0, n_vars - sizes[1] - sizes[2])
        role_df$group <- c(rep("Root", sizes[1]), rep("Mediator", sizes[2]), rep("Terminal", sizes[3]))
    }
    role_df
}

# ==============================================================================
# CAST Features (identical to 03_run)
# ==============================================================================
build_cast_features <- function(X_full_sc, all_vars, cast_vars, strong_edges, ate_results, boot_str = NULL) {
    X_base <- X_full_sc[, all_vars, drop = FALSE]
    p <- length(all_vars)

    # 1. ATE weighting
    ate_weights <- rep(1.0, p)
    names(ate_weights) <- all_vars
    for (v in all_vars) {
        idx <- which(ate_results$variable == v)
        if (length(idx) > 0 && isTRUE(ate_results$significant[idx[1]])) {
            coef_val <- ate_results$coef[idx[1]]
            if (is.finite(coef_val)) ate_weights[v] <- 1.0 + abs(coef_val)
        }
    }
    X_weighted <- X_base
    for (v in all_vars) X_weighted[[v]] <- X_weighted[[v]] * ate_weights[v]

    # 2. DAG Interactions
    interaction_cols <- list()
    edge_names <- c()
    if (nrow(strong_edges) > 0 && length(cast_vars) > 0) {
        for (k in 1:nrow(strong_edges)) {
            from_v <- strong_edges$from[k]
            to_v   <- strong_edges$to[k]
            if (from_v %in% cast_vars && to_v %in% cast_vars &&
                from_v %in% all_vars && to_v %in% all_vars) {
                col_name <- paste0("int_", from_v, "_", to_v)
                edge_weight <- strong_edges$strength[k]
                interaction_cols[[col_name]] <- X_base[[from_v]] * X_base[[to_v]] * edge_weight
                edge_names <- c(edge_names, col_name)
            }
        }
    }
    if (length(interaction_cols) > 0) {
        X_out <- cbind(X_weighted, as.data.frame(interaction_cols))
    } else {
        X_out <- X_weighted
    }
    list(data = X_out, n_base = p, n_interactions = length(interaction_cols),
         n_total = ncol(X_out), ate_weights = ate_weights, interaction_names = edge_names)
}

# ==============================================================================
# CI-MLP Architecture (identical to 03_run)
# ==============================================================================
CI_MLP <- nn_module("CI_MLP",
    initialize = function(n_input, hidden = 64, dropout = 0.2) {
        self$net <- nn_sequential(
            nn_linear(n_input, hidden), nn_layer_norm(hidden), nn_silu(), nn_dropout(dropout),
            nn_linear(hidden, hidden), nn_layer_norm(hidden), nn_silu(), nn_dropout(dropout),
            nn_linear(hidden, hidden), nn_layer_norm(hidden), nn_silu(), nn_dropout(dropout),
            nn_linear(hidden, as.integer(hidden %/% 2)), nn_layer_norm(as.integer(hidden %/% 2)), nn_silu(), nn_dropout(dropout * 0.5),
            nn_linear(as.integer(hidden %/% 2), 1L)
        )
    },
    forward = function(x) self$net(x)
)

focal_loss <- function(logits, targets, alpha = 0.25, gamma = 2.0) {
    bce <- nn_bce_with_logits_loss(reduction = "none")(logits, targets)
    probs <- torch_sigmoid(logits)
    pt <- targets * probs + (1 - targets) * (1 - probs)
    focal_weight <- alpha * (1 - pt)^gamma
    (focal_weight * bce)$mean()
}

flat_dataset <- dataset("FlatDS",
    initialize = function(X, y) {
        self$x <- torch_tensor(as.matrix(X), dtype = torch_float())
        self$y <- torch_tensor(y, dtype = torch_float())$unsqueeze(2)
    },
    .getitem = function(i) list(x = self$x[i, ], y = self$y[i, ]),
    .length = function() self$y$size(1)
)

train_nn <- function(model, train_dl, val_pred_fn, y_val_vec,
                     epochs = 200, lr = 1e-3, wd = 1e-4, patience = 30,
                     warmup_epochs = 10, focal_alpha = 0.25) {
    optimizer <- optim_adamw(model$parameters, lr = lr, weight_decay = wd)
    best_auc <- 0; best_state <- NULL; no_imp <- 0
    for (epoch in seq_len(epochs)) {
        current_lr <- if (epoch <= warmup_epochs) lr * epoch / warmup_epochs else
            1e-5 + 0.5 * (lr - 1e-5) * (1 + cos(pi * (epoch - warmup_epochs) / (epochs - warmup_epochs)))
        for (pg in optimizer$param_groups) pg$lr <- current_lr
        model$train()
        coro::loop(for (batch in train_dl) {
            optimizer$zero_grad()
            logits <- model(batch$x)
            loss <- focal_loss(logits, batch$y, alpha = focal_alpha, gamma = 2.0)
            if (is.nan(loss$item())) next
            loss$backward()
            nn_utils_clip_grad_norm_(model$parameters, max_norm = 1.0)
            optimizer$step()
        })
        model$eval()
        with_no_grad({ vp <- val_pred_fn(model) })
        va <- if (any(is.nan(vp))) 0 else tryCatch(
            as.numeric(pROC::auc(pROC::roc(y_val_vec, vp, quiet = TRUE))), error = function(e) 0)
        if (va > best_auc + 1e-4) {
            best_auc <- va; best_state <- lapply(model$state_dict(), function(p) p$clone()); no_imp <- 0
        } else { no_imp <- no_imp + 1 }
        if (no_imp >= patience) break
    }
    if (!is.null(best_state)) model$load_state_dict(best_state)
    list(model = model, best_val_auc = best_auc)
}

# ==============================================================================
# Prediction helper: predict with NN on arbitrary data
# ==============================================================================
predict_nn <- function(model, X_df) {
    model$eval()
    with_no_grad({
        xt <- torch_tensor(as.matrix(X_df), dtype = torch_float())
        pred <- as.numeric(torch_sigmoid(model(xt))$squeeze()$cpu())
    })
    pred
}

# ==============================================================================
# Traditional SDM train + predict-on-grid helper
# ==============================================================================
train_and_predict_sdm <- function(sdm_name, X_tr_raw, y_tr_raw, X_pred_grid) {
    switch(sdm_name,
        "RF" = {
            set.seed(42)
            m <- ranger::ranger(presence ~ ., data = cbind(presence = as.factor(y_tr_raw), X_tr_raw),
                                num.trees = 300, num.threads = 8, probability = TRUE, seed = 42)
            predict(m, data = X_pred_grid)$predictions[, "1"]
        },
        "Maxent" = {
            mx <- maxnet::maxnet(p = y_tr_raw, data = X_tr_raw,
                                 maxnet.formula(p = y_tr_raw, data = X_tr_raw))
            as.numeric(predict(mx, X_pred_grid, type = "logistic"))
        },
        "BRT" = {
            set.seed(42)
            brt <- gbm::gbm(presence ~ ., data = cbind(presence = y_tr_raw, X_tr_raw),
                             distribution = "bernoulli", n.trees = 500, interaction.depth = 5,
                             shrinkage = 0.01, cv.folds = 5, n.cores = 8, verbose = FALSE)
            bt <- gbm::gbm.perf(brt, method = "cv", plot.it = FALSE)
            predict(brt, X_pred_grid, n.trees = bt, type = "response")
        }
    )
}

# ==============================================================================
# Load full spatial grid (all hexagons for prediction)
# ==============================================================================
cat("  Loading full spatial grid...\n")
env_grid <- fread(env_file)
grid_hids <- env_grid$HID
grid_lons <- env_grid$lon
grid_lats <- env_grid$lat
cat(sprintf("  Full grid: %d hexagons\n", nrow(env_grid)))

# ==============================================================================
# Main Loop: For each species, train all models → predict on full grid
# ==============================================================================
sp_files <- list.files(data_dir, pattern = "^CAST_.*_screened\\.csv$", full.names = TRUE)
cat(sprintf("  Found %d species datasets.\n\n", length(sp_files)))

# Combined long-format collector
all_predictions_long <- data.frame()

for (sp_idx in seq_along(sp_files)) {
    f <- sp_files[sp_idx]
    sp_name_raw <- gsub("CAST_|_Res9_screened\\.csv$", "", basename(f))
    sp <- sp_name_raw

    # Check if already done
    pred_file <- file.path(out_dir, paste0("pred_", sp, ".csv"))
    if (file.exists(pred_file)) {
        cat(sprintf("  [%d/%d] %s — [Skip: already done]\n", sp_idx, length(sp_files), sp))
        next
    }

    cat(sprintf("  ═══ [%d/%d] Species: %s ═══\n", sp_idx, length(sp_files), sp))

    # ---- Load species data ----
    sp_df <- fread(f)
    meta_cols <- c("HID", "lon", "lat", "species", "sid", "family", "category", "presence", "fraction")
    env_cols <- setdiff(names(sp_df), meta_cols)

    # ---- Same train/test split as 03_run ----
    set.seed(42)
    train_idx <- sample(1:nrow(sp_df), size = round(0.7 * nrow(sp_df)))
    train_data <- sp_df[train_idx, ]
    y_train_all <- train_data$presence

    # ---- Full spatial grid env data (for prediction) ----
    X_grid_raw <- as.data.frame(env_grid[, ..env_cols, drop = FALSE])
    X_grid_raw[is.na(X_grid_raw)] <- 0

    # ---- Scaling parameters from training data ----
    X_train_full <- as.data.frame(train_data[, ..env_cols, drop = FALSE])
    X_train_full[is.na(X_train_full)] <- 0
    X_means <- colMeans(X_train_full, na.rm = TRUE)
    X_sds   <- apply(X_train_full, 2, sd, na.rm = TRUE)
    X_sds[X_sds < 1e-10] <- 1

    X_train_full_sc <- as.data.frame(scale(X_train_full, center = X_means, scale = X_sds))
    X_train_full_sc[is.na(X_train_full_sc)] <- 0
    X_grid_sc <- as.data.frame(scale(X_grid_raw, center = X_means, scale = X_sds))
    X_grid_sc[is.na(X_grid_sc)] <- 0

    selected_vars <- env_cols

    # ==================================================================
    # Step 2: DAG Learning (same as 03_run)
    # ==================================================================
    env_for_dag <- train_data[, ..selected_vars, drop = FALSE]
    env_for_dag_df <- as.data.frame(env_for_dag)
    for (col in names(env_for_dag_df)) env_for_dag_df[[col]] <- as.numeric(env_for_dag_df[[col]])
    env_for_dag_df <- na.omit(env_for_dag_df)
    if (nrow(env_for_dag_df) < 10) { cat("    Skip: too few complete cases.\n"); next }
    if (nrow(env_for_dag_df) > 8000) env_for_dag_df <- env_for_dag_df[sample(nrow(env_for_dag_df), 8000), ]

    set.seed(42)
    boot_str <- bnlearn::boot.strength(env_for_dag_df, R = 100, algorithm = "hc",
                                        algorithm.args = list(score = "bic-g"))
    strong_edges <- boot_str %>% filter(strength >= 0.7, direction >= 0.6)

    n_possible <- length(selected_vars) * (length(selected_vars) - 1) / 2
    dag_density <- nrow(strong_edges) / max(n_possible, 1)

    node_outdeg <- strong_edges %>% group_by(from) %>% summarise(out_degree = n(), .groups = "drop")

    # ==================================================================
    # Step 3: ATE Estimation (same as 03_run)
    # ==================================================================
    Y_full <- train_data$presence
    X_full <- as.data.frame(train_data[, ..selected_vars, drop = FALSE])
    X_full[is.na(X_full)] <- 0

    ate_results <- data.frame()
    for (v in selected_vars) {
        T_bin <- as.integer(X_full[[v]] > median(X_full[[v]], na.rm = TRUE))
        W <- X_full[, setdiff(selected_vars, v), drop = FALSE]
        tryCatch({
            set.seed(42)
            res <- dml_ate(Y = Y_full, T_var = T_bin, W = W, K = 2, num_trees = 200)
            ate_results <- rbind(ate_results, data.frame(
                variable = v, coef = res$ate, se = res$se,
                p_value = res$p_value, significant = res$significant,
                stringsAsFactors = FALSE))
        }, error = function(e) {})
    }

    # ==================================================================
    # Step 4: Adaptive CAST Screening (same as 03_run)
    # ==================================================================
    set.seed(42)
    rf_imp <- ranger::ranger(presence ~ ., data = cbind(presence = as.factor(Y_full), X_full),
                             num.trees = 300, num.threads = 8, importance = "permutation",
                             verbose = FALSE)$variable.importance
    dag_quality <- 1 - dag_density
    ate_sig_ratio <- if (nrow(ate_results) > 0) sum(ate_results$significant) / nrow(ate_results) else 0
    w_dag <- 0.15 + 0.15 * dag_quality
    w_ate <- 0.25 + 0.25 * ate_sig_ratio
    w_imp <- 1 - w_dag - w_ate

    screening_df <- data.frame(variable = selected_vars, stringsAsFactors = FALSE) %>%
        left_join(node_outdeg %>% rename(variable = from), by = "variable") %>%
        left_join(if (nrow(ate_results) > 0) ate_results %>% select(variable, coef, p_value, significant)
                  else data.frame(variable = character(0), coef = numeric(0),
                                  p_value = numeric(0), significant = logical(0)),
                  by = "variable") %>%
        mutate(out_degree = replace_na(out_degree, 0), abs_ate = abs(replace_na(coef, 0)),
               p_val = replace_na(p_value, 1), sig = replace_na(significant, FALSE),
               importance = rf_imp[variable],
               score_dag = normalize01(out_degree), score_ate_raw = normalize01(abs_ate),
               ate_penalty = pmin(1.0, -log10(pmax(p_val, 1e-10)) / 3),
               score_ate = score_ate_raw * ate_penalty, score_imp = normalize01(importance),
               score_total = w_dag * score_dag + w_ate * score_ate + w_imp * score_imp) %>%
        arrange(desc(score_total))

    min_keep_n <- max(5L, ceiling(length(selected_vars) * 0.5))
    if (length(unique(screening_df$score_total)) >= 2) {
        km <- kmeans(screening_df$score_total, centers = 2, nstart = 10)
        high_cluster <- which.max(km$centers)
        cast_vars_km <- screening_df$variable[km$cluster == high_cluster]
    } else { cast_vars_km <- screening_df$variable }
    cast_vars <- if (length(cast_vars_km) < min_keep_n)
        screening_df$variable[1:min(min_keep_n, nrow(screening_df))] else cast_vars_km

    cat(sprintf("    DAG:%d(d=%.2f) | ATE:%d/%d | screened:%d\n",
                nrow(strong_edges), dag_density,
                sum(ate_results$significant), length(selected_vars), length(cast_vars)))

    # ==================================================================
    # Build features for CAST and MLP_ATE on full grid
    # ==================================================================
    # CAST features (full causal: ATE weighting + DAG interactions)
    cast_train_info <- build_cast_features(X_train_full_sc, selected_vars, cast_vars,
                                            strong_edges, ate_results, boot_str)
    cast_grid_info  <- build_cast_features(X_grid_sc, selected_vars, cast_vars,
                                            strong_edges, ate_results, boot_str)

    # MLP_ATE features (ATE weighting only, no interactions)
    ate_train_info <- build_cast_features(X_train_full_sc, selected_vars, character(0),
                                           strong_edges, ate_results, boot_str)
    ate_grid_info  <- build_cast_features(X_grid_sc, selected_vars, character(0),
                                           strong_edges, ate_results, boot_str)

    # ==================================================================
    # Train/Val split for NN models
    # ==================================================================
    set.seed(123)
    pos_idx <- which(y_train_all == 1)
    neg_idx <- which(y_train_all == 0)
    val_pos <- sample(pos_idx, round(0.2 * length(pos_idx)))
    val_neg <- sample(neg_idx, round(0.2 * length(neg_idx)))
    val_idx <- c(val_pos, val_neg)
    y_val <- y_train_all[val_idx]
    y_tr  <- y_train_all[-val_idx]
    focal_alpha <- 1 - mean(y_tr)
    batch_size <- min(128L, max(32L, as.integer(length(y_tr) / 100)))

    hidden_size_cast <- max(32L, min(128L, as.integer(cast_train_info$n_total * 4)))
    hidden_size_full <- max(32L, min(128L, as.integer(length(selected_vars) * 8)))

    # ---- Initialize prediction output ----
    pred_out <- data.frame(
        HID = grid_hids, lon = grid_lons, lat = grid_lats,
        stringsAsFactors = FALSE
    )

    # Add presence info: mark hexagons where species is known present/absent
    sp_presence <- sp_df[, .(HID, presence)]
    pred_out <- pred_out %>% left_join(sp_presence, by = "HID")

    # ==================================================================
    # Model 1: CAST (CI-MLP with causal feature engineering)
    # ==================================================================
    cat("    Training CAST...")
    tryCatch({
        torch_manual_seed(SEED); set.seed(SEED)
        X_tr_cast <- cast_train_info$data[-val_idx, ]
        X_val_cast <- cast_train_info$data[val_idx, ]
        ds <- flat_dataset(X_tr_cast, y_tr)
        dl <- dataloader(ds, batch_size = batch_size, shuffle = TRUE, drop_last = TRUE)
        vt <- torch_tensor(as.matrix(X_val_cast), dtype = torch_float())
        m <- CI_MLP(ncol(X_tr_cast), hidden_size_cast, 0.2)
        res <- train_nn(m, dl, function(m) as.numeric(torch_sigmoid(m(vt))$squeeze()$cpu()),
                        y_val, epochs = 200, patience = 20, focal_alpha = focal_alpha)
        pred_out$HSS_CAST <- predict_nn(res$model, cast_grid_info$data)
        cat(sprintf(" val_AUC=%.4f ✓\n", res$best_val_auc))
    }, error = function(e) {
        pred_out$HSS_CAST <<- NA_real_
        cat(sprintf(" FAILED: %s\n", e$message))
    })

    # ==================================================================
    # Model 2: MLP_ATE (CI-MLP with ATE weighting only, no interactions)
    # ==================================================================
    cat("    Training MLP_ATE...")
    tryCatch({
        torch_manual_seed(SEED); set.seed(SEED)
        X_tr_ate <- ate_train_info$data[-val_idx, ]
        X_val_ate <- ate_train_info$data[val_idx, ]
        ds <- flat_dataset(X_tr_ate, y_tr)
        dl <- dataloader(ds, batch_size = batch_size, shuffle = TRUE, drop_last = TRUE)
        vt <- torch_tensor(as.matrix(X_val_ate), dtype = torch_float())
        m <- CI_MLP(ncol(X_tr_ate), hidden_size_full, 0.2)
        res <- train_nn(m, dl, function(m) as.numeric(torch_sigmoid(m(vt))$squeeze()$cpu()),
                        y_val, epochs = 200, patience = 20, focal_alpha = focal_alpha)
        pred_out$HSS_MLP_ATE <- predict_nn(res$model, ate_grid_info$data)
        cat(sprintf(" val_AUC=%.4f ✓\n", res$best_val_auc))
    }, error = function(e) {
        pred_out$HSS_MLP_ATE <<- NA_real_
        cat(sprintf(" FAILED: %s\n", e$message))
    })

    # ==================================================================
    # Model 3: MLP (vanilla, no causal features)
    # ==================================================================
    cat("    Training MLP...")
    tryCatch({
        torch_manual_seed(SEED); set.seed(SEED)
        X_tr_full <- X_train_full_sc[-val_idx, ]
        X_val_full <- X_train_full_sc[val_idx, ]
        ds <- flat_dataset(X_tr_full, y_tr)
        dl <- dataloader(ds, batch_size = batch_size, shuffle = TRUE, drop_last = TRUE)
        vt <- torch_tensor(as.matrix(X_val_full), dtype = torch_float())
        m <- CI_MLP(ncol(X_tr_full), hidden_size_full, 0.2)
        res <- train_nn(m, dl, function(m) as.numeric(torch_sigmoid(m(vt))$squeeze()$cpu()),
                        y_val, epochs = 200, patience = 20, focal_alpha = focal_alpha)
        pred_out$HSS_MLP <- predict_nn(res$model, X_grid_sc)
        cat(sprintf(" val_AUC=%.4f ✓\n", res$best_val_auc))
    }, error = function(e) {
        pred_out$HSS_MLP <<- NA_real_
        cat(sprintf(" FAILED: %s\n", e$message))
    })

    # ==================================================================
    # Model 4: RF (Random Forest)
    # ==================================================================
    cat("    Training RF...")
    tryCatch({
        pred_out$HSS_RF <- train_and_predict_sdm("RF", X_train_full, y_train_all, X_grid_raw)
        cat(" ✓\n")
    }, error = function(e) {
        pred_out$HSS_RF <<- NA_real_
        cat(sprintf(" FAILED: %s\n", e$message))
    })

    # ==================================================================
    # Model 5: Maxent
    # ==================================================================
    cat("    Training Maxent...")
    tryCatch({
        pred_out$HSS_Maxent <- train_and_predict_sdm("Maxent", X_train_full, y_train_all, X_grid_raw)
        cat(" ✓\n")
    }, error = function(e) {
        pred_out$HSS_Maxent <<- NA_real_
        cat(sprintf(" FAILED: %s\n", e$message))
    })

    # ==================================================================
    # Model 6: BRT (Boosted Regression Trees)
    # ==================================================================
    cat("    Training BRT...")
    tryCatch({
        pred_out$HSS_BRT <- train_and_predict_sdm("BRT", X_train_full, y_train_all, X_grid_raw)
        cat(" ✓\n")
    }, error = function(e) {
        pred_out$HSS_BRT <<- NA_real_
        cat(sprintf(" FAILED: %s\n", e$message))
    })

    # ==================================================================
    # Save per-species wide-format prediction file
    # ==================================================================
    write.csv(pred_out, pred_file, row.names = FALSE)
    cat(sprintf("    → Saved: %s (%d hexagons × %d models)\n\n", basename(pred_file), nrow(pred_out), 6))

    # Collect into long format for combined file
    model_cols <- grep("^HSS_", names(pred_out), value = TRUE)
    for (mc in model_cols) {
        model_name <- gsub("^HSS_", "", mc)
        sp_long <- data.frame(
            region  = REGION,
            species = sp,
            model   = model_name,
            HID     = pred_out$HID,
            lon     = pred_out$lon,
            lat     = pred_out$lat,
            presence = pred_out$presence,
            hss     = pred_out[[mc]],
            stringsAsFactors = FALSE
        )
        all_predictions_long <- rbind(all_predictions_long, sp_long)
    }

    # Save combined checkpoint
    write.csv(all_predictions_long,
              file.path(out_dir, "all_spatial_predictions.csv"),
              row.names = FALSE)

    gc()
}

cat("\n======================================================================\n")
cat("  All Spatial Predictions Complete!\n")
cat(sprintf("  Output directory: %s\n", out_dir))
cat(sprintf("  Per-species files: pred_{species}.csv (wide format)\n"))
cat(sprintf("  Combined file: all_spatial_predictions.csv (long format)\n"))
cat("======================================================================\n")
