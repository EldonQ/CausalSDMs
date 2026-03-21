# 05_generate_multiseed_predictions.R
# ==============================================================================
# CAST — Multi-Seed Spatial Predictions for Uncertainty Analysis
#
# Purpose:
#   For each species, retrain CAST, MLP_ATE, and MLP with multiple seeds
#   to generate prediction uncertainty estimates. Traditional SDMs (RF, Maxent,
#   BRT) are deterministic given a fixed seed and are skipped.
#
# Output per species:
#   pred_{species}_multiseed.csv with columns:
#     HID, lon, lat, presence,
#     HSS_CAST_s42, HSS_CAST_s71, HSS_CAST_s103, HSS_CAST_mean, HSS_CAST_sd,
#     HSS_MLP_s42, HSS_MLP_s71, HSS_MLP_s103, HSS_MLP_mean, HSS_MLP_sd,
#     (same for MLP_ATE)
#
# This script complements 04_generate_spatial_predictions.R by providing
# within-model uncertainty from stochastic NN training.
#
# Reference: GNN-SDM (Wu et al., 2025) BMU sensitivity + 95% CI
# ==============================================================================

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

# ---- Dependencies ----
pkgs <- c("tidyverse", "data.table", "bnlearn", "pROC", "caret",
          "ranger", "maxnet", "gbm", "torch", "grf")
for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
}
if (!torch_is_installed()) torch::install_torch()

device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
cat(sprintf("  Computing device: %s\n", as.character(device)))

# ---- Configuration ----
REGION <- "China_Res9"
SEEDS  <- c(42, 71, 103)  # Multiple seeds for uncertainty

data_dir <- "E:/CausalSDMs/outputs/EcoISEA3H/Res9/CAST_ready/species_data_screened"
env_file <- "E:/CausalSDMs/outputs/EcoISEA3H/Res9/CAST_ready/China_EnvData_Res9_Screened.csv"
out_dir  <- "E:/CausalSDMs/output/case2_eco/spatial_predictions"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Only process these species (set to NULL for all)
TARGET_SPECIES <- c("Rhinopithecus_roxellana", "Ovis_ammon", "Macaca_mulatta")

# ==============================================================================
# Utility Functions (identical to 04)
# ==============================================================================
normalize01 <- function(x) {
    r <- range(x, na.rm = TRUE)
    if (r[2] - r[1] < 1e-10) return(rep(0.5, length(x)))
    (x - r[1]) / (r[2] - r[1])
}

dml_ate <- function(Y, T_var, W, K = 2, num_trees = 300) {
    n <- length(Y)
    folds <- sample(rep(1:K, length.out = n))
    y_res <- numeric(n); t_res <- numeric(n)
    for (k in 1:K) {
        train_idx <- which(folds != k); test_idx <- which(folds == k)
        W_train <- W[train_idx, , drop = FALSE]; W_test <- W[test_idx, , drop = FALSE]
        rf_y <- ranger::ranger(y ~ ., data = cbind(y = Y[train_idx], W_train), num.trees = num_trees, verbose = FALSE)
        y_res[test_idx] <- Y[test_idx] - predict(rf_y, data = W_test)$predictions
        rf_t <- ranger::ranger(y ~ ., data = cbind(y = T_var[train_idx], W_train), num.trees = num_trees, verbose = FALSE)
        t_res[test_idx] <- T_var[test_idx] - predict(rf_t, data = W_test)$predictions
    }
    ate <- sum(t_res * y_res) / sum(t_res^2)
    residuals <- y_res - ate * t_res
    se <- sqrt(mean(residuals^2 * t_res^2) / (mean(t_res^2)^2) / n)
    p_value <- 2 * pnorm(-abs(ate / se))
    list(ate = ate, se = se, p_value = p_value, significant = p_value < 0.05)
}

build_cast_features <- function(X_full_sc, all_vars, cast_vars, strong_edges, ate_results, boot_str = NULL) {
    X_base <- X_full_sc[, all_vars, drop = FALSE]; p <- length(all_vars)
    ate_weights <- rep(1.0, p); names(ate_weights) <- all_vars
    for (v in all_vars) {
        idx <- which(ate_results$variable == v)
        if (length(idx) > 0 && isTRUE(ate_results$significant[idx[1]])) {
            coef_val <- ate_results$coef[idx[1]]
            if (is.finite(coef_val)) ate_weights[v] <- 1.0 + abs(coef_val)
        }
    }
    X_weighted <- X_base
    for (v in all_vars) X_weighted[[v]] <- X_weighted[[v]] * ate_weights[v]
    interaction_cols <- list(); edge_names <- c()
    if (nrow(strong_edges) > 0 && length(cast_vars) > 0) {
        for (k in 1:nrow(strong_edges)) {
            from_v <- strong_edges$from[k]; to_v <- strong_edges$to[k]
            if (from_v %in% cast_vars && to_v %in% cast_vars && from_v %in% all_vars && to_v %in% all_vars) {
                col_name <- paste0("int_", from_v, "_", to_v)
                interaction_cols[[col_name]] <- X_base[[from_v]] * X_base[[to_v]] * strong_edges$strength[k]
                edge_names <- c(edge_names, col_name)
            }
        }
    }
    X_out <- if (length(interaction_cols) > 0) cbind(X_weighted, as.data.frame(interaction_cols)) else X_weighted
    list(data = X_out, n_base = p, n_interactions = length(interaction_cols),
         n_total = ncol(X_out), ate_weights = ate_weights, interaction_names = edge_names)
}

# ==============================================================================
# CI-MLP + Training (identical to 04)
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
    (alpha * (1 - pt)^gamma * bce)$mean()
}

flat_dataset <- dataset("FlatDS",
    initialize = function(X, y) {
        self$x <- torch_tensor(as.matrix(X), dtype = torch_float(), device = device)
        self$y <- torch_tensor(y, dtype = torch_float(), device = device)$unsqueeze(2)
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

predict_nn <- function(model, X_df) {
    model$eval()
    with_no_grad({
        xt <- torch_tensor(as.matrix(X_df), dtype = torch_float(), device = device)
        pred <- as.numeric(torch_sigmoid(model(xt))$squeeze()$cpu())
    })
    pred
}

# ==============================================================================
# Main Loop
# ==============================================================================
cat("======================================================================\n")
cat("  Multi-Seed Predictions for Uncertainty Analysis\n")
cat("======================================================================\n")

env_grid <- fread(env_file)
grid_hids <- env_grid$HID; grid_lons <- env_grid$lon; grid_lats <- env_grid$lat
cat(sprintf("  Full grid: %d hexagons\n", nrow(env_grid)))

sp_files <- list.files(data_dir, pattern = "^CAST_.*_screened\\.csv$", full.names = TRUE)

# Filter to target species if specified
if (!is.null(TARGET_SPECIES)) {
    sp_files <- sp_files[sapply(sp_files, function(f) {
        sp_raw <- gsub("CAST_|_Res9_screened\\.csv$", "", basename(f))
        sp_raw %in% TARGET_SPECIES
    })]
}
cat(sprintf("  Species to process: %d\n\n", length(sp_files)))

for (sp_idx in seq_along(sp_files)) {
    f <- sp_files[sp_idx]
    sp <- gsub("CAST_|_Res9_screened\\.csv$", "", basename(f))
    ms_file <- file.path(out_dir, paste0("pred_", sp, "_multiseed.csv"))

    if (file.exists(ms_file)) {
        cat(sprintf("  [%d/%d] %s — [Skip: already done]\n", sp_idx, length(sp_files), sp))
        next
    }

    cat(sprintf("  ═══ [%d/%d] Species: %s ═══\n", sp_idx, length(sp_files), sp))

    sp_df <- fread(f)
    meta_cols <- c("HID", "lon", "lat", "species", "sid", "family", "category", "presence", "fraction")
    env_cols <- setdiff(names(sp_df), meta_cols)

    set.seed(42)
    train_idx <- sample(1:nrow(sp_df), size = round(0.7 * nrow(sp_df)))
    train_data <- sp_df[train_idx, ]
    y_train_all <- train_data$presence

    X_grid_raw <- as.data.frame(env_grid[, ..env_cols, drop = FALSE])
    X_grid_raw[is.na(X_grid_raw)] <- 0
    X_train_full <- as.data.frame(train_data[, ..env_cols, drop = FALSE])
    X_train_full[is.na(X_train_full)] <- 0
    X_means <- colMeans(X_train_full, na.rm = TRUE)
    X_sds <- apply(X_train_full, 2, sd, na.rm = TRUE); X_sds[X_sds < 1e-10] <- 1
    X_train_full_sc <- as.data.frame(scale(X_train_full, center = X_means, scale = X_sds))
    X_train_full_sc[is.na(X_train_full_sc)] <- 0
    X_grid_sc <- as.data.frame(scale(X_grid_raw, center = X_means, scale = X_sds))
    X_grid_sc[is.na(X_grid_sc)] <- 0
    selected_vars <- env_cols

    # ---- DAG + ATE (run once, same for all seeds) ----
    env_for_dag <- train_data[, c(selected_vars, "presence"), with = FALSE]
    env_for_dag_df <- as.data.frame(env_for_dag)
    for (col in names(env_for_dag_df)) env_for_dag_df[[col]] <- as.numeric(env_for_dag_df[[col]])
    env_for_dag_df <- na.omit(env_for_dag_df)
    if (nrow(env_for_dag_df) < 10) { cat("    Skip: too few cases\n"); next }
    if (nrow(env_for_dag_df) > 8000) env_for_dag_df <- env_for_dag_df[sample(nrow(env_for_dag_df), 8000), ]

    set.seed(42)
    boot_str <- bnlearn::boot.strength(env_for_dag_df, R = 100, algorithm = "hc",
                                         algorithm.args = list(score = "bic-g"))
    strong_edges <- boot_str %>% filter(strength >= 0.7, direction >= 0.6)
    strong_env_edges <- strong_edges %>% filter(from != "presence" & to != "presence")
    dag_density <- nrow(strong_env_edges) / max(length(selected_vars) * (length(selected_vars) - 1) / 2, 1)
    node_outdeg <- strong_env_edges %>% group_by(from) %>% summarise(out_degree = n(), .groups = "drop")

    Y_full <- train_data$presence
    X_full <- as.data.frame(train_data[, ..selected_vars, drop = FALSE]); X_full[is.na(X_full)] <- 0
    ate_results <- data.frame()
    for (v in selected_vars) {
        T_bin <- as.integer(X_full[[v]] > median(X_full[[v]], na.rm = TRUE))
        W <- X_full[, setdiff(selected_vars, v), drop = FALSE]
        tryCatch({ set.seed(42)
            res <- dml_ate(Y = Y_full, T_var = T_bin, W = W, K = 2, num_trees = 200)
            ate_results <- rbind(ate_results, data.frame(variable = v, coef = res$ate, se = res$se,
                p_value = res$p_value, significant = res$significant, stringsAsFactors = FALSE))
        }, error = function(e) {})
    }

    # Screening
    set.seed(42)
    rf_imp <- ranger::ranger(presence ~ ., data = cbind(presence = as.factor(Y_full), X_full),
                              num.trees = 300, num.threads = 8, importance = "permutation", verbose = FALSE)$variable.importance
    dag_quality <- 1 - dag_density
    ate_sig_ratio <- if (nrow(ate_results) > 0) sum(ate_results$significant) / nrow(ate_results) else 0
    w_dag <- 0.15 + 0.15 * dag_quality; w_ate <- 0.25 + 0.25 * ate_sig_ratio; w_imp <- 1 - w_dag - w_ate
    screening_df <- data.frame(variable = selected_vars, stringsAsFactors = FALSE) %>%
        left_join(node_outdeg %>% rename(variable = from), by = "variable") %>%
        left_join(if (nrow(ate_results) > 0) ate_results %>% select(variable, coef, p_value, significant)
                  else data.frame(variable = character(0), coef = numeric(0), p_value = numeric(0), significant = logical(0)), by = "variable") %>%
        mutate(out_degree = replace_na(out_degree, 0), abs_ate = abs(replace_na(coef, 0)),
               p_val = replace_na(p_value, 1), sig = replace_na(significant, FALSE), importance = rf_imp[variable],
               score_dag = normalize01(out_degree), score_ate_raw = normalize01(abs_ate),
               ate_penalty = pmin(1.0, -log10(pmax(p_val, 1e-10)) / 3),
               score_ate = score_ate_raw * ate_penalty, score_imp = normalize01(importance),
               score_total = w_dag * score_dag + w_ate * score_ate + w_imp * score_imp) %>%
        arrange(desc(score_total))
    min_keep_n <- max(5L, ceiling(length(selected_vars) * 0.5))
    if (length(unique(screening_df$score_total)) >= 2) {
        km <- kmeans(screening_df$score_total, centers = 2, nstart = 10)
        cast_vars_km <- screening_df$variable[km$cluster == which.max(km$centers)]
    } else { cast_vars_km <- screening_df$variable }
    cast_vars <- if (length(cast_vars_km) < min_keep_n)
        screening_df$variable[1:min(min_keep_n, nrow(screening_df))] else cast_vars_km

    # Build features
    cast_train_info <- build_cast_features(X_train_full_sc, selected_vars, cast_vars, strong_env_edges, ate_results, boot_str)
    cast_grid_info  <- build_cast_features(X_grid_sc, selected_vars, cast_vars, strong_env_edges, ate_results, boot_str)
    ate_train_info  <- build_cast_features(X_train_full_sc, selected_vars, character(0), strong_env_edges, ate_results, boot_str)
    ate_grid_info   <- build_cast_features(X_grid_sc, selected_vars, character(0), strong_env_edges, ate_results, boot_str)

    set.seed(123)
    pos_idx <- which(y_train_all == 1); neg_idx <- which(y_train_all == 0)
    val_pos <- sample(pos_idx, round(0.2 * length(pos_idx)))
    val_neg <- sample(neg_idx, round(0.2 * length(neg_idx)))
    val_idx <- c(val_pos, val_neg)
    y_val <- y_train_all[val_idx]; y_tr <- y_train_all[-val_idx]
    focal_alpha <- 1 - mean(y_tr)
    batch_size <- min(128L, max(32L, as.integer(length(y_tr) / 100)))

    hidden_size_cast <- max(32L, min(128L, as.integer(cast_train_info$n_total * 4)))
    hidden_size_full <- max(32L, min(128L, as.integer(length(selected_vars) * 8)))

    # ---- Initialize output ----
    pred_out <- data.frame(HID = grid_hids, lon = grid_lons, lat = grid_lats, stringsAsFactors = FALSE)
    sp_presence <- sp_df[, .(HID, presence)]
    pred_out <- pred_out %>% left_join(sp_presence, by = "HID")

    cat(sprintf("    DAG:%d(d=%.2f) | ATE:%d/%d | screened:%d\n",
                nrow(strong_env_edges), dag_density, sum(ate_results$significant), length(selected_vars), length(cast_vars)))

    # ==================================================================
    # Multi-seed training loop
    # ==================================================================
    nn_configs <- list(
        list(name = "CAST",    X_tr = cast_train_info$data[-val_idx, ], X_val = cast_train_info$data[val_idx, ], X_grid = cast_grid_info$data,  hidden = hidden_size_cast),
        list(name = "MLP_ATE", X_tr = ate_train_info$data[-val_idx, ],  X_val = ate_train_info$data[val_idx, ],  X_grid = ate_grid_info$data,   hidden = hidden_size_full),
        list(name = "MLP",     X_tr = X_train_full_sc[-val_idx, ],      X_val = X_train_full_sc[val_idx, ],      X_grid = X_grid_sc,            hidden = hidden_size_full)
    )

    for (cfg in nn_configs) {
        cat(sprintf("    Training %s (3 seeds)...", cfg$name))
        seed_preds <- list()

        for (si in seq_along(SEEDS)) {
            s <- SEEDS[si]
            tryCatch({
                torch_manual_seed(s); set.seed(s)
                ds <- flat_dataset(cfg$X_tr, y_tr)
                dl <- dataloader(ds, batch_size = batch_size, shuffle = TRUE, drop_last = TRUE)
                vt <- torch_tensor(as.matrix(cfg$X_val), dtype = torch_float(), device = device)
                m <- CI_MLP(ncol(cfg$X_tr), cfg$hidden, 0.2)
                m$to(device = device)
                res <- train_nn(m, dl, function(m) as.numeric(torch_sigmoid(m(vt))$squeeze()$cpu()),
                                y_val, epochs = 200, patience = 20, focal_alpha = focal_alpha)
                seed_preds[[si]] <- predict_nn(res$model, cfg$X_grid)
            }, error = function(e) {
                seed_preds[[si]] <<- rep(NA_real_, nrow(cfg$X_grid))
            })
        }

        # Save per-seed and aggregate columns
        for (si in seq_along(SEEDS)) {
            col_name <- paste0("HSS_", cfg$name, "_s", SEEDS[si])
            pred_out[[col_name]] <- seed_preds[[si]]
        }
        seed_mat <- do.call(cbind, seed_preds)
        pred_out[[paste0("HSS_", cfg$name, "_mean")]] <- rowMeans(seed_mat, na.rm = TRUE)
        pred_out[[paste0("HSS_", cfg$name, "_sd")]]   <- apply(seed_mat, 1, sd, na.rm = TRUE)
        cat(" ✓\n")
    }

    # Save
    write.csv(pred_out, ms_file, row.names = FALSE)
    cat(sprintf("    → Saved: %s (%d hexagons × %d seeds × 3 models)\n\n",
                basename(ms_file), nrow(pred_out), length(SEEDS)))
    gc()
    if (cuda_is_available()) torch::cuda_empty_cache()
}

cat("\n======================================================================\n")
cat("  Multi-Seed Predictions Complete!\n")
cat(sprintf("  Output directory: %s\n", out_dir))
cat("======================================================================\n")
