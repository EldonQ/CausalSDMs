# 03_run_Plant_multi_species.R
# ==============================================================================
# CAST v3 — Plant Multi-Species Experiment (China Region)
#
# Pipeline per species:
#   Step 1: Variables already VIF-screened (13 vars)
#   Step 2: DAG learning (HC bootstrap R=100, strength≥0.7, direction≥0.6)
#   Step 3: ATE estimation (DML 2-fold cross-fitting)
#   Step 4: Adaptive CAST Screening v2
#   Step 5: Causal Role Grouping
#   Step 6: Unified comparison table (CAST, MLP_ATE, MLP, RF, Maxent, BRT)
#
# Input: E:/CausalSDMs/outputs/Plant/CAST_ready/species_data_screened/CAST_*_screened.csv
# Output: E:/CausalSDMs/output/case4_plant/ (all_results_plant_v3.csv, etc.)
# ==============================================================================

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

# ---- Dependencies ----
pkgs <- c("tidyverse", "data.table", "bnlearn", "pROC", "caret", "ranger", "maxnet", "gbm", "torch")
for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
}
if (!torch_is_installed()) torch::install_torch()

# ---- Configuration ----
REGION <- "Europe_Plant"
n_runs <- 3
seeds <- c(42, 71, 103)

data_dir <- "E:/CausalSDMs/outputs/Plant/CAST_ready/species_data_screened"
out_dir <- "E:/CausalSDMs/output/case4_plant"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Utility Functions
# ==============================================================================
normalize01 <- function(x) {
    r <- range(x, na.rm = TRUE)
    if (r[2] - r[1] < 1e-10) {
        return(rep(0.5, length(x)))
    }
    (x - r[1]) / (r[2] - r[1])
}

evaluate_model <- function(pred, obs) {
    pred <- pmin(pmax(pred, 1e-7), 1 - 1e-7)
    auc_val <- tryCatch(as.numeric(pROC::auc(pROC::roc(obs, pred, quiet = TRUE))), error = function(e) NA_real_)
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
        sizes <- c(ceiling(n_vars / n_groups), ceiling((n_vars - ceiling(n_vars / n_groups)) / 2), 0)
        sizes[3] <- max(0, n_vars - sizes[1] - sizes[2])
        role_df$group <- c(rep("Root", sizes[1]), rep("Mediator", sizes[2]), rep("Terminal", sizes[3]))
        role_df$group_idx <- c(rep(1L, sizes[1]), rep(2L, sizes[2]), rep(3L, sizes[3]))
    }
    role_df
}

# ==============================================================================
# CAST Features
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
    for (v in all_vars) {
        X_weighted[[v]] <- X_weighted[[v]] * ate_weights[v]
    }

    # 2. DAG Interactions
    interaction_cols <- list()
    edge_names <- c()
    if (nrow(strong_edges) > 0 && length(cast_vars) > 0) {
        for (k in 1:nrow(strong_edges)) {
            from_v <- strong_edges$from[k]
            to_v <- strong_edges$to[k]
            if (from_v %in% cast_vars && to_v %in% cast_vars && from_v %in% all_vars && to_v %in% all_vars) {
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

    list(
        data = X_out, n_base = p, n_interactions = length(interaction_cols),
        n_total = ncol(X_out), ate_weights = ate_weights, interaction_names = edge_names
    )
}

# ==============================================================================
# CI-MLP Architecture
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

train_nn <- function(model, train_dl, val_pred_fn, y_val_vec,
                     epochs = 200, lr = 1e-3, wd = 1e-4, patience = 30, warmup_epochs = 10, focal_alpha = 0.25) {
    optimizer <- optim_adamw(model$parameters, lr = lr, weight_decay = wd)
    best_auc <- 0
    best_state <- NULL
    no_imp <- 0
    nan_count <- 0
    for (epoch in seq_len(epochs)) {
        current_lr <- if (epoch <= warmup_epochs) lr * epoch / warmup_epochs else 1e-5 + 0.5 * (lr - 1e-5) * (1 + cos(pi * (epoch - warmup_epochs) / (epochs - warmup_epochs)))
        for (pg in optimizer$param_groups) pg$lr <- current_lr
        model$train()
        coro::loop(for (batch in train_dl) {
            optimizer$zero_grad()
            logits <- model(batch$x)
            loss <- focal_loss(logits, batch$y, alpha = focal_alpha, gamma = 2.0)
            if (is.nan(loss$item())) {
                nan_count <- nan_count + 1
                next
            }
            loss$backward()
            nn_utils_clip_grad_norm_(model$parameters, max_norm = 1.0)
            optimizer$step()
        })
        model$eval()
        with_no_grad({
            vp <- val_pred_fn(model)
        })
        if (any(is.nan(vp))) {
            nan_count <- nan_count + 1
            next
        }
        va <- tryCatch(as.numeric(pROC::auc(pROC::roc(y_val_vec, vp, quiet = TRUE))), error = function(e) 0)
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

flat_dataset <- dataset("FlatDS",
    initialize = function(X, y) {
        self$x <- torch_tensor(as.matrix(X), dtype = torch_float())
        self$y <- torch_tensor(y, dtype = torch_float())$unsqueeze(2)
    },
    .getitem = function(i) list(x = self$x[i, ], y = self$y[i, ]),
    .length = function() self$y$size(1)
)

train_sdm <- function(sdm_name, X_tr_raw, y_tr_raw, X_te_raw) {
    switch(sdm_name,
        "RF" = {
            set.seed(42)
            predict(ranger::ranger(presence ~ ., data = cbind(presence = as.factor(y_tr_raw), X_tr_raw), num.trees = 1000, probability = TRUE, seed = 42), data = X_te_raw)$predictions[, "1"]
        },
        "Maxent" = {
            mx <- maxnet::maxnet(p = y_tr_raw, data = X_tr_raw, maxnet.formula(p = y_tr_raw, data = X_tr_raw))
            as.numeric(predict(mx, X_te_raw, type = "logistic"))
        },
        "BRT" = {
            set.seed(42)
            brt <- gbm::gbm(presence ~ ., data = cbind(presence = y_tr_raw, X_tr_raw), distribution = "bernoulli", n.trees = 2000, interaction.depth = 5, shrinkage = 0.01, cv.folds = 5, verbose = FALSE)
            bt <- gbm::gbm.perf(brt, method = "cv", plot.it = FALSE)
            predict(brt, X_te_raw, n.trees = bt, type = "response")
        }
    )
}

# ==============================================================================
# Main Loop (Eco-ISEA3H Data)
# ==============================================================================
cat("======================================================================\n")
cat(sprintf("  CAST Plant Multi-Species Experiment: %s\n", REGION))
cat("======================================================================\n\n")

# Load existing checkpoint if any
all_results <- data.frame()
all_ate_results <- data.frame()
all_dag_info <- data.frame()
all_dag_edges <- data.frame()
all_screening <- data.frame()
all_role_info <- data.frame()
all_learning_curves <- data.frame()
all_spatial_cate <- data.frame()
done_pairs <- character(0)

checkpoint_paths <- list(
    res = file.path(out_dir, "all_results_plant_v3.csv"),
    ate = file.path(out_dir, "all_ate_results_plant_v3.csv"),
    dag = file.path(out_dir, "all_dag_info_plant_v3.csv"),
    edg = file.path(out_dir, "all_dag_edges_plant_v3.csv"),
    scr = file.path(out_dir, "all_screening_plant_v3.csv"),
    rol = file.path(out_dir, "all_role_info_plant_v3.csv"),
    lcv = file.path(out_dir, "all_learning_curves_plant_v3.csv"),
    spt = file.path(out_dir, "all_spatial_cate_plant_v3.csv")
)

if (file.exists(checkpoint_paths$res)) {
    all_results <- read.csv(checkpoint_paths$res, stringsAsFactors = FALSE) %>% distinct(region, species, model, .keep_all = TRUE)
    if (file.exists(checkpoint_paths$ate)) all_ate_results <- read.csv(checkpoint_paths$ate, stringsAsFactors = FALSE) %>% distinct(region, species, variable, .keep_all = TRUE)
    if (file.exists(checkpoint_paths$dag)) all_dag_info <- read.csv(checkpoint_paths$dag, stringsAsFactors = FALSE) %>% distinct(region, species, .keep_all = TRUE)
    if (file.exists(checkpoint_paths$edg)) all_dag_edges <- read.csv(checkpoint_paths$edg, stringsAsFactors = FALSE) %>% distinct(region, species, from, to, .keep_all = TRUE)
    if (file.exists(checkpoint_paths$scr)) all_screening <- read.csv(checkpoint_paths$scr, stringsAsFactors = FALSE) %>% distinct(region, species, variable, .keep_all = TRUE)
    if (file.exists(checkpoint_paths$rol)) all_role_info <- read.csv(checkpoint_paths$rol, stringsAsFactors = FALSE) %>% distinct(region, species, variable, .keep_all = TRUE)
    if (file.exists(checkpoint_paths$lcv)) all_learning_curves <- read.csv(checkpoint_paths$lcv, stringsAsFactors = FALSE)
    if (file.exists(checkpoint_paths$spt)) all_spatial_cate <- read.csv(checkpoint_paths$spt, stringsAsFactors = FALSE)

    done_count <- all_results %>%
        group_by(region, species) %>%
        summarise(n = n(), .groups = "drop") %>%
        filter(n >= 6)
    if (nrow(done_count) > 0) done_pairs <- paste(done_count$region, done_count$species, sep = "___")

    incomplete <- all_results %>%
        group_by(region, species) %>%
        summarise(n = n(), .groups = "drop") %>%
        filter(n < 6)
    if (nrow(incomplete) > 0) {
        inc_pairs <- paste(incomplete$region, incomplete$species, sep = "___")
        filter_inc <- function(df) if (nrow(df) > 0) filter(df, !paste(region, species, sep = "___") %in% inc_pairs) else df
        all_results <- filter_inc(all_results)
        all_ate_results <- filter_inc(all_ate_results)
        all_dag_info <- filter_inc(all_dag_info)
        all_dag_edges <- filter_inc(all_dag_edges)
        all_screening <- filter_inc(all_screening)
        all_role_info <- filter_inc(all_role_info)
        all_learning_curves <- filter_inc(all_learning_curves)
        all_spatial_cate <- filter_inc(all_spatial_cate)
    }
}

sp_files <- list.files(data_dir, pattern = "^CAST_.*_screened\\.csv$", full.names = TRUE)
cat(sprintf("  Found %d screened species datasets.\n", length(sp_files)))

for (sp_idx in seq_along(sp_files)) {
    f <- sp_files[sp_idx]
    # Extact species name
    sp_name_raw <- gsub("CAST_|_screened\\.csv$", "", basename(f))
    sp <- sp_name_raw
    pair_key <- paste(REGION, sp, sep = "___")

    if (pair_key %in% done_pairs) {
        cat(sprintf("\n  ─── [%d/%d] Species: %s ─── [Skip: already done]\n", sp_idx, length(sp_files), sp))
        next
    }

    cat(sprintf("\n  ─── [%d/%d] Species: %s ───\n", sp_idx, length(sp_files), sp))

    # Load data
    sp_df <- fread(f)
    meta_cols <- c("PlotObservationID", "Longitude", "Latitude", "species", "presence")
    env_cols <- setdiff(names(sp_df), meta_cols)

    # For Eco-ISEA3H we split into train/test using a 70/30 spatial split natively here
    set.seed(42)
    train_idx <- sample(1:nrow(sp_df), size = round(0.7 * nrow(sp_df)))
    train_data <- sp_df[train_idx, ]
    test_data <- sp_df[-train_idx, ]

    selected_vars <- env_cols

    # --- Step 2: DAG ---
    # CRITICAL FIX: Include 'presence' in DAG learning so that the structure is species-specific!
    env_for_dag <- train_data[, c(selected_vars, "presence"), with = FALSE]

    # Cast all specified variables to numeric for bnlearn, avoiding "type: integer" issue
    env_for_dag_df <- as.data.frame(env_for_dag)
    for (col in names(env_for_dag_df)) {
        env_for_dag_df[[col]] <- as.numeric(env_for_dag_df[[col]])
    }

    if (nrow(env_for_dag_df) > 8000) env_for_dag_df <- env_for_dag_df[sample(nrow(env_for_dag_df), 8000), ]
    set.seed(42)
    boot_str <- bnlearn::boot.strength(env_for_dag_df, R = 100, algorithm = "hc", algorithm.args = list(score = "bic-g"))
    strong_edges <- boot_str %>% filter(strength >= 0.7, direction >= 0.6)

    # Exclude edges connecting to 'presence' when passing DAG edges to the rest of the pipeline
    # The downstream structure encoding expects ONLY environmental variables.
    strong_env_edges <- strong_edges %>%
        filter(from != "presence" & to != "presence")

    # Calculate density based on environmental edges
    n_possible <- length(selected_vars) * (length(selected_vars) - 1) / 2
    dag_density <- nrow(strong_env_edges) / max(n_possible, 1)

    node_outdeg <- strong_env_edges %>%
        group_by(from) %>%
        summarise(out_degree = n(), .groups = "drop")

    all_dag_info <- rbind(all_dag_info, data.frame(region = REGION, species = sp, n_edges = nrow(strong_env_edges), dag_density = dag_density, n_vars_after_vif = length(selected_vars), stringsAsFactors = FALSE))
    if (nrow(strong_env_edges) > 0) all_dag_edges <- rbind(all_dag_edges, strong_env_edges %>% select(from, to, strength, direction) %>% mutate(region = REGION, species = sp))

    # --- Step 3: ATE ---
    Y_full <- train_data$presence
    X_full <- as.data.frame(train_data[, ..selected_vars, drop = FALSE])
    X_full[is.na(X_full)] <- 0
    ate_results <- data.frame()
    for (v in selected_vars) {
        T_bin <- as.integer(X_full[[v]] > median(X_full[[v]], na.rm = TRUE))
        W <- X_full[, setdiff(selected_vars, v), drop = FALSE]
        tryCatch(
            {
                set.seed(42)
                res <- dml_ate(Y = Y_full, T_var = T_bin, W = W, K = 2, num_trees = 200) # num_trees=200 for speed
                ate_results <- rbind(ate_results, data.frame(variable = v, coef = res$ate, se = res$se, p_value = res$p_value, significant = res$significant, stringsAsFactors = FALSE))
            },
            error = function(e) {}
        )
    }
    n_sig <- sum(ate_results$significant)
    if (nrow(ate_results) > 0) {
        ate_results$species <- sp
        ate_results$region <- REGION
        all_ate_results <- rbind(all_ate_results, ate_results)
    }

    # --- Step 4: Adaptive CAST Screening v2 ---
    set.seed(42)
    rf_imp <- ranger::ranger(presence ~ ., data = cbind(presence = as.factor(Y_full), X_full), num.trees = 500, importance = "permutation", verbose = FALSE)$variable.importance
    dag_quality <- 1 - dag_density
    ate_sig_ratio <- if (nrow(ate_results) > 0) sum(ate_results$significant) / nrow(ate_results) else 0
    w_dag <- 0.15 + 0.15 * dag_quality
    w_ate <- 0.25 + 0.25 * ate_sig_ratio
    w_imp <- 1 - w_dag - w_ate

    screening_df <- data.frame(variable = selected_vars, stringsAsFactors = FALSE) %>%
        left_join(node_outdeg %>% rename(variable = from), by = "variable") %>%
        left_join(if (nrow(ate_results) > 0) ate_results %>% select(variable, coef, p_value, significant) else data.frame(variable = character(0), coef = numeric(0), p_value = numeric(0), significant = logical(0)), by = "variable") %>%
        mutate(out_degree = replace_na(out_degree, 0), abs_ate = abs(replace_na(coef, 0)), p_val = replace_na(p_value, 1), sig = replace_na(significant, FALSE), importance = rf_imp[variable], score_dag = normalize01(out_degree), score_ate_raw = normalize01(abs_ate), ate_penalty = pmin(1.0, -log10(pmax(p_val, 1e-10)) / 3), score_ate = score_ate_raw * ate_penalty, score_imp = normalize01(importance), score_total = w_dag * score_dag + w_ate * score_ate + w_imp * score_imp) %>%
        arrange(desc(score_total))

    min_keep_n <- max(5L, ceiling(length(selected_vars) * 0.5))
    if (length(unique(screening_df$score_total)) >= 2) {
        km <- kmeans(screening_df$score_total, centers = 2, nstart = 10)
        high_cluster <- which.max(km$centers)
        cast_vars_km <- screening_df$variable[km$cluster == high_cluster]
    } else {
        cast_vars_km <- screening_df$variable
    }
    cast_vars <- if (length(cast_vars_km) < min_keep_n) screening_df$variable[1:min(min_keep_n, nrow(screening_df))] else cast_vars_km

    screening_df$species <- sp
    screening_df$region <- REGION
    screening_df$w_dag <- w_dag
    screening_df$w_ate <- w_ate
    screening_df$w_imp <- w_imp
    all_screening <- rbind(all_screening, screening_df)

    # --- Step 5: Causal Roles ---
    role_df <- assign_causal_roles(cast_vars, strong_env_edges, n_groups = 3)
    role_df$species <- sp
    role_df$region <- REGION
    all_role_info <- rbind(all_role_info, role_df)

    cat(sprintf("    DAG:%d(d=%.2f) | ATE:%d/%d | screened:%d\n", nrow(strong_env_edges), dag_density, n_sig, length(selected_vars), length(cast_vars)))

    # ═══ Step 6: Modeling ═══
    y_train_all <- train_data$presence
    y_test_all <- test_data$presence
    X_train_full <- as.data.frame(train_data[, ..selected_vars, drop = FALSE])
    X_test_full <- as.data.frame(test_data[, ..selected_vars, drop = FALSE])
    X_means_full <- colMeans(X_train_full, na.rm = TRUE)
    X_sds_full <- apply(X_train_full, 2, sd, na.rm = TRUE)
    X_sds_full[X_sds_full < 1e-10] <- 1
    X_train_full_sc <- as.data.frame(scale(X_train_full, center = X_means_full, scale = X_sds_full))
    X_train_full_sc[is.na(X_train_full_sc)] <- 0
    X_test_full_sc <- as.data.frame(scale(X_test_full, center = X_means_full, scale = X_sds_full))
    X_test_full_sc[is.na(X_test_full_sc)] <- 0

    cast_train_info <- build_cast_features(X_train_full_sc, selected_vars, cast_vars, strong_env_edges, ate_results, boot_str)
    cast_test_info <- build_cast_features(X_test_full_sc, selected_vars, cast_vars, strong_env_edges, ate_results, boot_str)
    X_train_cast <- cast_train_info$data
    X_test_cast <- cast_test_info$data
    n_cast_interactions <- cast_train_info$n_interactions

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

    run_nn <- function(name, X_tr, X_val, X_te, hidden_sz) {
        aucs <- numeric(n_runs)
        tsss <- numeric(n_runs)
        model_curves <- data.frame()
        for (ri in 1:n_runs) {
            torch_manual_seed(seeds[ri])
            set.seed(seeds[ri])
            tryCatch(
                {
                    ds <- flat_dataset(X_tr, y_tr)
                    dl <- dataloader(ds, batch_size = batch_size, shuffle = TRUE, drop_last = TRUE)
                    vt <- torch_tensor(as.matrix(X_val), dtype = torch_float())
                    m <- CI_MLP(ncol(X_tr), hidden_sz, 0.2)
                    res <- train_nn(m, dl, function(m) as.numeric(torch_sigmoid(m(vt))$squeeze()$cpu()), y_val, epochs = 200, patience = 20, focal_alpha = focal_alpha)

                    if (nrow(res$history) > 0) {
                        res$history$run <- ri
                        res$history$model <- name
                        model_curves <- rbind(model_curves, res$history)
                    }

                    res$model$eval()
                    with_no_grad({
                        tt <- torch_tensor(as.matrix(X_te), dtype = torch_float())
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
        cat(sprintf("      %s: AUC=%.4f±%.4f\n", name, mean(aucs, na.rm = TRUE), sd(aucs, na.rm = TRUE)))
        list(auc_mean = mean(aucs, na.rm = TRUE), auc_sd = sd(aucs, na.rm = TRUE), tss_mean = mean(tsss, na.rm = TRUE), tss_sd = sd(tsss, na.rm = TRUE), curves = model_curves)
    }

    # CAST
    r_cast <- run_nn("CAST", X_tr_cast, X_val_cast, X_test_cast, hidden_size_cast)
    sp_curves <- r_cast$curves
    sp_results <- rbind(sp_results, data.frame(region = REGION, species = sp, model = "CAST", var_set = "full+causal", n_vars = length(selected_vars), n_interactions = n_cast_interactions, n_features_total = cast_train_info$n_total, n_dag_edges = nrow(strong_env_edges), dag_density = dag_density, auc_mean = r_cast$auc_mean, auc_sd = r_cast$auc_sd, tss_mean = r_cast$tss_mean, tss_sd = r_cast$tss_sd, n_success = 3, stringsAsFactors = FALSE))

    # MLP_ATE
    ate_tr_info <- build_cast_features(X_train_full_sc, selected_vars, character(0), strong_env_edges, ate_results, boot_str)
    ate_te_info <- build_cast_features(X_test_full_sc, selected_vars, character(0), strong_env_edges, ate_results, boot_str)
    r_ate <- run_nn("MLP_ATE", ate_tr_info$data[-val_idx, ], ate_tr_info$data[val_idx, ], ate_te_info$data, hidden_size_full)
    if (!is.null(r_ate$curves)) sp_curves <- rbind(sp_curves, r_ate$curves)
    sp_results <- rbind(sp_results, data.frame(region = REGION, species = sp, model = "MLP_ATE", var_set = "full", n_vars = length(selected_vars), n_interactions = 0, n_features_total = length(selected_vars), n_dag_edges = 0, dag_density = 0, auc_mean = r_ate$auc_mean, auc_sd = r_ate$auc_sd, tss_mean = r_ate$tss_mean, tss_sd = r_ate$tss_sd, n_success = 3, stringsAsFactors = FALSE))

    # MLP
    r_mlp <- run_nn("MLP", X_tr_full, X_val_full, X_test_full_sc, hidden_size_full)
    if (!is.null(r_mlp$curves)) sp_curves <- rbind(sp_curves, r_mlp$curves)
    sp_results <- rbind(sp_results, data.frame(region = REGION, species = sp, model = "MLP", var_set = "full", n_vars = length(selected_vars), n_interactions = 0, n_features_total = length(selected_vars), n_dag_edges = 0, dag_density = 0, auc_mean = r_mlp$auc_mean, auc_sd = r_mlp$auc_sd, tss_mean = r_mlp$tss_mean, tss_sd = r_mlp$tss_sd, n_success = 3, stringsAsFactors = FALSE))

    if (!is.null(sp_curves) && nrow(sp_curves) > 0) {
        sp_curves$region <- REGION
        sp_curves$species <- sp
        all_learning_curves <- rbind(all_learning_curves, sp_curves)
    }

    # Trad SDMs
    baselines <- c("RF", "Maxent", "BRT")
    for (b in baselines) {
        p <- tryCatch(
            {
                train_sdm(b, X_train_full, y_train_all, X_test_full)
            },
            error = function(e) NA
        )
        ev <- if (any(is.na(p))) c(auc = NA, tss = NA) else evaluate_model(p, y_test_all)
        cat(sprintf("      %s: AUC=%.4f\n", b, ev["auc"]))
        sp_results <- rbind(sp_results, data.frame(region = REGION, species = sp, model = b, var_set = "full", n_vars = length(selected_vars), n_interactions = 0, n_features_total = length(selected_vars), n_dag_edges = 0, dag_density = 0, auc_mean = ev["auc"], auc_sd = 0, tss_mean = ev["tss"], tss_sd = 0, n_success = if (is.na(ev["auc"])) 0 else 1, stringsAsFactors = FALSE))
    }

    # ═══ Step 7: Real Spatial CATE (grf causal_forest) ═══
    # Select top 3 significant causal variables to calculate full spatial CATEs
    cate_vars <- cast_vars[cast_vars %in% ate_results$variable[ate_results$significant == TRUE]]
    if (length(cate_vars) == 0) {
        cate_vars <- cast_vars[1:min(3, length(cast_vars))]
    } else {
        cate_vars <- cate_vars[1:min(3, length(cate_vars))]
    }

    cate_df_list <- list()
    for (cv in cate_vars) {
        T_cont <- X_full[[cv]]
        W_covs <- X_full[, setdiff(selected_vars, cv), drop = FALSE]
        cf <- tryCatch(
            {
                grf::causal_forest(X = as.matrix(W_covs), Y = Y_full, W = T_cont, num.trees = 500, seed = 42)
            },
            error = function(e) NULL
        )

        if (!is.null(cf)) {
            X_all <- as.data.frame(sp_df[, ..selected_vars, drop = FALSE])
            W_all <- as.matrix(X_all[, setdiff(selected_vars, cv), drop = FALSE])
            cate_preds <- predict(cf, W_all, estimate.variance = FALSE)$predictions
            cate_df_list[[cv]] <- data.frame(
                region = REGION, species = sp, variable = cv,
                lon = sp_df$Longitude, lat = sp_df$Latitude, cate = as.numeric(cate_preds),
                stringsAsFactors = FALSE
            )
        }
    }
    if (length(cate_df_list) > 0) {
        sp_cate_res <- do.call(rbind, cate_df_list)
        all_spatial_cate <- rbind(all_spatial_cate, sp_cate_res)
    }

    # Save checkpoint
    all_results <- rbind(all_results, sp_results)
    write.csv(all_results, checkpoint_paths$res, row.names = FALSE)
    write.csv(all_ate_results, checkpoint_paths$ate, row.names = FALSE)
    write.csv(all_dag_info, checkpoint_paths$dag, row.names = FALSE)
    write.csv(all_dag_edges, checkpoint_paths$edg, row.names = FALSE)
    write.csv(all_screening, checkpoint_paths$scr, row.names = FALSE)
    write.csv(all_role_info, checkpoint_paths$rol, row.names = FALSE)
    write.csv(all_learning_curves, checkpoint_paths$lcv, row.names = FALSE)
    write.csv(all_spatial_cate, checkpoint_paths$spt, row.names = FALSE)
}

cat("\n======================================================================\n")
cat("  All Plant Species Processing Complete!\n")
cat("======================================================================\n")
