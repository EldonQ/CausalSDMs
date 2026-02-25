################################################################################
# 01_cast_pipeline.R
# CAST v3 Single-Species Pipeline — disdat benchmark (any region, any species)
#
# This script is a single-species version of 02_multi_species_experiment.R.
# All functions, architectures, and logic are IDENTICAL to 02.
#
# ► Set `TARGET_REGION` and `TARGET_SPECIES` below, then source() this script.
#
# ══════════════════════════════════════════════════════════════════════════════
# AVAILABLE REGIONS AND SPECIES (disdat benchmark, Elith et al. 2020)
# ══════════════════════════════════════════════════════════════════════════════
#
# ── Region: AWT (Australian Wet Tropics) ──
#   Taxa: birds, bats, reptiles
#   Env vars: 13 (bc01–bc06, cti, li_meso, mi_ann, radnwet, radnwint, rainann, ruession)
#   Species (30):
#     awt01, awt02, awt03, awt04, awt05, awt06, awt07, awt08, awt09, awt10,
#     awt11, awt12, awt13, awt14, awt15, awt16, awt17, awt18, awt19, awt20,
#     awt21, awt22, awt23, awt24, awt25, awt26, awt27, awt28, awt29, awt30
#   Usage:
#     TARGET_REGION  <- "AWT"
#     TARGET_SPECIES <- "awt01"
#
# ── Region: CAN (Ontario, Canada) ──
#   Taxa: birds
#   Env vars: 11 (alt, asp, cti, dtr, frs, pre, slp, srad, tmax, tmin, tmp)
#   Species (30):
#     can01, can02, can03, can04, can05, can06, can07, can08, can09, can10,
#     can11, can12, can13, can14, can15, can16, can17, can18, can19, can20,
#     can21, can22, can23, can24, can25, can26, can27, can28, can29, can30
#   Usage:
#     TARGET_REGION  <- "CAN"
#     TARGET_SPECIES <- "can01"
#
# ── Region: NSW (New South Wales, Australia) ──
#   Taxa: vascular plants
#   Env vars: 11 (cti, disturb, mi, mi_run, rainann, raindq, ruession, soild, soiln, tmax, tmin)
#   Species (30):
#     nsw01, nsw02, nsw03, nsw04, nsw05, nsw06, nsw07, nsw08, nsw09, nsw10,
#     nsw11, nsw12, nsw13, nsw14, nsw15, nsw16, nsw17, nsw18, nsw19, nsw20,
#     nsw21, nsw22, nsw23, nsw24, nsw25, nsw26, nsw27, nsw28, nsw29, nsw30
#   Usage:
#     TARGET_REGION  <- "NSW"
#     TARGET_SPECIES <- "nsw01"
#
# ── Region: NZ (New Zealand) ──
#   Taxa: vascular plants
#   Env vars: 14 (age, deficit, dem, gdd, mat, mas, r2pet, rain, slope, solar, toxicite, vpd)
#   Species (52):
#     nz01, nz02, nz03, ..., nz52
#   Usage:
#     TARGET_REGION  <- "NZ"
#     TARGET_SPECIES <- "nz01"
#
# ── Region: SA (South America) ──
#   Taxa: vascular plants
#   Env vars: 6 (bio5, bio6, bio16, bio17, cti, disturession)
#   Species (30):
#     sa01, sa02, sa03, sa04, sa05, sa06, sa07, sa08, sa09, sa10,
#     sa11, sa12, sa13, sa14, sa15, sa16, sa17, sa18, sa19, sa20,
#     sa21, sa22, sa23, sa24, sa25, sa26, sa27, sa28, sa29, sa30
#   Usage:
#     TARGET_REGION  <- "SA"
#     TARGET_SPECIES <- "sa01"
#
# ── Region: SWI (Switzerland) ──
#   Taxa: vascular plants
#   Env vars: 13 (bcc, calc, ccc, ddeg, nutri, pday, precyy, sfroyy, slope, sradyy, swb, tavecc, topo)
#   Species (25 with PO≥200):
#     swi01, swi02, swi03, swi04, swi05, swi06, swi07, swi08, swi09, swi10,
#     swi11, swi12, swi13, swi14, swi15, swi16, swi17, swi18, swi19, swi20,
#     swi21, swi22, swi23, swi24, swi25
#   Usage:
#     TARGET_REGION  <- "SWI"
#     TARGET_SPECIES <- "swi06"
#
# ══════════════════════════════════════════════════════════════════════════════

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

# ╔══════════════════════════════════════════════════════════════════════════╗
# ║  ► CHANGE THESE TWO PARAMETERS ◄                                       ║
# ╚══════════════════════════════════════════════════════════════════════════╝
TARGET_REGION <- "SWI" # One of: "AWT", "CAN", "NSW", "NZ", "SA", "SWI"
TARGET_SPECIES <- "swi06" # Species ID matching the region prefix

# ---- Dependencies ----
pkgs <- c("tidyverse", "bnlearn", "pROC", "caret", "ranger", "maxnet", "gbm", "mgcv", "torch")
for (pkg in pkgs) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
}
if (!torch_is_installed()) torch::install_torch()

dir.create(sprintf("output/case2/%s/single_species", TARGET_REGION),
    recursive = TRUE, showWarnings = FALSE
)

# ==============================================================================
# Utility Functions  (identical to 02_multi_species_experiment.R)
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
            ceiling((n_vars - ceiling(n_vars / n_groups)) / 2), 0
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
# CI-MLP Feature Engineering (identical to 02)
# ==============================================================================
build_cimlp_features <- function(X, cast_vars, adj_matrix, ate_results,
                                 boot_str = NULL) {
    X_base <- X[, cast_vars, drop = FALSE]
    p <- length(cast_vars)

    ate_weights <- rep(1.0, p)
    names(ate_weights) <- cast_vars
    for (v in cast_vars) {
        idx <- which(ate_results$variable == v)
        if (length(idx) > 0 && ate_results$significant[idx[1]]) {
            ate_weights[v] <- 1.0 + abs(ate_results$coef[idx[1]])
        }
    }
    X_weighted <- X_base
    for (v in cast_vars) {
        X_weighted[[v]] <- X_weighted[[v]] * ate_weights[v]
    }

    interaction_cols <- list()
    edge_names <- c()
    for (i in 1:p) {
        for (j in 1:p) {
            if (i != j && adj_matrix[cast_vars[i], cast_vars[j]] > 0) {
                col_name <- paste0("int_", cast_vars[i], "_", cast_vars[j])
                edge_weight <- 1.0
                if (!is.null(boot_str)) {
                    edge_idx <- which(boot_str$from == cast_vars[i] &
                        boot_str$to == cast_vars[j])
                    if (length(edge_idx) > 0) {
                        edge_weight <- boot_str$strength[edge_idx[1]]
                    }
                }
                interaction_cols[[col_name]] <- X_base[[cast_vars[i]]] *
                    X_base[[cast_vars[j]]] * edge_weight
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
# CI-MLP & FlatNN Architectures (identical to 02)
# ==============================================================================
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

FlatNN <- nn_module("FlatNN",
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
# Training & Dataset (identical to 02)
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

flat_dataset <- dataset("FlatDS",
    initialize = function(X, y) {
        self$x <- torch_tensor(as.matrix(X), dtype = torch_float())
        self$y <- torch_tensor(y, dtype = torch_float())$unsqueeze(2)
    },
    .getitem = function(i) list(x = self$x[i, ], y = self$y[i, ]),
    .length = function() self$y$size(1)
)

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
# Single-Species Pipeline
# ==============================================================================
sp <- TARGET_SPECIES
region <- TARGET_REGION
n_runs <- 5
seeds <- c(42, 71, 103, 137, 251)

cat("======================================================================\n")
cat(sprintf("  CAST v3 Single-Species Pipeline: %s [%s]\n", sp, region))
cat("======================================================================\n\n")

# Load data — try region-specific path first, then legacy SWI path
train_file <- sprintf("output/case2/%s/train_data_%s.csv", region, sp)
test_file <- sprintf("output/case2/%s/test_data_%s.csv", region, sp)
if (!file.exists(train_file)) {
    # Fallback for old SWI data layout
    train_file <- sprintf("output/case2/train_data_%s.csv", sp)
    test_file <- sprintf("output/case2/test_data_%s.csv", sp)
}
stopifnot(file.exists(train_file), file.exists(test_file))

train_raw <- read.csv(train_file, stringsAsFactors = FALSE)
test_raw <- read.csv(test_file, stringsAsFactors = FALSE)
env_cols <- setdiff(names(train_raw), "presence")
for (col in env_cols) {
    train_raw[[col]] <- as.numeric(train_raw[[col]])
    test_raw[[col]] <- as.numeric(test_raw[[col]])
}

cat(sprintf(
    "  Train: %d (P=%d, BG=%d, prev=%.3f)\n",
    nrow(train_raw), sum(train_raw$presence == 1), sum(train_raw$presence == 0),
    mean(train_raw$presence == 1)
))
cat(sprintf(
    "  Test:  %d (prev=%.3f)\n",
    nrow(test_raw), mean(test_raw$presence == 1)
))

# ═══ Step 1: VIF ═══
cat("\n═══ Step 1: VIF Filtering ═══\n")
X_env <- train_raw[, env_cols, drop = FALSE]
sds <- apply(X_env, 2, sd, na.rm = TRUE)
X_env <- X_env[, sds > 1e-10, drop = FALSE]
n_start <- ncol(X_env)
repeat {
    vifs <- calc_vif(X_env)
    if (max(vifs, na.rm = TRUE) <= 10 || ncol(X_env) <= 3) break
    X_env <- X_env[, setdiff(names(X_env), names(which.max(vifs))), drop = FALSE]
}
selected_vars <- names(X_env)
cat(sprintf("  %d → %d variables (threshold=10)\n", n_start, length(selected_vars)))

train_data <- cbind(presence = train_raw$presence, train_raw[, selected_vars, drop = FALSE])
test_data <- cbind(presence = test_raw$presence, test_raw[, selected_vars, drop = FALSE])
train_data <- train_data[complete.cases(train_data), ]
test_data <- test_data[complete.cases(test_data), ]

# ═══ Step 2: DAG ═══
cat("\n═══ Step 2: DAG Learning (HC, R=200) ═══\n")
env_for_dag <- train_data[, selected_vars, drop = FALSE]
if (nrow(env_for_dag) > 8000) {
    set.seed(42)
    env_for_dag <- env_for_dag[sample(nrow(env_for_dag), 8000), ]
}
set.seed(42)
boot_str <- bnlearn::boot.strength(env_for_dag,
    R = 200, algorithm = "hc",
    algorithm.args = list(score = "bic-g")
)
strong_edges <- boot_str %>% filter(strength >= 0.7, direction >= 0.6)
n_possible <- length(selected_vars) * (length(selected_vars) - 1) / 2
dag_density <- nrow(strong_edges) / max(n_possible, 1)
node_outdeg <- strong_edges %>%
    group_by(from) %>%
    summarise(out_degree = n(), .groups = "drop")
cat(sprintf("  %d strong edges (density=%.2f)\n", nrow(strong_edges), dag_density))

# ═══ Step 3: ATE ═══
cat("\n═══ Step 3: ATE Estimation (DML, 2-fold CF) ═══\n")
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
            cat(sprintf(
                "  %-8s ATE=%+.4f SE=%.4f p=%.4f %s\n",
                v, res$ate, res$se, res$p_value, ifelse(res$significant, " *", "")
            ))
        },
        error = function(e) cat(sprintf("  %-8s FAILED: %s\n", v, e$message))
    )
}
n_sig <- sum(ate_results$significant)
cat(sprintf("  %d/%d significant (p<0.05)\n", n_sig, nrow(ate_results)))

# ═══ Step 4: Adaptive CAST Screening v2 ═══
cat("\n═══ Step 4: Adaptive CAST Screening v2 ═══\n")
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
cat(sprintf(
    "  Adaptive weights: w_dag=%.3f w_ate=%.3f w_imp=%.3f\n",
    w_dag, w_ate, w_imp
))

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

min_keep_n <- max(5L, ceiling(length(selected_vars) * 0.4))
if (length(unique(screening_df$score_total)) >= 2) {
    km <- kmeans(screening_df$score_total, centers = 2, nstart = 10)
    high_cluster <- which.max(km$centers)
    cast_vars_km <- screening_df$variable[km$cluster == high_cluster]
    threshold_method <- sprintf("kmeans(high=%d)", length(cast_vars_km))
} else {
    cast_vars_km <- screening_df$variable
    threshold_method <- "all-equal"
}
if (length(cast_vars_km) < min_keep_n) {
    cast_vars <- screening_df$variable[1:min(min_keep_n, nrow(screening_df))]
    threshold_method <- paste0(threshold_method, "+floor(", min_keep_n, ")")
} else {
    cast_vars <- cast_vars_km
}

cat(sprintf(
    "  %d → %d CAST variables (%s)\n",
    length(selected_vars), length(cast_vars), threshold_method
))

cat("  Screening scores:\n")
print(screening_df %>% select(variable, score_dag, score_ate, score_imp, score_total))

write.csv(screening_df, sprintf(
    "output/case2/%s/single_species/screening_%s.csv",
    region, sp
), row.names = FALSE)
write.csv(ate_results, sprintf(
    "output/case2/%s/single_species/ate_%s.csv",
    region, sp
), row.names = FALSE)

# ═══ Step 5: Causal Role Grouping ═══
cat("\n═══ Step 5: Causal Role Grouping ═══\n")
role_df <- assign_causal_roles(cast_vars, strong_edges, n_groups = 3)
adj_matrix <- build_adj_matrix(cast_vars, strong_edges)
n_dag_edges_cast <- sum(adj_matrix)
hidden_size <- max(32L, min(128L, as.integer(length(cast_vars) * 8)))

print(role_df %>% select(variable, out, inp, role_score, group))
cat(sprintf("  DAG edges among CAST vars: %d\n", n_dag_edges_cast))
cat(sprintf("  Hidden size: %d\n", hidden_size))

write.csv(role_df, sprintf(
    "output/case2/%s/single_species/roles_%s.csv",
    region, sp
), row.names = FALSE)

# ═══ Step 6: Model Training ═══
cat("\n═══ Step 6: Model Training & Evaluation ═══\n")
cat(sprintf("  A: Full-var baselines (%d vars)\n", length(selected_vars)))
cat(sprintf("  B: CAST-screened baselines (%d vars)\n", length(cast_vars)))
cat(sprintf(
    "  C: CI-MLP (%d vars + %d DAG interaction features)\n",
    length(cast_vars), n_dag_edges_cast
))

y_train_all <- train_data$presence
y_test_all <- test_data$presence

# Scale data
X_train_cast <- train_data[, cast_vars, drop = FALSE]
X_test_cast <- test_data[, cast_vars, drop = FALSE]
X_means_cast <- colMeans(X_train_cast, na.rm = TRUE)
X_sds_cast <- apply(X_train_cast, 2, sd, na.rm = TRUE)
X_sds_cast[X_sds_cast < 1e-10] <- 1
X_train_cast_sc <- as.data.frame(scale(X_train_cast, center = X_means_cast, scale = X_sds_cast))
X_test_cast_sc <- as.data.frame(scale(X_test_cast, center = X_means_cast, scale = X_sds_cast))
X_train_cast_sc[is.na(X_train_cast_sc)] <- 0
X_test_cast_sc[is.na(X_test_cast_sc)] <- 0

X_train_full <- train_data[, selected_vars, drop = FALSE]
X_test_full <- test_data[, selected_vars, drop = FALSE]
X_means_full <- colMeans(X_train_full, na.rm = TRUE)
X_sds_full <- apply(X_train_full, 2, sd, na.rm = TRUE)
X_sds_full[X_sds_full < 1e-10] <- 1
X_train_full_sc <- as.data.frame(scale(X_train_full, center = X_means_full, scale = X_sds_full))
X_test_full_sc <- as.data.frame(scale(X_test_full, center = X_means_full, scale = X_sds_full))
X_train_full_sc[is.na(X_train_full_sc)] <- 0
X_test_full_sc[is.na(X_test_full_sc)] <- 0

# CI-MLP features
ci_train_info <- build_cimlp_features(
    X_train_cast_sc, cast_vars,
    adj_matrix, ate_results, boot_str
)
ci_test_info <- build_cimlp_features(
    X_test_cast_sc, cast_vars,
    adj_matrix, ate_results, boot_str
)
X_train_ci <- ci_train_info$data
X_test_ci <- ci_test_info$data

cat(sprintf(
    "  CI-MLP features: %d base + %d interactions = %d total\n",
    ci_train_info$n_base, ci_train_info$n_interactions, ci_train_info$n_total
))
if (ci_train_info$n_interactions > 0) {
    cat(sprintf(
        "  Interactions: %s\n",
        paste(ci_train_info$interaction_names, collapse = ", ")
    ))
}

# Val split
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
hidden_size_full <- max(32L, min(128L, as.integer(length(selected_vars) * 8)))

X_tr_cast <- X_train_cast_sc[-val_idx, ]
X_val_cast <- X_train_cast_sc[val_idx, ]
X_tr_full <- X_train_full_sc[-val_idx, ]
X_val_full <- X_train_full_sc[val_idx, ]
X_tr_ci <- X_train_ci[-val_idx, ]
X_val_ci <- X_train_ci[val_idx, ]

sp_results <- data.frame()

# ---- Group C: CI-MLP ----
cat("\n  --- Group C: CI-MLP (CAST vars + DAG interactions) ---\n")
{
    aucs <- numeric(n_runs)
    tsss <- numeric(n_runs)
    for (ri in 1:n_runs) {
        torch_manual_seed(seeds[ri])
        set.seed(seeds[ri])
        tryCatch(
            {
                ds <- flat_dataset(X_tr_ci, y_tr)
                dl <- dataloader(ds, batch_size = batch_size, shuffle = TRUE, drop_last = TRUE)
                vt <- torch_tensor(as.matrix(X_val_ci), dtype = torch_float())
                m <- CI_MLP(ncol(X_tr_ci), hidden_size, 0.2)
                res <- train_nn(m, dl,
                    function(m) as.numeric(torch_sigmoid(m(vt))$squeeze()$cpu()),
                    y_val,
                    epochs = 200, patience = 30, focal_alpha = focal_alpha
                )
                res$model$eval()
                with_no_grad({
                    tt <- torch_tensor(as.matrix(X_test_ci), dtype = torch_float())
                    pred <- as.numeric(torch_sigmoid(res$model(tt))$squeeze()$cpu())
                })
                ev <- evaluate_model(pred, y_test_all)
                aucs[ri] <- ev["auc"]
                tsss[ri] <- ev["tss"]
                cat(sprintf("    CI-MLP run %d: AUC=%.4f TSS=%.4f\n", ri, ev["auc"], ev["tss"]))
            },
            error = function(e) {
                cat(sprintf("    CI-MLP run %d FAILED: %s\n", ri, e$message))
                aucs[ri] <<- NA
                tsss[ri] <<- NA
            }
        )
    }
    sp_results <- rbind(sp_results, data.frame(
        region = region, species = sp, model = "CI_MLP", var_set = "cast",
        n_vars = length(cast_vars), n_interactions = ci_train_info$n_interactions,
        n_features_total = ci_train_info$n_total,
        n_dag_edges = n_dag_edges_cast, dag_density = dag_density,
        auc_mean = mean(aucs, na.rm = TRUE), auc_sd = sd(aucs, na.rm = TRUE),
        tss_mean = mean(tsss, na.rm = TRUE), tss_sd = sd(tsss, na.rm = TRUE),
        n_success = sum(!is.na(aucs)), stringsAsFactors = FALSE
    ))
}

# ---- Group B: FlatNN_cast ----
cat("\n  --- Group B: FlatNN_cast (CAST vars, no structure) ---\n")
{
    aucs <- numeric(n_runs)
    tsss <- numeric(n_runs)
    for (ri in 1:n_runs) {
        torch_manual_seed(seeds[ri])
        set.seed(seeds[ri])
        tryCatch(
            {
                ds <- flat_dataset(X_tr_cast, y_tr)
                dl <- dataloader(ds, batch_size = batch_size, shuffle = TRUE, drop_last = TRUE)
                vt <- torch_tensor(as.matrix(X_val_cast), dtype = torch_float())
                m <- FlatNN(length(cast_vars), hidden_size, 0.2)
                res <- train_nn(m, dl,
                    function(m) as.numeric(torch_sigmoid(m(vt))$squeeze()$cpu()),
                    y_val,
                    epochs = 200, patience = 30, focal_alpha = focal_alpha
                )
                res$model$eval()
                with_no_grad({
                    tt <- torch_tensor(as.matrix(X_test_cast_sc), dtype = torch_float())
                    pred <- as.numeric(torch_sigmoid(res$model(tt))$squeeze()$cpu())
                })
                ev <- evaluate_model(pred, y_test_all)
                aucs[ri] <- ev["auc"]
                tsss[ri] <- ev["tss"]
                cat(sprintf("    FlatNN_cast run %d: AUC=%.4f TSS=%.4f\n", ri, ev["auc"], ev["tss"]))
            },
            error = function(e) {
                cat(sprintf("    FlatNN_cast run %d FAILED: %s\n", ri, e$message))
                aucs[ri] <<- NA
                tsss[ri] <<- NA
            }
        )
    }
    sp_results <- rbind(sp_results, data.frame(
        region = region, species = sp, model = "FlatNN_cast", var_set = "cast",
        n_vars = length(cast_vars), n_interactions = 0L,
        n_features_total = length(cast_vars),
        n_dag_edges = n_dag_edges_cast, dag_density = dag_density,
        auc_mean = mean(aucs, na.rm = TRUE), auc_sd = sd(aucs, na.rm = TRUE),
        tss_mean = mean(tsss, na.rm = TRUE), tss_sd = sd(tsss, na.rm = TRUE),
        n_success = sum(!is.na(aucs)), stringsAsFactors = FALSE
    ))
}

# ---- Group A: FlatNN_full ----
cat("\n  --- Group A: FlatNN_full (all post-VIF vars) ---\n")
{
    aucs <- numeric(n_runs)
    tsss <- numeric(n_runs)
    for (ri in 1:n_runs) {
        torch_manual_seed(seeds[ri])
        set.seed(seeds[ri])
        tryCatch(
            {
                ds <- flat_dataset(X_tr_full, y_tr)
                dl <- dataloader(ds, batch_size = batch_size, shuffle = TRUE, drop_last = TRUE)
                vt <- torch_tensor(as.matrix(X_val_full), dtype = torch_float())
                m <- FlatNN(length(selected_vars), hidden_size_full, 0.2)
                res <- train_nn(m, dl,
                    function(m) as.numeric(torch_sigmoid(m(vt))$squeeze()$cpu()),
                    y_val,
                    epochs = 200, patience = 30, focal_alpha = focal_alpha
                )
                res$model$eval()
                with_no_grad({
                    tt <- torch_tensor(as.matrix(X_test_full_sc), dtype = torch_float())
                    pred <- as.numeric(torch_sigmoid(res$model(tt))$squeeze()$cpu())
                })
                ev <- evaluate_model(pred, y_test_all)
                aucs[ri] <- ev["auc"]
                tsss[ri] <- ev["tss"]
                cat(sprintf("    FlatNN_full run %d: AUC=%.4f TSS=%.4f\n", ri, ev["auc"], ev["tss"]))
            },
            error = function(e) {
                cat(sprintf("    FlatNN_full run %d FAILED: %s\n", ri, e$message))
                aucs[ri] <<- NA
                tsss[ri] <<- NA
            }
        )
    }
    sp_results <- rbind(sp_results, data.frame(
        region = region, species = sp, model = "FlatNN_full", var_set = "full",
        n_vars = length(selected_vars), n_interactions = 0L,
        n_features_total = length(selected_vars),
        n_dag_edges = NA, dag_density = dag_density,
        auc_mean = mean(aucs, na.rm = TRUE), auc_sd = sd(aucs, na.rm = TRUE),
        tss_mean = mean(tsss, na.rm = TRUE), tss_sd = sd(tsss, na.rm = TRUE),
        n_success = sum(!is.na(aucs)), stringsAsFactors = FALSE
    ))
}

# ---- Traditional SDMs ----
cat("\n  --- Traditional SDM Baselines ---\n")
y_tr_raw <- train_data$presence
for (sdm_name in c("RF", "Maxent", "BRT", "GAM")) {
    for (vset in c("cast", "full")) {
        vars_use <- if (vset == "cast") cast_vars else selected_vars
        model_label <- paste0(sdm_name, "_", vset)
        tryCatch(
            {
                pred <- train_sdm(sdm_name, train_data[, vars_use, drop = FALSE],
                    y_tr_raw, test_data[, vars_use, drop = FALSE],
                    gam_vars = vars_use
                )
                ev <- evaluate_model(pred, y_test_all)
                sp_results <- rbind(sp_results, data.frame(
                    region = region, species = sp, model = model_label,
                    var_set = vset, n_vars = length(vars_use),
                    n_interactions = 0L, n_features_total = length(vars_use),
                    n_dag_edges = NA, dag_density = dag_density,
                    auc_mean = ev["auc"], auc_sd = 0,
                    tss_mean = ev["tss"], tss_sd = 0,
                    n_success = 1, stringsAsFactors = FALSE
                ))
                cat(sprintf("    %s: AUC=%.4f\n", model_label, ev["auc"]))
            },
            error = function(e) {
                cat(sprintf("    %s FAILED: %s\n", model_label, substr(e$message, 1, 60)))
            }
        )
    }
}

# Save results
write.csv(sp_results, sprintf(
    "output/case2/%s/single_species/results_%s.csv",
    region, sp
), row.names = FALSE)

# Also save to legacy path for SWI backward compat
if (region == "SWI") {
    dir.create("output/case2/single_species", recursive = TRUE, showWarnings = FALSE)
    write.csv(sp_results, sprintf("output/case2/single_species/results_%s.csv", sp),
        row.names = FALSE
    )
    write.csv(screening_df, sprintf("output/case2/single_species/screening_%s.csv", sp),
        row.names = FALSE
    )
    write.csv(ate_results, sprintf("output/case2/single_species/ate_%s.csv", sp),
        row.names = FALSE
    )
    write.csv(role_df, sprintf("output/case2/single_species/roles_%s.csv", sp),
        row.names = FALSE
    )
}

# ═══ Summary ═══
cat("\n╔═══════════════════════════════════════════════════════════════════════════╗\n")
cat(sprintf(
    "║  CAST v3 Pipeline Complete: %s [%s]                                    ║\n",
    sp, region
))
cat("╠═══════════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf(
    "║  VIF:       %d → %d vars                                                ║\n",
    n_start, length(selected_vars)
))
cat(sprintf(
    "║  DAG:       %d edges (density=%.2f)                                     ║\n",
    nrow(strong_edges), dag_density
))
cat(sprintf(
    "║  ATE:       %d/%d significant                                            ║\n",
    n_sig, nrow(ate_results)
))
cat(sprintf(
    "║  Screening: %d → %d CAST vars                                           ║\n",
    length(selected_vars), length(cast_vars)
))
cat(sprintf(
    "║  CI-MLP:    %d base + %d interactions = %d features                      ║\n",
    ci_train_info$n_base, ci_train_info$n_interactions, ci_train_info$n_total
))
cat("╠═══════════════════════════════════════════════════════════════════════════╣\n")

for (grp in c("full", "cast")) {
    label <- if (grp == "full") "A: Full-variable baselines" else "B: CAST-screened baselines"
    cat(sprintf("║  %s%s║\n", label, strrep(" ", 75 - 4 - nchar(label))))
    subset <- sp_results %>%
        filter(var_set == grp, model != "CI_MLP") %>%
        arrange(desc(auc_mean))
    for (i in 1:nrow(subset)) {
        r <- subset[i, ]
        cat(sprintf(
            "║    %-14s AUC=%.4f±%.4f  TSS=%.4f±%.4f  (%dv)             ║\n",
            r$model, r$auc_mean, r$auc_sd, r$tss_mean, r$tss_sd, r$n_vars
        ))
    }
}
cat("║  C: CI-MLP (CAST + DAG interactions)                                      ║\n")
ci_row <- sp_results %>% filter(model == "CI_MLP")
cat(sprintf(
    "║  ★ %-14s AUC=%.4f±%.4f  TSS=%.4f±%.4f  (%dv+%d)           ║\n",
    ci_row$model, ci_row$auc_mean, ci_row$auc_sd,
    ci_row$tss_mean, ci_row$tss_sd, ci_row$n_vars, ci_row$n_interactions
))

# Key comparison
flatnn_cast_auc <- (sp_results %>% filter(model == "FlatNN_cast"))$auc_mean
ci_auc <- ci_row$auc_mean
delta <- ci_auc - flatnn_cast_auc
cat("║                                                                           ║\n")
cat(sprintf(
    "║  Structure effect: CI-MLP vs FlatNN_cast = %+.4f AUC%s            ║\n",
    delta, ifelse(delta > 0, " ✓", "")
))
cat("╚═══════════════════════════════════════════════════════════════════════════╝\n")
