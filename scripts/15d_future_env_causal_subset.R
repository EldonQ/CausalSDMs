#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 15d_future_env_causal_subset.R
# 功能说明: 在 19+4 未来变量体系下，使用因果核心变量 ∩ (bioc01-19 + 4 静态)
#          的交集重训 4 个 SDM，并进行未来情景投影，与 full 19+4 模型对比
# 输入:  output/15b_causal_retraining/core_drivers_selection.csv
#       output/15_future_env_19plus4/training_data_19plus4.csv
#       output/15_future_env_19plus4/future_bioc_china_SSP*.tif
# 输出: output/15_future_env_19plus4/models_causal_19plus4/*.rds
#       output/15_future_env_19plus4/evaluation_causal_19plus4.csv
#       output/15_future_env_19plus4/rasters/*_causal_19plus4_*.tif
#       figures/15_future_env_19plus4/*_causal_19plus4_*.png
# ============================================================================

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

packages <- c(
  "raster", "sf", "tidyverse", "maxnet", "randomForest", "mgcv",
  "nnet", "pROC", "sysfonts", "showtext", "ggplot2", "viridis"
)
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

try(
  {
    sysfonts::font_add(
      family = "Arial",
      regular = "C:/Windows/Fonts/arial.ttf",
      bold = "C:/Windows/Fonts/arialbd.ttf",
      italic = "C:/Windows/Fonts/ariali.ttf",
      bolditalic = "C:/Windows/Fonts/arialbi.ttf"
    )
    showtext::showtext_opts(dpi = 2400)
    showtext::showtext_auto(enable = TRUE)
  },
  silent = TRUE
)

source("scripts/visualization/viz_utils.R")

dir.create("output/15_future_env_19plus4/models_causal_19plus4", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("因果 ∩ 19+4 变量的未来预测\n")
cat("======================================\n\n")

ssp_scenarios <- c("SSP126", "SSP245", "SSP370", "SSP585")

# -----------------------------
# 步骤1: 读取因果核心变量与 19+4 训练数据，求交集
# -----------------------------

cat("步骤 1/4: 读取因果核心变量与 19+4 训练数据...\n")

core_path <- "output/15b_causal_retraining/core_drivers_selection.csv"
train_19p4_path <- "output/15_future_env_19plus4/training_data_19plus4.csv"

if (!file.exists(core_path)) stop("未找到因果核心变量列表: ", core_path)
if (!file.exists(train_19p4_path)) stop("未找到 19+4 训练数据: ", train_19p4_path)

core_sel <- read.csv(core_path, stringsAsFactors = FALSE)
core_drivers <- unique(core_sel$variable)

train_df <- read.csv(train_19p4_path)

future_var_set <- setdiff(colnames(train_df), c("id", "species", "lon", "lat", "presence"))

# 静态变量直接按名称相交
static_candidates <- c("dem_avg", "slope_avg", "geology_total", "soil_avg_01")
static_vars <- intersect(core_drivers, intersect(static_candidates, future_var_set))

# 气候变量：将 core_drivers 中的 hydro_wavg_XX 按编号映射到 bioXX
hydro_vars <- core_drivers[grepl("^hydro_wavg_0*[0-9]+$", core_drivers)]
hydro_idx <- suppressWarnings(as.integer(sub("^hydro_wavg_0*([0-9]+)$", "\\1", hydro_vars)))
hydro_idx <- hydro_idx[!is.na(hydro_idx)]
mapped_bio <- sprintf("bio%02d", hydro_idx)
mapped_bio <- mapped_bio[mapped_bio %in% future_var_set]

vars_causal_future <- unique(c(static_vars, mapped_bio))

if (length(vars_causal_future) < 3) {
  warning("因果核心 ∩ 19+4（映射后）变量数少于 3 个，结果仅供探索性参考")
}

cat("  ✓ 因果核心变量数: ", length(core_drivers), "\n", sep = "")
cat("  ✓ 19+4 变量数: ", length(future_var_set), "\n", sep = "")
cat("  ✓ 因果映射后用于未来预测的变量数: ", length(vars_causal_future), "\n", sep = "")
if (length(static_vars) > 0) {
  cat("    静态变量 (来自因果核心 ∩ 19+4): ", paste(static_vars, collapse = ", "), "\n", sep = "")
}
if (length(mapped_bio) > 0) {
  cat("    气候变量 (由 hydro_wavg_XX → bioXX 映射): ", paste(mapped_bio, collapse = ", "), "\n", sep = "")
}
if (length(vars_causal_future) > 0) {
  cat("    最终用于简化模型的变量列表:\n")
  for (i in seq_along(vars_causal_future)) {
    cat(sprintf("    %2d. %s\n", i, vars_causal_future[i]))
  }
}

# -----------------------------
# 步骤2: 在交集变量集上重训 4 个 SDM
# -----------------------------

cat("\n步骤 2/4: 在交集变量集上重训 4 个 SDM...\n")

train_df <- train_df[stats::complete.cases(train_df[, c("presence", vars_causal_future)]), ]
set.seed(20251121)
pres_idx <- which(train_df$presence == 1)
back_idx <- which(train_df$presence == 0)
train_idx <- c(
  sample(pres_idx, round(0.8 * length(pres_idx))),
  sample(back_idx, round(0.8 * length(back_idx)))
)
test_idx <- setdiff(seq_len(nrow(train_df)), train_idx)

X_train <- train_df[train_idx, vars_causal_future, drop = FALSE]
X_test <- train_df[test_idx, vars_causal_future, drop = FALSE]
y_train <- train_df$presence[train_idx]
y_test <- train_df$presence[test_idx]

models_causal <- list()

# Maxnet
cat("  - Maxnet (causal 19+4 交集)...\n")
mx_model <- maxnet::maxnet(p = y_train, data = X_train, f = maxnet::maxnet.formula(y_train, X_train))
models_causal$Maxnet <- mx_model

# RF
cat("  - RF (causal 19+4 交集)...\n")
rf_model <- randomForest::randomForest(x = X_train, y = factor(y_train), ntree = 500)
models_causal$RF <- rf_model

# GAM
cat("  - GAM (causal 19+4 交集)...\n")
form_terms <- paste0("s(", vars_causal_future, ", k=5)", collapse = " + ")
form <- as.formula(paste0("presence ~ ", form_terms, " + s(lon,lat, k=10)"))
gam_model <- mgcv::gam(form, data = cbind(train_df[train_idx, ], X_train), family = binomial(link = "logit"))
models_causal$GAM <- gam_model

# NN
cat("  - NN (causal 19+4 交集)...\n")
mu <- sapply(X_train, mean, na.rm = TRUE)
sdv <- sapply(X_train, sd, na.rm = TRUE)
sdv[sdv == 0 | is.na(sdv)] <- 1
X_train_s <- as.data.frame(sweep(sweep(as.matrix(X_train), 2, mu, "-"), 2, sdv, "/"))
X_test_s <- as.data.frame(sweep(sweep(as.matrix(X_test), 2, mu, "-"), 2, sdv, "/"))
size_hidden <- max(3, floor(length(vars_causal_future) / 2))

tmp_nn <- nnet::nnet(
  x = X_train_s, y = y_train, size = size_hidden, linout = FALSE,
  rang = 0.1, decay = 5e-4, maxit = 500, trace = FALSE
)
models_causal$NN <- list(model = tmp_nn, mean = mu, sd = sdv, vars = vars_causal_future)

saveRDS(models_causal$Maxnet, "output/15_future_env_19plus4/models_causal_19plus4/maxnet_causal_19plus4.rds")
saveRDS(models_causal$RF, "output/15_future_env_19plus4/models_causal_19plus4/rf_causal_19plus4.rds")
saveRDS(models_causal$GAM, "output/15_future_env_19plus4/models_causal_19plus4/gam_causal_19plus4.rds")
saveRDS(models_causal$NN, "output/15_future_env_19plus4/models_causal_19plus4/nn_causal_19plus4.rds")

# 评估
cat("\n  - 评估简化模型性能...\n")

eval_list <- list()

# Maxnet
pred_mx <- as.numeric(predict(models_causal$Maxnet, X_test, type = "logistic"))
auc_mx <- as.numeric(pROC::auc(pROC::roc(y_test, pred_mx, quiet = TRUE)))

# RF
pred_rf <- as.numeric(predict(models_causal$RF, newdata = X_test, type = "prob")[, "1"])
auc_rf <- as.numeric(pROC::auc(pROC::roc(y_test, pred_rf, quiet = TRUE)))

# GAM
pred_gam <- as.numeric(predict(models_causal$GAM, newdata = cbind(train_df[test_idx, ], X_test), type = "response"))
auc_gam <- as.numeric(pROC::auc(pROC::roc(y_test, pred_gam, quiet = TRUE)))

# NN
pred_nn <- as.numeric(nnet:::predict.nnet(models_causal$NN$model, as.matrix(X_test_s), type = "raw"))
auc_nn <- as.numeric(pROC::auc(pROC::roc(y_test, pred_nn, quiet = TRUE)))

eval_all <- dplyr::bind_rows(
  data.frame(model = "Maxnet", AUC = auc_mx),
  data.frame(model = "RF", AUC = auc_rf),
  data.frame(model = "GAM", AUC = auc_gam),
  data.frame(model = "NN", AUC = auc_nn)
)
write.csv(eval_all, "output/15_future_env_19plus4/evaluation_causal_19plus4.csv", row.names = FALSE)

cat("  ✓ 简化模型训练与评估完成\n\n")

# -----------------------------
# 步骤3: 在 19+4 未来环境上进行因果子集投影
# -----------------------------

cat("步骤 3/4: 在 19+4 未来环境上进行因果子集投影...\n")

future_bioc_list <- list()
for (ssp in ssp_scenarios) {
  tif_path <- file.path("output/15_future_env_19plus4", paste0("future_bioc_china_", ssp, ".tif"))
  if (!file.exists(tif_path)) {
    warning("未找到未来 bioc tif: ", tif_path)
    next
  }
  future_bioc_list[[ssp]] <- raster::brick(tif_path)
}

fa <- raster::brick("earthenvstreams_china/flow_acc.tif")[[2]]
fa_vals <- raster::getValues(fa)
fa_vals[fa_vals <= 0] <- NA
river_mask <- raster::setValues(fa, fa_vals)
rm(fa_vals)

# hydro_ord2_6 <- viz_load_hydrorivers_ord2_6(crs_target = raster::crs(fa)) # Removed

make_predict_fun <- function(model_name, model_obj) {
  if (model_name == "Maxnet") {
    return(function(m, df) {
      as.numeric(predict(m, df, type = "logistic"))
    })
  }
  if (model_name == "RF") {
    return(function(m, df) {
      as.numeric(predict(m, newdata = df, type = "prob")[, "1"])
    })
  }
  if (model_name == "GAM") {
    return(function(m, df) {
      as.numeric(predict(m, newdata = df, type = "response"))
    })
  }
  if (model_name == "NN") {
    return(function(m, df) {
      mu <- m$mean
      sdv <- m$sd
      mod <- m$model
      vars <- m$vars
      sdv[sdv == 0 | is.na(sdv)] <- 1
      x <- as.matrix(df[, vars, drop = FALSE])
      x <- sweep(x, 2, mu[vars], "-")
      x <- sweep(x, 2, sdv[vars], "/")
      as.numeric(nnet:::predict.nnet(mod, x, type = "raw"))
    })
  }
}

build_future_env_causal <- function(bioc_stack) {
  names(bioc_stack) <- sprintf("bio%02d", 1:raster::nlayers(bioc_stack))
  base_ref <- bioc_stack[["bio01"]]
  static_list <- list()
  if (file.exists("earthenvstreams_china/elevation.tif")) {
    dem_br <- raster::brick("earthenvstreams_china/elevation.tif")[[4]]
    dem_br <- suppressWarnings(raster::projectRaster(dem_br, base_ref, method = "bilinear"))
    names(dem_br) <- "dem_avg"
    static_list[[length(static_list) + 1]] <- dem_br
  }
  if (file.exists("earthenvstreams_china/slope.tif")) {
    slope_br <- raster::brick("earthenvstreams_china/slope.tif")[[4]] / 100
    slope_br <- suppressWarnings(raster::projectRaster(slope_br, base_ref, method = "bilinear"))
    names(slope_br) <- "slope_avg"
    static_list[[length(static_list) + 1]] <- slope_br
  }
  if (file.exists("earthenvstreams_china/geology_weighted_sum.tif")) {
    geo_total <- raster::calc(raster::brick("earthenvstreams_china/geology_weighted_sum.tif"), fun = function(x) {
      if (all(is.na(x))) NA else sum(x, na.rm = TRUE)
    })
    geo_total <- suppressWarnings(raster::projectRaster(geo_total, base_ref, method = "bilinear"))
    names(geo_total) <- "geology_total"
    static_list[[length(static_list) + 1]] <- geo_total
  }
  if (file.exists("earthenvstreams_china/soil_average.tif")) {
    soil_soc <- raster::brick("earthenvstreams_china/soil_average.tif")[[1]]
    soil_soc <- suppressWarnings(raster::projectRaster(soil_soc, base_ref, method = "bilinear"))
    names(soil_soc) <- "soil_avg_01"
    static_list[[length(static_list) + 1]] <- soil_soc
  }
  stk <- raster::stack(c(bioc_stack[[1:19]], static_list))
  lon_r <- raster::init(stk[[1]], fun = "x")
  names(lon_r) <- "lon"
  lat_r <- raster::init(stk[[1]], fun = "y")
  names(lat_r) <- "lat"
  stk <- raster::addLayer(stk, lon_r, lat_r)
  return(stk)
}

all_summaries <- list()

for (ssp in names(future_bioc_list)) {
  bioc_china <- future_bioc_list[[ssp]]
  env_full <- build_future_env_causal(bioc_china)
  # 仅保留交集变量 + lon/lat
  vars_use <- unique(c(vars_causal_future, "lon", "lat"))
  env_stk <- env_full[[vars_use]]

  out_dir_ras <- file.path("output/15_future_env_19plus4/rasters", paste0(ssp, "_causal"))
  out_dir_fig <- file.path("figures/15_future_env_19plus4", paste0(ssp, "_causal"))
  dir.create(out_dir_ras, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)

  summary_rows <- list()
  for (mn in names(models_causal)) {
    cat("  -> 情景 ", ssp, " | 模型 ", mn, " (causal) ...\n", sep = "")
    mdl <- models_causal[[mn]]
    pred_fun <- make_predict_fun(mn, mdl)
    tif_path <- file.path(out_dir_ras, paste0("pred_", tolower(mn), "_causal_19plus4.tif"))
    if (file.exists(tif_path)) {
      try(
        {
          file.remove(tif_path)
        },
        silent = TRUE
      )
    }
    pred_r <- raster::predict(env_stk,
      model = mdl, fun = pred_fun, filename = tif_path,
      overwrite = TRUE, progress = "text"
    )
    pred_r <- raster::clamp(pred_r, lower = 0, upper = 1, useValues = TRUE)
    river_mask_ref <- suppressWarnings(raster::projectRaster(river_mask, pred_r, method = "ngb"))
    pred_r_river <- raster::mask(pred_r, river_mask_ref)
    tif_mask_path <- file.path(out_dir_ras, paste0("pred_", tolower(mn), "_causal_19plus4_river.tif"))
    if (file.exists(tif_mask_path)) {
      try(
        {
          file.remove(tif_mask_path)
        },
        silent = TRUE
      )
    }
    raster::writeRaster(pred_r_river, tif_mask_path, overwrite = TRUE)

    vals <- raster::getValues(pred_r_river)
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        scenario = ssp, model = mn, n_pixels_river = length(vals),
        mean = mean(vals), sd = sd(vals), min = min(vals), max = max(vals),
        p10 = as.numeric(quantile(vals, 0.1)),
        p50 = as.numeric(quantile(vals, 0.5)),
        p90 = as.numeric(quantile(vals, 0.9))
      )
    }

    out_base <- file.path(out_dir_fig, paste0("prediction_", tolower(mn), "_causal_19plus4"))
    # pred_r_river_thick <- viz_thicken_river_raster(pred_r_river, w_size = 3) # Removed
    viz_save_raster_map(
      r = pred_r_river, out_base = out_base,
      title = paste0(ssp, " - ", mn, " (causal 19+4 subset)"),
      palette = "magma", q_limits = c(0.01, 0.99),
      china_path = "earthenvstreams_china/china_boundary.shp",
      width_in = 8, height_in = 6
      # hydrorivers_sf = hydro_ord2_6 # Removed
    )
  }

  if (length(summary_rows) > 0) {
    summary_df <- dplyr::bind_rows(summary_rows)
    write.csv(summary_df, file.path(out_dir_ras, "prediction_summary_causal_19plus4.csv"), row.names = FALSE)
    summary_df$scenario <- factor(summary_df$scenario, levels = c("SSP126", "SSP245", "SSP370", "SSP585"))
    all_summaries[[ssp]] <- summary_df
  }
}

if (length(all_summaries) > 0) {
  trends <- dplyr::bind_rows(all_summaries)
  out_trend_csv <- "output/15_future_env_19plus4/prediction_trends_all_models_causal_19plus4.csv"
  out_trend_png <- "figures/15_future_env_19plus4/habitat_trends_all_models_causal_19plus4.png"
  write.csv(trends, out_trend_csv, row.names = FALSE)
  p_trend <- ggplot(trends, aes(x = scenario, y = mean, group = model, color = model)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.4) +
    labs(
      title = "Habitat Suitability Trend (Causal 19+4 Subset)",
      x = "Scenario", y = "Mean Predicted Suitability", color = "Model"
    ) +
    viz_theme_nature(base_size = 8, title_size = 9)
  ggsave(out_trend_png, p_trend, width = 4.8, height = 3.2, dpi = 2400, bg = "transparent")
}

cat("\n======================================\n")
cat("因果 ∩ 19+4 未来预测完成\n")
cat("======================================\n\n")
