#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 15c_future_env_19plus4.R
# 功能说明: 使用 WorldClim 当前 bioc01-19 + 4 个静态变量重训 4 个 SDM 模型，
#          并在 CMIP6 SSP 情景下基于 19+4 变量进行未来河网适生度预测
# 变量集合: bioc01-bioc19 (WorldClim 当前 + 未来), dem_avg, slope_avg,
#          geology_total, soil_avg_01（静态 earthenvstreams_china）
# 输出目录: output/15_future_env_19plus4
#          figures/15_future_env_19plus4
# 作者: Nature级别科研项目
# 日期: 2025-11-21
# ============================================================================

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

packages <- c(
  "raster", "terra", "sf", "tidyverse", "maxnet", "randomForest",
  "mgcv", "nnet", "pROC", "sysfonts", "showtext", "ggplot2", "viridis"
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

dir.create("output/15_future_env_19plus4", showWarnings = FALSE, recursive = TRUE)
dir.create("output/15_future_env_19plus4/models", showWarnings = FALSE, recursive = TRUE)
dir.create("output/15_future_env_19plus4/rasters", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/15_future_env_19plus4", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("19+4 变量下的未来河网适生度预测\n")
cat("======================================\n\n")

# -----------------------------
# 配置
# -----------------------------

ssp_scenarios <- c("SSP126", "SSP245", "SSP370", "SSP585")
current_bioc_path <- "earthenvstreams_china/hydroclim_weighted_average+sum.tif"
future_bioc_template <- "E:/WorldClim/Future/%s/wc2.1_30s_bioc_BCC-CSM2-MR_%s_2041-2060.tif"

china <- sf::st_read("data-main/vector/china.shp", quiet = TRUE)
bbox <- sf::st_bbox(china)
bbox_ext <- bbox + c(-2, -2, 2, 2)

# -----------------------------
# 步骤1: 构建当前 19+4 训练数据
# -----------------------------

cat("步骤 1/4: 构建当前 19+4 训练数据...\n")

base_df <- read.csv("output/04_collinearity/collinearity_removed.csv")
if (!all(c("id", "lon", "lat", "presence") %in% names(base_df))) {
  stop("缺少 id/lon/lat/presence 于 output/04_collinearity/collinearity_removed.csv")
}

if (!file.exists(current_bioc_path)) {
  stop("未找到当前 WorldClim bioc tif: ", current_bioc_path)
}

bioc_current <- raster::brick(current_bioc_path)
if (raster::nlayers(bioc_current) < 19) {
  stop("当前 bioc 图层少于 19 个，请检查 current_bioc_path")
}

names(bioc_current) <- sprintf("bio%02d", 1:raster::nlayers(bioc_current))

static_list <- list()
if (file.exists("earthenvstreams_china/elevation.tif")) {
  dem_br <- raster::brick("earthenvstreams_china/elevation.tif")[[4]]
  names(dem_br) <- "dem_avg"
  static_list[[length(static_list) + 1]] <- dem_br
}
if (file.exists("earthenvstreams_china/slope.tif")) {
  slope_br <- raster::brick("earthenvstreams_china/slope.tif")[[4]] / 100
  names(slope_br) <- "slope_avg"
  static_list[[length(static_list) + 1]] <- slope_br
}
if (file.exists("earthenvstreams_china/geology_weighted_sum.tif")) {
  geo_total <- raster::calc(raster::brick("earthenvstreams_china/geology_weighted_sum.tif"), fun = function(x) {
    if (all(is.na(x))) NA else sum(x, na.rm = TRUE)
  })
  names(geo_total) <- "geology_total"
  static_list[[length(static_list) + 1]] <- geo_total
}
if (file.exists("earthenvstreams_china/soil_average.tif")) {
  soil_soc <- raster::brick("earthenvstreams_china/soil_average.tif")[[1]]
  names(soil_soc) <- "soil_avg_01"
  static_list[[length(static_list) + 1]] <- soil_soc
}

env_stack_current <- raster::stack(c(bioc_current[[1:19]], static_list))

coords <- base_df[, c("lon", "lat")]
spvals <- raster::extract(env_stack_current, coords)
train_df <- cbind(base_df[, c("id", "species", "lon", "lat", "presence")], spvals)
train_df <- train_df[stats::complete.cases(train_df), ]

sel_vars <- colnames(train_df)[!(colnames(train_df) %in% c("id", "species", "lon", "lat", "presence"))]
write.csv(data.frame(variable = sel_vars), "output/15_future_env_19plus4/selected_variables_19plus4.csv", row.names = FALSE)
write.csv(train_df, "output/15_future_env_19plus4/training_data_19plus4.csv", row.names = FALSE)

cat("  ✓ 变量数: ", length(sel_vars), "\n", sep = "")
cat("  ✓ 样本数: ", nrow(train_df), "\n\n", sep = "")

# -----------------------------
# 步骤2: 19+4 上重训 4 个 SDM
# -----------------------------

cat("步骤 2/4: 19+4 上重训 4 个 SDM...\n")

set.seed(20251121)
pres_idx <- which(train_df$presence == 1)
back_idx <- which(train_df$presence == 0)
train_idx <- c(
  sample(pres_idx, round(0.8 * length(pres_idx))),
  sample(back_idx, round(0.8 * length(back_idx)))
)
test_idx <- setdiff(seq_len(nrow(train_df)), train_idx)

X_train <- train_df[train_idx, sel_vars, drop = FALSE]
X_test <- train_df[test_idx, sel_vars, drop = FALSE]
y_train <- train_df$presence[train_idx]
y_test <- train_df$presence[test_idx]

# Maxnet
cat("  - Maxnet...\n")
mx_model <- maxnet::maxnet(p = y_train, data = X_train, f = maxnet::maxnet.formula(y_train, X_train))
saveRDS(mx_model, "output/15_future_env_19plus4/models/maxnet_19plus4.rds")
pred_test_mx <- as.numeric(predict(mx_model, X_test, type = "logistic"))
auc_mx <- as.numeric(pROC::auc(pROC::roc(y_test, pred_test_mx, quiet = TRUE)))

# RF
cat("  - RF...\n")
rf_model <- randomForest::randomForest(x = X_train, y = factor(y_train), ntree = 500)
saveRDS(rf_model, "output/15_future_env_19plus4/models/rf_19plus4.rds")
pred_test_rf <- as.numeric(predict(rf_model, newdata = X_test, type = "prob")[, "1"])
auc_rf <- as.numeric(pROC::auc(pROC::roc(y_test, pred_test_rf, quiet = TRUE)))

# GAM
cat("  - GAM...\n")
form_terms <- paste0("s(", sel_vars, ")", collapse = " + ")
form <- as.formula(paste0("presence ~ ", form_terms, " + s(lon,lat)"))
gam_model <- mgcv::gam(form, data = cbind(train_df[train_idx, ], X_train), family = binomial(link = "logit"))
saveRDS(gam_model, "output/15_future_env_19plus4/models/gam_19plus4.rds")
pred_test_gam <- as.numeric(predict(gam_model, newdata = cbind(train_df[test_idx, ], X_test), type = "response"))
auc_gam <- as.numeric(pROC::auc(pROC::roc(y_test, pred_test_gam, quiet = TRUE)))

# NN
cat("  - NN...\n")
mu <- sapply(X_train, mean, na.rm = TRUE)
sdv <- sapply(X_train, sd, na.rm = TRUE)
sdv[sdv == 0 | is.na(sdv)] <- 1
X_train_s <- as.data.frame(sweep(sweep(as.matrix(X_train), 2, mu, "-"), 2, sdv, "/"))
X_test_s <- as.data.frame(sweep(sweep(as.matrix(X_test), 2, mu, "-"), 2, sdv, "/"))
size_hidden <- max(3, floor(length(sel_vars) / 2))
nn_model <- nnet::nnet(
  x = X_train_s, y = y_train, size = size_hidden, linout = FALSE,
  rang = 0.1, decay = 5e-4, maxit = 500, trace = FALSE
)
saveRDS(
  list(model = nn_model, mean = mu, sd = sdv, vars = sel_vars),
  "output/15_future_env_19plus4/models/nn_19plus4.rds"
)
pred_test_nn <- as.numeric(nnet:::predict.nnet(nn_model, as.matrix(X_test_s), type = "raw"))
auc_nn <- as.numeric(pROC::auc(pROC::roc(y_test, pred_test_nn, quiet = TRUE)))


eval_all <- dplyr::bind_rows(
  data.frame(model = "Maxnet", AUC = auc_mx),
  data.frame(model = "RF", AUC = auc_rf),
  data.frame(model = "GAM", AUC = auc_gam),
  data.frame(model = "NN", AUC = auc_nn)
)
write.csv(eval_all, "output/15_future_env_19plus4/evaluation_summary_19plus4.csv", row.names = FALSE)

cat("  ✓ 重训完成\n\n")

# -----------------------------
# 步骤3: 未来 bioc01-19 预处理
# -----------------------------

cat("步骤 3/4: 处理未来 bioc01-19...\n")

future_bioc_list <- list()
stats_list <- list()

for (ssp in ssp_scenarios) {
  cat("  处理情景: ", ssp, "\n", sep = "")
  bioc_file <- sprintf(future_bioc_template, ssp, tolower(ssp))
  if (!file.exists(bioc_file)) {
    cat("    ✗ 文件不存在: ", bioc_file, "\n", sep = "")
    next
  }
  bioc_raster <- raster::brick(bioc_file)
  bioc_china <- raster::crop(bioc_raster, raster::extent(bbox_ext))
  bioc_china <- raster::mask(bioc_china, china)
  names(bioc_china) <- sprintf("bio%02d", 1:raster::nlayers(bioc_china))
  out_tif <- file.path("output/15_future_env_19plus4", paste0("future_bioc_china_", ssp, ".tif"))
  raster::writeRaster(bioc_china, out_tif, overwrite = TRUE)
  future_bioc_list[[ssp]] <- bioc_china

  for (i in 1:raster::nlayers(bioc_china)) {
    vals <- raster::getValues(bioc_china[[i]])
    vals <- vals[!is.na(vals)]
    stats_list[[length(stats_list) + 1]] <- data.frame(
      scenario = ssp,
      bioc = sprintf("bio%02d", i),
      mean = mean(vals),
      sd = sd(vals),
      min = min(vals),
      max = max(vals)
    )
  }

  rm(bioc_raster)
  gc(verbose = FALSE)
}

if (length(stats_list) > 0) {
  stats_df <- dplyr::bind_rows(stats_list)
  write.csv(stats_df, "output/15_future_env_19plus4/future_bioc_statistics_19plus4.csv", row.names = FALSE)
}

cat("  ✓ 未来 bioc 处理完成\n\n")

# -----------------------------
# 步骤4: 19+4 变量下的未来河网适生度预测
# -----------------------------

cat("步骤 4/4: 19+4 变量下的未来河网适生度预测...\n")

models_path <- "output/15_future_env_19plus4/models"
mdl_files <- c(
  Maxnet = file.path(models_path, "maxnet_19plus4.rds"),
  RF     = file.path(models_path, "rf_19plus4.rds"),
  GAM    = file.path(models_path, "gam_19plus4.rds"),
  NN     = file.path(models_path, "nn_19plus4.rds")
)
miss_m <- names(mdl_files)[!file.exists(unlist(mdl_files))]
if (length(miss_m) > 0) stop(paste0("缺少模型: ", paste(miss_m, collapse = ", ")))

fa <- raster::brick("earthenvstreams_china/flow_acc.tif")[[2]]
fa_vals <- raster::getValues(fa)
fa_vals[fa_vals <= 0] <- NA
river_mask <- raster::setValues(fa, fa_vals)
rm(fa_vals)

# hydro_ord2_6 <- viz_load_hydrorivers_ord2_6(crs_target = raster::crs(fa)) # Removed: using raster thickening

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

build_future_env_19plus4 <- function(bioc_stack) {
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
  env_stk <- build_future_env_19plus4(bioc_china)
  out_dir_ras <- file.path("output/15_future_env_19plus4/rasters", ssp)
  out_dir_fig <- file.path("figures/15_future_env_19plus4", ssp)
  dir.create(out_dir_ras, showWarnings = FALSE, recursive = TRUE)
  dir.create(out_dir_fig, showWarnings = FALSE, recursive = TRUE)

  summary_rows <- list()
  for (mn in names(mdl_files)) {
    cat("  -> 情景 ", ssp, " | 模型 ", mn, " ...\n", sep = "")
    mdl <- readRDS(mdl_files[[mn]])
    pred_fun <- make_predict_fun(mn, mdl)
    tif_path <- file.path(out_dir_ras, paste0("pred_", tolower(mn), "_19plus4.tif"))
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
    tif_mask_path <- file.path(out_dir_ras, paste0("pred_", tolower(mn), "_19plus4_river.tif"))
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

    out_base <- file.path(out_dir_fig, paste0("prediction_", tolower(mn), "_19plus4"))
    # pred_r_river_thick <- viz_thicken_river_raster(pred_r_river, w_size = 3) # Removed
    viz_save_raster_map(
      r = pred_r_river, out_base = out_base,
      title = paste0(ssp, " - ", mn, " (19+4)"),
      palette = "magma", q_limits = c(0.01, 0.99),
      china_path = "earthenvstreams_china/china_boundary.shp",
      width_in = 8, height_in = 6
      # hydrorivers_sf = hydro_ord2_6 # Removed
    )
  }

  if (length(summary_rows) > 0) {
    summary_df <- dplyr::bind_rows(summary_rows)
    write.csv(summary_df, file.path(out_dir_ras, "prediction_summary_19plus4.csv"), row.names = FALSE)
    summary_df$scenario <- factor(summary_df$scenario, levels = c("SSP126", "SSP245", "SSP370", "SSP585"))
    all_summaries[[ssp]] <- summary_df
  }
}

if (length(all_summaries) > 0) {
  trends <- dplyr::bind_rows(all_summaries)
  out_trend_csv <- "output/15_future_env_19plus4/prediction_trends_all_models_19plus4.csv"
  out_trend_png <- "figures/15_future_env_19plus4/habitat_trends_all_models_19plus4.png"
  write.csv(trends, out_trend_csv, row.names = FALSE)
  p_trend <- ggplot(trends, aes(x = scenario, y = mean, group = model, color = model)) +
    geom_line(linewidth = 0.7) +
    geom_point(size = 1.4) +
    labs(
      title = "Habitat Suitability Trend across SSP Scenarios (19+4)",
      x = "Scenario", y = "Mean Predicted Suitability", color = "Model"
    ) +
    viz_theme_nature(base_size = 8, title_size = 9)
  ggsave(out_trend_png, p_trend, width = 4.8, height = 3.2, dpi = 2400, bg = "transparent")
}

cat("\n======================================\n")
cat("19+4 未来预测完成\n")
cat("======================================\n\n")
