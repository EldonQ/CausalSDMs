#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 11c_shap_contrib_maps.R
# 功能说明: 基于已训练模型与真实环境栅格，计算局部SHAP贡献并在河网内输出空间图；
#          针对 Maxnet / RF / GAM / NN 采用 fastshap 模型无关接口逐块计算SHAP。
# 输入文件: output/05_model_maxnet/model.rds, output/06_model_rf/model.rds,
#          output/07_model_gam/model.rds, output/05b_model_nn/model.rds
#          earthenvstreams_china/*.tif, selected_variables.csv, extracted_variables.csv
# 输出文件: output/11_prediction_maps/rasters/shap_{model}_{var}.tif
#          figures/11_prediction_maps/shap_{model}_{var}.png
#          output/11_prediction_maps/shap_maps_summary.csv
# 作者: Nature级别科研项目
# 日期: 2025-10-24
# ==============================================================================

# 初始化
rm(list = ls())
gc()
setwd("E:/SDM01")

packages <- c("tidyverse", "raster", "fastshap", "nnet", "randomForest", "maxnet", "mgcv", "sf", "rnaturalearth", "viridis", "sysfonts", "showtext", "terra", "svglite", "ggplot2")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

dir.create("output/11_prediction_maps/rasters", showWarnings = FALSE, recursive = TRUE)

# 字体与制图
try({
  sysfonts::font_add(
    family = "Arial",
    regular = "C:/Windows/Fonts/arial.ttf",
    bold = "C:/Windows/Fonts/arialbd.ttf",
    italic = "C:/Windows/Fonts/ariali.ttf",
    bolditalic = "C:/Windows/Fonts/arialbi.ttf"
  )
  showtext::showtext_opts(dpi = 1200)
  showtext::showtext_auto(enable = TRUE)
}, silent = TRUE)

# 统一绘图工具（Nature风格/PNG+SVG/Arial）
source("scripts/visualization/viz_utils.R")

cat("\n======================================\n")
cat("局部SHAP贡献空间图 (真实数据，河网掩膜)\n")
cat("======================================\n\n")

# 变量与环境栅格
sel_vars <- read.csv("output/04_collinearity/selected_variables.csv", stringsAsFactors = FALSE)$variable
var_map <- read.csv("output/02_env_extraction/extracted_variables.csv", stringsAsFactors = FALSE)
var_map <- var_map[var_map$variable %in% sel_vars, c("variable", "file", "band")]

build_env_stack <- function(var_map_df, base_dir = "earthenvstreams_china") {
  groups <- split(var_map_df, var_map_df$file)
  stk_list <- list()
  for(fn in names(groups)) {
    g <- groups[[fn]]
    r <- raster::brick(file.path(base_dir, fn))
    r_sel <- r[[g$band]]
    names(r_sel) <- g$variable
    stk_list[[length(stk_list) + 1]] <- r_sel
    rm(r, r_sel)
    gc(verbose = FALSE)
  }
  stk <- raster::stack(stk_list)
  stk <- stk[[sel_vars]]
  return(stk)
}

env_stack <- build_env_stack(var_map)

# 若后续模型（如GAM）包含 s(lon,lat)，需要在环境栅格中加入经纬度图层
lon_r <- raster::init(env_stack[[1]], fun = 'x') ; names(lon_r) <- 'lon'
lat_r <- raster::init(env_stack[[1]], fun = 'y') ; names(lat_r) <- 'lat'
env_stack <- raster::addLayer(env_stack, lon_r, lat_r)

# 河网掩膜
fa <- raster::brick("earthenvstreams_china/flow_acc.tif")[[2]]
fa_vals <- raster::getValues(fa) ; fa_vals[fa_vals <= 0] <- NA
river_mask <- raster::setValues(fa, fa_vals)
rm(fa_vals)

# 中国边界
china <- rnaturalearth::ne_countries(country = "China", scale = "medium", returnclass = "sf")

# 模型与预测封装
model_paths <- c(
  Maxnet = "output/05_model_maxnet/model.rds",
  RF = "output/06_model_rf/model.rds",
  GAM = "output/07_model_gam/model.rds",
  NN = "output/05b_model_nn/model.rds"
)

make_pred_fun <- function(model_name, model_obj) {
  if(model_name == "Maxnet") {
    return(function(object, newdata) { as.numeric(predict(object, newdata, type = "logistic")) })
  }
  if(model_name == "RF") {
    return(function(object, newdata) { as.numeric(predict(object, newdata = newdata, type = "prob")[, "1"]) })
  }
  if(model_name == "GAM") {
    return(function(object, newdata) { as.numeric(predict(object, newdata = newdata, type = "response")) })
  }
  if(model_name == "NN") {
    return(function(object, newdata) {
      mu <- object$mean; sdv <- object$sd; mod <- object$model; vars <- object$vars
      sdv[sdv == 0 | is.na(sdv)] <- 1
      x <- as.matrix(newdata[, vars, drop = FALSE])
      x <- sweep(x, 2, mu[vars], "-")
      x <- sweep(x, 2, sdv[vars], "/")
      as.numeric(nnet:::predict.nnet(mod, x, type = "raw"))
    })
  }
}

available_models <- names(model_paths)[file.exists(model_paths)]
if(length(available_models) == 0) stop("未发现已训练模型，无法计算SHAP贡献图")

# 从变量重要性中选取Top变量，控制图件数量
imp_path <- "output/09_variable_importance/importance_summary.csv"
if(file.exists(imp_path)) {
  imp_df <- read.csv(imp_path)
  top_vars <- imp_df %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(mean_imp = mean(importance_normalized, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(mean_imp)) %>%
    dplyr::pull(variable)
  top_vars <- intersect(top_vars, sel_vars)
  if(length(top_vars) > 12) top_vars <- top_vars[1:12]
} else {
  top_vars <- sel_vars[seq_len(min(12, length(sel_vars)))]
}

# 分块设置
bs <- raster::blockSize(env_stack)
summary_rows <- list()

for(mn in available_models) {
  cat("模型 ", mn, " 的 SHAP 空间映射...\n", sep = "")
  mdl <- readRDS(model_paths[[mn]])
  pred_fun <- make_pred_fun(mn, mdl)

  # 注意：fastshap::explain 需要一批样本；这里对每个块内像元使用 X_block 直接求SHAP
  # 计算成本高：可适当减少 nsim 或 top_vars 数量
  for(v in top_vars) {
    cat("  - 变量 ", v, " ...\n", sep = "")
    out_path <- file.path("output/11_prediction_maps/rasters", paste0("shap_", tolower(mn), "_", gsub("[^A-Za-z0-9_]+", "_", v), ".tif"))
    if(file.exists(out_path)) { try({ file.remove(out_path) }, silent = TRUE) }
    out_r <- raster::raster(env_stack, layer = 1)
    out_r <- raster::setValues(out_r, NA_real_)

    wr <- raster::writeStart(out_r, filename = out_path, overwrite = TRUE)

    for(i in seq_len(bs$n)) {
      X_block <- raster::getValues(env_stack, row = bs$row[i], nrows = bs$nrows[i])
      X_df <- as.data.frame(X_block)
      if(nrow(X_df) == 0) {
        wr <- raster::writeValues(wr, rep(NA_real_, 0), bs$row[i])
        next
      }
      # 删除包含NA的行，SHAP在NA时将无法计算；对应像元结果记为NA
      na_rows <- !stats::complete.cases(X_df)
      X_ok <- X_df[!na_rows, , drop = FALSE]

      if(nrow(X_ok) > 0) {
        # fastshap 仅对使用到的变量计算；确保变量集合对齐
        vars <- if(mn == "NN") mdl$vars else sel_vars
        vars <- intersect(vars, colnames(X_ok))
        # 若为GAM，预测需要 lon/lat；将其并入 X_used，但 feature_names 仍仅对 vars 求SHAP
        if(mn == "GAM") {
          extra_xy <- intersect(c("lon","lat"), colnames(X_ok))
          union_vars <- unique(c(vars, extra_xy))
          X_used <- X_ok[, union_vars, drop = FALSE]
        } else {
          X_used <- X_ok[, vars, drop = FALSE]
        }

        shap_mat <- try({
          fastshap::explain(
            object = mdl,
            X = X_used,
            pred_wrapper = pred_fun,
            feature_names = vars,
            nsim = 32,
            adjust = TRUE
          )
        }, silent = TRUE)

        if(inherits(shap_mat, "try-error")) {
          vec <- rep(NA_real_, nrow(X_df))
        } else {
          shap_df <- as.data.frame(shap_mat)
          vec <- rep(NA_real_, nrow(X_df))
          if(v %in% colnames(shap_df)) {
            vec[!na_rows] <- shap_df[[v]]
          }
        }
      } else {
        vec <- rep(NA_real_, nrow(X_df))
      }

      wr <- raster::writeValues(wr, vec, bs$row[i])
    }
    out_r <- raster::writeStop(wr)

    # 掩膜到河网
    out_riv <- raster::mask(out_r, river_mask)
    raster::writeRaster(out_riv, out_path, overwrite = TRUE)

    # 统计摘要
    vals <- raster::getValues(out_riv)
    vals <- vals[is.finite(vals)]
    if(length(vals) > 0) {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        model = mn, variable = v,
        n_pixels_river = length(vals),
        mean = mean(vals), sd = sd(vals),
        min = min(vals), max = max(vals),
        p10 = as.numeric(stats::quantile(vals, 0.1)),
        p50 = as.numeric(stats::quantile(vals, 0.5)),
        p90 = as.numeric(stats::quantile(vals, 0.9))
      )
    } else {
      summary_rows[[length(summary_rows) + 1]] <- data.frame(
        model = mn, variable = v,
        n_pixels_river = 0, mean = NA, sd = NA, min = NA, max = NA, p10 = NA, p50 = NA, p90 = NA
      )
    }

    # 统一出图风格（PNG+SVG）
    out_base <- file.path("figures/11_prediction_maps", paste0("shap_", tolower(mn), "_", gsub("[^A-Za-z0-9_]+", "_", v)))
    viz_save_raster_map(r = out_riv, out_base = out_base,
                        title = paste("Local SHAP:", mn, "-", v),
                        palette = "magma", q_limits = c(0.01, 0.99),
                        china_path = "earthenvstreams_china/china_boundary.shp",
                        width_in = 8, height_in = 6)

    
  }
}

# 汇总表
if(length(summary_rows) > 0) {
  shap_sum <- dplyr::bind_rows(summary_rows)
  write.csv(shap_sum, "output/11_prediction_maps/shap_maps_summary.csv", row.names = FALSE)
}

cat("\n完成：局部SHAP贡献空间图与统计输出\n\n")
