#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 15a_retrain_future_vars.R
# 功能说明: 基于仅“未来可获得”的气候与静态变量，使用当前真实数据重训4个模型
# 变量集合: temp_ann_mean(由当前 monthly_tmax/tmin 平均)、prec_ann_sum(由当前 monthly_prec_sum 汇总)
#          + elevation, slope, geology_weighted_sum, soil_average（静态）
# 输出: output/15_future_env/training_data_future_vars.csv
#      output/15_future_env/selected_variables_future.csv
#      output/15_future_env/models/{maxnet,nn,rf,gam}.rds
#      output/15_future_env/evaluation_*.csv
# 图件: figures/15_future_env/model_eval_*.png（Arial, 1200dpi）
# 作者: Nature级别科研项目
# 日期: 2025-10-24
# ============================================================================

# 初始化
rm(list = ls())
gc()
setwd("E:/SDM01")

packages <- c("tidyverse", "raster", "sf", "maxnet", "randomForest", "nnet", "mgcv", "pROC", "sysfonts", "showtext")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

dir.create("output/15_future_env", showWarnings = FALSE, recursive = TRUE)
dir.create("output/15_future_env/models", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/15_future_env", showWarnings = FALSE, recursive = TRUE)

# 字体
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

cat("\n======================================\n")
cat("未来可用变量 重筛+重训\n")
cat("======================================\n\n")

# 统一可视化工具（主题字号保持一致，便于后续可能的图件扩展）
source("scripts/visualization/viz_utils.R")

# 读取样本（仅使用真实出现/背景样本）
base_df <- read.csv("output/04_collinearity/collinearity_removed.csv")
if(!all(c("id","lon","lat","presence") %in% names(base_df))) {
  stop("缺少关键列 id/lon/lat/presence 于 output/04_collinearity/collinearity_removed.csv")
}

# 构建当前“未来可用变量”栅格栈（基于 earthenvstreams_china）
# 1) 年平均温度 temp_ann_mean = (mean(monthly_tmax_average) + mean(monthly_tmin_average)) / 2
# 2) 年降水量 prec_ann_sum = sum(monthly_prec_sum)
# 3) 静态: elevation.tif, slope.tif, geology_weighted_sum.tif, soil_average.tif

get_band_mean <- function(br) { raster::calc(br, fun = function(x){ if(all(is.na(x))) NA else mean(x, na.rm = TRUE) }) }
get_band_sum  <- function(br) { raster::calc(br, fun = function(x){ if(all(is.na(x))) NA else sum(x, na.rm = TRUE) }) }

# 读取月度栅格（12波段）
tx_br <- raster::brick("earthenvstreams_china/monthly_tmax_average.tif") / 10  # ℃×10 → ℃
tn_br <- raster::brick("earthenvstreams_china/monthly_tmin_average.tif") / 10
pr_br <- raster::brick("earthenvstreams_china/monthly_prec_sum.tif")

# 计算年度汇总
ann_tmax <- get_band_mean(tx_br)
ann_tmin <- get_band_mean(tn_br)
ann_temp_mean <- (ann_tmax + ann_tmin) / 2; names(ann_temp_mean) <- "temp_ann_mean"
ann_prec_sum <- get_band_sum(pr_br); names(ann_prec_sum) <- "prec_ann_sum"

# 静态层（选取与47变量体系一致的代表性band，并进行单位换算）
static_list <- list()

# 平均高程 dem_avg (band 4)
if(file.exists("earthenvstreams_china/elevation.tif")) {
  dem_br <- raster::brick("earthenvstreams_china/elevation.tif")[[4]]
  names(dem_br) <- "dem_avg"
  static_list[[length(static_list)+1]] <- dem_br
}

# 平均坡度 slope_avg (band 4，度×100 → 度)
if(file.exists("earthenvstreams_china/slope.tif")) {
  slope_br <- raster::brick("earthenvstreams_china/slope.tif")[[4]] / 100
  names(slope_br) <- "slope_avg"
  static_list[[length(static_list)+1]] <- slope_br
}

# 地质：所有 geology weighted sum band 求和（静态特征），保持与未来可用一致
if(file.exists("earthenvstreams_china/geology_weighted_sum.tif")) {
  geo_total <- raster::calc(raster::brick("earthenvstreams_china/geology_weighted_sum.tif"), fun = function(x) {
    if(all(is.na(x))) NA else sum(x, na.rm = TRUE)
  })
  names(geo_total) <- "geology_total"
  static_list[[length(static_list)+1]] <- geo_total
}

# 土壤：soil_average.tif 第1波段（soil_avg_01：SOC），单位保持 g/kg
if(file.exists("earthenvstreams_china/soil_average.tif")) {
  soil_soc <- raster::brick("earthenvstreams_china/soil_average.tif")[[1]]
  names(soil_soc) <- "soil_avg_01"
  static_list[[length(static_list)+1]] <- soil_soc
}

# 合栈
env_stack <- raster::stack(c(list(ann_temp_mean, ann_prec_sum), static_list))
# 补充 lon/lat（供GAM使用）
lon_r <- raster::init(env_stack[[1]], fun = 'x'); names(lon_r) <- 'lon'
lat_r <- raster::init(env_stack[[1]], fun = 'y'); names(lat_r) <- 'lat'

# 提取样本环境
coords <- base_df[, c("lon","lat")]
spvals <- raster::extract(env_stack, coords)
train_df <- cbind(base_df[, c("id","species","lon","lat","presence")], spvals)

# 去除含NA的样本
train_df <- train_df[stats::complete.cases(train_df), ]

# 选用变量清单（未来场景可获得的6项）
sel_vars <- c("temp_ann_mean","prec_ann_sum","dem_avg","slope_avg","geology_total","soil_avg_01")
sel_vars <- intersect(sel_vars, colnames(train_df))
write.csv(data.frame(variable = sel_vars), "output/15_future_env/selected_variables_future.csv", row.names = FALSE)
write.csv(train_df, "output/15_future_env/training_data_future_vars.csv", row.names = FALSE)

cat("  ✓ 变量数: ", length(sel_vars), "\n", sep = "")
cat("  ✓ 样本数: ", nrow(train_df), "\n\n", sep = "")

# 分层划分
set.seed(20251024)
pres_idx <- which(train_df$presence == 1)
back_idx <- which(train_df$presence == 0)
train_idx <- c(sample(pres_idx, round(0.8 * length(pres_idx))), sample(back_idx, round(0.8 * length(back_idx))))
test_idx <- setdiff(seq_len(nrow(train_df)), train_idx)

X_train <- train_df[train_idx, sel_vars, drop = FALSE]
X_test  <- train_df[test_idx,  sel_vars, drop = FALSE]
y_train <- train_df$presence[train_idx]
y_test  <- train_df$presence[test_idx]

# 训练 Maxnet
cat("训练 Maxnet...\n")
mx_model <- maxnet::maxnet(p = y_train, data = X_train, f = maxnet.formula(y_train, X_train))
saveRDS(mx_model, "output/15_future_env/models/maxnet.rds")
pred_test_mx <- as.numeric(predict(mx_model, X_test, type = "logistic"))
auc_mx <- as.numeric(pROC::auc(pROC::roc(y_test, pred_test_mx, quiet = TRUE)))
write.csv(data.frame(model = "Maxnet", AUC = auc_mx), "output/15_future_env/evaluation_maxnet.csv", row.names = FALSE)

# 训练 RF
cat("训练 RF...\n")
rf_model <- randomForest::randomForest(x = X_train, y = factor(y_train), ntree = 500)
saveRDS(rf_model, "output/15_future_env/models/rf.rds")
pred_test_rf <- as.numeric(predict(rf_model, newdata = X_test, type = "prob")[, "1"])
auc_rf <- as.numeric(pROC::auc(pROC::roc(y_test, pred_test_rf, quiet = TRUE)))
write.csv(data.frame(model = "RF", AUC = auc_rf), "output/15_future_env/evaluation_rf.csv", row.names = FALSE)

# 训练 GAM（含 s(lon,lat)）
cat("训练 GAM...\n")
form_terms <- paste0("s(", sel_vars, ")", collapse = " + ")
form <- as.formula(paste0("presence ~ ", form_terms, " + s(lon,lat)") )
gam_model <- mgcv::gam(form, data = cbind(train_df[train_idx, ], X_train), family = binomial(link = "logit"))
saveRDS(gam_model, "output/15_future_env/models/gam.rds")
pred_test_gam <- as.numeric(predict(gam_model, newdata = cbind(train_df[test_idx, ], X_test), type = "response"))
auc_gam <- as.numeric(pROC::auc(pROC::roc(y_test, pred_test_gam, quiet = TRUE)))
write.csv(data.frame(model = "GAM", AUC = auc_gam), "output/15_future_env/evaluation_gam.csv", row.names = FALSE)

# 训练 NN（标准化）
cat("训练 NN...\n")
mu <- sapply(X_train, mean, na.rm = TRUE)
sdv <- sapply(X_train, sd, na.rm = TRUE); sdv[sdv == 0 | is.na(sdv)] <- 1
X_train_s <- as.data.frame(sweep(sweep(as.matrix(X_train), 2, mu, "-"), 2, sdv, "/"))
X_test_s  <- as.data.frame(sweep(sweep(as.matrix(X_test),  2, mu, "-"), 2, sdv, "/"))
size_hidden <- max(3, floor(length(sel_vars) / 2))
nn_model <- nnet::nnet(x = X_train_s, y = y_train, size = size_hidden, linout = FALSE, rang = 0.1, decay = 5e-4, maxit = 500, trace = FALSE)
saveRDS(list(model = nn_model, mean = mu, sd = sdv, vars = sel_vars), "output/15_future_env/models/nn.rds")
pred_test_nn <- as.numeric(nnet:::predict.nnet(nn_model, as.matrix(X_test_s), type = "raw"))
auc_nn <- as.numeric(pROC::auc(pROC::roc(y_test, pred_test_nn, quiet = TRUE)))
write.csv(data.frame(model = "NN", AUC = auc_nn), "output/15_future_env/evaluation_nn.csv", row.names = FALSE)

# 汇总评估
eval_all <- dplyr::bind_rows(
  read.csv("output/15_future_env/evaluation_maxnet.csv"),
  read.csv("output/15_future_env/evaluation_rf.csv"),
  read.csv("output/15_future_env/evaluation_gam.csv"),
  read.csv("output/15_future_env/evaluation_nn.csv")
)
write.csv(eval_all, "output/15_future_env/evaluation_summary.csv", row.names = FALSE)

cat("\n======================================\n")
cat("重训完成\n")
cat("======================================\n\n")

cat("✓ 变量清单: output/15_future_env/selected_variables_future.csv\n")
cat("✓ 训练数据: output/15_future_env/training_data_future_vars.csv\n")
cat("✓ 模型: output/15_future_env/models/*.rds\n")
cat("✓ 评估: output/15_future_env/evaluation_summary.csv\n\n")
