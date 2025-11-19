#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 03_response_curves_clean.R
# 功能说明: 重新绘制10_response_curves.R的所有结果（简洁版：保留坐标轴、无数字、透明画布、无网格）
#           包含：GAM响应曲线、ALE曲线（多模型）、所有统计输出
# 参考脚本: scripts/10_response_curves.R（完整功能复刻，仅改变绘图样式）
# 输入文件: output/07_model_gam/model.rds
#          output/09_variable_importance/importance_summary.csv
#          output/04_collinearity/collinearity_removed.csv
#          output/05_model_maxnet/model.rds
#          output/06_model_rf/model.rds
#          output/05b_model_nn/model.rds
# 输出文件: Drawing/output/response_curves/individual/*.png
#          Drawing/output/response_curves/response_curves_top10.png
#          Drawing/output/response_curves/ale/*.png
#          Drawing/output/response_curves/ale/*.csv
# 作者: Nature级别科研项目
# 日期: 2025-11-09
# ==============================================================================

# -----------------------------
# 0. 初始化环境
# -----------------------------
rm(list = ls())
gc()
setwd("E:/SDM01")

# -----------------------------
# 1. 加载必要的包（缺则自动安装）
# -----------------------------
packages <- c(
  "tidyverse",     # 数据整理
  "mgcv",          # GAM模型
  "ggplot2",       # 绘图
  "cowplot",       # 图形组合
  "iml",           # ALE曲线
  "nnet",          # NN模型
  "randomForest",  # RF模型
  "maxnet"         # Maxnet模型
)
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 统一可视化工具（Nature风格、Arial、高dpi）
source("scripts/visualization/viz_utils.R")

# -----------------------------
# 2. 目录准备
# -----------------------------
dir.create("Drawing/output/response_curves", showWarnings = FALSE, recursive = TRUE)
dir.create("Drawing/output/response_curves/individual", showWarnings = FALSE, recursive = TRUE)
dir.create("Drawing/output/response_curves/ale", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("响应曲线重绘（简洁版：保留轴线、无数字、透明背景）\n")
cat("======================================\n\n")

# -----------------------------
# 3. 读取模型和变量重要性
# -----------------------------
cat("步骤 1/4: 读取模型和数据...\n")

gam_model <- readRDS("output/07_model_gam/model.rds")
var_importance <- read.csv("output/09_variable_importance/importance_summary.csv") %>%
  filter(model == "GAM", variable != "lon,lat")

# 选择Top 10变量
top_vars <- var_importance %>%
  arrange(desc(importance_normalized)) %>%
  head(10) %>%
  pull(variable)

cat("  ✓ GAM模型已加载\n")
cat("  ✓ Top 10变量: ", length(top_vars), "\n", sep = "")

# 读取训练数据（用于概率尺度响应曲线：固定其他变量为中位数/最常值）
train_df <- read.csv("output/04_collinearity/collinearity_removed.csv")
exclude_cols <- c("id", "species", "source", "presence", "presence.1")
all_predictors <- setdiff(names(train_df), exclude_cols)

# 构建基准观测：数值型取中位数，类别型取众数
base_row <- as.list(train_df[1, all_predictors, drop = TRUE])
for (nm in all_predictors) {
  v <- train_df[[nm]]
  if (is.numeric(v)) {
    base_row[[nm]] <- stats::median(v, na.rm = TRUE)
  } else {
    lv <- names(sort(table(v), decreasing = TRUE))[1]
    base_row[[nm]] <- if (is.na(lv)) NA else lv
  }
}
base_row <- as.data.frame(base_row, stringsAsFactors = FALSE)

# -----------------------------
# 4. 绘制单个响应曲线（简洁版）
# -----------------------------
cat("\n步骤 2/4: 绘制单变量响应曲线（简洁版）...\n")

plot_list <- list()
curves_data <- list()

for (i in seq_along(top_vars)) {
  var <- top_vars[i]
  cat("  - ", var, "\n", sep = "")
  
  vx <- train_df[[var]]
  vx <- vx[is.finite(as.numeric(vx))]
  if (length(vx) == 0) next
  
  rng <- stats::quantile(as.numeric(vx), probs = c(0.01, 0.99), na.rm = TRUE)
  x_seq <- seq(rng[1], rng[2], length.out = 200)
  
  newd <- base_row[rep(1, length(x_seq)), , drop = FALSE]
  newd[[var]] <- x_seq
  
  pred <- as.numeric(predict(gam_model, newdata = newd, type = "response"))
  dfp <- data.frame(x = x_seq, y = pred)
  
  # 保存数据用于统计
  curves_data[[length(curves_data) + 1]] <- data.frame(
    variable = var,
    rank = i,
    x_min = rng[1],
    x_max = rng[2],
    y_min = min(pred, na.rm = TRUE),
    y_max = max(pred, na.rm = TRUE),
    y_mean = mean(pred, na.rm = TRUE)
  )
  
  # 简洁版主题：保留坐标轴线、无网格、透明画布、无坐标数字
  p <- ggplot(dfp, aes(x = x, y = y)) +
    geom_line(linewidth = 0.7, color = "black") +
    labs(title = paste0("Response: ", var), 
         x = var, 
         y = "Presence Probability") +
    coord_cartesian(ylim = c(0, 1)) +
    viz_theme_nature(base_size = 9) +
    theme(
      # 去掉背景网格
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      # 透明画布
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      # 去掉坐标轴数字，但保留坐标轴线
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      # 保留坐标轴线（简单的两条线）
      axis.line = element_line(color = "black", linewidth = 0.5),
      # 保留坐标轴标题
      axis.title.x = element_text(size = 9, family = "Arial"),
      axis.title.y = element_text(size = 9, family = "Arial"),
      plot.title = element_text(size = 10, family = "Arial", face = "bold")
    )
  
  # 保存单个图形（PNG 2400 dpi + SVG，透明背景）
  ggsave(
    filename = paste0("Drawing/output/response_curves/individual/", var, ".png"),
    plot = p, width = 2.4, height = 2.4, units = "in", dpi = 2400, bg = "transparent"
  )
  ggsave(
    filename = paste0("Drawing/output/response_curves/individual/", var, ".svg"),
    plot = p, width = 2.4, height = 2.4, bg = "transparent"
  )
  
  plot_list[[length(plot_list) + 1]] <- p
}

cat("  ✓ 单变量曲线: Drawing/output/response_curves/individual/\n")

# -----------------------------
# 5. 绘制组合图（Top 10）
# -----------------------------
cat("\n步骤 3/4: 绘制Top 10组合图...\n")

if (length(plot_list) > 0) {
  comb <- cowplot::plot_grid(plotlist = plot_list, ncol = 2, align = "hv")
  
  ggsave(
    "Drawing/output/response_curves/response_curves_top10.png",
    plot = comb, width = 4.8, height = 6, units = "in", dpi = 2400, bg = "transparent"
  )
  ggsave(
    "Drawing/output/response_curves/response_curves_top10.svg",
    plot = comb, width = 4.8, height = 6, bg = "transparent"
  )
  
  cat("  ✓ 组合图: Drawing/output/response_curves/response_curves_top10.png\n")
}

# -----------------------------
# 6. 输出统计表与日志
# -----------------------------
cat("\n步骤 4/4: 输出统计信息...\n")

if (length(curves_data) > 0) {
  curves_summary <- dplyr::bind_rows(curves_data)
  write.csv(curves_summary, "Drawing/output/response_curves/curves_summary.csv", row.names = FALSE)
  cat("  ✓ 曲线统计: Drawing/output/response_curves/curves_summary.csv\n")
}

# 保存变量列表
top_vars_df <- data.frame(
  rank = seq_along(top_vars),
  variable = top_vars,
  importance = var_importance %>%
    filter(variable %in% top_vars) %>%
    arrange(desc(importance_normalized)) %>%
    head(10) %>%
    pull(importance_normalized)
)
write.csv(top_vars_df, "Drawing/output/response_curves/top_variables.csv", row.names = FALSE)

# -----------------------------
# 7. ALE曲线计算（模型无关解释，多模型）
# -----------------------------
cat("\n步骤 5/5: 计算 ALE 曲线（多模型）...\n")

# 读取用于建模的数据集
model_data <- read.csv("output/04_collinearity/collinearity_removed.csv")
exclude_cols <- c("id", "species", "lon", "lat", "source", "presence", "presence.1")
env_vars <- setdiff(names(model_data), exclude_cols)
X_all <- model_data[, env_vars, drop = FALSE]
y_all <- model_data$presence

# 若GAM包含 s(lon,lat)，提供经纬度
has_lon <- "lon" %in% names(model_data)
has_lat <- "lat" %in% names(model_data)
lonlat_df <- if (has_lon && has_lat) model_data[, c("lon", "lat"), drop = FALSE] else NULL

# 读取各模型
model_paths <- c(
  Maxnet = "output/05_model_maxnet/model.rds",
  RF = "output/06_model_rf/model.rds",
  GAM = "output/07_model_gam/model.rds",
  NN = "output/05b_model_nn/model.rds"
)

available_models <- names(model_paths)[file.exists(model_paths)]

if (length(available_models) == 0) {
  cat("  ✗ 未发现已训练模型，跳过 ALE 计算\n")
} else {
  # 选择ALE变量（Top 15）
  imp_path <- "output/09_variable_importance/importance_summary.csv"
  if (file.exists(imp_path)) {
    imp_df <- read.csv(imp_path)
    top_from_imp <- imp_df %>%
      group_by(variable) %>%
      summarise(mean_imp = mean(importance_normalized, na.rm = TRUE), .groups = "drop") %>%
      arrange(desc(mean_imp)) %>%
      pull(variable)
  } else {
    top_from_imp <- env_vars
  }
  ale_vars <- intersect(top_from_imp, env_vars)
  if (length(ale_vars) > 15) ale_vars <- ale_vars[1:15]
  
  # 预测函数工厂
  make_pred_fun <- function(model_name, model_obj) {
    if (model_name == "Maxnet") {
      return(function(object, newdata) { 
        as.numeric(predict(object, newdata, type = "logistic")) 
      })
    }
    if (model_name == "RF") {
      return(function(object, newdata) { 
        as.numeric(predict(object, newdata = newdata, type = "prob")[, "1"]) 
      })
    }
    if (model_name == "GAM") {
      return(function(object, newdata) { 
        as.numeric(predict(object, newdata = newdata, type = "response")) 
      })
    }
    if (model_name == "NN") {
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
  
  # 计算并输出ALE
  ale_all <- list()
  for (mn in available_models) {
    cat("  -> 模型 ", mn, " 的 ALE ...\n", sep = "")
    mdl <- readRDS(model_paths[[mn]])
    pred_fun <- make_pred_fun(mn, mdl)
    
    # 若为GAM且存在经纬度，附加 lon/lat
    data_for_model <- X_all
    if (mn == "GAM" && !is.null(lonlat_df)) {
      for (cc in c("lon", "lat")) {
        if (!(cc %in% colnames(data_for_model)) && (cc %in% colnames(lonlat_df))) {
          data_for_model[[cc]] <- lonlat_df[[cc]]
        }
      }
    }
    
    predictor <- iml::Predictor$new(
      model = mdl,
      data = data_for_model,
      y = y_all,
      predict.function = pred_fun,
      class = NULL
    )
    
    for (v in ale_vars) {
      fe <- iml::FeatureEffect$new(predictor, feature = v, method = "ale", grid.size = 40)
      
      # 保存CSV
      res <- fe$results
      res$model <- mn
      res$variable <- v
      v_sanit <- gsub("[^A-Za-z0-9_]+", "_", v)
      write.csv(res, 
                file = paste0("Drawing/output/response_curves/ale/ale_", tolower(mn), "_", v_sanit, ".csv"), 
                row.names = FALSE)
      
      # 绘制ALE图（简洁版：保留轴线、无数字、透明背景）
      plt <- plot(fe)
      
      # 移除rug图层
      try({
        plt$layers <- Filter(function(ly) { !inherits(ly$geom, "GeomRug") }, plt$layers)
      }, silent = TRUE)
      
      # 应用简洁版主题
      plt <- plt + 
        labs(title = paste0("ALE - ", mn, ": ", v), x = v, y = "ALE of .y") +
        theme_minimal(base_size = 9, base_family = "Arial") +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          axis.line = element_line(color = "black", linewidth = 0.5),
          axis.title.x = element_text(size = 9, family = "Arial"),
          axis.title.y = element_text(size = 9, family = "Arial"),
          plot.title = element_text(size = 10, family = "Arial", face = "bold")
        )
      
      # 保存PNG（2400 dpi，透明背景）
      ggsave(
        filename = paste0("Drawing/output/response_curves/ale/ale_", tolower(mn), "_", v_sanit, ".png"),
        plot = plt, width = 2.4, height = 2.4, units = "in", dpi = 2400, bg = "transparent"
      )
      
      # 保存SVG
      ggsave(
        filename = paste0("Drawing/output/response_curves/ale/ale_", tolower(mn), "_", v_sanit, ".svg"),
        plot = plt, width = 2.4, height = 2.4, bg = "transparent"
      )
      
      ale_all[[length(ale_all) + 1]] <- res
      rm(fe)
      gc(verbose = FALSE)
    }
  }
  
  # 汇总ALE结果
  if (length(ale_all) > 0) {
    ale_df <- dplyr::bind_rows(ale_all)
    write.csv(ale_df, "Drawing/output/response_curves/ale/ale_summary.csv", row.names = FALSE)
    cat("  ✓ ALE 结果已保存至 Drawing/output/response_curves/ale/\n")
  } else {
    cat("  ✗ 未产出 ALE 结果\n")
  }
}

# -----------------------------
# 8. 输出日志
# -----------------------------
sink("Drawing/output/response_curves/processing_log.txt")
cat("响应曲线重绘日志（简洁版）\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
cat("绘图设置:\n")
cat("  - 坐标轴: 保留简单的两条线\n")
cat("  - 坐标轴数字: 已移除\n")
cat("  - 背景网格: 已移除\n")
cat("  - 画布背景: 透明\n")
cat("  - 输出分辨率: 2400 dpi\n\n")
cat("Top 10变量:\n")
print(top_vars_df)
if (length(curves_data) > 0) {
  cat("\n曲线统计:\n")
  print(curves_summary)
}
if (exists("ale_all") && length(ale_all) > 0) {
  cat("\nALE 计算完成:\n")
  cat("  - 模型数: ", length(available_models), "\n", sep = "")
  cat("  - 变量数: ", length(ale_vars), "\n", sep = "")
  cat("  - 总曲线数: ", length(ale_all), "\n", sep = "")
}
sink()

cat("\n======================================\n")
cat("响应曲线重绘完成（完整功能）\n")
cat("======================================\n\n")

cat(sprintf("✓ GAM响应曲线: %d 个变量\n", length(top_vars)))
cat(sprintf("✓ 单图输出: Drawing/output/response_curves/individual/\n"))
cat(sprintf("✓ 组合图: Drawing/output/response_curves/response_curves_top10.png\n"))
if (length(available_models) > 0) {
  cat(sprintf("✓ ALE 曲线: %d 个模型 × %d 个变量\n", length(available_models), length(ale_vars)))
  cat(sprintf("✓ ALE 输出: Drawing/output/response_curves/ale/\n"))
}
cat(sprintf("✓ 分辨率: 2400 dpi（超高质量）\n"))
cat(sprintf("✓ 样式: 保留坐标轴线，无数字，透明背景\n"))

cat("\n✓ 脚本执行完成!\n\n")

