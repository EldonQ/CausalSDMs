#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 09_variable_importance_viz.R
# 功能说明: 可视化并对比所有模型的变量重要性（小提琴图）
# 输入文件: output/05_model_maxnet/variable_importance.csv
#          output/06_model_rf/variable_importance.csv
#          output/07_model_gam/variable_importance.csv
# 输出文件: figures/09_variable_importance/importance_violin.png
#          figures/09_variable_importance/importance_heatmap.png
#          output/09_variable_importance/importance_summary.csv
# 作者: Nature级别科研项目
# 日期: 2025-10-20
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/SDM01")

# 加载必要的包
packages <- c("tidyverse", "ggplot2", "viridis", "reshape2", "pheatmap", "fastshap", "iml", "sysfonts", "showtext", "nnet", "randomForest", "maxnet", "mgcv", "patchwork")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 中文注释：注册并启用 Arial 字体，确保 PNG/SVG 输出嵌入或渲染为 Arial（期刊级别）
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

dir.create("output/09_variable_importance", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/09_variable_importance", showWarnings = FALSE, recursive = TRUE)
dir.create("output/09_variable_importance/shap", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/09_variable_importance/shap", showWarnings = FALSE, recursive = TRUE)

# 统一可视化工具与主题（Nature风格、Arial、1200dpi、viridis）
source("scripts/visualization/viz_utils.R")

cat("\n======================================\n")
cat("变量重要性可视化 (3个模型)\n")
cat("======================================\n\n")

# 1. 读取所有模型的变量重要性
cat("步骤 1/4: 读取变量重要性数据...\n")

# Maxnet
maxnet_imp <- read.csv("output/05_model_maxnet/variable_importance.csv") %>%
  mutate(model = "Maxnet") %>%
  select(model, variable, importance)

# RF
rf_imp <- read.csv("output/06_model_rf/variable_importance.csv") %>%
  mutate(model = "RF") %>%
  select(model, variable, importance)

# GAM
gam_imp <- read.csv("output/07_model_gam/variable_importance.csv") %>%
  mutate(model = "GAM") %>%
  select(model, variable, importance) %>%
  filter(variable != "lon,lat")  # 排除空间项

cat("  ✓ Maxnet: ", nrow(maxnet_imp), " 个变量\n", sep = "")
cat("  ✓ RF: ", nrow(rf_imp), " 个变量\n", sep = "")
cat("  ✓ GAM: ", nrow(gam_imp), " 个变量\n", sep = "")

# 合并所有模型
all_imp <- bind_rows(maxnet_imp, rf_imp, gam_imp)

# 2. 标准化变量重要性（按模型）
cat("\n步骤 2/4: 标准化变量重要性...\n")

all_imp <- all_imp %>%
  group_by(model) %>%
  mutate(importance_normalized = (importance - min(importance, na.rm = TRUE)) / 
         (max(importance, na.rm = TRUE) - min(importance, na.rm = TRUE) + 1e-10)) %>%
  ungroup()

# 识别在所有模型中都存在的变量
common_vars <- all_imp %>%
  group_by(variable) %>%
  summarise(n_models = n_distinct(model), .groups = "drop") %>%
  filter(n_models == 3) %>%
  pull(variable)

cat("  - 共同变量数: ", length(common_vars), "\n", sep = "")

# 计算平均重要性，选择Top N（避免拥挤，突出主图）
TOP_N <- 20
top_vars <- all_imp %>%
  filter(variable %in% common_vars) %>%
  group_by(variable) %>%
  summarise(mean_importance = mean(importance_normalized, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mean_importance)) %>%
  head(TOP_N) %>%
  pull(variable)

cat("  - 绘图变量数: ", length(top_vars), "\n", sep = "")

# 筛选数据
plot_data <- all_imp %>%
  filter(variable %in% top_vars)

# 3. 绘制分模型“棒棒图”（lollipop），提高可读性（每模型单图 + 汇总图）
cat("\n步骤 3/4: 绘制分模型 lollipop 图...\n")

# 调色（每模型一色，简洁克制）
model_colors <- c("Maxnet" = "#1F78B4", "RF" = "#33A02C", "GAM" = "#E31A1C")

make_lollipop <- function(df_all, model_name, color, top_n = TOP_N) {
  # 中文注释：按单模型排序取 Top N，使用 0→importance 的线段 + 端点强调
  dfm <- df_all %>%
    dplyr::filter(model == model_name) %>%
    dplyr::group_by(variable) %>%
    dplyr::summarise(imp = mean(importance_normalized, na.rm = TRUE), .groups = "drop") %>%
    dplyr::arrange(imp)
  if(nrow(dfm) == 0) return(NULL)
  if(nrow(dfm) > top_n) dfm <- dfm %>% dplyr::slice((nrow(dfm)-top_n+1):n())
  dfm$variable <- factor(dfm$variable, levels = dfm$variable)
  p <- ggplot(dfm, aes(x = imp, y = variable)) +
    geom_segment(aes(x = 0, xend = imp, y = variable, yend = variable),
                 color = color, linewidth = 0.6, alpha = 0.85) +
    geom_point(color = color, size = 1.6) +
    scale_x_continuous(limits = c(0, 1), expand = expansion(mult = c(0, 0.02))) +
    labs(title = paste0("Variable Importance - ", model_name), x = "Normalized Importance", y = NULL) +
    viz_theme_nature(base_size = 8) +
    theme(
      plot.title = element_text(size = 9, face = "bold"),
      axis.text.y = element_text(size = 6),
      axis.text.x = element_text(size = 6)
    )
  height_in <- max(3.0, 0.22 * nrow(dfm) + 1.0)
  out_base <- file.path("figures/09_variable_importance", paste0("lollipop_", tolower(model_name)))
  ggsave(paste0(out_base, ".png"), p, width = 4.2, height = height_in, dpi = 1200, bg = "white")
  ggsave(paste0(out_base, ".svg"), p, width = 4.2, height = height_in, bg = "white")
  cat("  ✓ ", model_name, " : ", out_base, ".png/.svg\n", sep = "")
  return(list(plot = p, height = height_in))
}

plot_list <- list()
for(mn in c("Maxnet", "RF", "GAM")) {
  if(mn %in% unique(plot_data$model)) {
    res <- make_lollipop(plot_data, mn, model_colors[[mn]], top_n = TOP_N)
    if(!is.null(res)) plot_list[[length(plot_list)+1]] <- res
  }
}

# 合并为一张“汇总图”（纵向拼接，主图突出、便于阅读对比）
if(length(plot_list) > 0) {
  library(patchwork)
  pp <- plot_list[[1]]$plot
  if(length(plot_list) > 1) {
    for(ii in 2:length(plot_list)) pp <- pp / plot_list[[ii]]$plot
  }
  total_h <- sum(vapply(plot_list, function(x) x$height, numeric(1)))
  ggsave("figures/09_variable_importance/importance_lollipop_by_model.png", pp,
         width = 4.6, height = total_h, dpi = 1200, bg = "white")
  ggsave("figures/09_variable_importance/importance_lollipop_by_model.svg", pp,
         width = 4.6, height = total_h, bg = "white")
  cat("  ✓ 汇总图: figures/09_variable_importance/importance_lollipop_by_model.png/.svg\n")
}

# 4. 绘制热图
cat("\n步骤 4/4: 绘制变量重要性热图...\n")

# 准备热图数据（改用 ggplot2 统一风格）
heatmap_long <- plot_data %>%
  select(model, variable, importance_normalized)

# 行顺序：按变量平均重要性排序（与小提琴一致）
row_order <- plot_data %>%
  group_by(variable) %>%
  summarise(mean_imp = mean(importance_normalized), .groups = "drop") %>%
  arrange(desc(mean_imp)) %>% pull(variable)
heatmap_long$variable <- factor(heatmap_long$variable, levels = rev(row_order))
heatmap_long$model <- factor(heatmap_long$model, levels = c("Maxnet","RF","GAM"))

p_heat <- ggplot(heatmap_long, aes(x = model, y = variable, fill = importance_normalized)) +
  geom_tile(color = "white", linewidth = 0.2) +
  scale_fill_viridis(name = "Normalized\nImportance", option = "C", limits = c(0,1)) +
  labs(title = "Variable Importance Heatmap", x = NULL, y = NULL) +
  viz_theme_nature(base_size = 8) +
  theme(
    plot.title = element_text(size = 9, face = "bold"),
    axis.text.x = element_text(size = 7),
    axis.text.y = element_text(size = 6)
  )

heatmap_height_in <- max(3.0, 0.20 * length(row_order) + 1.0)
ggsave("figures/09_variable_importance/importance_heatmap.png", p_heat,
       width = 3.6, height = heatmap_height_in, dpi = 1200, bg = "white")
ggsave("figures/09_variable_importance/importance_heatmap.svg", p_heat,
       width = 3.6, height = heatmap_height_in, bg = "white")

 

cat("  ✓ 热图: figures/09_variable_importance/importance_heatmap.png / .svg\n")

# ========================= 新增：SHAP 全局/局部分析 ========================= #
cat("\n[新增] 计算 SHAP（全局与局部）...\n")

# 中文注释：为了兼容各模型（Maxnet/RF/GAM/NN），我们采用 fastshap 的模型无关接口，
# 为每个能找到的模型分别计算：
# 1) 全局重要性（|SHAP|的均值）
# 2) 局部SHAP值（样本×变量），便于后续空间映射或案例解释
# 3) 依赖图（SHAP vs 变量值）

# 读取用于建模的样本数据
model_data <- read.csv("output/04_collinearity/collinearity_removed.csv")
exclude_cols <- c("id", "species", "lon", "lat", "source", "presence", "presence.1")
env_vars <- setdiff(names(model_data), exclude_cols)
X_all <- model_data[, env_vars, drop = FALSE]
y_all <- model_data$presence

# 采样以控制计算量（如数据量很大时）
set.seed(20251024)
sample_n <- min(2000, nrow(X_all))
idx <- sample(seq_len(nrow(X_all)), sample_n)
X_sample <- X_all[idx, , drop = FALSE]
y_sample <- y_all[idx]
# 中文注释：若后续模型（如GAM）包含空间项 s(lon,lat)，则需要为预测提供经纬度列
has_lon <- "lon" %in% names(model_data)
has_lat <- "lat" %in% names(model_data)
lonlat_sample <- if(has_lon && has_lat) model_data[idx, c("lon","lat"), drop = FALSE] else NULL

# 已训练模型路径
model_paths <- c(
  Maxnet = "output/05_model_maxnet/model.rds",
  RF = "output/06_model_rf/model.rds",
  GAM = "output/07_model_gam/model.rds",
  NN = "output/05b_model_nn/model.rds"
)

available_models <- names(model_paths)[file.exists(model_paths)]
if(length(available_models) == 0) {
  cat("  ✗ 未发现已训练模型，跳过 SHAP 计算\n")
} else {
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

  shap_summaries <- list()
  for(mn in available_models) {
    cat("  -> 模型 ", mn, " 的 SHAP ...\n", sep = "")
    mdl <- readRDS(model_paths[[mn]])
    pred_fun <- make_pred_fun(mn, mdl)

    # 确定变量集合
    vars <- if(mn == "NN") mdl$vars else env_vars
    vars <- intersect(vars, colnames(X_sample))
    X_used <- X_sample[, vars, drop = FALSE]
    # 若为GAM且存在经纬度，则附加 lon/lat 以满足 s(lon,lat) 预测需求
    if(mn == "GAM" && !is.null(lonlat_sample)) {
      # 避免重复列
      for(cc in c("lon","lat")) {
        if(!(cc %in% colnames(X_used)) && (cc %in% colnames(lonlat_sample))) {
          X_used[[cc]] <- lonlat_sample[[cc]]
        }
      }
    }

    # 计算 SHAP 值（使用调整以减缓特征相关性的影响）
    shap_vals <- fastshap::explain(
      object = mdl,
      X = X_used,
      pred_wrapper = pred_fun,
      nsim = 64,
      adjust = TRUE
    )
    # 中文注释：fastshap 返回矩阵；为便于按变量名索引，转为数据框
    shap_df <- as.data.frame(shap_vals)

    # 保存局部 SHAP 值
    shap_out <- cbind(data.frame(id = idx, presence = y_sample), shap_df)
    write.csv(shap_out, file = file.path("output/09_variable_importance/shap", paste0("shap_values_", tolower(mn), ".csv")), row.names = FALSE)

    # 计算全局重要性并保存
    # 中文注释：全局重要性中排除经纬度列，聚焦环境变量
    shap_df_env <- shap_df[, setdiff(colnames(shap_df), c("lon","lat")), drop = FALSE]
    shap_imp <- data.frame(variable = colnames(shap_df_env), importance = colMeans(abs(as.matrix(shap_df_env)), na.rm = TRUE)) %>%
      arrange(desc(importance))
    shap_imp$model <- mn
    write.csv(shap_imp, file = file.path("output/09_variable_importance/shap", paste0("shap_global_", tolower(mn), ".csv")), row.names = FALSE)

    # 绘制全局重要性条形图（英语标签，Arial，1200dpi）
    p_bar <- ggplot(shap_imp %>% head(30), aes(x = reorder(variable, importance), y = importance)) +
      geom_col(fill = "#4DAF4A") +
      coord_flip() +
      labs(title = paste0("SHAP Global Importance - ", mn), x = "Variable", y = "Mean |SHAP|") +
      theme_minimal(base_family = "Arial", base_size = 7) +
      theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
    ggsave(filename = file.path("figures/09_variable_importance/shap", paste0("shap_global_bar_", tolower(mn), ".png")), p_bar, width = 3, height = 4, dpi = 1200, device = png, type = "cairo-png")
    

    # 依赖图（Top 6）
    dep_vars <- shap_imp$variable[seq_len(min(6, nrow(shap_imp)))]
    for(v in dep_vars) {
      if(!(v %in% colnames(shap_df_env))) next
      df_dep <- data.frame(x = X_used[[v]], shap = shap_df_env[[v]])
      p_dep <- ggplot(df_dep, aes(x = x, y = shap)) +
        geom_point(alpha = 0.3, size = 0.3, color = "#377EB8") +
        geom_smooth(method = "loess", se = TRUE, color = "#E41A1C", size = 0.4) +
        labs(title = paste0("SHAP Dependence - ", mn, ": ", v), x = v, y = "SHAP value") +
        theme_minimal(base_family = "Arial", base_size = 7) +
        theme(plot.title = element_text(size = 9, face = "bold", hjust = 0.5))
      v_sanit <- gsub("[^A-Za-z0-9_]+", "_", v)
      ggsave(filename = file.path("figures/09_variable_importance/shap", paste0("shap_dependence_", tolower(mn), "_", v_sanit, ".png")), p_dep, width = 2.4, height = 2.4, dpi = 1200, device = png, type = "cairo-png")
      
    }

    shap_summaries[[length(shap_summaries) + 1]] <- shap_imp
    rm(shap_vals)
    gc(verbose = FALSE)
  }

  if(length(shap_summaries) > 0) {
    shap_all <- bind_rows(shap_summaries)
    write.csv(shap_all, file.path("output/09_variable_importance/shap", "shap_global_summary.csv"), row.names = FALSE)
    cat("  ✓ SHAP 结果已保存至 output/09_variable_importance/shap/\n")
  } else {
    cat("  ✗ 未产出 SHAP 结果\n")
  }
}

# 保存汇总数据
summary_data <- all_imp %>%
  filter(variable %in% common_vars) %>%
  select(model, variable, importance, importance_normalized) %>%
  arrange(model, desc(importance_normalized))

write.csv(summary_data, "output/09_variable_importance/importance_summary.csv", row.names = FALSE)

# 日志
sink("output/09_variable_importance/processing_log.txt")
cat("变量重要性对比日志\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
cat("模型数量: 3 (Maxnet, RF, GAM)\n")
cat("共同变量数: ", length(common_vars), "\n", sep = "")
cat("绘图变量数: ", length(top_vars), "\n\n", sep = "")
cat("Top 10变量 (按平均重要性):\n")
top10 <- all_imp %>%
  filter(variable %in% common_vars) %>%
  group_by(variable) %>%
  summarise(mean_importance = mean(importance_normalized), .groups = "drop") %>%
  arrange(desc(mean_importance)) %>%
  head(10)
print(top10)
sink()

cat("\n======================================\n")
cat("变量重要性可视化完成\n")
cat("======================================\n\n")

cat("Top 10最重要变量:\n")
for(i in seq_len(min(10, nrow(top10)))) {
  cat(sprintf("  %2d. %s (%.3f)\n", i, top10$variable[i], top10$mean_importance[i]))
}

cat("\n✓ 脚本执行完成!\n\n")
