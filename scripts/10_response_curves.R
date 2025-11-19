#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 10_response_curves.R
# 功能说明: 绘制GAM模型的环境变量响应曲线（偏依赖图）
# 方法: 基于GAM平滑项绘制单变量响应
# 输入文件: output/07_model_gam/model.rds
#          output/09_variable_importance/importance_summary.csv
# 输出文件: figures/10_response_curves/response_curves_top10.png
#          figures/10_response_curves/individual/*.png
# 作者: Nature级别科研项目
# 日期: 2025-10-20
# ==============================================================================

# 初始化环境
rm(list = ls())
gc()
setwd("E:/SDM01")

# 加载必要的包
packages <- c("tidyverse", "mgcv", "ggplot2", "gridExtra", "viridis", "iml", "sysfonts", "showtext", "nnet", "randomForest", "maxnet", "cowplot")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 中文注释：注册并启用 Arial 字体，确保 PDF/PNG 输出嵌入或渲染为 Arial（期刊级别）
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

dir.create("output/10_response_curves", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/10_response_curves", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/10_response_curves/individual", showWarnings = FALSE, recursive = TRUE)
dir.create("output/10_response_curves/ale", showWarnings = FALSE, recursive = TRUE)
dir.create("figures/10_response_curves/ale", showWarnings = FALSE, recursive = TRUE)

# 统一可视化工具（Nature风格、Arial、1200dpi）
source("scripts/visualization/viz_utils.R")

cat("\n======================================\n")
cat("GAM响应曲线绘制\n")
cat("======================================\n\n")

# 采用统一主题：viz_theme_nature（小标题/轴字号已统一收敛，防止溢出）

# 1. 读取模型和变量重要性
cat("步骤 1/3: 读取模型和数据...\n")

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
for(nm in all_predictors) {
  v <- train_df[[nm]]
  if(is.numeric(v)) {
    base_row[[nm]] <- stats::median(v, na.rm = TRUE)
  } else {
    lv <- names(sort(table(v), decreasing = TRUE))[1]
    base_row[[nm]] <- if(is.na(lv)) NA else lv
  }
}
base_row <- as.data.frame(base_row, stringsAsFactors = FALSE)

# 2. 绘制单个响应曲线（概率尺度，ggplot）
cat("\n步骤 2/3: 绘制单变量响应曲线...\n")

plot_list <- list()
for(i in seq_along(top_vars)) {
  var <- top_vars[i]
  cat("  - ", var, "\n", sep = "")

  vx <- train_df[[var]]
  vx <- vx[is.finite(as.numeric(vx))]
  if(length(vx) == 0) next
  rng <- stats::quantile(as.numeric(vx), probs = c(0.01, 0.99), na.rm = TRUE)
  x_seq <- seq(rng[1], rng[2], length.out = 200)

  newd <- base_row[rep(1, length(x_seq)), , drop = FALSE]
  newd[[var]] <- x_seq

  pred <- as.numeric(predict(gam_model, newdata = newd, type = "response"))
  dfp <- data.frame(x = x_seq, y = pred)

  p <- ggplot(dfp, aes(x = x, y = y)) +
    geom_line(linewidth = 0.6, color = "black") +
    labs(title = paste0("Response Curve: ", var), x = var, y = "Presence Probability") +
    coord_cartesian(ylim = c(0, 1)) +
    viz_theme_nature(base_size = 8, title_size = 9)

  ggsave(filename = paste0("figures/10_response_curves/individual/", var, ".png"),
         plot = p, width = 2.4, height = 2.4, units = "in", dpi = 1200, bg = "white")

  plot_list[[length(plot_list) + 1]] <- p
}

cat("  ✓ 单变量曲线: figures/10_response_curves/individual/\n")

# 3. 绘制组合图
cat("\n步骤 3/3: 绘制Top 10组合图...\n")

if(length(plot_list) > 0) {
  comb <- cowplot::plot_grid(plotlist = plot_list, ncol = 2, align = "hv")
  ggsave("figures/10_response_curves/response_curves_top10.png",
         plot = comb, width = 4.8, height = 6, units = "in", dpi = 1200, bg = "white")
}

cat("  ✓ 组合图: figures/10_response_curves/response_curves_top10.png\n")

# ========================= 新增：ALE曲线（模型无关解释） ========================= #
cat("\n[新增] 计算 ALE 曲线（模型无关解释）...\n")

# 中文注释：ALE 比 PDP 更稳健地处理相关特征，这里对 Top 变量与现有四类模型
#（Maxnet / RF / GAM / NN）分别计算 ALE，并分别输出高分辨率图与CSV。

# 读取用于建模的数据集，以获取自变量矩阵与响应
model_data <- read.csv("output/04_collinearity/collinearity_removed.csv")
exclude_cols <- c("id", "species", "lon", "lat", "source", "presence", "presence.1")
env_vars <- setdiff(names(model_data), exclude_cols)
X_all <- model_data[, env_vars, drop = FALSE]
y_all <- model_data$presence

# 中文注释：若后续模型（如GAM）包含 s(lon,lat)，需要为预测提供经纬度列
has_lon <- "lon" %in% names(model_data)
has_lat <- "lat" %in% names(model_data)
lonlat_df <- if(has_lon && has_lat) model_data[, c("lon","lat"), drop = FALSE] else NULL

# 读取各模型文件（存在则计算，不存在则跳过）
model_paths <- c(
  Maxnet = "output/05_model_maxnet/model.rds",
  RF = "output/06_model_rf/model.rds",
  GAM = "output/07_model_gam/model.rds",
  NN = "output/05b_model_nn/model.rds"
)

available_models <- names(model_paths)[file.exists(model_paths)]
if(length(available_models) == 0) {
  cat("  ✗ 未发现已训练模型，跳过 ALE 计算\n")
} else {
  # 依据重要性选择变量（与上文Top变量交集，若无则取出现频率最高的前15个）
  imp_path <- "output/09_variable_importance/importance_summary.csv"
  if(file.exists(imp_path)) {
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
  if(length(ale_vars) > 15) ale_vars <- ale_vars[1:15]  # 中文注释：限制绘制数量以控制运行时间

  # 工具：为不同模型提供预测函数（返回概率）
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

  # 计算并输出
  ale_all <- list()
  for(mn in available_models) {
    cat("  -> 模型 ", mn, " 的 ALE ...\n", sep = "")
    mdl <- readRDS(model_paths[[mn]])
    pred_fun <- make_pred_fun(mn, mdl)
    # 若为GAM且存在经纬度，则附加 lon/lat 以满足 s(lon,lat) 预测需求
    data_for_model <- X_all
    if(mn == "GAM" && !is.null(lonlat_df)) {
      # 避免重复列
      for(cc in c("lon","lat")) {
        if(!(cc %in% colnames(data_for_model)) && (cc %in% colnames(lonlat_df))) {
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

    for(v in ale_vars) {
      fe <- iml::FeatureEffect$new(predictor, feature = v, method = "ale", grid.size = 40)
      # 保存CSV（保持英文字段名）
      res <- fe$results
      # 添加模型与变量列
      res$model <- mn
      res$variable <- v
      # 文件名中清理变量名
      v_sanit <- gsub("[^A-Za-z0-9_]+", "_", v)
      write.csv(res, file = file.path("output/10_response_curves/ale", paste0("ale_", tolower(mn), "_", v_sanit, ".csv")), row.names = FALSE)

      # 单图绘制（PNG + PDF，Arial，1200dpi） —— 移除底部rug，统一风格
      plt <- plot(fe)
      try({
        plt$layers <- Filter(function(ly){ !inherits(ly$geom, "GeomRug") }, plt$layers)
      }, silent = TRUE)
      plt <- plt + ggplot2::labs(title = paste0("ALE - ", mn, ": ", v), x = v, y = "ALE of .y") + theme_nature()

      png(file.path("figures/10_response_curves/ale", paste0("ale_", tolower(mn), "_", v_sanit, ".png")),
          width = 2400, height = 2400, res = 1200, type = "cairo-png", family = "Arial")
      print(plt)
      dev.off()

      

      ale_all[[length(ale_all) + 1]] <- res
      rm(fe)
      gc(verbose = FALSE)
    }
  }

  if(length(ale_all) > 0) {
    ale_df <- dplyr::bind_rows(ale_all)
    write.csv(ale_df, "output/10_response_curves/ale/ale_summary.csv", row.names = FALSE)
    cat("  ✓ ALE 结果已保存至 output/10_response_curves/ale/\n")
  } else {
    cat("  ✗ 未产出 ALE 结果\n")
  }
}

# 日志
sink("output/10_response_curves/processing_log.txt")
cat("GAM响应曲线绘制日志\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n", sep = "")
cat("Top 10变量:\n")
print(data.frame(rank = seq_along(top_vars), variable = top_vars))
sink()

cat("\n======================================\n")
cat("响应曲线绘制完成\n")
cat("======================================\n\n")

cat("✓ 脚本执行完成!\n\n")
