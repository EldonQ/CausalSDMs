#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 11b_ensemble_prediction_map.R
# 功能说明: 生成2模型集成(Maxnet, RF)的物种分布预测地图 (Publication Quality)
# 规格: Arial字体, 2400 DPI, 透明背景, 无网格线, 透明图例, 无重叠
# 输入文件: output/11_prediction_maps/rasters/pred_*_river.tif
# 输出文件: figures/11_prediction_maps/combined_ensemble.png
# 作者: Nature级别科研项目
# 日期: 2025-12-17
# ==============================================================================

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

# 加载必要的包
packages <- c(
    "raster", "terra", "ggplot2", "sf", "rnaturalearth",
    "viridis", "sysfonts", "showtext", "svglite"
)
for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, dependencies = TRUE)
        library(pkg, character.only = TRUE)
    }
}

# 注册Arial字体
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

dir.create("figures/11_prediction_maps", showWarnings = FALSE, recursive = TRUE)

cat("\n======================================\n")
cat("集成预测地图生成 (Publication Quality)\n")
cat("======================================\n\n")

# 1. 读取两个模型的河网预测栅格
cat("步骤 1/4: 读取模型预测栅格...\n")

model_files <- c(
    "output/11_prediction_maps/rasters/pred_maxnet_river.tif",
    "output/11_prediction_maps/rasters/pred_rf_river.tif"
)

# 检查文件存在性
for (f in model_files) {
    if (!file.exists(f)) {
        stop("文件不存在: ", f, "\n请先运行 11_current_prediction_maps.R")
    }
}

r_maxnet <- raster(model_files[1])
r_rf <- raster(model_files[2])

cat("  ✓ 已读取2个模型预测栅格\n")

# 2. 计算等权重集成平均
cat("\n步骤 2/4: 计算集成预测 (等权重平均)...\n")

# 堆叠并计算均值
stack_all <- stack(r_maxnet, r_rf)
ensemble_mean <- calc(stack_all, fun = mean, na.rm = TRUE)

# 保存集成栅格
writeRaster(ensemble_mean,
    "output/11_prediction_maps/rasters/pred_ensemble_river.tif",
    overwrite = TRUE
)
cat("  ✓ 集成栅格已保存\n")

# 3. 准备绘图数据
cat("\n步骤 3/4: 准备绘图数据...\n")

# 转换为数据框
df <- as.data.frame(ensemble_mean, xy = TRUE)
colnames(df) <- c("lon", "lat", "prob")
df <- df[!is.na(df$prob), ]

# 限制到合理分位数范围避免极端值
q_low <- quantile(df$prob, 0.01, na.rm = TRUE)
q_high <- quantile(df$prob, 0.99, na.rm = TRUE)
df$prob_clamp <- pmax(pmin(df$prob, q_high), q_low)

# 读取中国边界
china <- ne_countries(country = "China", scale = "medium", returnclass = "sf")

cat("  ✓ 数据点数: ", nrow(df), "\n", sep = "")

# 4. 绘制出版质量地图
cat("\n步骤 4/4: 绘制出版质量地图...\n")

# Publication-quality theme
# Publication-quality theme
theme_publication <- theme_minimal(base_size = 12, base_family = "Arial") +
    theme(
        # 标题
        plot.title = element_text(
            face = "bold", size = 14, hjust = 0.5,
            margin = margin(b = 10)
        ),
        # 坐标轴
        axis.title = element_text(face = "bold", size = 10),
        axis.text = element_text(color = "black", size = 9),
        axis.ticks = element_line(color = "black", linewidth = 0.3),
        # 无网格线
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        # 边框
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
        # 透明背景
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA),
        # 图例 - 完全透明
        legend.background = element_rect(fill = NA, color = NA),
        legend.key = element_rect(fill = NA, color = NA),
        legend.box.background = element_rect(fill = NA, color = NA),
        legend.title = element_text(face = "bold", size = 10),
        legend.text = element_text(size = 9),
        legend.position = c(0.88, 0.20),
        legend.key.height = unit(0.5, "cm"),
        legend.key.width = unit(0.3, "cm"),
        legend.spacing.y = unit(0.1, "cm"),
        # 边距
        plot.margin = margin(t = 10, r = 10, b = 10, l = 10)
    )

# 创建地图
p <- ggplot() +
    # 河网预测热图 (使用Spectral色系: 蓝-黄-红)
    geom_point(
        data = df, aes(x = lon, y = lat, color = prob_clamp),
        size = 0.05, shape = 15, alpha = 0.8
    ) +
    # 颜色标度 - Spectral
    scale_color_gradientn(
        colors = rev(RColorBrewer::brewer.pal(11, "Spectral")),
        name = "Suitability",
        limits = c(0, 1),
        breaks = seq(0, 1, 0.2),
        guide = guide_colorbar(
            title.position = "top",
            title.hjust = 0.5,
            barwidth = 0.8,
            barheight = 4, # 缩短图例条
            frame.colour = NA, # 无边框
            ticks.colour = "black",
            ticks.linewidth = 0.5
        )
    ) +
    # 中国边界
    geom_sf(data = china, fill = NA, color = "gray20", linewidth = 0.4) +
    # 坐标设置
    coord_sf(xlim = c(73, 136), ylim = c(17, 55), expand = FALSE) +
    # 标签
    labs(
        title = "Ensemble Prediction",
        x = "Longitude",
        y = "Latitude"
    ) +
    # 主题
    theme_publication

# PNG输出 (2400 DPI, 透明背景)
png("figures/11_prediction_maps/combined_ensemble.png",
    width = 7, height = 5.5, units = "in", res = 2400,
    bg = "transparent", family = "Arial"
)
print(p)
dev.off()

# SVG输出
svglite("figures/11_prediction_maps/combined_ensemble.svg",
    width = 7, height = 5.5, bg = "transparent"
)
print(p)
dev.off()

cat("  ✓ PNG: figures/11_prediction_maps/combined_ensemble.png\n")
cat("  ✓ SVG: figures/11_prediction_maps/combined_ensemble.svg\n")

# 统计摘要
cat("\n======================================\n")
cat("集成预测统计\n")
cat("======================================\n")
cat(sprintf("  河网像元数: %d\n", nrow(df)))
cat(sprintf("  均值: %.3f\n", mean(df$prob)))
cat(sprintf("  标准差: %.3f\n", sd(df$prob)))
cat(sprintf("  范围: %.3f - %.3f\n", min(df$prob), max(df$prob)))

cat("\n✓ 脚本执行完成!\n\n")
