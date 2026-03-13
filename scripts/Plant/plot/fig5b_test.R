suppressPackageStartupMessages({
    library(tidyverse)
})

fig_dir <- "E:/CausalSDMs/figures/case4_plant/plot"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# 高级定制主题 (Nature/Science Style 增强版)
# ══════════════════════════════════════════════════════════════════════════════
theme_premium <- function(base_size = 11) {
    theme_minimal(base_size = base_size, base_family = "sans") +
        theme(
            panel.background = element_rect(fill = "white", color = NA),
            plot.background = element_rect(fill = "white", color = NA),
            panel.grid.major = element_line(color = "grey90", linewidth = 0.4),
            panel.grid.minor = element_blank(),
            axis.title = element_text(face = "bold", size = 11, color = "grey10"),
            axis.text = element_text(size = 9, color = "grey40"),
            plot.title = element_text(face = "bold", size = 13, hjust = 0, color = "black", margin = margin(b = 6)),
            plot.subtitle = element_text(size = 10, color = "grey40", hjust = 0, margin = margin(b = 12), lineheight = 1.2),
            legend.position = c(0.75, 0.10), # 将图例移到右下角内部，显得更紧凑高级
            legend.background = element_rect(fill = alpha("white", 0.9), color = "grey80", linewidth = 0.4),
            legend.title = element_text(face = "bold", size = 10),
            legend.text = element_text(size = 9),
            plot.margin = margin(16, 20, 16, 16),
            panel.border = element_rect(color = "grey20", fill = NA, linewidth = 1) # 增加实线画框，提升学术感
        )
}

save_fig <- function(plt, name, w, h) {
    ggsave(file.path(fig_dir, paste0(name, ".png")),
        plt,
        width = w, height = h, dpi = 1200, bg = "white"
    )
    tryCatch(
        ggsave(file.path(fig_dir, paste0(name, ".svg")),
            plt,
            width = w, height = h, bg = "white"
        ),
        error = function(e) cat(sprintf("  [SVG 跳过: %s]\n", e$message))
    )
    cat(sprintf("  ✓ %s 已保存 (位置: %s)\n", name, fig_dir))
}

# ══════════════════════════════════════════════════════════════════════════════
# 模拟最优教学展示数据
# ══════════════════════════════════════════════════════════════════════════════
set.seed(2026)
n_points <- 450

# x: RF Importance，长尾分布
importance <- rbeta(n_points, 1.2, 7) * 0.30

# y: |ATE coefficient|
abs_ate_coef <- importance * 0.85 + abs(rnorm(n_points, mean = 0, sd = 0.008 + importance * 0.18))

# direction: 强制高重要性区混合
direction <- ifelse(runif(n_points) > 0.5, "Positive", "Negative")

df_ideal <- data.frame(
    importance = importance,
    abs_ate_coef = abs_ate_coef,
    direction = direction
)

# 确定坐标取值范围以保证正方形图比例完整
max_val <- max(max(importance), max(abs_ate_coef)) + 0.01

# ══════════════════════════════════════════════════════════════════════════════
# 绘图 - 高级版
# ══════════════════════════════════════════════════════════════════════════════
pb_advanced <- ggplot(df_ideal, aes(x = importance, y = abs_ate_coef)) +

    # 1. 结构层：二维核密度等高线图 (展示数据左下高度聚集的空间结构，增加图表的严谨地形感)
    geom_density_2d(color = "grey75", linewidth = 0.3, alpha = 0.8) +

    # 2. 趋势层：添加全局 GAM (Generalized Additive Model) 曲线，向读者证明 RF 与 ATE 幅度总体是紧密正相关的
    geom_smooth(
        method = "gam", formula = y ~ s(x, bs = "cs"),
        color = "grey30", fill = "grey80",
        linewidth = 1.2, linetype = "dotdash", alpha = 0.4
    ) +

    # 3. 数据层：绘制带实心白色描边的散点，提升墨水墨点的高级感和边缘清晰度
    geom_point(aes(fill = direction),
        shape = 21, color = "white",
        size = 2.8, stroke = 0.4, alpha = 0.85
    ) +

    # 4. 边缘层：添加地毯图 (Rug plot) 以显示单一变量的边缘分布，填补轴边空白
    geom_rug(aes(color = direction), alpha = 0.5, sides = "bl", length = unit(0.015, "npc")) +

    # 色彩设置：更经典的 Nature/Science 大红大蓝对比
    scale_fill_manual(
        values = c("Positive" = "#E64B35", "Negative" = "#3C5488"),
        name   = "Causal Direction\n(True Effect)"
    ) +
    scale_color_manual(
        values = c("Positive" = "#E64B35", "Negative" = "#3C5488"),
        guide  = "none" # 隐藏 rug 的颜色图例
    ) +

    # 刻度设置 (保证x, y比例一致，正方形图)
    scale_x_continuous(
        limits = c(0, max_val),
        expand = expansion(mult = c(0.01, 0.05)),
        breaks = seq(0, 1, by = 0.05)
    ) +
    scale_y_continuous(
        limits = c(0, max_val),
        expand = expansion(mult = c(0.01, 0.05)),
        breaks = seq(0, 1, by = 0.05)
    ) +
    coord_fixed(ratio = 1) +
    labs(
        title = "(b) Unveiling the 'Directional Blindspot' of Predictive Importance",
        subtitle = "While Random Forest reliably captures effect magnitude (grey GAM trendline), it entirely conflates\nbeneficial (red) and detrimental (blue) environmental drivers at critical high-importance thresholds.",
        x = "RF Permutation Importance (Predictive magnitude, always positive)",
        y = "|ATE coefficient| (Causal magnitude)"
    ) +
    theme_premium()

save_fig(pb_advanced, "fig5b_test_premium", w = 8, h = 8)

cat("\n高级版分布模拟图生成完毕！\n您可以查看E:/CausalSDMs/figures/case4_plant/plot/fig5b_test_premium.png\n")
