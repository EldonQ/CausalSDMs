#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 01b_variable_prescreening.R
# 功能说明: 对原始栅格变量进行质量检查（Quality Check）
# 策略: 宽松检查，只移除明显有问题的变量，真正筛选在提取后进行
# 输入文件: earthenvstreams_china/*.tif (河流网络环境变量)
# 输出文件: output/01b_variable_prescreening/qualified_variables.csv
# 作者: Nature级别科研项目
# 日期: 2025-10-20
# ==============================================================================

# 清空环境
rm(list = ls())
gc()

# 设置工作目录
setwd("E:/SDM01")

# 加载必要的包（加入showtext/sysfonts以确保PDF/PNG均可使用Arial字体）
packages <- c("raster", "tidyverse", "sf", "showtext", "sysfonts")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 创建输出目录
if(!dir.exists("output")) dir.create("output")
if(!dir.exists("output/01b_variable_prescreening")) {
  dir.create("output/01b_variable_prescreening", recursive = TRUE)
}
if(!dir.exists("figures")) dir.create("figures")
if(!dir.exists("figures/01b_variable_prescreening")) {
  dir.create("figures/01b_variable_prescreening", recursive = TRUE)
}

## 字体配置：强制使用系统 Arial，并启用 showtext 以嵌入/描边字体
## 这样可避免 pdf() 报错 “没有 'Arial' 这样的系列”，并满足期刊字体要求
suppressWarnings({
  if(.Platform$OS.type == "windows") {
    ar_path <- "C:/Windows/Fonts/arial.ttf"
    ar_b_path <- "C:/Windows/Fonts/arialbd.ttf"
    ar_i_path <- "C:/Windows/Fonts/ariali.ttf"
    ar_bi_path <- "C:/Windows/Fonts/arialbi.ttf"
    if(file.exists(ar_path)) {
      sysfonts::font_add(family = "Arial",
                         regular = ar_path,
                         bold = if(file.exists(ar_b_path)) ar_b_path else ar_path,
                         italic = if(file.exists(ar_i_path)) ar_i_path else ar_path,
                         bolditalic = if(file.exists(ar_bi_path)) ar_bi_path else ar_path)
    }
  }
})
showtext::showtext_auto(enable = TRUE)
showtext::showtext_opts(dpi = 1200)

cat("======================================\n")
cat("环境变量质量检查\n")
cat("策略: 宽松检查，真正筛选在提取后\n")
cat("======================================\n\n")

# ------------------------------------------------------------------------------
# 0. 若存在手工白名单（47变量），则直接输出并结束（严格变量模式）
# ------------------------------------------------------------------------------
override_path <- "scripts/variables_selected_47.csv"
if (file.exists(override_path)) {
  cat("检测到变量白名单: ", override_path, "\n", sep = "")
  cat("启用严格变量模式：仅使用这47个变量，直接生成合格变量列表并结束。\n\n")
  qualified_vars_override <- read.csv(override_path, stringsAsFactors = FALSE)
  # 保存为标准输出文件，供下游脚本统一读取
  if(!dir.exists("output/01b_variable_prescreening")) {
    dir.create("output/01b_variable_prescreening", recursive = TRUE)
  }
  write.csv(qualified_vars_override,
            "output/01b_variable_prescreening/qualified_variables.csv",
            row.names = FALSE)
  cat("  ✓ 已保存: output/01b_variable_prescreening/qualified_variables.csv\n")
  # 为记录可重复性，另存一份变量清单到 figures 同名txt（便于审稿材料）
  if(!dir.exists("output/variable_lists")) dir.create("output/variable_lists", recursive = TRUE)
  write.csv(qualified_vars_override, "output/variable_lists/selected_variables_47.csv", row.names = FALSE)
  cat("  ✓ 已保存: output/variable_lists/selected_variables_47.csv\n\n")
  cat("严格变量模式完成，退出本脚本。\n")
  quit(save = "no", status = 0)
}

# ------------------------------------------------------------------------------
# 1. 读取中国边界（用于计算研究区域内的统计量）
# ------------------------------------------------------------------------------
cat("步骤 1/6: 读取研究区域...\n")

china_boundary <- st_read("earthenvstreams_china/china_boundary.shp", quiet = TRUE)
cat("  ✓ 已加载中国边界\n")

# ------------------------------------------------------------------------------
# 2. 定义所有可用的变量配置
# ------------------------------------------------------------------------------
cat("\n步骤 2/6: 定义变量配置...\n")

# 完整的变量配置（包含所有可能的变量）
variable_config <- list(
  # 1. 地形变量 - elevation.tif (4个bands)
  elevation = list(
    file = "elevation.tif",
    total_bands = 4,
    band_names = c("dem_min", "dem_max", "dem_range", "dem_avg"),
    category = "地形",
    scale = 1
  ),
  
  # 2. 坡度变量 - slope.tif (4个bands)
  slope = list(
    file = "slope.tif",
    total_bands = 4,
    band_names = c("slope_min", "slope_max", "slope_range", "slope_avg"),
    category = "地形",
    scale = 0.01  # 转换为度
  ),
  
  # 3. 水流累积 - flow_acc.tif (2个bands)
  flow = list(
    file = "flow_acc.tif",
    total_bands = 2,
    band_names = c("flow_length", "flow_acc"),
    category = "水文",
    scale = 1
  ),
  
  # 4. 水文气候 - hydroclim_average+sum.tif (19个bands)
  hydroclim_avg = list(
    file = "hydroclim_average+sum.tif",
    total_bands = 19,
    band_names = sprintf("hydro_avg_%02d", 1:19),
    category = "水文气候",
    scale = 1,
    temp_bands = 1:11,  # 需要除以10的bands
    prec_bands = 12:19  # 降水bands
  ),
  
  # 5. 水文气候加权 - hydroclim_weighted_average+sum.tif (19个bands)
  hydroclim_wavg = list(
    file = "hydroclim_weighted_average+sum.tif",
    total_bands = 19,
    band_names = sprintf("hydro_wavg_%02d", 1:19),
    category = "水文气候",
    scale = 1,
    temp_bands = 1:11,
    prec_bands = 12:19
  ),
  
  # 6. 月度最低温-平均 - monthly_tmin_average.tif (12个bands)
  tmin_avg = list(
    file = "monthly_tmin_average.tif",
    total_bands = 12,
    band_names = sprintf("tmin_avg_%02d", 1:12),
    category = "月度温度",
    scale = 0.1
  ),
  
  # 7. 月度最低温-加权 - monthly_tmin_weighted_average.tif (12个bands)
  tmin_wavg = list(
    file = "monthly_tmin_weighted_average.tif",
    total_bands = 12,
    band_names = sprintf("tmin_wavg_%02d", 1:12),
    category = "月度温度",
    scale = 0.1
  ),
  
  # 8. 月度最高温-平均 - monthly_tmax_average.tif (12个bands)
  tmax_avg = list(
    file = "monthly_tmax_average.tif",
    total_bands = 12,
    band_names = sprintf("tmax_avg_%02d", 1:12),
    category = "月度温度",
    scale = 0.1
  ),
  
  # 9. 月度最高温-加权 - monthly_tmax_weighted_average.tif (12个bands)
  tmax_wavg = list(
    file = "monthly_tmax_weighted_average.tif",
    total_bands = 12,
    band_names = sprintf("tmax_wavg_%02d", 1:12),
    category = "月度温度",
    scale = 0.1
  ),
  
  # 10. 月度降水-总和 - monthly_prec_sum.tif (12个bands)
  prec_sum = list(
    file = "monthly_prec_sum.tif",
    total_bands = 12,
    band_names = sprintf("prec_sum_%02d", 1:12),
    category = "月度降水",
    scale = 1
  ),
  
  # 11. 月度降水-加权总和 - monthly_prec_weighted_sum.tif (12个bands)
  prec_wsum = list(
    file = "monthly_prec_weighted_sum.tif",
    total_bands = 12,
    band_names = sprintf("prec_wsum_%02d", 1:12),
    category = "月度降水",
    scale = 1
  ),
  
  # 12. 土地覆盖-平均 - landcover_average.tif (12个bands)
  lc_avg = list(
    file = "landcover_average.tif",
    total_bands = 12,
    band_names = sprintf("lc_avg_%02d", 1:12),
    category = "土地覆盖",
    scale = 1
  ),
  
  # 13. 土地覆盖-最大 - landcover_maximum.tif (12个bands)
  lc_max = list(
    file = "landcover_maximum.tif",
    total_bands = 12,
    band_names = sprintf("lc_max_%02d", 1:12),
    category = "土地覆盖",
    scale = 1
  ),
  
  # 14. 土地覆盖-最小 - landcover_minimum.tif (12个bands)
  lc_min = list(
    file = "landcover_minimum.tif",
    total_bands = 12,
    band_names = sprintf("lc_min_%02d", 1:12),
    category = "土地覆盖",
    scale = 1
  ),
  
  # 15. 土地覆盖-范围 - landcover_range.tif (12个bands)
  lc_range = list(
    file = "landcover_range.tif",
    total_bands = 12,
    band_names = sprintf("lc_range_%02d", 1:12),
    category = "土地覆盖",
    scale = 1
  ),
  
  # 16. 土地覆盖-加权平均 - landcover_weighted_average.tif (12个bands)
  lc_wavg = list(
    file = "landcover_weighted_average.tif",
    total_bands = 12,
    band_names = sprintf("lc_wavg_%02d", 1:12),
    category = "土地覆盖",
    scale = 1
  ),
  
  # 17. 土壤-平均 - soil_average.tif (10个bands)
  soil_avg = list(
    file = "soil_average.tif",
    total_bands = 10,
    band_names = sprintf("soil_avg_%02d", 1:10),
    category = "土壤",
    scale = 1
  ),
  
  # 18. 土壤-最大 - soil_maximum.tif (10个bands)
  soil_max = list(
    file = "soil_maximum.tif",
    total_bands = 10,
    band_names = sprintf("soil_max_%02d", 1:10),
    category = "土壤",
    scale = 1
  ),
  
  # 19. 土壤-最小 - soil_minimum.tif (10个bands)
  soil_min = list(
    file = "soil_minimum.tif",
    total_bands = 10,
    band_names = sprintf("soil_min_%02d", 1:10),
    category = "土壤",
    scale = 1
  ),
  
  # 20. 土壤-范围 - soil_range.tif (10个bands)
  soil_range = list(
    file = "soil_range.tif",
    total_bands = 10,
    band_names = sprintf("soil_range_%02d", 1:10),
    category = "土壤",
    scale = 1
  ),
  
  # 21. 土壤-加权平均 - soil_weighted_average.tif (10个bands)
  soil_wavg = list(
    file = "soil_weighted_average.tif",
    total_bands = 10,
    band_names = sprintf("soil_wavg_%02d", 1:10),
    category = "土壤",
    scale = 1
  ),
  
  # 22. 地质-加权总和 - geology_weighted_sum.tif (92个bands)
  geology = list(
    file = "geology_weighted_sum.tif",
    total_bands = 92,
    band_names = sprintf("geo_wsum_%02d", 1:92),
    category = "地质",
    scale = 1
  )
)

# 计算总变量数
total_vars <- sum(sapply(variable_config, function(x) x$total_bands))
cat("  总变量数: ", total_vars, " 个\n", sep = "")
cat("  变量组数: ", length(variable_config), " 组\n", sep = "")

# ------------------------------------------------------------------------------
# 3. 计算每个变量的栅格统计量
# ------------------------------------------------------------------------------
cat("\n步骤 3/6: 计算栅格统计量（这可能需要一些时间）...\n")

env_dir <- "earthenvstreams_china"
var_stats <- data.frame()

for(var_group in names(variable_config)) {
  config <- variable_config[[var_group]]
  file_path <- file.path(env_dir, config$file)
  
  cat("\n  处理: ", config$file, "...\n", sep = "")
  
  if(!file.exists(file_path)) {
    cat("    ✗ 文件不存在，跳过\n")
    next
  }
  
  tryCatch({
    # 读取栅格数据
    r <- brick(file_path)
    
    # 检查实际的band数量
    actual_bands <- nlayers(r)
    if(actual_bands < config$total_bands) {
      cat("    ⚠ 警告: 配置为", config$total_bands, "个bands，实际只有", actual_bands, "个\n")
      config$total_bands <- actual_bands
    }
    
    # 对每个band进行统计
    for(i in 1:config$total_bands) {
      band_name <- config$band_names[i]
      cat("    分析: ", band_name, "... ", sep = "")
      
      # 提取单个band（增加错误处理）
      tryCatch({
        r_band <- r[[i]]
      }, error = function(e) {
        cat("✗ Band访问错误: ", e$message, "\n", sep = "")
        next
      })
      
      # 计算统计量（只采样有效像元）
      # 河流网络数据只在河流位置有值，采样时跳过NA
      sample_vals <- tryCatch({
        sampleRandom(r_band, size = 10000, na.rm = TRUE)
      }, error = function(e) {
        values(r_band)
      })
      
      # 清理NoData标记值
      sample_vals[sample_vals %in% c(-127, -999, -9999) | sample_vals < -1000] <- NA
      sample_vals <- na.omit(sample_vals)
      
      # 应用单位转换（如果需要）
      if(!is.null(config$temp_bands) && i %in% config$temp_bands) {
        sample_vals <- sample_vals / 10
      } else if(!is.null(config$scale) && config$scale != 1) {
        sample_vals <- sample_vals * config$scale
      }
      
      # 计算统计量
      # 注意：sample_vals已经是过滤后的有效值
      n_valid <- length(sample_vals)
      n_total <- ncell(r_band)  # 栅格总像元数
      
      # 对于河流网络数据，缺失值比例是正常的（非河流区域都是NA）
      # 这里我们关注的是有效值的质量，而不是缺失值比例
      missing_pct <- (n_total - n_valid) / n_total * 100
      
      if(n_valid >= 10) {  # 至少需要10个有效值才能统计
        var_mean <- mean(sample_vals, na.rm = TRUE)
        var_sd <- sd(sample_vals, na.rm = TRUE)
        var_min <- min(sample_vals, na.rm = TRUE)
        var_max <- max(sample_vals, na.rm = TRUE)
        var_range <- var_max - var_min
        
        # 计算变异系数（标准差/均值的绝对值，避免除以0）
        if(abs(var_mean) > 1e-10) {
          cv <- abs(var_sd / var_mean)
        } else {
          cv <- NA
        }
        
        # 计算零值比例（对于地质等分类变量）
        zero_pct <- sum(sample_vals == 0, na.rm = TRUE) / n_valid * 100
        
        # 计算标准化范围（max-min）/sd，判断是否有异常值
        if(!is.na(var_sd) && var_sd > 1e-10) {
          normalized_range <- var_range / var_sd
        } else {
          normalized_range <- NA
        }
        
      } else {
        # 有效值太少，标记为NA
        var_mean <- NA
        var_sd <- NA
        var_min <- NA
        var_max <- NA
        var_range <- NA
        cv <- NA
        zero_pct <- NA
        normalized_range <- NA
        cat("有效值不足 (", n_valid, ")\n", sep = "")
      }
      
      # 保存统计结果
      var_stats <- rbind(var_stats, data.frame(
        variable = band_name,
        file = config$file,
        band = i,
        category = config$category,
        n_valid = n_valid,
        n_total = n_total,
        missing_pct = round(missing_pct, 2),
        mean = round(var_mean, 4),
        sd = round(var_sd, 4),
        min = round(var_min, 4),
        max = round(var_max, 4),
        range = round(var_range, 4),
        cv = round(cv, 4),
        zero_pct = round(zero_pct, 2),
        normalized_range = round(normalized_range, 2),
        stringsAsFactors = FALSE
      ))
      
      cat("完成\n")
    }
    
    rm(r)
    gc(verbose = FALSE)
    
  }, error = function(e) {
    cat("    ✗ 错误: ", e$message, "\n", sep = "")
  })
}

cat("\n  ✓ 已分析 ", nrow(var_stats), " 个变量\n", sep = "")

# 保存完整统计结果
write.csv(var_stats,
          "output/01b_variable_prescreening/all_variables_stats.csv",
          row.names = FALSE)
cat("  ✓ 已保存: output/01b_variable_prescreening/all_variables_stats.csv\n")

# ------------------------------------------------------------------------------
# 4. 应用筛选标准
# ------------------------------------------------------------------------------
cat("\n步骤 4/6: 应用筛选标准...\n")

# 定义质量检查标准（极宽松，只移除明显问题）
# 策略：宁可多留，不可错删。真正的变量筛选在提取后进行（04_collinearity_analysis.R）
cat("\n  质量检查标准（宽松模式）:\n")
cat("    1. 有效值数量 >= 50 (数据量不能太少)\n")
cat("    2. 标准差 > 0 (不能是完全常数)\n")
cat("    3. 零值比例 < 99% (不能几乎为空)\n")
cat("    ※ 真正的变量筛选在提取到发生点后进行\n\n")

# 应用宽松的质量检查标准
var_stats$passed_qc <- (
  var_stats$n_valid >= 50 &                              # 至少50个有效值
  !is.na(var_stats$sd) & var_stats$sd > 0 &             # 有变异
  (is.na(var_stats$zero_pct) | var_stats$zero_pct < 99) # 不是空数据
)

# 统计检查结果
n_total <- nrow(var_stats)
n_passed <- sum(var_stats$passed_qc, na.rm = TRUE)
n_removed <- n_total - n_passed

cat("  质量检查结果:\n")
cat("    - 原始变量数: ", n_total, "\n", sep = "")
cat("    - 通过检查: ", n_passed, " (", round(n_passed/n_total*100, 1), "%)\n", sep = "")
cat("    - 存在问题: ", n_removed, " (", round(n_removed/n_total*100, 1), "%)\n\n", sep = "")

# 简化的问题统计
if(n_removed > 0) {
  cat("  数据质量问题:\n")
  cat("    - 有效值不足 (<50): ", sum(var_stats$n_valid < 50, na.rm = TRUE), "\n", sep = "")
  cat("    - 完全常数 (sd=0): ", sum(!is.na(var_stats$sd) & var_stats$sd == 0, na.rm = TRUE), "\n", sep = "")
  cat("    - 几乎为空 (>99%零值): ", sum(!is.na(var_stats$zero_pct) & var_stats$zero_pct >= 99, na.rm = TRUE), "\n\n", sep = "")
}

# 按类别统计
cat("  按类别统计:\n")
category_summary <- var_stats %>%
  group_by(category) %>%
  summarise(
    total = n(),
    passed = sum(passed_qc, na.rm = TRUE),
    removed = total - passed
  ) %>%
  arrange(desc(passed))

print(as.data.frame(category_summary))

# ------------------------------------------------------------------------------
# 5. 生成变量列表
# ------------------------------------------------------------------------------
cat("\n步骤 5/6: 生成变量列表...\n")

# 通过质量检查的变量
qualified_vars <- var_stats %>%
  filter(passed_qc) %>%
  arrange(category, variable) %>%
  select(variable, category, file, band, n_valid, mean, sd, cv, zero_pct)

write.csv(qualified_vars,
          "output/01b_variable_prescreening/qualified_variables.csv",
          row.names = FALSE)
cat("  ✓ 已保存: output/01b_variable_prescreening/qualified_variables.csv\n")

# 有问题的变量（如果有）
if(n_removed > 0) {
  problem_vars <- var_stats %>%
    filter(!passed_qc) %>%
    arrange(category, variable) %>%
    select(variable, category, n_valid, sd, zero_pct)
  
  write.csv(problem_vars,
            "output/01b_variable_prescreening/problem_variables.csv",
            row.names = FALSE)
  cat("  ✓ 已保存: output/01b_variable_prescreening/problem_variables.csv\n")
}

# 保存完整统计结果
write.csv(var_stats,
          "output/01b_variable_prescreening/all_variables_stats.csv",
          row.names = FALSE)
cat("  ✓ 已保存: output/01b_variable_prescreening/all_variables_stats.csv\n")

# ------------------------------------------------------------------------------
# 6. 生成可视化报告（简化）
# ------------------------------------------------------------------------------
cat("\n步骤 6/6: 生成可视化报告...\n")

# 图1: 按类别统计（顺序输出PNG与PDF，避免设备叠加；同时加大画布防止边距错误）
category_matrix <- as.matrix(category_summary[, c("passed", "removed")])
rownames(category_matrix) <- category_summary$category

# PNG版本（1200dpi，4.8x3.6英寸）
png("figures/01b_variable_prescreening/category_summary.png", 
    width = 4800, height = 3600, res = 1200, type = "cairo")
par(mar = c(7, 6, 3, 2), family = "Arial")  # 调整边距并使用Arial
barplot(t(category_matrix), beside = TRUE,
        col = c("#2E7D32", "#D32F2F"),
        legend.text = c("Qualified", "Problem"),
        args.legend = list(x = "topright", bty = "n", cex = 0.8),
        ylab = "Number of Variables",
        main = "Quality Check Results by Category",
        las = 2, cex.names = 0.85)
dev.off()

 
cat("  ✓ 已保存: figures/01b_variable_prescreening/category_summary.png\n")

# 图2: 有效值数量 vs 标准差（如果有问题变量）
if(n_removed > 0) {
  # PNG
  png("figures/01b_variable_prescreening/quality_scatter.png",
      width = 4800, height = 3600, res = 1200, type = "cairo")
  par(mar = c(5.5, 6, 3, 2), family = "Arial")
  plot(var_stats$n_valid, log10(var_stats$sd + 1e-10),
       col = ifelse(var_stats$passed_qc, "#2E7D32", "#D32F2F"),
       pch = 16, cex = 0.6,
       xlab = "Number of Valid Values", 
       ylab = "log10(Standard Deviation)",
       main = "Variable Quality Overview")
  abline(v = 50, col = "red", lty = 2)
  legend("bottomright", 
         legend = c("Qualified", "Problem"),
         col = c("#2E7D32", "#D32F2F"), pch = 16, bty = "n", cex = 0.8)
  dev.off()

  
  cat("  ✓ 已保存: figures/01b_variable_prescreening/quality_scatter.png\n")
}

# ------------------------------------------------------------------------------
# 7. 生成处理报告
# ------------------------------------------------------------------------------
cat("\n生成处理报告...\n")

sink("output/01b_variable_prescreening/quality_check_report.txt")
cat("=======================================================\n")
cat("环境变量质量检查报告\n")
cat("=======================================================\n")
cat("处理时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("1. 检查策略\n")
cat("-----------------------------------------------------------\n")
cat("本脚本对原始栅格变量进行质量检查（Quality Check），\n")
cat("只移除明显有数据质量问题的变量。\n\n")
cat("策略：宁可多留，不可错删\n")
cat("原因：河流网络数据的真正筛选需要在提取到发生点后进行\n\n")

cat("2. 检查标准（宽松模式）\n")
cat("-----------------------------------------------------------\n")
cat("  • 有效值数量 >= 50\n")
cat("  • 标准差 > 0 (非完全常数)\n")
cat("  • 零值比例 < 99% (非空数据)\n\n")

cat("3. 检查结果\n")
cat("-----------------------------------------------------------\n")
cat("原始变量数:", n_total, "\n")
cat("通过检查:", n_passed, " (", round(n_passed/n_total*100, 1), "%)\n", sep = "")
cat("存在问题:", n_removed, " (", round(n_removed/n_total*100, 1), "%)\n\n", sep = "")

cat("4. 按类别统计\n")
cat("-----------------------------------------------------------\n")
print(as.data.frame(category_summary))
cat("\n")

if(n_removed > 0) {
  cat("5. 数据质量问题\n")
  cat("-----------------------------------------------------------\n")
  cat("有效值不足 (<50):", sum(var_stats$n_valid < 50, na.rm = TRUE), "个\n")
  cat("完全常数 (sd=0):", sum(!is.na(var_stats$sd) & var_stats$sd == 0, na.rm = TRUE), "个\n")
  cat("几乎为空 (>99%零值):", sum(!is.na(var_stats$zero_pct) & var_stats$zero_pct >= 99, na.rm = TRUE), "个\n\n")
}

cat("6. 下一步工作流程\n")
cat("-----------------------------------------------------------\n")
cat("通过质量检查的", n_passed, "个变量将进入下一步处理：\n\n", sep = "")
cat("  02_env_extraction_and_cleaning.R\n")
cat("    ↓ 提取并清洗：出现点/背景点环境值\n")
cat("  04_collinearity_analysis.R ★严格变量模式已启用时跳过剔除★\n")
cat("    ↓ 基于实际数据进行变量选择\n")
cat("    • 相关性分析 (|r| > 0.7)\n")
cat("    • VIF分析 (VIF > 10)\n")
cat("    • 生态学意义评估\n")
cat("    ↓\n")
cat("  最终建模变量集\n\n")

cat("7. 重要说明\n")
cat("-----------------------------------------------------------\n")
cat("河流网络环境变量的特点：\n")
cat("  • 只在河流位置有值（非河流区域为NoData）\n")
cat("  • 全局统计不代表研究区域特征\n")
cat("  • 真正的变量筛选必须在提取到发生点后进行\n\n")
cat("因此，本脚本只做基础质量检查，不做严格筛选。\n\n")

cat("=======================================================\n")
cat("报告结束\n")
cat("=======================================================\n")
sink()

cat("  ✓ 已保存: output/01b_variable_prescreening/quality_check_report.txt\n")

# ------------------------------------------------------------------------------
# 完成
# ------------------------------------------------------------------------------
cat("\n======================================\n")
cat("质量检查完成！\n")
cat("======================================\n")
cat("通过检查: ", n_passed, " / ", n_total, " (", round(n_passed/n_total*100, 1), "%)\n\n", sep = "")

cat("输出文件:\n")
cat("  • qualified_variables.csv - 通过质量检查的变量\n")
if(n_removed > 0) {
  cat("  • problem_variables.csv - 存在问题的变量\n")
}
cat("  • all_variables_stats.csv - 完整统计信息\n")
cat("  • quality_check_report.txt - 详细报告\n\n")

cat("下一步: 运行 02_env_extraction_and_cleaning.R 提取并清洗环境变量\n")
cat("变量筛选将在 04_collinearity_analysis.R 中进行（严格变量模式将直接保留白名单）\n\n")

cat("脚本执行完成!\n")

