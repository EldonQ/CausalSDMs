#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 02_env_extraction_and_cleaning.R
# 功能说明: 从中国区域裁剪的环境变量提取数据并处理缺失值
# 策略: 使用白名单变量提取，严格执行单位转换和缺失值填补
# 输入文件:
#   - output/01_data_preparation/species_occurrence_cleaned.csv
#   - output/01b_variable_prescreening/qualified_variables.csv
# 输出文件:
#   - output/02_env_extraction/occurrence_with_env_complete.csv
#   - output/02_env_extraction/extracted_variables.csv
# 作者: Nature级别科研项目
# 日期: 2025-10-20
# ==============================================================================

# 清空环境
rm(list = ls())
gc()

# 设置工作目录
setwd("E:/CausalSDMs")

# 加载必要的包
packages <- c("raster", "sp", "tidyverse")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 创建输出目录
if (!dir.exists("output/02_env_extraction")) {
  dir.create("output/02_env_extraction", recursive = TRUE)
}

cat("======================================\n")
cat("环境变量提取与数据清洗\n")
cat("策略: 依据tif元数据进行单位转换，并确保输出数据完整无缺失\n")
cat("======================================\n\n")

# ------------------------------------------------------------------------------
# 1. 读取清洗后的出现点数据
# ------------------------------------------------------------------------------
cat("步骤 1/7: 读取清洗后的出现点数据...\n")

occurrences <- read.csv("output/01_data_preparation/species_occurrence_cleaned.csv")
cat("  - 记录数: ", nrow(occurrences), "\n", sep = "")
cat("  - 物种数: ", length(unique(occurrences$species)), "\n", sep = "")

# 创建空间点对象
coordinates(occurrences) <- ~ lon + lat
proj4string(occurrences) <- CRS("+proj=longlat +datum=WGS84")
occ_sp <- occurrences

# 准备环境数据框
env_data <- data.frame(
  id = 1:length(occ_sp),
  species = occ_sp$species,
  lon = coordinates(occ_sp)[, 1],
  lat = coordinates(occ_sp)[, 2],
  source = occ_sp$source
)

# ------------------------------------------------------------------------------
# 2. 读取初筛后的合格变量列表
# ------------------------------------------------------------------------------
cat("\n步骤 2/7: 读取初筛后的合格变量列表...\n")

qualified_vars <- read.csv("output/01b_variable_prescreening/qualified_variables.csv")
cat("  - 合格变量数: ", nrow(qualified_vars), "\n", sep = "")

# 按文件组织变量
vars_by_file <- qualified_vars %>%
  arrange(file, band) %>%
  group_by(file) %>%
  summarise(
    bands = list(band),
    var_names = list(variable),
    category = first(category),
    .groups = "drop"
  )

# ------------------------------------------------------------------------------
# 3. 提取环境变量
# ------------------------------------------------------------------------------
cat("\n步骤 3/7: 提取环境变量...\n")

env_dir <- "earthenvstreams_china"
all_var_names <- c()
extraction_log <- list()

for (i in 1:nrow(vars_by_file)) {
  file_name <- vars_by_file$file[i]
  file_path <- file.path(env_dir, file_name)
  bands <- vars_by_file$bands[[i]]
  var_names <- vars_by_file$var_names[[i]]

  cat("  [", i, "/", nrow(vars_by_file), "] 处理: ", file_name, " (", length(bands), " 个变量)...\n", sep = "")

  tryCatch(
    {
      r <- brick(file_path)
      r_selected <- r[[bands]]
      
      # 提取值 (增加5km缓冲以减少边缘NA)
      values <- raster::extract(r_selected, occ_sp, buffer = 5000, fun = mean, na.rm = TRUE)
      if (is.null(dim(values))) values <- matrix(values, ncol = 1)

      # 转换NoData值为NA
      values[values == -127] <- NA
      values[values == -999] <- NA
      values[values == -9999] <- NA
      values[values < -1000] <- NA

      # === 单位转换逻辑 (CRITICAL) ===
      # 1. 温度 (Bioclim 1-11): 原值通常放大10倍存储 -> 除以10还原为℃
      #    hydro_wavg_01 ~ 11 对应 BIO1 ~ BIO11 (温度相关)
      #    hydro_wavg_12 ~ 19 对应 BIO12 ~ BIO19 (降水相关，通常为mm，不需转换)
      #    tmin_, tmax_ 也通常放大10倍
      temp_vars <- grepl("^(tmin_|tmax_|hydro_wavg_0[1-9]|hydro_wavg_1[01])", var_names)
      if (any(temp_vars)) {
        cols <- which(temp_vars)
        values[, cols] <- values[, cols] / 10
      }

      # 2. 坡度 (Slope): 原值通常放大100倍存储 -> 除以100还原为度或百分比
      slope_vars <- grepl("^slope_", var_names)
      if (any(slope_vars)) {
        cols <- which(slope_vars)
        values[, cols] <- values[, cols] / 100
      }

      # 3. 土壤 pH: soil_wavg_02 通常是 pH * 10 -> 除以10还原
      ph_vars <- var_names %in% c("soil_wavg_02")
      if (any(ph_vars)) {
        cols <- which(ph_vars)
        values[, cols] <- values[, cols] / 10
      }
      
      # 设置列名
      colnames(values) <- var_names
      env_data <- cbind(env_data, values)
      all_var_names <- c(all_var_names, var_names)
      
      rm(r, r_selected, values); gc(verbose = FALSE)
      cat("    ✓ 成功提取 (单位已转换)\n")
    },
    error = function(e) {
      cat("    ✗ 错误: ", e$message, "\n", sep = "")
    }
  )
}

# ------------------------------------------------------------------------------
# 4. 数据质量检查与数值范围核对
# ------------------------------------------------------------------------------
cat("\n步骤 4/7: 数据质量与数值范围检查...\n")

# 打印数值范围摘要，供用户核对单位是否正确
cat("  - 关键变量数值范围预览 (请核对单位):\n")
summary_stats <- function(x) {
  sprintf("Min: %.1f, Median: %.1f, Max: %.1f, Mean: %.1f", 
          min(x, na.rm=T), median(x, na.rm=T), max(x, na.rm=T), mean(x, na.rm=T))
}

# 检查几个代表性变量
check_vars <- c("hydro_wavg_01", "hydro_wavg_12", "slope_mean", "soil_wavg_02", "elevation")
for (v in check_vars) {
  if (v %in% names(env_data)) {
    cat("    * ", v, ": ", summary_stats(env_data[[v]]), "\n", sep = "")
  }
}

# 缺失值统计
missing_counts <- colSums(is.na(env_data[, all_var_names]))
missing_pct <- missing_counts / nrow(env_data) * 100
env_data$n_missing <- rowSums(is.na(env_data[, all_var_names]))
env_data$missing_pct <- round(env_data$n_missing / length(all_var_names) * 100, 1)

cat("\n  - 记录级缺失情况:\n")
cat("    * 完整记录: ", sum(env_data$missing_pct == 0), "\n", sep = "")
cat("    * 严重缺失(>10%): ", sum(env_data$missing_pct >= 10), "\n", sep = "")

# ------------------------------------------------------------------------------
# 5. 处理缺失数据 (确保完整性)
# ------------------------------------------------------------------------------
cat("\n步骤 5/7: 处理缺失数据 (确保数据完整)...\n")

# 策略: 
# 1. 剔除缺失严重的记录 (>= 10% 变量缺失)
# 2. 对剩余记录进行中位数插补，填补少量NA

# 1. 剔除
env_clean <- env_data %>% filter(missing_pct < 10)
n_removed <- nrow(env_data) - nrow(env_clean)
cat("  - 移除缺失严重记录: ", n_removed, " 条\n", sep = "")

# 2. 插补
vars_to_impute <- all_var_names[colSums(is.na(env_clean[, all_var_names])) > 0]
if (length(vars_to_impute) > 0) {
  cat("  - 对 ", length(vars_to_impute), " 个含缺失变量进行中位数插补...\n", sep = "")
  for (var in vars_to_impute) {
    median_val <- median(env_clean[[var]], na.rm = TRUE)
    env_clean[[var]][is.na(env_clean[[var]])] <- median_val
  }
  cat("  ✓ 插补完成\n")
} else {
  cat("  - 无需插补\n")
}

# 最终完整性检查
if (anyNA(env_clean[, all_var_names])) {
  stop("处理后仍存在缺失值，请检查插补逻辑！")
} else {
  cat("  ✓ 数据完整性校验通过 (所有记录无缺失值)\n")
}

# 清理辅助列
env_clean$n_missing <- NULL
env_clean$missing_pct <- NULL

# ------------------------------------------------------------------------------
# 6. 保存结果
# ------------------------------------------------------------------------------
cat("\n步骤 6/7: 保存结果...\n")

write.csv(env_clean, "output/02_env_extraction/occurrence_with_env_complete.csv", row.names = FALSE)
cat("  ✓ 已保存: output/02_env_extraction/occurrence_with_env_complete.csv\n")

# 保存变量元数据
var_info <- data.frame(index = 1:length(all_var_names), variable = all_var_names, stringsAsFactors = FALSE) %>%
  left_join(qualified_vars[, c("variable", "category", "file", "band")], by = "variable")
write.csv(var_info, "output/02_env_extraction/extracted_variables.csv", row.names = FALSE)
cat("  ✓ 已保存: output/02_env_extraction/extracted_variables.csv\n")

# ------------------------------------------------------------------------------
# 7. 生成详细日志 (含数值范围)
# ------------------------------------------------------------------------------
sink("output/02_env_extraction/processing_log.txt")
cat("环境变量提取与数据清洗日志\nUpdated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("=== 1. 数据概况 ===\n")
cat("原始记录: ", nrow(env_data), "\n")
cat("清洗后记录: ", nrow(env_clean), " (保留率 ", round(nrow(env_clean)/nrow(env_data)*100, 1), "%)\n")
cat("变量数量: ", length(all_var_names), "\n\n")

cat("=== 2. 数值范围检查 (核对单位) ===\n")
cat("变量名 | Min | Median | Max | Mean | 单位说明\n")
cat("--- | --- | --- | --- | --- | ---\n")
# 遍历所有变量打印统计
for (v in all_var_names) {
  vals <- env_clean[[v]]
  unit_note <- ""
  if (grepl("wavg_0[1-9]|wavg_1[01]|tmin|tmax", v)) unit_note <- "Temperature (℃)"
  if (grepl("slope", v)) unit_note <- "Slope (deg/%)"
  if (grepl("soil_wavg_02", v)) unit_note <- "pH"
  
  cat(sprintf("%s | %.2f | %.2f | %.2f | %.2f | %s\n", 
              v, min(vals), median(vals), max(vals), mean(vals), unit_note))
}

cat("\n=== 3. 处理操作 ===\n")
cat("- 单位转换: 温度/10, 坡度/100, pH/10\n")
cat("- 缺失处理: 移除缺失率>=10%的行，其余进行中位数插补\n")
sink()

cat("\n步骤 7/7: 处理日志已生成 -> output/02_env_extraction/processing_log.txt\n")
cat("\n脚本执行完成！\n")
