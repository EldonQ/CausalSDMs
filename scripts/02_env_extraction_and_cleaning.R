#!/usr/bin/env Rscript
# ==============================================================================
# 脚本名称: 02_env_extraction_and_cleaning.R
# 功能说明: 从中国区域裁剪的环境变量提取数据并处理缺失值
# 策略: 使用白名单的47个变量，一步完成提取和清洗
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
setwd("E:/SDM01")

# 加载必要的包
packages <- c("raster", "sp", "tidyverse")
for(pkg in packages) {
  if(!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 创建输出目录
if(!dir.exists("output/02_env_extraction")) {
  dir.create("output/02_env_extraction", recursive = TRUE)
}

cat("======================================\n")
cat("环境变量提取与数据清洗\n")
cat("策略: 使用白名单/初筛后的合格变量（当前应为47个）\n")
cat("======================================\n\n")

# ------------------------------------------------------------------------------
# 1. 读取清洗后的出现点数据
# ------------------------------------------------------------------------------
cat("步骤 1/7: 读取清洗后的出现点数据...\n")

occurrences <- read.csv("output/01_data_preparation/species_occurrence_cleaned.csv")
cat("  - 记录数: ", nrow(occurrences), "\n", sep = "")
cat("  - 物种数: ", length(unique(occurrences$species)), "\n", sep = "")
cat("  - 来源分布:\n")
print(table(occurrences$source))

# 创建空间点对象
coordinates(occurrences) <- ~lon+lat
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

# 按类别统计
cat("  - 变量类别分布:\n")
category_counts <- table(qualified_vars$category)
for(cat_name in names(category_counts)) {
  cat("    * ", cat_name, ": ", category_counts[cat_name], " 个\n", sep = "")
}

# 按文件组织变量（提高提取效率）
vars_by_file <- qualified_vars %>%
  arrange(file, band) %>%
  group_by(file) %>%
  summarise(
    bands = list(band),
    var_names = list(variable),
    category = first(category),
    .groups = "drop"
  )

cat("  - 需要读取的文件数: ", nrow(vars_by_file), "\n\n", sep = "")

# ------------------------------------------------------------------------------
# 3. 提取环境变量
# ------------------------------------------------------------------------------
cat("步骤 3/7: 提取环境变量...\n")

env_dir <- "earthenvstreams_china"
all_var_names <- c()
extraction_log <- list()

for(i in 1:nrow(vars_by_file)) {
  file_name <- vars_by_file$file[i]
  file_path <- file.path(env_dir, file_name)
  bands <- vars_by_file$bands[[i]]
  var_names <- vars_by_file$var_names[[i]]
  
  cat("  [", i, "/", nrow(vars_by_file), "] 处理: ", file_name, " (", 
      length(bands), " 个变量)...\n", sep = "")
  
  tryCatch({
    # 读取.tif文件
    r <- brick(file_path)
    
    # 提取指定的波段
    r_selected <- r[[bands]]
    
    # 提取值
    values <- raster::extract(r_selected, occ_sp)
    if(is.null(dim(values))) {
      values <- matrix(values, ncol = 1)
    }
    
    # 转换NoData值为NA
    values[values == -127] <- NA
    values[values == -999] <- NA
    values[values == -9999] <- NA
    values[values < -1000] <- NA
    
    # 应用单位转换（根据变量类型）
    # 温度变量除以10（包括加权Bioclim的温度类：hydro_wavg_01~11）
    temp_vars <- grepl("^(tmin_|tmax_|hydro_wavg_0[1-9]|hydro_wavg_1[01])", var_names)
    if(any(temp_vars)) {
      temp_cols <- which(temp_vars)
      values[, temp_cols] <- values[, temp_cols] / 10
    }
    
    # 坡度变量除以100
    slope_vars <- grepl("^slope_", var_names)
    if(any(slope_vars)) {
      slope_cols <- which(slope_vars)
      values[, slope_cols] <- values[, slope_cols] / 100
    }

    # 土壤 pH 变量除以10（soil_wavg_02 为 pH×10 存储）
    ph_vars <- var_names %in% c("soil_wavg_02")
    if(any(ph_vars)) {
      ph_cols <- which(ph_vars)
      values[, ph_cols] <- values[, ph_cols] / 10
    }
    
    # 设置列名
    colnames(values) <- var_names
    
    # 添加到数据框
    env_data <- cbind(env_data, values)
    all_var_names <- c(all_var_names, var_names)
    
    # 记录提取结果
    extraction_log[[file_name]] <- data.frame(
      variable = var_names,
      band = bands,
      n_extracted = nrow(values),
      n_missing = colSums(is.na(values)),
      stringsAsFactors = FALSE
    )
    
    cat("    ✓ 成功提取 ", length(var_names), " 个变量\n", sep = "")
    
    rm(r, r_selected, values)
    gc(verbose = FALSE)
    
  }, error = function(e) {
    cat("    ✗ 错误: ", e$message, "\n", sep = "")
  })
}

cat("\n  总提取变量数: ", length(all_var_names), "\n", sep = "")

# ------------------------------------------------------------------------------
# 4. 数据质量检查
# ------------------------------------------------------------------------------
cat("\n步骤 4/7: 数据质量检查...\n")

# 检查缺失值
missing_counts <- colSums(is.na(env_data[, all_var_names]))
missing_pct <- missing_counts / nrow(env_data) * 100

cat("  - 变量级缺失值统计:\n")
cat("    * 完全无缺失的变量: ", sum(missing_counts == 0), " 个\n", sep = "")
cat("    * 缺失<10%的变量: ", sum(missing_pct < 10), " 个\n", sep = "")
cat("    * 缺失10-30%的变量: ", sum(missing_pct >= 10 & missing_pct < 30), " 个\n", sep = "")
cat("    * 缺失>30%的变量: ", sum(missing_pct >= 30), " 个\n", sep = "")

# 计算每个记录的缺失率
env_data$n_missing <- rowSums(is.na(env_data[, all_var_names]))
env_data$missing_pct <- round(env_data$n_missing / length(all_var_names) * 100, 1)

cat("\n  - 记录级缺失值统计:\n")
cat("    * 完整(0%缺失): ", sum(env_data$missing_pct == 0), " 个点\n", sep = "")
cat("    * 接近完整(<10%缺失): ", sum(env_data$missing_pct > 0 & env_data$missing_pct < 10), " 个点\n", sep = "")
cat("    * 部分缺失(10-50%): ", sum(env_data$missing_pct >= 10 & env_data$missing_pct < 50), " 个点\n", sep = "")
cat("    * 严重缺失(>50%): ", sum(env_data$missing_pct >= 50), " 个点\n", sep = "")

# ------------------------------------------------------------------------------
# 5. 处理缺失数据
# ------------------------------------------------------------------------------
cat("\n步骤 5/7: 处理缺失数据...\n")

# 策略1：移除缺失率≥10%的记录
env_clean <- env_data %>%
  filter(missing_pct < 10)

cat("  - 筛选后记录数: ", nrow(env_clean), " / ", nrow(env_data), 
    " (保留 ", round(nrow(env_clean) / nrow(env_data) * 100, 1), "%)\n", sep = "")

# 策略2：对于缺失<10%的记录，进行中位数插补
vars_to_impute <- all_var_names[colSums(is.na(env_clean[, all_var_names])) > 0]

if(length(vars_to_impute) > 0) {
  cat("  - 需要插补的变量数: ", length(vars_to_impute), "\n", sep = "")
  
  for(var in vars_to_impute) {
    n_missing <- sum(is.na(env_clean[[var]]))
    median_val <- median(env_clean[[var]], na.rm = TRUE)
    env_clean[[var]][is.na(env_clean[[var]])] <- median_val
  }
  cat("  ✓ 插补完成（使用中位数）\n")
} else {
  cat("  - 无需插补，数据已完整\n")
}

# 验证完整性
complete_check <- complete.cases(env_clean[, all_var_names])
cat("  - 完整记录验证: ", sum(complete_check), " / ", nrow(env_clean), 
    ifelse(sum(complete_check) == nrow(env_clean), " ✓", " ✗"), "\n", sep = "")

# 移除辅助列
env_clean$n_missing <- NULL
env_clean$missing_pct <- NULL

# ------------------------------------------------------------------------------
# 6. 保存结果
# ------------------------------------------------------------------------------
cat("\n步骤 6/7: 保存结果...\n")

# 保存完整的清洗后数据
write.csv(env_clean,
          "output/02_env_extraction/occurrence_with_env_complete.csv",
          row.names = FALSE)
cat("  ✓ 已保存: output/02_env_extraction/occurrence_with_env_complete.csv\n")

# 保存变量信息
var_info <- data.frame(
  index = 1:length(all_var_names),
  variable = all_var_names,
  stringsAsFactors = FALSE
) %>%
  left_join(qualified_vars[, c("variable", "category", "file", "band")], 
            by = "variable")

write.csv(var_info,
          "output/02_env_extraction/extracted_variables.csv",
          row.names = FALSE)
cat("  ✓ 已保存: output/02_env_extraction/extracted_variables.csv\n")

# 保存物种统计
species_counts <- env_clean %>%
  group_by(species) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(desc(n))

write.csv(species_counts,
          "output/02_env_extraction/species_counts_clean.csv",
          row.names = FALSE)
cat("  ✓ 已保存: output/02_env_extraction/species_counts_clean.csv\n")

# ------------------------------------------------------------------------------
# 7. 生成处理日志
# ------------------------------------------------------------------------------
cat("\n步骤 7/7: 生成处理日志...\n")

sink("output/02_env_extraction/processing_log.txt")
cat("环境变量提取与数据清洗日志\n")
cat("处理时间:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

cat("=== 数据概览 ===\n")
cat("原始出现点数量: ", nrow(env_data), "\n", sep = "")
cat("清洗后记录数量: ", nrow(env_clean), "\n", sep = "")
cat("保留率: ", round(nrow(env_clean) / nrow(env_data) * 100, 1), "%\n", sep = "")
cat("提取变量数: ", length(all_var_names), "\n", sep = "")
cat("数据完整率: 100%\n\n")

cat("=== 处理策略 ===\n")
cat("1. 使用白名单的47个变量（四组：G1/G2/G3/G4）\n")
cat("2. 单位统一：温度(hydro_wavg_01–11)÷10；坡度÷100；土壤pH(soil_wavg_02)÷10\n")
cat("3. 移除记录级缺失率≥10%的出现点\n")
cat("4. 对剩余<10%缺失的记录进行中位数插补\n")
cat("5. 确保所有记录完整无缺失\n\n")

cat("=== 变量类别分布 ===\n")
var_category_counts <- table(var_info$category)
print(as.data.frame(var_category_counts))

cat("\n=== 空间范围 ===\n")
cat("经度范围: ", range(env_clean$lon), "\n", sep = "")
cat("纬度范围: ", range(env_clean$lat), "\n\n", sep = "")

cat("=== 物种统计 ===\n")
cat("总物种数: ", nrow(species_counts), "\n", sep = "")
cat("记录数≥5的物种: ", sum(species_counts$n >= 5), " 个\n", sep = "")
cat("记录数≥3的物种: ", sum(species_counts$n >= 3), " 个\n", sep = "")
cat("记录数=1的物种: ", sum(species_counts$n == 1), " 个\n\n", sep = "")

cat("前20个物种记录数:\n")
print(head(species_counts, 20))

sink()

cat("  ✓ 已保存: output/02_env_extraction/processing_log.txt\n")

# ------------------------------------------------------------------------------
# 摘要
# ------------------------------------------------------------------------------
cat("\n======================================\n")
cat("环境变量提取与清洗完成\n")
cat("======================================\n")
cat("提取变量数: ", length(all_var_names), "\n", sep = "")
cat("原始记录: ", nrow(env_data), " 个点\n", sep = "")
cat("清洗后记录: ", nrow(env_clean), " 个点 (保留 ",
    round(nrow(env_clean) / nrow(env_data) * 100, 1), "%)\n", sep = "")
cat("数据完整率: 100%\n")
cat("物种数: ", length(unique(env_clean$species)), "\n\n", sep = "")

cat("⚠ 重要说明:\n")
cat("  - 使用了白名单的", length(all_var_names), "个变量（应为47个）\n", sep = "")
cat("  - 移除了", nrow(env_data) - nrow(env_clean), "个缺失严重的记录\n", sep = "")
cat("  - 单位已统一：温度÷10、坡度÷100、土壤pH÷10；所有保留记录均已完整插补\n")
cat("  - 可直接用于后续共线性分析和建模\n\n")

cat("脚本执行完成！\n")

