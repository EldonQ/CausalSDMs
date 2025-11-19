################################################################################
# 脚本名称：00_crop_china_boundary_V2.R
# 功能描述：将earthenvstreams文件夹中的所有nc环境变量裁剪到中国区域（改进版）
# 使用省级行政区合并为中国整体边界
# 日期：2025-10-11
################################################################################

# 清空环境
rm(list = ls())
gc()

# 设置工作目录
setwd("E:/SDM01")

# 加载必需的R包
library(terra)
library(sf)

# 日志函数
log_message <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
}

################################################################################
# 第一步：读取并合并省级行政区为中国边界
################################################################################

log_message("读取省级行政区数据...")

# 尝试多个省级文件
province_files <- c(
  "data-main/vector/2019全国行政区划/省.shp",
  "data-main/vector/china-all/Province/省级行政区.shp"
)

china_boundary <- NULL
for (pf in province_files) {
  if (file.exists(pf)) {
    log_message(paste0("使用文件: ", basename(pf)))
    china_boundary <- vect(pf)
    break
  }
}

if (is.null(china_boundary)) {
  stop("错误：找不到省级行政区文件！")
}

log_message(paste0("读取到 ", nrow(china_boundary), " 个省级行政区"))

# 合并所有省份为单个中国边界
log_message("合并所有省份为中国整体边界...")
china_boundary <- aggregate(china_boundary, dissolve = TRUE)
log_message("✓ 合并完成")

# 确保WGS84坐标系
china_boundary <- project(china_boundary, "EPSG:4326")
log_message("✓ 坐标系已设置为WGS84")

################################################################################
# 第二步：创建输出文件夹
################################################################################

output_dir <- "earthenvstreams_china"
if (dir.exists(output_dir)) {
  log_message("删除旧的输出文件夹...")
  unlink(output_dir, recursive = TRUE)
}
dir.create(output_dir, recursive = TRUE)
log_message(paste0("✓ 创建输出文件夹：", output_dir))

################################################################################
# 第三步：获取所有nc文件列表
################################################################################

nc_files <- list.files("earthenvstreams", pattern = "\\.nc$", full.names = TRUE)
log_message(paste0("找到 ", length(nc_files), " 个nc文件待处理"))

################################################################################
# 第四步：批量裁剪nc文件
################################################################################

log_message("开始批量裁剪...")

# 处理记录
processing_log <- data.frame(
  file_name = character(),
  n_layers = integer(),
  original_size_mb = numeric(),
  cropped_size_mb = numeric(),
  processing_time_sec = numeric(),
  status = character(),
  stringsAsFactors = FALSE
)

# 循环处理每个nc文件
for (nc_file in nc_files) {
  
  file_name <- basename(nc_file)
  log_message(paste0("\n[", which(nc_files == nc_file), "/", length(nc_files), "] ", file_name))
  
  start_time <- Sys.time()
  
  tryCatch({
    
    # 读取nc文件
    log_message("  读取...")
    nc_rast <- rast(nc_file)
    n_layers <- nlyr(nc_rast)
    original_size <- file.info(nc_file)$size / 1024^2
    
    # 第一步：裁剪到边界范围（矩形）
    log_message("  裁剪...")
    nc_cropped <- crop(nc_rast, china_boundary)
    
    # 第二步：掩膜保留边界内的数据
    log_message("  掩膜...")
    nc_cropped <- mask(nc_cropped, china_boundary)
    
    # 检查数据
    test_values <- values(nc_cropped[[1]], na.rm = TRUE)
    n_valid <- sum(!is.na(test_values))
    log_message(paste0("  有效像元数: ", n_valid))
    
    # 保存
    output_file <- file.path(output_dir, gsub("\\.nc$", ".tif", file_name))
    log_message("  保存...")
    writeRaster(nc_cropped, output_file, overwrite = TRUE, gdal = c("COMPRESS=LZW"))
    
    cropped_size <- file.info(output_file)$size / 1024^2
    processing_time <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    
    log_message(paste0("  ✓ 完成 (", round(processing_time, 1), "秒, ",
                      round(original_size, 1), "→", round(cropped_size, 1), " MB)"))
    
    # 记录
    processing_log <- rbind(processing_log, data.frame(
      file_name = file_name,
      n_layers = n_layers,
      original_size_mb = round(original_size, 2),
      cropped_size_mb = round(cropped_size, 2),
      processing_time_sec = round(processing_time, 2),
      status = "Success",
      stringsAsFactors = FALSE
    ))
    
    # 清理内存
    rm(nc_rast, nc_cropped)
    gc()
    
  }, error = function(e) {
    log_message(paste0("  ✗ 错误: ", e$message))
    processing_log <<- rbind(processing_log, data.frame(
      file_name = file_name,
      n_layers = NA,
      original_size_mb = NA,
      cropped_size_mb = NA,
      processing_time_sec = NA,
      status = paste0("Error: ", e$message),
      stringsAsFactors = FALSE
    ))
  })
}

################################################################################
# 第五步：保存处理记录和边界
################################################################################

# 保存处理日志
write.csv(processing_log, file.path(output_dir, "cropping_log.csv"), row.names = FALSE)

# 保存中国边界
writeVector(china_boundary, file.path(output_dir, "china_boundary.shp"), overwrite = TRUE)

# 生成摘要
log_message("\n========================================")
log_message("处理完成！")
log_message(paste0("成功：", sum(processing_log$status == "Success"), " / 总数：", nrow(processing_log)))
log_message(paste0("原始总大小：", round(sum(processing_log$original_size_mb, na.rm = TRUE) / 1024, 2), " GB"))
log_message(paste0("裁剪后总大小：", round(sum(processing_log$cropped_size_mb, na.rm = TRUE) / 1024, 2), " GB"))
log_message(paste0("总处理时间：", round(sum(processing_log$processing_time_sec, na.rm = TRUE) / 60, 2), " 分钟"))
log_message("========================================")

log_message("脚本执行完毕！")


