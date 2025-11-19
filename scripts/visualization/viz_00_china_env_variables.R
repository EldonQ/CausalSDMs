################################################################################
# 脚本名称：viz_00_china_env_variables.R
# 功能描述：可视化裁剪到中国区域的环境变量（契合02脚本的变量提取）
# 作者：Nature级别科研项目
# 日期：2025-10-11
# 注意：只绘制02_env_extraction.R脚本中实际使用的波段
################################################################################

# 清空环境
rm(list = ls())
gc()

# 设置工作目录
setwd("E:/SDM01")

# 加载必需的R包（包含svglite以输出SVG矢量图）
required_packages <- c("terra", "sf", "tidyverse", "RColorBrewer", "viridis", "cowplot", "scales", "svglite")

for (pkg in required_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE)
    library(pkg, character.only = TRUE)
  }
}

# 日志函数
log_message <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
}

################################################################################
# 第一步：设置绘图参数（Nature期刊标准）
################################################################################

log_message("设置Nature期刊绘图标准...")

# Nature期刊要求
DPI <- 1200  # 分辨率
FONT_FAMILY <- "Arial"  # 无衬线字体
BASE_SIZE <- 12  # 基础字体大小

# 设置ggplot2主题（Nature风格）
theme_nature <- function() {
  theme_minimal(base_size = BASE_SIZE, base_family = FONT_FAMILY) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5),
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 12, face = "bold"),
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      legend.position = "right",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(size = 0.3, color = "gray90"),
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )
}

# 可选：仅绘制指定变量组（例如只看河网flow）
# 将 ENABLED_GROUPS 设为 NULL 则绘制全部
ENABLED_GROUPS <- NULL

# 创建输出文件夹
output_dir <- "figures/00_china_env_variables"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# 数据文件夹
data_dir <- "earthenvstreams_china"

################################################################################
# 第二步：读取中国边界
################################################################################

log_message("读取中国边界...")

china_boundary <- vect(file.path(data_dir, "china_boundary.shp"))
china_sf <- st_as_sf(china_boundary)

log_message("边界读取完成")
log_message("出图范围: 根据中国边界自动确定（最佳效果）")

################################################################################
# 第三步：定义变量配置（覆盖所有文件的全部波段，自动匹配单位与命名）
################################################################################

log_message("定义变量配置（覆盖全部文件与波段）...")

# 土地覆盖12类（与EarthEnv CLC一致，图注使用英文，单位为%）
lc_classes <- c(
  "Evergreen/deciduous needleleaf trees", "Evergreen broadleaf trees",
  "Deciduous broadleaf trees", "Mixed/other trees",
  "Shrubs", "Herbaceous vegetation",
  "Cultivated & managed vegetation", "Regularly flooded shrub/herbaceous",
  "Urban/built-up", "Snow/ice",
  "Barren/sparse vegetation", "Open water"
)

# Hydroclim 19个指标的英文名称（与层顺序一致，单位见注释）
hydro_names <- c(
  "Annual Mean Upstream Temperature (°C)",                 # 1  (÷10)
  "Mean Upstream Diurnal Range (°C)",                      # 2  (÷10)
  "Upstream Isothermality (index, ×100)",                  # 3  (原始×100)
  "Upstream Temperature Seasonality (index, ×100)",        # 4  (原始×100)
  "Max Upstream Temp of Warmest Month (°C)",               # 5  (÷10)
  "Min Upstream Temp of Coldest Month (°C)",               # 6  (÷10)
  "Upstream Temperature Annual Range (°C)",                # 7  (÷10)
  "Mean Upstream Temp of Wettest Quarter (°C)",            # 8  (÷10)
  "Mean Upstream Temp of Driest Quarter (°C)",             # 9  (÷10)
  "Mean Upstream Temp of Warmest Quarter (°C)",            # 10 (÷10)
  "Mean Upstream Temp of Coldest Quarter (°C)",            # 11 (÷10)
  "Annual Upstream Precipitation (mm)",                    # 12
  "Upstream Precipitation of Wettest Month (mm)",          # 13
  "Upstream Precipitation of Driest Month (mm)",           # 14
  "Upstream Precipitation Seasonality (index, ×100)",      # 15 (原始×100)
  "Upstream Precipitation of Wettest Quarter (mm)",        # 16
  "Upstream Precipitation of Driest Quarter (mm)",         # 17
  "Upstream Precipitation of Warmest Quarter (mm)",        # 18
  "Upstream Precipitation of Coldest Quarter (mm)"         # 19
)

# 月份英文缩写
month_names <- month.abb

# 构建完整变量配置（全部文件×全部波段）
variable_config <- list(
  # 地形
  elevation = list(
    file = "elevation.tif",
    bands = 1:4,
    names = c("DEM Minimum (m)", "DEM Maximum (m)", "DEM Range (m)", "DEM Average (m)"),
    group = "Topography",
    color_scheme = "terrain.colors"
  ),
  # 坡度
  slope = list(
    file = "slope.tif",
    bands = 1:4,
    names = c("Slope Minimum (degrees)", "Slope Maximum (degrees)", "Slope Range (degrees)", "Slope Average (degrees)"),
    group = "Topography",
    color_scheme = "YlOrRd"
  ),
  # 水文（河网长度计数/集水区像元计数）
  flow = list(
    file = "flow_acc.tif",
    bands = 1:2,
    names = c("Flow Length (cell count)", "Flow Accumulation (cell count)"),
    group = "Hydrology",
    color_scheme = "Blues"
  ),
  # 上游Bioclim（平均/总和）
  hydroclim_avg = list(
    file = "hydroclim_average+sum.tif",
    bands = 1:19,
    names = hydro_names,
    group = "Hydroclimate",
    color_scheme = "RdYlBu"
  ),
  # 上游Bioclim（加权）
  hydroclim_wavg = list(
    file = "hydroclim_weighted_average+sum.tif",
    bands = 1:19,
    names = paste0(hydro_names, " (weighted)"),
    group = "Hydroclimate",
    color_scheme = "RdYlBu"
  ),
  # 月尺度气候（上游平均）
  tmin_avg = list(
    file = "monthly_tmin_average.tif",
    bands = 1:12,
    names = sprintf("Tmin %s (°C)", month_names),
    group = "Monthly Temperature",
    color_scheme = "RdBu"
  ),
  tmax_avg = list(
    file = "monthly_tmax_average.tif",
    bands = 1:12,
    names = sprintf("Tmax %s (°C)", month_names),
    group = "Monthly Temperature",
    color_scheme = "RdYlBu"
  ),
  prec_sum = list(
    file = "monthly_prec_sum.tif",
    bands = 1:12,
    names = sprintf("Precipitation %s (mm)", month_names),
    group = "Monthly Precipitation",
    color_scheme = "YlGnBu"
  ),
  # 月尺度气候（上游加权）
  tmin_wavg = list(
    file = "monthly_tmin_weighted_average.tif",
    bands = 1:12,
    names = sprintf("Tmin %s (°C, weighted)", month_names),
    group = "Monthly Temperature",
    color_scheme = "RdBu"
  ),
  tmax_wavg = list(
    file = "monthly_tmax_weighted_average.tif",
    bands = 1:12,
    names = sprintf("Tmax %s (°C, weighted)", month_names),
    group = "Monthly Temperature",
    color_scheme = "RdYlBu"
  ),
  prec_wsum = list(
    file = "monthly_prec_weighted_sum.tif",
    bands = 1:12,
    names = sprintf("Precipitation %s (mm, weighted sum)", month_names),
    group = "Monthly Precipitation",
    color_scheme = "YlGnBu"
  ),
  # 土地覆盖（5种统计）
  landcover_min = list(
    file = "landcover_minimum.tif",
    bands = 1:12,
    names = paste0(lc_classes, " (minimum, %)")
    , group = "Land Cover", color_scheme = "YlGn"
  ),
  landcover_max = list(
    file = "landcover_maximum.tif",
    bands = 1:12,
    names = paste0(lc_classes, " (maximum, %)")
    , group = "Land Cover", color_scheme = "YlGn"
  ),
  landcover_range = list(
    file = "landcover_range.tif",
    bands = 1:12,
    names = paste0(lc_classes, " (range, %)")
    , group = "Land Cover", color_scheme = "YlGn"
  ),
  landcover_avg = list(
    file = "landcover_average.tif",
    bands = 1:12,
    names = paste0(lc_classes, " (average, %)")
    , group = "Land Cover", color_scheme = "YlGn"
  ),
  landcover_wavg = list(
    file = "landcover_weighted_average.tif",
    bands = 1:12,
    names = paste0(lc_classes, " (weighted average, %)")
    , group = "Land Cover", color_scheme = "YlGn"
  ),
  # 土壤（5种统计）
  soil_min = list(
    file = "soil_minimum.tif",
    bands = 1:10,
    names = c("Soil organic carbon (g/kg)", "Soil pH (unit)", "Sand (%)", "Silt (%)",
              "Clay (%)", ">2mm coarse fragments (%)", "Cation exchange capacity (cmol/kg)",
              "Bulk density (kg/m³)", "Depth to bedrock (cm)", "Probability of R horizon (%)"),
    group = "Soil", color_scheme = "BrBG"
  ),
  soil_max = list(
    file = "soil_maximum.tif",
    bands = 1:10,
    names = c("Soil organic carbon (g/kg)", "Soil pH (unit)", "Sand (%)", "Silt (%)",
              "Clay (%)", ">2mm coarse fragments (%)", "Cation exchange capacity (cmol/kg)",
              "Bulk density (kg/m³)", "Depth to bedrock (cm)", "Probability of R horizon (%)"),
    group = "Soil", color_scheme = "BrBG"
  ),
  soil_range = list(
    file = "soil_range.tif",
    bands = 1:10,
    names = c("Soil organic carbon (g/kg)", "Soil pH (unit)", "Sand (%)", "Silt (%)",
              "Clay (%)", ">2mm coarse fragments (%)", "Cation exchange capacity (cmol/kg)",
              "Bulk density (kg/m³)", "Depth to bedrock (cm)", "Probability of R horizon (%)"),
    group = "Soil", color_scheme = "BrBG"
  ),
  soil_avg = list(
    file = "soil_average.tif",
    bands = 1:10,
    names = c("Soil organic carbon (g/kg)", "Soil pH (unit)", "Sand (%)", "Silt (%)",
              "Clay (%)", ">2mm coarse fragments (%)", "Cation exchange capacity (cmol/kg)",
              "Bulk density (kg/m³)", "Depth to bedrock (cm)", "Probability of R horizon (%)"),
    group = "Soil", color_scheme = "BrBG"
  ),
  soil_wavg = list(
    file = "soil_weighted_average.tif",
    bands = 1:10,
    names = c("Soil organic carbon (g/kg)", "Soil pH (unit)", "Sand (%)", "Silt (%)",
              "Clay (%)", ">2mm coarse fragments (%)", "Cation exchange capacity (cmol/kg)",
              "Bulk density (kg/m³)", "Depth to bedrock (cm)", "Probability of R horizon (%)"),
    group = "Soil", color_scheme = "BrBG"
  ),
  # 地质（加权计数）
  geology_wsum = list(
    file = "geology_weighted_sum.tif",
    bands = 1:92,
    names = sprintf("Geology weighted count %02d", 1:92),
    group = "Geology", color_scheme = "Spectral"
  ),
  # 质量控制（用于可视化与屏蔽）
  quality_control = list(
    file = "quality_control.tif",
    bands = 1:2,
    names = c("QC Missing cells (flag)", "QC Cells removed (flag)"),
    group = "Quality Control", color_scheme = "Greys"
  )
)

# 若设置了 ENABLED_GROUPS，则仅保留该组对应的配置
if (!is.null(ENABLED_GROUPS)) {
  keep_idx <- sapply(variable_config, function(x) {
    isTRUE(x$group %in% ENABLED_GROUPS)
  })
  variable_config <- variable_config[keep_idx]
  log_message(paste0("本次仅绘制组: ", paste(unique(sapply(variable_config, function(x) x$group)), collapse=", "))) 
}

# 统计总图数
total_plots <- sum(sapply(variable_config, function(x) length(x$bands)))
log_message(paste0("总共需要绘制 ", total_plots, " 个图"))

################################################################################
# 第四步：批量绘图函数
################################################################################

# 绘制单个栅格图层的函数（自动缩放单位 + 质量控制屏蔽 + 双格式导出）
plot_raster_layer <- function(raster_layer, layer_name, color_scheme, china_sf,
                              output_file, file_name, band_idx,
                              qc_rast = NULL, width = 8, height = 6) {
  
  tryCatch({
    
    # -------------------------
    # 质量控制：对齐并屏蔽被剔除像元（cells_removed==1），统计缺失填补数
    # -------------------------
    qc_note <- ""
    if (!is.null(qc_rast)) {
      # 重采样质量控制栅格以匹配当前图层分辨率/范围
      qc_aligned <- tryCatch({
        resample(qc_rast, raster_layer, method = "near")
      }, error = function(e) {
        qc_rast
      })
      cells_missing <- qc_aligned[[1]]  # missing_cells
      cells_removed <- qc_aligned[[2]]  # cells_removed
      # 屏蔽被剔除像元
      raster_layer <- mask(raster_layer, cells_removed, maskvalues = 1)
      # 统计缺失填补与剔除数量
      missing_count <- tryCatch({
        as.numeric(global(ifel(cells_missing == 1, 1, 0), "sum", na.rm = TRUE)[1,1])
      }, error = function(e) NA)
      removed_count <- tryCatch({
        as.numeric(global(ifel(cells_removed == 1, 1, 0), "sum", na.rm = TRUE)[1,1])
      }, error = function(e) NA)
      qc_note <- paste0("QC: filled=", ifelse(is.na(missing_count), "NA", missing_count),
                        ", removed=", ifelse(is.na(removed_count), "NA", removed_count))
    }

    # -------------------------
    # 单位缩放：依据文件与波段类型进行数值缩放（仅用于显示，不改写源数据）
    # -------------------------
    scale_raster_layer <- function(x, file_name, band_idx) {
      # 坡度：度×100 -> 度
      if (file_name == "slope.tif") return(x / 100)
      # 月最低/最高温（平均/加权）：℃×10 -> ℃
      if (file_name %in% c("monthly_tmin_average.tif", "monthly_tmax_average.tif",
                           "monthly_tmin_weighted_average.tif", "monthly_tmax_weighted_average.tif")) {
        return(x / 10)
      }
      # 土壤pH：pH×10 -> pH
      if (file_name %in% c("soil_minimum.tif", "soil_maximum.tif", "soil_range.tif",
                           "soil_average.tif", "soil_weighted_average.tif")) {
        if (band_idx == 2) return(x / 10) else return(x)
      }
      # 上游Bioclim：温度类波段（1,2,5~11）℃×10 -> ℃；03/04/15为指数(×100)不缩放；降水类不缩放
      if (file_name %in% c("hydroclim_average+sum.tif", "hydroclim_weighted_average+sum.tif")) {
        if (band_idx %in% c(1,2,5,6,7,8,9,10,11)) return(x / 10) else return(x)
      }
      # 其他：不缩放
      return(x)
    }
    raster_layer <- scale_raster_layer(raster_layer, file_name, band_idx)

    # 转换为数据框（用于ggplot2）
    raster_df <- as.data.frame(raster_layer, xy = TRUE, na.rm = TRUE)
    colnames(raster_df)[3] <- "value"
    
    # 如果数据为空，跳过
    if (nrow(raster_df) == 0) {
      log_message(paste0("  - 跳过（无数据）: ", layer_name))
      return(FALSE)
    }
    
    # 数据质量检查
    n_valid <- sum(!is.na(raster_df$value))
    data_range <- range(raster_df$value, na.rm = TRUE)
    log_message(paste0("    有效像元: ", n_valid, 
                      ", 值域: [", round(data_range[1], 2), ", ", 
                      round(data_range[2], 2), "]"))
    
    # 处理极端值（使用分位数裁剪）
    q01 <- quantile(raster_df$value, 0.01, na.rm = TRUE)
    q99 <- quantile(raster_df$value, 0.99, na.rm = TRUE)
    
    # 选择配色方案
    if (color_scheme == "terrain.colors") {
      colors <- terrain.colors(100)
    } else if (color_scheme %in% c("Blues", "Greens", "Reds", "Oranges", 
                                     "YlOrRd", "YlGnBu", "YlGn", "RdYlBu", 
                                     "RdBu", "BrBG", "Spectral")) {
      # 对于RColorBrewer调色板，根据数据类型调整方向
      if (color_scheme %in% c("RdBu", "RdYlBu", "BrBG")) {
        colors <- brewer.pal(11, color_scheme)  # 保持原始方向（冷暖色）
      } else {
        colors <- brewer.pal(9, color_scheme)
      }
    } else {
      colors <- viridis(100, option = "viridis")
    }
    
    # 创建ggplot图（图题使用英文；图注加入QC说明）
    p <- ggplot() +
      geom_raster(data = raster_df, aes(x = x, y = y, fill = value)) +
      geom_sf(data = china_sf, fill = NA, color = "black", size = 0.3, alpha = 0.8) +
      scale_fill_gradientn(
        colors = colors, 
        na.value = "transparent",
        name = "",
        limits = c(q01, q99),
        oob = scales::squish  # 将超出范围的值压缩到边界
      ) +
      coord_sf(expand = FALSE) +  # 根据中国边界自动确定范围
      labs(
        title = layer_name,
        subtitle = qc_note,
        x = "Longitude (°E)",
        y = "Latitude (°N)"
      ) +
      theme_nature() +
      theme(
        legend.key.width = unit(0.5, "cm"),
        legend.key.height = unit(1.5, "cm")
      )
    
    # 保存PNG（1200 dpi）
    ggsave(
      filename = paste0(output_file, ".png"),
      plot = p,
      width = width,
      height = height,
      dpi = DPI,
      units = "in",
      bg = "white"
    )
    # 保存SVG（矢量图）
    ggsave(
      filename = paste0(output_file, ".svg"),
      plot = p,
      width = width,
      height = height,
      units = "in",
      bg = "white",
      device = "svg"
    )

    return(TRUE)
    
  }, error = function(e) {
    log_message(paste0("  - 错误: ", e$message))
    return(FALSE)
  })
}

################################################################################
# 第五步：批量处理所有变量
################################################################################

log_message("开始批量绘图...")

# 创建处理记录
plot_log <- data.frame(
  file = character(),
  band = integer(),
  layer_name = character(),
  status = character(),
  stringsAsFactors = FALSE
)

plot_counter <- 0

# 预加载质量控制栅格（如存在）
qc_path <- file.path(data_dir, "quality_control.tif")
qc_rast <- NULL
if (file.exists(qc_path)) {
  log_message("加载质量控制图层 quality_control.tif 用于屏蔽与标注...")
  qc_rast <- rast(qc_path)
}

# 循环处理每组变量
for (var_group in names(variable_config)) {
  
  config <- variable_config[[var_group]]
  file_path <- file.path(data_dir, config$file)
  
  log_message(paste0("\n处理变量组: ", var_group, " (", config$file, ")"))
  
  if (!file.exists(file_path)) {
    log_message(paste0("  - 文件不存在，跳过"))
    next
  }
  
  tryCatch({
    
    # 读取栅格数据
    raster_data <- rast(file_path)
    
    # 循环处理每个指定的波段
    for (i in seq_along(config$bands)) {
      
      band_idx <- config$bands[i]
      layer_name <- config$names[i]
      
      # 提取指定波段
      raster_layer <- raster_data[[band_idx]]
      
      # 生成输出文件名（安全的文件名）
      safe_name <- gsub("[^A-Za-z0-9_]", "_", layer_name)
      safe_name <- gsub("_{2,}", "_", safe_name)  # 移除多余下划线
      output_file <- file.path(output_dir, 
                               paste0(sprintf("%02d", plot_counter + 1), "_", 
                                     var_group, "_", safe_name))
      
      # 绘制图形
      log_message(paste0("  [", plot_counter + 1, "/", total_plots, "] ", layer_name))
      
      success <- plot_raster_layer(
        raster_layer = raster_layer,
        layer_name = layer_name,
        color_scheme = config$color_scheme,
        china_sf = china_sf,
        output_file = output_file,
        file_name = config$file,
        band_idx = band_idx,
        qc_rast = qc_rast
      )
      
      # 记录处理结果
      plot_log <- rbind(plot_log, data.frame(
        file = config$file,
        band = band_idx,
        layer_name = layer_name,
        status = ifelse(success, "Success", "Failed"),
        stringsAsFactors = FALSE
      ))
      
      if (success) {
        plot_counter <- plot_counter + 1
      }
      
      # 清理内存
      rm(raster_layer)
      gc()
    }
    
    # 清理内存
    rm(raster_data)
    gc()
    
  }, error = function(e) {
    log_message(paste0("  - 处理错误: ", e$message))
  })
}

################################################################################
# 第六步：保存处理记录
################################################################################

log_message("\n保存处理记录...")

# 保存绘图日志
log_file <- file.path(output_dir, "plotting_log.csv")
write.csv(plot_log, log_file, row.names = FALSE)
log_message(paste0("绘图日志已保存：", log_file))

# 保存变量清单
variable_summary <- plot_log %>%
  group_by(file) %>%
  summarise(
    total_bands = n(),
    success_count = sum(status == "Success"),
    failed_count = sum(status != "Success")
  )

summary_file <- file.path(output_dir, "variable_summary.csv")
write.csv(variable_summary, summary_file, row.names = FALSE)
log_message(paste0("变量汇总已保存：", summary_file))

################################################################################
# 第七步：生成摘要报告
################################################################################

log_message("========================================")
log_message("可视化完成摘要：")
log_message(paste0("- 总图数：", nrow(plot_log)))
log_message(paste0("- 成功：", sum(plot_log$status == "Success")))
log_message(paste0("- 失败：", sum(plot_log$status != "Success")))
log_message(paste0("- 输出文件夹：", output_dir))
log_message(paste0("- 图片格式：PNG (1200 dpi)"))
log_message(paste0("- 字体：Arial"))
log_message("========================================")

log_message("脚本执行完毕！")

