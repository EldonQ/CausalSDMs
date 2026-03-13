# ══════════════════════════════════════════════════════════════════════════════
# 绘制 CAST 模型所用 32 种哺乳动物的提取分布点地图
# ══════════════════════════════════════════════════════════════════════════════
suppressPackageStartupMessages({
    library(sf)
    library(dplyr)
    library(data.table)
    library(ggplot2)
    library(purrr)
    library(viridis)
})

# 1. 路径配置
data_dir <- "E:/CausalSDMs/outputs/EcoISEA3H/Res9/CAST_ready/species_data_screened"
china_shp_file <- "E:/CausalSDMs/data-main/vector/china.shp"
fig_dir <- "E:/CausalSDMs/figures/EcoISEA3H/Res9/plot"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# 2. 读取底图
cat("Loading China boundary...\n")
sf_use_s2(FALSE)
china_sf <- st_read(china_shp_file, quiet = TRUE)
if (is.na(st_crs(china_sf)) || st_crs(china_sf)$epsg != 4326) {
    china_sf <- st_transform(china_sf, 4326)
}
china_sf <- st_make_valid(china_sf)

# 3. 扫描并合并 32 个物种的筛查后分布数据
cat("Scanning screened species CSVs...\n")
csv_files <- list.files(data_dir, pattern = "\\.csv$", full.names = TRUE)

# 我们只提取有分布的点 (presence == 1) 或按 fraction 进行渲染
# 为了高效，提取 lon, lat, species, presence, fraction 字段
all_species_list <- map(csv_files, function(f) {
    dt <- fread(f, select = c("lon", "lat", "species", "presence", "fraction"))
    # 只保留有物种存在的像元，用于作图
    dt <- dt[presence == 1]
    return(dt)
})
all_species_dt <- rbindlist(all_species_list)

# 将物种名称中的 "_" 替换为空格，更适合展示
all_species_dt[, species := gsub("_", " ", species)]

# 按物种包含像元数量排序（确保画出的面板有序）
sp_counts <- all_species_dt[, .N, by = species][order(-N)]
all_species_dt$species <- factor(all_species_dt$species, levels = sp_counts$species)

cat(sprintf("Succesfully loaded %d occurrence points across %d species.\n", 
            nrow(all_species_dt), length(csv_files)))

# 4. 绘制所有 32 物种的多面板图 (Faceted Map)
cat("Plotting the faceted distribution map...\n")

p_all <- ggplot() +
    # 底图
    geom_sf(data = china_sf, fill = "#f4f6f7", color = "#bdc3c7", linewidth = 0.2) +
    # 分布点 (以 fraction 值映射颜色)
    geom_point(data = all_species_dt, 
               aes(x = lon, y = lat, color = fraction), 
               size = 0.3, alpha = 0.8) +
    scale_color_viridis_c(option = "turbo", name = "Occurrence\nFraction", limits = c(0, 1)) +
    facet_wrap(~species, ncol = 6) +
    theme_void(base_size = 11) +
    theme(
        strip.text = element_text(face = "italic", size = 9, margin = margin(b=3)),
        legend.position = "bottom",
        legend.key.width = unit(2, "cm"),
        plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b=10, t=10)),
        plot.subtitle = element_text(size = 11, hjust=0.5, color="grey30", margin=margin(b=15)),
        plot.background = element_rect(fill = "white", color = NA),
        panel.spacing = unit(0.5, "lines")
    ) +
    labs(
        title = "Geographic Distribution of 32 Mammal Species in China",
        subtitle = "Data source: EcoISEA3H (Resolution 9) • CAST Screened Data",
        x = NULL, y = NULL
    )

# 5. 保存
out_file <- file.path(fig_dir, "CAST_32Species_Distribution_Map.png")
ggsave(out_file, p_all, width = 14, height = 16, dpi = 400, bg = "white")

cat(sprintf("✓ Plot saved to: %s\n", out_file))
