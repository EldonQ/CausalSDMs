# 01_prepare_CAST_envdata.R
# ============================================================
# Eco-ISEA3H | Resolution 9 | China Region
# Prepare environmental data for CAST modeling pipeline
#
# This script:
#   1. Loads annotated species list from 00_extract (with filtering)
#   2. Loads ALL environmental data layers from Eco-ISEA3H
#   3. Subsets to China HIDs and merges into a single env table
#   4. For each qualifying species, creates a model-ready dataset:
#      HID, lon, lat, presence (0/1), fraction, env variables
#   5. Saves per-species CSVs and a summary
#
# ── Configurable parameters ──────────────────────────────
MIN_CELLS <- 100 # Minimum grid-cell coverage to include a species
FRAC_THRES <- 0.0 # Fraction threshold to convert to presence (> this)
# ─────────────────────────────────────────────────────────
# ============================================================

library(data.table)
library(dplyr)
library(stringr)
library(ggplot2)

# ── Paths ────────────────────────────────────────────────
base_dir <- "E:/CausalSDMs/Eco-ISEA3H/data/ISEA3H09"
out_dir_00 <- "E:/CausalSDMs/outputs/EcoISEA3H/Res9"
out_dir <- "E:/CausalSDMs/outputs/EcoISEA3H/Res9/CAST_ready"
fig_dir <- "E:/CausalSDMs/figures/EcoISEA3H/Res9"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("============================================================\n")
cat("  01_prepare_CAST_envdata.R\n")
cat(sprintf("  MIN_CELLS = %d  |  FRAC_THRES = %.2f\n", MIN_CELLS, FRAC_THRES))
cat("============================================================\n")

# ============================================================
# STEP 1: Load species list & filter
# ============================================================
cat("\nStep 1: Loading annotated species list...\n")

sp_file <- file.path(out_dir_00, "China_Species_List_Res9_Annotated.csv")
sp_all <- fread(sp_file)
cat(sprintf("  Total species: %d\n", nrow(sp_all)))

# Filter species by minimum coverage
sp_selected <- sp_all[n_cells >= MIN_CELLS]
cat(sprintf(
    "  After filter (n_cells >= %d): %d species\n",
    MIN_CELLS, nrow(sp_selected)
))
print(sp_selected[, .(Family, SID, scientific_name, category, n_cells, mean_frac)])

# Load centroids
coord_file <- file.path(out_dir_00, "China_Centroids_Res9.csv")
china_coords <- fread(coord_file)
china_hids <- china_coords$HID
cat(sprintf("  China HIDs: %d\n", length(china_hids)))

# ============================================================
# STEP 2: Load all environmental variables
# ============================================================
cat("\nStep 2: Loading environmental variables...\n")

# ── 2a. WorldClim BIO variables (19 bioclimatic) ─────────
cat("  [1/7] WorldClim BIO1-19...\n")
bio_file <- file.path(
    base_dir, "WorldClim30AS_V02",
    "ISEA3H09_WorldClim30AS_V02_BIO_Mean.txt"
)
bio_df <- fread(bio_file)
# Rename: BIO01_Mean -> bio01, etc. (cleaner for modeling)
bio_cols <- grep("^BIO", names(bio_df), value = TRUE)
new_bio <- tolower(str_replace(bio_cols, "_Mean$", ""))
setnames(bio_df, bio_cols, new_bio)
bio_df <- bio_df[HID %in% china_hids]
cat(sprintf("    %d variables, %d rows\n", length(new_bio), nrow(bio_df)))

# ── 2b. Elevation (SRTM) ────────────────────────────────
cat("  [2/7] Elevation (SRTM)...\n")
elev_file <- file.path(
    base_dir, "SRTM30PLUS_V11",
    "ISEA3H09_SRTM30PLUS_V11_Elevation_Mean.txt"
)
elev_df <- fread(elev_file)
setnames(elev_df, "Elevation_Mean", "elevation")
elev_df <- elev_df[HID %in% china_hids]
cat(sprintf("    %d rows\n", nrow(elev_df)))

# ── 2c. ENVIREM Climate indices ──────────────────────────
cat("  [3/7] ENVIREM Climate...\n")
envclim_file <- file.path(
    base_dir, "ENVIREM30AS_V01",
    "ISEA3H09_ENVIREM30AS_V01_Climate_Mean.txt"
)
envclim_df <- fread(envclim_file)
# Clean names: AnnualPET_Mean -> annualpet, etc.
ec_cols <- setdiff(names(envclim_df), "HID")
new_ec <- tolower(str_replace(ec_cols, "_Mean$", ""))
setnames(envclim_df, ec_cols, new_ec)
envclim_df <- envclim_df[HID %in% china_hids]
cat(sprintf("    %d variables, %d rows\n", length(new_ec), nrow(envclim_df)))

# ── 2d. ENVIREM Topography ───────────────────────────────
cat("  [4/7] ENVIREM Topography...\n")
topo_file <- file.path(
    base_dir, "ENVIREM30AS_V01",
    "ISEA3H09_ENVIREM30AS_V01_Topography_Mean.txt"
)
topo_df <- fread(topo_file)
tc_cols <- setdiff(names(topo_df), "HID")
new_tc <- tolower(str_replace(tc_cols, "_Mean$", ""))
setnames(topo_df, tc_cols, new_tc)
topo_df <- topo_df[HID %in% china_hids]
cat(sprintf("    %d variables, %d rows\n", length(new_tc), nrow(topo_df)))

# ── 2e. MODIS Vegetation Cover (Tree, NonTree, Bare) ─────
cat("  [5/7] MODIS VCF...\n")
vcf_file <- file.path(
    base_dir, "MOD44B_V06",
    "ISEA3H09_MOD44B_V06_Y2018_VCF_Mean.txt"
)
vcf_df <- fread(vcf_file)
vc_cols <- setdiff(names(vcf_df), "HID")
new_vc <- tolower(str_replace(vc_cols, "_Mean$", ""))
setnames(vcf_df, vc_cols, new_vc)
vcf_df <- vcf_df[HID %in% china_hids]
cat(sprintf("    %d variables, %d rows\n", length(new_vc), nrow(vcf_df)))

# ── 2f. ECMWF ERA-40 Climate Extremes (ETCCDI) ──────────
cat("  [6/7] ECMWF ETCCDI climate extremes...\n")
etccdi_file <- file.path(
    base_dir, "ECMWF_ERA40",
    "ISEA3H09_ECMWF_ERA40_Y1958_Y2001_ETCCDI_IDW1N10.txt"
)
etccdi_df <- fread(etccdi_file)
et_cols <- setdiff(names(etccdi_df), "HID")
new_et <- paste0("etccdi_", tolower(str_replace(et_cols, "_IDW1N10$", "")))
setnames(etccdi_df, et_cols, new_et)
etccdi_df <- etccdi_df[HID %in% china_hids]
cat(sprintf("    %d variables, %d rows\n", length(new_et), nrow(etccdi_df)))

# ── 2g. MODIS Land Cover IGBP Fractions (16 continuous) ──
# IGBP class codes:
#   01=Evergreen Needleleaf Forest   02=Evergreen Broadleaf Forest
#   03=Deciduous Needleleaf Forest   04=Deciduous Broadleaf Forest
#   05=Mixed Forest                  06=Closed Shrublands
#   07=Open Shrublands               08=Woody Savannas
#   09=Savannas                      10=Grasslands
#   11=Permanent Wetlands            12=Croplands
#   13=Urban and Built-Up            14=Cropland/Natural Vegetation Mosaic
#   15=Snow and Ice                  16=Barren or Sparsely Vegetated
cat("  [7/7] MODIS Land Cover (IGBP Fractions)...\n")
lc_file <- file.path(
    base_dir, "MCD12Q1_V06",
    "ISEA3H09_MCD12Q1_V06_Y2018_IGBP_Fractions.txt"
)
# Check if 2018 exists; fallback to latest available
if (!file.exists(lc_file)) {
    lc_files <- list.files(file.path(base_dir, "MCD12Q1_V06"),
        pattern = "Fractions\\.txt$", full.names = TRUE
    )
    lc_file <- tail(sort(lc_files), 1)
}
lc_df <- fread(lc_file)
lc_cols <- grep("^IGBP_", names(lc_df), value = TRUE)
# Rename: IGBP_01_Fraction -> igbp_01, etc.
new_lc <- tolower(str_replace(lc_cols, "^IGBP_(\\d+)_Fraction$", "igbp_\\1"))
setnames(lc_df, lc_cols, new_lc)
lc_df <- lc_df[HID %in% china_hids]
cat(sprintf("    %d variables, %d rows\n", length(new_lc), nrow(lc_df)))

# ============================================================
# STEP 3: Merge all env data by HID
# ============================================================
cat("\nStep 3: Merging all environmental layers...\n")

env_master <- china_coords # starts with HID, lon, lat
env_master <- merge(env_master, bio_df, by = "HID", all.x = TRUE)
env_master <- merge(env_master, elev_df, by = "HID", all.x = TRUE)
env_master <- merge(env_master, envclim_df, by = "HID", all.x = TRUE)
env_master <- merge(env_master, topo_df, by = "HID", all.x = TRUE)
env_master <- merge(env_master, vcf_df, by = "HID", all.x = TRUE)
env_master <- merge(env_master, etccdi_df, by = "HID", all.x = TRUE)
env_master <- merge(env_master, lc_df, by = "HID", all.x = TRUE)

env_cols <- setdiff(names(env_master), c("HID", "lon", "lat"))
cat(sprintf(
    "  -> %d HIDs x %d environmental variables\n",
    nrow(env_master), length(env_cols)
))

# Replace NoData sentinels with NA
# WorldClim: -100, ENVIREM: -1000, VCF: -1
nodata_vals <- c(-100, -1000, -1)
for (col in env_cols) {
    if (is.numeric(env_master[[col]])) {
        env_master[get(col) %in% nodata_vals, (col) := NA]
    }
}

# Report NA coverage
na_pct <- sapply(env_master[, ..env_cols], function(x) round(100 * mean(is.na(x)), 1))
cat("  Variables with >10% NA:\n")
high_na <- na_pct[na_pct > 10]
if (length(high_na) > 0) {
    for (nm in names(high_na)) cat(sprintf("    %-35s %5.1f%%\n", nm, high_na[nm]))
} else {
    cat("    None — good coverage!\n")
}

# Save master env table
fwrite(env_master, file.path(out_dir, "China_EnvData_Res9_Master.csv"))
cat(sprintf(
    "  -> China_EnvData_Res9_Master.csv (%d x %d)\n",
    nrow(env_master), ncol(env_master)
))

# ============================================================
# STEP 4: Create per-species CAST-ready datasets
# ============================================================
cat("\nStep 4: Building per-species CAST-ready datasets...\n")
cat(sprintf(
    "  Generating for %d species (n_cells >= %d)...\n",
    nrow(sp_selected), MIN_CELLS
))

# Species data directory
sp_data_dir <- file.path(out_dir, "species_data")
dir.create(sp_data_dir, recursive = TRUE, showWarnings = FALSE)

iucn_dir <- file.path(base_dir, "IUCNRL_V201901")

sp_summary <- list()

for (i in seq_len(nrow(sp_selected))) {
    sp <- sp_selected[i]
    sp_name <- ifelse(!is.na(sp$scientific_name), sp$scientific_name,
        paste0("SID_", sp$SID)
    )
    sp_label <- gsub(" ", "_", sp_name)

    cat(sprintf(
        "  [%2d/%d] %-30s (SID %d, %d cells)...\n",
        i, nrow(sp_selected), sp_name, sp$SID, sp$n_cells
    ))

    # Read fraction data from family file
    fam_file <- list.files(iucn_dir, pattern = sp$Family, full.names = TRUE)
    if (length(fam_file) == 0) {
        cat("    [!] Family file not found, skipping.\n")
        next
    }

    sp_frac <- fread(fam_file[1], select = c("HID", sp$SID_col))
    setnames(sp_frac, sp$SID_col, "fraction")

    # Create full dataset: all China HIDs with env data
    sp_df <- copy(env_master)
    sp_df <- merge(sp_df, sp_frac, by = "HID", all.x = TRUE)
    sp_df[is.na(fraction), fraction := 0]
    sp_df[, presence := as.integer(fraction > FRAC_THRES)]

    # Add species metadata
    sp_df[, `:=`(
        species  = sp_name,
        sid      = sp$SID,
        family   = sp$Family,
        category = sp$category
    )]

    # Reorder columns
    meta_cols <- c(
        "HID", "lon", "lat", "species", "sid", "family", "category",
        "presence", "fraction"
    )
    setcolorder(sp_df, c(meta_cols, setdiff(names(sp_df), meta_cols)))

    # Remove rows with all NA in env variables (ocean/no-data cells)
    n_before <- nrow(sp_df)
    env_numeric <- intersect(env_cols, names(sp_df))
    sp_df <- sp_df[rowSums(!is.na(sp_df[, ..env_numeric])) > 0]
    n_after <- nrow(sp_df)

    # Save
    out_file <- file.path(
        sp_data_dir,
        sprintf("CAST_%s_Res9.csv", sp_label)
    )
    fwrite(sp_df, out_file)

    # Record summary
    sp_summary[[i]] <- data.frame(
        species = sp_name,
        SID = sp$SID,
        family = sp$Family,
        category = sp$category,
        n_hids = n_after,
        n_presence = sum(sp_df$presence),
        n_absence = sum(sp_df$presence == 0),
        prevalence = round(mean(sp_df$presence), 4),
        n_env_vars = length(env_numeric),
        n_env_na_rows = n_before - n_after,
        file = basename(out_file),
        stringsAsFactors = FALSE
    )
}

# ============================================================
# STEP 5: Save summary and generate overview figure
# ============================================================
cat("\nStep 5: Saving summary...\n")

sp_summary_df <- bind_rows(sp_summary)
fwrite(sp_summary_df, file.path(out_dir, "CAST_Species_Summary.csv"))

cat("\n  === CAST-Ready Species Summary ===\n")
cat(sprintf(
    "  %-30s %-6s %8s %8s %10s\n",
    "Species", "IUCN", "Present", "Absent", "Prevalence"
))
cat(paste(rep("-", 75), collapse = ""), "\n")
tmp <- sp_summary_df[order(-sp_summary_df$n_presence), ]
for (j in seq_len(nrow(tmp))) {
    cat(sprintf(
        "  %-30s %-6s %8d %8d %9.1f%%\n",
        tmp$species[j], tmp$category[j],
        tmp$n_presence[j], tmp$n_absence[j],
        tmp$prevalence[j] * 100
    ))
}

# Prevalence overview figure
p_prev <- ggplot(
    sp_summary_df,
    aes(
        x = reorder(species, prevalence), y = prevalence,
        fill = category
    )
) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("%.1f%%", prevalence * 100)),
        hjust = -0.1, size = 3
    ) +
    coord_flip() +
    scale_y_continuous(labels = scales::percent, limits = c(0, 1)) +
    theme_minimal(base_size = 11) +
    labs(
        title = "Species Prevalence in China (CAST-Ready Dataset)",
        subtitle = sprintf(
            "n_cells >= %d | %d species | %d env variables",
            MIN_CELLS, nrow(sp_summary_df), sp_summary_df$n_env_vars[1]
        ),
        x = NULL, y = "Prevalence (proportion of cells with presence)",
        fill = "IUCN Status"
    )
ggsave(file.path(fig_dir, "05_CAST_Species_Prevalence.png"), p_prev,
    width = 10, height = 8, dpi = 300
)
cat("  -> 05_CAST_Species_Prevalence.png\n")

# ── Env variable correlation heatmap ─────────────────────
cat("  Generating env correlation heatmap...\n")

# Select only numeric env columns for correlation
env_numeric_data <- env_master[, ..env_cols]
env_numeric_data <- env_numeric_data[, .SD, .SDcols = sapply(env_numeric_data, is.numeric)]
# Remove columns with >30% NA
good_cols <- names(which(colMeans(is.na(env_numeric_data)) < 0.3))
cor_mat <- cor(env_numeric_data[, ..good_cols], use = "pairwise.complete.obs")

# Simple heatmap via ggplot (top-triangle pairs with high correlation)
high_cor <- which(abs(cor_mat) > 0.85 & upper.tri(cor_mat), arr.ind = TRUE)
if (nrow(high_cor) > 0) {
    high_cor_df <- data.frame(
        var1 = good_cols[high_cor[, 1]],
        var2 = good_cols[high_cor[, 2]],
        r = cor_mat[high_cor]
    ) %>% arrange(desc(abs(r)))

    fwrite(high_cor_df, file.path(out_dir, "Env_HighCorrelations_r85.csv"))
    cat(sprintf("  -> %d variable pairs with |r| > 0.85 saved.\n", nrow(high_cor_df)))
} else {
    cat("  No variable pairs with |r| > 0.85.\n")
}

# ── Final summary ────────────────────────────────────────
cat("\n====================================================\n")
cat("=== CAST Environment Data Preparation Complete! ===\n")
cat("====================================================\n")
cat(sprintf("  Species qualifying:  %d (n_cells >= %d)\n", nrow(sp_summary_df), MIN_CELLS))
cat(sprintf("  Env variables:       %d\n", length(env_cols)))
cat(sprintf("  Master env table:    %s\n", file.path(out_dir, "China_EnvData_Res9_Master.csv")))
cat(sprintf("  Per-species dir:     %s\n", sp_data_dir))
cat(sprintf("  Summary:             %s\n", file.path(out_dir, "CAST_Species_Summary.csv")))
cat("====================================================\n")
cat("\nReady for CAST pipeline! Each species CSV contains:\n")
cat("  HID, lon, lat, species, sid, family, category,\n")
cat("  presence (0/1), fraction, + all env variables\n")
cat("====================================================\n")
