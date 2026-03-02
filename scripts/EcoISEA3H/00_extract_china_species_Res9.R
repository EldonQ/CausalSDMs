# 00_extract_china_species_Res9.R
# ============================================================
# Eco-ISEA3H | Resolution 9 | China Region Data Extraction
#
# Steps:
#   1. Load China boundary & find HIDs within China (cached)
#   2. Scan all 30 IUCNRL species files for China occurrences
#   3. Map IUCN SID -> species names via iucnredlist API v4
#      - Stage A: assessments_by_taxonomy() per family
#      - Stage B: assessment_data_many() for full data
#      - Extract: $taxon$scientific_name, $red_list_category$code
#   3b. Attach lon/lat centroids for CAST modeling
#   4. Save annotated species table, summaries & figures
#
# Requires:
#   install.packages(c("sf", "dplyr", "data.table", "ggplot2",
#                      "stringr", "purrr"))
#   devtools::install_github("IUCN-UK/iucnredlist")
# ============================================================

library(sf)
library(dplyr)
library(data.table)
library(ggplot2)
library(stringr)
library(purrr)

# ── Output directories ──────────────────────────────────────
out_dir <- "E:/CausalSDMs/outputs/EcoISEA3H/Res9"
fig_dir <- "E:/CausalSDMs/figures/EcoISEA3H/Res9"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# STEP 1: China HIDs (cached after first run)
# ============================================================
hid_cache <- file.path(out_dir, "china_hids_Res9.rds")
centroids_file <- "E:/CausalSDMs/Eco-ISEA3H/data/Spatial/Text/Centroids_ISEA3H09_Geodetic_V_WGS84.txt"

# Always load centroids (needed for lon/lat later)
cat("Step 1: Loading ISEA3H09 centroids...\n")
centroids_df <- fread(centroids_file)

if (file.exists(hid_cache)) {
    cat("  Loading cached China HIDs...\n")
    china_hids <- readRDS(hid_cache)
} else {
    cat("  Loading China boundary...\n")
    sf_use_s2(FALSE)
    china_sf <- st_read("E:/CausalSDMs/data-main/vector/china.shp", quiet = TRUE)
    if (is.na(st_crs(china_sf)) || st_crs(china_sf)$epsg != 4326) {
        china_sf <- st_transform(china_sf, 4326)
    }
    china_sf <- st_make_valid(china_sf)

    cat("  Intersecting centroids with China boundary...\n")
    centroids_sf <- st_as_sf(centroids_df, coords = c("X", "Y"), crs = 4326)
    idx <- lengths(st_intersects(centroids_sf, china_sf)) > 0
    china_hids <- centroids_sf$HID[idx]
    china_hids_sf <- centroids_sf[idx, ]

    saveRDS(china_hids, hid_cache)
    cat(sprintf("  -> Found %d HIDs in China. Cached.\n", length(china_hids)))

    # Verification map
    p <- ggplot() +
        geom_sf(data = china_sf, fill = "#dce9f5", color = "grey40", linewidth = 0.3) +
        geom_sf(data = china_hids_sf, color = "#e74c3c", size = 0.15, alpha = 0.4) +
        theme_minimal(base_size = 12) +
        labs(
            title    = "Eco-ISEA3H Resolution 9: HIDs in China",
            subtitle = sprintf("n = %d hexagonal grid cells | cell area ~ 2,591 km2", length(china_hids)),
            caption  = "Source: Eco-ISEA3H; Boundary: data-main/vector/china.shp"
        )
    ggsave(file.path(fig_dir, "01_China_HIDs_Verification.png"), p, width = 9, height = 7, dpi = 300)
    cat("  -> Verification map saved.\n")
}

cat(sprintf("  China HIDs: %d\n", length(china_hids)))

# ============================================================
# STEP 2: Scan IUCNRL files for species present in China
# ============================================================
cat("\nStep 2: Scanning IUCNRL species files...\n")

iucn_dir <- "E:/CausalSDMs/Eco-ISEA3H/data/ISEA3H09/IUCNRL_V201901"
species_files <- list.files(iucn_dir, pattern = "\\.txt$", full.names = TRUE)

get_family <- function(f) {
    str_match(basename(f), "V201901_([A-Za-z]+)_Fractions")[, 2]
}

results <- list()

for (file in species_files) {
    family <- get_family(file)
    cat(sprintf("  Processing %-20s ...", family))
    sp_df <- fread(file)
    sp_china <- sp_df[HID %in% china_hids]

    if (nrow(sp_china) > 0) {
        sp_cols <- setdiff(names(sp_china), "HID")
        if (length(sp_cols) > 0) {
            presence_n <- sapply(sp_china[, ..sp_cols], function(v) sum(v > 0, na.rm = TRUE))
            mean_frac <- sapply(sp_china[, ..sp_cols], function(v) mean(v[v > 0], na.rm = TRUE))
            present_mask <- presence_n > 0
            if (any(present_mask)) {
                results[[family]] <- data.frame(
                    Family = family,
                    SID_col = names(presence_n)[present_mask],
                    n_cells = presence_n[present_mask],
                    mean_frac = round(mean_frac[present_mask], 4),
                    stringsAsFactors = FALSE
                )
                cat(sprintf(" %d species found\n", sum(present_mask)))
            } else {
                cat(" 0 species found\n")
            }
        }
    } else {
        cat(" 0 HIDs matched\n")
    }
}

china_sp_raw <- bind_rows(results)
china_sp_raw <- china_sp_raw %>%
    mutate(SID = as.integer(str_remove(str_remove(SID_col, "^SID0*"), "_Fraction$")))

cat(sprintf(
    "\n  -> Total species in China: %d across %d families\n",
    nrow(china_sp_raw), n_distinct(china_sp_raw$Family)
))

# ============================================================
# STEP 3: Map SID -> Species name via iucnredlist (IUCN API v4)
# ============================================================
# Full workflow (verified):
#   Stage A: assessments_by_taxonomy(api, level="family", name=X)
#            -> minimal tibble with sis_taxon_id, assessment_id
#   Stage B: assessment_data_many(api, assessment_ids)
#            -> list of raw API responses
#   Extract: from each response element:
#            $taxon$scientific_name  (NOT $taxonomy!)
#            $sis_taxon_id           (integer)
#            $red_list_category$code (tibble, take [1])
# ============================================================
cat("\nStep 3: Looking up species names from IUCN Red List API v4...\n")

# API token
Sys.setenv(IUCN_REDLIST_KEY = "AnM3rT3DyqYDEX9FDvWhVXUCJfmTQrDEJF4K")
iucn_key <- Sys.getenv("IUCN_REDLIST_KEY")

# Check if cached lookup exists (to avoid re-downloading 500s of data)
lookup_cache <- file.path(out_dir, "iucn_lookup_cache.rds")

if (file.exists(lookup_cache)) {
    cat("  Loading cached IUCN lookup table...\n")
    lookup_df <- readRDS(lookup_cache)
    cat(sprintf("  -> %d species in cache.\n", nrow(lookup_df)))
} else if (nchar(iucn_key) == 0 || !requireNamespace("iucnredlist", quietly = TRUE)) {
    cat("  !! iucnredlist not installed or API key missing.\n")
    cat("  !! Run: devtools::install_github('IUCN-UK/iucnredlist')\n")
    lookup_df <- NULL
} else {
    library(iucnredlist)

    cat("  Initialising API connection...\n")
    api <- init_api(iucn_key)

    # Stage A: get assessment IDs per family
    china_families <- unique(china_sp_raw$Family)
    cat(sprintf("  Stage A: Querying %d families...\n", length(china_families)))

    family_assessments <- list()
    for (fam in china_families) {
        cat(sprintf("    -> %s\n", fam))
        tryCatch(
            {
                res <- assessments_by_taxonomy(api,
                    level = "family", name = fam,
                    latest = TRUE, wait_time = 0.5
                )
                if (!is.null(res) && nrow(res) > 0) {
                    res$Family <- fam
                    family_assessments[[fam]] <- res
                }
            },
            error = function(e) {
                cat(sprintf("    [!] %s failed: %s\n", fam, conditionMessage(e)))
            }
        )
    }

    all_assessments <- bind_rows(family_assessments)
    cat(sprintf("  Stage A done: %d assessments.\n", nrow(all_assessments)))

    # Stage B: fetch full assessment data
    cat("  Stage B: Fetching full assessment data (this may take ~8 min)...\n")
    full_data <- assessment_data_many(api, assessment_ids = all_assessments$assessment_id)

    # Extract taxonomy using simple for loop (verified working)
    cat("  Extracting taxonomy from full_data...\n")
    lookup_list <- list()
    for (i in seq_along(full_data)) {
        tryCatch(
            {
                asmt <- full_data[[i]]
                tx <- asmt$taxon
                rl_cat <- asmt$red_list_category
                lookup_list[[i]] <- data.frame(
                    sid = as.integer(asmt$sis_taxon_id),
                    scientific_name = as.character(tx$scientific_name[1]),
                    order_name = as.character(tx$order_name[1]),
                    class_name = as.character(tx$class_name[1]),
                    family_api = as.character(tx$family_name[1]),
                    category = as.character(rl_cat$code[1]),
                    stringsAsFactors = FALSE
                )
            },
            error = function(e) {}
        )
    }

    lookup_df <- bind_rows(lookup_list)
    lookup_df <- lookup_df[!is.na(lookup_df$sid), ]
    lookup_df <- lookup_df[!duplicated(lookup_df$sid), ]

    # Cache the lookup table so we never re-download
    saveRDS(lookup_df, lookup_cache)
    cat(sprintf("  -> Extracted & cached %d unique species.\n", nrow(lookup_df)))
}

# Join to China species table
if (!is.null(lookup_df)) {
    china_sp_annotated <- china_sp_raw %>%
        left_join(lookup_df, by = c("SID" = "sid")) %>%
        select(
            Family, SID, SID_col, scientific_name, category,
            order_name, class_name, n_cells, mean_frac
        ) %>%
        arrange(Family, desc(n_cells))
} else {
    china_sp_annotated <- china_sp_raw %>%
        mutate(
            scientific_name = NA_character_, category = NA_character_,
            order_name = NA_character_, class_name = NA_character_
        ) %>%
        select(
            Family, SID, SID_col, scientific_name, category,
            order_name, class_name, n_cells, mean_frac
        ) %>%
        arrange(Family, desc(n_cells))
}

n_named <- sum(!is.na(china_sp_annotated$scientific_name))
cat(sprintf(
    "\n  -> Annotated %d / %d species with scientific names.\n",
    n_named, nrow(china_sp_annotated)
))

# ============================================================
# STEP 3b: Attach lon/lat from centroids (for CAST modeling)
# ============================================================
cat("\nStep 3b: Preparing lon/lat coordinates for CAST...\n")

china_centroids <- copy(centroids_df[centroids_df$HID %in% china_hids, ])
setnames(china_centroids, c("X", "Y"), c("lon", "lat"), skip_absent = TRUE)

cat(sprintf("  -> %d China HIDs with lon/lat ready.\n", nrow(china_centroids)))

# ============================================================
# STEP 4: Save all outputs & figures
# ============================================================
cat("\nStep 4: Saving outputs...\n")

# Annotated species list
write.csv(china_sp_annotated,
    file.path(out_dir, "China_Species_List_Res9_Annotated.csv"),
    row.names = FALSE
)
cat("  -> China_Species_List_Res9_Annotated.csv\n")

# China centroids
fwrite(china_centroids, file.path(out_dir, "China_Centroids_Res9.csv"))
cat("  -> China_Centroids_Res9.csv\n")

# Family summary
family_summary <- china_sp_annotated %>%
    group_by(Family) %>%
    summarise(
        n_species    = n(),
        max_coverage = max(n_cells),
        avg_coverage = round(mean(n_cells), 1)
    ) %>%
    arrange(desc(n_species))

write.csv(family_summary,
    file.path(out_dir, "China_Species_Family_Summary_Res9.csv"),
    row.names = FALSE
)
cat("  -> China_Species_Family_Summary_Res9.csv\n")

# ── Figure: species count per family ─────────────────────
p2 <- ggplot(family_summary, aes(x = reorder(Family, n_species), y = n_species)) +
    geom_col(fill = "#2980b9", width = 0.7) +
    geom_text(aes(label = n_species), hjust = -0.2, size = 3.5) +
    coord_flip() +
    theme_minimal(base_size = 12) +
    labs(
        title = "Species Count per Family in China (IUCNRL Res 9)",
        subtitle = paste0("Total: ", nrow(china_sp_annotated), " species"),
        x = NULL, y = "Number of species"
    )
ggsave(file.path(fig_dir, "02_China_Species_per_Family.png"), p2,
    width = 9, height = 6, dpi = 300
)
cat("  -> 02_China_Species_per_Family.png\n")

# ── Figure: top-40 species by coverage ───────────────────
top40 <- china_sp_annotated %>%
    slice_max(n_cells, n = 40) %>%
    mutate(label = ifelse(!is.na(scientific_name), scientific_name, SID_col))

p3 <- ggplot(top40, aes(x = reorder(label, n_cells), y = n_cells, fill = Family)) +
    geom_col(width = 0.75) +
    coord_flip() +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom") +
    labs(
        title = "Top 40 Species by Grid-Cell Coverage in China (Res 9)",
        subtitle = "n_cells = number of HIDs where species fraction > 0",
        x = NULL, y = "Number of HIDs with presence"
    )
ggsave(file.path(fig_dir, "03_China_Top40_Species_Coverage.png"), p3,
    width = 10, height = 12, dpi = 300
)
cat("  -> 03_China_Top40_Species_Coverage.png\n")

# ── Figure: Top-10 species spatial distribution in China ──
cat("\nStep 5: Generating Top-10 species distribution maps...\n")

# Select top 10 species by coverage
top10 <- china_sp_annotated %>%
    slice_max(n_cells, n = 10) %>%
    mutate(label = ifelse(!is.na(scientific_name), scientific_name, SID_col))

# Load China boundary for background (from cache or file)
sf_use_s2(FALSE)
china_shp_file <- "E:/CausalSDMs/data-main/vector/china.shp"
china_sf <- st_read(china_shp_file, quiet = TRUE)
if (is.na(st_crs(china_sf)) || st_crs(china_sf)$epsg != 4326) {
    china_sf <- st_transform(china_sf, 4326)
}
china_sf <- st_make_valid(china_sf)

# For each top-10 species, extract fraction data from original files
# and join with centroids for lon/lat
cat("  Reading fraction data for top-10 species...\n")

sp_map_data <- list()
for (i in seq_len(nrow(top10))) {
    sp_family <- top10$Family[i]
    sp_col <- top10$SID_col[i]
    sp_label <- top10$label[i]

    # Find the species file
    sp_file <- list.files(iucn_dir, pattern = sp_family, full.names = TRUE)
    if (length(sp_file) == 0) next

    dt <- fread(sp_file[1], select = c("HID", sp_col))
    # Filter to China HIDs with presence
    dt <- dt[HID %in% china_hids & get(sp_col) > 0]
    setnames(dt, sp_col, "fraction")

    # Join lon/lat from china_centroids
    dt <- merge(dt, china_centroids, by = "HID", all.x = TRUE)
    dt$species <- sp_label

    sp_map_data[[i]] <- dt
    cat(sprintf("    %2d. %-30s %d points\n", i, sp_label, nrow(dt)))
}

sp_map_all <- rbindlist(sp_map_data)
# Order species factor by coverage (descending)
sp_map_all$species <- factor(sp_map_all$species,
    levels = top10$label
)

# Create faceted map
p4 <- ggplot() +
    geom_sf(
        data = china_sf, fill = "#f5f5f5", color = "grey60",
        linewidth = 0.2
    ) +
    geom_point(
        data = sp_map_all,
        aes(x = lon, y = lat, color = fraction),
        size = 0.4, alpha = 0.7
    ) +
    scale_color_viridis_c(
        option = "turbo", name = "Fraction",
        limits = c(0, 1)
    ) +
    facet_wrap(~species, ncol = 2) +
    theme_minimal(base_size = 9) +
    theme(
        strip.text       = element_text(face = "italic", size = 9),
        legend.position  = "bottom",
        legend.key.width = unit(1.5, "cm"),
        panel.grid       = element_blank(),
        axis.text        = element_text(size = 6)
    ) +
    labs(
        title = "Top 10 Species Distribution in China (Eco-ISEA3H Res 9)",
        subtitle = "Color = species occurrence fraction within each hexagonal grid cell",
        x = "Longitude", y = "Latitude"
    )

ggsave(file.path(fig_dir, "04_China_Top10_Species_Distribution.png"), p4,
    width = 10, height = 22, dpi = 300
)
cat("  -> 04_China_Top10_Species_Distribution.png\n")

# ── Final summary ────────────────────────────────────────
cat("\n================================================\n")
cat("=== Eco-ISEA3H China Res9 Extraction Complete ===\n")
cat("================================================\n")
cat(sprintf("  HIDs (grids):    %d  (lon/lat attached)\n", length(china_hids)))
cat(sprintf("  Species total:   %d\n", nrow(china_sp_annotated)))
cat(sprintf("  Species named:   %d\n", n_named))
cat(sprintf("  Families:        %d\n", nrow(family_summary)))
cat(sprintf("  Outputs dir:     %s\n", out_dir))
cat(sprintf("  Figures dir:     %s\n", fig_dir))
cat("================================================\n")
