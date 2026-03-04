# 01_prepare_CAST_envdata.R
# ============================================================
# Plant Dataset - Environment Preparation
# Reference: E:\CausalSDMs\scripts\EcoISEA3H\01_prepare_CAST_envdata.R
# ============================================================

library(data.table)
library(dplyr)
library(reticulate)

out_dir_base <- "E:/CausalSDMs/outputs/Plant"
out_dir <- file.path(out_dir_base, "CAST_ready")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

data_dir <- "E:/CausalSDMs/maskSDMData/data"

cat("============================================================\n")
cat("  01_prepare_CAST_envdata.R (Plant Dataset)\n")
cat("============================================================\n")

# 1. Load data
cat("Step 1: Loading environmental datasets...\n")
loc <- fread(file.path(data_dir, "location_data.csv"))
wc <- fread(file.path(data_dir, "worldclim_data.csv"))
sg <- fread(file.path(data_dir, "soilgrid_data.csv"))
topo <- fread(file.path(data_dir, "topographic_data.csv"))
hum <- fread(file.path(data_dir, "human_data.csv"))

clean_df <- function(df) {
    if ("V1" %in% names(df)) df[, `V1` := NULL]
    if ("Unnamed: 0" %in% names(df)) df[, `Unnamed: 0` := NULL]
    if ("Unnamed: 0.1" %in% names(df)) df[, `Unnamed: 0.1` := NULL]
    return(df)
}
loc <- clean_df(loc)
wc <- clean_df(wc)
sg <- clean_df(sg)
topo <- clean_df(topo)
hum <- clean_df(hum)

# 2. Merge Environment Master
cat("Step 2: Merging into Master Environment table...\n")
env_master <- loc
env_master <- merge(env_master, wc, by = "PlotObservationID", all.x = TRUE)
env_master <- merge(env_master, sg, by = "PlotObservationID", all.x = TRUE)
env_master <- merge(env_master, topo, by = "PlotObservationID", all.x = TRUE)
env_master <- merge(env_master, hum, by = "PlotObservationID", all.x = TRUE)

cat("  -> Master env dim:", dim(env_master)[1], "x", dim(env_master)[2], "\n")

# Handling NA values
init_rows <- nrow(env_master)
env_master <- na.omit(env_master)
cat(sprintf("  -> After NA omit: %d plots remaining out of %d\n", nrow(env_master), init_rows))

fwrite(env_master, file.path(out_dir, "Plant_EnvData_Master.csv"))

# 3. Read species occurrences
cat("\nStep 3: Processing Species Occurrences via Numpy...\n")
use_python("C:/Users/LQ/.conda/envs/lab/python.exe", required = FALSE)
np <- import("numpy", convert = FALSE)
occ_file <- file.path(data_dir, "species_occurrences.npy")
# Numpy boolean arrays load nicely through reticulate
occ <- np$load(occ_file, mmap_mode = "r")

sp_names <- fread(file.path(data_dir, "species_names.csv"))

cat("  Calculating species prevalence... (this may take a moment)\n")
col_sums <- py_to_r(occ$sum(axis = 0L))

sp_names$n_presence <- col_sums

# Select top 30 most prevalent species to keep experiments tractable
sp_selected <- sp_names[order(-n_presence)][1:30]
cat(sprintf("  Selecting top %d most prevalent species.\n", nrow(sp_selected)))
print(sp_selected[, c("Species Name", "n_presence")])

sp_dir <- file.path(out_dir, "species_data")
dir.create(sp_dir, showWarnings = FALSE)

sp_summary <- list()

cat("\nStep 4: Writing per-species CAST-ready CSVs...\n")

loc_orig <- fread(file.path(data_dir, "location_data.csv"))
loc_orig$orig_idx <- 1:nrow(loc_orig)
env_master <- merge(env_master, loc_orig[, .(PlotObservationID, orig_idx)], by = "PlotObservationID", all.x = TRUE)

for (i in seq_len(nrow(sp_selected))) {
    nm <- sp_selected$`Species Name`[i]
    idx <- sp_selected$Index[i]

    cat(sprintf(
        "  [%2d/30] %-25s (Index: %d, Presences: %d)...\n",
        i, nm, idx, sp_selected$n_presence[i]
    ))

    # Python index is 0-based
    occ_vec <- py_to_r(occ[, as.integer(idx)])

    # Map occurrence directly using orig_idx which refers to row index in species_occurrences.npy
    sp_df <- copy(env_master)
    sp_df$presence <- as.integer(occ_vec[sp_df$orig_idx])
    sp_df$species <- nm
    sp_df[, orig_idx := NULL]

    setcolorder(sp_df, c("PlotObservationID", "Longitude", "Latitude", "species", "presence"))

    sp_name_clean <- gsub("[^A-Za-z0-9]", "_", nm)
    out_file <- file.path(sp_dir, sprintf("CAST_%s.csv", sp_name_clean))

    fwrite(sp_df, out_file)

    sp_summary[[i]] <- data.frame(
        species = nm,
        index = idx,
        n_presence = sum(sp_df$presence),
        n_absence = sum(sp_df$presence == 0),
        file = basename(out_file)
    )
}

smry_df <- bind_rows(sp_summary)
fwrite(smry_df, file.path(out_dir, "CAST_Species_Summary.csv"))
cat("============================================================\n")
cat(" Data preparation complete! Proceed to VIF screening.\n")
cat("============================================================\n")
