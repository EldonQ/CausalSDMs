################################################################################
# 00_data_preparation.R
# CAST v3 — Data preparation for ALL disdat regions (6 regions, 226 species)
#
# Dataset: Elith et al. (2020) SDM Benchmark Dataset
#   Paper: Elith J. et al. (2020) Presence-only and presence-absence data for
#         comparing species distribution modeling methods. Biodiversity
#         Informatics, 15(2), 69-80.
#   R pkg: disdat (CRAN)
#   Install: install.packages("disdat") → ready to use
#
# 6 Regions:
#   AWT — Australian Wet Tropics   (birds, bats, reptiles)
#   CAN — Ontario, Canada          (birds)
#   NSW — New South Wales           (plants)
#   NZ  — New Zealand               (plants)
#   SA  — South America              (plants)
#   SWI — Switzerland                (plants)
#
# Output per region:
#   output/case2/{REGION}/train_data_{sp}.csv
#   output/case2/{REGION}/test_data_{sp}.csv
#   output/case2/{REGION}/species_summary.csv
#   output/case2/{REGION}/env_variable_stats.csv
# Global:
#   output/case2/all_regions_summary.csv
################################################################################

rm(list = ls())
gc()
setwd("E:/CausalSDMs")

# ---- Dependencies ----
if (!require(disdat, quietly = TRUE)) {
    install.packages("disdat", repos = "https://cloud.r-project.org")
}
library(disdat)
library(tidyverse)

# ---- Configuration ----
# Regions to process (all 6 disdat regions)
# AWT和NSW有多个group，需分组加载PA数据
regions <- c("AWT", "CAN", "NSW", "NZ", "SA", "SWI")
region_groups <- list(
    AWT = c("bird", "plant"),
    CAN = NULL,
    NSW = c("ba", "db", "nb", "ot", "ou", "rt", "ru", "sr"),
    NZ  = NULL,
    SA  = NULL,
    SWI = NULL
)
min_po <- 200

cat("======================================================================\n")
cat("  CAST v3 Data Preparation — disdat 6 Regions\n")
cat("  Elith et al. (2020) SDM Benchmark (226 species)\n")
cat("======================================================================\n\n")

all_regions_summary <- data.frame()

for (region in regions) {
    cat(sprintf("\n═══ Region: %s ═══\n", region))

    dir.create(sprintf("output/case2/%s", region), recursive = TRUE, showWarnings = FALSE)
    dir.create(sprintf("figures/case2/%s", region), recursive = TRUE, showWarnings = FALSE)

    # ---- 1. Load disdat data ----
    po <- disdat::disPo(region)
    bg <- disdat::disBg(region)

    # PA数据：有group的区域需逐group加载后合并
    groups <- region_groups[[region]]
    if (is.null(groups)) {
        pa <- disdat::disPa(region)
        env <- disdat::disEnv(region)
    } else {
        pa_list <- lapply(groups, function(g) {
            tryCatch(disdat::disPa(region, group = g),
                error = function(e) NULL)
        })
        pa_list <- pa_list[!sapply(pa_list, is.null)]
        # 合并所有group的PA（取共有列的并集）
        if (length(pa_list) > 0) {
            all_cols <- unique(unlist(lapply(pa_list, names)))
            pa_list <- lapply(pa_list, function(df) {
                missing <- setdiff(all_cols, names(df))
                for (m in missing) df[[m]] <- NA
                df[, all_cols]
            })
            pa <- do.call(rbind, pa_list)
        } else {
            cat(sprintf("  ⚠ No PA data loaded for %s → skip\n", region))
            next
        }
        env_list <- lapply(groups, function(g) {
            tryCatch(disdat::disEnv(region, group = g),
                error = function(e) NULL)
        })
        env_list <- env_list[!sapply(env_list, is.null)]
        if (length(env_list) > 0) {
            env <- do.call(rbind, env_list)
            env <- env[!duplicated(env), ]
        } else {
            env <- data.frame()
        }
    }

    cat(sprintf(
        "  PO: %d rows | BG: %d rows | PA: %d rows\n",
        nrow(po), nrow(bg), nrow(pa)
    ))

    # ---- 2. Identify env vars and species ----
    meta_cols <- c("siteid", "spid", "x", "y", "occ", "group")
    env_cols <- setdiff(names(po), meta_cols)
    # Also remove env columns not present in bg
    env_cols <- intersect(env_cols, names(bg))
    cat(sprintf(
        "  Env vars (%d): %s\n", length(env_cols),
        paste(env_cols, collapse = ", ")
    ))

    # Species from PO
    species_list <- sort(unique(po$spid))
    cat(sprintf("  Species (PO): %d\n", length(species_list)))

    # PA species columns (pattern varies: swi01, awt01, can01, etc.)
    region_prefix <- tolower(region)
    sp_cols_pa <- grep(paste0("^", region_prefix), names(pa), value = TRUE)
    cat(sprintf("  PA species columns: %d\n", length(sp_cols_pa)))

    # ---- 3. Merge PA with environmental data ----
    merge_keys <- intersect(names(pa), names(env))
    pa_full <- merge(pa, env, by = merge_keys, all.x = TRUE)

    # Verify env columns present
    env_in_pa <- intersect(env_cols, names(pa_full))
    if (length(env_in_pa) < length(env_cols)) {
        missing <- setdiff(env_cols, env_in_pa)
        cat(sprintf("  ⚠ Missing env in PA: %s\n", paste(missing, collapse = ", ")))
        env_cols <- env_in_pa # Use only available vars
    }

    # ---- 4. Build train/test for each species ----
    species_counts <- po %>%
        group_by(spid) %>%
        summarise(n_po = n(), .groups = "drop") %>%
        arrange(desc(n_po))

    viable_species <- species_counts %>%
        filter(n_po >= min_po) %>%
        pull(spid)

    # Background data (shared)
    bg_data <- bg[, env_cols, drop = FALSE] %>% mutate(presence = 0)

    region_summary <- data.frame()

    for (sp in viable_species) {
        # Training: PO + BG
        sp_po <- po %>%
            filter(spid == sp) %>%
            select(all_of(env_cols)) %>%
            mutate(presence = 1)
        train <- bind_rows(sp_po, bg_data) %>% drop_na()

        # Testing: independent PA
        if (sp %in% names(pa_full)) {
            test <- pa_full %>%
                select(all_of(c(sp, env_cols))) %>%
                rename(presence = !!sp) %>%
                drop_na()
        } else {
            cat(sprintf("    ⚠ %s not in PA → skipping\n", sp))
            next
        }

        if (nrow(test) < 10 || sum(test$presence == 1) < 5) {
            cat(sprintf("    ⚠ %s too few test records → skipping\n", sp))
            next
        }

        write.csv(train, sprintf("output/case2/%s/train_data_%s.csv", region, sp),
            row.names = FALSE
        )
        write.csv(test, sprintf("output/case2/%s/test_data_%s.csv", region, sp),
            row.names = FALSE
        )

        region_summary <- rbind(region_summary, data.frame(
            region = region,
            species = sp,
            n_train = nrow(train),
            n_presence_train = sum(train$presence == 1),
            n_background = sum(train$presence == 0),
            n_test = nrow(test),
            n_presence_test = sum(test$presence == 1),
            n_absence_test = sum(test$presence == 0),
            prevalence_test = round(mean(test$presence == 1), 3),
            n_env_vars = length(env_cols),
            stringsAsFactors = FALSE
        ))
    }

    write.csv(region_summary,
        sprintf("output/case2/%s/species_summary.csv", region),
        row.names = FALSE
    )

    cat(sprintf(
        "  ✓ %d species prepared (of %d with PO≥%d)\n",
        nrow(region_summary), length(viable_species), min_po
    ))

    # Save env stats
    env_stats <- bg[, env_cols, drop = FALSE] %>%
        pivot_longer(everything(), names_to = "variable", values_to = "value") %>%
        group_by(variable) %>%
        summarise(
            mean = round(mean(value, na.rm = TRUE), 2),
            sd = round(sd(value, na.rm = TRUE), 2),
            min = min(value, na.rm = TRUE),
            max = max(value, na.rm = TRUE),
            .groups = "drop"
        )
    env_stats$region <- region
    write.csv(env_stats,
        sprintf("output/case2/%s/env_variable_stats.csv", region),
        row.names = FALSE
    )

    all_regions_summary <- rbind(all_regions_summary, region_summary)
}

# Also create SWI symlinks for backward compatibility
# (01_cast_pipeline.R used to read from output/case2/ directly)
swi_species <- all_regions_summary %>% filter(region == "SWI")
for (i in 1:nrow(swi_species)) {
    sp <- swi_species$species[i]
    src_train <- sprintf("output/case2/SWI/train_data_%s.csv", sp)
    dst_train <- sprintf("output/case2/train_data_%s.csv", sp)
    if (file.exists(src_train) && !file.exists(dst_train)) {
        file.copy(src_train, dst_train)
    }
    src_test <- sprintf("output/case2/SWI/test_data_%s.csv", sp)
    dst_test <- sprintf("output/case2/test_data_%s.csv", sp)
    if (file.exists(src_test) && !file.exists(dst_test)) {
        file.copy(src_test, dst_test)
    }
}
swi_summary_src <- "output/case2/SWI/species_summary.csv"
swi_summary_dst <- "output/case2/species_summary.csv"
if (file.exists(swi_summary_src) && !file.exists(swi_summary_dst)) {
    file.copy(swi_summary_src, swi_summary_dst)
}

write.csv(all_regions_summary, "output/case2/all_regions_summary.csv", row.names = FALSE)

# ==============================================================================
# Summary
# ==============================================================================
cat("\n======================================================================\n")
cat("  Data Preparation Complete!\n")
cat("======================================================================\n")

region_counts <- all_regions_summary %>%
    group_by(region) %>%
    summarise(n = n(), .groups = "drop")

for (i in 1:nrow(region_counts)) {
    cat(sprintf("  %s: %d species\n", region_counts$region[i], region_counts$n[i]))
}
cat(sprintf("  ────────────────\n"))
cat(sprintf(
    "  Total: %d species across %d regions\n",
    nrow(all_regions_summary), length(unique(all_regions_summary$region))
))

cat("\n  Output structure:\n")
for (r in regions) {
    n <- sum(all_regions_summary$region == r)
    cat(sprintf("    output/case2/%s/  (%d species)\n", r, n))
}
cat("    output/case2/all_regions_summary.csv\n")
cat("\n  Next: source('scripts/case2/01_cast_pipeline.R')  # single species\n")
cat("   or:  source('scripts/case2/02_multi_species_experiment.R')  # all regions\n")
