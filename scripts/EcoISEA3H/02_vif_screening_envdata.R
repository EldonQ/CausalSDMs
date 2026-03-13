# 02_vif_screening_envdata.R
# ============================================================
# Eco-ISEA3H | Resolution 9 | Global Region
# VIF-based variable screening before CAST modeling
#
# Strategy (3-stage):
#   Stage 1: Expert-guided pre-filter — remove obvious redundancies
#   Stage 2: Iterative VIF elimination (threshold = 10)
#   Stage 3: Rebuild per-species CAST-ready datasets with clean vars
#
# Input:  outputs/EcoISEA3H/Res9/CAST_ready/Global_EnvData_Res9_Master.csv
# Output: outputs/EcoISEA3H/Res9/CAST_ready/Global_EnvData_Res9_Screened.csv
#         outputs/EcoISEA3H/Res9/CAST_ready/species_data_screened/
# ============================================================

library(data.table)
library(dplyr)
library(car) # vif()
library(ggplot2)
library(corrplot)

# ── Configurable ─────────────────────────────────────────
VIF_THRESHOLD <- 10 # Threshold 5 follows Zuur et al. (2010) to ensure robust causal DAG learning and structure identifiability.
MIN_CELLS <- 200 # Same as 01_prepare
# ─────────────────────────────────────────────────────────

out_dir <- "E:/CausalSDMs/outputs/EcoISEA3H/Global_Res9/CAST_ready"
fig_dir <- "E:/CausalSDMs/figures/EcoISEA3H/Global_Res9"
sp_dir_in <- file.path(out_dir, "species_data")
sp_dir_out <- file.path(out_dir, "species_data_screened")
dir.create(sp_dir_out, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat("  02_vif_screening_envdata.R\n")
cat(sprintf("  VIF threshold = %d\n", VIF_THRESHOLD))
cat("================================================================\n")

# ============================================================
# STEP 1: Load and expert pre-filter
# ============================================================
cat("\nStep 1: Loading env data & expert pre-filter...\n")

env <- fread(file.path(out_dir, "Global_EnvData_Res9_Master.csv"))
meta_cols <- c("HID", "lon", "lat")
all_env <- setdiff(names(env), meta_cols)
cat(sprintf("  Starting variables: %d\n", length(all_env)))

# ── Stage 1a: Remove categorical (not suitable for VIF) ──
cat_cols <- all_env[!sapply(env[, ..all_env], is.numeric)]
cat(sprintf(
    "  Categorical (set aside): %s\n",
    ifelse(length(cat_cols) > 0, paste(cat_cols, collapse = ", "), "none")
))

numeric_env <- setdiff(all_env, cat_cols)

# ── Stage 1b: Expert redundancy removal ──────────────────
# Rule: when two variables measure essentially the same thing,
# keep the one that is more standard / ecologically interpretable
expert_remove <- c(
    # --- Temperature redundancy ---
    # growingdegdays0/5 = accumulated temperature, highly redundant with bio01
    "growingdegdays0", "growingdegdays5",
    # thermicityindex = (bio06 + bio11 + bio01)*10, pure linear combo
    "thermicityindex",
    # bio08/bio09 (temp of wettest/driest quarter) = bio subset with noise
    "bio08", "bio09",
    # bio10/bio11 (temp of warmest/coldest quarter) - keep bio05/bio06 instead
    "bio10", "bio11",

    # --- ETCCDI temperature extremes redundant with BIO ---
    # etccdi_tnn (min of Tmin) ≈ bio06; etccdi_txx (max of Tmax) ≈ bio05
    "etccdi_tnn", "etccdi_txx", "etccdi_tnx", "etccdi_txn",
    # etccdi_fd (frost days) ≈ bio06; etccdi_id (icing days) ≈ bio06
    "etccdi_id",
    # etccdi_su (summer days), etccdi_tr (tropical nights) ≈ bio05
    "etccdi_su", "etccdi_tr",
    # etccdi_tn10p/tn90p/tx10p/tx90p = percentile indices, secondary
    "etccdi_tn10p", "etccdi_tn90p", "etccdi_tx10p", "etccdi_tx90p",
    # etccdi_wsdi/csdi (warm/cold spell) = duration metrics, secondary
    "etccdi_wsdi", "etccdi_csdi",
    # etccdi_dtr (diurnal T range) ≈ bio02
    "etccdi_dtr",

    # --- Precipitation redundancy ---
    # etccdi_prcptot ≈ bio12; etccdi_r1mm ≈ bio12
    "etccdi_prcptot", "etccdi_r1mm",
    # etccdi_r10mm/r20mm = threshold indices, keep rx1day/rx5day instead
    "etccdi_r10mm", "etccdi_r20mm",
    # etccdi_r95p/r99p = extreme precip, redundant with rx1day
    "etccdi_r95p", "etccdi_r99p",

    # --- ENVIREM redundancy ---
    # petcoldestquarter/petwarmestquarter/petwettestquarter/petdriestquarter
    # = quarterly PET, keep annualpet + petseasonality instead
    "petcoldestquarter", "petwarmestquarter", "petwettestquarter", "petdriestquarter",
    # monthcountbytemp10 ≈ growingdegdays / bio01
    "monthcountbytemp10",
    # maxtempcoldestest / mintempwarmest = extremes already in bio05/06
    "maxtempcoldestest", "mintempwarmest",
    # continentality ≈ bio07 (annual temperature range)
    "continentality"
)

# Only remove what actually exists
expert_remove <- intersect(expert_remove, numeric_env)
cat(sprintf("  Expert pre-filter removes: %d variables\n", length(expert_remove)))
cat("    Removed: ", paste(expert_remove, collapse = ", "), "\n")

remaining <- setdiff(numeric_env, expert_remove)
cat(sprintf("  Remaining for VIF: %d variables\n", length(remaining)))

# ============================================================
# STEP 2: Iterative VIF elimination
# ============================================================
cat("\nStep 2: Iterative VIF screening (threshold = %d)...\n", VIF_THRESHOLD)

# Prepare clean data (no NAs for lm)
env_clean <- env[, c("HID", remaining), with = FALSE]
env_clean <- na.omit(env_clean)
cat(sprintf("  Complete cases: %d / %d\n", nrow(env_clean), nrow(env)))

# Iterative VIF removal
current_vars <- remaining
iteration <- 0
vif_log <- list()

repeat {
    iteration <- iteration + 1

    # Fit dummy regression to get VIF
    fmla <- as.formula(paste("HID ~", paste(current_vars, collapse = " + ")))

    # First check for perfect collinearity (aliased coefficients)
    mod <- lm(fmla, data = env_clean)
    al <- alias(mod)
    if (!is.null(al$Complete)) {
        # Extract the first aliased variable name
        aliased_var <- rownames(al$Complete)[1]
        cat(sprintf("  [!] Perfect collinearity detected. Removing: %s\n", aliased_var))
        current_vars <- setdiff(current_vars, aliased_var)
        next
    }

    tryCatch(
        {
            vif_vals <- vif(mod)
        },
        error = function(e) {
            cat(sprintf("  [!] VIF error at iteration %d: %s\n", iteration, e$message))
            # Try removing near-constant columns
            sds <- sapply(env_clean[, ..current_vars], sd, na.rm = TRUE)
            if (any(sds < 1e-10)) {
                bad <- names(sds[sds < 1e-10])
                cat(sprintf("  Removing near-constant: %s\n", paste(bad, collapse = ", ")))
                current_vars <<- setdiff(current_vars, bad)
            }
            return(NULL)
        }
    )

    if (!exists("vif_vals") || is.null(vif_vals)) next

    max_vif <- max(vif_vals)
    max_var <- names(which.max(vif_vals))

    vif_log[[iteration]] <- data.frame(
        iteration = iteration,
        n_vars = length(current_vars),
        max_vif = round(max_vif, 2),
        removed = max_var,
        stringsAsFactors = FALSE
    )

    if (max_vif <= VIF_THRESHOLD) {
        cat(sprintf(
            "\n  ✅ Converged at iteration %d: all VIF <= %d (%d variables remain)\n",
            iteration, VIF_THRESHOLD, length(current_vars)
        ))
        break
    }

    cat(sprintf(
        "  Iter %2d: %2d vars | max VIF = %8.1f -> remove %-30s\n",
        iteration, length(current_vars), max_vif, max_var
    ))
    current_vars <- setdiff(current_vars, max_var)

    if (length(current_vars) < 5) {
        cat("  [!] Too few variables remaining. Stopping.\n")
        break
    }
}

final_vars <- current_vars
cat(sprintf("\n  Final variable set: %d variables\n", length(final_vars)))

# Add back categorical if present
if (length(cat_cols) > 0) {
    final_vars_full <- c(final_vars, cat_cols)
} else {
    final_vars_full <- final_vars
}

# Print final VIF values
cat("\n  Final VIF values:\n")
fmla_f <- as.formula(paste("HID ~", paste(final_vars, collapse = " + ")))
mod_f <- lm(fmla_f, data = env_clean)
vif_f <- sort(vif(mod_f), decreasing = TRUE)
for (nm in names(vif_f)) {
    cat(sprintf("    %-40s VIF = %6.2f\n", nm, vif_f[nm]))
}

# Save VIF screening log
vif_log_df <- bind_rows(vif_log)
fwrite(vif_log_df, file.path(out_dir, "VIF_Screening_Log.csv"))

# ============================================================
# STEP 3: Group retained variables by ecological dimension
# ============================================================
cat("\nStep 3: Ecological dimension grouping...\n")

classify_var <- function(v) {
    if (grepl("^bio0[1-7]$|^bio1[01]$|temp|tnn|txx|fd|gsl|su", v, ignore.case = TRUE)) {
        return("Temperature")
    }
    if (grepl("^bio1[2-9]$|prec|rain|rx|cdd|cwd|sdii|r1|r10|r20|r95|r99", v, ignore.case = TRUE)) {
        return("Precipitation")
    }
    if (grepl("pet|arid|moisture|ember", v, ignore.case = TRUE)) {
        return("Water-Energy")
    }
    if (grepl("elev|tri|topo|slope", v, ignore.case = TRUE)) {
        return("Topography")
    }
    if (grepl("tree|bare|nontree|landcover|vcf", v, ignore.case = TRUE)) {
        return("Vegetation/Land")
    }
    if (grepl("etccdi", v, ignore.case = TRUE)) {
        return("Climate Extremes")
    }
    return("Other")
}

var_groups <- data.frame(
    variable = final_vars_full,
    dimension = sapply(final_vars_full, classify_var),
    stringsAsFactors = FALSE
) %>% arrange(dimension, variable)

cat("\n  Retained variables by ecological dimension:\n")
for (dim in unique(var_groups$dimension)) {
    vars_in_dim <- var_groups$variable[var_groups$dimension == dim]
    cat(sprintf(
        "    [%s] (%d): %s\n", dim, length(vars_in_dim),
        paste(vars_in_dim, collapse = ", ")
    ))
}

fwrite(var_groups, file.path(out_dir, "VIF_Final_Variables.csv"))

# ============================================================
# STEP 4: Save screened env data & rebuild species datasets
# ============================================================
cat("\nStep 4: Saving screened data...\n")

# Screened master env table
env_screened <- env[, c("HID", "lon", "lat", final_vars_full), with = FALSE]
fwrite(env_screened, file.path(out_dir, "Global_EnvData_Res9_Screened.csv"))
cat(sprintf(
    "  -> Global_EnvData_Res9_Screened.csv (%d x %d)\n",
    nrow(env_screened), ncol(env_screened)
))

# Rebuild per-species datasets
cat("  Rebuilding per-species CAST-ready datasets...\n")

sp_summary <- fread(file.path(out_dir, "CAST_Species_Summary.csv"))

for (i in seq_len(nrow(sp_summary))) {
    old_file <- file.path(sp_dir_in, sp_summary$file[i])
    if (!file.exists(old_file)) next

    sp_df <- fread(old_file)
    # Keep meta + screened variables
    keep_cols <- c(
        "HID", "lon", "lat", "species", "sid", "family", "category",
        "presence", "fraction", final_vars_full
    )
    keep_cols <- intersect(keep_cols, names(sp_df))
    sp_df <- sp_df[, ..keep_cols]

    new_file <- file.path(
        sp_dir_out,
        gsub("\\.csv$", "_screened.csv", sp_summary$file[i])
    )
    fwrite(sp_df, new_file)
}

cat(sprintf("  -> %d species rebuilt in %s\n", nrow(sp_summary), sp_dir_out))

# ============================================================
# STEP 5: Visualization
# ============================================================
cat("\nStep 5: Generating figures...\n")

# 5a. Correlation heatmap of final variables
cor_mat <- cor(env_clean[, ..final_vars], use = "pairwise.complete.obs")

png(file.path(fig_dir, "06_VIF_Screened_Correlation.png"),
    width = 1200, height = 1000, res = 150
)
corrplot(cor_mat,
    method = "color", type = "upper",
    tl.cex = 0.6, tl.col = "black",
    col = colorRampPalette(c("#2166AC", "white", "#B2182B"))(200),
    title = sprintf(
        "Env Variable Correlations After VIF Screening (%d vars)",
        length(final_vars)
    ),
    mar = c(0, 0, 2, 0)
)
dev.off()
cat("  -> 06_VIF_Screened_Correlation.png\n")

# 5b. VIF reduction trajectory
if (nrow(vif_log_df) > 1) {
    p_vif <- ggplot(vif_log_df, aes(x = iteration, y = max_vif)) +
        geom_line(color = "#2980B9", linewidth = 1) +
        geom_point(color = "#2980B9", size = 2) +
        geom_hline(
            yintercept = VIF_THRESHOLD, linetype = "dashed",
            color = "#E74C3C", linewidth = 0.6
        ) +
        annotate("text",
            x = max(vif_log_df$iteration) * 0.7, y = VIF_THRESHOLD + 5,
            label = sprintf("Threshold = %d", VIF_THRESHOLD),
            color = "#E74C3C", fontface = "bold", size = 3.5
        ) +
        scale_y_log10() +
        theme_minimal(base_size = 11) +
        labs(
            title = "VIF Screening Convergence",
            subtitle = sprintf(
                "%d -> %d variables",
                vif_log_df$n_vars[1], length(final_vars)
            ),
            x = "Iteration", y = "Max VIF (log scale)"
        )
    ggsave(file.path(fig_dir, "07_VIF_Convergence.png"), p_vif,
        width = 8, height = 5, dpi = 300
    )
    cat("  -> 07_VIF_Convergence.png\n")
}

# 5c. Before/After variable count comparison
before_after <- data.frame(
    stage = factor(c("Original (69)", "Expert filter", "VIF screened"),
        levels = c("Original (69)", "Expert filter", "VIF screened")
    ),
    n_vars = c(length(numeric_env), length(remaining), length(final_vars))
)

p_ba <- ggplot(before_after, aes(x = stage, y = n_vars, fill = stage)) +
    geom_col(width = 0.6) +
    geom_text(aes(label = n_vars), vjust = -0.5, fontface = "bold", size = 4) +
    scale_fill_manual(values = c("#E74C3C", "#F39C12", "#27AE60"), guide = "none") +
    theme_minimal(base_size = 12) +
    labs(
        title = "Variable Screening Pipeline",
        subtitle = sprintf(
            "69 → %d → %d variables (%.0f%% reduction)",
            length(remaining), length(final_vars),
            100 * (1 - length(final_vars) / length(numeric_env))
        ),
        x = NULL, y = "Number of variables"
    )
ggsave(file.path(fig_dir, "08_Variable_Screening_Pipeline.png"), p_ba,
    width = 7, height = 5, dpi = 300
)
cat("  -> 08_Variable_Screening_Pipeline.png\n")

# ── Final summary ────────────────────────────────────────
cat("\n================================================================\n")
cat("=== VIF Screening Complete ===\n")
cat("================================================================\n")
cat(sprintf("  Original variables:   %d\n", length(numeric_env)))
cat(sprintf("  Expert removed:       %d\n", length(expert_remove)))
cat(sprintf("  VIF removed:          %d\n", length(remaining) - length(final_vars)))
cat(sprintf(
    "  Final variables:      %d (+ %d categorical)\n",
    length(final_vars), length(cat_cols)
))
cat(sprintf(
    "  Max VIF:              %.2f (threshold = %d)\n",
    max(vif_f), VIF_THRESHOLD
))
cat(sprintf(
    "  Screened env table:   %s\n",
    file.path(out_dir, "Global_EnvData_Res9_Screened.csv")
))
cat(sprintf("  Species data dir:     %s\n", sp_dir_out))
cat("================================================================\n")
