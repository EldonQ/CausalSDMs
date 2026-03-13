# 02_vif_screening_envdata.R
# ============================================================
# Plant Dataset - VIF Screening
# Reference: E:\CausalSDMs\scripts\EcoISEA3H\02_vif_screening_envdata.R
# ============================================================
library(data.table)
library(dplyr)
library(car)

VIF_THRESHOLD <- 10
out_dir <- "E:/CausalSDMs/outputs/Plant/CAST_ready"
fig_dir <- "E:/CausalSDMs/figures/Plant/plot"
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

cat("================================================================\n")
cat("  02_vif_screening_envdata.R (Plant)\n")
cat("================================================================\n")

env <- fread(file.path(out_dir, "Plant_EnvData_Master.csv"))
meta_cols <- c("PlotObservationID", "Longitude", "Latitude")
all_env <- setdiff(names(env), meta_cols)

cat(sprintf("  Starting variables: %d\n", length(all_env)))

numeric_env <- all_env[sapply(env[, ..all_env], is.numeric)]

# Expert pre-filter (drop some highly redundant bio vars to speed up VIF)
expert_remove <- c(
    "bio_8", "bio_9", "bio_10", "bio_11", "bio_12",
    "bio_13", "bio_14", "bio_16", "bio_17"
)
expert_remove <- intersect(expert_remove, numeric_env)
remaining <- setdiff(numeric_env, expert_remove)

cat("\nStep 1: Iterative VIF screening (Threshold = 5)...\n")
env_clean <- env[, ..remaining]
env_clean <- na.omit(env_clean)

current_vars <- remaining
iteration <- 0

# Random subset to speed up lm
set.seed(42)
if (nrow(env_clean) > 8000) {
    env_clean_samp <- env_clean[sample(nrow(env_clean), 8000), ]
} else {
    env_clean_samp <- env_clean
}

env_clean_samp$DUMMY <- rnorm(nrow(env_clean_samp))

repeat {
    iteration <- iteration + 1
    fmla <- as.formula(paste("DUMMY ~", paste(current_vars, collapse = " + ")))

    mod <- lm(fmla, data = env_clean_samp)
    al <- alias(mod)
    if (!is.null(al$Complete)) {
        aliased_var <- rownames(al$Complete)[1]
        cat(sprintf("  [!] Perfect collinearity: removing %s\n", aliased_var))
        current_vars <- setdiff(current_vars, aliased_var)
        next
    }

    vif_vals <- tryCatch(
        {
            vif(mod)
        },
        error = function(e) NULL
    )
    if (is.null(vif_vals)) break

    max_vif <- max(vif_vals)
    max_var <- names(which.max(vif_vals))

    if (max_vif <= VIF_THRESHOLD) {
        cat(sprintf("  ✅ Converged: Max VIF = %.2f\n", max_vif))
        break
    }

    cat(sprintf("  Iter %2d: remove %-20s (VIF=%.1f)\n", iteration, max_var, max_vif))
    current_vars <- setdiff(current_vars, max_var)
    if (length(current_vars) < 5) break
}

final_vars <- current_vars
cat(sprintf("\n  Final variable set: %d variables\n", length(final_vars)))
print(final_vars)

fwrite(data.frame(variable = final_vars), file.path(out_dir, "VIF_Final_Variables.csv"))

cat("\nStep 2: Saving screened datasets...\n")
env_screened <- env[, c("PlotObservationID", "Longitude", "Latitude", final_vars), with = FALSE]
fwrite(env_screened, file.path(out_dir, "Plant_EnvData_Screened.csv"))

sp_dir_in <- file.path(out_dir, "species_data")
sp_dir_out <- file.path(out_dir, "species_data_screened")
dir.create(sp_dir_out, showWarnings = FALSE)

sp_smry <- fread(file.path(out_dir, "CAST_Species_Summary.csv"))
for (i in 1:nrow(sp_smry)) {
    fin <- file.path(sp_dir_in, sp_smry$file[i])
    if (file.exists(fin)) {
        sp_df <- fread(fin)
        keep <- c("PlotObservationID", "Longitude", "Latitude", "species", "presence", final_vars)
        sp_df <- sp_df[, ..keep]
        fout <- file.path(sp_dir_out, gsub("\\.csv", "_screened.csv", sp_smry$file[i]))
        fwrite(sp_df, fout)
    }
}
cat("================================================================\n")
cat(" Done. Proceed to 03_run_Plant_multi_species.R\n")
cat("================================================================\n")
