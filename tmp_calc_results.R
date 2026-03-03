library(dplyr)
library(tidyr)

d <- read.csv("E:/CausalSDMs/output/case2_eco/all_results_v3.csv", stringsAsFactors = FALSE) %>%
    filter(!is.na(auc_mean))

screen <- read.csv("E:/CausalSDMs/output/case2_eco/all_screening_v3.csv", stringsAsFactors = FALSE)

sp_meta <- read.csv("E:/CausalSDMs/outputs/EcoISEA3H/Res9/CAST_ready/CAST_Species_Summary.csv", stringsAsFactors = FALSE) %>%
    mutate(species = gsub(" ", "_", species))

d <- d %>% left_join(sp_meta %>% select(species, family), by = "species")
screen <- screen %>% left_join(sp_meta %>% select(species, family), by = "species")

# 3.1 & Table 2:
# Groups A and B: RF, MaxEnt, BRT, GAM do not exist in case2_eco based on previous scripts. Let's see what models actually exist.
cat("Available models:\n")
print(unique(d$model))

cat("\nTable 2 proxy (Mean AUC and TSS):\n")
tab2 <- d %>%
    group_by(model) %>%
    summarise(
        mean_auc = mean(auc_mean), sd_auc = sd(auc_mean),
        mean_tss = mean(tss_mean, na.rm = TRUE), sd_tss = sd(tss_mean, na.rm = TRUE)
    )
print(as.data.frame(tab2))

# Screening effect (A vs B)
# if models like RF_full, RF_cast exist, we compare them.
# Or if `model` is "CAST", "MLP", "MLP_ATE", "RF", "Maxent", "BRT". What are their variable sets?
# Let's see how `model` is defined.
# If they only have CAST, MLP_ATE, MLP, RF, BRT, Maxent, we can compare MLP vs CAST.

cat("\n3.2 Species count:\n")
cat("Total species in all_results_v3.csv: ", length(unique(d$species)), "\n")

paired_wide <- d %>%
    filter(model %in% c("MLP", "CAST")) %>%
    select(species, model, auc_mean) %>%
    pivot_wider(names_from = model, values_from = auc_mean, values_fn = max) %>%
    filter(!is.na(CAST), !is.na(MLP)) %>%
    mutate(delta = CAST - MLP)
n_cast_wins <- sum(paired_wide$delta > 0)
n_cast_ties <- sum(paired_wide$delta == 0)
n_total <- nrow(paired_wide)
cat(sprintf("\n3.2 CI-MLP vs MLP_cast wins: %d (%.1f%%)\n", n_cast_wins, 100 * n_cast_wins / n_total))
cat(sprintf("Average AUC advantage: %.4f\n", mean(paired_wide$delta)))

# Variables screening percentage
cast_screened <- screen %>%
    group_by(species) %>%
    summarise(n_total = n(), .groups = "drop")
cast_nvars <- d %>%
    filter(model == "CAST") %>%
    select(species, n_vars, family)
var_red <- cast_screened %>%
    inner_join(cast_nvars, by = "species") %>%
    mutate(red_pct = 1 - n_vars / n_total)

cat(sprintf("\n3.3 Variable reduction average: %.1f%%\n", mean(var_red$red_pct, na.rm = TRUE) * 100))

# DAG density processing
dag_info <- read.csv("E:/CausalSDMs/output/case2_eco/all_dag_info_v3.csv", stringsAsFactors = FALSE)
cat(sprintf("\n3.4 DAG Density average: %.2f\n", mean(dag_info$density, na.rm = TRUE)))
cat("Correlation Density vs AUC gain:\n")
dag_auc <- paired_wide %>% left_join(dag_info, by = "species")
cor_res <- cor.test(dag_auc$density, dag_auc$delta)
print(cor_res)

# Decomposing CAST (Screening vs Structure)
# In section 3.5, it mentions "A vs B" screening effect, which is FlatNN_full vs FlatNN_cast.
# Do we have MLP_full in the dataset? Let's check model names:
cat("Unique models:", paste(unique(d$model), collapse = ", "), "\n")

cat("\nTable 3 per Family / Region ? In Eco we have families instead of regions.\n")
tab3 <- var_red %>%
    group_by(family) %>%
    summarise(
        n_sp = n(),
        vif_vars = mean(n_total),
        cast_vars = mean(n_vars),
        red_pct = mean(red_pct) * 100
    )
print(as.data.frame(tab3))
