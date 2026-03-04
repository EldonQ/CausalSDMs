import os

in_file = "E:/CausalSDMs/scripts/EcoISEA3H/03_run_Eco_multi_species.R"
out_file = "E:/CausalSDMs/scripts/Plant/03_run_Plant_multi_species.R"

with open(in_file, "r", encoding="utf-8") as f:
    code = f.read()

# Replacements
code = code.replace('REGION <- "China_Res9"', 'REGION <- "Europe_Plant"')
code = code.replace('data_dir <- "E:/CausalSDMs/outputs/EcoISEA3H/Res9/CAST_ready/species_data_screened"', 
                    'data_dir <- "E:/CausalSDMs/outputs/Plant/CAST_ready/species_data_screened"')
code = code.replace('out_dir <- "E:/CausalSDMs/output/case2_eco"', 
                    'out_dir <- "E:/CausalSDMs/output/case4_plant"')

# Checkpoints
code = code.replace('file.path(out_dir, "all_results_v3.csv")', 'file.path(out_dir, "all_results_plant.csv")')
code = code.replace('file.path(out_dir, "all_ate_results_v3.csv")', 'file.path(out_dir, "all_ate_results_plant.csv")')
code = code.replace('file.path(out_dir, "all_dag_info_v3.csv")', 'file.path(out_dir, "all_dag_info_plant.csv")')
code = code.replace('file.path(out_dir, "all_dag_edges_v3.csv")', 'file.path(out_dir, "all_dag_edges_plant.csv")')
code = code.replace('file.path(out_dir, "all_screening_v3.csv")', 'file.path(out_dir, "all_screening_plant.csv")')
code = code.replace('file.path(out_dir, "all_role_info_v3.csv")', 'file.path(out_dir, "all_role_info_plant.csv")')

# Meta columns
code = code.replace('meta_cols <- c("HID", "lon", "lat", "species", "sid", "family", "category", "presence", "fraction")', 
                    'meta_cols <- c("PlotObservationID", "Longitude", "Latitude", "species", "presence")')

# Cosmetic naming
code = code.replace('Eco-ISEA3H Multi-Species Experiment', 'Plant Multi-Species Experiment')
code = code.replace('03_run_EcoISEA3H_multi_species.R', '03_run_Plant_multi_species.R')
code = code.replace('EcoISEA3H', 'Plant')

# Fix the split from "f <- sp_files[sp_idx]" to "CAST_.*_Res9_screened.csv"
code = code.replace('sp_name_raw <- gsub("CAST_|_Res9_screened\\\\.csv$", "", basename(f))', 
                    'sp_name_raw <- gsub("CAST_|_screened\\\\.csv$", "", basename(f))')

with open(out_file, "w", encoding="utf-8") as f:
    f.write(code)

print("Created 03_run_Plant_multi_species.R successfully.")
