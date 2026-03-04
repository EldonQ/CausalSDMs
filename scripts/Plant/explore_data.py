import pandas as pd
import numpy as np
import os

data_dir = "E:/CausalSDMs/maskSDMData/data"

print("--- CSV Files ---")
csv_files = ["location_data.csv", "species_names.csv", "cover_data.csv", "worldclim_data.csv", "height_data.csv", "human_data.csv", "metadata_data.csv", "soilgrid_data.csv", "topographic_data.csv"]
for f in csv_files:
    try:
        df = pd.read_csv(os.path.join(data_dir, f))
        print(f"{f}: shape={df.shape}, cols={list(df.columns)[:5]}...")
    except Exception as e:
        print(f"Error reading {f}: {e}")

print("\n--- NPY Files ---")
npy_files = ["europe_map_coordinates.npy", "europe_map_human.npy", "europe_map_soilgrids.npy", "europe_map_topography.npy", "europe_map_worldclim.npy", "species_occurrences.npy"]
for f in npy_files:
    try:
        arr = np.load(os.path.join(data_dir, f), mmap_mode='r')
        print(f"{f}: shape={arr.shape}, dtype={arr.dtype}")
    except Exception as e:
        print(f"Error reading {f}: {e}")
