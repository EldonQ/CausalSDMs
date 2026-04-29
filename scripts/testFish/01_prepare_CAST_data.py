from __future__ import annotations

"""
01_prepare_CAST_data.py
=======================
testFish | Lota lota | RiverATLAS v10 | Global

Build a CAST-ready dataset for Lota lota species distribution modeling.

Workflow (mirrors EcoISEA3H/01_prepare_CAST_envdata.R):
  1. Load species occurrence data (Lota_lota.csv)
  2. Load ALL RiverATLAS v10 shapefiles (global coverage)
     - Extract reach geometry and all environmental variables
     - 10 regional shapefiles: af, ar, as, au, eu, gr, na, sa_north, sa_south, si
  3. Spatial join: match species points to nearest river reaches
  4. Presence / pseudo-absence construction
  5. Merge reach-level environmental variables
  6. Save CAST-ready master dataset + per-species CSV

References
---------
RiverATLAS v10 — HydroATLAS: Lehner & Grill (2013)
  https://www.hydrosheds.org/page/hydroatlas
  Creative Commons CC-BY 4.0
"""

import os
import sys
from pathlib import Path

# Fix GDAL DLL loading on Windows (conda envs store GDAL DLLs in Library/bin)
dll_dir = os.path.join(os.path.dirname(sys.executable), "Library", "bin")
if os.path.isdir(dll_dir) and dll_dir not in os.environ.get("PATH", ""):
    os.add_dll_directory(dll_dir)

import geopandas as gpd
import numpy as np
import pandas as pd
from shapely.geometry import Point
from sklearn.neighbors import BallTree

# ── Paths ──────────────────────────────────────────────────────────────────────
SCRIPT_DIR = Path(__file__).parent.resolve()
DATA_DIR   = SCRIPT_DIR                          # Lota_lota.csv lives here
RIVER_DIR  = DATA_DIR / "RiverATLAS_Data_v10_shp" / "RiverATLAS_v10_shp"
OUT_DIR    = SCRIPT_DIR / "output"
FIG_DIR    = SCRIPT_DIR / "figures"
OUT_DIR.mkdir(exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)

# RiverATLAS shapefiles (regional, together cover the globe)
SHP_FILES = {
    "af":       RIVER_DIR / "RiverATLAS_v10_af.shp",   # Africa
    "ar":       RIVER_DIR / "RiverATLAS_v10_ar.shp",   # Arctic
    "as":       RIVER_DIR / "RiverATLAS_v10_as.shp",   # Asia
    "au":       RIVER_DIR / "RiverATLAS_v10_au.shp",   # Australasia
    "eu":       RIVER_DIR / "RiverATLAS_v10_eu.shp",   # Europe
    "gr":       RIVER_DIR / "RiverATLAS_v10_gr.shp",   # Greenland
    "na":       RIVER_DIR / "RiverATLAS_v10_na.shp",   # North America
    "sa_north": RIVER_DIR / "RiverATLAS_v10_sa_north.shp",  # South America North
    "sa_south": RIVER_DIR / "RiverATLAS_v10_sa_south.shp",  # South America South
    "si":       RIVER_DIR / "RiverATLAS_v10_si.shp",   # South Asia
}

# ── Configurable parameters ────────────────────────────────────────────────────
MAX_JOIN_DIST_KM = 50       # Max distance (km) for point→reach spatial join
                         # Reaches > 50 km from any point carry no info
PSEUDOABS_RATIO = 3         # N pseudo-absences per presence
RANDOM_SEED = 42

np.random.seed(RANDOM_SEED)


# ══════════════════════════════════════════════════════════════════════════════
# HELPER: load a shapefile, extract non-spatial columns, return clean GeoDataFrame
# ══════════════════════════════════════════════════════════════════════════════
def load_shapefile(shp_path: Path) -> gpd.GeoDataFrame:
    """Load a RiverATLAS shapefile and return a clean GeoDataFrame."""
    print(f"    Loading {shp_path.name} ...", end=" ", flush=True)
    gdf = gpd.read_file(str(shp_path))
    # RiverATLAS stores geometry as river-line or reach-centroid
    # We keep geometry so we can compute distances
    n = len(gdf)
    print(f"{n:,} reaches loaded")
    return gdf


def identify_env_columns(gdf: gpd.GeoDataFrame) -> list:
    """
    RiverATLAS v10 columns include administrative IDs, geometry, and
    attribute variables.  Return the list of env-attribute column names
    (excluding metadata, geometry, and administrative fields).
    """
    META_EXCLUDE = {
        "geometry", "Hylak_id", "NextDownID", "MainRiver", "DistDown",
        "DistSea",  "CatchStorM", "VolR", "UpCumulS", "Shape_Leng", "Shape_Area",
        "geometry_type", "gid",
    }
    cols = [c for c in gdf.columns
            if c not in META_EXCLUDE
            and not c.lower().startswith("shape")
            and not c.lower().startswith("geom")]
    return cols


# ══════════════════════════════════════════════════════════════════════════════
# STEP 1: Load species occurrence data
# ══════════════════════════════════════════════════════════════════════════════
def load_species_data(csv_path: Path) -> pd.DataFrame:
    """Load and validate the species occurrence CSV."""
    print(f"\n[Step 1] Loading species data: {csv_path}")
    df = pd.read_csv(csv_path)
    required = {"taxon", "longitude", "latitude", "occurrenceStatus"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CSV缺少必要列: {missing}")

    # Keep only presence records for now
    df_pres = df[df["occurrenceStatus"].str.upper() == "PRESENT"].copy()
    df_pres = df_pres.dropna(subset=["longitude", "latitude"])

    # Remove obvious ocean outliers (lon/lat on land)
    df_pres = df_pres[
        (df_pres["longitude"].between(-180, 180)) &
        (df_pres["latitude"].between(-90, 90))
    ].reset_index(drop=True)

    print(f"    Total records : {len(df):,}")
    print(f"    Present records: {len(df_pres):,}")
    return df_pres


# ══════════════════════════════════════════════════════════════════════════════
# STEP 2: Load ALL RiverATLAS shapefiles and concatenate
# ══════════════════════════════════════════════════════════════════════════════
def load_global_riveratlas() -> tuple[gpd.GeoDataFrame, list[str]]:
    """
    Load all 10 RiverATLAS regional shapefiles, concatenate into one
    global GeoDataFrame, and return a list of env-var column names.
    """
    print(f"\n[Step 2] Loading RiverATLAS v10 global river network")
    print(f"    Regions: {', '.join(SHP_FILES.keys())}")

    all_gdfs = []
    env_cols_global = set()

    for region, shp_path in SHP_FILES.items():
        if not shp_path.exists():
            print(f"    [WARN] {shp_path.name} not found — skipping {region}")
            continue

        gdf = load_shapefile(shp_path)
        env_cols = identify_env_columns(gdf)
        env_cols_global.update(env_cols)

        # Tag with region for diagnostics
        gdf["src_region"] = region
        all_gdfs.append(gdf)

    if not all_gdfs:
        raise RuntimeError("No RiverATLAS shapefiles could be loaded!")

    global_gdf = gpd.GeoDataFrame(pd.concat(all_gdfs, ignore_index=True), crs=all_gdfs[0].crs)
    env_cols_sorted = sorted(env_cols_global)

    print(f"\n    Global reaches: {len(global_gdf):,}")
    print(f"    Environmental columns: {len(env_cols_sorted)}")
    print(f"    CRS: {global_gdf.crs}")
    return global_gdf, env_cols_sorted


# ══════════════════════════════════════════════════════════════════════════════
# STEP 3: Spatial join — match species points to nearest river reaches
# ══════════════════════════════════════════════════════════════════════════════
def spatial_join_points_to_reaches(
    species_df: pd.DataFrame,
    river_gdf:  gpd.GeoDataFrame,
    max_dist_km: float = MAX_JOIN_DIST_KM,
) -> pd.DataFrame:
    """
    For each species point, find the nearest river reach.
    Uses BallTree on lat/lon (great-circle distance) for speed.

    Returns a DataFrame with columns:
      lon, lat, reach_idx, distance_km, + env vars from the matched reach
    """
    print(f"\n[Step 3] Spatial join — species points → river reaches (max_dist={max_dist_km} km)")

    # Convert species points to a GeoDataFrame in WGS84
    species_gdf = gpd.GeoDataFrame(
        species_df.copy(),
        geometry=[Point(lon, lat) for lon, lat in
                  zip(species_df["longitude"], species_df["latitude"])],
        crs="EPSG:4326"
    )

    # Project river reaches to a metric CRS for accurate buffering
    # Use EPSG:4326 + degree→km conversion handled by BallTree
    # For the join itself, work in metres on a PlateCarree-projected CRS
    river_proj = river_gdf.to_crs("EPSG:4326")   # long/lat but same crs as points
    # Actually, use a metric CRS centred near the study area (global: eqc)
    # EPSG:6933 = Cylindrical Equal Area (global, metric)
    river_eqc = river_gdf.to_crs("EPSG:6933")
    species_eqc = species_gdf.to_crs("EPSG:6933")

    reach_centroids = np.radians(
        np.array([(g.centroid.x, g.centroid.y) for g in river_eqc.geometry])
    )
    point_coords = np.radians(
        np.array([(g.centroid.x, g.centroid.y) for g in species_eqc.geometry])
    )

    # Earth radius in km
    EARTH_R = 6371.0

    print(f"    Building BallTree for {len(reach_centroids):,} reaches ...")
    tree = BallTree(reach_centroids, metric="haversine")

    max_dist_rad = max_dist_km / EARTH_R
    print(f"    Querying nearest reaches (k=1, radius={max_dist_km} km) ...")
    dist_rad, idx_rad = tree.query(point_coords, k=1, return_distance=True)
    dist_km = (dist_rad * EARTH_R).flatten()

    matched = idx_rad.flatten()
    kept = dist_km <= max_dist_km

    n_matched = kept.sum()
    n_unmatched = (~kept).sum()
    print(f"    Matched : {n_matched:,} points ({100*n_matched/len(kept):.1f}%)")
    print(f"    Unmatched (> {max_dist_km} km): {n_unmatched:,}")

    # Build result table
    result = species_df.copy()
    result["reach_idx"] = matched
    result["dist_km"] = np.round(dist_km, 3)
    result["matched"] = kept

    return result


def attach_env_data(
    joined_df:  pd.DataFrame,
    river_gdf:   gpd.GeoDataFrame,
    env_cols:    list[str],
) -> pd.DataFrame:
    """
    Attach environmental variables from the matched reach to each species point.
    """
    print(f"\n[Step 3b] Attaching {len(env_cols)} environmental variables ...")

    # Reset index so reach_idx maps correctly
    river_reset = river_gdf.reset_index(drop=True)

    # Filter to matched points
    matched_df = joined_df[joined_df["matched"]].copy()
    n_m = len(matched_df)

    # Pull env data for matched reaches
    reach_indices = matched_df["reach_idx"].values
    env_data = river_reset.loc[reach_indices, env_cols].reset_index(drop=True)

    # Replace NoData sentinel (-9999) with NA
    for col in env_cols:
        if pd.api.types.is_numeric_dtype(env_data[col]):
            env_data[col] = env_data[col].replace(-9999, np.nan)

    for col in env_cols:
        joined_df.loc[joined_df["matched"], col] = env_data[col].values

    # Report NA coverage
    na_pct = joined_df[env_cols].isna().mean() * 100
    high_na = na_pct[na_pct > 10].sort_values(ascending=False)
    print(f"    Variables with >10% NA after join:")
    if len(high_na) > 0:
        for nm, val in high_na.head(10).items():
            print(f"      {nm:<35s} {val:5.1f}%")
    else:
        print("      None — good coverage!")

    return joined_df


# ══════════════════════════════════════════════════════════════════════════════
# STEP 4: Presence / Pseudo-absence construction
# ══════════════════════════════════════════════════════════════════════════════
def build_presence_pseudoabsence(
    matched_df: pd.DataFrame,
    river_gdf:  gpd.GeoDataFrame,
    env_cols:   list[str],
    ratio:      int = PSEUDOABS_RATIO,
) -> pd.DataFrame:
    """
    - Presence: species points matched to a reach (dist_km <= threshold)
    - Pseudo-absence: randomly sample reaches that are:
        (a) far from any presence point (> max_dist_km), and
        (b) NOT matched to any presence point
    """
    print(f"\n[Step 4] Building presence / pseudo-absence dataset (ratio = {ratio}:1)")

    # ── Presence ──────────────────────────────────────────────────────────────
    pres_df = matched_df[matched_df["matched"]].copy()
    pres_df["presence"] = 1

    n_pres = len(pres_df)
    n_pseudo = n_pres * ratio
    print(f"    Presence points: {n_pres:,}")

    # ── Pseudo-absence ────────────────────────────────────────────────────────
    # Reaches matched to presence
    matched_reach_idx = set(pres_df["reach_idx"].values)

    # Project to metric CRS
    river_eqc = river_gdf.to_crs("EPSG:6933")

    # Presence reach centroids (in metres)
    pres_centroids = np.array([
        (river_eqc.loc[i].geometry.centroid.x,
         river_eqc.loc[i].geometry.centroid.y)
        for i in pres_df["reach_idx"]
    ])

    # Build tree over presence centroids
    pres_tree = BallTree(pres_centroids)

    # Candidate pseudo-absence reaches: those NOT in matched set
    candidate_mask = ~river_eqc.index.isin(matched_reach_idx)
    candidate_idx = river_eqc.index[candidate_mask].tolist()
    n_candidates = len(candidate_idx)
    print(f"    Candidate pseudo-absence reaches: {n_candidates:,}")

    if n_candidates == 0:
        print("    [WARN] No candidate pseudo-absence reaches found.")
        return pres_df

    # Sample candidates (up to n_pseudo random reaches)
    n_to_sample = min(n_pseudo, n_candidates)
    sampled_idx = np.random.choice(candidate_idx, size=n_to_sample, replace=False)

    # Filter to those far enough (> MAX_JOIN_DIST_KM * 2) from any presence reach
    sampled_coords = np.array([
        (river_eqc.loc[i].geometry.centroid.x,
         river_eqc.loc[i].geometry.centroid.y)
        for i in sampled_idx
    ])
    # Convert to radians for haversine
    sampled_coords_rad = np.radians(
        gpd.GeoSeries(
            [Point(x, y) for x, y in sampled_coords],
            crs="EPSG:6933"
        ).to_crs("EPSG:4326").toadians()
    )
    pres_centroids_rad = np.radians(
        gpd.GeoSeries(
            [Point(x, y) for x, y in pres_centroids],
            crs="EPSG:6933"
        ).to_crs("EPSG:4326").toadians()
    )
    dist_rad, _ = pres_tree.query(sampled_coords_rad, k=1, return_distance=True)
    min_dist_km = (dist_rad.flatten() * 6371.0)

    # Keep only reaches > 2 * max_dist from any presence
    safe_mask = min_dist_km > (MAX_JOIN_DIST_KM * 2)
    final_pseudo_idx = sampled_idx[safe_mask]

    print(f"    After distance filter (> 2×{MAX_JOIN_DIST_KM} km from presence): {len(final_pseudo_idx):,}")

    # Cap to desired ratio
    final_pseudo_idx = final_pseudo_idx[:n_pres * ratio]

    # Build pseudo-absence DataFrame
    pseudo_river = river_eqc.loc[final_pseudo_idx].reset_index(drop=True)
    pseudo_lonlat = (
        pseudo_river
        .to_crs("EPSG:4326")
        .geometry.centroid
    )

    pseudo_df = pd.DataFrame({
        "longitude": pseudo_lonlat.x.values,
        "latitude":  pseudo_lonlat.y.values,
        "presence":  0,
        "reach_idx": final_pseudo_idx,
        "dist_km":   np.nan,
        "matched":   False,
    })
    for col in env_cols:
        pseudo_df[col] = river_eqc.loc[final_pseudo_idx, col].reset_index(drop=True).values
        # Replace NoData sentinel
        if pd.api.types.is_numeric_dtype(pseudo_df[col].dtype):
            pseudo_df[col] = pseudo_df[col].replace(-9999, np.nan)

    # Add species metadata
    for col in ["taxon"]:
        if col in pres_df.columns:
            pseudo_df[col] = pres_df[col].iloc[0]

    # Combine
    all_cols = [c for c in pres_df.columns if c in pseudo_df.columns]
    result = pd.concat([pres_df[all_cols], pseudo_df[all_cols]], ignore_index=True)
    result["reach_idx"] = result["reach_idx"].astype(int)

    print(f"\n    Final dataset: {len(result):,} rows")
    print(f"      Presence:       {int(result['presence'].sum()):,}")
    print(f"      Pseudo-absence: {int((result['presence'] == 0).sum()):,}")
    print(f"      Prevalence:      {result['presence'].mean():.3f}")

    return result


# ══════════════════════════════════════════════════════════════════════════════
# STEP 5: Save outputs
# ══════════════════════════════════════════════════════════════════════════════
def save_outputs(
    cast_df:   pd.DataFrame,
    env_cols:   list[str],
    species_df: pd.DataFrame,
):
    """Save all output files."""
    print(f"\n[Step 5] Saving outputs ...")

    # ── 5a. Master CAST-ready dataset ─────────────────────────────────────────
    meta_cols = ["longitude", "latitude", "taxon", "presence",
                 "reach_idx", "dist_km", "matched"]
    out_cols  = meta_cols + [c for c in env_cols if c in cast_df.columns]
    out_df = cast_df[[c for c in out_cols if c in cast_df.columns]].copy()

    master_path = OUT_DIR / "Lota_lota_CAST_Ready.csv"
    out_df.to_csv(master_path, index=False, float_format="%.6g")
    print(f"    Master CSV : {master_path}")
    print(f"    Shape: {out_df.shape}")

    # ── 5b. Summary statistics ────────────────────────────────────────────────
    summary = {
        "species":          species_df["taxon"].iloc[0],
        "n_presence":      int(out_df["presence"].sum()),
        "n_pseudoabsence": int((out_df["presence"] == 0).sum()),
        "n_matched":        int(out_df["matched"].sum()),
        "n_unmatched":      int((~out_df["matched"]).sum()),
        "prevalence":       round(float(out_df["presence"].mean()), 4),
        "n_env_vars":       len([c for c in env_cols if c in out_df.columns]),
        "max_join_dist_km": MAX_JOIN_DIST_KM,
        "pseudoabs_ratio":  PSEUDOABS_RATIO,
    }
    summary_df = pd.DataFrame([summary])
    summary_path = OUT_DIR / "CAST_Data_Summary.csv"
    summary_df.to_csv(summary_path, index=False)
    print(f"    Summary    : {summary_path}")

    # ── 5c. NA coverage table ─────────────────────────────────────────────────
    if env_cols:
        env_data = out_df[[c for c in env_cols if c in out_df.columns]]
        na_cov = pd.DataFrame({
            "variable":     env_data.columns,
            "n_total":      len(env_data),
            "n_missing":    env_data.isna().sum().values,
            "pct_missing":  (env_data.isna().mean() * 100).round(2).values,
        })
        na_cov_path = OUT_DIR / "Env_NA_Coverage.csv"
        na_cov.to_csv(na_cov_path, index=False)
        print(f"    NA coverage: {na_cov_path}")

    return summary


# ══════════════════════════════════════════════════════════════════════════════
# STEP 6: Quick diagnostic plots
# ══════════════════════════════════════════════════════════════════════════════
def plot_distribution_check(cast_df: pd.DataFrame, species_df: pd.DataFrame):
    """Plot presence/pseudo-absence distribution map."""
    try:
        import matplotlib as mpl
        import matplotlib.pyplot as plt
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature

        mpl.use("Agg")
        mpl.rcParams["font.family"] = "DejaVu Sans"

        pres = cast_df[cast_df["presence"] == 1]
        abse = cast_df[cast_df["presence"] == 0]

        fig = plt.figure(figsize=(14, 7), facecolor="white")
        ax  = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())

        ax.add_feature(cfeature.OCEAN,    facecolor="#EAF6FF", zorder=0)
        ax.add_feature(cfeature.LAND,     facecolor="#F5F5F5", zorder=1)
        ax.add_feature(cfeature.COASTLINE, edgecolor="#4A4A4A", linewidth=0.4, zorder=2)
        ax.add_feature(cfeature.BORDERS,  edgecolor="#AAAAAA", linewidth=0.15, zorder=2)

        # Pseudo-absence (grey)
        if len(abse) > 0:
            ax.scatter(abse["longitude"], abse["latitude"],
                       c="lightgrey", s=2, alpha=0.4,
                       transform=ccrs.PlateCarree(), zorder=4, label="Pseudo-absence")
        # Presence (red)
        ax.scatter(pres["longitude"], pres["latitude"],
                   c="#E74C3C", s=6, alpha=0.75,
                   transform=ccrs.PlateCarree(), zorder=5, label="Presence")

        ax.set_global()
        ax.set_title(f"Lota lota — CAST-ready dataset\n"
                     f"Presence={len(pres):,} | Pseudo-absence={len(abse):,}",
                     fontsize=11)
        ax.legend(loc="lower left", fontsize=8)

        out_path = FIG_DIR / "CAST_Data_Distribution_Map.png"
        fig.savefig(out_path, dpi=150, bbox_inches="tight", facecolor="white")
        plt.close(fig)
        print(f"    Distribution map: {out_path}")

    except Exception as e:
        print(f"    [WARN] Could not generate plot: {e}")


# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════
def main():
    print("=" * 65)
    print("  01_prepare_CAST_data.py  |  Lota lota × RiverATLAS v10")
    print("=" * 65)

    # ── Species data ────────────────────────────────────────────────────────────
    csv_file = DATA_DIR / "Lota_lota.csv"
    species_df = load_species_data(csv_file)

    # ── RiverATLAS global network ───────────────────────────────────────────────
    river_gdf, env_cols = load_global_riveratlas()

    # ── Spatial join ────────────────────────────────────────────────────────────
    joined = spatial_join_points_to_reaches(species_df, river_gdf)

    # ── Attach env vars ────────────────────────────────────────────────────────
    joined_with_env = attach_env_data(joined, river_gdf, env_cols)

    # ── Presence / pseudo-absence ──────────────────────────────────────────────
    cast_df = build_presence_pseudoabsence(joined_with_env, river_gdf, env_cols)

    # ── Save ───────────────────────────────────────────────────────────────────
    summary = save_outputs(cast_df, env_cols, species_df)

    # ── Diagnostic plot ─────────────────────────────────────────────────────────
    plot_distribution_check(cast_df, species_df)

    # ── Final report ────────────────────────────────────────────────────────────
    print("\n" + "=" * 65)
    print("=== CAST Data Preparation Complete ===")
    print("=" * 65)
    print(f"  Species          : {summary['species']}")
    print(f"  Presence points  : {summary['n_presence']:,}")
    print(f"  Pseudo-absences  : {summary['n_pseudoabsence']:,}")
    print(f"  Prevalence       : {summary['prevalence']:.4f}")
    print(f"  Env variables    : {summary['n_env_vars']}")
    print(f"  Max join dist    : {summary['max_join_dist_km']} km")
    print(f"  Output dir       : {OUT_DIR}")
    print("=" * 65)


if __name__ == "__main__":
    main()
