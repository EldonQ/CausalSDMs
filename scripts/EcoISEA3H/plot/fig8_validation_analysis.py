#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig8_validation_analysis.py
================================================================
CAST Paper — Prediction Validation Analysis (Fig 8 Suite)

Three validation components:
  A. Prediction Uncertainty Maps
     - Ensemble SD: standard deviation of HSS across all 6 models
     - If multi-seed data available: within-model SD across seeds
  B. MESS Environmental Extrapolation Detection
     - Identifies grid cells where environmental conditions fall
       outside the range of training data (Elith et al., 2010)
  C. Inter-model Spatial Consistency Metrics
     - Pairwise: Cosine similarity, Warren's I, Pearson's r
     - Reference: GNN-SDM (Wu et al., 2025)

Run:
  cd E:/CausalSDMs && python scripts/EcoISEA3H/plot/fig8_validation_analysis.py

Dependencies:
  pip install cartopy geopandas scipy matplotlib numpy pandas
================================================================
"""

import os, sys, warnings
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.path import Path as MplPath
from scipy.interpolate import griddata
from scipy.spatial.distance import cosine as cosine_dist
from scipy.stats import pearsonr
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import geopandas as gpd
warnings.filterwarnings("ignore")

# ==============================================================================
# ■ USER CONFIGURATION
# ==============================================================================
BASE_DIR      = "E:/CausalSDMs"
os.chdir(BASE_DIR)

PRED_DIR      = "output/case2_eco/spatial_predictions"
TRAIN_DIR     = "outputs/EcoISEA3H/Res9/CAST_ready/species_data_screened"
ENV_GRID_FILE = "outputs/EcoISEA3H/Res9/CAST_ready/China_EnvData_Res9_Screened.csv"
CHINA_SHP     = "plot-function-main/data/china.shp"
DASH_SHP      = "plot-function-main/data/dashline.shp"
OUT_DIR       = "figures/case2_eco/validation"

# Species to analyze (None = all available)
TARGET_SPECIES = [
    "Rhinopithecus_roxellana",
    "Ovis_ammon",
    "Macaca_mulatta",
]

MODELS = ["CAST", "MLP_ATE", "MLP", "RF", "Maxent", "BRT"]
META_COLS = {"HID", "lon", "lat", "species", "sid", "family", "category",
             "presence", "fraction"}

FIG_DPI       = 2400
CHINA_EXTENT  = [73.5, 135, 18, 53.5]
INTERP_RES    = 0.06
INTERP_METHOD = "nearest"
DISPLAY_RES   = 0.02

os.makedirs(OUT_DIR, exist_ok=True)

# ==============================================================================
# ■ Colormap Definitions
# ==============================================================================
def make_hss_cmap():
    colors = ["#0d0887", "#3b049a", "#7201a8", "#a52c60",
              "#d44842", "#ed7953", "#f0f921"]
    return LinearSegmentedColormap.from_list("viridis_habitat", colors, N=256)

def make_uncertainty_cmap():
    """Yellow-Orange-Red for uncertainty (SD)"""
    colors = ["#ffffcc", "#ffeda0", "#fed976", "#feb24c",
              "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026"]
    return LinearSegmentedColormap.from_list("uncertainty", colors, N=256)

def make_mess_cmap():
    """Diverging: brown(extrapolation) - white(0) - green(interpolation)"""
    colors = ["#8c510a", "#bf812d", "#dfc27d", "#f6e8c3",
              "#f5f5f5",
              "#c7eae5", "#80cdc1", "#35978f", "#01665e"]
    return LinearSegmentedColormap.from_list("mess_brg", colors, N=256)

def make_agreement_cmap():
    """Sequential for model agreement count"""
    colors = ["#f7fbff", "#deebf7", "#c6dbef", "#9ecae1",
              "#6baed6", "#4292c6", "#2171b5", "#084594"]
    return LinearSegmentedColormap.from_list("agreement", colors, N=256)

HSS_CMAP       = make_hss_cmap()
UNCERT_CMAP    = make_uncertainty_cmap()
MESS_CMAP      = make_mess_cmap()
AGREE_CMAP     = make_agreement_cmap()

# ==============================================================================
# ■ Cartopy & Interpolation Tools (shared with fig6/fig7)
# ==============================================================================
lon_min, lon_max = CHINA_EXTENT[0], CHINA_EXTENT[1]
lat_min, lat_max = CHINA_EXTENT[2], CHINA_EXTENT[3]
grid_lon = np.arange(lon_min, lon_max + INTERP_RES * 0.5, INTERP_RES)
grid_lat = np.arange(lat_min, lat_max + INTERP_RES * 0.5, INTERP_RES)

if DISPLAY_RES and DISPLAY_RES > INTERP_RES:
    step = max(1, int(round(DISPLAY_RES / INTERP_RES)))
    plot_lon = grid_lon[::step]
    plot_lat = grid_lat[::step]
else:
    step = 1
    plot_lon = grid_lon
    plot_lat = grid_lat

_china_geom_cache = None
_boundary_cache = {}

def get_china_geometry():
    global _china_geom_cache
    effective_buffer = 0.5 * INTERP_RES
    if _china_geom_cache is not None:
        buf_used, geom = _china_geom_cache
        if buf_used == effective_buffer:
            return geom
    if not os.path.exists(CHINA_SHP):
        return None
    gdf = gpd.read_file(CHINA_SHP)
    if gdf is None or len(gdf) == 0:
        return None
    union = gdf.unary_union
    geom = union.buffer(effective_buffer, resolution=2) if effective_buffer else union
    _china_geom_cache = (effective_buffer, geom)
    return geom

def get_boundary_geoms(shp_path):
    global _boundary_cache
    if shp_path not in _boundary_cache:
        if not os.path.exists(shp_path):
            _boundary_cache[shp_path] = []
        else:
            reader = shpreader.Reader(shp_path)
            _boundary_cache[shp_path] = list(reader.geometries())
            reader.close()
    return _boundary_cache[shp_path]

def _points_inside_china(points_xy, china_geom):
    from shapely.geometry import MultiPolygon
    if china_geom is None:
        return np.ones(len(points_xy), dtype=bool)
    polys = list(china_geom.geoms) if isinstance(china_geom, MultiPolygon) else [china_geom]
    inside = np.zeros(len(points_xy), dtype=bool)
    for poly in polys:
        if poly.is_empty or poly.exterior is None:
            continue
        path = MplPath(np.array(poly.exterior.xy).T)
        inside |= path.contains_points(points_xy)
    return inside

def build_china_mask(grid_lon_1d, grid_lat_1d, china_geom):
    nlat, nlon = len(grid_lat_1d), len(grid_lon_1d)
    if china_geom is None:
        return np.ones((nlat, nlon), dtype=bool)
    try:
        center_lon = (grid_lon_1d[:-1] + grid_lon_1d[1:]) * 0.5
        center_lat = (grid_lat_1d[:-1] + grid_lat_1d[1:]) * 0.5
        lon_c, lat_c = np.meshgrid(center_lon, center_lat)
        pts = np.column_stack([lon_c.ravel(), lat_c.ravel()])
        inside_c = _points_inside_china(pts, china_geom).reshape(nlat - 1, nlon - 1)
        mask = np.zeros((nlat, nlon), dtype=bool)
        mask[:nlat - 1, :nlon - 1] = inside_c
        mask[nlat - 1, :] = mask[nlat - 2, :]
        mask[:, nlon - 1] = mask[:, nlon - 2]
        return mask
    except Exception:
        return np.ones((nlat, nlon), dtype=bool)

def add_china_boundary_cached(ax, geoms, **kwargs):
    if geoms:
        ax.add_geometries(geoms, ccrs.PlateCarree(), **kwargs)

def setup_china_axes_subplot(ax, geoms_china, geoms_dash):
    ax.set_extent(CHINA_EXTENT, crs=ccrs.PlateCarree())
    ax.set_aspect("auto")
    if geoms_china:
        add_china_boundary_cached(ax, geoms_china, ec="black", fc="none", linewidth=0.5)
    if geoms_dash:
        add_china_boundary_cached(ax, geoms_dash, ec="black", fc="none", linewidth=0.6)
    ax.patch.set_facecolor("white")
    ax.set_xticks([])
    ax.set_yticks([])
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1)

def compute_one_grid(lons, lats, vals, china_mask):
    if len(lons) < 30:
        return None
    points = np.column_stack([lons, lats])
    glon2d, glat2d = np.meshgrid(grid_lon, grid_lat)
    g = griddata(points, vals, (glon2d, glat2d), method=INTERP_METHOD)
    g_near = griddata(points, vals, (glon2d, glat2d), method="nearest")
    g = np.where(np.isnan(g), g_near, g)
    out = np.array(g, dtype=float, copy=True)
    out[~china_mask] = np.nan
    return out

# Pre-compute
china_geom = get_china_geometry()
china_mask = build_china_mask(grid_lon, grid_lat, china_geom)
geoms_china = get_boundary_geoms(CHINA_SHP)
geoms_dash = get_boundary_geoms(DASH_SHP)

def format_species_name(sp):
    parts = sp.split("_")
    if len(parts) >= 2:
        return f"$\\it{{{parts[0]}}}$  $\\it{{{parts[1]}}}$"
    return sp

def plot_interpolated_map(ax, lon, lat, values, cmap, vmin, vmax,
                          title="", show_cbar=True, cbar_label="",
                          diverging=False, extend="neither"):
    g_out = compute_one_grid(lon, lat, values, china_mask)
    if g_out is None:
        g_out = np.full((len(grid_lat), len(grid_lon)), np.nan)
    plot_vals = g_out[::step, ::step] if step > 1 else g_out

    if diverging and vmin < 0 < vmax:
        norm_i = TwoSlopeNorm(vmin=vmin, vcenter=0, vmax=vmax)
    else:
        norm_i = mcolors.Normalize(vmin=vmin, vmax=vmax)

    cmap_copy = cmap.copy() if hasattr(cmap, "copy") else plt.get_cmap(cmap).copy()
    cmap_copy.set_bad("none")

    mesh = ax.pcolormesh(plot_lon, plot_lat, plot_vals,
                         transform=ccrs.PlateCarree(),
                         cmap=cmap_copy, norm=norm_i,
                         shading="auto", rasterized=True, zorder=1)
    setup_china_axes_subplot(ax, geoms_china, geoms_dash)
    if title:
        ax.set_title(title, fontsize=10, fontweight="bold", color="black", pad=4,
                     path_effects=[pe.withStroke(linewidth=2, foreground="white")])
    if show_cbar:
        cax = ax.inset_axes([0.05, -0.06, 0.9, 0.04])
        cbar = plt.colorbar(mesh, cax=cax, orientation="horizontal", extend=extend)
        cbar.set_label(cbar_label, fontsize=7, color="black", labelpad=1)
        cbar.ax.tick_params(labelsize=6, colors="black", pad=1)
        cbar.outline.set_edgecolor("black")
        cbar.outline.set_linewidth(0.5)
    return mesh

# ==============================================================================
# ■ COMPONENT A: Prediction Uncertainty
# ==============================================================================
def compute_ensemble_stats(pred_df, models=None):
    """Compute ensemble mean, SD, and model agreement across available models."""
    if models is None:
        models = MODELS
    hss_cols = [f"HSS_{m}" for m in models if f"HSS_{m}" in pred_df.columns]
    if len(hss_cols) < 2:
        return None
    hss_matrix = pred_df[hss_cols].values.astype(float)
    ens_mean = np.nanmean(hss_matrix, axis=1)
    ens_sd   = np.nanstd(hss_matrix, axis=1, ddof=1)
    # Model agreement: how many models predict presence (HSS >= 0.5)
    n_agree  = np.nansum(hss_matrix >= 0.5, axis=1)
    # 95% CI width (using t-distribution for small n)
    from scipy.stats import t as t_dist
    n_models = np.sum(~np.isnan(hss_matrix), axis=1)
    t_val = t_dist.ppf(0.975, df=np.maximum(n_models - 1, 1))
    ci_width = 2 * t_val * ens_sd / np.sqrt(np.maximum(n_models, 1))
    return {
        "ensemble_mean": ens_mean,
        "ensemble_sd": ens_sd,
        "ci95_width": ci_width,
        "n_agree": n_agree,
        "n_models": len(hss_cols),
    }

def plot_uncertainty_panel(species_name, pred_df, out_path):
    """Generate 2×2 uncertainty panel: CAST HSS, Ensemble SD, 95% CI Width, Agreement."""
    stats = compute_ensemble_stats(pred_df)
    if stats is None:
        print(f"    ⚠ No ensemble stats for {species_name}")
        return
    lon, lat = pred_df["lon"].values, pred_df["lat"].values

    fig = plt.figure(figsize=(16, 12), facecolor="white")
    fig.suptitle(f"Prediction Uncertainty — {format_species_name(species_name)}",
                 fontsize=14, fontweight="bold", color="black", y=0.97)
    gs = GridSpec(2, 2, figure=fig, hspace=0.30, wspace=0.12,
                  top=0.92, bottom=0.06, left=0.03, right=0.97)

    # Panel 1: CAST HSS (reference)
    ax0 = fig.add_subplot(gs[0, 0], projection=ccrs.PlateCarree())
    hss_cast = pred_df["HSS_CAST"].values if "HSS_CAST" in pred_df.columns \
               else stats["ensemble_mean"]
    vmax = max(np.nanmax(hss_cast), 1.0)
    plot_interpolated_map(ax0, lon, lat, hss_cast, HSS_CMAP, 0, vmax,
                          title="(a) CAST Prediction (HSS)",
                          cbar_label="Habitat Suitability Score")

    # Panel 2: Ensemble SD
    ax1 = fig.add_subplot(gs[0, 1], projection=ccrs.PlateCarree())
    sd_max = np.nanpercentile(stats["ensemble_sd"], 99)
    plot_interpolated_map(ax1, lon, lat, stats["ensemble_sd"],
                          UNCERT_CMAP, 0, max(sd_max, 0.01),
                          title=f"(b) Ensemble SD ({stats['n_models']} models)",
                          cbar_label="Standard Deviation")

    # Panel 3: 95% CI Width
    ax2 = fig.add_subplot(gs[1, 0], projection=ccrs.PlateCarree())
    ci_max = np.nanpercentile(stats["ci95_width"], 99)
    plot_interpolated_map(ax2, lon, lat, stats["ci95_width"],
                          UNCERT_CMAP, 0, max(ci_max, 0.01),
                          title="(c) 95% Confidence Interval Width",
                          cbar_label="CI Width")

    # Panel 4: Model Agreement
    ax3 = fig.add_subplot(gs[1, 1], projection=ccrs.PlateCarree())
    plot_interpolated_map(ax3, lon, lat, stats["n_agree"],
                          AGREE_CMAP, 0, stats["n_models"],
                          title="(d) Model Agreement (HSS ≥ 0.5)",
                          cbar_label="Number of Models Agreeing")

    fig.savefig(out_path, dpi=FIG_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"    ✓ Uncertainty: {os.path.basename(out_path)}")


# ==============================================================================
# ■ COMPONENT B: MESS Environmental Extrapolation
# ==============================================================================
def compute_mess(training_env, prediction_env, env_cols):
    """
    Compute MESS (Multivariate Environmental Similarity Surface).
    Elith et al. (2010) Methods in Ecology and Evolution.

    For each variable i at each prediction location j:
      p = % of training points with value ≤ prediction value
      If p=0:   f = (val - min_train) / range × 100  (negative if below min)
      If 0<p≤50: f = 2 × p
      If 50≤p<100: f = 2 × (100 - p)
      If p=100: f = (max_train - val) / range × 100  (negative if above max)

    MESS_j = min(f_ij) across all variables
    MESS < 0 → extrapolation (at least one variable outside training range)
    """
    n_pred = len(prediction_env)
    n_vars = len(env_cols)
    similarity = np.full((n_pred, n_vars), np.nan)

    for i, col in enumerate(env_cols):
        train_vals = training_env[col].dropna().values
        pred_vals  = prediction_env[col].values

        if len(train_vals) < 5:
            similarity[:, i] = 100.0
            continue

        min_val   = np.min(train_vals)
        max_val   = np.max(train_vals)
        range_val = max_val - min_val

        if range_val < 1e-10:
            # Constant variable → anything different is extrapolation
            similarity[:, i] = np.where(
                np.abs(pred_vals - min_val) < 1e-10, 100.0, -100.0
            )
            continue

        # Vectorized percentile: % of training values ≤ each pred value
        sorted_train = np.sort(train_vals)
        p = np.searchsorted(sorted_train, pred_vals, side="right") \
            / len(train_vals) * 100.0

        f = np.zeros(n_pred)

        # Case 1: p == 0 → prediction below all training values
        mask0 = p == 0
        f[mask0] = (pred_vals[mask0] - min_val) / range_val * 100

        # Case 2: 0 < p ≤ 50
        mask_low = (p > 0) & (p <= 50)
        f[mask_low] = 2 * p[mask_low]

        # Case 3: 50 < p < 100
        mask_high = (p > 50) & (p < 100)
        f[mask_high] = 2 * (100 - p[mask_high])

        # Case 4: p == 100 → prediction above all training values
        mask100 = p == 100
        f[mask100] = (max_val - pred_vals[mask100]) / range_val * 100

        similarity[:, i] = f

    # MESS = minimum similarity across all variables
    mess_vals = np.nanmin(similarity, axis=1)
    # Most dissimilar variable (the one causing the lowest similarity)
    most_dissimilar_idx = np.nanargmin(similarity, axis=1)
    most_dissimilar_var = [env_cols[idx] for idx in most_dissimilar_idx]

    return mess_vals, most_dissimilar_var, similarity


def plot_mess_panel(species_name, pred_df, mess_vals, most_dissimilar_var,
                    env_cols, out_path):
    """Generate MESS panel: MESS map + Extrapolation binary mask."""
    lon, lat = pred_df["lon"].values, pred_df["lat"].values

    fig = plt.figure(figsize=(16, 6.5), facecolor="white")
    fig.suptitle(f"MESS Extrapolation Analysis — {format_species_name(species_name)}",
                 fontsize=14, fontweight="bold", color="black", y=0.97)
    gs = GridSpec(1, 2, figure=fig, wspace=0.12,
                  top=0.90, bottom=0.08, left=0.03, right=0.97)

    # Panel 1: MESS values (diverging)
    ax0 = fig.add_subplot(gs[0, 0], projection=ccrs.PlateCarree())
    mess_min = np.nanpercentile(mess_vals, 1)
    mess_max = np.nanpercentile(mess_vals, 99)
    abs_lim = max(abs(mess_min), abs(mess_max), 1.0)
    plot_interpolated_map(ax0, lon, lat, mess_vals,
                          MESS_CMAP, -abs_lim, abs_lim,
                          title="(a) MESS Values",
                          cbar_label="MESS (< 0 = extrapolation)",
                          diverging=True, extend="both")

    # Panel 2: Extrapolation mask (binary)
    ax1 = fig.add_subplot(gs[1, 0] if gs.nrows > 1 else gs[0, 1],
                          projection=ccrs.PlateCarree())
    extrap_mask = (mess_vals < 0).astype(float)
    binary_cmap = LinearSegmentedColormap.from_list(
        "extrap", ["#2ca02c", "#d62728"], N=2)
    plot_interpolated_map(ax1, lon, lat, extrap_mask,
                          binary_cmap, 0, 1,
                          title="(b) Extrapolation Risk (MESS < 0)",
                          cbar_label="Green=Interpolation  Red=Extrapolation")

    # Stats annotation
    n_total = len(mess_vals)
    n_extrap = np.sum(mess_vals < 0)
    pct_extrap = n_extrap / n_total * 100
    stats_text = (f"Total cells: {n_total:,}\n"
                  f"Extrapolation: {n_extrap:,} ({pct_extrap:.1f}%)\n"
                  f"Interpolation: {n_total - n_extrap:,} ({100 - pct_extrap:.1f}%)")
    fig.text(0.50, 0.02, stats_text, ha="center", va="bottom",
             fontsize=9, color="black", fontfamily="monospace",
             bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow",
                       edgecolor="gray", alpha=0.9))

    fig.savefig(out_path, dpi=FIG_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"    ✓ MESS: {os.path.basename(out_path)}")


# ==============================================================================
# ■ COMPONENT C: Inter-model Spatial Consistency
# ==============================================================================
def warrens_i(p1, p2):
    """
    Warren's I niche overlap metric (Warren et al. 2008).
    I = Σ √(p1_norm × p2_norm) where p_norm sums to 1.
    Range: [0, 1], 1 = perfect overlap.
    """
    valid = ~(np.isnan(p1) | np.isnan(p2)) & (p1 >= 0) & (p2 >= 0)
    a, b = p1[valid], p2[valid]
    sa, sb = a.sum(), b.sum()
    if sa < 1e-10 or sb < 1e-10:
        return np.nan
    a_norm = a / sa
    b_norm = b / sb
    return float(np.sum(np.sqrt(a_norm * b_norm)))

def schoeners_d(p1, p2):
    """
    Schoener's D niche overlap (Schoener 1968).
    D = 1 - 0.5 × Σ |p1_norm - p2_norm|
    Range: [0, 1], 1 = identical.
    """
    valid = ~(np.isnan(p1) | np.isnan(p2)) & (p1 >= 0) & (p2 >= 0)
    a, b = p1[valid], p2[valid]
    sa, sb = a.sum(), b.sum()
    if sa < 1e-10 or sb < 1e-10:
        return np.nan
    return float(1.0 - 0.5 * np.sum(np.abs(a / sa - b / sb)))

def cosine_similarity(p1, p2):
    """Cosine similarity between two suitability vectors."""
    valid = ~(np.isnan(p1) | np.isnan(p2))
    a, b = p1[valid], p2[valid]
    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    if norm_a < 1e-10 or norm_b < 1e-10:
        return np.nan
    return float(np.dot(a, b) / (norm_a * norm_b))

def pearson_corr(p1, p2):
    """Pearson correlation between two suitability vectors."""
    valid = ~(np.isnan(p1) | np.isnan(p2))
    a, b = p1[valid], p2[valid]
    if len(a) < 10:
        return np.nan, np.nan
    r, p = pearsonr(a, b)
    return float(r), float(p)

def compute_pairwise_metrics(pred_df, models=None):
    """Compute all pairwise consistency metrics."""
    if models is None:
        models = MODELS
    avail = [m for m in models if f"HSS_{m}" in pred_df.columns]
    n = len(avail)
    results = []

    for i in range(n):
        for j in range(i + 1, n):
            p1 = pred_df[f"HSS_{avail[i]}"].values.astype(float)
            p2 = pred_df[f"HSS_{avail[j]}"].values.astype(float)

            cos = cosine_similarity(p1, p2)
            wi  = warrens_i(p1, p2)
            sd  = schoeners_d(p1, p2)
            r, pval = pearson_corr(p1, p2)

            results.append({
                "model_a": avail[i], "model_b": avail[j],
                "cosine_sim": cos, "warrens_I": wi,
                "schoeners_D": sd, "pearson_r": r, "pearson_p": pval,
            })
    return pd.DataFrame(results), avail

def plot_consistency_heatmaps(species_name, metrics_df, model_list, out_path):
    """Generate heatmap matrices for consistency metrics."""
    n = len(model_list)
    metric_names = [
        ("cosine_sim", "Cosine Similarity"),
        ("warrens_I", "Warren's I"),
        ("pearson_r", "Pearson's r"),
    ]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.5), facecolor="white")
    fig.suptitle(f"Inter-model Spatial Consistency — {format_species_name(species_name)}",
                 fontsize=14, fontweight="bold", color="black", y=1.02)

    for ax_idx, (metric_key, metric_label) in enumerate(metric_names):
        mat = np.ones((n, n))
        for _, row in metrics_df.iterrows():
            i = model_list.index(row["model_a"])
            j = model_list.index(row["model_b"])
            val = row[metric_key]
            if np.isnan(val):
                val = 0
            mat[i, j] = val
            mat[j, i] = val

        ax = axes[ax_idx]
        im = ax.imshow(mat, cmap="YlGnBu", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xticklabels(model_list, fontsize=8, rotation=45, ha="right")
        ax.set_yticklabels(model_list, fontsize=8)
        ax.set_title(metric_label, fontsize=11, fontweight="bold", pad=8)

        # Annotate values
        for i in range(n):
            for j in range(n):
                text_color = "white" if mat[i, j] > 0.7 else "black"
                ax.text(j, i, f"{mat[i, j]:.3f}",
                        ha="center", va="center", fontsize=8,
                        color=text_color, fontweight="bold")
        plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)

    fig.tight_layout()
    fig.savefig(out_path, dpi=FIG_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"    ✓ Consistency: {os.path.basename(out_path)}")


# ==============================================================================
# ■ COMPONENT D: Prediction + MESS Credibility Zone Overlay
# ==============================================================================
def make_zone_cmap():
    """3-zone: red (extrapolation) → amber (buffer) → green (reliable)"""
    colors = ["#d62728", "#ff7f0e", "#2ca02c"]
    return LinearSegmentedColormap.from_list("zone", colors, N=3)


def plot_credibility_map(species_name, pred_df, mess_vals, env_grid_df,
                        train_df, out_path):
    """
    Generate CAST prediction heatmap with MESS credibility zone overlay.

    Three zones:
      Zone 1 (green):  Within known distribution range
      Zone 2 (amber):  Buffer/adjacent area, MESS ≥ 0 (environmental interpolation)
      Zone 3 (red):    MESS < 0 (environmental extrapolation — lower confidence)
    """
    lon = pred_df["lon"].values
    lat = pred_df["lat"].values
    hss_cast = pred_df["HSS_CAST"].values if "HSS_CAST" in pred_df.columns \
               else pred_df["HSS_MLP"].values

    # ── Identify known presence range ──
    if train_df is not None and "presence" in train_df.columns:
        pres = train_df[train_df["presence"] == 1]
        pres_lons = pres["lon"].values
        pres_lats = pres["lat"].values
    else:
        pres_lons = np.array([])
        pres_lats = np.array([])

    # ── Build zone classification for each grid cell ──
    # Zone values: 1=known range, 2=buffer(MESS≥0), 3=extrapolation(MESS<0)
    # Use env_grid coordinates (same order as mess_vals)
    grid_lons_env = env_grid_df["lon"].values
    grid_lats_env = env_grid_df["lat"].values

    # Classify each grid cell
    zone_vals = np.full(len(mess_vals), 3.0)  # default: extrapolation
    zone_vals[mess_vals >= 0] = 2.0            # interpolation (buffer)

    # Mark cells within known presence range (convex hull buffer)
    if len(pres_lons) >= 3:
        from scipy.spatial import ConvexHull
        try:
            pts_pres = np.column_stack([pres_lons, pres_lats])
            hull = ConvexHull(pts_pres)
            hull_path = MplPath(pts_pres[hull.vertices])
            grid_pts = np.column_stack([grid_lons_env, grid_lats_env])
            # Zone 1: within convex hull of presence points (with ~50km buffer ≈ 0.5°)
            from shapely.geometry import Polygon as ShapelyPolygon
            hull_poly = ShapelyPolygon(pts_pres[hull.vertices])
            hull_buffered = hull_poly.buffer(0.5)  # ~50km at mid-latitudes
            hull_path_buf = MplPath(
                np.array(hull_buffered.exterior.xy).T
            )
            in_range = hull_path_buf.contains_points(grid_pts)
            zone_vals[in_range & (mess_vals >= 0)] = 1.0
        except Exception:
            pass

    # ── Create figure: 1×2 (HSS + HSS with zones) ──
    fig = plt.figure(figsize=(17, 7.5), facecolor="white")
    fig.suptitle(
        f"Prediction Credibility Zones — {format_species_name(species_name)}",
        fontsize=14, fontweight="bold", color="black", y=0.97,
    )
    gs = GridSpec(1, 2, figure=fig, wspace=0.10,
                  top=0.90, bottom=0.10, left=0.03, right=0.97)

    # ── Panel 1: CAST HSS prediction (clean, like fig7) ──
    ax0 = fig.add_subplot(gs[0, 0], projection=ccrs.PlateCarree())
    vmax = max(np.nanmax(hss_cast), 1.0)
    plot_interpolated_map(ax0, lon, lat, hss_cast, HSS_CMAP, 0, vmax,
                          title="(a) CAST Prediction (HSS)",
                          cbar_label="Habitat Suitability Score")

    # ── Panel 2: HSS + MESS credibility contour overlay ──
    ax1 = fig.add_subplot(gs[0, 1], projection=ccrs.PlateCarree())

    # Base layer: HSS heatmap (same as panel 1)
    hss_grid = compute_one_grid(lon, lat, hss_cast, china_mask)
    if hss_grid is None:
        hss_grid = np.full((len(grid_lat), len(grid_lon)), np.nan)
    plot_hss = hss_grid[::step, ::step] if step > 1 else hss_grid
    hss_cmap = HSS_CMAP.copy()
    hss_cmap.set_bad("none")
    norm_hss = mcolors.Normalize(vmin=0, vmax=vmax)
    ax1.pcolormesh(plot_lon, plot_lat, plot_hss,
                   transform=ccrs.PlateCarree(),
                   cmap=hss_cmap, norm=norm_hss,
                   shading="auto", rasterized=True, zorder=1)

    # Overlay: MESS=0 contour line (interpolation/extrapolation boundary)
    mess_grid = compute_one_grid(grid_lons_env, grid_lats_env, mess_vals, china_mask)
    if mess_grid is not None:
        plot_mess = mess_grid[::step, ::step] if step > 1 else mess_grid
        lon2d, lat2d = np.meshgrid(plot_lon, plot_lat)
        # MESS = 0 contour (thick red line = extrapolation boundary)
        c0 = ax1.contour(lon2d, lat2d, plot_mess,
                         levels=[0], colors=["#d62728"],
                         linewidths=[1.8], linestyles=["solid"],
                         transform=ccrs.PlateCarree(), zorder=4)
        ax1.clabel(c0, fmt="MESS=0", fontsize=6, colors=["#d62728"])

    # Overlay: Zone 1 boundary (known range convex hull)
    if len(pres_lons) >= 3:
        try:
            hull_poly_plot = ShapelyPolygon(pts_pres[hull.vertices])
            hull_buf_plot = hull_poly_plot.buffer(0.5)
            from matplotlib.patches import Polygon as MplPolygon
            from shapely.geometry import MultiPolygon as ShapelyMultiPolygon
            polys = [hull_buf_plot] if not isinstance(hull_buf_plot, ShapelyMultiPolygon) \
                    else list(hull_buf_plot.geoms)
            for poly in polys:
                coords = np.array(poly.exterior.xy).T
                ax1.plot(coords[:, 0], coords[:, 1],
                         color="#2ca02c", linewidth=1.5, linestyle="--",
                         transform=ccrs.PlateCarree(), zorder=5,
                         label="Known range (~50km buffer)")
        except Exception:
            pass

    setup_china_axes_subplot(ax1, geoms_china, geoms_dash)
    ax1.set_title("(b) HSS + Credibility Zones", fontsize=10,
                  fontweight="bold", color="black", pad=4,
                  path_effects=[pe.withStroke(linewidth=2, foreground="white")])

    # HSS colorbar
    cax1 = ax1.inset_axes([0.05, -0.06, 0.9, 0.04])
    sm = plt.cm.ScalarMappable(cmap=hss_cmap, norm=norm_hss)
    cbar1 = plt.colorbar(sm, cax=cax1, orientation="horizontal")
    cbar1.set_label("HSS", fontsize=7, labelpad=1)
    cbar1.ax.tick_params(labelsize=6, pad=1)
    cbar1.outline.set_edgecolor("black")
    cbar1.outline.set_linewidth(0.5)

    # Legend for contour lines
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color="#2ca02c", lw=1.5, ls="--",
               label="Zone 1: Known range (high confidence)"),
        Line2D([0], [0], color="#d62728", lw=1.8, ls="-",
               label="MESS=0: Interpolation/Extrapolation boundary"),
    ]
    ax1.legend(handles=legend_elements, loc="lower left",
              fontsize=7, framealpha=0.9, edgecolor="gray",
              fancybox=True, bbox_to_anchor=(0.01, 0.01))

    # Stats annotation
    n_total = len(zone_vals)
    n_z1 = int(np.sum(zone_vals == 1))
    n_z2 = int(np.sum(zone_vals == 2))
    n_z3 = int(np.sum(zone_vals == 3))
    stats_text = (
        f"Zone 1 (known range): {n_z1:,} ({n_z1/n_total*100:.1f}%)\n"
        f"Zone 2 (interpolation): {n_z2:,} ({n_z2/n_total*100:.1f}%)\n"
        f"Zone 3 (extrapolation): {n_z3:,} ({n_z3/n_total*100:.1f}%)"
    )
    fig.text(0.50, 0.02, stats_text, ha="center", va="bottom",
             fontsize=8, color="black", fontfamily="monospace",
             bbox=dict(boxstyle="round,pad=0.4", facecolor="lightyellow",
                       edgecolor="gray", alpha=0.9))

    fig.savefig(out_path, dpi=FIG_DPI, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(f"    ✓ Credibility: {os.path.basename(out_path)}")


# ==============================================================================
# ■ MAIN
# ==============================================================================
def main():
    print("=" * 60)
    print("  CAST Fig 8 — Prediction Validation Analysis")
    print("=" * 60)

    # ── Load environmental grid data (for MESS) ──
    print("\n  Loading environmental grid data...")
    env_grid = pd.read_csv(ENV_GRID_FILE)
    env_cols = [c for c in env_grid.columns if c not in META_COLS]
    print(f"  Env grid: {len(env_grid):,} cells, {len(env_cols)} variables")
    print(f"  Variables: {env_cols}")

    # ── Discover species ──
    pred_files = [f for f in os.listdir(PRED_DIR) if f.startswith("pred_") and f.endswith(".csv")]
    available_species = [f.replace("pred_", "").replace(".csv", "") for f in pred_files]
    if TARGET_SPECIES:
        species_list = [sp for sp in TARGET_SPECIES if sp in available_species]
    else:
        species_list = available_species
    print(f"  Species to analyze: {len(species_list)}")

    # ── Summary collectors ──
    all_consistency = []
    all_mess_summary = []

    for sp_idx, sp in enumerate(species_list):
        print(f"\n  ─── [{sp_idx + 1}/{len(species_list)}] {sp} ───")

        # Load prediction data
        pred_path = os.path.join(PRED_DIR, f"pred_{sp}.csv")
        pred_df = pd.read_csv(pred_path)
        print(f"    Predictions: {len(pred_df):,} cells")

        # Load training data
        train_path = os.path.join(TRAIN_DIR, f"CAST_{sp}_Res9_screened.csv")
        if not os.path.exists(train_path):
            print(f"    ⚠ Training data not found: {train_path}")
            train_df = None
        else:
            train_df = pd.read_csv(train_path)
            print(f"    Training data: {len(train_df):,} cells")

        # ──────────────────────────────────────────────────────────
        # A. Prediction Uncertainty
        # ──────────────────────────────────────────────────────────
        print("    [A] Computing prediction uncertainty...")
        plot_uncertainty_panel(
            sp, pred_df,
            os.path.join(OUT_DIR, f"fig8_uncertainty_{sp}.png"),
        )

        # ──────────────────────────────────────────────────────────
        # B. MESS Extrapolation Detection
        # ──────────────────────────────────────────────────────────
        if train_df is not None:
            print("    [B] Computing MESS extrapolation analysis...")
            # Use presence locations only as the reference
            train_env = train_df[train_df["presence"] == 1][env_cols].copy()
            if len(train_env) < 10:
                # Fall back to all training data
                train_env = train_df[env_cols].copy()

            pred_env = env_grid[env_cols].copy()
            mess_vals, most_dissim, _ = compute_mess(train_env, pred_env, env_cols)

            # Save MESS to CSV
            mess_out = pd.DataFrame({
                "HID": env_grid["HID"], "lon": env_grid["lon"], "lat": env_grid["lat"],
                "MESS": mess_vals, "most_dissimilar_var": most_dissim,
            })
            mess_csv_path = os.path.join(OUT_DIR, f"mess_{sp}.csv")
            mess_out.to_csv(mess_csv_path, index=False)

            plot_mess_panel(
                sp, pred_df, mess_vals, most_dissim, env_cols,
                os.path.join(OUT_DIR, f"fig8_mess_{sp}.png"),
            )

            # Summary stats
            n_extrap = int(np.sum(mess_vals < 0))
            all_mess_summary.append({
                "species": sp,
                "n_total": len(mess_vals),
                "n_extrapolation": n_extrap,
                "pct_extrapolation": n_extrap / len(mess_vals) * 100,
                "mess_min": float(np.nanmin(mess_vals)),
                "mess_median": float(np.nanmedian(mess_vals)),
                "mess_mean": float(np.nanmean(mess_vals)),
            })
            # ──────────────────────────────────────────────────────────
            # D. Credibility Zone Map (HSS + MESS overlay)
            # ──────────────────────────────────────────────────────────
            print("    [D] Generating credibility zone map...")
            plot_credibility_map(
                sp, pred_df, mess_vals, env_grid, train_df,
                os.path.join(OUT_DIR, f"fig8_credibility_{sp}.png"),
            )
        else:
            print("    [B] Skipped MESS (no training data)")
            mess_vals = None

        # ──────────────────────────────────────────────────────────
        # C. Inter-model Consistency
        # ──────────────────────────────────────────────────────────
        print("    [C] Computing inter-model consistency metrics...")
        metrics_df, model_list = compute_pairwise_metrics(pred_df)
        if not metrics_df.empty:
            metrics_df["species"] = sp
            all_consistency.append(metrics_df)
            plot_consistency_heatmaps(
                sp, metrics_df, model_list,
                os.path.join(OUT_DIR, f"fig8_consistency_{sp}.png"),
            )

    # ── Save summary tables ──
    if all_consistency:
        cons_all = pd.concat(all_consistency, ignore_index=True)
        cons_path = os.path.join(OUT_DIR, "fig8_consistency_summary.csv")
        cons_all.to_csv(cons_path, index=False)
        print(f"\n  ✓ Consistency summary: {cons_path}")

        # Print summary
        print("\n  ── Inter-model Consistency Summary ──")
        for sp in species_list:
            sp_df = cons_all[cons_all["species"] == sp]
            if sp_df.empty:
                continue
            print(f"  {sp}:")
            print(f"    Cosine sim:  mean={sp_df['cosine_sim'].mean():.4f} "
                  f"(range: {sp_df['cosine_sim'].min():.4f}–{sp_df['cosine_sim'].max():.4f})")
            print(f"    Warren's I:  mean={sp_df['warrens_I'].mean():.4f} "
                  f"(range: {sp_df['warrens_I'].min():.4f}–{sp_df['warrens_I'].max():.4f})")
            print(f"    Pearson r:   mean={sp_df['pearson_r'].mean():.4f} "
                  f"(range: {sp_df['pearson_r'].min():.4f}–{sp_df['pearson_r'].max():.4f})")

    if all_mess_summary:
        mess_sum = pd.DataFrame(all_mess_summary)
        mess_path = os.path.join(OUT_DIR, "fig8_mess_summary.csv")
        mess_sum.to_csv(mess_path, index=False)
        print(f"\n  ✓ MESS summary: {mess_path}")
        print("\n  ── MESS Extrapolation Summary ──")
        for _, row in mess_sum.iterrows():
            print(f"  {row['species']}: {row['pct_extrapolation']:.1f}% extrapolation "
                  f"(MESS range: {row['mess_min']:.1f} to {row['mess_median']:.1f} median)")

    print(f"\n{'=' * 60}")
    print(f"  All validation figures saved to: {OUT_DIR}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
