#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
fig7_spatial_prediction_maps.py
================================================================
CAST Paper — Species Distribution Prediction Maps (Fig 7)

Three-layer progressive design:
  Layer 1: HSS prediction maps for selected models
  Layer 2: Contrastive difference maps (CAST vs baseline)
  Layer 3: CATE overlay (optional, from existing fig6 data)

KEY FEATURE: Fully configurable — choose any species, any models to compare.
This version uses cartopy + griddata interpolation to match fig6 styling.
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.patheffects as pe
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.path import Path as MplPath
from scipy.interpolate import griddata
import cartopy.crs as ccrs
import cartopy.io.shapereader as shpreader
import geopandas as gpd
import warnings
warnings.filterwarnings("ignore")

# ==============================================================================
# ■ USER CONFIGURATION
# ==============================================================================
PRED_DIR  = "E:/CausalSDMs/output/case2_eco/spatial_predictions"
CATE_FILE = "E:/CausalSDMs/output/case2_eco/all_spatial_cate_v3.csv"
OUT_DIR   = "E:/CausalSDMs/figures/case2_eco"
CHINA_SHP  = "plot-function-main/data/china.shp"
DASH_SHP   = "plot-function-main/data/dashline.shp"

SPECIES_LIST = [
    "Rhinopithecus_roxellana",
    "Ovis_ammon",
    "Macaca_mulatta"
]
MODELS_TO_COMPARE = ["CAST", "MLP", "RF", "Maxent", "BRT"]
LAYOUT_MODE = "full"

FIG_DPI = 300
CHINA_EXTENT  = [73.5, 135, 18, 53.5]
INTERP_RES    = 0.06
INTERP_METHOD = "nearest"
DISPLAY_RES   = 0.02
HEX_SIZE = 1.8 

# Grid setup for interpolation
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

def make_hss_cmap():
    colors = ["#0d0887", "#3b049a", "#7201a8", "#a52c60", "#d44842", "#ed7953", "#f0f921"]
    return LinearSegmentedColormap.from_list("viridis_habitat", colors, N=256)

def make_diff_cmap():
    colors = ["#2166ac", "#4393c3", "#92c5de", "#d1e5f0", "#f7f7f7", "#fddbc7", "#f4a582", "#d6604d", "#b2182b"]
    return LinearSegmentedColormap.from_list("diff_rb", colors, N=256)

HSS_CMAP  = make_hss_cmap()
DIFF_CMAP = make_diff_cmap()

# ==============================================================================
# ■ Cartopy & Interpolation Tools (from fig6)
# ==============================================================================
_china_geom_cache = None
_boundary_cache = {}

def get_china_geometry():
    global _china_geom_cache
    effective_buffer = 0.5 * INTERP_RES
    if _china_geom_cache is not None:
        buf_used, geom = _china_geom_cache
        if buf_used == effective_buffer:
            return geom
    if not os.path.exists(CHINA_SHP): return None
    gdf = gpd.read_file(CHINA_SHP)
    if gdf is None or len(gdf) == 0: return None
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
    if china_geom is None: return np.ones(len(points_xy), dtype=bool)
    if isinstance(china_geom, MultiPolygon): polys = list(china_geom.geoms)
    else: polys = [china_geom]
    inside = np.zeros(len(points_xy), dtype=bool)
    for poly in polys:
        if poly.is_empty or poly.exterior is None: continue
        path = MplPath(np.array(poly.exterior.xy).T)
        inside |= path.contains_points(points_xy)
    return inside

def build_china_mask(grid_lon_1d, grid_lat_1d, china_geom):
    nlat, nlon = len(grid_lat_1d), len(grid_lon_1d)
    if china_geom is None: return np.ones((nlat, nlon), dtype=bool)
    try:
        center_lon = (grid_lon_1d[:-1] + grid_lon_1d[1:]) * 0.5
        center_lat = (grid_lat_1d[:-1] + grid_lat_1d[1:]) * 0.5
        lon_c, lat_c = np.meshgrid(center_lon, center_lat)
        points = np.column_stack([lon_c.ravel(), lat_c.ravel()])
        inside_c = _points_inside_china(points, china_geom).reshape(nlat - 1, nlon - 1)
        mask = np.zeros((nlat, nlon), dtype=bool)
        mask[: nlat - 1, : nlon - 1] = inside_c
        mask[nlat - 1, :] = mask[nlat - 2, :]
        mask[:, nlon - 1] = mask[:, nlon - 2]
        return mask
    except Exception:
        return np.ones((nlat, nlon), dtype=bool)

def add_china_boundary_cached(ax, geoms, **kwargs):
    if not geoms: return
    ax.add_geometries(geoms, ccrs.PlateCarree(), **kwargs)

def setup_china_axes_subplot(ax, geoms_china, geoms_dash):
    ax.set_extent(CHINA_EXTENT, crs=ccrs.PlateCarree())
    ax.set_aspect('auto')
    if geoms_china: add_china_boundary_cached(ax, geoms_china, ec="black", fc="none", linewidth=0.5)
    if geoms_dash: add_china_boundary_cached(ax, geoms_dash, ec="black", fc="none", linewidth=0.6)
    ax.patch.set_facecolor("none")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(True)
    ax.patch.set_edgecolor("black")
    ax.patch.set_linewidth(1)
    ax.patch.set_facecolor("white")
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("black")
        spine.set_linewidth(1)

def compute_one_grid(lons, lats, cate_vals, china_mask):
    if len(lons) < 30: return None
    points = np.column_stack([lons, lats])
    grid_lon_2d, grid_lat_2d = np.meshgrid(grid_lon, grid_lat)
    grid_cate = griddata(points, cate_vals, (grid_lon_2d, grid_lat_2d), method=INTERP_METHOD)
    grid_cate_nearest = griddata(points, cate_vals, (grid_lon_2d, grid_lat_2d), method="nearest")
    grid_cate = np.where(np.isnan(grid_cate), grid_cate_nearest, grid_cate)
    out = np.array(grid_cate, dtype=float, copy=True)
    out[~china_mask] = np.nan
    return out

china_geom = get_china_geometry()
china_mask = build_china_mask(grid_lon, grid_lat, china_geom)
geoms_china = get_boundary_geoms(CHINA_SHP)
geoms_dash = get_boundary_geoms(DASH_SHP)

# ==============================================================================
def load_species_predictions(species_name):
    f = os.path.join(PRED_DIR, f"pred_{species_name}.csv")
    if not os.path.exists(f): raise FileNotFoundError(f)
    return pd.read_csv(f)

def load_cate_data(species_name):
    if not os.path.exists(CATE_FILE): return None
    df = pd.read_csv(CATE_FILE)
    sp = df[df["species"] == species_name.replace("_", " ")]
    if sp.empty: sp = df[df["species"] == species_name]
    return sp if not sp.empty else None

def format_species_name(sp):
    parts = sp.split("_")
    if len(parts) >= 2: return f"$\\it{{{parts[0]}}}$  $\\it{{{parts[1]}}}$"
    return sp

def plot_interpolated_map(ax, lon, lat, values, cmap, vmin, vmax, title="", show_cbar=True, cbar_label=""):
    g_out = compute_one_grid(lon, lat, values, china_mask)
    if g_out is None: g_out = np.full((len(grid_lat), len(grid_lon)), np.nan)
    
    plot_vals = g_out[::step, ::step] if step > 1 else g_out
    norm_i = mcolors.Normalize(vmin=vmin, vmax=vmax)
    
    cmap_copy = plt.get_cmap(cmap).copy() if isinstance(cmap, str) else cmap.copy()
    cmap_copy.set_bad("none")
    
    mesh = ax.pcolormesh(plot_lon, plot_lat, plot_vals, transform=ccrs.PlateCarree(),
                         cmap=cmap_copy, norm=norm_i, shading="auto", rasterized=True, zorder=1)
    
    setup_china_axes_subplot(ax, geoms_china, geoms_dash)
    
    if title:
        ax.set_title(title, fontsize=9, fontweight='bold', color='black', pad=4,
                     path_effects=[pe.withStroke(linewidth=2, foreground='white')])
    if show_cbar:
        # Instead of shrinking map layout, use an inset cbar below
        cax = ax.inset_axes([0.05, -0.05, 0.9, 0.04])
        cbar = plt.colorbar(mesh, cax=cax, orientation='horizontal')
        cbar.set_label(cbar_label, fontsize=7, color='black', labelpad=1)
        cbar.ax.tick_params(labelsize=6, colors='black', pad=1)
        cbar.outline.set_edgecolor('black')
        cbar.outline.set_linewidth(0.5)

def plot_species_full(species_name, df, cate_df=None):
    primary = MODELS_TO_COMPARE[0]
    baselines = MODELS_TO_COMPARE[1:]
    avail = [b for b in baselines if f"HSS_{b}" in df.columns]
    if not avail: return None
    
    n_rows = len(avail)
    n_cols = 4 if (cate_df is not None and not cate_df.empty) else 3
    
    fig = plt.figure(figsize=(n_cols * 4.2, n_rows * 3.2 + 0.8), facecolor='white')
    fig.suptitle(f"Species Distribution Predictions — {format_species_name(species_name)}",
                 fontsize=14, fontweight='bold', color='black', y=0.98)
    gs = GridSpec(n_rows, n_cols, figure=fig, hspace=0.35, wspace=0.15, top=0.92, bottom=0.08, left=0.03, right=0.97)
    
    lon, lat = df["lon"].values, df["lat"].values
    hss_primary = df[f"HSS_{primary}"].values
    hss_vmax = max(np.nanmax(hss_primary), 1.0)
    
    for row_i, baseline in enumerate(avail):
        hss_b = df[f"HSS_{baseline}"].values
        d_val = hss_primary - hss_b
        
        ax0 = fig.add_subplot(gs[row_i, 0], projection=ccrs.PlateCarree())
        plot_interpolated_map(ax0, lon, lat, hss_primary, HSS_CMAP, 0, hss_vmax,
                              title=f"{primary}" if row_i==0 else "", show_cbar=(row_i==n_rows-1), cbar_label="HSS")
        ax0.text(-0.08, 0.5, f"vs {baseline}", transform=ax0.transAxes, fontsize=8, color='black',
                 rotation=90, va='center', ha='center', fontweight='bold')
                 
        ax1 = fig.add_subplot(gs[row_i, 1], projection=ccrs.PlateCarree())
        plot_interpolated_map(ax1, lon, lat, hss_b, HSS_CMAP, 0, hss_vmax,
                              title=baseline if row_i==0 else "", show_cbar=(row_i==n_rows-1), cbar_label="HSS")
                              
        ax2 = fig.add_subplot(gs[row_i, 2], projection=ccrs.PlateCarree())
        d_abs = max(abs(np.nanmin(d_val)), abs(np.nanmax(d_val)), 0.01)
        plot_interpolated_map(ax2, lon, lat, d_val, DIFF_CMAP, -d_abs, d_abs,
                              title=(f"{primary} − {baseline}" if row_i==0 else ""), show_cbar=(row_i==n_rows-1), cbar_label="ΔHSS")
                              
        if n_cols == 4:
            ax3 = fig.add_subplot(gs[row_i, 3], projection=ccrs.PlateCarree())
            top_var = cate_df["variable"].value_counts().index[0]
            sub = cate_df[cate_df["variable"] == top_var]
            if row_i == 0:
                c_abs = max(abs(sub["cate"].min()), abs(sub["cate"].max()), 0.01)
                plot_interpolated_map(ax3, sub["lon"].values, sub["lat"].values, sub["cate"].values, DIFF_CMAP, -c_abs, c_abs,
                                      title=f"CATE: {top_var}", show_cbar=True, cbar_label="CATE")
            else:
                ax3.set_visible(False)
    return fig

def plot_multi_species_comparison(species_list, model_a, model_b):
    n_sp = len(species_list)
    fig = plt.figure(figsize=(13.5, n_sp * 3.2 + 1.0), facecolor='white')
    fig.suptitle(f"Habitat Suitability: {model_a} vs {model_b}", fontsize=14, fontweight='bold', color='black', y=0.98)
    gs = GridSpec(n_sp, 3, figure=fig, hspace=0.35, wspace=0.12, top=0.93, bottom=0.08, left=0.06, right=0.97)

    for row_i, sp in enumerate(species_list):
        try: df = load_species_predictions(sp)
        except Exception: continue
        
        ca, cb = f"HSS_{model_a}", f"HSS_{model_b}"
        if ca not in df.columns or cb not in df.columns: continue
        
        lon, lat, ha, hb = df["lon"].values, df["lat"].values, df[ca].values, df[cb].values
        dv = ha - hb
        vmax = max(np.nanmax(ha), np.nanmax(hb), 1.0)
        d_abs = max(abs(np.nanmin(dv)), abs(np.nanmax(dv)), 0.01)
        
        ax0 = fig.add_subplot(gs[row_i, 0], projection=ccrs.PlateCarree())
        plot_interpolated_map(ax0, lon, lat, ha, HSS_CMAP, 0, vmax, title=(model_a if row_i==0 else ""), show_cbar=(row_i==n_sp-1), cbar_label="HSS")
        ax0.text(-0.12, 0.5, format_species_name(sp), transform=ax0.transAxes, fontsize=9, color='black', rotation=90, va='center', ha='center', fontweight='bold', path_effects=[pe.withStroke(linewidth=2, foreground='white')])
        
        ax1 = fig.add_subplot(gs[row_i, 1], projection=ccrs.PlateCarree())
        plot_interpolated_map(ax1, lon, lat, hb, HSS_CMAP, 0, vmax, title=(model_b if row_i==0 else ""), show_cbar=(row_i==n_sp-1), cbar_label="HSS")
        
        ax2 = fig.add_subplot(gs[row_i, 2], projection=ccrs.PlateCarree())
        plot_interpolated_map(ax2, lon, lat, dv, DIFF_CMAP, -d_abs, d_abs, title=(f"{model_a} − {model_b}" if row_i==0 else ""), show_cbar=(row_i==n_sp-1), cbar_label="ΔHSS")
        
    return fig

if __name__ == "__main__":
    os.makedirs(OUT_DIR, exist_ok=True)
    print("=" * 60)
    print("  CAST Fig 7 — Species Distribution Prediction Maps")
    for sp in SPECIES_LIST:
        print(f"\n  Processing: {sp}")
        try: df = load_species_predictions(sp)
        except Exception as e:
            print("   ", e); continue
        
        fig = plot_species_full(sp, df, load_cate_data(sp))
        if fig is not None:
            out = os.path.join(OUT_DIR, f"fig7_{sp}_{LAYOUT_MODE}.png")
            fig.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor=fig.get_facecolor())
            plt.close(fig)
            print(f"    → Saved: {out}")
            
    print(f"\n  Generating multi-species comparison: {MODELS_TO_COMPARE[0]} vs baselines...")
    for baseline in MODELS_TO_COMPARE[1:]:
        fig = plot_multi_species_comparison(SPECIES_LIST, MODELS_TO_COMPARE[0], baseline)
        if fig is not None:
            out = os.path.join(OUT_DIR, f"fig7_comparison_{MODELS_TO_COMPARE[0]}_vs_{baseline}.png")
            fig.savefig(out, dpi=FIG_DPI, bbox_inches='tight', facecolor=fig.get_facecolor())
            plt.close(fig)
            print(f"    → Saved: {out}")
    print("\n  All figures generated!")
