#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig 6 组合图：3×3 排布，横向（列）标题为物种名称，每列为该物种的环境变量 CATE 图。

- 每张子图：去掉标题和图例，仅保留变量名称文字；四边黑色方框。
- 全图底部一个共享 colorbar；图例文字与效应范围在脚本顶部 CBAR_* / CATE_* 常量中自定义。
- 输出 1200 dpi PNG/SVG，英文标注，Arial 字体。

图例与效应范围说明、CATE 热图扩展玩法见：
  figures/case2_eco/tables/CATE_heatmap_visualization_notes.md

依赖：与 fig6_spatial_cate_maps.py 相同（cartopy, geopandas, scipy, matplotlib, pandas）。

运行：cd E:/CausalSDMs && python scripts/EcoISEA3H/plot/fig6_spatial_cate_maps_3x3.py
"""

import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.path import Path as MplPath
from matplotlib.ticker import FixedLocator, FixedFormatter
from scipy.interpolate import griddata

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import geopandas as gpd

warnings.filterwarnings("ignore")

# ══════════════════════════════════════════════════════════════════════════════
# 路径与参数（与 fig6 一致）
# ══════════════════════════════════════════════════════════════════════════════
BASE_DIR = "E:/CausalSDMs"
os.chdir(BASE_DIR)

CATE_CSV   = "output/case2_eco/all_spatial_cate_v3.csv"
CHINA_SHP  = "plot-function-main/data/china.shp"
DASH_SHP   = "plot-function-main/data/dashline.shp"
FIG_DIR    = "figures/case2_eco/plot/cate_maps"
OUT_NAME   = "fig6_cate_3x3_composite"

INTERP_RES    = 0.06
INTERP_METHOD = "nearest"
DISPLAY_RES   = 0.02   # 子图降采样加快绘制
DPI           = 2400
CHINA_EXTENT  = [73.5, 135, 18, 53.5]

# ══════════════════════════════════════════════════════════════════════════════
# 底部图例（colorbar）自定义
# ══════════════════════════════════════════════════════════════════════════════
# 效应范围模式：
#   "relative" = 每张子图按自身 CATE 范围归一化，亮度一致；色条仅表示方向（-1/0/1）
#   "global"   = 全图统一 CATE 刻度，便于跨图比较绝对值；色条显示实际效应范围
CATE_SCALE_MODE = "relative"   # "relative" | "global"

# 仅当 CATE_SCALE_MODE == "global" 时生效：
#   None = 自动取全 8 张图 |CATE| 的 98 分位数作为 ±lim
#   正数 = 手动指定 ±lim，如 0.08 表示色条范围 [-0.08, 0.08]
CATE_GLOBAL_LIM = 0.08

# 图例主标题（色条下方一行字）
CBAR_LABEL = "CATE"

# 色条刻度位置（必须是数值）。None = 默认 3 档
CBAR_TICKS = [-0.08, -0.04, 0, 0.04, 0.08]

# 色条刻度上显示的文字（与 CBAR_TICKS 数量一致）。None = 用数值自动转成字符串
CBAR_TICK_LABELS = ["-0.08", "-0.04", "0", "0.04", "0.08"]

# 若希望 global 模式下色条显示数值，可设为 True，脚本会用实际 ±lim 生成刻度标签
CBAR_SHOW_NUMERIC_IN_GLOBAL = True

# 变量英文显示名（仅用于子图内变量文字）
VAR_LABELS = {
    "aridityindexthornthwaite": "Aridity Index",
    "bio02":                    "Diurnal Range (Bio02)",
    "bio15":                    "Precip. Seasonality (Bio15)",
    "bio19":                    "Precip. Coldest Qtr (Bio19)",
    "elevation":                "Elevation",
    "etccdi_cwd":               "Consecutive Wet Days",
    "landcover_igbp":           "Land Cover (IGBP)",
    "maxtempcoldest":           "Tmax Coldest Month",
    "nontree":                  "Non-tree Vegetation",
    "topowet":                  "Topographic Wetness",
    "tri":                      "Terrain Ruggedness",
}

def get_var_label(v):
    return VAR_LABELS.get(v, v)

def fmt_species(s):
    return s.replace("_", " ")

# 3×3 排布定义：列 = 物种（横向标题），行 = 第 1/2/3 个变量；(row, col) -> (species, variable)
# 顺序：Macaca_mulatta, Ovis_ammon, Rhinopithecus_roxellana
GRID_SPEC = [
    # row0: 各物种第 1 个变量
    ("Macaca_mulatta", "bio19"),
    ("Ovis_ammon", "bio02"),
    ("Rhinopithecus_roxellana", "aridityindexthornthwaite"),
    # row1: 各物种第 2 个变量
    ("Macaca_mulatta", "etccdi_cwd"),
    ("Ovis_ammon", "elevation"),
    ("Rhinopithecus_roxellana", "maxtempcoldest"),
    # row2: 各物种第 3 个变量（最后一格无图）
    ("Macaca_mulatta", "maxtempcoldest"),
    ("Ovis_ammon", "maxtempcoldest"),
    None,  # (2,2) 空白
]

# ══════════════════════════════════════════════════════════════════════════════
# 地图工具（从 fig6 复用）
# ══════════════════════════════════════════════════════════════════════════════
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
    if isinstance(china_geom, MultiPolygon):
        polys = list(china_geom.geoms)
    else:
        polys = [china_geom]
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
    if not geoms:
        return
    ax.add_geometries(geoms, ccrs.PlateCarree(), **kwargs)

def setup_china_axes_subplot(ax, geoms_china, geoms_dash):
    """子图用：仅范围、省界、南海断续线，无海岸线河流湖泊，无标题。"""
    ax.set_extent(CHINA_EXTENT, crs=ccrs.PlateCarree())
    ax.set_aspect('auto')  # 这会允许 subplots_adjust 的 hspace 生效，但会轻微扭曲地图
    if geoms_china:
        add_china_boundary_cached(ax, geoms_china, ec="black", fc="none", linewidth=0.5)
    if geoms_dash:
        add_china_boundary_cached(ax, geoms_dash, ec="black", fc="none", linewidth=0.6)
    ax.patch.set_facecolor("none")
    # 子图保留四边黑色方框
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

# ══════════════════════════════════════════════════════════════════════════════
# 单格数据：插值 + 中国 mask，返回栅格与 98 分位（用于全局 colorbar）
# ══════════════════════════════════════════════════════════════════════════════
def compute_one_grid(lons, lats, cate_vals, grid_lon, grid_lat, china_mask):
    if len(lons) < 30:
        return None, None
    points = np.column_stack([lons, lats])
    grid_lon_2d, grid_lat_2d = np.meshgrid(grid_lon, grid_lat)
    grid_cate = griddata(points, cate_vals, (grid_lon_2d, grid_lat_2d), method=INTERP_METHOD)
    grid_cate_nearest = griddata(points, cate_vals, (grid_lon_2d, grid_lat_2d), method="nearest")
    grid_cate = np.where(np.isnan(grid_cate), grid_cate_nearest, grid_cate)
    out = np.array(grid_cate, dtype=float, copy=True)
    out[~china_mask] = np.nan
    valid = out[~np.isnan(out)]
    if len(valid) == 0:
        return None, None
    cate_lim = np.percentile(np.abs(valid), 98)
    if cate_lim < 1e-8:
        cate_lim = np.max(np.abs(valid)) + 1e-6
    return out, cate_lim

# ══════════════════════════════════════════════════════════════════════════════
# 主流程
# ══════════════════════════════════════════════════════════════════════════════
def main():
    print("Fig 6 组合图 3×3：读取数据与预计算...")
    cate_all = pd.read_csv(CATE_CSV)
    lon_min, lon_max = CHINA_EXTENT[0], CHINA_EXTENT[1]
    lat_min, lat_max = CHINA_EXTENT[2], CHINA_EXTENT[3]
    grid_lon = np.arange(lon_min, lon_max + INTERP_RES * 0.5, INTERP_RES)
    grid_lat = np.arange(lat_min, lat_max + INTERP_RES * 0.5, INTERP_RES)
    china_geom = get_china_geometry()
    china_mask = build_china_mask(grid_lon, grid_lat, china_geom)
    geoms_china = get_boundary_geoms(CHINA_SHP)
    geoms_dash = get_boundary_geoms(DASH_SHP)

    # 为 8 个任务计算栅格并收集全局 cate_lim
    tasks = [x for x in GRID_SPEC if x is not None]
    grids = []
    cate_lims = []
    for (sp, var) in tasks:
        mask = (cate_all["species"] == sp) & (cate_all["variable"] == var)
        sub = cate_all.loc[mask].dropna(subset=["lon", "lat", "cate"])
        g, lim = compute_one_grid(
            sub["lon"].values, sub["lat"].values, sub["cate"].values,
            grid_lon, grid_lat, china_mask
        )
        if g is None:
            grids.append(np.full((len(grid_lat), len(grid_lon)), np.nan))
            cate_lims.append(0.0)
        else:
            grids.append(g)
            cate_lims.append(lim)
    for k in range(len(cate_lims)):
        if cate_lims[k] < 1e-8:
            cate_lims[k] = 0.01

    # 全局刻度（global 模式）：用于跨图比较绝对值
    if CATE_SCALE_MODE == "global":
        global_lim = CATE_GLOBAL_LIM
        if global_lim is None or global_lim <= 0:
            global_lim = max(cate_lims) if cate_lims else 0.01
        if global_lim < 1e-8:
            global_lim = 0.01

    # 绘图降采样
    if DISPLAY_RES and DISPLAY_RES > INTERP_RES:
        step = max(1, int(round(DISPLAY_RES / INTERP_RES)))
        plot_lon = grid_lon[::step]
        plot_lat = grid_lat[::step]
    else:
        step = 1
        plot_lon = grid_lon
        plot_lat = grid_lat

    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad("none")

    # 使用 Arial 无衬线字体（Nature 风格）
    plt.rcParams["font.sans-serif"] = ["Arial", "DejaVu Sans"]
    plt.rcParams["axes.unicode_minus"] = False

    # 3×3 画布，列标题为物种名
    # hspace：行间距（单子图高度的比例），越小越紧凑；0.005 为极小行距，接近贴紧
    # wspace：列间距（单子图宽度的比例）
    fig, axes = plt.subplots(
        3, 3,
        subplot_kw={"projection": ccrs.PlateCarree()},
        figsize=(12, 10),
        # gridspec_kw={"hspace": 0.02, "wspace": 0.08}
    )
    species_col_titles = ["Macaca mulatta", "Ovis ammon", "Rhinopithecus roxellana"]
    for j in range(3):
        axes[0, j].set_title(species_col_titles[j], fontsize=11, fontweight="bold", fontfamily="sans-serif")

    idx = 0
    mesh_for_cbar = None  # 用于 global 模式 colorbar
    for i in range(3):
        for j in range(3):
            ax = axes[i, j]
            spec = GRID_SPEC[i * 3 + j]
            if spec is None:
                ax.set_visible(False)
                continue
            sp, var = spec
            grid_cate = grids[idx]
            lim = cate_lims[idx]
            idx += 1
            if CATE_SCALE_MODE == "global":
                norm_i = matplotlib.colors.TwoSlopeNorm(vmin=-global_lim, vcenter=0, vmax=global_lim)
            else:
                norm_i = matplotlib.colors.TwoSlopeNorm(vmin=-lim, vcenter=0, vmax=lim)
            if step > 1:
                plot_cate = grid_cate[::step, ::step]
            else:
                plot_cate = grid_cate
            mesh = ax.pcolormesh(
                plot_lon, plot_lat, plot_cate,
                transform=ccrs.PlateCarree(),
                cmap=cmap, norm=norm_i,
                shading="auto", rasterized=True, zorder=1,
            )
            mesh_for_cbar = mesh
            setup_china_axes_subplot(ax, geoms_china, geoms_dash)
            var_display = get_var_label(var)
            ax.text(0.02, 0.98, var_display, transform=ax.transAxes,
                    fontsize=9, va="top", ha="left", fontfamily="sans-serif",
                    bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8, edgecolor="none"))

    # 底部 colorbar：直接针对自定义刻度进行排布
    cbar_ax = fig.add_axes([0.22, 0.08, 0.56, 0.024]) # 稍微调高一点位置
    
    # 确定刻度与范围
    if CBAR_TICKS is not None:
        ticks = [float(x) for x in CBAR_TICKS]
        vmin, vmax = min(ticks), max(ticks)
    else:
        if CATE_SCALE_MODE == "global":
            vmin, vmax = -global_lim, global_lim
            ticks = [-global_lim, 0.0, global_lim]
        else:
            vmin, vmax = -1.0, 1.0
            ticks = [-1.0, 0.0, 1.0]

    # 标签对齐
    if CBAR_TICK_LABELS is not None:
        labels = list(CBAR_TICK_LABELS)
    else:
        labels = [f"{t:.3g}" for t in ticks]

    if len(labels) != len(ticks):
        raise ValueError(f"CBAR_TICK_LABELS 长度({len(labels)})必须与 CBAR_TICKS 长度({len(ticks)})一致。")

    # 创建 ScalarMappable，确保其范围与我们想要的颜色映射范围一致
    # 哪怕子图是 relative 的，图例也作为一个“统一标尺”显示
    norm_cbar = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    sm = plt.cm.ScalarMappable(norm=norm_cbar, cmap=cmap)
    sm.set_array([])

    # 绘制 colorbar，使用 extend='both'
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal", extend="both")
    
    # 【核⼼修复】显式设置 colorbar 轴的范围，保证 FixedLocator 的数值对应到正确像素位置
    cbar.ax.set_xlim(vmin, vmax)
    cbar.ax.xaxis.set_major_locator(FixedLocator(ticks))
    cbar.ax.xaxis.set_major_formatter(FixedFormatter(labels))

    cbar.set_label(CBAR_LABEL, fontsize=10, fontfamily="sans-serif", labelpad=5)
    nticks = len(ticks)
    cbar.ax.tick_params(labelsize=8)

    # 保存前再次收紧行/列间距
    # fig.subplots_adjust(hspace=0.04, wspace=0.08)
    fig.subplots_adjust(left=0.05, right=0.95, bottom=0.12, top=0.85, 
                        hspace=0.05, wspace=0.08)

    plt.savefig(
        os.path.join(FIG_DIR, f"{OUT_NAME}.png"),
        dpi=DPI, pad_inches=0.05, facecolor="white",
    )
    try:
        plt.savefig(
            os.path.join(FIG_DIR, f"{OUT_NAME}.svg"),
                pad_inches=0.05, facecolor="white",
        )
    except Exception:
        pass
    plt.close()
    print(f"已保存: {FIG_DIR}/{OUT_NAME}.png (及 .svg)")


if __name__ == "__main__":
    main()
