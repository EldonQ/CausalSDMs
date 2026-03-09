#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig 6: 空间 CATE 热力图 — 批量生成，单图单文件

采用 plotbook.ipynb 的渲染方案:
  - cartopy 地图投影 (PlateCarree)
  - scipy.interpolate.griddata 高分辨率插值 → 类 TIF 精度
  - pcolormesh 栅格渲染 (像素级平滑)
  - 中国省界 (china.shp) + 南海断续线 (dashline.shp)
  - 南海小地图 inset
  - 海岸线 / 河流 / 湖泊
  - 经纬网格 + 坐标标注

可配置参数 (脚本头部修改):
  TARGET_SPECIES : 指定物种 (None = 全部)
  TARGET_VARS    : 指定变量 (None = 全部已有 CATE 变量)
  INTERP_RES     : 插值分辨率 (度, 越小越精细)

输出: figures/case2_eco/plot/cate_maps/fig6_cate_{species}_{variable}.png/svg

依赖: pip install cartopy geopandas scipy matplotlib numpy pandas xarray

运行: cd E:/CausalSDMs && python scripts/EcoISEA3H/plot/fig6_spatial_cate_maps.py
"""

import os
import sys
import warnings

# 先设置后端再 import pyplot，避免 GUI/DLL 问题
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
from matplotlib.path import Path as MplPath
from matplotlib.patches import PathPatch
from scipy.interpolate import griddata

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
import geopandas as gpd

warnings.filterwarnings("ignore")

# ══════════════════════════════════════════════════════════════════════════════
# ★ 可配置参数
# ══════════════════════════════════════════════════════════════════════════════
# 绘制一个物种用 ["物种名"]；多个用 ["A", "B"]；None = 全部
TARGET_SPECIES = None        # 例: ["Alces_alces"] 或 ["Ovis_ammon", "Capra_sibirica"]
TARGET_VARS    = None          # None = 全部已有 CATE 变量; 或 ["elevation", "bio19"]
INTERP_RES     = 0.01        # 格点密度：越小点越多；0.06 约 6km，热图点更密
# scipy.griddata 仅支持三种: "nearest"=离散格点 | "linear"=线性平滑 | "cubic"=三次平滑
INTERP_METHOD  = "nearest"
MASK_BUFFER    = 0.03         # 裁剪时边界外扩(度)，避免省界旁出现白边
DPI            = 1200          # 输出分辨率

# ══════════════════════════════════════════════════════════════════════════════
# 路径配置
# ══════════════════════════════════════════════════════════════════════════════
BASE_DIR = "E:/CausalSDMs"
os.chdir(BASE_DIR)

CATE_CSV    = "output/case2_eco/all_spatial_cate_v3.csv"
ATE_CSV     = "output/case2_eco/all_ate_results_v3.csv"
CHINA_SHP   = "plot-function-main/data/china.shp"      # plotbook 目录下的中国省界
DASH_SHP    = "plot-function-main/data/dashline.shp"    # 南海断续线
FIG_DIR     = "figures/case2_eco/plot/cate_maps"
TBL_DIR     = "figures/case2_eco/tables"

os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TBL_DIR, exist_ok=True)

# ══════════════════════════════════════════════════════════════════════════════
# 变量英文显示名称
# ══════════════════════════════════════════════════════════════════════════════
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

# ══════════════════════════════════════════════════════════════════════════════
# 中国地图绘制工具 (参照 plotbook one_map_china / sub_china_map)
# ══════════════════════════════════════════════════════════════════════════════
# 中国全图范围 (与 plotbook plot.py one_map_china 一致: set_extent([70, 140, 15, 55]))
CHINA_EXTENT = [70, 140, 15, 55]
# 南海小地图范围
SCS_EXTENT   = [104.5, 125, 0, 26]

# 中国境内多边形 (用于裁剪 CATE 栅格，等价于 plotbook 的 da.salem.roi(shape=shpfile))
_china_geom = None

def get_china_geometry():
    """加载中国省界 shapefile 并合并为单一多边形，用于栅格裁剪。边界略外扩以消除白边。"""
    global _china_geom
    if _china_geom is not None:
        return _china_geom
    if not os.path.exists(CHINA_SHP):
        return None
    gdf = gpd.read_file(CHINA_SHP)
    if gdf is None or len(gdf) == 0:
        return None
    union = gdf.unary_union
    # 小幅外扩，使贴边格点被保留，避免省界旁出现空白
    _china_geom = union.buffer(MASK_BUFFER) if MASK_BUFFER else union
    return _china_geom

def mask_grid_to_china(grid_lon_1d, grid_lat_1d, grid_values_2d, china_geom):
    """
    将栅格裁剪到中国境内：境外像元置为 NaN。
    等价于 plotbook 中 da.salem.roi(shape=shpfile) 的效果。
    使用向量化 point-in-polygon 加速。
    """
    if china_geom is None:
        return grid_values_2d
    try:
        from shapely.geometry import Polygon, MultiPolygon
        if isinstance(china_geom, MultiPolygon):
            polys = list(china_geom.geoms)
        else:
            polys = [china_geom]
        out = np.array(grid_values_2d, dtype=float, copy=True)
        nlat, nlon = out.shape
        # 所有格点 (lon, lat)
        lon_2d, lat_2d = np.meshgrid(grid_lon_1d, grid_lat_1d)
        points = np.column_stack([lon_2d.ravel(), lat_2d.ravel()])
        inside = np.zeros(len(points), dtype=bool)
        for poly in polys:
            if poly.is_empty:
                continue
            ext = poly.exterior
            if ext is None:
                continue
            from matplotlib.path import Path
            path = Path(np.array(ext.xy).T)
            inside |= path.contains_points(points)
        out.ravel()[~inside] = np.nan
        return out
    except Exception:
        return grid_values_2d

def add_china_boundary(ax, shp_path, **kwargs):
    """叠加中国省界 shapefile"""
    if not os.path.exists(shp_path):
        return
    proj = ccrs.PlateCarree()
    reader = shpreader.Reader(shp_path)
    geometries = reader.geometries()
    ax.add_geometries(geometries, proj, **kwargs)
    reader.close()

def setup_china_axes(ax, add_rivers=True, add_lakes=True,
                     add_coastlines=True, add_gridlines=True):
    """配置中国主图的基础地理要素"""
    ax.set_extent(CHINA_EXTENT, crs=ccrs.PlateCarree())

    # 海岸线
    if add_coastlines:
        ax.coastlines(color="0.3", linewidth=0.6, zorder=3)

    # 河流
    if add_rivers:
        ax.add_feature(cfeature.RIVERS, linewidth=0.3, edgecolor="steelblue",
                        zorder=2)
    # 湖泊
    if add_lakes:
        ax.add_feature(cfeature.LAKES, facecolor="lightblue",
                        edgecolor="steelblue", linewidth=0.2, zorder=2)

    # 中国省界
    add_china_boundary(ax, CHINA_SHP, ec="black", fc="none", linewidth=0.6)
    # 南海断续线
    add_china_boundary(ax, DASH_SHP, ec="black", fc="none", linewidth=0.8)

    # 边框美化
    for spine in ax.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color("0.4")

    # 经纬网格
    if add_gridlines:
        gl = ax.gridlines(
            crs=ccrs.PlateCarree(), draw_labels=True,
            linewidth=0.3, color="grey", alpha=0.5, linestyle="--",
            xlocs=np.arange(70, 141, 10),
            ylocs=np.arange(15, 56, 10),
            x_inline=False, y_inline=False,
        )
        gl.top_labels = False
        gl.right_labels = False
        gl.xlabel_style = {"size": 7, "color": "0.3"}
        gl.ylabel_style = {"size": 7, "color": "0.3"}


def add_scs_inset(fig, ax_main, data_lons, data_lats, data_vals, levels, cmap, norm):
    """在主图右下角添加南海小地图 inset，小地图右下角贴合主图右下角"""
    pos = ax_main.get_position()
    inset_w = pos.width * 0.28
    inset_h = pos.height * 0.35
    # 小地图右下角 = 主图右下角：右对齐 pos.x1，底对齐 pos.y0
    inset_x = pos.x1 - inset_w
    inset_y = pos.y0

    ax_inset = fig.add_axes(
        [inset_x, inset_y, inset_w, inset_h],
        projection=ccrs.PlateCarree()
    )
    ax_inset.set_extent(SCS_EXTENT, crs=ccrs.PlateCarree())

    # 绘制同样的 pcolormesh 数据 (如果覆盖南海区域)
    if data_vals is not None and len(data_vals) > 0:
        ax_inset.pcolormesh(
            data_lons, data_lats, data_vals,
            transform=ccrs.PlateCarree(),
            cmap=cmap, norm=norm,
            shading="auto", rasterized=True
        )

    ax_inset.coastlines(color="0.3", linewidth=0.4)
    add_china_boundary(ax_inset, CHINA_SHP, ec="black", fc="none", linewidth=0.5)
    add_china_boundary(ax_inset, DASH_SHP, ec="black", fc="none", linewidth=0.6)

    ax_inset.gridlines(
        draw_labels=False, linewidth=0.1,
        color="gray", alpha=0.6, linestyle="--"
    )

    for spine in ax_inset.spines.values():
        spine.set_linewidth(0.5)
        spine.set_color("0.4")

    return ax_inset


# ══════════════════════════════════════════════════════════════════════════════
# 核心渲染函数: 将散点 CATE 插值为高精度栅格, 以 pcolormesh 渲染
# ══════════════════════════════════════════════════════════════════════════════
def render_cate_map(lons, lats, cate_vals, species_name, var_name, ate_coef):
    """
    单张 CATE 热力图:
    1. 在中国范围内构建规则网格 (与 plotbook 一致)
    2. scipy.griddata 插值
    3. 按 china.shp 裁剪到境内 (等价 plotbook 的 salem.roi(shape=shpfile))
    4. pcolormesh 仅绘制境内像元，避免越界
    """
    # ── 网格覆盖中国全境 (与 plotbook one_map_china set_extent([70,140,15,55]) 一致) ──
    lon_min, lon_max = CHINA_EXTENT[0], CHINA_EXTENT[1]
    lat_min, lat_max = CHINA_EXTENT[2], CHINA_EXTENT[3]
    grid_lon = np.arange(lon_min, lon_max + INTERP_RES * 0.5, INTERP_RES)
    grid_lat = np.arange(lat_min, lat_max + INTERP_RES * 0.5, INTERP_RES)
    grid_lon_2d, grid_lat_2d = np.meshgrid(grid_lon, grid_lat)

    # ── scipy.griddata 插值 ──────────────────────────────────────────────────
    points = np.column_stack([lons, lats])
    grid_cate = griddata(points, cate_vals, (grid_lon_2d, grid_lat_2d),
                         method=INTERP_METHOD)

    # ── 裁剪到中国境内 (境外置 NaN，与 plotbook damask=da.salem.roi(shape=shpfile) 等价) ──
    china_geom = get_china_geometry()
    grid_cate = mask_grid_to_china(grid_lon, grid_lat, grid_cate, china_geom)

    # ── 对称发散色标 (以 0 为中心，仅用境内有效值) ─────────────────────────────
    valid = grid_cate[~np.isnan(grid_cate)]
    if len(valid) == 0:
        return None
    cate_lim = np.percentile(np.abs(valid), 98)
    if cate_lim < 1e-8:
        cate_lim = np.max(np.abs(valid)) + 1e-6

    # Nature 风格发散配色: 蓝-白-红
    cmap = plt.cm.RdBu_r
    norm = matplotlib.colors.TwoSlopeNorm(vmin=-cate_lim, vcenter=0, vmax=cate_lim)

    # ── 创建画布 ─────────────────────────────────────────────────────────────
    fig = plt.figure(figsize=(10, 8), facecolor="white")
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

    # ── pcolormesh 渲染 (plotbook 核心方法) ───────────────────────────────────
    mesh = ax.pcolormesh(
        grid_lon, grid_lat, grid_cate,
        transform=ccrs.PlateCarree(),
        cmap=cmap, norm=norm,
        shading="auto",
        rasterized=True,          # 关键: 矢量化输出中保持栅格精细
        zorder=1,
    )

    # ── 叠加地理要素 ─────────────────────────────────────────────────────────
    setup_china_axes(ax)

    # ── colorbar ─────────────────────────────────────────────────────────────
    cbar = fig.colorbar(
        mesh, ax=ax, orientation="horizontal",
        fraction=0.046, pad=0.08, shrink=0.7,
        extend="both",
    )
    cbar.set_label("CATE (Conditional Average Treatment Effect)",
                   fontsize=9, labelpad=4, fontfamily="sans-serif")
    cbar.ax.tick_params(labelsize=8)

    # ── 标题（先设好，避免与 colorbar/小地图重叠）────────────────────────────
    var_display = get_var_label(var_name)
    sp_display  = fmt_species(species_name)
    ax.set_title(
        f"Spatial CATE: {var_display}",
        fontsize=14, fontweight="bold", pad=10, fontfamily="sans-serif"
    )
    subtitle = (f"{sp_display}  |  ATE = {ate_coef:.4f}  |  "
                f"n = {len(lons):,} grid cells  |  "
                f"Interp. {INTERP_RES}°")
    # 副标题置于主标题上方，留出间距避免与主标题重叠
    ax.text(
        0.5, 1.03, subtitle,
        transform=ax.transAxes, ha="center", va="bottom",
        fontsize=9, color="grey", fontstyle="italic",
        fontfamily="sans-serif"
    )

    # 先收紧边距，再按最终主图位置放置小地图，避免小地图错位；top 留足给标题
    plt.subplots_adjust(left=0.04, right=0.96, top=0.90, bottom=0.14)
    fig.canvas.draw()
    # 南海小地图：必须在 subplots_adjust 之后添加，位置才贴合主图右下角
    add_scs_inset(fig, ax, grid_lon, grid_lat, grid_cate, None, cmap, norm)
    return fig


# ══════════════════════════════════════════════════════════════════════════════
# 主流程: 读取数据 → 批量渲染
# ══════════════════════════════════════════════════════════════════════════════
def main():
    print("=" * 60)
    print("  Fig 6: 空间 CATE 热力图 (plotbook 渲染方案)")
    print("=" * 60)

    # ── 读取数据 ─────────────────────────────────────────────────────────────
    print("\n读取 CATE 数据...")
    cate_all = pd.read_csv(CATE_CSV)
    ate_all  = pd.read_csv(ATE_CSV)
    ate_all["coef"] = pd.to_numeric(ate_all["coef"], errors="coerce")

    print(f"  CATE: {len(cate_all):,} 条, "
          f"{cate_all['species'].nunique()} 物种, "
          f"{cate_all['variable'].nunique()} 变量")

    # ── 筛选绘图任务 ─────────────────────────────────────────────────────────
    combos = cate_all[["species", "variable"]].drop_duplicates()

    if TARGET_SPECIES is not None:
        combos = combos[combos["species"].isin(TARGET_SPECIES)]
    if TARGET_VARS is not None:
        combos = combos[combos["variable"].isin(TARGET_VARS)]

    # 合并 ATE 信息
    combos = combos.merge(
        ate_all[["species", "variable", "coef"]],
        on=["species", "variable"], how="left"
    )
    combos["coef"] = combos["coef"].fillna(0)
    combos = combos.sort_values(["species", "variable"]).reset_index(drop=True)

    n_tasks = len(combos)
    print(f"  待绘制: {n_tasks} 张 CATE 热力图\n")

    # 保存任务清单
    combos.to_csv(os.path.join(TBL_DIR, "fig6_cate_plot_tasks.csv"), index=False)

    # ── 批量绘制 ─────────────────────────────────────────────────────────────
    success = 0
    for i, row in combos.iterrows():
        sp  = row["species"]
        var = row["variable"]
        ate = row["coef"]

        print(f"[{i+1}/{n_tasks}] {sp} x {var} (ATE={ate:.4f})")

        # 提取该物种×变量的 CATE 数据
        mask = (cate_all["species"] == sp) & (cate_all["variable"] == var)
        sub = cate_all.loc[mask].dropna(subset=["lon", "lat", "cate"])

        if len(sub) < 30:
            print(f"  ⚠ 数据不足 ({len(sub)} 点), 跳过")
            continue

        try:
            fig = render_cate_map(
                sub["lon"].values, sub["lat"].values, sub["cate"].values,
                sp, var, ate
            )

            if fig is None:
                print("  ⚠ 插值后无有效数据, 跳过")
                continue

            fname = f"fig6_cate_{sp}_{var}"
            fig.savefig(
                os.path.join(FIG_DIR, f"{fname}.png"),
                dpi=DPI, bbox_inches="tight", pad_inches=0.02, facecolor="white"
            )
            # SVG 版本
            try:
                fig.savefig(
                    os.path.join(FIG_DIR, f"{fname}.svg"),
                    bbox_inches="tight", pad_inches=0.02, facecolor="white"
                )
            except Exception:
                pass

            plt.close(fig)
            success += 1
            print(f"  ✓ 已保存")

        except Exception as e:
            print(f"  ✗ 错误: {e}")
            plt.close("all")

    print(f"\n{'=' * 60}")
    print(f"  完成: {success}/{n_tasks} 张 CATE 热力图已保存至 {FIG_DIR}")
    print(f"{'=' * 60}")


if __name__ == "__main__":
    main()
