#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig 6: 空间 CATE 热力图 — Plant 欧洲案例，与 MaskSDM-MEE-main 一致

研究范围: 欧洲 bbox 经度 -10~31°、纬度 36~56°（与 MaskSDM-MEE generate_map_data 一致）。
绘制方式: 不做插值；用 Natural Earth 110m 国家面做陆地面 mask，只保留陆地上的 CATE 点，
  散点绘制（GeoDataFrame + Point + gdf.plot），颜色不溢出海洋，边界由 shp 约束。

输出: figures/case4_plant/plot/cate_maps/fig6_cate_{species}_{variable}.png/svg
依赖: geopandas matplotlib numpy pandas shapely
运行: cd E:/CausalSDMs && python scripts/Plant/plot/fig6_spatial_cate_maps.py
"""

import os
import warnings
import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import geopandas as gpd
from shapely.geometry import Point

warnings.filterwarnings("ignore")

# ══════════════════════════════════════════════════════════════════════════════
# ★ 可配置参数
# ══════════════════════════════════════════════════════════════════════════════
TARGET_SPECIES = None
TARGET_VARS    = None
MARKER_SIZE    = 1.2    # 陆地点散点大小，与 MaskSDM-MEE predictions_maps 一致（约 1）
DPI            = 1200
N_WORKERS      = 8
BACKGROUND_COLOR = "#EAEAEA"

# ══════════════════════════════════════════════════════════════════════════════
# 路径配置
# ══════════════════════════════════════════════════════════════════════════════
BASE_DIR = "E:/CausalSDMs"
os.chdir(BASE_DIR)

# Plant 案例输出 (03_run_Plant_multi_species.R)
CATE_CSV    = "output/case4_plant/all_spatial_cate_plant.csv"
ATE_CSV     = "output/case4_plant/all_ate_results_plant.csv"
FIG_DIR     = "figures/case4_plant/plot/cate_maps"
TBL_DIR     = "figures/case4_plant/tables"
# 欧洲研究范围: 与 MaskSDM-MEE-main 完全一致
EUROPE_EXTENT = [-10.0, 31.0, 36.0, 56.0]   # lon_min, lon_max, lat_min, lat_max
WORLD_SHP   = "110m_cultural/ne_110m_admin_0_countries.shp"

os.makedirs(FIG_DIR, exist_ok=True)
os.makedirs(TBL_DIR, exist_ok=True)

# ══════════════════════════════════════════════════════════════════════════════
# 变量英文显示名称
# ══════════════════════════════════════════════════════════════════════════════
VAR_LABELS = {
    "aridityindexthornthwaite": "Aridity Index",
    "bio02": "Diurnal Range", "bio_2": "Mean Diurnal Range",
    "bio15": "Precip. Seasonality", "bio_15": "Precip. Seasonality",
    "bio19": "Precip. Coldest Qtr", "bio_19": "Precip. Coldest Qtr",
    "bio03": "Isothermality", "bio_3": "Isothermality",
    "bio18": "Precip. Warmest Qtr", "bio_18": "Precip. Warmest Qtr",
    "elevation": "Elevation", "Elevation": "Elevation",
    "etccdi_cwd": "Consecutive Wet Days",
    "landcover_igbp": "Land Cover (IGBP)",
    "maxtempcoldest": "Tmax Coldest Month",
    "nontree": "Non-tree Vegetation",
    "topowet": "Topographic Wetness",
    "tri": "Terrain Ruggedness",
    "Slope": "Slope", "Aspect": "Aspect",
    "ORCDRC": "Soil Organic C", "PHIHOX": "Soil pH",
    "CECSOL": "Soil CEC", "CLYPPT": "Clay Content",
    "SLTPPT": "Silt Content", "BDTICM": "Bulk Density",
    "Lights2009": "Night Lights", "Built2009": "Built-up",
    "Croplands2005": "Croplands", "Pasture2009": "Pasture",
}

def get_var_label(v):
    return VAR_LABELS.get(v, v)

def fmt_species(s):
    return s.replace("_", " ")

# ══════════════════════════════════════════════════════════════════════════════
# 欧洲 CATE 图：与 MaskSDM-MEE 一致 — 陆地面 mask + 散点，无插值
# ══════════════════════════════════════════════════════════════════════════════
_world_gdf_cache = None
_land_union_cache = None

def _get_world():
    """欧洲底图：本地 Natural Earth 110m 国家面 (110m_cultural)。"""
    global _world_gdf_cache
    if _world_gdf_cache is None:
        path = os.path.abspath(os.path.join(BASE_DIR, WORLD_SHP))
        if not os.path.exists(path):
            raise FileNotFoundError(f"World shapefile not found: {path}")
        _world_gdf_cache = gpd.read_file(path)
    return _world_gdf_cache

def _get_land_union():
    """陆地面（国家面 union），用于只保留陆地上的 CATE 点，边界由 shp 约束。"""
    global _land_union_cache
    if _land_union_cache is None:
        _land_union_cache = _get_world().unary_union
    return _land_union_cache

def _mask_to_land(lons, lats, cate_vals):
    """只保留在陆地面内的点，与 MaskSDM-MEE generate_map_data 中 within(world.unary_union) 一致。向量化。"""
    land = _get_land_union()
    if land is None:
        return lons, lats, cate_vals
    pts = gpd.GeoSeries([Point(x, y) for x, y in zip(lons, lats)], crs="EPSG:4326")
    keep = pts.within(land).values
    return np.asarray(lons)[keep], np.asarray(lats)[keep], np.asarray(cate_vals)[keep]

def render_cate_map_europe_masksdm_style(lons, lats, cate_vals, species_name, var_name, ate_coef):
    """
    欧洲单张 CATE 图：底图 world + 仅陆地点散点（无插值），与 MaskSDM-MEE predictions_maps 一致。
    颜色不溢出海洋，边界由 shp 约束。
    """
    xmin, xmax = EUROPE_EXTENT[0], EUROPE_EXTENT[1]
    ymin, ymax = EUROPE_EXTENT[2], EUROPE_EXTENT[3]

    # 只保留陆地上的点，不做任何插值
    lons_land, lats_land, cate_land = _mask_to_land(
        np.asarray(lons, dtype=float),
        np.asarray(lats, dtype=float),
        np.asarray(cate_vals, dtype=float),
    )
    if len(cate_land) < 10:
        return None

    valid = cate_land[~np.isnan(cate_land)]
    if len(valid) == 0:
        return None
    cate_lim = np.percentile(np.abs(valid), 98)
    if cate_lim < 1e-8:
        cate_lim = np.max(np.abs(valid)) + 1e-6
    norm = mcolors.TwoSlopeNorm(vmin=-cate_lim, vcenter=0, vmax=cate_lim)
    cmap = plt.cm.RdBu_r

    fig, ax = plt.subplots(1, 1, figsize=(10, 8), facecolor="white")
    world = _get_world()
    world.plot(ax=ax, color=BACKGROUND_COLOR, edgecolor="none", zorder=0)
    # 散点绘制陆地点，与 MaskSDM-MEE 一致：GeoDataFrame + Point + gdf.plot(column=..., markersize=...)
    gdf = gpd.GeoDataFrame(
        {"cate": cate_land},
        geometry=[Point(lon, lat) for lon, lat in zip(lons_land, lats_land)],
        crs="EPSG:4326",
    )
    gdf.plot(ax=ax, column="cate", cmap=cmap, norm=norm, markersize=MARKER_SIZE, legend=False, zorder=1)
    world.plot(ax=ax, facecolor="none", edgecolor="0.35", linewidth=0.4, zorder=2)
    ax.set_xlim((xmin, xmax))
    ax.set_ylim((ymin, ymax))
    ax.set_aspect("equal")
    ax.grid(False)
    ax.set_xticks([])
    ax.set_yticks([])

    # 为散点图添加 colorbar（用 mappable 从 norm 构造）
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])
    cbar_ax = fig.add_axes([0.25, 0.08, 0.5, 0.02])
    cbar = fig.colorbar(sm, cax=cbar_ax, orientation="horizontal")
    cbar.set_label("CATE (Conditional Average Treatment Effect)", fontsize=9, fontfamily="sans-serif")
    cbar.ax.tick_params(labelsize=8)

    var_display = get_var_label(var_name)
    sp_display = fmt_species(species_name)
    ax.set_title(f"Spatial CATE: {var_display}", fontsize=14, fontweight="bold", pad=10, fontfamily="sans-serif")
    fig.text(0.5, 0.96, f"{sp_display}  |  ATE = {ate_coef:.4f}  |  n = {len(lons_land):,} land points (no interpolation)",
             ha="center", va="bottom", fontsize=9, color="grey", style="italic", fontfamily="sans-serif")
    plt.subplots_adjust(left=0.02, right=0.98, top=0.90, bottom=0.14)
    return fig


def render_cate_map(lons, lats, cate_vals, species_name, var_name, ate_coef):
    """Plant 欧洲案例：陆地面 mask + 散点绘制，无插值（与 MaskSDM-MEE 一致）。"""
    return render_cate_map_europe_masksdm_style(
        lons, lats, cate_vals, species_name, var_name, ate_coef,
    )


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

    # ── 欧洲范围（无插值、无预计算网格）────────────────────────────────────
    print(f"  范围: {EUROPE_EXTENT} (陆地点散点，无插值)\n")
    fig_dir_abs = os.path.abspath(FIG_DIR)
    n_workers = max(0, int(N_WORKERS))

    # ── 批量绘制（单进程或多进程）────────────────────────────────────────────
    success = 0
    if n_workers <= 1:
        for i, row in combos.iterrows():
            sp  = row["species"]
            var = row["variable"]
            ate = row["coef"]
            print(f"[{i+1}/{n_tasks}] {sp} x {var} (ATE={ate:.4f})")
            mask = (cate_all["species"] == sp) & (cate_all["variable"] == var)
            sub = cate_all.loc[mask].dropna(subset=["lon", "lat", "cate"])
            if len(sub) < 30:
                print(f"  ⚠ 数据不足 ({len(sub)} 点), 跳过")
                continue
            try:
                fig = render_cate_map(
                    sub["lon"].values, sub["lat"].values, sub["cate"].values,
                    sp, var, ate,
                )
                if fig is None:
                    print("  ⚠ 陆地点不足或无有效数据, 跳过")
                    continue
                fname = f"fig6_cate_{sp}_{var}"
                fig.savefig(
                    os.path.join(FIG_DIR, f"{fname}.png"),
                    dpi=DPI, bbox_inches="tight", pad_inches=0.02, transparent=True,
                )
                try:
                    fig.savefig(
                        os.path.join(FIG_DIR, f"{fname}.svg"),
                        bbox_inches="tight", pad_inches=0.02, transparent=True,
                    )
                except Exception:
                    pass
                plt.close(fig)
                success += 1
                print(f"  ✓ 已保存")
            except Exception as e:
                print(f"  ✗ 错误: {e}")
                plt.close("all")
    else:
        from concurrent.futures import ProcessPoolExecutor, as_completed
        task_args = []
        for _, row in combos.iterrows():
            sp, var, ate = row["species"], row["variable"], row["coef"]
            mask = (cate_all["species"] == sp) & (cate_all["variable"] == var)
            sub = cate_all.loc[mask].dropna(subset=["lon", "lat", "cate"])
            if len(sub) < 30:
                continue
            task_args.append((
                sp, var, ate,
                sub["lon"].values.copy(), sub["lat"].values.copy(), sub["cate"].values.copy(),
                fig_dir_abs, DPI,
            ))
        print(f"  使用 {n_workers} 个进程并行渲染 {len(task_args)} 张图...\n")
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(_render_one_task, a): a for a in task_args}
            for i, future in enumerate(as_completed(futures)):
                ok, sp, var, err = future.result()
                if ok:
                    success += 1
                    print(f"  ✓ [{i+1}/{len(task_args)}] {sp} x {var}")
                else:
                    print(f"  ✗ [{i+1}/{len(task_args)}] {sp} x {var}  {err or ''}")

    print(f"\n{'=' * 60}")
    print(f"  完成: {success}/{n_tasks} 张 CATE 热力图已保存至 {FIG_DIR}")
    print(f"{'=' * 60}")


def _render_one_task(args):
    """
    单张图渲染任务，供多进程调用。
    args: (sp, var, ate, lons, lats, cate_vals, fig_dir, dpi).
    返回 (success, sp, var, err_msg or None)。
    """
    (sp, var, ate, lons, lats, cate_vals, fig_dir, dpi) = args
    try:
        fig = render_cate_map(lons, lats, cate_vals, sp, var, ate)
        if fig is None:
            return (False, sp, var, "无有效数据")
        fname = f"fig6_cate_{sp}_{var}"
        fig.savefig(
            os.path.join(fig_dir, f"{fname}.png"),
            dpi=dpi, bbox_inches="tight", pad_inches=0.02, transparent=True,
        )
        try:
            fig.savefig(
                os.path.join(fig_dir, f"{fname}.svg"),
                bbox_inches="tight", pad_inches=0.02, transparent=True,
            )
        except Exception:
            pass
        plt.close(fig)
        return (True, sp, var, None)
    except Exception as e:
        plt.close("all")
        return (False, sp, var, str(e))


if __name__ == "__main__":
    main()
