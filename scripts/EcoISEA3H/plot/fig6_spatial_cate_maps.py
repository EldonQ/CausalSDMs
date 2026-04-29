#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig 6: 空间 CATE 热力图 — 批量生成，单图单文件

采用 plotbook.ipynb 的渲染方案:
  - cartopy 地图投影 (PlateCarree)
  - scipy.interpolate.griddata 高分辨率插值 → 类 TIF 精度
  - pcolormesh 栅格渲染 (像素级平滑)
  - 中国省界 (china.shp) + 南海断续线 (dashline.shp)
  - 可选：南海小地图 inset（SHOW_SCS_INSET，默认关闭）
  - 海岸线 / 河流 / 湖泊
  - 经纬网格 + 坐标标注
  - CATE 最终显示范围：中国省界 mask + 各 SDM 的 HSS 适宜区（CAST/BRT/Maxent/RF 分别出图）

可配置参数 (脚本头部修改):
  TARGET_SPECIES : 指定物种 (None = 全部)
  TARGET_VARS    : 指定变量 (None = 全部已有 CATE 变量)
  INTERP_RES     : 插值分辨率 (度, 越小越精细)
  USE_HSS_MASK   : 是否按 HSS 阈值截断（需 spatial_predictions/pred_*.csv）
  HSS_CONF_THRESHOLD : HSS 生境置信下限（格点 HSS < 该值则不绘制 CATE）

输出: figures/case2_eco/plot/cate_maps/fig6_cate_{species}_{variable}_maskHSS_{CAST|BRT|Maxent|RF}.png/svg

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
import geopandas as gpd

warnings.filterwarnings("ignore")
    
# ══════════════════════════════════════════════════════════════════════════════
# ★ 可配置参数
# ══════════════════════════════════════════════════════════════════════════════
# 绘制一个物种用 ["物种名"]；多个用 ["A", "B"]；None = 全部
# Rhinopithecus roxellana (Cercopithecidae)
# Ovis ammon (Bovidae)
# Macaca mulatta (Cercopithecidae)

TARGET_SPECIES = ["Rhinopithecus_roxellana", "Ovis_ammon", "Macaca_mulatta"]   # 例: ["Alces_alces"] 或 ["Ovis_ammon", "Capra_sibirica"]
TARGET_VARS    = None          # None = 全部已有 CATE 变量; 或 ["elevation", "bio19"]
INTERP_RES     = 0.01       # 格点密度：越小点越多；0.06 约 6km，热图点更密
# scipy.griddata 仅支持三种: "nearest"=最快 | "linear"=折中 | "cubic"=最慢最平滑
INTERP_METHOD  = "cubic"
MASK_BUFFER    = 0.0         # 裁剪时边界外扩(度)；不再外扩，通过 nearest 填充彻底消除白边，避免热图移除边界
DPI            = 1200       # 输出分辨率
SHOW_SCS_INSET = False      # True=主图右下角绘制南海小地图；False=不绘制

# ★ HSS 生境掩膜（与 fig7/fig8 的 pred_{species}.csv 一致）
USE_HSS_MASK       = True
PRED_DIR           = "output/case2_eco/spatial_predictions"
# 列名 HSS_CAST / HSS_BRT / HSS_Maxent / HSS_RF（Maxent 即 MaxNet 类输出）
HSS_MODELS         = ["CAST", "BRT", "Maxent", "RF"]
HSS_CONF_THRESHOLD = 0.1   # 仅在该 HSS 以上绘制 CATE；可按文献或 TSS 最优点调整
HSS_INTERP_METHOD  = "cubic"  # HSS 栅格化用 linear，避免 nearest 把适宜区外扩到凸包外

# ★ 速度优化
N_WORKERS      = 8         # 并行进程数，0 或 1=单进程，4/8 等=多进程加速
FAST_PLOT      = True       # True=低分辨率海岸线、不画河流湖泊，出图更快
DISPLAY_RES    = 0.0        # 若设为浮点(如 0.02)，绘图用此分辨率降采样，加快 pcolormesh

# ══════════════════════════════════════════════════════════════════════════════
# 路径配置
# ══════════════════════════════════════════════════════════════════════════════
BASE_DIR = "E:/CausalSDMs"
os.chdir(BASE_DIR)

CATE_CSV    = "output/case2_eco/all_spatial_cate_v3.csv"
ATE_CSV     = "output/case2_eco/all_ate_results_v3.csv"
CHINA_SHP   = "refPackage/plot-function-main/data/china.shp"      # plotbook 目录下的中国省界
DASH_SHP    = "refPackage/plot-function-main/data/dashline.shp"    # 南海断续线
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
# 仅中国区范围（收紧为大陆+南海，不包含蒙古/俄罗斯/印度等周边国家视野）
CHINA_EXTENT = [73.5, 135, 18, 53.5]
# 南海小地图范围
SCS_EXTENT   = [104.5, 125, 0, 26]

# 中国境内多边形 (用于裁剪 CATE 栅格)；缓存 key 含 MASK_BUFFER/INTERP_RES 以便参数变更后重算
_china_geom_cache = None  # (buffer_used, geometry)
# 省界/南海 shapefile 几何缓存，避免每张图重复读盘
_boundary_cache = {}

def _normalize_china_polygonal(geom):
    """将省界 union 结果规范为面几何，供 contains_xy 使用（含 GeometryCollection 拆出多边形）。"""
    if geom is None or getattr(geom, "is_empty", True):
        return None
    gt = geom.geom_type
    if gt in ("Polygon", "MultiPolygon"):
        return geom
    if gt == "GeometryCollection":
        from shapely.ops import unary_union
        parts = [g for g in geom.geoms if g.geom_type in ("Polygon", "MultiPolygon")]
        if not parts:
            return None
        return unary_union(parts)
    return geom

def get_china_geometry():
    """
    加载中国省界 shapefile 并合并为单一多边形，用于栅格裁剪。
    外扩 = MASK_BUFFER + 半格(0.5*INTERP_RES)，使 pcolormesh 格心贴省界时仍被保留，消除边缘空白。
    """
    global _china_geom_cache
    # 外扩量：用户设置的 MASK_BUFFER + 半格补偿（避免格心刚好在边界外导致白边）
    effective_buffer = (MASK_BUFFER or 0) + 0.5 * INTERP_RES
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
    geom = _normalize_china_polygonal(geom)
    _china_geom_cache = (effective_buffer, geom)
    return geom

def get_boundary_geoms(shp_path):
    """缓存 shapefile 几何列表，避免每张图重复读盘。返回 list 供 add_geometries 使用。"""
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
    """对 (N,2) 的 points_xy 做境内判断；优先 Shapely contains_xy（含孔洞/多部件），失败则回退到外轮廓。"""
    if china_geom is None:
        return np.ones(len(points_xy), dtype=bool)
    geom = _normalize_china_polygonal(china_geom)
    if geom is None:
        return np.zeros(len(points_xy), dtype=bool)
    try:
        from shapely import contains_xy
        x = np.ascontiguousarray(points_xy[:, 0], dtype=float)
        y = np.ascontiguousarray(points_xy[:, 1], dtype=float)
        return np.asarray(contains_xy(geom, x, y), dtype=bool)
    except Exception:
        from shapely.geometry import MultiPolygon
        if isinstance(geom, MultiPolygon):
            polys = list(geom.geoms)
        else:
            polys = [geom]
        inside = np.zeros(len(points_xy), dtype=bool)
        for poly in polys:
            if poly.is_empty or poly.exterior is None:
                continue
            path = MplPath(np.array(poly.exterior.xy).T)
            inside |= path.contains_points(points_xy)
        return inside

def mask_grid_to_china(grid_lon_1d, grid_lat_1d, grid_values_2d, china_geom):
    """
    将栅格裁剪到中国境内：境外像元置为 NaN。
    按 pcolormesh 的格元中心做 point-in-polygon，避免省界内侧出现白边。
    """
    if china_geom is None:
        return grid_values_2d
    try:
        nlat, nlon = grid_values_2d.shape
        # 格元中心（与 pcolormesh 的 cell 一一对应）
        center_lon = (grid_lon_1d[:-1] + grid_lon_1d[1:]) * 0.5
        center_lat = (grid_lat_1d[:-1] + grid_lat_1d[1:]) * 0.5
        lon_c, lat_c = np.meshgrid(center_lon, center_lat)
        points = np.column_stack([lon_c.ravel(), lat_c.ravel()])
        inside_c = _points_inside_china(points, china_geom).reshape(nlat - 1, nlon - 1)
        # 扩成 (nlat, nlon)，与 grid_values_2d 对齐：最后一列/行沿用前一格
        mask = np.zeros((nlat, nlon), dtype=bool)
        mask[: nlat - 1, : nlon - 1] = inside_c
        mask[nlat - 1, :] = mask[nlat - 2, :]
        mask[:, nlon - 1] = mask[:, nlon - 2]
        out = np.array(grid_values_2d, dtype=float, copy=True)
        out[~mask] = np.nan
        return out
    except Exception:
        return grid_values_2d

def build_china_mask(grid_lon_1d, grid_lat_1d, china_geom):
    """
    一次性计算栅格在中国境内的布尔 mask（True=境内）。
    按 pcolormesh 的格元中心判断，与 mask_grid_to_china 一致，消除边缘空白。
    """
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

def find_pred_csv(species_name, pred_dir):
    """定位 spatial_predictions 下的 pred_{species}.csv（兼容空格/下划线命名）。"""
    if not pred_dir or not os.path.isdir(pred_dir):
        return None
    s = (species_name or "").strip()
    candidates = [
        os.path.join(pred_dir, f"pred_{s}.csv"),
        os.path.join(pred_dir, f"pred_{s.replace(' ', '_')}.csv"),
        os.path.join(pred_dir, f"pred_{s.replace(' ', '')}.csv"),
    ]
    for p in candidates:
        if os.path.isfile(p):
            return p
    return None


def hss_column_for_model(model):
    """pred 表中的列名（与 fig7/fig8 一致）；MaxNet 类输出在表里记为 HSS_Maxent。"""
    key = (model or "").strip()
    aliases = {
        "CAST": "HSS_CAST",
        "BRT": "HSS_BRT",
        "Maxent": "HSS_Maxent",
        "MAXENT": "HSS_Maxent",
        "maxent": "HSS_Maxent",
        "MaxNet": "HSS_Maxent",
        "maxnet": "HSS_Maxent",
        "RF": "HSS_RF",
    }
    if key in aliases:
        return aliases[key]
    return f"HSS_{key}"


def build_hss_suitable_mask_on_grid(pred_df, hss_col, threshold, grid_lon, grid_lat, interp_method):
    """
    将 pred 表中的 HSS 插值到与 CATE 相同的规则网格，返回与 china_mask 同形状的 bool：
    True 表示该格点 HSS 插值有效且 >= threshold（生境置信/适宜区内）。
    使用 linear/cubic 等时不再用 nearest 回填，避免适宜区沿凸包错误外延。
    """
    if pred_df is None or len(pred_df) == 0 or hss_col not in pred_df.columns:
        return None
    lon = pd.to_numeric(pred_df.get("lon"), errors="coerce")
    lat = pd.to_numeric(pred_df.get("lat"), errors="coerce")
    hss = pd.to_numeric(pred_df[hss_col], errors="coerce")
    ok = lon.notna() & lat.notna() & hss.notna()
    lon = lon[ok].to_numpy(dtype=float)
    lat = lat[ok].to_numpy(dtype=float)
    hss = hss[ok].to_numpy(dtype=float)
    if len(hss) < 30:
        return None
    grid_lon_2d, grid_lat_2d = np.meshgrid(grid_lon, grid_lat)
    points = np.column_stack([lon, lat])
    try:
        g = griddata(points, hss, (grid_lon_2d, grid_lat_2d), method=interp_method)
    except Exception:
        g = griddata(points, hss, (grid_lon_2d, grid_lat_2d), method="linear")
    suitable = np.isfinite(g) & (g >= float(threshold))
    return suitable


def precompute_hss_plot_masks(species_names, grid_lon, grid_lat, china_mask, pred_dir, models, threshold, interp_method):
    """
    返回 dict[(species, model)] = 与 china_mask 同形状的 bool（中国 ∧ HSS 适宜区）。
    """
    masks = {}
    for sp in species_names:
        path = find_pred_csv(sp, pred_dir)
        if path is None:
            print(f"  [HSS] 未找到预测文件: pred_{sp}.csv （已查 PRED_DIR={pred_dir}）")
            continue
        try:
            pdf = pd.read_csv(path)
        except Exception as e:
            print(f"  [HSS] 读取失败 {path}: {e}")
            continue
        for model in models:
            col = hss_column_for_model(model)
            if col not in pdf.columns:
                alt = [c for c in pdf.columns if str(c).upper().startswith("HSS_")]
                print(f"  [HSS] {sp} 无列 {col}，可用 HSS 列: {alt[:8]}{'...' if len(alt) > 8 else ''}")
                continue
            sm = build_hss_suitable_mask_on_grid(pdf, col, threshold, grid_lon, grid_lat, interp_method)
            if sm is None:
                continue
            if sm.shape != china_mask.shape:
                print(f"  [HSS] {sp} {model} 掩膜形状 {sm.shape} 与网格 {china_mask.shape} 不一致，跳过。")
                continue
            masks[(sp, model)] = china_mask & sm.astype(bool)
    return masks


def add_china_boundary(ax, shp_path, **kwargs):
    """叠加中国省界 shapefile（未缓存时读盘）"""
    if not os.path.exists(shp_path):
        return
    proj = ccrs.PlateCarree()
    reader = shpreader.Reader(shp_path)
    geometries = reader.geometries()
    ax.add_geometries(geometries, proj, **kwargs)
    reader.close()

def add_china_boundary_cached(ax, geoms, **kwargs):
    """叠加省界/南海，使用已缓存的几何列表，避免重复读盘"""
    if not geoms:
        return
    ax.add_geometries(geoms, ccrs.PlateCarree(), **kwargs)

def setup_china_axes(ax, add_rivers=True, add_lakes=True,
                     add_coastlines=True, add_gridlines=True,
                     geoms_china=None, geoms_dash=None, fast_plot=False):
    """
    配置中国主图：仅中国区范围 + 省界 + 南海断续线；保留河流、湖泊；不画全球海岸线（避免周边国家轮廓）。
    无经纬度、无边框、无网格线。
    """
    ax.set_extent(CHINA_EXTENT, crs=ccrs.PlateCarree())

    # 不画全球海岸线，仅用省界+南海断续线勾勒中国，实现“仅保留中国区”
    if add_coastlines:
        res = "110m" if fast_plot else "50m"
        try:
            ax.coastlines(resolution=res, color="0.3", linewidth=0.6, zorder=3)
        except Exception:
            ax.coastlines(color="0.3", linewidth=0.6, zorder=3)

    if not fast_plot:
        if add_rivers:
            ax.add_feature(cfeature.RIVERS, linewidth=0.3, edgecolor="steelblue", zorder=2)
        if add_lakes:
            ax.add_feature(cfeature.LAKES, facecolor="lightblue",
                           edgecolor="steelblue", linewidth=0.2, zorder=2)

    # 中国省界、南海断续线（优先用缓存）
    if geoms_china is not None:
        add_china_boundary_cached(ax, geoms_china, ec="black", fc="none", linewidth=0.6)
    else:
        add_china_boundary(ax, CHINA_SHP, ec="black", fc="none", linewidth=0.6)
    if geoms_dash is not None:
        add_china_boundary_cached(ax, geoms_dash, ec="black", fc="none", linewidth=0.8)
    else:
        add_china_boundary(ax, DASH_SHP, ec="black", fc="none", linewidth=0.8)

    # 除中国区外不绘制：主图背景透明（境外区域透明）
    ax.patch.set_facecolor("none")
    ax.patch.set_edgecolor("white")
    ax.patch.set_linewidth(0)
    # 去掉经纬度、网格线、地图黑框
    ax.set_frame_on(False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    _outline = getattr(ax, "outline_patch", None)
    if _outline is not None:
        _outline.set_visible(False)
        _outline.set_edgecolor("none")
    _bg = getattr(ax, "background_patch", None)
    if _bg is not None:
        _bg.set_visible(False)


def add_scs_inset(fig, ax_main, data_lons, data_lats, data_vals, levels, cmap, norm,
                  geoms_china=None, geoms_dash=None):
    """在主图右下角添加南海小地图 inset，带黑框，位置适当下移。"""
    pos = ax_main.get_position()
    inset_w = pos.width * 0.28
    inset_h = pos.height * 0.35
    inset_x = pos.x1 - inset_w
    # 适当下移：相对主图底部再向下偏移一截
    inset_y = pos.y0 - 0.04

    ax_inset = fig.add_axes(
        [inset_x, inset_y, inset_w, inset_h],
        projection=ccrs.PlateCarree()
    )
    ax_inset.set_extent(SCS_EXTENT, crs=ccrs.PlateCarree())

    if data_vals is not None and len(data_vals) > 0:
        ax_inset.pcolormesh(
            data_lons, data_lats, data_vals,
            transform=ccrs.PlateCarree(),
            cmap=cmap, norm=norm,
            shading="auto", rasterized=True
        )

    ax_inset.coastlines(color="0.3", linewidth=0.4)
    if geoms_china is not None:
        add_china_boundary_cached(ax_inset, geoms_china, ec="black", fc="none", linewidth=0.5)
    else:
        add_china_boundary(ax_inset, CHINA_SHP, ec="black", fc="none", linewidth=0.5)
    if geoms_dash is not None:
        add_china_boundary_cached(ax_inset, geoms_dash, ec="black", fc="none", linewidth=0.6)
    else:
        add_china_boundary(ax_inset, DASH_SHP, ec="black", fc="none", linewidth=0.6)

    # 南海小图保留黑色矩形边框
    ax_inset.set_frame_on(True)
    ax_inset.patch.set_edgecolor("black")
    ax_inset.patch.set_linewidth(1)
    ax_inset.patch.set_facecolor("white")

    return ax_inset


# ══════════════════════════════════════════════════════════════════════════════
# 核心渲染函数: 将散点 CATE 插值为高精度栅格, 以 pcolormesh 渲染
# ══════════════════════════════════════════════════════════════════════════════
def render_cate_map(lons, lats, cate_vals, species_name, var_name, ate_coef,
                    grid_lon=None, grid_lat=None, china_mask=None,
                    geoms_china=None, geoms_dash=None, display_res=None,
                    plot_hss_model=None, plot_hss_threshold=None):
    """
    单张 CATE 热力图。可选传入预计算的 grid_lon, grid_lat, china_mask 与边界几何以加速。
    china_mask: 最终为 True 的格点才保留 CATE（通常 = 中国 ∧ HSS 适宜区）。
    display_res: 若设置，绘图时降采样到此分辨率以加快 pcolormesh。
    plot_hss_model / plot_hss_threshold: 仅在副标题标注当前 HSS 掩膜来源与阈值。
    """
    lon_min, lon_max = CHINA_EXTENT[0], CHINA_EXTENT[1]
    lat_min, lat_max = CHINA_EXTENT[2], CHINA_EXTENT[3]
    use_precomputed = grid_lon is not None and grid_lat is not None and china_mask is not None

    if use_precomputed:
        grid_lon_2d, grid_lat_2d = np.meshgrid(grid_lon, grid_lat)
    else:
        grid_lon = np.arange(lon_min, lon_max + INTERP_RES * 0.5, INTERP_RES)
        grid_lat = np.arange(lat_min, lat_max + INTERP_RES * 0.5, INTERP_RES)
        grid_lon_2d, grid_lat_2d = np.meshgrid(grid_lon, grid_lat)

    # ── scipy.griddata 插值 ──────────────────────────────────────────────────
    points = np.column_stack([lons, lats])
    grid_cate = griddata(points, cate_vals, (grid_lon_2d, grid_lat_2d),
                         method=INTERP_METHOD)
                         
    # ★ 核心修复：用 nearest 填充 linear 不能覆盖的边界地带 (Convex Hull之外的数据)
    grid_cate_nearest = griddata(points, cate_vals, (grid_lon_2d, grid_lat_2d), method='nearest')
    grid_cate = np.where(np.isnan(grid_cate), grid_cate_nearest, grid_cate)

    # ── 裁剪到中国境内 ───────────────────────────────────────────────────────
    if use_precomputed:
        out = np.array(grid_cate, dtype=float, copy=True)
        out[~china_mask] = np.nan
        grid_cate = out
    else:
        china_geom = get_china_geometry()
        grid_cate = mask_grid_to_china(grid_lon, grid_lat, grid_cate, china_geom)

    # ── 对称发散色标 (以 0 为中心，仅用境内有效值) ─────────────────────────────
    valid = grid_cate[~np.isnan(grid_cate)]
    if len(valid) == 0:
        return None
    cate_lim = np.percentile(np.abs(valid), 98)
    if cate_lim < 1e-8:
        cate_lim = np.max(np.abs(valid)) + 1e-6

    # 可选：绘图降采样，减少 pcolormesh 像元数以加速
    plot_lon, plot_lat, plot_cate = grid_lon, grid_lat, grid_cate
    if display_res is not None and display_res > INTERP_RES:
        step = max(1, int(round(display_res / INTERP_RES)))
        plot_lon = grid_lon[::step]
        plot_lat = grid_lat[::step]
        plot_cate = grid_cate[::step, ::step]

    # Nature 风格发散配色: 蓝-白-红；境外（NaN）设为透明
    cmap = plt.cm.RdBu_r.copy()
    cmap.set_bad("none")
    norm = matplotlib.colors.TwoSlopeNorm(vmin=-cate_lim, vcenter=0, vmax=cate_lim)

    # ── 创建画布（除中国区外透明背景）────────────────────────────────────────
    fig = plt.figure(figsize=(10, 8), facecolor="none")
    fig.patch.set_edgecolor("none")
    fig.patch.set_linewidth(0)
    ax = fig.add_subplot(111, projection=ccrs.PlateCarree())

    # ── pcolormesh 渲染：仅中国区有颜色，境外透明 ─────────────────────────────
    mesh = ax.pcolormesh(
        plot_lon, plot_lat, plot_cate,
        transform=ccrs.PlateCarree(),
        cmap=cmap, norm=norm,
        shading="auto",
        rasterized=True,
        zorder=1,
    )

    # ── 叠加地理要素：仅中国区省界+南海断续线，不绘制境外（无海岸线/河流/湖泊），主图背景透明 ──
    setup_china_axes(ax,
                     add_coastlines=False, add_rivers=False, add_lakes=False, add_gridlines=False,
                     geoms_china=geoms_china, geoms_dash=geoms_dash, fast_plot=FAST_PLOT)

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
    if plot_hss_model is not None and plot_hss_threshold is not None:
        subtitle += f"  |  Mask: HSS-{plot_hss_model} ≥ {plot_hss_threshold:g}"
    # 副标题置于主标题上方，留出间距避免与主标题重叠
    ax.text(
        0.5, 1.03, subtitle,
        transform=ax.transAxes, ha="center", va="bottom",
        fontsize=9, color="grey", fontstyle="italic",
        fontfamily="sans-serif"
    )

    # 先收紧边距，再按最终主图位置放置小地图，避免小地图错位；top 留足给标题
    plt.subplots_adjust(left=0.04, right=0.96, top=0.90, bottom=0.14)
    if SHOW_SCS_INSET:
        fig.canvas.draw()
        # 南海小地图：必须在 subplots_adjust 之后添加，位置才贴合主图右下角
        add_scs_inset(fig, ax, plot_lon, plot_lat, plot_cate, None, cmap, norm,
                      geoms_china=geoms_china, geoms_dash=geoms_dash)
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
    if USE_HSS_MASK:
        print(f"  待绘制: {n_tasks} 组 (物种×变量)；每组按 HSS 模型各出一张（{HSS_MODELS}）\n")
    else:
        print(f"  待绘制: {n_tasks} 张 CATE 热力图（仅国界 mask）\n")

    # 保存任务清单
    combos.to_csv(os.path.join(TBL_DIR, "fig6_cate_plot_tasks.csv"), index=False)

    # ── 预计算：网格 + 中国 mask + 省界几何（全任务共用，显著加速）────────────
    print("预计算网格与边界（一次性）...")
    lon_min, lon_max = CHINA_EXTENT[0], CHINA_EXTENT[1]
    lat_min, lat_max = CHINA_EXTENT[2], CHINA_EXTENT[3]
    grid_lon = np.arange(lon_min, lon_max + INTERP_RES * 0.5, INTERP_RES)
    grid_lat = np.arange(lat_min, lat_max + INTERP_RES * 0.5, INTERP_RES)
    china_geom = get_china_geometry()
    if china_geom is None:
        print("  ⚠ 未加载中国省界几何 (CHINA_SHP)，无法按国界裁剪；请检查路径。")
    china_mask = build_china_mask(grid_lon, grid_lat, china_geom)
    geoms_china = get_boundary_geoms(CHINA_SHP)
    geoms_dash = get_boundary_geoms(DASH_SHP)
    print(f"  网格: {len(grid_lon)} x {len(grid_lat)}, 境内像元: {np.sum(china_mask):,}")

    hss_plot_masks = {}
    if USE_HSS_MASK:
        pred_dir = PRED_DIR if os.path.isabs(PRED_DIR) else os.path.join(BASE_DIR, PRED_DIR)
        sp_unique = combos["species"].drop_duplicates().tolist()
        hss_plot_masks = precompute_hss_plot_masks(
            sp_unique, grid_lon, grid_lat, china_mask, pred_dir,
            HSS_MODELS, HSS_CONF_THRESHOLD, HSS_INTERP_METHOD,
        )
        print(f"  [HSS] 已构建 (物种, 模型) 掩膜: {len(hss_plot_masks)} 个，阈值 ≥ {HSS_CONF_THRESHOLD}\n")
    else:
        print()

    fig_dir_abs = os.path.abspath(FIG_DIR)
    n_workers = max(0, int(N_WORKERS))

    def _save_cate_figure(fig, fname_base):
        fig.savefig(
            os.path.join(FIG_DIR, f"{fname_base}.png"),
            dpi=DPI, bbox_inches="tight", pad_inches=0.02, transparent=True,
        )
        try:
            fig.savefig(
                os.path.join(FIG_DIR, f"{fname_base}.svg"),
                bbox_inches="tight", pad_inches=0.02, transparent=True,
            )
        except Exception:
            pass
        plt.close(fig)

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
            if USE_HSS_MASK:
                for model in HSS_MODELS:
                    key = (sp, model)
                    if key not in hss_plot_masks:
                        print(f"  ⚠ 无 HSS 掩膜 ({sp}, {model})，跳过该模型")
                        continue
                    plot_mask = hss_plot_masks[key]
                    if not np.any(plot_mask):
                        print(f"  ⚠ HSS 掩膜全空 ({sp}, {model})，跳过")
                        continue
                    try:
                        fig = render_cate_map(
                            sub["lon"].values, sub["lat"].values, sub["cate"].values,
                            sp, var, ate,
                            grid_lon=grid_lon, grid_lat=grid_lat, china_mask=plot_mask,
                            geoms_china=geoms_china, geoms_dash=geoms_dash,
                            display_res=DISPLAY_RES,
                            plot_hss_model=model, plot_hss_threshold=HSS_CONF_THRESHOLD,
                        )
                        if fig is None:
                            print(f"  ⚠ [{model}] 插值后无有效数据, 跳过")
                            continue
                        fname = f"fig6_cate_{sp}_{var}_maskHSS_{model}"
                        _save_cate_figure(fig, fname)
                        success += 1
                        print(f"  ✓ 已保存 [{model}] → {fname}.png")
                    except Exception as e:
                        print(f"  ✗ [{model}] 错误: {e}")
                        plt.close("all")
            else:
                try:
                    fig = render_cate_map(
                        sub["lon"].values, sub["lat"].values, sub["cate"].values,
                        sp, var, ate,
                        grid_lon=grid_lon, grid_lat=grid_lat, china_mask=china_mask,
                        geoms_china=geoms_china, geoms_dash=geoms_dash,
                        display_res=DISPLAY_RES,
                    )
                    if fig is None:
                        print("  ⚠ 插值后无有效数据, 跳过")
                        continue
                    fname = f"fig6_cate_{sp}_{var}"
                    _save_cate_figure(fig, fname)
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
            if USE_HSS_MASK:
                for model in HSS_MODELS:
                    key = (sp, model)
                    if key not in hss_plot_masks:
                        continue
                    plot_mask = hss_plot_masks[key]
                    if not np.any(plot_mask):
                        continue
                    task_args.append((
                        sp, var, ate,
                        sub["lon"].values.copy(), sub["lat"].values.copy(), sub["cate"].values.copy(),
                        grid_lon, grid_lat, np.asarray(plot_mask, dtype=bool),
                        geoms_china, geoms_dash,
                        fig_dir_abs, DPI, model, float(HSS_CONF_THRESHOLD),
                    ))
            else:
                task_args.append((
                    sp, var, ate,
                    sub["lon"].values.copy(), sub["lat"].values.copy(), sub["cate"].values.copy(),
                    grid_lon, grid_lat, np.asarray(china_mask, dtype=bool),
                    geoms_china, geoms_dash,
                    fig_dir_abs, DPI, None, None,
                ))
        print(f"  使用 {n_workers} 个进程并行渲染 {len(task_args)} 张图...\n")
        with ProcessPoolExecutor(max_workers=n_workers) as executor:
            futures = {executor.submit(_render_one_task, a): a for a in task_args}
            for i, future in enumerate(as_completed(futures)):
                ok, sp, var, model_tag, err = future.result()
                if ok:
                    success += 1
                    tag = f" maskHSS_{model_tag}" if model_tag else ""
                    print(f"  ✓ [{i+1}/{len(task_args)}] {sp} x {var}{tag}")
                else:
                    print(f"  ✗ [{i+1}/{len(task_args)}] {sp} x {var}  {err or ''}")

    print(f"\n{'=' * 60}")
    print(f"  完成: {success} 张 CATE 热力图已保存至 {FIG_DIR}")
    print(f"{'=' * 60}")


def _render_one_task(args):
    """
    单张图渲染任务，供多进程调用。args: (sp, var, ate, lons, lats, cate_vals,
    grid_lon, grid_lat, plot_mask, geoms_china, geoms_dash, fig_dir, dpi,
    hss_model, hss_thr)。
    返回 (success, sp, var, hss_model or None, err_msg or None)。
    """
    (sp, var, ate, lons, lats, cate_vals, grid_lon, grid_lat, plot_mask,
     geoms_china, geoms_dash, fig_dir, dpi, hss_model, hss_thr) = args
    try:
        fig = render_cate_map(
            lons, lats, cate_vals, sp, var, ate,
            grid_lon=grid_lon, grid_lat=grid_lat, china_mask=plot_mask,
            geoms_china=geoms_china, geoms_dash=geoms_dash,
            display_res=DISPLAY_RES,
            plot_hss_model=hss_model, plot_hss_threshold=hss_thr,
        )
        if fig is None:
            return (False, sp, var, hss_model, "无有效数据")
        fname = (f"fig6_cate_{sp}_{var}_maskHSS_{hss_model}"
                 if hss_model else f"fig6_cate_{sp}_{var}")
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
        return (True, sp, var, hss_model, None)
    except Exception as e:
        plt.close("all")
        return (False, sp, var, hss_model, str(e))


if __name__ == "__main__":
    main()
