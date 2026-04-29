"""
plot_Lota_lota_distribution.py
绘制 Lota lota（江鳕）的全球分布点图
底图使用 refPackage/plot-function-main/data/countries.shp
"""
import sys
from pathlib import Path

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import cartopy.io.shapereader as shpreader
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# 路径配置
SCRIPT_DIR = Path(__file__).parent.resolve()
DATA_DIR   = SCRIPT_DIR                  # Lota_lota.csv 所在目录
REF_DIR    = SCRIPT_DIR.parent / "refPackage" / "plot-function-main" / "data"

# ----------------------------- 样式设置 -----------------------------
mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["axes.linewidth"] = 0.5
mpl.rcParams["axes.spines.bottom"] = True
mpl.rcParams["axes.spines.left"] = True
mpl.rcParams["axes.spines.top"] = True
mpl.rcParams["axes.spines.right"] = True

# 配色
LAND_COLOR  = "#F5F5F5"
COAST_COLOR = "#4A4A4A"
OCEAN_COLOR = "#EAF6FF"
POINT_COLOR = "#1A1A1A"
POINT_EDGE  = "#000000"


# ----------------------------- 读取数据 -----------------------------
def load_occurrence_data(csv_path: Path) -> pd.DataFrame:
    """读取物种分布CSV"""
    df = pd.read_csv(csv_path)
    required = {"taxon", "longitude", "latitude", "occurrenceStatus"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"CSV缺少必要列: {missing}")
    return df


def get_map_geoms(shp_path: Path, name_attr: str = None):
    """读取shapefile要素"""
    reader = shpreader.Reader(str(shp_path))
    if name_attr:
        geoms = [rec.geometry for rec in reader.records() if rec.attributes.get(name_attr)]
    else:
        geoms = [rec.geometry for rec in reader.records()]
    reader.close()
    return geoms


# ----------------------------- 核心绘图 -----------------------------
def plot_distribution_points(
    df: pd.DataFrame,
    output_path: Path,
    title: str = "Lota lota (Eurasian Burbot) — Global Distribution",
    figsize: tuple = (14, 7),
    point_size: float = 3,
    coastline_color: str = COAST_COLOR,
    land_color: str = LAND_COLOR,
    ocean_color: str = OCEAN_COLOR,
    point_color: str = POINT_COLOR,
    point_edge_color: str = POINT_EDGE,
    projection=ccrs.Robinson(),
    dpi: int = 150,
):
    """
    在全球底图上绘制物种分布点

    Parameters
    ----------
    df : pd.DataFrame
        含 longitude, latitude 列的DataFrame
    output_path : Path
        图片输出路径
    title : str
        图标题
    figsize : tuple
        图形尺寸（英寸）
    point_size : float
        散点大小
    coastline_color : str
        海岸线颜色
    land_color : str
        陆地颜色
    ocean_color : str
        海洋颜色
    point_color : str
        分布点颜色
    point_edge_color : str
        分布点描边颜色
    projection : ccrs.Projection
        地图投影
    dpi : int
        输出分辨率
    """
    # 筛选有效坐标
    df_valid = df.dropna(subset=["longitude", "latitude"])
    lon = df_valid["longitude"].values
    lat = df_valid["latitude"].values

    # 统计信息
    n_total = len(df)
    n_valid = len(df_valid)

    # ----------------------------- 创建画布 -----------------------------
    fig = plt.figure(figsize=figsize, facecolor="white")
    ax = fig.add_subplot(1, 1, 1, projection=projection)

    # 海洋背景
    ax.add_feature(
        cfeature.OCEAN,
        facecolor=OCEAN_COLOR,
        zorder=0,
    )

    # 陆地填充
    countries_shp = REF_DIR / "countries.shp"
    if countries_shp.exists():
        land_geoms = get_map_geoms(countries_shp)
        ax.add_geometries(
            land_geoms,
            crs=ccrs.PlateCarree(),
            facecolor=land_color,
            edgecolor=coastline_color,
            linewidth=0.3,
            zorder=1,
        )
    else:
        ax.add_feature(cfeature.LAND, facecolor=land_color, zorder=1)

    # 绘制分布点
    scatter = ax.scatter(
        lon,
        lat,
        c=point_color,
        s=point_size,
        alpha=0.75,
        edgecolors=point_edge_color,
        linewidths=0.15,
        transform=ccrs.PlateCarree(),
        zorder=5,
        label=f"Lota lota (n={n_valid})",
    )

    # 海岸线（叠加在最上层）
    ax.coastlines(
        color=coastline_color,
        linewidth=0.4,
        zorder=3,
    )

    # 国家边界（更细）
    ax.add_feature(
        cfeature.BORDERS,
        linestyle="-",
        edgecolor="#888888",
        linewidth=0.15,
        zorder=2,
    )

    ax.set_global()

    # ----------------------------- 标题与标注 -----------------------------
    fig.suptitle(
        title,
        fontsize=13,
        fontweight="bold",
        y=0.97,
        color="#1A1A1A",
    )

    ax.set_title(
        f"n = {n_valid:,}  records",
        fontsize=9,
        pad=6,
        color="#555555",
        style="italic",
    )

    # 图例
    legend = ax.legend(
        loc="lower left",
        fontsize=8.5,
        frameon=True,
        framealpha=0.85,
        edgecolor="#CCCCCC",
        fancybox=False,
        borderpad=0.4,
        labelspacing=0.3,
    )
    legend.get_frame().set_linewidth(0.5)

    # 来源标注
    fig.text(
        0.99, 0.01,
        "Data: GBIF / iNaturalist",
        ha="right",
        va="bottom",
        fontsize=7,
        color="#888888",
        style="italic",
    )

    # ----------------------------- 保存 -----------------------------
    fig.savefig(
        output_path,
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
        edgecolor="none",
    )
    plt.close(fig)
    print(f"[OK] 图片已保存至: {output_path}")
    print(f"     有效记录: {n_valid:,} / {n_total:,}")


# ----------------------------- 主程序 -----------------------------
if __name__ == "__main__":
    csv_file = DATA_DIR / "Lota_lota.csv"
    out_file = DATA_DIR / "Lota_lota_distribution.png"

    if not csv_file.exists():
        print(f"[ERROR] 数据文件不存在: {csv_file}")
        sys.exit(1)

    print(f"[INFO] 读取数据: {csv_file}")
    df = load_occurrence_data(csv_file)
    print(f"[INFO] 共读取 {len(df):,} 条记录")

    plot_distribution_points(df, out_file, dpi=150)
