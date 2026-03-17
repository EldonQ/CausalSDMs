#!/usr/bin/env python3
"""
Fig 1: CAST Framework Architecture
黑色边框 · Font Awesome 图标 · 结构可视化 · 最少文字
参考风格: eco01 (黑色边框) + MaskSDM 论文 (图标驱动)
输出: figures/case2_eco/plot/
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle
from matplotlib.font_manager import FontProperties
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.offsetbox import OffsetImage, AnnotationBbox
import numpy as np
import os
from io import BytesIO
import inspect
import urllib.parse
import urllib.request

# ─── pyconify 图标（优先） ───
# 说明：
# - 本脚本会优先尝试使用 pyconify 从 Iconify 拉取 PNG 图标并嵌入 matplotlib。
# - 若当前环境未安装 pyconify / Pillow 或无法联网，将自动回退到 Font Awesome 文本图标
#   或 eco01 风格的手绘图标（保证脚本仍可运行并出图）。
PYCONIFY_OK = False
PIL_OK = False
try:
    import pyconify  # type: ignore
    PYCONIFY_OK = True
    print('[pyconify] available')
except Exception:
    print('[pyconify] not installed — fallback to FA')

try:
    from PIL import Image  # type: ignore
    PIL_OK = True
except Exception:
    PIL_OK = False

try:
    # matplotlib 自带的 PNG reader（无 Pillow 时兜底）
    from matplotlib._png import read_png  # type: ignore
except Exception:
    read_png = None

# SVG -> PNG 转换（可选依赖）
CAIROSVG_OK = False
try:
    import cairosvg  # type: ignore
    CAIROSVG_OK = True
    print('[cairosvg] available (svg->png enabled)')
except Exception:
    CAIROSVG_OK = False
    if PYCONIFY_OK:
        print('[cairosvg] not installed — pyconify icons will fallback unless you install cairosvg')

# ─── 路径配置 ───
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..', '..'))
FIG_DIR = os.path.join(PROJECT_ROOT, 'figures', 'case2_eco', 'plot')
os.makedirs(FIG_DIR, exist_ok=True)

# ─── Font Awesome 6 加载 ───
FA = None
try:
    import fontawesomefree as _fam
    _fp = os.path.join(
        os.path.dirname(_fam.__file__),
        'static', 'fontawesomefree', 'otfs',
        'Font Awesome 6 Free-Solid-900.otf')
    if os.path.exists(_fp):
        FA = FontProperties(fname=_fp)
        print('[FA] Font Awesome 6 loaded')
    else:
        print(f'[FA] otf not found at {_fp}')
except ImportError:
    print('[FA] fontawesomefree not installed — using fallback')

# ─── 全局字体 ───
plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'figure.facecolor': 'white',
})

# ─── 颜色方案 (参考 eco01: 黑色边框 + 内部彩色) ───
BK = 'black'
GY = '#777777'
LG = '#CCCCCC'
WH = 'white'

C_BLUE   = '#4472C4'
C_RED    = '#C0504D'
C_GREEN  = '#548235'
C_ORANGE = '#ED7D31'
C_TEAL   = '#2B8C8C'
C_BROWN  = '#A0522D'


# ════════════════════════════════════════════════════════════════
# 基础绘图 (完全参考 eco01 的 box / arr)
# ════════════════════════════════════════════════════════════════

def box(ax, cx, cy, w, h, lw=1.0):
    """纯黑圆角边框白底矩形"""
    ax.add_patch(FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.06",
        facecolor=WH, edgecolor=BK, linewidth=lw, zorder=2))


def arr(ax, x1, y1, x2, y2, lw=1.0, color=BK, pad=0.12):
    """
    箭头（不与方框边框重叠）。

    pad: 端点内缩距离（数据坐标），会沿箭头方向从两端各缩短 pad。
         对阶段间水平箭头，0.10–0.15 通常足够避免压到黑色边框。
    """
    dx, dy = (x2 - x1), (y2 - y1)
    d = np.hypot(dx, dy)
    if d < 1e-9:
        return
    ux, uy = dx / d, dy / d
    xa, ya = x1 + ux * pad, y1 + uy * pad
    xb, yb = x2 - ux * pad, y2 - uy * pad
    ax.annotate(
        '',
        xy=(xb, yb),
        xytext=(xa, ya),
        arrowprops=dict(
            arrowstyle='->',
            color=color,
            lw=lw,
            shrinkA=0,
            shrinkB=0,
            mutation_scale=10,
        ),
        zorder=10
    )


def fa_icon(ax, x, y, code, size=14, color=BK):
    """渲染 Font Awesome 图标，返回是否成功"""
    if FA is not None:
        ax.text(x, y, code, fontproperties=FA,
                fontsize=size, color=color,
                ha='center', va='center', zorder=8)
        return True
    return False


_ICON_CACHE = {}


def _pyconify_get_png_bytes(icon: str, color: str = None, size: int = 64):
    """
    从 pyconify 拉取 PNG 字节。

    兼容不同 pyconify 版本的接口：优先尝试 `pyconify.png(...)`，其次尝试
    `pyconify.get(...)` / `pyconify.icon(...)` 等常见命名。
    """
    key = (icon, color, size)
    if key in _ICON_CACHE:
        return _ICON_CACHE[key]

    if not PYCONIFY_OK:
        return None

    try:
        # pyconify 在当前环境提供的是 svg(key=...) -> bytes
        if not hasattr(pyconify, "svg"):
            return None
        # 注意：你当前环境的签名为 svg(*key, color=..., height=..., width=...)，key 需位置参数传入
        svg_bytes = pyconify.svg(icon, color=color, height=size)
        if not svg_bytes:
            return None
        if not CAIROSVG_OK:
            return None
        png_bytes = cairosvg.svg2png(bytestring=svg_bytes, output_width=size, output_height=size)
        if png_bytes:
            _ICON_CACHE[key] = png_bytes
            return png_bytes
    except Exception as e:
        print(f"[pyconify-svg] failed to fetch/convert {icon}: {e}")
        return None


def iconify_icon(ax, x, y, icon: str, size_px: int = 26, color: str = BK, alpha: float = 1.0):
    """
    用 pyconify 的 Iconify 图标渲染到 matplotlib。
    成功返回 True；失败返回 False（调用方可回退到 FA 或手绘）。
    """
    if not PYCONIFY_OK:
        return False

    # 拉取稍大一点的 PNG，避免缩放锯齿
    png_bytes = _pyconify_get_png_bytes(icon=icon, color=color, size=max(64, int(size_px * 3)))
    if not png_bytes:
        return False

    try:
        if PIL_OK:
            im = Image.open(BytesIO(png_bytes)).convert("RGBA")
            rgba = np.asarray(im)
            base = max(im.size)
        else:
            # 无 Pillow：用 matplotlib._png.read_png 读取
            if read_png is None:
                return False
            rgba = read_png(BytesIO(png_bytes))
            base = max(rgba.shape[0], rgba.shape[1])

        # matplotlib 的 zoom 是相对原图像素；用 size_px 控制显示大小
        zoom = size_px / float(base)
        oi = OffsetImage(rgba, zoom=zoom, resample=True)
        oi.set_alpha(alpha)
        ab = AnnotationBbox(oi, (x, y), frameon=False, box_alignment=(0.5, 0.5), zorder=9)
        ax.add_artist(ab)
        return True
    except Exception as e:
        print(f"[pyconify] render failed for {icon}: {e}")
        return False


# ════════════════════════════════════════════════════════════════
# 方法可视化 (自定义绘制)
# ════════════════════════════════════════════════════════════════

def draw_dag(ax, cx, cy, s=1.0):
    """DAG 因果结构示意"""
    nd = {0: (cx - 0.72 * s, cy + 0.42 * s),
          1: (cx - 0.72 * s, cy - 0.42 * s),
          2: (cx,             cy + 0.15 * s),
          3: (cx,             cy - 0.28 * s),
          4: (cx + 0.72 * s,  cy)}
    nc = [C_BLUE, C_TEAL, C_ORANGE, C_GREEN, C_RED]
    edges = [(0, 2), (1, 2), (1, 3), (2, 4), (3, 4)]
    r = 0.13 * s
    for a, b in edges:
        xa, ya = nd[a]; xb, yb = nd[b]
        dx, dy = xb - xa, yb - ya
        d = np.hypot(dx, dy)
        off = r + 0.02 * s
        ax.annotate(
            '', xy=(xb - dx / d * off, yb - dy / d * off),
            xytext=(xa + dx / d * off, ya + dy / d * off),
            arrowprops=dict(arrowstyle='->', color='#555555',
                            lw=0.7 * s, shrinkA=0, shrinkB=0,
                            connectionstyle='arc3,rad=0.06'),
            zorder=6)
    for i, (x, y) in nd.items():
        ax.add_patch(Circle((x, y), r,
                            fc=nc[i], ec=WH, lw=0.5, zorder=7))


def draw_forest(ax, cx, cy, s=1.0):
    """森林图 — ATE 效应量示意"""
    effs = [0.35, -0.20, 0.50, 0.08, -0.38]
    ax.plot([cx, cx], [cy - 0.52 * s, cy + 0.52 * s],
            color=LG, lw=0.5, ls='--', zorder=3)
    for i, e in enumerate(effs):
        yy = cy + (2 - i) * 0.21 * s
        ci = abs(e) * 0.38 + 0.05
        col = C_RED if e > 0 else C_BLUE
        ax.plot([cx + (e - ci) * s, cx + (e + ci) * s], [yy, yy],
                color=col, lw=0.65, zorder=4, solid_capstyle='round')
        ax.plot(cx + e * s, yy, 'D', color=col, ms=2.6 * s, zorder=5)


def draw_screen_feat(ax, cx, cy, s=1.0):
    """自适应筛选 + 特征工程示意"""
    # ── 上半: 三准则评分条 ──
    top_y = cy + 0.4 * s
    bxs = [cx - 0.5 * s, cx, cx + 0.5 * s]
    bcl = [C_BLUE, C_RED, C_GREEN]
    bht = [0.38, 0.26, 0.32]
    bw = 0.28 * s
    for bx, col, bh in zip(bxs, bcl, bht):
        ax.add_patch(Rectangle(
            (bx - bw / 2, top_y), bw, bh * s * 0.55,
            fc=col, ec='none', zorder=5, alpha=0.75))
    # 汇聚箭头
    for bx in bxs:
        arr(ax, bx, top_y - 0.03 * s, cx, top_y - 0.22 * s,
            lw=0.35, color=GY, pad=0.05 * s)
    # 合成条
    comp_y = top_y - 0.38 * s
    ax.add_patch(Rectangle(
        (cx - 0.32 * s, comp_y), 0.64 * s, 0.12 * s,
        fc=C_ORANGE, ec=BK, lw=0.3, zorder=5, alpha=0.7))

    # ── 连接箭头 ──
    arr(ax, cx, comp_y - 0.02 * s,
        cx, cy - 0.22 * s, lw=0.4, color=GY, pad=0.05 * s)

    # ── 下半: 特征矩阵增强 ──
    my = cy - 0.5 * s
    rows, co, cn = 3, 2, 4
    cw, ch = 0.09 * s, 0.07 * s
    ox = cx - 0.55 * s
    for i in range(rows):
        for j in range(co):
            ax.add_patch(Rectangle(
                (ox + j * cw, my + i * ch), cw * 0.88, ch * 0.85,
                fc=C_BLUE, ec=WH, lw=0.08, zorder=5, alpha=0.55))
    arr(ax, ox + co * cw + 0.05 * s, my + rows * ch / 2,
        cx + 0.08 * s, my + rows * ch / 2,
        lw=0.4, color=GY, pad=0.05 * s)
    nx = cx + 0.14 * s
    fcl = [C_BLUE, C_BLUE, C_ORANGE, C_GREEN]
    np.random.seed(7)
    for i in range(rows):
        for j in range(cn):
            ax.add_patch(Rectangle(
                (nx + j * cw, my + i * ch), cw * 0.88, ch * 0.85,
                fc=fcl[j], ec=WH, lw=0.08, zorder=5,
                alpha=0.4 + 0.35 * np.random.random()))


def draw_nn(ax, cx, cy, s=1.0):
    """CI-MLP 神经网络示意"""
    layers = [3, 5, 5, 3, 1]
    lx = [cx + (i - 2) * 0.34 * s for i in range(5)]
    pos = []
    for x, n in zip(lx, layers):
        pos.append([(x, cy + (j - (n - 1) / 2) * 0.14 * s)
                    for j in range(n)])
    for li in range(len(layers) - 1):
        for x1, y1 in pos[li]:
            for x2, y2 in pos[li + 1]:
                ax.plot([x1, x2], [y1, y2],
                        color='#D0D8E8', lw=0.12, zorder=3)
    for ps in pos:
        for x, y in ps:
            ax.add_patch(Circle((x, y), 0.04 * s,
                                fc=C_TEAL, ec=WH, lw=0.12, zorder=5))


def _mini_method_box(ax, x, y, w, h, title, icon_key=None, icon_color=GY, title_color=BK):
    """CI-MLP 框内的小模块（方法图子块）。"""
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h,
        boxstyle="round,pad=0.05",
        facecolor="white",
        edgecolor="#222222",
        linewidth=0.7,
        zorder=4
    ))
    # 图标（能渲染则用 pyconify，否则留空不破坏布局）
    if icon_key is not None:
        iconify_icon(ax, x + 0.18, y + h - 0.18, icon_key, size_px=18, color=icon_color, alpha=0.95)
    ax.text(x + w / 2, y + h - 0.22, title,
            ha="center", va="center", fontsize=7.6, color=title_color, fontweight="bold", zorder=6)


def draw_ci_mlp_panel(ax, cx, cy, s=1.0):
    """
    CI-MLP（偏方法图）：
    - 上层：ATE 加权输入 + DAG 引导交互特征构造
    - 下层：更干净的 MLP 学习器示意
    """
    # 坐标与尺度（以原 draw_nn 的中心点为参考）
    top_band_y = cy + 0.58 * s
    top_band_h = 0.72 * s
    left_x = cx - 0.95 * s
    total_w = 1.90 * s

    # ── 上层：两个方法子块 ──────────────────────────────────────────
    gap = 0.10 * s
    mb_w = (total_w - gap) / 2
    mb_h = top_band_h
    mb_y = top_band_y - mb_h / 2

    # ATE weighting
    x1 = left_x
    _mini_method_box(
        ax, x1, mb_y, mb_w, mb_h,
        title="ATE weighting",
        icon_key="mdi:scale-balance",
        icon_color=C_RED
    )
    ax.text(x1 + mb_w / 2, mb_y + 0.32 * s,
            r"$\tilde{x}_j = w_j x_j$",
            ha="center", va="center", fontsize=7.2, color=BK, zorder=7)
    bx = x1 + 0.18
    by = mb_y + 0.12 * s
    for k, a in enumerate([0.55, 0.35, 0.22]):
        yy = by + k * 0.10 * s
        ax.add_patch(Rectangle((bx, yy), 0.22, 0.03 * s, fc="#D9D9D9", ec="none", zorder=6))
        ax.add_patch(Rectangle((bx + 0.30, yy), 0.22 * (0.6 + a), 0.03 * s,
                               fc=C_RED, ec="none", alpha=0.55, zorder=6))
        arr(ax, bx + 0.24, yy + 0.015 * s, bx + 0.30, yy + 0.015 * s,
            lw=0.6, color="#999999", pad=0.02)

    # DAG-guided interactions
    x2 = left_x + mb_w + gap
    _mini_method_box(
        ax, x2, mb_y, mb_w, mb_h,
        title="DAG interactions",
        icon_key="mdi:graph-outline",
        icon_color=C_BLUE
    )
    ax.text(x2 + mb_w / 2, mb_y + 0.32 * s,
            r"$x_i x_j\;(i\!\to\!j\in DAG)$",
            ha="center", va="center", fontsize=7.0, color=BK, zorder=7)
    vby = mb_y + 0.12 * s
    ax.add_patch(Rectangle((x2 + 0.16, vby + 0.06 * s), 0.20, 0.03 * s,
                           fc=C_BLUE, ec="none", alpha=0.55, zorder=6))
    ax.add_patch(Rectangle((x2 + 0.16, vby - 0.02 * s), 0.20, 0.03 * s,
                           fc=C_GREEN, ec="none", alpha=0.55, zorder=6))
    ax.text(x2 + 0.44, vby + 0.02 * s, "×", ha="center", va="center",
            fontsize=9, color="#666666", zorder=7)
    ax.add_patch(Rectangle((x2 + 0.54, vby + 0.01 * s), 0.24, 0.06 * s,
                           fc=C_ORANGE, ec=BK, lw=0.25, alpha=0.55, zorder=6))
    ax.add_patch(Circle((x2 + 0.20, vby + 0.22 * s), 0.015 * s,
                        fc=C_BLUE, ec=WH, lw=0.2, zorder=7))
    ax.add_patch(Circle((x2 + 0.32, vby + 0.16 * s), 0.015 * s,
                        fc=C_ORANGE, ec=WH, lw=0.2, zorder=7))
    arr(ax, x2 + 0.20, vby + 0.22 * s, x2 + 0.32, vby + 0.16 * s,
        lw=0.6, color="#777777", pad=0.02)

    # ── 下层：MLP（更干净的连线）────────────────────────────────────
    mlp_cy = cy - 0.12 * s
    layers2 = [4, 6, 5, 3, 1]
    lx2 = [cx + (i - 2) * 0.30 * s for i in range(len(layers2))]
    pos2 = []
    for x, n in zip(lx2, layers2):
        pos2.append([(x, mlp_cy + (j - (n - 1) / 2) * 0.12 * s) for j in range(n)])
    for li in range(len(layers2) - 1):
        left = pos2[li]
        right = pos2[li + 1]
        for i_idx, (x1p, y1p) in enumerate(left):
            for j_idx, (x2p, y2p) in enumerate(right):
                if (i_idx + j_idx) % 2 == 0:
                    ax.plot([x1p, x2p], [y1p, y2p], color="#D6DFEF", lw=0.10, zorder=3)
    for ps in pos2:
        for x, y in ps:
            ax.add_patch(Circle((x, y), 0.035 * s, fc=C_TEAL, ec=WH, lw=0.12, zorder=5))

    # 上层 → MLP 的汇聚箭头
    arr(ax, cx - 0.40 * s, mb_y - 0.04 * s, cx - 0.25 * s, mlp_cy + 0.18 * s,
        lw=0.7, color="#666666", pad=0.04)
    arr(ax, cx + 0.40 * s, mb_y - 0.04 * s, cx + 0.25 * s, mlp_cy + 0.18 * s,
        lw=0.7, color="#666666", pad=0.04)


def draw_heatmap(ax, cx, cy, s=1.0):
    """物种分布预测热力图"""
    cmap = LinearSegmentedColormap.from_list(
        'pred', [C_BLUE, '#F0F0F0', C_RED], N=64)
    n = 5
    cell = 0.1 * s
    np.random.seed(42)
    v = np.random.random((n, n))
    for i in range(n):
        for j in range(n):
            ax.add_patch(Rectangle(
                (cx - n * cell / 2 + j * cell,
                 cy - n * cell / 2 + i * cell),
                cell * 0.92, cell * 0.92,
                fc=cmap(v[i, j]), ec=WH, lw=0.1, zorder=5))


# ════════════════════════════════════════════════════════════════
# 输入数据图标 (FA 优先，fallback 到 eco01 风格)
# ════════════════════════════════════════════════════════════════

def draw_input_column(ax, cx, cy_top, s=1.0):
    """左侧输入变量图标列"""
    items = [
        # (FA_code, label, color, fallback, iconify_id)
        ('\uf0c2', 'Climate',    C_BLUE,   'wave',     'mdi:weather-partly-cloudy'),
        ('\uf6fc', 'Topography', C_BROWN,  'mountain', 'mdi:terrain'),
        ('\uf043', 'Hydrology',  C_ORANGE, 'lines',    'mdi:water'),
        ('\uf4d8', 'Land cover', C_GREEN,  'tree',     'mdi:pine-tree'),
        ('\uf3c5', 'Occurrence', C_RED,    'dot',      'mdi:map-marker'),
    ]
    gap = 0.48 * s
    for i, (fa_code, name, col, fallback, iconify_id) in enumerate(items):
        y = cy_top - i * gap
        drawn = iconify_icon(ax, cx, y, iconify_id, size_px=int(22 * s), color=col, alpha=0.98)
        if not drawn:
            drawn = fa_icon(ax, cx, y, fa_code, size=11, color=col)
        if not drawn:
            # eco01 风格回退
            if fallback == 'wave':
                t = np.linspace(cx - 0.16 * s, cx + 0.16 * s, 25)
                ax.plot(t, y + np.sin(np.linspace(0, 2.5 * np.pi, 25)) * 0.05 * s,
                        color=col, lw=1.2, zorder=5)
            elif fallback == 'mountain':
                px = [cx - 0.17 * s, cx - 0.04 * s, cx + 0.04 * s,
                      cx + 0.13 * s, cx + 0.17 * s]
                py = [y - 0.05 * s, y + 0.07 * s, y,
                      y + 0.05 * s, y - 0.05 * s]
                ax.plot(px, py, color=col, lw=1.2, zorder=5)
            elif fallback == 'lines':
                for dy in [0.03, 0, -0.03]:
                    ax.plot([cx - 0.12 * s, cx + 0.12 * s],
                            [y + dy * s] * 2,
                            color=col, lw=1.0, zorder=5)
            elif fallback == 'tree':
                ax.plot([cx, cx], [y - 0.05 * s, y + 0.01 * s],
                        color=col, lw=1.2, zorder=5)
                ax.add_patch(Circle((cx, y + 0.05 * s), 0.05 * s,
                                    fc=col, ec=WH, lw=0.3, zorder=6))
            elif fallback == 'dot':
                ax.plot(cx, y, 'o', color=col, ms=5, zorder=5)
        ax.text(cx + 0.26 * s, y, name,
                ha='left', va='center', fontsize=7.5, color=BK)


# ════════════════════════════════════════════════════════════════
# 主函数
# ════════════════════════════════════════════════════════════════

def main():
    fig, ax = plt.subplots(1, 1, figsize=(15, 6), dpi=300)
    ax.set_xlim(0, 15)
    ax.set_ylim(0, 6)
    ax.axis('off')

    # ── 标题 ──
    ax.text(7.5, 5.65, 'CAST',
            ha='center', fontsize=20, fontweight='bold', color=BK)
    ax.text(7.5, 5.3,
            'Causally-Adaptive Species distribution modeling Technique',
            ha='center', fontsize=10, color=GY)

    # ── 布局常量 ──
    MY  = 2.65              # 中线 (箭头 y)
    BH  = 3.55              # 框体高度
    BT  = MY + BH / 2       # 4.425
    BB  = MY - BH / 2       # 0.875

    # 阶段框 (cx, width)
    S = [(2.8, 2.3), (5.55, 2.3), (8.65, 2.85), (11.65, 2.3)]

    stg_names    = ['Structure\nDiscovery',
                    'Effect\nEstimation',
                    'Screening &\nEngineering',
                    'CI-MLP']
    method_names = ['Bootstrap HC',
                    'Double ML',
                    'Multi-criteria',
                    'Focal · Ensemble']
    output_names = ['Consensus DAG',
                    'ATE · CATE',
                    'Causal features',
                    'SDM Prediction']
    # 右上角装饰图标（pyconify 优先，fallback 到 FA）
    fa_corner = ['\uf0e8', '\uf24e', '\uf0b0', '\uf5dc']
    iconify_corner = ['mdi:graph-outline', 'mdi:chart-line', 'mdi:filter-variant', 'mdi:brain']

    # ── 绘制4个阶段框 (黑色边框) ──
    for i, ((cx, w), sn, mn, on) in enumerate(
            zip(S, stg_names, method_names, output_names)):
        box(ax, cx, MY, w, BH, lw=1.0)

        # 右上角装饰图标（浅灰）
        ok = iconify_icon(
            ax,
            cx + w / 2 - 0.28,
            BT - 0.28,
            iconify_corner[i],
            size_px=19,
            color=LG,
            alpha=0.9
        )
        if not ok:
            fa_icon(ax, cx + w / 2 - 0.28, BT - 0.28,
                    fa_corner[i], size=10, color=LG)

        # 阶段名称
        ax.text(cx, BT - 0.42, sn,
                ha='center', va='center', fontsize=10,
                fontweight='bold', color=BK, linespacing=1.15)

        # 方法名称
        ax.text(cx, BT - 1.0, mn,
                ha='center', fontsize=8, color=GY)

        # 输出标签
        ax.text(cx, BB + 0.2, on,
                ha='center', fontsize=7.5, color=BK)

    # ── 方法可视化 (各阶段核心图标) ──
    iy = MY - 0.18
    draw_dag(ax, S[0][0], iy, s=0.88)
    draw_forest(ax, S[1][0], iy, s=0.85)
    draw_screen_feat(ax, S[2][0], iy, s=0.88)
    draw_ci_mlp_panel(ax, S[3][0], iy, s=1.0)

    # ── 输入数据列 ──
    inp_x = 0.55
    draw_input_column(ax, inp_x, 3.55, s=0.85)
    ax.text(inp_x, 1.15, 'VIF-screened',
            ha='center', fontsize=7, color=GY, style='italic')

    # 输入 → Stage 1
    arr(ax, inp_x + 0.62, MY, S[0][0] - S[0][1] / 2, MY, lw=1.0, pad=0.14)

    # ── 阶段间箭头 + 数据流标注 ──
    flow_labels = ['DAG', 'ATE', 'Features']
    for i in range(3):
        x1 = S[i][0] + S[i][1] / 2
        x2 = S[i + 1][0] - S[i + 1][1] / 2
        arr(ax, x1, MY, x2, MY, lw=1.0, pad=0.14)
        ax.text((x1 + x2) / 2, MY + 0.22, flow_labels[i],
                ha='center', fontsize=6.5, color=GY, style='italic')

    # ── 输出 (右侧: SDM预测 + 空间CATE + 模型评估) ──
    out_x = 13.9
    arr(ax, S[3][0] + S[3][1] / 2, MY, out_x - 0.55, MY, lw=1.0, pad=0.14)

    # SDM 预测热力图
    draw_heatmap(ax, out_x, MY + 0.55, s=0.78)
    ax.text(out_x, MY + 0.05, 'SDM Prediction',
            ha='center', fontsize=7.5, color=BK)

    # CATE 空间图
    draw_heatmap(ax, out_x, MY - 0.6, s=0.6)
    ax.text(out_x, MY - 1.05, 'Spatial CATE',
            ha='center', fontsize=7, color=GY)

    # 评估指标
    ax.text(out_x, MY - 1.42, 'AUC · TSS',
            ha='center', fontsize=7, color=GY, style='italic')

    # ── 保存：仅 PNG（按你的要求不输出 SVG） ──
    fp = os.path.join(FIG_DIR, 'fig1_cast_framework_architecture_pyconify.png')
    fig.savefig(
        fp,
        dpi=1200,
        format='png',
        bbox_inches='tight',
        facecolor='white',
        edgecolor='none'
    )
    print(f'[OK] {os.path.basename(fp)}')

    plt.close()


if __name__ == '__main__':
    main()
