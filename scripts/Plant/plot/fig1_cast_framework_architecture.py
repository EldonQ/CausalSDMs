#!/usr/bin/env python3
"""
Fig 1: CAST Framework Architecture
黑色边框 · Font Awesome 图标 · 结构可视化 · 最少文字
参考风格: eco01 (黑色边框) + MaskSDM 论文 (图标驱动)
输出: figures/case4_plant/plot/
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle
from matplotlib.font_manager import FontProperties
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import os

# ─── 路径配置 ───
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..', '..'))
FIG_DIR = os.path.join(PROJECT_ROOT, 'figures', 'case4_plant', 'plot')
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


def arr(ax, x1, y1, x2, y2, lw=1.0, color=BK):
    """箭头 (参考 eco01)"""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                shrinkA=1, shrinkB=1), zorder=10)


def fa_icon(ax, x, y, code, size=14, color=BK):
    """渲染 Font Awesome 图标，返回是否成功"""
    if FA is not None:
        ax.text(x, y, code, fontproperties=FA,
                fontsize=size, color=color,
                ha='center', va='center', zorder=8)
        return True
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
            lw=0.35, color=GY)
    # 合成条
    comp_y = top_y - 0.38 * s
    ax.add_patch(Rectangle(
        (cx - 0.32 * s, comp_y), 0.64 * s, 0.12 * s,
        fc=C_ORANGE, ec=BK, lw=0.3, zorder=5, alpha=0.7))

    # ── 连接箭头 ──
    arr(ax, cx, comp_y - 0.02 * s,
        cx, cy - 0.22 * s, lw=0.4, color=GY)

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
        lw=0.4, color=GY)
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
        ('\uf0c2', 'Climate',    C_BLUE,   'wave'),
        ('\uf6fc', 'Topography', C_BROWN,  'mountain'),
        ('\uf043', 'Hydrology',  C_ORANGE, 'lines'),
        ('\uf4d8', 'Land cover', C_GREEN,  'tree'),
        ('\uf3c5', 'Occurrence', C_RED,    'dot'),
    ]
    gap = 0.48 * s
    for i, (fa_code, name, col, fallback) in enumerate(items):
        y = cy_top - i * gap
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
    # 右上角装饰 FA 图标
    fa_corner = ['\uf0e8', '\uf24e', '\uf0b0', '\uf5dc']

    # ── 绘制4个阶段框 (黑色边框) ──
    for i, ((cx, w), sn, mn, on) in enumerate(
            zip(S, stg_names, method_names, output_names)):
        box(ax, cx, MY, w, BH, lw=1.0)

        # FA 装饰图标 (右上角，浅灰)
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
    draw_nn(ax, S[3][0], iy, s=1.0)

    # ── 输入数据列 ──
    inp_x = 0.55
    draw_input_column(ax, inp_x, 3.55, s=0.85)
    ax.text(inp_x, 1.15, 'VIF-screened',
            ha='center', fontsize=7, color=GY, style='italic')

    # 输入 → Stage 1
    arr(ax, inp_x + 0.62, MY, S[0][0] - S[0][1] / 2, MY, lw=1.0)

    # ── 阶段间箭头 + 数据流标注 ──
    flow_labels = ['DAG', 'ATE', 'Features']
    for i in range(3):
        x1 = S[i][0] + S[i][1] / 2
        x2 = S[i + 1][0] - S[i + 1][1] / 2
        arr(ax, x1, MY, x2, MY, lw=1.0)
        ax.text((x1 + x2) / 2, MY + 0.22, flow_labels[i],
                ha='center', fontsize=6.5, color=GY, style='italic')

    # ── 输出 (右侧: SDM预测 + 空间CATE + 模型评估) ──
    out_x = 13.9
    arr(ax, S[3][0] + S[3][1] / 2, MY, out_x - 0.55, MY, lw=1.0)

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

    # ── 保存 ──
    for ext in ['png', 'svg']:
        fp = os.path.join(FIG_DIR,
                          f'fig1_cast_framework_architecture.{ext}')
        dpi_val = 1200 if ext == 'png' else 150
        fig.savefig(fp, dpi=dpi_val, format=ext,
                    bbox_inches='tight',
                    facecolor='white', edgecolor='none')
        print(f'[OK] {os.path.basename(fp)}')

    plt.close()


if __name__ == '__main__':
    main()
