#!/usr/bin/env python3
"""
CAST核心方法技术路线图
黑色文字 · 黑色边框 · 无填充 · 图标/符号彩色 · 无Output框
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import os

plt.rcParams.update({
    'font.family': 'Arial',
    'font.size': 10,
    'figure.facecolor': 'white',
})

BK = 'black'
GY = '#777777'
LG = '#CCCCCC'

C_BLUE   = '#4472C4'
C_RED    = '#C0504D'
C_GREEN  = '#548235'
C_BROWN  = '#A0522D'
C_ORANGE = '#ED7D31'
C_TEAL   = '#2B8C8C'
C_PURPLE = '#7B68AD'

DIV_CMAP = LinearSegmentedColormap.from_list(
    'div', [C_BLUE, '#F0F0F0', C_RED], N=256)


def box(ax, cx, cy, w, h, lw=1.0):
    """白底黑框圆角矩形"""
    ax.add_patch(FancyBboxPatch(
        (cx - w / 2, cy - h / 2), w, h,
        boxstyle="round,pad=0.06",
        facecolor='white', edgecolor=BK, linewidth=lw, zorder=2))


def arr(ax, x1, y1, x2, y2, lw=1.0, color=BK):
    """箭头"""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color=color, lw=lw,
                                shrinkA=1, shrinkB=1), zorder=10)


# ────────────────────────── DAG (彩色节点) ──────────────────────────
def draw_dag(ax, cx, cy, s=1.0):
    nodes = {
        1: (cx - 1.05 * s, cy + 0.48 * s),
        2: (cx - 1.05 * s, cy - 0.48 * s),
        3: (cx,            cy + 0.22 * s),
        4: (cx,            cy - 0.30 * s),
        5: (cx + 1.05 * s, cy),
    }
    ncolors = {1: C_BLUE, 2: C_TEAL, 3: C_ORANGE, 4: C_GREEN, 5: C_RED}
    for a, b in [(1, 3), (2, 3), (2, 4), (3, 5), (4, 5)]:
        xa, ya = nodes[a]; xb, yb = nodes[b]
        ax.annotate('', xy=(xb, yb), xytext=(xa, ya),
                    arrowprops=dict(arrowstyle='->', color='#555555', lw=0.8,
                                    shrinkA=5.5 * s, shrinkB=5.5 * s,
                                    connectionstyle='arc3,rad=0.06'), zorder=6)
    for idx, (x, y) in nodes.items():
        ax.add_patch(Circle((x, y), 0.16 * s,
                            fc=ncolors[idx], ec='white', lw=0.6, zorder=7))


# ────────────────────────── ATE 森林图 (红正蓝负) ──────────────────
def draw_forest(ax, cx, cy, s=1.0):
    effs = [0.35, -0.20, 0.50, 0.08, -0.38]
    ax.plot([cx, cx], [cy - 0.65 * s, cy + 0.65 * s],
            color=LG, lw=0.5, ls='--', zorder=3)
    for i, e in enumerate(effs):
        yy = cy + (2 - i) * 0.26 * s
        ci = abs(e) * 0.5
        col = C_RED if e > 0 else C_BLUE
        ax.plot([cx + (e - ci) * s, cx + (e + ci) * s], [yy, yy],
                color=col, lw=0.7, zorder=4, solid_capstyle='round')
        ax.plot(cx + e * s, yy, 'o', color=col, ms=2.8 * s, zorder=5)


# ────────────────────────── CATE 热力图 (蓝白红) ──────────────────
def draw_heatmap(ax, cx, cy, s=1.0):
    grid = np.array([
        [ 0.7,  0.9,  0.5,  0.1, -0.2],
        [ 0.3,  0.8,  1.0,  0.6,  0.0],
        [ 0.0,  0.3,  0.7,  0.9,  0.5],
        [-0.2,  0.0,  0.4,  0.6,  0.8]])
    cw, ch = 0.17 * s, 0.13 * s
    for i in range(4):
        for j in range(5):
            v = (grid[i, j] + 0.4) / 1.4
            rgba = DIV_CMAP(v)
            ax.add_patch(plt.Rectangle(
                (cx - 2.5 * cw + j * cw, cy - 2 * ch + i * ch),
                cw, ch, fc=rgba, ec='white', lw=0.25, zorder=4))


# ────────────────────────── 神经网络 (蓝色节点) ──────────────────
def draw_nn(ax, cx, cy, s=1.0):
    layers = [4, 6, 6, 4, 1]
    lx = [cx + (i - 2) * 0.44 * s for i in range(5)]
    positions = []
    for x, n in zip(lx, layers):
        positions.append([(x, cy + (j - (n - 1) / 2) * 0.16 * s) for j in range(n)])
    for li in range(len(layers) - 1):
        for x1, y1 in positions[li]:
            for x2, y2 in positions[li + 1]:
                ax.plot([x1, x2], [y1, y2], color='#D0D8E8', lw=0.14, zorder=3)
    for ps in positions:
        for x, y in ps:
            ax.add_patch(Circle((x, y), 0.046 * s,
                                fc=C_BLUE, ec='white', lw=0.15, zorder=5))


# ────────────────────────── 输入图标 (各类别彩色) ──────────────────
def draw_icons(ax, cx, cy, s=1.0):
    gap = 1.5 * s
    xs = [cx - 1.5 * gap, cx - 0.5 * gap, cx + 0.5 * gap, cx + 1.5 * gap]
    names  = ['Climate', 'Topography', 'Soil', 'Land cover']
    colors = [C_BLUE,     C_BROWN,      C_ORANGE, C_GREEN]
    for x, name, col in zip(xs, names, colors):
        if name == 'Climate':
            t = np.linspace(x - 0.24 * s, x + 0.24 * s, 30)
            ax.plot(t, cy + np.sin(np.linspace(0, 2.5 * np.pi, 30)) * 0.08 * s,
                    color=col, lw=1.2, zorder=5)
        elif name == 'Topography':
            px = [x - 0.26 * s, x - 0.06 * s, x + 0.06 * s, x + 0.19 * s, x + 0.26 * s]
            py = [cy - 0.09 * s, cy + 0.11 * s, cy - 0.01 * s, cy + 0.08 * s, cy - 0.09 * s]
            ax.plot(px, py, color=col, lw=1.2, zorder=5)
        elif name == 'Soil':
            for dy in [0.05, 0, -0.05]:
                ax.plot([x - 0.2 * s, x + 0.2 * s], [cy + dy * s] * 2,
                        color=col, lw=1.2, zorder=5)
        elif name == 'Land cover':
            ax.plot([x, x], [cy - 0.09 * s, cy + 0.01 * s],
                    color=col, lw=1.2, zorder=5)
            ax.add_patch(Circle((x, cy + 0.07 * s), 0.08 * s,
                                fc=col, ec='white', lw=0.4, zorder=6))
        ax.text(x, cy - 0.20 * s, name, ha='center', va='top',
                fontsize=7, color=BK)


# ════════════════════════════════════════════════════════════════
# 主函数
# ════════════════════════════════════════════════════════════════
def main():
    os.makedirs('output/case2', exist_ok=True)

    fig, ax = plt.subplots(1, 1, figsize=(10, 12), dpi=300)
    ax.set_xlim(0, 10)
    ax.set_ylim(3.2, 15.7)
    ax.axis('off')

    CX = 5.0
    W = 7.2
    GAP = 0.18
    LX = CX - W / 2

    iy  = 14.8;  ih  = 1.25
    s1y = 12.5;  s1h = 1.85
    s2y = 9.55;  s2h = 2.5
    s3y = 5.35;  s3h = 3.2

    # ═══════════════ INPUT ═══════════════
    box(ax, CX, iy, W, ih, lw=1.2)
    ax.text(CX, iy + 0.28, 'Environmental Variables',
            ha='center', fontsize=15, fontweight='bold', color=BK)
    ax.text(CX, iy - 0.02, 'p raw variables  ·  n observations',
            ha='center', fontsize=10, color=GY)
    draw_icons(ax, CX, iy - 0.38, s=0.72)

    # ──── VIF arrow ────
    arr(ax, CX, iy - ih / 2 - GAP, CX, s1y + s1h / 2 + GAP)
    mid_v = (iy - ih / 2 + s1y + s1h / 2) / 2
    ax.text(CX + 0.2, mid_v, 'VIF filtering',
            ha='left', fontsize=8.5, color=GY)

    # ═══════════════ STAGE 1 ═══════════════
    box(ax, CX, s1y, W, s1h, lw=1.2)
    ax.text(LX + 0.35, s1y + 0.58, 'Stage 1', ha='left',
            fontsize=13, fontweight='bold', color=BK)
    ax.text(LX + 1.75, s1y + 0.58, 'Causal Structure Learning', ha='left',
            fontsize=13, color=BK)
    ax.text(LX + 0.45, s1y + 0.08, 'Bootstrap Hill-Climbing',
            ha='left', fontsize=10, color=GY)
    ax.text(LX + 0.45, s1y - 0.20, 'B = 200  ·  BIC-Gaussian',
            ha='left', fontsize=10, color=GY)

    draw_dag(ax, CX + W / 2 - 1.5, s1y - 0.02, s=0.78)

    ax.text(CX, s1y - 0.70, 'Consensus DAG  ·  out-degree  ·  edge strength',
            ha='center', fontsize=9, color=BK)

    arr(ax, CX, s1y - s1h / 2 - GAP, CX, s2y + s2h / 2 + GAP)

    # ═══════════════ STAGE 2 ═══════════════
    box(ax, CX, s2y, W, s2h, lw=1.2)
    ax.text(LX + 0.35, s2y + 0.98, 'Stage 2', ha='left',
            fontsize=13, fontweight='bold', color=BK)
    ax.text(LX + 1.75, s2y + 0.98, 'Causal Effect Estimation', ha='left',
            fontsize=13, color=BK)

    ate_cx  = CX - W / 4
    cate_cx = CX + W / 4

    ax.text(ate_cx, s2y + 0.48, 'ATE (Global)', ha='center',
            fontsize=11, fontweight='bold', color=BK)
    ax.text(ate_cx, s2y + 0.14, 'Double Machine Learning', ha='center',
            fontsize=8, color=GY)
    draw_forest(ax, ate_cx, s2y - 0.38, s=0.72)

    ax.plot([CX, CX], [s2y - 0.88, s2y + 0.68],
            color=LG, lw=0.5, ls=':', zorder=3)

    ax.text(cate_cx, s2y + 0.48, 'CATE (Spatial)', ha='center',
            fontsize=11, fontweight='bold', color=BK)
    ax.text(cate_cx, s2y + 0.14, 'Causal Forest  ·  2000 trees', ha='center',
            fontsize=8, color=GY)
    draw_heatmap(ax, cate_cx, s2y - 0.42, s=0.88)

    ax.text(CX, s2y - 1.05,
            'Global ATE per variable  +  Spatial CATE per location',
            ha='center', fontsize=9, color=BK)

    arr(ax, CX, s2y - s2h / 2 - GAP, CX, s3y + s3h / 2 + GAP)

    # ═══════════════ STAGE 3 ═══════════════
    box(ax, CX, s3y, W, s3h, lw=1.2)
    ax.text(LX + 0.35, s3y + 1.33, 'Stage 3', ha='left',
            fontsize=13, fontweight='bold', color=BK)
    ax.text(LX + 1.75, s3y + 1.33, 'Adaptive Screening  +  CI-MLP', ha='left',
            fontsize=13, color=BK)

    crit_y  = s3y + 0.58
    crit_xs = [LX + 0.95, LX + 2.1, LX + 3.25]
    crit_colors = [C_TEAL, C_RED, C_GREEN]
    for i, (lbl, col) in enumerate(zip(
            ['DAG\nOut-degree', 'ATE\nEffect size', 'RF\nImportance'],
            crit_colors)):
        box(ax, crit_xs[i], crit_y, 1.05, 0.58, lw=0.6)
        ax.text(crit_xs[i], crit_y, lbl, ha='center', va='center',
                fontsize=7, color=BK, linespacing=1.2)
        ax.plot(crit_xs[i], crit_y + 0.35, 's', color=col,
                ms=4, zorder=8, clip_on=False)

    merge_cx = crit_xs[1]
    for bx in crit_xs:
        arr(ax, bx, crit_y - 0.32, merge_cx, crit_y - 0.68, lw=0.5)

    score_y = s3y - 0.38
    box(ax, merge_cx, score_y, 3.2, 0.36, lw=0.6)
    ax.text(merge_cx, score_y, 'Adaptive weighted composite score',
            ha='center', va='center', fontsize=7.5, color=BK)

    ax.text(LX + 0.6, s3y - 0.92, 'ATE-weighted features',
            ha='left', fontsize=7.5, color=GY)
    ax.text(LX + 0.6, s3y - 1.18, '+ DAG interaction features',
            ha='left', fontsize=7.5, color=GY)

    nn_cx = CX + W / 2 - 1.35
    arr(ax, merge_cx + 1.65, score_y, nn_cx - 1.05, score_y, lw=0.8)

    draw_nn(ax, nn_cx, s3y - 0.30, s=0.95)
    ax.text(nn_cx, s3y - 1.08, 'CI-MLP', ha='center',
            fontsize=11, fontweight='bold', color=BK)
    ax.text(nn_cx, s3y - 1.36, 'Focal loss  ·  5 seeds', ha='center',
            fontsize=7, color=GY)

    # ──── 保存 ────
    fig.savefig('output/case2/CAST_technical_roadmap.png', dpi=1200,
                bbox_inches='tight', facecolor='white', edgecolor='none')
    print('[OK] PNG saved')
    fig.savefig('output/case2/CAST_technical_roadmap.svg', format='svg',
                bbox_inches='tight', facecolor='white', edgecolor='none')
    print('[OK] SVG saved')
    plt.close()


if __name__ == '__main__':
    main()
