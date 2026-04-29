#!/usr/bin/env python3
"""
Extract and save individual panel icons from the CAST technical roadmap figure.
Each panel is saved as a separate PNG at 2400 DPI with transparent background.

Run from project root: python scripts/EcoISEA3H/plot/fig1_panels_export.py
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..', '..'))
OUT_DIR = os.path.join(PROJECT_ROOT, 'figures', 'case2_eco', 'panels')
os.makedirs(OUT_DIR, exist_ok=True)

BK = 'black'
GY = '#777777'
LG = "#000000"

C_BLUE   = '#4472C4'
C_RED    = '#C0504D'
C_GREEN  = '#548235'
C_BROWN  = '#A0522D'
C_ORANGE = '#ED7D31'
C_TEAL   = '#2B8C8C'

DIV_CMAP = LinearSegmentedColormap.from_list(
    'div', [C_BLUE, '#F0F0F0', C_RED], N=256)


# ---------------------------------------------------------------------------
# Helper drawing functions (identical to the original script)
# ---------------------------------------------------------------------------

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


def draw_forest(ax, cx, cy, s=1.0, center_color=LG, center_lw=1.4, center_ls='--'):
    """
    在指定的 matplotlib 轴对象上绘制一棵"森林"风格的水平条带图。

    参数:
        ax  : matplotlib.axes.Axes，绘图区域
        cx  : float，中心 x 坐标
        cy  : float，中心 y 坐标
        s   : float，整体缩放因子，控制所有尺寸（长度、粗细、点大小等）
    """
    # 五个偏移量，代表每个条带相对于中心的水平位置（正值在右，负值在左）
    effs = [0.45, -0.20, 0.50, 0.20, -0.38]

    # 绘制一条从 (cx, cy-0.65*s) 到 (cx, cy+0.65*s) 的竖直虚线
    # center_color 控制线的颜色，可传入 'red'、'#FF0000' 等
    # center_lw 控制线的粗细，数值越大线越粗
    # center_ls 控制线型：'--' 虚线、'-' 实线、':' 点线、'-.' 点划线
    ax.plot([cx, cx], [cy - 0.65 * s, cy + 0.65 * s],
            color=center_color, lw=center_lw, ls=center_ls, zorder=3)

    # 循环绘制五个水平条带及末端的圆点
    for i, e in enumerate(effs):
        # 计算当前条带的 y 坐标：从上到下依次为 cy+0.52*s, cy+0.26*s, cy, cy-0.26*s, cy-0.52*s
        yy = cy + (2 - i) * 0.26 * s

        # 条带的半宽：基于偏移量 e 的绝对值，再利用缩放因子 s 控制整体长度
        ci = abs(e) * 0.6   # 该值决定了条带的水平长度（条带范围为 [e-ci, e+ci]）
        # 条带颜色：正偏移 e>0 用红色，负偏移用蓝色（C_RED, C_BLUE 需提前定义）
        col = C_RED if e > 0 else C_BLUE

        # 绘制水平条带（线段）
        # 线段端点：从 cx + (e-ci)*s 到 cx + (e+ci)*s，y 坐标均为 yy
        # lw=0.7 控制线条粗细，可增大数值使条带更粗，或乘以 s 实现随缩放变化
        # solid_capstyle='round' 使线段端点圆润（避免方形截断），改为 'butt' 可得到平头
        ax.plot([cx + (e - ci) * s, cx + (e + ci) * s], [yy, yy],
                color=col, lw=3.0, zorder=4, solid_capstyle='round')

        # 绘制条带末端的圆点（位于条带中央，即偏移量 e 处）
        # ms=2.8*s 控制圆点的大小（markersize），可通过修改该系数调整
        # 'o' 表示圆形，可改为 's'（方形）、'^'（三角形）等改变形状
        # 颜色 col 与条带相同，若需不同颜色可单独设置
        # 注意：圆点没有描边（edgecolor 默认与 facecolor 相同），如需描边需添加 edgecolor='k', linewidth=0.5 等参数
        ax.plot(cx + e * s, yy, 'o', color=col, ms=7.5 * s, zorder=5)

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


def draw_nn(ax, cx, cy, s=1.0):
    # ==================== 网络结构配置 ====================
    # layers: 定义每一层的神经元数量
    # 当前为 5 层：输入层4个 → 隐藏层1(6个) → 隐藏层2(6个) → 隐藏层3(4个) → 输出层1个
    # 【可调】修改列表中的数字即可改变各层神经元数量
    layers = [4, 6, 6, 4, 1]

    # ==================== 计算各层的水平位置 ====================
    # lx: 每一层圆心的 x 坐标
    # range(5) 产生 0,1,2,3,4，减去 2 后以 cx 为中心对称分布
    # 0.44 * s 是相邻两层之间的水平间距
    # 【可调】修改 0.44 可改变层与层之间的水平距离
    lx = [cx + (i - 2) * 0.4 * s for i in range(5)]

    # ==================== 计算每个神经元的坐标 ====================
    # positions: 二维列表，positions[li][j] 表示第 li 层第 j 个神经元的 (x, y) 坐标
    positions = []
    for x, n in zip(lx, layers):
        # 对于每一层，计算该层 n 个神经元的垂直分布
        # (n - 1) / 2 使该层神经元以 cy 为中心上下对称排列
        # 0.16 * s 是同一层内相邻神经元的垂直间距
        # 【可调】修改 0.16 可改变同一层内神经元之间的垂直距离
        positions.append([(x, cy + (j - (n - 1) / 2) * 0.20 * s) for j in range(n)])

    # ==================== 绘制神经元之间的连线 ====================
    # 遍历每一对相邻层
    for li in range(len(layers) - 1):
        # 遍历当前层的每一个神经元
        for x1, y1 in positions[li]:
            # 遍历下一层的每一个神经元
            for x2, y2 in positions[li + 1]:
                # 绘制从 (x1,y1) 到 (x2,y2) 的连线
                ax.plot([x1, x2], [y1, y2],
                        color="#000000",   # 【连线颜色】当前为近黑色，可改为 'red'、'#FF0000' 等
                        lw=0.8,           # 【连线粗细】线宽（line width），数值越大线越粗
                        zorder=3)          # 图层顺序，zorder 越大越在上层显示

    # ==================== 绘制神经元节点（圆形） ====================
    for ps in positions:          # 遍历每一层
        for x, y in ps:           # 遍历该层的每一个神经元坐标
            ax.add_patch(Circle((x, y),
                                0.066 * s,      # 【神经元大小】圆的半径，数值越大圆越大
                                fc=C_BLUE,      # 【填充颜色】facecolor，神经元内部颜色（C_BLUE 需外部定义，如 'skyblue'）
                                ec='blue',      # 【描边颜色】edgecolor，圆圈的边框颜色
                                lw=0.2,        # 【描边粗细】边框线宽，数值越大描边越粗
                                zorder=5))      # 图层顺序，确保圆覆盖在连线上方


def draw_icons(ax, cx, cy, s=1.0):
    gap = 1.5 * s
    xs = [cx - 1.5 * gap, cx - 0.5 * gap, cx + 0.5 * gap, cx + 1.5 * gap]
    names  = ['Climate', 'Topography', 'Hydrology', 'Land cover']
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
        elif name == 'Hydrology':
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


# ---------------------------------------------------------------------------
# Save a panel to PNG with transparent background
# ---------------------------------------------------------------------------

def save_panel(draw_fn, filename, figsize=(3, 3), s=1.0, extra_ax_calls=None):
    """Generic helper to draw a panel and save as PNG."""
    fig, ax = plt.subplots(figsize=figsize)
    ax.set_aspect('equal')
    ax.axis('off')
    if extra_ax_calls:
        extra_ax_calls(ax)
    draw_fn(ax, 0, 0, s=s)
    fig.savefig(os.path.join(OUT_DIR, filename),
                dpi=2400, bbox_inches='tight',
                facecolor='none', edgecolor='none')
    plt.close(fig)
    print(f'[OK] {filename}')


# ---------------------------------------------------------------------------
# Individual panel savers
# ---------------------------------------------------------------------------

def save_dag():
    """Causal DAG panel (Stage 1)."""
    save_panel(draw_dag, 'panel_dag.png', figsize=(2.5, 2.2), s=0.85)


def save_forest():
    """Forest plot (Stage 2 - ATE)."""
    def extra(ax):
        ax.set_xlim(-0.65, 0.65)
        ax.set_ylim(-0.80, 0.75)
    save_panel(draw_forest, 'panel_forest.png', figsize=(2.2, 2.8), s=0.72, extra_ax_calls=extra)


def save_heatmap():
    """Heatmap / spatial proxy panel (Stage 2 - CATE)."""
    def extra(ax):
        ax.set_xlim(-0.55, 0.55)
        ax.set_ylim(-0.35, 0.35)
    save_panel(draw_heatmap, 'panel_heatmap.png', figsize=(2.5, 2.2), s=0.88, extra_ax_calls=extra)


def save_nn():
    """Neural network diagram (Stage 3 - CI-MLP)."""
    def extra(ax):
        ax.set_xlim(-1.0, 1.0)
        ax.set_ylim(-0.60, 0.60)
    save_panel(draw_nn, 'panel_nn.png', figsize=(3.5, 2.2), s=0.95, extra_ax_calls=extra)


def save_icons():
    """Environmental variable icons (Climate, Topography, Hydrology, Land cover)."""
    fig, ax = plt.subplots(figsize=(7.5, 2.0))
    ax.set_aspect('equal')
    ax.axis('off')
    draw_icons(ax, 0, 0, s=0.72)
    fig.savefig(os.path.join(OUT_DIR, 'panel_icons.png'),
                dpi=2400, bbox_inches='tight',
                facecolor='none', edgecolor='none')
    plt.close(fig)
    print('[OK] panel_icons.png')


def save_criteria():
    """Stage 3 criteria boxes (DAG Out-degree, ATE Effect size, RF Importance)."""
    fig, ax = plt.subplots(figsize=(3.5, 1.5))
    ax.set_aspect('equal')
    ax.axis('off')
    crit_xs = [-1.15, 0.0, 1.15]
    crit_y  = 0.0
    crit_colors = [C_TEAL, C_RED, C_GREEN]
    labels = ['DAG\nOut-degree', 'ATE\nEffect size', 'RF\nImportance']
    for i, (bx, col, lbl) in enumerate(zip(crit_xs, crit_colors, labels)):
        ax.add_patch(FancyBboxPatch(
            (bx - 0.52, crit_y - 0.29), 1.04, 0.58,
            boxstyle="round,pad=0.06",
            facecolor='white', edgecolor=BK, linewidth=0.6, zorder=2))
        ax.text(bx, crit_y, lbl, ha='center', va='center',
                fontsize=7, color=BK, linespacing=1.2)
        ax.plot(bx, crit_y + 0.35, 's', color=col, ms=4, zorder=8, clip_on=False)
    fig.savefig(os.path.join(OUT_DIR, 'panel_criteria.png'),
                dpi=2400, bbox_inches='tight',
                facecolor='none', edgecolor='none')
    plt.close(fig)
    print('[OK] panel_criteria.png')


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    save_dag()
    save_forest()
    save_heatmap()
    save_nn()
    save_icons()
    save_criteria()
    print(f'\nAll panels saved to: {OUT_DIR}')


if __name__ == '__main__':
    main()
