#!/usr/bin/env python3
"""
Fig 1: CAST Methodological Framework Architecture.
Generates an elegant, left-to-right architecture diagram for the CAST framework.
Run from project root: python scripts/EcoISEA3H/plot/fig1_cast_framework_architecture.py
"""

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Circle, Rectangle, Polygon
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.patches as patches
import numpy as np
import os

# Paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.normpath(os.path.join(SCRIPT_DIR, '..', '..', '..'))
# Note: Saving to a general methodological path since this diagram is dataset-agnostic
FIG_DIR = os.path.join(PROJECT_ROOT, 'figures', 'methodology', 'plot')
os.makedirs(FIG_DIR, exist_ok=True)

# Aesthetic configuration
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'font.size': 10,
    'figure.facecolor': 'white',
})

# Colors (Maintaining the original scientific palette)
BK = '#222222'; GY = '#666666'; LG = '#EAEAEA'
C_BLUE = '#4472C4'; C_RED = '#C0504D'; C_GREEN = '#548235'
C_BROWN = '#A0522D'; C_ORANGE = '#ED7D31'; C_TEAL = '#2B8C8C'; C_PURPLE = '#7030A0'

# Panel background colors (Subtle pastels)
BG_DATA = '#F2F6FA'; BG_CAST = '#F4FAF4'; BG_CIMLP = '#FCF4F4'

def rbox(ax, x, y, w, h, fc='white', ec=BK, lw=1.0, pad=0.1, zorder=2, alpha=1.0, ls='-'):
    """Draw a rounded bounding box."""
    box = FancyBboxPatch((x, y), w, h, boxstyle=f"round,pad={pad}", 
                         facecolor=fc, edgecolor=ec, linewidth=lw, 
                         zorder=zorder, alpha=alpha, linestyle=ls)
    ax.add_patch(box)
    return box

def d_arr(ax, x1, y1, x2, y2, lw=1.5, color=GY, head_width=0.12):
    """Draw a straight directed arrow."""
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle=f"->,head_length={head_width*2.5},head_width={head_width}", 
                                color=color, lw=lw, shrinkA=0, shrinkB=0), zorder=3)

def iso_rect(ax, cx, cy, w, h, offset_y, fc, alpha=0.8):
    """Draw an isometric (slanted) rectangle to represent a map layer."""
    pts = np.array([
        [cx, cy + h/2 + offset_y],
        [cx + w/2, cy + offset_y],
        [cx, cy - h/2 + offset_y],
        [cx - w/2, cy + offset_y]
    ])
    ax.add_patch(Polygon(pts, closed=True, facecolor=fc, edgecolor='w', lw=0.5, alpha=alpha, zorder=4+int(offset_y*10)))

def draw_data_stack(ax, cx, cy, w, h):
    """Draw a stack of raster/environmental maps."""
    colors = [C_BLUE, C_TEAL, C_ORANGE, C_GREEN]
    for i, c in enumerate(colors):
        iso_rect(ax, cx, cy, w, h, i*0.4, c, alpha=0.5)
    
    # Simulate presence/absence points on top
    top_y = cy + 3*0.4
    np.random.seed(42)
    for _ in range(12):
        px = cx + np.random.uniform(-w/3, w/3)
        py = top_y + np.random.uniform(-h/4, h/4)
        ax.plot(px, py, 'o', color=C_RED, ms=3.5, zorder=12, markeredgecolor='white', markeredgewidth=0.5)

def draw_dag(ax, cx, cy, s=1.0):
    """Draw a Bayesian Network / DAG representation."""
    nodes = {
        1: (cx - 0.7*s, cy + 0.4*s),
        2: (cx - 0.7*s, cy - 0.4*s),
        3: (cx,         cy + 0.1*s),
        4: (cx,         cy - 0.5*s),
        5: (cx + 0.7*s, cy)
    }
    nc = {1: C_BLUE, 2: C_TEAL, 3: C_ORANGE, 4: C_GREEN, 5: C_RED}
    for a, b in [(1, 3), (2, 3), (2, 4), (3, 5), (4, 5)]:
        xa, ya = nodes[a]; xb, yb = nodes[b]
        ax.annotate('', xy=(xb, yb), xytext=(xa, ya),
                    arrowprops=dict(arrowstyle='->', color=GY, lw=1.2,
                                    shrinkA=7*s, shrinkB=7*s, connectionstyle='arc3,rad=0.08'), zorder=5)
    for i, (x, y) in nodes.items():
        ax.add_patch(Circle((x, y), 0.14*s, fc=nc[i], ec='white', lw=0.5, zorder=6))

def draw_tokens(ax, cx, cy, labels, colors, w=0.6, h=0.25, dy=0.35):
    """Draw vertically stacked token representations."""
    for i, (lbl, col) in enumerate(zip(labels, colors)):
        y = cy - i*dy
        ax.add_patch(Rectangle((cx - w/2, y - h/2), w, h, fc=col, ec='white', lw=0.8, alpha=0.85, zorder=4))
        if lbl:
            ax.text(cx, y, lbl, ha='center', va='center', color='white', fontsize=8, fontweight='bold', zorder=5)

def draw_nn_layers(ax, cx, cy, layers, w=0.45, s=0.35):
    """Draw a standard multi-layer perceptron."""
    lx = [cx + i*w for i in range(len(layers))]
    positions = []
    for x, n in zip(lx, layers):
        positions.append([(x, cy + (j - (n - 1) / 2) * s) for j in range(n)])
    for li in range(len(layers) - 1):
        for x1, y1 in positions[li]:
            for x2, y2 in positions[li + 1]:
                ax.plot([x1, x2], [y1, y2], color='#D0D8E8', lw=0.6, zorder=3)
    for i, ps in enumerate(positions):
        for x, y in ps:
            fc = C_BLUE if i < len(layers)-1 else C_RED
            ax.add_patch(Circle((x, y), 0.05, fc=fc, ec='white', lw=0.3, zorder=5))

def main():
    fig, ax = plt.subplots(1, 1, figsize=(16, 7.8), dpi=300)
    ax.set_xlim(0, 16)
    ax.set_ylim(0, 7.8)
    ax.axis('off')

    # --------------------------------------------------------------------------
    # (a) Input Data Panel
    # --------------------------------------------------------------------------
    rbox(ax, 0.4, 0.5, 2.2, 6.4, fc=BG_DATA, ec=C_BLUE, lw=1.2, pad=0.0)
    ax.text(1.5, 6.5, "(a) Observational Data", ha='center', fontsize=13, fontweight='bold', color=C_BLUE)
    
    # Env Grids
    ax.text(1.5, 5.7, "Environmental Grids", ha='center', fontsize=11, fontweight='bold')
    draw_data_stack(ax, 1.5, 4.0, 1.6, 1.1)
    ax.text(1.5, 3.4, r"$X_{env} \in \mathbb{R}^{N \times P_{all}}$", ha='center', fontsize=10)
    
    # Presences
    ax.text(1.5, 2.6, "Species Occurrences", ha='center', fontsize=11, fontweight='bold')
    rbox(ax, 0.8, 1.3, 1.4, 1.0, fc='white', pad=0)
    ax.plot([0.9, 2.1], [2.0, 2.0], color=GY, lw=0.6)
    ax.text(1.5, 2.1, "Presence (1)", ha='center', va='bottom', fontsize=9, color=C_RED)
    ax.text(1.5, 1.5, "Pseudo-Absence (0)", ha='center', va='bottom', fontsize=9, color=C_BLUE)
    ax.text(1.5, 1.0, r"$Y_{obs} \in \{0, 1\}^N$", ha='center', fontsize=10)

    # Transition Arrow
    d_arr(ax, 2.6, 3.8, 3.2, 3.8, lw=3, color=C_BLUE, head_width=0.2) 

    # --------------------------------------------------------------------------
    # (b) CAST Core Engine Panel
    # --------------------------------------------------------------------------
    rbox(ax, 3.2, 0.5, 5.5, 6.4, fc=BG_CAST, ec=C_GREEN, lw=1.2, pad=0.0)
    ax.text(5.95, 6.5, "(b) Causal Analysis & Screening (CAST)", ha='center', fontsize=13, fontweight='bold', color=C_GREEN)
    
    # Sub-module 1: Structure
    rbox(ax, 3.5, 4.6, 2.3, 1.8, fc='white', ec=C_GREEN, pad=0)
    ax.text(4.65, 6.1, "1. Structure Learning", ha='center', fontsize=11, fontweight='bold', color=BK)
    draw_dag(ax, 4.65, 5.3, s=1.1)
    ax.text(4.65, 4.75, "DAG Density & Edge Strength", ha='center', fontsize=8, color=GY)

    # Sub-module 2: Effect
    rbox(ax, 3.5, 2.0, 2.3, 1.8, fc='white', ec=C_GREEN, pad=0)
    ax.text(4.65, 3.5, "2. Effect Estimation", ha='center', fontsize=11, fontweight='bold', color=BK)
    ax.plot([4.65, 4.65], [2.3, 3.0], color=LG, lw=0.8, ls='--', zorder=3)
    effs = [0.2, -0.3, 0.05, 0.4]
    for i, e in enumerate(effs):
        y = 2.9 - i*0.18
        col = C_RED if e>0 else C_BLUE
        ci = abs(e)*0.3 + 0.05
        ax.plot([4.65+(e-ci), 4.65+(e+ci)], [y, y], color=col, lw=1, zorder=4, solid_capstyle='round')
        ax.plot(4.65+e, y, 'o', color=col, ms=4, zorder=5)
    ax.text(4.65, 2.15, "Global Average Treatment Effect", ha='center', fontsize=8, color=GY)

    # Sub-module 3: Screening
    rbox(ax, 6.4, 2.4, 2.0, 3.6, fc='white', ec=C_GREEN, pad=0)
    ax.text(7.4, 5.7, "3. Adaptive\nScreening", ha='center', fontsize=11, fontweight='bold', color=BK)
    
    rbox(ax, 6.6, 4.7, 1.6, 0.4, fc=LG, ec=GY, pad=0)
    ax.text(7.4, 4.9, "DAG Topology Rank", ha='center', va='center', fontsize=8.5)
    rbox(ax, 6.6, 4.1, 1.6, 0.4, fc=LG, ec=GY, pad=0)
    ax.text(7.4, 4.3, "ATE Significance", ha='center', va='center', fontsize=8.5)
    rbox(ax, 6.6, 3.5, 1.6, 0.4, fc=LG, ec=GY, pad=0)
    ax.text(7.4, 3.7, "RF Importance", ha='center', va='center', fontsize=8.5)
    
    # Internal routes in CAST
    d_arr(ax, 5.8, 5.5, 6.4, 4.9)
    d_arr(ax, 5.8, 2.9, 6.4, 4.3)
    
    ax.plot([7.4, 7.4], [3.5, 2.9], color=GY, lw=1.2)
    ax.text(7.4, 2.7, "Causal Subset", ha='center', va='center', fontsize=10, fontweight='bold', color=C_GREEN)
    ax.text(7.4, 2.48, r"$X_C \in \mathbb{R}^{N \times P_{causal}}$", ha='center', va='center', fontsize=9, color=GY)
    
    # Transition Arrow
    d_arr(ax, 8.7, 3.8, 9.3, 3.8, lw=3, color=C_GREEN, head_width=0.2) 

    # --------------------------------------------------------------------------
    # (c) CI-MLP Panel
    # --------------------------------------------------------------------------
    rbox(ax, 9.3, 0.5, 4.6, 6.4, fc=BG_CIMLP, ec=C_RED, lw=1.2, pad=0.0)
    ax.text(11.6, 6.5, "(c) Causal-Informed NN (CI-MLP)", ha='center', fontsize=13, fontweight='bold', color=C_RED)
    
    # Tokens column
    ax.text(10.1, 5.9, "Subset Features", ha='center', fontsize=9, fontweight='bold')
    colors = [C_BLUE, C_TEAL, C_ORANGE, C_BROWN]
    draw_tokens(ax, 9.8, 5.3, [r"$x_1$", r"$x_2$", r"$x_3$", r"$x_{k}$"], colors, w=0.5, h=0.25, dy=0.38)
    # The 'Missing/Screened out' ones depicted optionally (MaskSDM style hint)
    ax.text(9.8, 3.8, "[...]", ha='center', fontsize=14, color=GY)
    
    # Structural Ops / Injection from CAST
    ax.text(10.7, 5.9, "Structural Injection", ha='center', fontsize=9, fontweight='bold')
    
    style = "Simple, tail_width=0.5, head_width=4, head_length=8"
    kw = dict(arrowstyle=style, color=C_GREEN, alpha=0.3)
    # Inject ATE to weights
    ax.add_patch(patches.FancyArrowPatch((5.8, 2.1), (10.7, 4.8), connectionstyle="arc3,rad=0.15", **kw, zorder=1))
    rbox(ax, 10.4, 4.6, 0.6, 0.4, fc='white', ec=C_RED, pad=0)
    ax.text(10.7, 4.8, r"$\otimes \mathbf{W}_{ATE}$", ha='center', va='center', fontsize=9, fontweight='bold')
    d_arr(ax, 10.1, 4.8, 10.4, 4.8)
    
    # Inject DAG to interactions
    ax.add_patch(patches.FancyArrowPatch((5.8, 5.7), (10.7, 3.6), connectionstyle="arc3,rad=-0.15", **kw, zorder=1))
    rbox(ax, 10.4, 3.4, 0.6, 0.4, fc='white', ec=C_RED, pad=0)
    ax.text(10.7, 3.6, r"$\oplus \mathbf{X}_{int}$", ha='center', va='center', fontsize=9, fontweight='bold')
    
    # Embeddings / Concat (CISO style inspiration)
    ax.text(11.6, 5.9, "Causal Tensor", ha='center', fontsize=9, fontweight='bold')
    emb_colors = [C_BLUE, C_TEAL, C_ORANGE, C_BROWN, C_PURPLE, C_PURPLE]
    draw_tokens(ax, 11.6, 5.3, ["","","","","",""], emb_colors, w=0.5, h=0.18, dy=0.25)
    
    d_arr(ax, 11.0, 4.8, 11.3, 5.0)
    d_arr(ax, 11.0, 3.6, 11.3, 4.0)
    
    # Neural Net
    ax.text(12.8, 5.9, "Deep MLP Architecture", ha='center', fontsize=9, fontweight='bold')
    draw_nn_layers(ax, 12.2, 4.7, [6, 8, 4, 1], w=0.45, s=0.35)
    ax.text(12.9, 3.4, "Focal Loss", ha='center', fontsize=10, fontweight='bold', color=C_RED)
    
    # Output score
    d_arr(ax, 13.55, 4.7, 14.1, 4.7, lw=2, color=C_BROWN)

    # --------------------------------------------------------------------------
    # (d) Output Panel
    # --------------------------------------------------------------------------
    rbox(ax, 14.3, 0.5, 1.6, 6.4, fc='white', ec=C_BROWN, lw=1.2, pad=0.0, ls='--')
    ax.text(15.1, 6.5, "(d) Output", ha='center', fontsize=13, fontweight='bold', color=C_BROWN)
    
    ax.text(15.1, 5.4, "Habitat\nSuitability", ha='center', fontsize=10, fontweight='bold')
    iso_rect(ax, 15.1, 4.7, 1.3, 0.9, 0, C_RED, alpha=0.7)
    
    ax.text(15.1, 2.9, "Spatial\nCATE Maps", ha='center', fontsize=10, fontweight='bold')
    iso_rect(ax, 15.1, 2.2, 1.3, 0.9, 0, C_PURPLE, alpha=0.7)
    
    # Connect MLP to HS
    d_arr(ax, 15.1, 4.2, 15.1, 3.6, lw=1.5, color=GY)

    # Final rendering
    fig.savefig(os.path.join(FIG_DIR, 'fig1_cast_framework_architecture.png'), dpi=1200, bbox_inches='tight')
    fig.savefig(os.path.join(FIG_DIR, 'fig1_cast_framework_architecture.svg'), format='svg', bbox_inches='tight')
    print('[OK] fig1_cast_framework_architecture.png saved to methodology/plot.')
    plt.close()

if __name__ == '__main__':
    main()
