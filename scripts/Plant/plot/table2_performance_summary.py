"""
Table 2 & Table S1: Model Performance Comparison — Plant
Booktabs-style academic tables (Nature/MEE journal format)
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from scipy.stats import wilcoxon

os.chdir("E:/CausalSDMs")
fig_dir = "figures/case4_plant/plot"
os.makedirs(fig_dir, exist_ok=True)

# ── Load data ────────────────────────────────────────────────────────────────
results = pd.read_csv("output/case4_plant/all_results_v3.csv")
dag_info = pd.read_csv("output/case4_plant/all_dag_info_v3.csv")
ate_data = pd.read_csv("output/case4_plant/all_ate_results_v3.csv")
sp_meta = pd.read_csv(
    "outputs/Plant/Res9/CAST_ready/CAST_Species_Summary.csv"
)
sp_meta["species"] = sp_meta["species"].str.replace(" ", "_")

n_species = results["species"].nunique()
model_order = ["CAST", "MLP_ATE", "MLP", "RF", "BRT", "Maxent"]

# ══════════════════════════════════════════════════════════════════════════════
# Compute all statistics
# ══════════════════════════════════════════════════════════════════════════════
summary_rows = []
for model in model_order:
    mdf = results[results["model"] == model].copy()
    if mdf.empty:
        continue
    auc = mdf["auc_mean"].dropna()
    tss = mdf["tss_mean"].dropna()
    row = {
        "model": model,
        "auc_mean": auc.mean(), "auc_sd": auc.std(),
        "auc_min": auc.min(), "auc_max": auc.max(),
        "tss_mean": tss.mean(), "tss_sd": tss.std(),
        "tss_min": tss.min(), "tss_max": tss.max(),
        "n_vars": mdf["n_vars"].mean(),
        "n_features": mdf["n_features_total"].mean(),
        "n_interactions": mdf["n_interactions"].mean() if model == "CAST" else 0,
        "dag_density": mdf["dag_density"].mean() if model == "CAST" else 0,
        "n": len(auc),
    }
    summary_rows.append(row)
summary_df = pd.DataFrame(summary_rows)

# Pairwise vs CAST
cast_sp = results[results["model"] == "CAST"][["species", "auc_mean", "tss_mean"]].copy()
cast_sp.columns = ["species", "cast_auc", "cast_tss"]

comp_rows = []
for model in model_order[1:]:
    mdf = results[results["model"] == model][["species", "auc_mean", "tss_mean"]].copy()
    mdf.columns = ["species", "other_auc", "other_tss"]
    mg = pd.merge(cast_sp, mdf, on="species", how="inner")
    da = mg["cast_auc"] - mg["other_auc"]
    dt = mg["cast_tss"] - mg["other_tss"]
    w, t, l = (da > 1e-4).sum(), ((da >= -1e-4) & (da <= 1e-4)).sum(), (da < -1e-4).sum()
    try:
        ws, wp = wilcoxon(mg["cast_auc"], mg["other_auc"])
    except Exception:
        ws, wp = np.nan, np.nan
    comp_rows.append({
        "model": model, "dauc_mean": da.mean(), "dauc_sd": da.std(),
        "dtss_mean": dt.mean(),
        "win": w, "tie": t, "loss": l,
        "win_pct": w / len(da) * 100,
        "wilcox_p": wp, "n": len(mg)
    })
comp_df = pd.DataFrame(comp_rows)

# ══════════════════════════════════════════════════════════════════════════════
# Booktabs-style rendering helper
# ══════════════════════════════════════════════════════════════════════════════
def draw_booktabs_table(ax, header_rows, data_rows, group_labels=None,
                        col_widths=None, bold_cols=None, bold_best=None,
                        title=None, footnote=None):
    """
    Draw a clean booktabs-style table.
    header_rows: list of lists (can be multi-level)
    data_rows: list of dicts with 'cells', 'group' (optional), 'bold_mask' (optional)
    """
    ax.axis("off")

    n_cols = len(header_rows[0])
    n_data = len(data_rows)
    n_header = len(header_rows)

    # Layout params
    row_h = 0.038
    header_h = 0.042
    x_start = 0.04
    x_end = 0.96
    total_w = x_end - x_start

    if col_widths is None:
        col_widths = [total_w / n_cols] * n_cols

    # Normalize col_widths to total_w
    cw_sum = sum(col_widths)
    col_widths = [w / cw_sum * total_w for w in col_widths]

    # Colors
    c_bg_header = "#F2F2F2"
    c_rule_thick = "#333333"
    c_rule_thin = "#CCCCCC"
    c_text = "#333333"
    c_text_light = "#666666"
    c_bg_group = "#F7F7F7"

    # Font
    font_props = {"family": "sans-serif", "size": 9}
    font_header = {"family": "sans-serif", "size": 9, "weight": "bold"}
    font_bold = {"family": "sans-serif", "size": 9, "weight": "bold"}
    font_group = {"family": "sans-serif", "size": 9, "weight": "bold"}

    # Calculate y positions (top to bottom)
    if title:
        y_title = 0.97
    y_top = 0.92
    y_after_header = y_top - n_header * header_h
    y_bottom = y_after_header - n_data * row_h
    if footnote:
        y_footnote = y_bottom - 0.03

    # ── Top rule ─────────────────────────────────────────────────────────
    ax.plot([x_start, x_end], [y_top, y_top],
            color=c_rule_thick, linewidth=1.5, clip_on=False,
            transform=ax.transAxes)

    # ── Header background ────────────────────────────────────────────────
    from matplotlib.patches import FancyBboxPatch
    header_rect = plt.Rectangle(
        (x_start, y_after_header), total_w, n_header * header_h,
        facecolor=c_bg_header, edgecolor="none",
        transform=ax.transAxes, clip_on=False
    )
    ax.add_patch(header_rect)

    # ── Header text ──────────────────────────────────────────────────────
    for hi, hrow in enumerate(header_rows):
        y_h = y_top - (hi + 0.5) * header_h
        x_pos = x_start
        for ci, txt in enumerate(hrow):
            x_center = x_pos + col_widths[ci] / 2
            ax.text(x_center, y_h, txt,
                    ha="center", va="center", color=c_text,
                    fontdict=font_header, transform=ax.transAxes)
            x_pos += col_widths[ci]

    # ── Mid rule (after header) ──────────────────────────────────────────
    ax.plot([x_start, x_end], [y_after_header, y_after_header],
            color=c_rule_thick, linewidth=1.0, clip_on=False,
            transform=ax.transAxes)

    # ── Data rows ────────────────────────────────────────────────────────
    prev_group = None
    for ri, drow in enumerate(data_rows):
        y_row = y_after_header - (ri + 0.5) * row_h
        y_row_top = y_after_header - ri * row_h
        y_row_bot = y_after_header - (ri + 1) * row_h

        cells = drow["cells"]
        group = drow.get("group", None)
        is_group_header = drow.get("is_group_header", False)
        bold_mask = drow.get("bold_mask", [False] * n_cols)

        # Group separator line
        if group and group != prev_group and not is_group_header:
            ax.plot([x_start, x_end], [y_row_top, y_row_top],
                    color=c_rule_thin, linewidth=0.5, clip_on=False,
                    transform=ax.transAxes)

        # Group header row
        if is_group_header:
            ax.plot([x_start, x_end], [y_row_top, y_row_top],
                    color=c_rule_thin, linewidth=0.7, clip_on=False,
                    transform=ax.transAxes)
            # Group header background
            grp_rect = plt.Rectangle(
                (x_start, y_row_bot), total_w, row_h,
                facecolor=c_bg_group, edgecolor="none",
                transform=ax.transAxes, clip_on=False
            )
            ax.add_patch(grp_rect)
            ax.text(x_start + 0.01, y_row, cells[0],
                    ha="left", va="center", color=c_text,
                    fontdict=font_group, transform=ax.transAxes)
            prev_group = group
            continue

        # Cell text
        x_pos = x_start
        for ci, txt in enumerate(cells):
            x_center = x_pos + col_widths[ci] / 2
            fd = font_bold if bold_mask[ci] else font_props
            tc = c_text if bold_mask[ci] else c_text_light
            # First column (model name) slightly left-aligned with indent
            if ci == 0:
                ax.text(x_pos + 0.015, y_row, txt,
                        ha="left", va="center", color=c_text,
                        fontdict=font_props, transform=ax.transAxes)
            else:
                ax.text(x_center, y_row, txt,
                        ha="center", va="center", color=tc,
                        fontdict=fd, transform=ax.transAxes)
            x_pos += col_widths[ci]

        prev_group = group

    # ── Bottom rule ──────────────────────────────────────────────────────
    ax.plot([x_start, x_end],
            [y_after_header - n_data * row_h, y_after_header - n_data * row_h],
            color=c_rule_thick, linewidth=1.5, clip_on=False,
            transform=ax.transAxes)

    # ── Title ────────────────────────────────────────────────────────────
    if title:
        ax.text(0.5, y_title, title,
                ha="center", va="center", color=c_text,
                fontdict={"family": "sans-serif", "size": 12, "weight": "bold"},
                transform=ax.transAxes)

    # ── Footnote ─────────────────────────────────────────────────────────
    if footnote:
        ax.text(x_start, y_footnote, footnote,
                ha="left", va="top", color="#888888",
                fontdict={"family": "sans-serif", "size": 7},
                transform=ax.transAxes, wrap=True)

    return y_bottom


# ══════════════════════════════════════════════════════════════════════════════
# TABLE 2: Main performance summary
# ══════════════════════════════════════════════════════════════════════════════
model_display = {
    "CAST": "CAST", "MLP_ATE": "MLP+ATE", "MLP": "MLP",
    "RF": "RF", "BRT": "BRT", "Maxent": "Maxent"
}
model_group_map = {
    "CAST": "Causal (proposed)",
    "MLP_ATE": "Ablation", "MLP": "Ablation",
    "RF": "Traditional SDM", "BRT": "Traditional SDM", "Maxent": "Traditional SDM"
}

# Find best values for bolding
best_auc = summary_df["auc_mean"].max()
best_tss = summary_df["tss_mean"].max()

header = [["Model", "Mechanism", "Test AUC", "Test TSS", "ΔAUC vs MLP", "W / T / L", "p-value"]]
col_w = [1.5, 1.6, 1.4, 1.4, 1.1, 1.0, 0.9]

mechanism_map = {
    "CAST": "DAG-guided Sparse",
    "MLP_ATE": "ATE-Weighted", 
    "MLP": "Flat Dense",
    "RF": "Tree Ensemble", 
    "BRT": "Tree Ensemble", 
    "Maxent": "Maxent Baseline"
}

# Build ΔAUC vs MLP for CAST
cast_vs_mlp = comp_df[comp_df["model"] == "MLP"]
cast_dauc_mlp = cast_vs_mlp["dauc_mean"].values[0] if len(cast_vs_mlp) > 0 else 0

data_rows = []
prev_grp = None
for _, row in summary_df.iterrows():
    m = row["model"]
    grp = model_group_map[m]

    # Insert group header
    if grp != prev_grp:
        data_rows.append({
            "cells": [grp] + [""] * 6,
            "group": grp,
            "is_group_header": True
        })
        prev_grp = grp

    # AUC and TSS
    auc_str = f"{row['auc_mean']:.4f}±{row['auc_sd']:.4f}"
    tss_str = f"{row['tss_mean']:.4f}±{row['tss_sd']:.4f}"
    mechanism = mechanism_map.get(m, "—")

    # ΔAUC vs MLP
    if m == "CAST":
        dauc_str = f"+{cast_dauc_mlp:.4f}" if cast_dauc_mlp > 0 else f"{cast_dauc_mlp:.4f}"
    elif m == "MLP":
        dauc_str = "— (ref)"
    else:
        cr = comp_df[comp_df["model"] == m]
        if len(cr) > 0:
            mlp_auc = summary_df[summary_df["model"] == "MLP"]["auc_mean"].values[0]
            d = row["auc_mean"] - mlp_auc
            dauc_str = f"{d:+.4f}"
        else:
            dauc_str = "—"

    # Win/Tie/Loss vs CAST
    cr = comp_df[comp_df["model"] == m]
    if len(cr) > 0:
        c = cr.iloc[0]
        wtl_str = f"{c['win']:.0f}/{c['tie']:.0f}/{c['loss']:.0f}"
        wp = c["wilcox_p"]
        if np.isnan(wp):
            p_str = "n/a"
        elif wp < 0.001:
            p_str = "< 0.001"
        else:
            p_str = f"{wp:.3f}"
    else:
        wtl_str = "—"
        p_str = "—"

    # Bold mask
    bold = [False] * 7
    if abs(row["auc_mean"] - best_auc) < 1e-6:
        bold[2] = True
    if abs(row["tss_mean"] - best_tss) < 1e-6:
        bold[3] = True

    data_rows.append({
        "cells": [f"  {model_display[m]}", mechanism, auc_str, tss_str,
                  dauc_str, wtl_str, p_str],
        "group": grp,
        "bold_mask": bold
    })

fig1, ax1 = plt.subplots(figsize=(14, 6))
draw_booktabs_table(
    ax1, header, data_rows, col_widths=col_w,
    title=f"Table 2  Model performance comparison across {n_species} plant species",
    footnote=(
        "Notes: All models share the same 11 foundational eco-hydroclimatic drivers. "
        "AUC and TSS are reported as mean±sd across species. "
        "Bold indicates the best-performing models. "
        "W/T/L = Win/Tie/Loss comparison against CAST (AUC). "
        "p-value from paired Wilcoxon signed-rank test."
    )
)
plt.savefig(os.path.join(fig_dir, "table2_performance_summary.png"),
            dpi=1200, bbox_inches="tight", facecolor="white")
plt.close()
print("✓ Saved table2_performance_summary.png")

# Save CSV
summary_df.to_csv(os.path.join(fig_dir, "table2_performance_summary.csv"),
                   index=False, float_format="%.6f")
comp_df.to_csv(os.path.join(fig_dir, "table2_pairwise_comparisons.csv"),
               index=False, float_format="%.6f")
print("✓ Saved CSV files")

# ══════════════════════════════════════════════════════════════════════════════
# TABLE S1: Per-species detail
# ══════════════════════════════════════════════════════════════════════════════
print("\n═══ Generating Table S1 (Per-Species Detail) ═══")

ate_summary = ate_data.groupby("species").agg(
    n_sig=("significant", "sum"),
).reset_index()

wide = results.pivot_table(
    index="species", columns="model", values="auc_mean", aggfunc="first"
).reset_index()
wide = wide.merge(sp_meta[["species", "family"]], on="species", how="left")
wide = wide.merge(dag_info[["species", "n_edges", "dag_density"]], on="species", how="left")
wide = wide.merge(ate_summary, on="species", how="left")

if "CAST" in wide.columns and "MLP" in wide.columns:
    wide["delta_mlp"] = wide["CAST"] - wide["MLP"]
if "CAST" in wide.columns and "RF" in wide.columns:
    wide["delta_rf"] = wide["CAST"] - wide["RF"]

model_cols_present = [c for c in model_order if c in wide.columns]
wide["best"] = wide[model_cols_present].idxmax(axis=1)
wide = wide.sort_values(["family", "species"])

# Header
s1_header = [["Species", "Family", "CAST", "MLP+ATE", "MLP", "RF", "BRT", "Maxent",
              "ΔAUC\nvs MLP", "ΔAUC\nvs RF", "DAG\ndensity", "#DAG\nedges", "#sig\nATEs", "Best"]]
s1_colw = [2.0, 1.2, 1.0, 1.0, 1.0, 0.8, 0.8, 0.8, 0.9, 0.9, 0.8, 0.7, 0.6, 0.8]

s1_rows = []
prev_fam = None
for _, r in wide.iterrows():
    fam = r.get("family", "")
    if fam != prev_fam:
        s1_rows.append({
            "cells": [fam if fam else "Unknown"] + [""] * 13,
            "group": fam, "is_group_header": True
        })
        prev_fam = fam

    sp_name = r["species"].replace("_", " ")

    # Find best AUC for this species
    sp_aucs = {m: r.get(m, np.nan) for m in model_cols_present}
    sp_best_val = max([v for v in sp_aucs.values() if pd.notna(v)], default=0)

    cells = [f"  {sp_name}", ""]
    bold = [False, False]

    for m in model_order:
        v = r.get(m, np.nan)
        if pd.notna(v):
            cells.append(f"{v:.4f}")
            bold.append(abs(v - sp_best_val) < 1e-6)
        else:
            cells.append("—")
            bold.append(False)

    # ΔAUC vs MLP, vs RF
    dm = r.get("delta_mlp", np.nan)
    dr = r.get("delta_rf", np.nan)
    cells.append(f"{dm:+.4f}" if pd.notna(dm) else "—")
    bold.append(False)
    cells.append(f"{dr:+.4f}" if pd.notna(dr) else "—")
    bold.append(False)

    # DAG density, edges, sig ATEs
    dd = r.get("dag_density", np.nan)
    cells.append(f"{dd:.3f}" if pd.notna(dd) else "—")
    bold.append(False)
    ne = r.get("n_edges", np.nan)
    cells.append(f"{int(ne)}" if pd.notna(ne) else "—")
    bold.append(False)
    ns = r.get("n_sig", np.nan)
    cells.append(f"{int(ns)}" if pd.notna(ns) else "—")
    bold.append(False)

    # Best model
    cells.append(r.get("best", "—"))
    bold.append(False)

    s1_rows.append({
        "cells": cells, "group": fam, "bold_mask": bold
    })

n_total_rows = len(s1_rows)
fig_h = max(10, n_total_rows * 0.35 + 4)
fig2, ax2 = plt.subplots(figsize=(18, fig_h))

draw_booktabs_table(
    ax2, s1_header, s1_rows, col_widths=s1_colw,
    title=f"Table S1  Per-species performance detail ({n_species} species, Plant)",
    footnote=(
        "Notes: Bold = best AUC for each species. Species grouped by taxonomic family. "
        "ΔAUC = CAST − competitor (positive = CAST advantage). "
        "#sig ATEs = number of environmental variables with significant causal effect (p < 0.05). "
        "DAG density = proportion of possible edges retained in consensus causal graph."
    )
)

plt.savefig(os.path.join(fig_dir, "tableS1_per_species_detail.png"),
            dpi=1200, bbox_inches="tight", facecolor="white")
plt.close()
print("✓ Saved tableS1_per_species_detail.png")

wide.to_csv(os.path.join(fig_dir, "tableS1_per_species_detail.csv"),
            index=False, float_format="%.4f")
print("✓ Saved tableS1_per_species_detail.csv")

# ══════════════════════════════════════════════════════════════════════════════
# Generating Markdown equivalents
# ══════════════════════════════════════════════════════════════════════════════
print("\n═══ Generating Markdown Tables ═══")

def to_md_table(header_row, d_rows):
    h = [str(x).replace('\n', '<br>') for x in header_row]
    res = ["| " + " | ".join(h) + " |", "|" + "|".join(["---"]*len(h)) + "|"]
    for r in d_rows:
        is_grp = r.get("is_group_header", False)
        cells = r["cells"]
        mask = r.get("bold_mask", [False]*len(cells))
        fmt = []
        for i, c in enumerate(cells):
            v = str(c).replace('\n', '<br>').strip()
            if is_grp and i == 0:
                fmt.append(f"**{v}**")
            elif i < len(mask) and mask[i] and v != "":
                fmt.append(f"**{v}**")
            else:
                fmt.append(v)
        res.append("| " + " | ".join(fmt) + " |")
    return "\n".join(res)

try:
    with open(os.path.join(fig_dir, "table2_performance_summary.md"), "w", encoding="utf-8") as f:
        f.write(f"### Table 2: Model performance comparison across {n_species} plant species (Plant, China)\n\n")
        f.write(to_md_table(header[0], data_rows))
        f.write("\n\n*Notes: n = 32 species. AUC = Area Under ROC Curve; TSS = True Skill Statistic. ")
        f.write("Bold = best across all models. Win/Tie/Loss = species-level comparison vs CAST (AUC). ")
        f.write("Wilcoxon p = paired signed-rank test (CAST vs. competitor). ")
        f.write("ΔAUC vs MLP = mean AUC difference relative to Flat MLP baseline. ")
        f.write("CAST uses causal DAG-guided interaction features; MLP+ATE uses ATE-weighted inputs only.*\n")
    print("✓ Saved table2_performance_summary.md")

    with open(os.path.join(fig_dir, "tableS1_per_species_detail.md"), "w", encoding="utf-8") as f:
        f.write(f"### Table S1: Per-species performance detail ({n_species} species, Plant)\n\n")
        f.write(to_md_table(s1_header[0], s1_rows))
        f.write("\n\n*Notes: Bold = best AUC for each species. Species grouped by taxonomic family. ")
        f.write("ΔAUC = CAST − competitor (positive = CAST advantage). ")
        f.write("#sig ATEs = number of environmental variables with significant causal effect (p < 0.05). ")
        f.write("DAG density = proportion of possible edges retained in consensus causal graph.*\n")
    print("✓ Saved tableS1_per_species_detail.md")
except Exception as e:
    print(f"Error generating markdown: {e}")

print("\n══════════ Table generation complete ══════════")


