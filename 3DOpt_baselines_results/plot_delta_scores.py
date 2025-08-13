#!/usr/bin/env python3
"""
plot_deltas.py
Show per-task differences in 3DOpt scores: ChemGE minus Random Sampler,
for each 3D generator (CCDC, OBabel, RDKit).
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import re
import textwrap
import matplotlib as mpl

# Font/style
mpl.rcParams.update({
    'font.size': 10,
    'font.family': 'Arial'
})

# ---- File names (adjust if yours differ) ----
R_FILES = {
    'CCDC':   'Rnd_ccdc.csv',
    'OBabel': 'Rnd_obabel.csv',
    'RDKit':  'Rnd_rdkit.csv',
    # 'CSD entries': 'Rnd_ccdc_entries.csv',  # no ChemGE pair, so we omit
}
C_FILES = {
    'CCDC':   'ChemGE_ccdc.csv',
    'OBabel': 'ChemGE_obabel.csv',
    'RDKit':  'ChemGE_rdkit.csv',
}

colors = {
    'CCDC':   '#1f77b4',
    'OBabel': '#ff7f0e',
    'RDKit':  '#2ca02c',
}

def load_scores(path):
    """Return (df_per_task, agg_value_or_nan). Drops Aggregate from per-task."""
    df = pd.read_csv(path)
    if 'task' not in df.columns or '3DOpt_Score' not in df.columns:
        raise ValueError(f"{path} missing 'task' or '3DOpt_Score'.")

    # Aggregate
    mask_agg = df['task'].astype(str).str.strip().str.lower().eq('aggregate')
    agg = float(df.loc[mask_agg, '3DOpt_Score'].iloc[0]) if mask_agg.any() else np.nan

    # Per-task only
    df = df[~mask_agg].copy()
    return df[['task', '3DOpt_Score']].reset_index(drop=True), agg

def clean_task_label(t):
    return re.sub(r'^\d+_', '', str(t))

# Load and align
pairs = {}
bench_deltas = {}  # per generator (ChemGE - Random)
tasks = None

for gen, rfile in R_FILES.items():
    cfile = C_FILES.get(gen)
    if not cfile:
        continue
    df_r, agg_r = load_scores(rfile)
    df_c, agg_c = load_scores(cfile)

    # Align on task (inner join)
    merged = pd.merge(df_r, df_c, on='task', suffixes=('_R', '_C'))
    # Sort by numeric index prefix in task (e.g., "10_ABEHAU")
    def task_idx(t):
        try: return int(str(t).split('_', 1)[0])
        except: return 10**9
    merged = merged.sort_values(key=lambda s: s.map(task_idx), by='task')
    merged['delta'] = merged['3DOpt_Score_C'] - merged['3DOpt_Score_R']

    pairs[gen] = merged
    bench_deltas[gen] = (agg_c - agg_r) if (np.isfinite(agg_c) and np.isfinite(agg_r)) else np.nan

    if tasks is None:
        tasks = merged['task'].values

# Plot
if not pairs:
    raise RuntimeError("No matching ChemGE/Random pairs found to plot.")

n_tasks = len(tasks)
gens = list(pairs.keys())
n_g = len(gens)

bar_width = 0.7 / n_g
x = np.arange(n_tasks)

fig, ax = plt.subplots(figsize=(6.93, 2.75))
fig.subplots_adjust(bottom=0.30)

for i, gen in enumerate(gens):
    y = pairs[gen]['delta'].values
    ax.bar(
        x + i * bar_width, y, width=bar_width,
        label=gen, color=colors.get(gen, 'gray'),
        edgecolor='black', linewidth=0.4
    )

# Axes
ax.set_ylabel("Δ3DOpt Score (ChemGE − RanSam)", fontsize=8, fontname='Arial')
ax.set_xlabel("Task", fontsize=10, fontname='Arial')

ax.set_xticks(x + bar_width * (n_g - 1) / 2)
ax.set_xticklabels(
    [textwrap.fill(clean_task_label(t), 12) for t in tasks],
    rotation=45, ha="right", rotation_mode="anchor",
    fontsize=8, fontname='Arial'
)

# Zero line & grid
ax.axhline(0.0, color='black', linewidth=0.6)
ax.grid(axis='y', linestyle='--', linewidth=0.4, alpha=0.7)
ax.set_axisbelow(True)

# Y limits (auto with small padding)
ymin, ymax = ax.get_ylim()
pad = (ymax - ymin) * 0.05 if ymax > ymin else 0.05
ax.set_ylim(ymin - pad, ymax + pad)

# Legend
leg = ax.legend(
    loc='lower left',
    fontsize=7, prop={'family': 'Arial', 'size': 7},
    frameon=True, framealpha=1.0, edgecolor='black'
)
leg.get_frame().set_linewidth(0.4)

# Spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# X limits buffer
group_width = bar_width * n_g
buffer = bar_width * 1.2
left = 0 - buffer
right = n_tasks - 1 + group_width + buffer - bar_width
ax.set_xlim(left, right)

# Save
outpath = Path("Delta_3DOpt_scores.svg")
fig.tight_layout(pad=0.5)
fig.savefig(outpath, format='svg')
plt.show()
print(f"✅ Saved {outpath.resolve()}")
