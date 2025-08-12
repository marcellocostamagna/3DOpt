#!/usr/bin/env python3
"""
plotting_ChemGE.py
Plot grouped ChemGE 3DOpt scores with consistent styles and layout.
"""

import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import re
import textwrap

# ---- Set global font to Arial size 10 ----
import matplotlib
matplotlib.rcParams.update({
    'font.size': 10,
    'font.family': 'Arial'
})

# ------------------------------------------------------------
# 1) Define files and labels (or get from CLI)
# ------------------------------------------------------------
file_labels = [
    ('ChemGE_ccdc', 'CCDC'),
    ('ChemGE_obabel', 'OBabel'),
    ('ChemGE_rdkit', 'RDKit'),
]
if len(sys.argv) > 1:
    file_labels = [(Path(f).stem, Path(f).stem.replace("_", " ")) for f in sys.argv[1:]]

colors = ['#1f77b4', '#ff7f0e', '#2ca02c']  # CCDC, OBabel, RDKit

# ------------------------------------------------------------
# 2) Load and align data
# ------------------------------------------------------------
dfs = []
labels = []
bench_totals = {}  # label -> aggregate 3DOpt benchmark score

for fname, label in file_labels:
    df = pd.read_csv(fname + ".csv")
    if 'task' not in df.columns or '3DOpt_Score' not in df.columns:
        raise ValueError(f"{fname}.csv must contain 'task' and '3DOpt_Score' columns.")

    # Capture aggregate (benchmark score) if present
    agg_mask = df['task'].astype(str).str.strip().str.lower().eq('aggregate')
    if agg_mask.any():
        bench_totals[label] = float(df.loc[agg_mask, '3DOpt_Score'].iloc[0])
    else:
        bench_totals[label] = np.nan

    # Drop aggregate row for plotting
    df = df[~agg_mask].copy()

    dfs.append(df[['task', '3DOpt_Score']].reset_index(drop=True))
    labels.append(label)

# Align tasks based on the first dataframe
tasks = dfs[0]['task'].values
n_tasks = len(tasks)
n_methods = len(dfs)

# Build score matrix
all_scores = np.zeros((n_methods, n_tasks))
for i, df in enumerate(dfs):
    all_scores[i] = df['3DOpt_Score'].values

# ------------------------------------------------------------
# 3) Plotting setup
# ------------------------------------------------------------
bar_width = 0.7 / n_methods
x = np.arange(n_tasks)

fig, ax = plt.subplots(figsize=(6.93, 2.75))
fig.subplots_adjust(bottom=0.30)

for i in range(n_methods):
    ax.bar(
        [p + i * bar_width for p in x],
        all_scores[i],
        width=bar_width,
        label=labels[i],
        color=colors[i % len(colors)],
        edgecolor='black',
        linewidth=0.4
    )

# ------------------------------------------------------------
# 4) Axis styling and labels
# ------------------------------------------------------------
def clean_label(t):
    return re.sub(r'^\d+_', '', t)

ax.set_ylabel("3DOpt Score", fontsize=10, fontname='Arial')
ax.set_xlabel("Task", fontsize=10, fontname='Arial')

# Center ticks
ax.set_xticks([p + bar_width * (n_methods - 1) / 2 for p in x])
ax.set_xticklabels(
    [textwrap.fill(clean_label(t), 12) for t in tasks],
    rotation=45, ha="right", rotation_mode="anchor",
    fontsize=8, fontname='Arial'
)
ax.tick_params(axis='both', labelsize=8)

# Y axis 0–1 with clean ticks
ax.set_ylim(0, 1)
ax.yaxis.set_major_formatter(
    plt.FuncFormatter(lambda y, _: f"{y:.1f}".rstrip('0').rstrip('.') if '.' in f"{y:.1f}" else f"{int(y)}")
)
for label in ax.get_yticklabels():
    label.set_fontsize(10)
    label.set_fontname('Arial')

# Legend (kept lower right per your script)
legend = ax.legend(
    loc='lower left',
    fontsize=7,
    title_fontsize=7,
    prop={'family': 'Arial', 'size': 7},
    framealpha=1.0, frameon=True, edgecolor='black'
)
legend.get_frame().set_linewidth(0.4)

# Horizontal grid lines for readability
ax.grid(axis='y', linestyle='--', linewidth=0.4, alpha=0.7)
ax.set_axisbelow(True)

# Hide top and right spines
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Tight X limits with small buffer
group_width = bar_width * n_methods
buffer = bar_width * 1.2
left = 0 - buffer
right = n_tasks - 1 + group_width + buffer - bar_width
ax.set_xlim(left, right)

# ------------------------------------------------------------
# 5) Benchmark-score box 
# ------------------------------------------------------------
lines = []
for lab in labels:
    val = bench_totals.get(lab, np.nan)
    if not np.isnan(val):
        lines.append(f"{lab}: {val:.1f}")

box_text = (r"$\mathbf{3DOpt\ Benchmark\ Scores}$" + "\n" + "\n".join(lines)) if lines \
           else r"$\mathbf{3DOpt\ Benchmark\ Scores}$\n—"

ax.text(
    0.98, 0.05,  # slightly above the x-axis
    box_text,
    transform=ax.transAxes,
    ha='right', va='bottom',
    fontsize=8, fontname='Arial',
    bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black', linewidth=0.4)
)

# ------------------------------------------------------------
# 6) Save the figure
# ------------------------------------------------------------
fig.tight_layout(pad=0.5)
output_path = Path("3DOpt_ChemGE_Scores.svg")
fig.savefig(output_path, format='svg')
plt.show()
print(f"✅ Plot saved to {output_path.resolve()}")
