#!/usr/bin/env python3

"""
Plot a per-target breakdown of CAE matching results from a summary CSV,
showing the proportion of matches, distorted matches, and no matches
as a stacked bar chart. Adapted for CAE analysis pipelines.

Usage:
    python plot_cae_matching_breakdown.py summary.csv [output.svg]

"""

import sys
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import re

# Set all fonts globally to Arial, size 10
import matplotlib
matplotlib.rcParams.update({
    'font.size': 10,
    'font.family': 'Arial'
})

# CLI args
csv_path = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("cae_summary_0_5.csv")
out_path = Path(sys.argv[2]) if len(sys.argv) > 2 else csv_path.with_name(csv_path.stem + "_breakdown.svg")
csv_name = Path(csv_path).name

match = re.search(r'(\d+)_(\d+)\.csv$', csv_name)
if match:
    thr_str = f"{match.group(1)}.{match.group(2)}"
else:
    thr_str = "?"

# 1. Load
df = pd.read_csv(csv_path)

# 2. Compute totals and percentages
df['total'] = df['matched'] + df['distorted'] + df['no_match']
df['pct_matched'] = df['matched'] / df['total'] * 100
df['pct_distorted'] = df['distorted'] / df['total'] * 100
df['pct_no_match'] = df['no_match'] / df['total'] * 100

# 3. Coverage stats
overall_perfect = df['matched'].sum() / df['total'].sum() * 100
overall_total = (df['matched'].sum() + df['distorted'].sum()) / df['total'].sum() * 100

# 4. Plot setup
targets = df['target']
n = len(df)
ind = np.arange(n)
bar_width = 0.8

# Increase left/right margin by 1 bar each side
xbuffer = 1.0
fig, ax = plt.subplots(figsize=(6.93, 3.3))
fig.subplots_adjust(bottom=0.29)

colors = {
    'matched': '#1f8836',      # strong green
    'distorted': '#ffd700',    # gold/yellow
    'no_match': '#e15759',     # reddish
}

# Stacked bars
bars1 = ax.bar(ind, df['pct_matched'], bar_width, label='Match', color=colors['matched'], edgecolor='black', linewidth=0.4)
bars2 = ax.bar(ind, df['pct_distorted'], bar_width, bottom=df['pct_matched'], label='Distorted match', color=colors['distorted'], edgecolor='black', linewidth=0.4)
bars3 = ax.bar(ind, df['pct_no_match'], bar_width, bottom=df['pct_matched'] + df['pct_distorted'], label='No match', color=colors['no_match'], edgecolor='black', linewidth=0.4)

# Labels & ticks
ax.set_ylabel('Connected Atom Environments (%)', fontsize=10, fontname='Arial')
ax.set_xlabel('Target', fontsize=10, fontname='Arial')
# ax.set_title(f'CAE Matching For Similarity Threshold {thr_str}', fontsize=10, fontname='Arial', fontweight='bold')
ax.set_xticks(ind)
ax.set_xticklabels(targets, rotation=45, ha='right', rotation_mode='anchor', fontsize=8, fontname='Arial')
ax.tick_params(axis='both', labelsize=8)
ax.set_ylim(0, 103)  # a bit above 100%

# Y ticks formatting
ax.yaxis.set_major_formatter(
    plt.FuncFormatter(lambda x, _: f"{x:.1f}".rstrip('0').rstrip('.') if '.' in f"{x:.1f}" else f"{int(x)}")
)
for label in ax.get_yticklabels():
    label.set_fontsize(10)
    label.set_fontname('Arial')

# More x space left/right
group_width = bar_width
ax.set_xlim(-xbuffer, n - 1 + group_width + xbuffer - bar_width)

# Legend in bottom left *inside axes*
leg = ax.legend(
    title="CAE classification",
    loc='lower left',
    bbox_to_anchor=(0.02, 0),
    frameon=True,
    fontsize=8,
    title_fontsize=8,
    prop={'family': 'Arial', 'size': 8},
    framealpha=1
)
leg.get_frame().set_edgecolor('black')
leg.get_frame().set_linewidth(0.4)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

# Coverage box in bottom right inside axes (above the axis)
textstr = (
    f"Overall coverage: {overall_perfect:.1f}%\n"
    f"Overall total coverage:   {overall_total:.1f}%"
)
ax.text(
    0.965, 0.04, textstr,
    transform=ax.transAxes,
    fontsize=10,
    va='bottom', ha='right',
    fontname='Arial',
    bbox=dict(
        boxstyle="round,pad=0.3",
        facecolor="white",
        edgecolor="black",
        linewidth=0.4,
        alpha=1)
)

fig.tight_layout(pad=0.7)
fig.savefig(out_path, format='svg')
print(f"✅  Plot saved to {out_path.resolve()}")

plt.show()
