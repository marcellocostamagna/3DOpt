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

# 2. Read task grouping and group order files
task_grouping_df = pd.read_csv("../task_grouping/task_grouping")
group_order_df = pd.read_csv("../task_grouping/group_order")

# Create task-to-group_id mapping (task -> integer group ID)
# Strip numeric prefix from task names to match with refcodes in task_grouping
def get_refcode(task_name):
    """Extract refcode by removing numeric prefix (e.g., '1_ABAHIW' -> 'ABAHIW')"""
    return re.sub(r'^\d+_', '', str(task_name))

task_to_group_id = dict(zip(task_grouping_df['task'], task_grouping_df['group']))

# Create group_id-to-info mappings
group_id_to_name = dict(zip(group_order_df['group'], group_order_df['name']))
group_id_to_position = dict(zip(group_order_df['group'], group_order_df['position']))

# Get original targets
original_targets = df['target'].tolist()

# Assign each target to its group ID and position
target_info = []
for target in original_targets:
    refcode = get_refcode(target)
    group_id = task_to_group_id.get(refcode, -1)  # -1 for unknown targets
    position = group_id_to_position.get(group_id, 999)  # Unknown groups go to end
    target_info.append((target, group_id, position))

# Sort targets by group position, maintaining original order within groups
target_info_sorted = sorted(target_info, key=lambda x: (x[2], original_targets.index(x[0])))
reordered_targets = [t[0] for t in target_info_sorted]

# Create reordering index
reorder_idx = [original_targets.index(target) for target in reordered_targets]

# Reorder dataframe
df = df.iloc[reorder_idx].reset_index(drop=True)

# Calculate group spans for visual indicators
task_groups = []
current_group_id = None
group_start = 0
for idx, target in enumerate(df['target']):
    refcode = get_refcode(target)
    group_id = task_to_group_id.get(refcode, -1)
    if group_id != current_group_id:
        if current_group_id is not None:
            # Get the group name for display
            group_name = group_id_to_name.get(current_group_id, 'Unknown')
            task_groups.append((group_start, idx - 1, group_name))
        current_group_id = group_id
        group_start = idx
# Add the last group
if current_group_id is not None:
    group_name = group_id_to_name.get(current_group_id, 'Unknown')
    task_groups.append((group_start, len(df) - 1, group_name))

# 3. Compute totals and percentages
df['total'] = df['matched'] + df['distorted'] + df['no_match']
df['pct_matched'] = df['matched'] / df['total'] * 100
df['pct_distorted'] = df['distorted'] / df['total'] * 100
df['pct_no_match'] = df['no_match'] / df['total'] * 100

# 4. Coverage stats
overall_perfect = df['matched'].sum() / df['total'].sum() * 100
overall_total = (df['matched'].sum() + df['distorted'].sum()) / df['total'].sum() * 100

# 5. Plot setup
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

# Add visual group indicators above the plot
y_pos = 1.08  # Position above the plot (in axis coordinates)
line_y = 1.04  # Y position for the bracket lines
extra_x_buffer = 0.4

for start_idx, end_idx, group_name in task_groups:
    # Calculate x positions for the group span (centered on bars)
    x_start = start_idx - extra_x_buffer
    x_end = end_idx + extra_x_buffer
    x_center = (x_start + x_end) / 2
    
    # Draw horizontal line for the group
    ax.plot([x_start, x_end], [line_y, line_y], 
            color='black', linewidth=0.8, 
            transform=ax.get_xaxis_transform(), clip_on=False)
    
    # Add vertical ticks at the ends
    tick_height = 0.015
    ax.plot([x_start, x_start], [line_y - tick_height, line_y + tick_height], 
            color='black', linewidth=0.8, 
            transform=ax.get_xaxis_transform(), clip_on=False)
    ax.plot([x_end, x_end], [line_y - tick_height, line_y + tick_height], 
            color='black', linewidth=0.8, 
            transform=ax.get_xaxis_transform(), clip_on=False)
    
    # Add group label
    ax.text(x_center, y_pos, group_name, 
            ha='center', va='bottom', fontsize=9, fontname='Arial',
            transform=ax.get_xaxis_transform())

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
