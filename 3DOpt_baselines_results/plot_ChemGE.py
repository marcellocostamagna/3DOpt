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

# ------------------------------------------------------------
# 3) Read task grouping and group order files
# ------------------------------------------------------------
task_grouping_df = pd.read_csv("../task_grouping/task_grouping")
group_order_df = pd.read_csv("../task_grouping/group_order")

# Create task-to-group_id mapping (task -> integer group ID)
# Strip numeric prefix from task_grouping to match with CSV task names
def get_refcode(task_name):
    """Extract refcode by removing numeric prefix (e.g., '1_ABAHIW' -> 'ABAHIW')"""
    return re.sub(r'^\d+_', '', str(task_name))

task_to_group_id = dict(zip(task_grouping_df['task'], task_grouping_df['group']))

# Create group_id-to-info mappings
group_id_to_name = dict(zip(group_order_df['group'], group_order_df['name']))
group_id_to_position = dict(zip(group_order_df['group'], group_order_df['position']))

# Get original tasks from first dataframe
original_tasks = dfs[0]['task'].values

# Assign each task to its group ID and position
task_info = []
for task in original_tasks:
    refcode = get_refcode(task)
    group_id = task_to_group_id.get(refcode, -1)  # -1 for unknown tasks
    position = group_id_to_position.get(group_id, 999)  # Unknown groups go to end
    task_info.append((task, group_id, position))

# Sort tasks by group position, maintaining original order within groups
task_info_sorted = sorted(task_info, key=lambda x: (x[2], list(original_tasks).index(x[0])))
reordered_tasks = [t[0] for t in task_info_sorted]

# Create reordering index
reorder_idx = [list(original_tasks).index(task) for task in reordered_tasks]

# Reorder all dataframes
for i in range(len(dfs)):
    dfs[i] = dfs[i].iloc[reorder_idx].reset_index(drop=True)

tasks = reordered_tasks
n_tasks = len(tasks)
n_methods = len(dfs)

# Build score matrix with reordered data
all_scores = np.zeros((n_methods, n_tasks))
for i, df in enumerate(dfs):
    all_scores[i] = df['3DOpt_Score'].values

# Calculate group spans for visual indicators
task_groups = []
current_group_id = None
group_start = 0
for idx, task in enumerate(tasks):
    refcode = get_refcode(task)
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
    task_groups.append((group_start, len(tasks) - 1, group_name))

# ------------------------------------------------------------
# 4) Plotting setup
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
# 5) Axis styling and labels
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
# 6) Add visual group indicators above the plot
# ------------------------------------------------------------
y_pos = 1.08  # Position above the plot (in axis coordinates)
line_y = 1.04  # Y position for the bracket lines
extra_x_buffer = 0.4

for start_idx, end_idx, group_name in task_groups:
    # Calculate x positions for the group span
    x_start = start_idx + bar_width * (n_methods - 1) / 2 - extra_x_buffer
    x_end = end_idx + bar_width * (n_methods - 1) / 2 + extra_x_buffer
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

# ------------------------------------------------------------
# 7) Benchmark-score box 
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
# 8) Save the figure
# ------------------------------------------------------------
fig.tight_layout(pad=0.5)
output_path = Path("3DOpt_ChemGE_Scores.svg")
fig.savefig(output_path, format='svg')
plt.show()
print(f"✅ Plot saved to {output_path.resolve()}")
