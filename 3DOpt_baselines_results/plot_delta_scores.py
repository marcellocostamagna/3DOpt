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

# ------------------------------------------------------------
# Read task grouping and group order files
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

# Get original tasks
original_tasks = tasks

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

# Reorder all pairs dataframes
for gen in pairs:
    pairs[gen] = pairs[gen].iloc[reorder_idx].reset_index(drop=True)

tasks = reordered_tasks

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
ax.set_ylabel(r"$\Delta$$_{Task\ Score}$ (ChemGE − RanSam)", 
              fontsize=8, fontname='Arial')


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

# ------------------------------------------------------------
# Add visual group indicators above the plot
# ------------------------------------------------------------
y_pos = 1.08  # Position above the plot (in axis coordinates)
line_y = 1.04  # Y position for the bracket lines
extra_x_buffer = 0.4

for start_idx, end_idx, group_name in task_groups:
    # Calculate x positions for the group span
    x_start = start_idx + bar_width * (n_g - 1) / 2 - extra_x_buffer
    x_end = end_idx + bar_width * (n_g - 1) / 2 + extra_x_buffer
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

# Save
outpath = Path("Delta_3DOpt_scores.svg")
fig.tight_layout(pad=0.5)
fig.savefig(outpath, format='svg')
plt.show()
print(f"✅ Saved {outpath.resolve()}")
