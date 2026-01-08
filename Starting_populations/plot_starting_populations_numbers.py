import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import textwrap
import re
import matplotlib

matplotlib.rcParams.update({
    'font.size': 10,
    'font.family': 'Arial'
})

csv_path = Path("Starting_populations.csv")
output_fig = Path("Starting_populations")

df = pd.read_csv(csv_path)
df.set_index("target", inplace=True)

numeric_cols = [c for c in df.columns if df[c].dtype.kind in "if"]
df = df[numeric_cols]

# ------------------------------------------------------------
# Read task grouping and group order files
# ------------------------------------------------------------
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

# Get original tasks
original_tasks = df.index.tolist()

# Assign each task to its group ID and position
task_info = []
for task in original_tasks:
    refcode = get_refcode(task)
    group_id = task_to_group_id.get(refcode, -1)  # -1 for unknown tasks
    position = group_id_to_position.get(group_id, 999)  # Unknown groups go to end
    task_info.append((task, group_id, position))

# Sort tasks by group position, maintaining original order within groups
task_info_sorted = sorted(task_info, key=lambda x: (x[2], original_tasks.index(x[0])))
reordered_tasks = [t[0] for t in task_info_sorted]

# Reorder dataframe
df = df.reindex(reordered_tasks)

# Calculate group spans for visual indicators
task_groups = []
current_group_id = None
group_start = 0
for idx, task in enumerate(df.index):
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
    task_groups.append((group_start, len(df) - 1, group_name))

scale_factor = 1e5
df_scaled = df / scale_factor

num_targets = len(df_scaled)
num_cols = len(df_scaled.columns)
bar_width = 0.8 / num_cols         
x = range(num_targets)              

fig, ax = plt.subplots(figsize=(6.93, 4.5))
fig.subplots_adjust(bottom=0.30)

colors = ['#1f77b4', '#ff7f0e', '#2ca02c']

for i, col in enumerate(df_scaled.columns):
    bars = ax.bar(
        [p + i * bar_width for p in x],
        df_scaled[col],
        width=bar_width,
        label=str(col),
        color=colors[i % len(colors)],
        edgecolor='black',
        linewidth=0.4
    )

ax.set_ylabel(r"Number of molecules ($\times 10^5$)", fontsize=10, fontname='Arial')
ax.set_xlabel("Target", fontsize=10, fontname='Arial')

def clean_label(t):
    return re.sub(r'^\d+_', '', t)

ax.set_xticks([p + bar_width * (num_cols - 1) / 2 for p in x])
ax.set_xticklabels(
    [textwrap.fill(clean_label(t), 12) for t in df_scaled.index],
    rotation=45,
    ha="right",
    rotation_mode="anchor",
    fontsize=8,
    fontname='Arial'
)
ax.tick_params(axis='both', labelsize=8)

ax.yaxis.set_major_formatter(
    plt.FuncFormatter(lambda x, _: f"{x:.1f}".rstrip('0').rstrip('.') if '.' in f"{x:.1f}" else f"{int(x)}")
)
for label in ax.get_yticklabels():
    label.set_fontsize(10)
    label.set_fontname('Arial')

ax.legend(
    title="Similarity threshold",
    loc='upper right', 
    fontsize=7,
    title_fontsize=7,
    prop={'family': 'Arial', 'size': 7},
)

ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)

group_width = bar_width * num_cols
buffer = bar_width * 1.2 
left = 0 - buffer
right = num_targets - 1 + group_width + buffer - bar_width
ax.set_xlim(left, right)

# ------------------------------------------------------------
# Add visual group indicators above the plot
# ------------------------------------------------------------
y_pos = 1.08  # Position above the plot (in axis coordinates)
line_y = 1.04  # Y position for the bracket lines
extra_x_buffer = 0.4

for start_idx, end_idx, group_name in task_groups:
    # Calculate x positions for the group span
    x_start = start_idx + bar_width * (num_cols - 1) / 2 - extra_x_buffer
    x_end = end_idx + bar_width * (num_cols - 1) / 2 + extra_x_buffer
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

fig.tight_layout(pad=0.5)
svg_path = output_fig.with_suffix('.svg')
fig.savefig(svg_path, format='svg')
print(f"✅  Plot saved to {svg_path.resolve()}")
