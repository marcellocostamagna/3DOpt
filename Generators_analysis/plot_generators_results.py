#!/usr/bin/env python3

"""
Plots benchmarking results (success rate, similarity, runtime) for 3D molecular generators.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import textwrap
import re
from pathlib import Path

# ---- Set global font to Arial size 10 ----
import matplotlib
matplotlib.rcParams.update({
    'font.size': 10,
    'font.family': 'Arial'
})

RESULTS_DIR = Path("3D_Generators_analysis")  

success_csv = RESULTS_DIR / "success_rates_per_target.csv"
agg_csv     = RESULTS_DIR / "average_scores_and_times.csv"

PLOT_OUT = Path("./Generators_plots") 
PLOT_OUT.mkdir(exist_ok=True)

method_labels = ['CCDC', 'OBabel', 'RDKit']
csv_methods   = ['ccdc', 'obabel', 'rdkit']
colors        = ['#1f77b4', '#ff7f0e', '#2ca02c'] 

markers = ['o', 's', '^']

hist_df = pd.read_csv(success_csv)
aggdf = pd.read_csv(agg_csv)
targets = aggdf['Target'].tolist()
hist_df = hist_df.set_index('Target').reindex(targets).reset_index()

# ------------------------------------------------------------
# Read task grouping and group order files
# ------------------------------------------------------------
task_grouping_df = pd.read_csv("../task_grouping/task_grouping")
group_order_df = pd.read_csv("../task_grouping/group_order")

# Create task-to-group_id mapping (task -> integer group ID)
task_to_group_id = dict(zip(task_grouping_df['task'], task_grouping_df['group']))

# Create group_id-to-info mappings
group_id_to_name = dict(zip(group_order_df['group'], group_order_df['name']))
group_id_to_position = dict(zip(group_order_df['group'], group_order_df['position']))

# Get original tasks
original_tasks = targets

# Assign each task to its group ID and position
task_info = []
for task in original_tasks:
    group_id = task_to_group_id.get(task, -1)  # -1 for unknown tasks
    position = group_id_to_position.get(group_id, 999)  # Unknown groups go to end
    task_info.append((task, group_id, position))

# Sort tasks by group position, maintaining original order within groups
task_info_sorted = sorted(task_info, key=lambda x: (x[2], list(original_tasks).index(x[0])))
reordered_tasks = [t[0] for t in task_info_sorted]

# Create reordering index
reorder_idx = [list(original_tasks).index(task) for task in reordered_tasks]

# Reorder dataframes
aggdf = aggdf.iloc[reorder_idx].reset_index(drop=True)
hist_df = hist_df.iloc[reorder_idx].reset_index(drop=True)

targets = reordered_tasks
n_tasks = len(targets)
n_methods = len(method_labels)
x = np.arange(n_tasks)

# Calculate group spans for visual indicators
task_groups = []
current_group_id = None
group_start = 0
for idx, task in enumerate(targets):
    group_id = task_to_group_id.get(task, -1)
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
    task_groups.append((group_start, len(targets) - 1, group_name))

def clean_label(t):
    return re.sub(r'^\d+_', '', t)

bar_width = 0.7 / n_methods
group_width = bar_width * n_methods
buffer = bar_width * 1.2
left = 0 - buffer
right = n_tasks - 1 + group_width + buffer - bar_width

def plot_success_bar(ax, x, hist_df, bar_width, method_labels, csv_methods, colors, targets, show_xlabel=True, show_legend=True, add_group_indicators=False):
    for i in range(len(csv_methods)):
        ax.bar(
            [p + i * bar_width for p in x],
            hist_df[csv_methods[i]].values,
            width=bar_width,
            label=method_labels[i],
            color=colors[i],
            edgecolor='black',
            linewidth=0.4
        )
    ax.set_ylabel("Success Rate", fontsize=10, fontname='Arial', labelpad=3)
    if show_xlabel:
        ax.set_xlabel("Task", fontsize=10, fontname='Arial', labelpad=3)
        ax.set_xticks([p + bar_width for p in x])
        ax.set_xticklabels(
            [textwrap.fill(clean_label(t), 12) for t in targets],
            rotation=45, ha="right", rotation_mode="anchor",
            fontsize=8, fontname='Arial'
        )
    else:
        ax.set_xticks([])
        ax.set_xticklabels([])
    ax.tick_params(axis='both', labelsize=8)
    ax.set_ylim(0, 1.08)
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: f"{x:.1f}".rstrip('0').rstrip('.') if '.' in f"{x:.1f}" else f"{int(x)}")
    )
    for label in ax.get_yticklabels():
        label.set_fontsize(10)
        label.set_fontname('Arial')
    if show_legend:
        legend = ax.legend(
            loc='lower right',
            fontsize=7, prop={'family': 'Arial', 'size': 7},
            framealpha=1.0, frameon=True, edgecolor='black'
        )
        legend.get_frame().set_linewidth(0.4)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(left, right)
    
    # Add group indicators if requested
    if add_group_indicators:
        y_pos = 1.08
        line_y = 1.04
        extra_x_buffer = 0.4
        for start_idx, end_idx, group_name in task_groups:
            x_start = start_idx + bar_width - extra_x_buffer
            x_end = end_idx + bar_width + extra_x_buffer
            x_center = (x_start + x_end) / 2
            ax.plot([x_start, x_end], [line_y, line_y], 
                    color='black', linewidth=0.8, 
                    transform=ax.get_xaxis_transform(), clip_on=False)
            tick_height = 0.015
            ax.plot([x_start, x_start], [line_y - tick_height, line_y + tick_height], 
                    color='black', linewidth=0.8, 
                    transform=ax.get_xaxis_transform(), clip_on=False)
            ax.plot([x_end, x_end], [line_y - tick_height, line_y + tick_height], 
                    color='black', linewidth=0.8, 
                    transform=ax.get_xaxis_transform(), clip_on=False)
            ax.text(x_center, y_pos, group_name, 
                    ha='center', va='bottom', fontsize=9, fontname='Arial',
                    transform=ax.get_xaxis_transform())

def plot_similarity_scatter(ax, x, aggdf, method_labels, csv_methods, colors, markers, targets, show_xlabel=True, show_legend=True, add_group_indicators=False, x_offset=0):
    handles = []
    for i, m in enumerate(csv_methods):
        y = aggdf[f"{m}_sim_mean"].values
        sc = ax.scatter(
            x + x_offset, y, label=method_labels[i],
            color=colors[i], marker=markers[i], s=25,
            edgecolor='black', linewidth=0.5
        )
        handles.append(sc)
    ax.set_ylabel("HSR Similarity", fontsize=10, fontname='Arial', labelpad=3)
    if show_xlabel:
        ax.set_xlabel("Task", fontsize=10, fontname='Arial', labelpad=3)
        ax.set_xticks(x)
        ax.set_xticklabels(
            [textwrap.fill(clean_label(t), 12) for t in targets],
            rotation=45, ha="right", rotation_mode="anchor",
            fontsize=8, fontname='Arial'
        )
    else:
        ax.set_xticks([])
        ax.set_xticklabels([])
    ax.tick_params(axis='both', labelsize=8)
    ax.set_ylim(0, 1)
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: f"{x:.1f}" if isinstance(x, float) else str(int(x)))
    )
    for label in ax.get_yticklabels():
        label.set_fontsize(10)
        label.set_fontname('Arial')
    if show_legend:
        legend = ax.legend(
            handles, method_labels,
            loc='lower right',
            fontsize=7, prop={'family': 'Arial', 'size': 7},
            framealpha=1.0, frameon=True, edgecolor='black'
        )
        legend.get_frame().set_linewidth(0.4)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(left, right)
    
    # Add group indicators if requested
    if add_group_indicators:
        y_pos = 1.08
        line_y = 1.04
        extra_x_buffer = 0.4
        for start_idx, end_idx, group_name in task_groups:
            x_start = start_idx + x_offset - extra_x_buffer
            x_end = end_idx + x_offset + extra_x_buffer
            x_center = (x_start + x_end) / 2
            ax.plot([x_start, x_end], [line_y, line_y], 
                    color='black', linewidth=0.8, 
                    transform=ax.get_xaxis_transform(), clip_on=False)
            tick_height = 0.015
            ax.plot([x_start, x_start], [line_y - tick_height, line_y + tick_height], 
                    color='black', linewidth=0.8, 
                    transform=ax.get_xaxis_transform(), clip_on=False)
            ax.plot([x_end, x_end], [line_y - tick_height, line_y + tick_height], 
                    color='black', linewidth=0.8, 
                    transform=ax.get_xaxis_transform(), clip_on=False)
            ax.text(x_center, y_pos, group_name, 
                    ha='center', va='bottom', fontsize=9, fontname='Arial',
                    transform=ax.get_xaxis_transform())

def plot_time_scatter(ax, x, aggdf, method_labels, csv_methods, colors, markers, targets, show_xlabel=True, show_legend=True, add_group_indicators=False, x_offset=0):
    handles = []
    for i, m in enumerate(csv_methods):
        y = aggdf[f"{m}_time_mean"].values
        sc = ax.scatter(
            x + x_offset, y, label=method_labels[i],
            color=colors[i], marker=markers[i], s=25,
            edgecolor='black', linewidth=0.5
        )
        handles.append(sc)
    ax.set_ylabel("Runtime (s)", fontsize=10, fontname='Arial', labelpad=3)
    if show_xlabel:
        ax.set_xlabel("Task", fontsize=10, fontname='Arial', labelpad=3)
        ax.set_xticks(x)
        ax.set_xticklabels(
            [textwrap.fill(clean_label(t), 12) for t in targets],
            rotation=45, ha="right", rotation_mode="anchor",
            fontsize=8, fontname='Arial'
        )
    else:
        ax.set_xticks([])
        ax.set_xticklabels([])
    ax.tick_params(axis='both', labelsize=8)
    ax.set_ylim(bottom=0)
    ax.yaxis.set_major_formatter(
        plt.FuncFormatter(lambda x, _: f"{x:.1f}" if isinstance(x, float) else str(int(x)))
    )
    for label in ax.get_yticklabels():
        label.set_fontsize(10)
        label.set_fontname('Arial')
    if show_legend:
        legend = ax.legend(
            handles, method_labels,
            loc='upper left',
            fontsize=7, prop={'family': 'Arial', 'size': 7},
            framealpha=1.0, frameon=True, edgecolor='black'
        )
        legend.get_frame().set_linewidth(0.4)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_xlim(left, right)
    
    # Add group indicators if requested
    if add_group_indicators:
        y_pos = 1.08
        line_y = 1.04
        extra_x_buffer = 0.4
        for start_idx, end_idx, group_name in task_groups:
            x_start = start_idx + x_offset - extra_x_buffer
            x_end = end_idx + x_offset + extra_x_buffer
            x_center = (x_start + x_end) / 2
            ax.plot([x_start, x_end], [line_y, line_y], 
                    color='black', linewidth=0.8, 
                    transform=ax.get_xaxis_transform(), clip_on=False)
            tick_height = 0.015
            ax.plot([x_start, x_start], [line_y - tick_height, line_y + tick_height], 
                    color='black', linewidth=0.8, 
                    transform=ax.get_xaxis_transform(), clip_on=False)
            ax.plot([x_end, x_end], [line_y - tick_height, line_y + tick_height], 
                    color='black', linewidth=0.8, 
                    transform=ax.get_xaxis_transform(), clip_on=False)
            ax.text(x_center, y_pos, group_name, 
                    ha='center', va='bottom', fontsize=9, fontname='Arial',
                    transform=ax.get_xaxis_transform())

# --- 1. Success rate grouped bar plot ---
fig, ax = plt.subplots(figsize=(6.93, 2))
fig.subplots_adjust(left=0.09, right=0.99, bottom=0.30)
plot_success_bar(ax, x, hist_df, bar_width, method_labels, csv_methods, colors, targets, show_xlabel=True, show_legend=True, add_group_indicators=True)
fig.tight_layout(pad=0.5)
fig.savefig(PLOT_OUT / "success_rate_grouped_bar.svg", format='svg')
plt.close(fig)
print(f"✅ Success rate plot saved to {PLOT_OUT / 'success_rate_grouped_bar.svg'}")

# --- 2. Similarity scatter plot ---
fig, ax = plt.subplots(figsize=(6.93, 2))
fig.subplots_adjust(left=0.09, right=0.99, bottom=0.30)
plot_similarity_scatter(ax, x, aggdf, method_labels, csv_methods, colors, markers, targets, show_xlabel=True, show_legend=True, add_group_indicators=True)
fig.tight_layout(pad=0.5)
fig.savefig(PLOT_OUT / "similarity_scatter_grouped.svg", format='svg')
plt.close(fig)
print(f"✅ Similarity scatter plot saved to {PLOT_OUT / 'similarity_scatter_grouped.svg'}")

# --- 3. Generation time scatter plot ---
fig, ax = plt.subplots(figsize=(6.93, 2))
fig.subplots_adjust(left=0.09, right=0.99, bottom=0.30)
plot_time_scatter(ax, x, aggdf, method_labels, csv_methods, colors, markers, targets, show_xlabel=True, show_legend=True, add_group_indicators=True)
fig.tight_layout(pad=0.5)
fig.savefig(PLOT_OUT / "time_scatter_grouped.svg", format='svg')
plt.close(fig)
print(f"✅ Generation time scatter plot saved to {PLOT_OUT / 'time_scatter_grouped.svg'}")

# --- 4. Stacked panel plot ---
fig, axes = plt.subplots(
    nrows=3, ncols=1, sharex=True,
    figsize=(6.93, 6),
    gridspec_kw=dict(hspace=0.07, bottom=0.19, top=0.94, left=0.09, right=0.99)
)
# Top: Similarity (no xlabels, legend in lower right, with group indicators at top)
# Use x_offset=bar_width to align scatter points with bar chart ticks
plot_similarity_scatter(axes[0], x, aggdf, method_labels, csv_methods, colors, markers, targets, show_xlabel=False, show_legend=True, add_group_indicators=True, x_offset=bar_width)
# Middle: Times (no xlabels, legend in upper left, no group indicators)
# Use x_offset=bar_width to align scatter points with bar chart ticks
plot_time_scatter(axes[1], x, aggdf, method_labels, csv_methods, colors, markers, targets, show_xlabel=False, show_legend=True, add_group_indicators=False, x_offset=bar_width)
# Bottom: Success rates (xlabels, legend in lower right, no group indicators)
plot_success_bar(axes[2], x, hist_df, bar_width, method_labels, csv_methods, colors, targets, show_xlabel=True, show_legend=True, add_group_indicators=False)
fig.align_ylabels(axes)
fig.tight_layout(pad=0.5)
fig.savefig(PLOT_OUT / "stacked_generators_panel.svg", format='svg', bbox_inches='tight')
plt.close(fig)
print(f"✅ Stacked panel plot saved to {PLOT_OUT / 'stacked_generators_panel.svg'}")
