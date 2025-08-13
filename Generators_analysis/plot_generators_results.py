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
n_tasks = len(targets)
n_methods = len(method_labels)
x = np.arange(n_tasks)

def clean_label(t):
    return re.sub(r'^\d+_', '', t)

bar_width = 0.7 / n_methods
group_width = bar_width * n_methods
buffer = bar_width * 1.2
left = 0 - buffer
right = n_tasks - 1 + group_width + buffer - bar_width

def plot_success_bar(ax, x, hist_df, bar_width, method_labels, csv_methods, colors, targets, show_xlabel=True, show_legend=True):
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

def plot_similarity_scatter(ax, x, aggdf, method_labels, csv_methods, colors, markers, targets, show_xlabel=True, show_legend=True):
    handles = []
    for i, m in enumerate(csv_methods):
        y = aggdf[f"{m}_sim_mean"].values
        sc = ax.scatter(
            x, y, label=method_labels[i],
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

def plot_time_scatter(ax, x, aggdf, method_labels, csv_methods, colors, markers, targets, show_xlabel=True, show_legend=True):
    handles = []
    for i, m in enumerate(csv_methods):
        y = aggdf[f"{m}_time_mean"].values
        sc = ax.scatter(
            x, y, label=method_labels[i],
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

# --- 1. Success rate grouped bar plot ---
fig, ax = plt.subplots(figsize=(6.93, 2))
fig.subplots_adjust(left=0.09, right=0.99, bottom=0.30)
plot_success_bar(ax, x, hist_df, bar_width, method_labels, csv_methods, colors, targets, show_xlabel=True, show_legend=True)
fig.tight_layout(pad=0.5)
fig.savefig(PLOT_OUT / "success_rate_grouped_bar.svg", format='svg')
plt.close(fig)
print(f"✅ Success rate plot saved to {PLOT_OUT / 'success_rate_grouped_bar.svg'}")

# --- 2. Similarity scatter plot ---
fig, ax = plt.subplots(figsize=(6.93, 2))
fig.subplots_adjust(left=0.09, right=0.99, bottom=0.30)
plot_similarity_scatter(ax, x, aggdf, method_labels, csv_methods, colors, markers, targets, show_xlabel=True, show_legend=True)
fig.tight_layout(pad=0.5)
fig.savefig(PLOT_OUT / "similarity_scatter_grouped.svg", format='svg')
plt.close(fig)
print(f"✅ Similarity scatter plot saved to {PLOT_OUT / 'similarity_scatter_grouped.svg'}")

# --- 3. Generation time scatter plot ---
fig, ax = plt.subplots(figsize=(6.93, 2))
fig.subplots_adjust(left=0.09, right=0.99, bottom=0.30)
plot_time_scatter(ax, x, aggdf, method_labels, csv_methods, colors, markers, targets, show_xlabel=True, show_legend=True)
fig.tight_layout(pad=0.5)
fig.savefig(PLOT_OUT / "time_scatter_grouped.svg", format='svg')
plt.close(fig)
print(f"✅ Generation time scatter plot saved to {PLOT_OUT / 'time_scatter_grouped.svg'}")

# --- 4. Stacked panel plot ---
fig, axes = plt.subplots(
    nrows=3, ncols=1, sharex=True,
    figsize=(6.93, 6),
    gridspec_kw=dict(hspace=0.07, bottom=0.19, top=0.99, left=0.09, right=0.99)
)
# Top: Similarity (no xlabels, legend in lower right)
plot_similarity_scatter(axes[0], x, aggdf, method_labels, csv_methods, colors, markers, targets, show_xlabel=False, show_legend=True)
# Middle: Times (no xlabels, legend in upper left)
plot_time_scatter(axes[1], x, aggdf, method_labels, csv_methods, colors, markers, targets, show_xlabel=False, show_legend=True)
# Bottom: Histogram (xlabels, legend in lower right)
plot_success_bar(axes[2], x, hist_df, bar_width, method_labels, csv_methods, colors, targets, show_xlabel=True, show_legend=True)
fig.align_ylabels(axes)
fig.tight_layout(pad=0.5)
fig.savefig(PLOT_OUT / "stacked_generators_panel.svg", format='svg')
plt.close(fig)
print(f"✅ Stacked panel plot saved to {PLOT_OUT / 'stacked_generators_panel.svg'}")
