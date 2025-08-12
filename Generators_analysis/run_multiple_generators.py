#!/usr/bin/env python3

"""
Run and aggregate multiple rounds of 3D structure generation benchmarks.
- Launches N independent runs of 'run_single_generators_test.py' on a set of targets and methods.
- Aggregates per-run and per-target results (success rate, average similarity, generation time).
- Outputs summary CSVs and scatter plots for benchmarking generator performance.

Results are saved in BASE_OUT, including:
- success_rates_per_target.csv
- success_rates_per_run.csv
- overall_success_rates.csv
- average_scores_and_times.csv
- SVG plots
"""

import os
import subprocess
import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt

# === CONFIGURATION ===
N_RUNS = 10  # Adjust as needed
BASE_OUT = "3D_Generators_analysis"
SINGLE_SCRIPT = "run_single_generators.py"
TARGETS = [
    'ABAHIW', 'ABAKIZ', 'ABADOX', 'ABABIP', 'GASQOK', 'ABEKIE', 'NIWPUE01',
    'ABEKIF', 'APUFEX', 'ABEHAU', 'TITTUO', 'EGEYOG', 'ABOBUP', 'XIDTOW',
    'ACNCOB10', 'TACXUQ', 'ACAZFE', 'NIVHEJ', 'ADUPAS', 'DAJLAC', 'OFOWIS',
    'CATSUL', 'HESMUQ01', 'GUDQOL', 'ABEVAG', 'AKOQOH', 'ADARUT', 'AFECIA',
    'ACOVUL', 'AFIXEV', 'ABAYAF', 'RULJAM'
]
methods = ["ccdc", "rdkit", "obabel"]
method_labels = ['CCDC', 'RDKit', 'OBabel']
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
markers = ['o', 's', '^']

os.makedirs(BASE_OUT, exist_ok=True)

# --- Use only RDKit-safe random seeds ---
seeds = np.random.randint(0, 2_147_483_647, size=N_RUNS, dtype=np.int32)
print(f"Selected seeds: {seeds}")

# --- 1. Run all single jobs ---
run_dirs = []
for i in range(1, N_RUNS+1):
    out_dir = os.path.join(BASE_OUT, f"run_{i}")
    run_dirs.append(out_dir)
    seed = int(seeds[i-1])
    print(f"Seed for run {i}: {seed} (type: {type(seed)})")
    cmd = ["python3", SINGLE_SCRIPT, out_dir, str(seed)]
    print(f"Run {i}/{N_RUNS}: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

print("All runs complete. Aggregating results...")

# --- 2. Gather results ---
all_tables = []
for out_dir in run_dirs:
    # Look for a CSV file starting with scores_and_times
    found = glob.glob(os.path.join(out_dir, "scores_and_times*.csv"))
    if found:
        csv_path = found[0]
        df = pd.read_csv(csv_path)
        df["run"] = os.path.basename(out_dir)
        all_tables.append(df)
    else:
        print(f"Warning: missing scores_and_times*.csv in {out_dir}")

if not all_tables:
    print("No results found.")
    exit()

alldf = pd.concat(all_tables, ignore_index=True)

# --- 3a. Compute success rate per target (across runs) ---
target_rows = []
for target in TARGETS:
    row = {"Target": target}
    for m in methods:
        n_success = (alldf.loc[alldf["Target"] == target, f"{m}_ok"] == "✓").sum()
        n_total = N_RUNS
        row[m] = n_success / n_total
    target_rows.append(row)

success_target_df = pd.DataFrame(target_rows)
success_target_csv = os.path.join(BASE_OUT, "success_rates_per_target.csv")
success_target_df.to_csv(success_target_csv, index=False)
print(f"Success rates per target saved to {success_target_csv}")

# --- 3b. Compute success rates per run ---
success_rate_rows = []
for i, out_dir in enumerate(run_dirs):
    df = alldf[alldf['run'] == os.path.basename(out_dir)]
    row = {"run": os.path.basename(out_dir)}
    for m in methods:
        n_success = (df[f"{m}_ok"] == "✓").sum()
        row[f"{m}_success_count"] = n_success
        row[f"{m}_success_rate"] = n_success / len(TARGETS)
    success_rate_rows.append(row)
success_rate_df = pd.DataFrame(success_rate_rows)
success_rate_csv = os.path.join(BASE_OUT, "success_rates_per_run.csv")
success_rate_df.to_csv(success_rate_csv, index=False)
print(f"Success rates per run saved to {success_rate_csv}")

# --- 4. Compute average success rate per method over all runs ---
overall_success = {}
for m in methods:
    n_success = (alldf[f"{m}_ok"] == "✓").sum()
    n_total = len(TARGETS) * N_RUNS
    overall_success[m] = (n_success, n_total, n_success / n_total)
overall_success_df = pd.DataFrame(
    {
        "method": list(overall_success.keys()),
        "n_success": [overall_success[m][0] for m in methods],
        "n_total": [overall_success[m][1] for m in methods],
        "success_rate": [overall_success[m][2] for m in methods],
    }
)
overall_success_csv = os.path.join(BASE_OUT, "overall_success_rates.csv")
overall_success_df.to_csv(overall_success_csv, index=False)
print(f"Overall success rates saved to {overall_success_csv}")

# --- 5. Compute mean & std for time/similarity, per method & target ---
aggrows = []
for target in TARGETS:
    row = {"Target": target}
    for m in methods:
        vals_sim = pd.to_numeric(alldf.loc[alldf["Target"] == target, f"{m}_sim"], errors='coerce')
        vals_time = pd.to_numeric(alldf.loc[alldf["Target"] == target, f"{m}_time"], errors='coerce')
        row[f"{m}_sim_mean"] = vals_sim.mean()
        row[f"{m}_sim_std"] = vals_sim.std()
        row[f"{m}_time_mean"] = vals_time.mean()
        row[f"{m}_time_std"] = vals_time.std()
    aggrows.append(row)
aggdf = pd.DataFrame(aggrows)
aggdf.to_csv(os.path.join(BASE_OUT, "average_scores_and_times.csv"), index=False)
print(f"Averaged scores and times saved to {os.path.join(BASE_OUT, 'average_scores_and_times.csv')}")

# --- 6. PLOTTING ---

x = np.arange(len(TARGETS))

# --- Similarity SCATTER ---
fig, ax = plt.subplots(figsize=(7, 4))
for i, m in enumerate(methods):
    y = aggdf[f"{m}_sim_mean"].values
    ax.scatter(
        x, y, label=method_labels[i],
        color=colors[i], marker=markers[i], s=60,
        edgecolor='black', linewidth=0.7
    )
ax.set_ylabel("Avg Similarity", fontsize=10, fontname='Arial')
ax.set_xlabel("Target", fontsize=10, fontname='Arial')
ax.set_title("Average 3D Similarity per Target", fontsize=10, fontname='Arial', fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(TARGETS, rotation=45, ha="right", fontsize=8, fontname='Arial')
ax.tick_params(axis='both', labelsize=8)
ax.set_ylim(0, 1)
ax.legend(loc='lower right', fontsize=7, frameon=True, edgecolor='black')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
fig.tight_layout(pad=0.5)
sim_plot_path = os.path.join(BASE_OUT, "average_similarities_scatter.svg")
fig.savefig(sim_plot_path, format='svg')
plt.close(fig)
print(f"✅ Similarity scatter plot saved to {sim_plot_path}")

# --- Times SCATTER ---
fig, ax = plt.subplots(figsize=(7, 4))
for i, m in enumerate(methods):
    y = aggdf[f"{m}_time_mean"].values
    ax.scatter(
        x, y, label=method_labels[i],
        color=colors[i], marker=markers[i], s=60,
        edgecolor='black', linewidth=0.7
    )
ax.set_ylabel("Avg Time (s)", fontsize=10, fontname='Arial')
ax.set_xlabel("Target", fontsize=10, fontname='Arial')
ax.set_title("Average Generation Time per Target", fontsize=10, fontname='Arial', fontweight='bold')
ax.set_xticks(x)
ax.set_xticklabels(TARGETS, rotation=45, ha="right", fontsize=8, fontname='Arial')
ax.tick_params(axis='both', labelsize=8)
ax.set_ylim(bottom=0)
ax.legend(loc='upper right', fontsize=7, frameon=True, edgecolor='black')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
fig.tight_layout(pad=0.5)
time_plot_path = os.path.join(BASE_OUT, "average_times_scatter.svg")
fig.savefig(time_plot_path, format='svg')
plt.close(fig)
print(f"✅ Time scatter plot saved to {time_plot_path}")

print("\n🎉 ALL DONE. Results and plots saved in", BASE_OUT)
