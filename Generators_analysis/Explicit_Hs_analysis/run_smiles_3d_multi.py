"""
Run multiple SMILES→3D benchmarks and aggregate results.

- Launches N runs of run_smiles_3d_single.py with different seeds.
- Aggregates similarities and success rates for implicit vs explicit-H SMILES.
- Outputs in BASE_OUT:
    * similarity_comparison_mean.csv / std.csv
    * success_rates_per_target.csv / per_run.csv / overall_success_rates.csv
    * average_similarities_scatter_CSD.svg / _Explicit.svg
    * atom_counts_comparison_run1.csv (from run_1)

Usage:
    python run_smiles_3d_multi.py
    (set N_RUNS, BASE_OUT, SINGLE_SCRIPT at top of file)
"""

import os
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# === CONFIGURATION ===
N_RUNS = 100
BASE_OUT = "Implicit_vs_explicit_SMILES"
SINGLE_SCRIPT = "run_smiles_3d_single.py"  

TARGETS = [
    'ABAHIW', 'ABAKIZ', 'ABADOX', 'ABABIP', 'GASQOK', 'ABEKIE', 'NIWPUE01',
    'ABEKIF', 'APUFEX', 'ABEHAU', 'TITTUO', 'EGEYOG', 'ABOBUP', 'XIDTOW',
    'ACNCOB10', 'TACXUQ', 'ACAZFE', 'NIVHEJ', 'ADUPAS', 'DAJLAC', 'OFOWIS',
    'CATSUL', 'HESMUQ01', 'GUDQOL', 'ABEVAG', 'AKOQOH', 'ADARUT', 'AFECIA',
    'ACOVUL', 'AFIXEV', 'ABAYAF', 'RULJAM'
]
TARGET_LABELS = [f"{i}_{rc}" for i, rc in enumerate(TARGETS, 1)]

SOURCES = ["CSD_SMILES", "RDKit_SMILES_EXPL"]
SOURCE_LABELS = {"CSD_SMILES": "CSD", "RDKit_SMILES_EXPL": "Explicit"}
METHODS = ["ccdc", "rdkit", "obabel"]
METHOD_LABELS = {"ccdc": "CCDC", "rdkit": "RDKit", "obabel": "OBabel"}

# Plot styling
colors = ['#1f77b4', '#ff7f0e', '#2ca02c']
markers = ['o', 's', '^']


def read_similarity_table(csv_path, target_labels, ref_cols=None):
    """
    Read one run's similarity_comparison.csv as a numeric DataFrame.
    - Tries MultiIndex header; falls back to flat header.
    - Converts "N.A." to NaN via to_numeric(errors='coerce').
    - Reindexes rows to target_labels; reindexes columns to ref_cols if given.
    """
    if not os.path.exists(csv_path):
        return None
    try:
        df = pd.read_csv(csv_path, header=[0, 1], index_col=0)
    except Exception:
        df = pd.read_csv(csv_path, index_col=0)
    df = df.apply(pd.to_numeric, errors="coerce")
    df = df.reindex(target_labels)
    if ref_cols is not None:
        try:
            df = df.reindex(columns=ref_cols)
        except Exception:
            pass
    return df


def main():
    os.makedirs(BASE_OUT, exist_ok=True)

    # --- Use only RDKit-safe random seeds (int32) ---
    seeds = np.random.randint(0, 2_147_483_647, size=N_RUNS, dtype=np.int32)
    print(f"Selected seeds: {seeds.tolist()}")

    # --- 1. Run all single jobs ---
    run_dirs = []
    for i in range(1, N_RUNS + 1):
        out_dir = os.path.join(BASE_OUT, f"run_{i}")
        run_dirs.append(out_dir)
        seed = int(seeds[i - 1])
        cmd = ["python3", SINGLE_SCRIPT, "--outdir", out_dir, "--seed", str(seed)]
        print(f"Run {i}/{N_RUNS}: {' '.join(cmd)}")
        subprocess.run(cmd, check=True)

    print("All runs complete. Aggregating results...\n")

    # --- 2. Gather per-run 'similarity_comparison.csv' ---
    sim_tables = []
    present_run_dirs = []
    ref_cols = None
    for out_dir in run_dirs:
        csv_path = os.path.join(out_dir, "similarity_comparison.csv")
        df = read_similarity_table(csv_path, TARGET_LABELS, ref_cols=None)
        if df is not None:
            if ref_cols is None:
                ref_cols = df.columns  
            df = df.reindex(columns=ref_cols)
            sim_tables.append(df)
            present_run_dirs.append(out_dir)
        else:
            print(f"Warning: missing similarity_comparison.csv in {out_dir}")

    if not sim_tables:
        print("No similarity tables found. Nothing to aggregate.")
        raise SystemExit(2)

    # --- 3. Compute mean & std similarity across runs ---
    sim_stack = np.stack([df.values for df in sim_tables], axis=0)  # (runs, targets, columns)
    sim_mean = pd.DataFrame(sim_stack.mean(axis=0), index=TARGET_LABELS, columns=ref_cols)
    if sim_stack.shape[0] > 1:
        sim_std_vals = sim_stack.std(axis=0, ddof=1)
    else:
        sim_std_vals = np.zeros_like(sim_stack.mean(axis=0))
    sim_std = pd.DataFrame(sim_std_vals, index=TARGET_LABELS, columns=ref_cols)

    mean_csv = os.path.join(BASE_OUT, "similarity_comparison_mean.csv")
    std_csv  = os.path.join(BASE_OUT, "similarity_comparison_std.csv")
    sim_mean.to_csv(mean_csv)
    sim_std.to_csv(std_csv)
    print(f"Saved: {mean_csv}")
    print(f"Saved: {std_csv}")

    # --- 4. SUCCESS metrics reconstructed from similarity tables ---
    # Success = cell is not NaN
    success_tables = [(~df.isna()).astype(int) for df in sim_tables]

    # 4a) Success rate per target & (source, method): fraction of runs with non-NaN
    rows = []
    for r_idx, tgt in enumerate(TARGET_LABELS):
        row = {"Target": tgt}
        for src in SOURCES:
            for m in METHODS:
                col = (src, m) if isinstance(ref_cols, pd.MultiIndex) else f"{src}_{m}"
                try:
                    col_idx = list(ref_cols).index(col)
                except ValueError:
                    # skip if column not present in the data
                    continue
                # collect across runs
                successes = [tbl.iloc[r_idx, col_idx] for tbl in success_tables]
                n_total = len(successes)
                rate = (np.sum(successes) / n_total) if n_total else np.nan
                row[f"{SOURCE_LABELS[src]}_{m}"] = rate
        rows.append(row)
    success_per_target = pd.DataFrame(rows)
    success_per_target_csv = os.path.join(BASE_OUT, "success_rates_per_target.csv")
    success_per_target.to_csv(success_per_target_csv, index=False)
    print(f"Saved: {success_per_target_csv}")

    # 4b) Success rates per run (counts & rates)
    rows = []
    for out_dir, tbl in zip(present_run_dirs, success_tables):
        run_tag = os.path.basename(out_dir)
        row = {"run": run_tag}
        for src in SOURCES:
            for m in METHODS:
                col = (src, m) if isinstance(ref_cols, pd.MultiIndex) else f"{src}_{m}"
                if col not in tbl.columns:
                    continue
                # count across all targets for this run
                col_successes = int(tbl[col].sum())
                row[f"{SOURCE_LABELS[src]}_{m}_count"] = col_successes
                row[f"{SOURCE_LABELS[src]}_{m}_rate"] = col_successes / len(TARGETS)
        rows.append(row)
    success_per_run = pd.DataFrame(rows)
    success_per_run_csv = os.path.join(BASE_OUT, "success_rates_per_run.csv")
    success_per_run.to_csv(success_per_run_csv, index=False)
    print(f"Saved: {success_per_run_csv}")

    # 4c) Overall success across all runs/targets
    rows = []
    # total runs considered (some runs may be missing; use present_run_dirs)
    n_runs_present = len(present_run_dirs)
    for src in SOURCES:
        for m in METHODS:
            col = (src, m) if isinstance(ref_cols, pd.MultiIndex) else f"{src}_{m}"
            # sum successes across runs and targets
            n_success = 0
            n_total = 0
            for tbl in success_tables:
                if col in tbl.columns:
                    n_success += int(tbl[col].sum())
                    n_total += tbl.shape[0]  # one attempt per target
            rate = (n_success / n_total) if n_total else np.nan
            rows.append({
                "SMILES_Source": SOURCE_LABELS[src],
                "Method": METHOD_LABELS[m],
                "n_success": int(n_success),
                "n_total": int(n_total),
                "success_rate": rate
            })
    overall_success = pd.DataFrame(rows)
    overall_success_csv = os.path.join(BASE_OUT, "overall_success_rates.csv")
    overall_success.to_csv(overall_success_csv, index=False)
    print(f"Saved: {overall_success_csv}")

    # --- 5. Optional: copy atom-count comparison from run_1 ---
    run1_atom_counts = os.path.join(run_dirs[0], "atom_counts_comparison.csv")
    if os.path.exists(run1_atom_counts):
        dst = os.path.join(BASE_OUT, "atom_counts_comparison_run1.csv")
        pd.read_csv(run1_atom_counts).to_csv(dst, index=False)
        print(f"Saved: {dst}")
    else:
        print("Note: atom_counts_comparison.csv not found in run_1 (nothing copied).")

    # --- 6. Plotting from MEAN similarities ---
    mean_tidy = sim_mean.copy()
    if isinstance(mean_tidy.columns, pd.MultiIndex):
        mean_tidy.columns = mean_tidy.columns.map(lambda t: f"{t[0]}_{t[1]}")
    else:
        mean_tidy.columns = mean_tidy.columns.astype(str)

    x = np.arange(len(mean_tidy.index))

    def plot_group(source_key, title, filename, ylim=(0, 1)):
        cols = [c for c in mean_tidy.columns if c.startswith(source_key + "_")]
        ordered_cols, ordered_labels = [], []
        for m in METHODS:
            col = f"{source_key}_{m}"
            if col in cols:
                ordered_cols.append(col)
                ordered_labels.append(METHOD_LABELS.get(m, m))
        if not ordered_cols:
            print(f"Skipping plot '{title}': no columns found for {source_key}")
            return

        fig, ax = plt.subplots(figsize=(8, 4.5))
        for i, (lab, col) in enumerate(zip(ordered_labels, ordered_cols)):
            y = mean_tidy[col].values
            ax.scatter(
                x, y, label=lab,
                color=colors[i % len(colors)],
                marker=markers[i % len(markers)],
                s=60, edgecolor='black', linewidth=0.7
            )
        ax.set_ylabel("Avg Similarity", fontsize=10)
        ax.set_xlabel("Target", fontsize=10)
        src_label = SOURCE_LABELS.get(source_key, source_key)
        ax.set_title(title.replace("{SRC}", src_label), fontsize=11, fontweight='bold')
        ax.set_xticks(x)
        ax.set_xticklabels(TARGET_LABELS, rotation=45, ha="right", fontsize=8)
        ax.tick_params(axis='both', labelsize=8)
        if ylim:
            ax.set_ylim(*ylim)
        ax.legend(loc='best', fontsize=8, frameon=True, edgecolor='black')
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        fig.tight_layout(pad=0.7)
        path = os.path.join(BASE_OUT, filename)
        fig.savefig(path, format='svg')
        plt.close(fig)
        print(f"✅ Plot saved to {path}")

    plot_group("CSD_SMILES",
               title="Average 3D Similarity per Target — {SRC} SMILES",
               filename="average_similarities_scatter_CSD.svg",
               ylim=(0, 1))

    plot_group("RDKit_SMILES_EXPL",
               title="Average 3D Similarity per Target — {SRC} SMILES",
               filename="average_similarities_scatter_Explicit.svg",
               ylim=(0, 1))

    print("\n🎉 ALL DONE. Results and plots saved in", BASE_OUT)


if __name__ == "__main__":
    main()
