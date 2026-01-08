#!/usr/bin/env python3
import os
import re
import glob
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.ticker import ScalarFormatter
from matplotlib.ticker import FuncFormatter

# ---------------- CONFIG ----------------
INPUT_FOLDERS = [
    "Starting_populations_0_3",
    "Starting_populations_0_4",
    "Starting_populations_0_5",
]
PER_TARGET_BINS = 30
COMBINED_BINS = 50
PARENT_OUT = Path("Similarity_distributions")
COLORS = ['#1f77b4', '#ff7f0e', '#2ca02c'] 
# ----------------------------------------

# ---- Match the plotting style you like ----
import matplotlib
matplotlib.rcParams.update({
    'font.size': 10,
    'font.family': 'Arial'
})

def infer_suffix_from_dir(path: str):
    """Extract trailing underscore style threshold from folder name: 'Starting_populations_0_5' -> '0_5'."""
    name = os.path.basename(os.path.normpath(path))
    m = re.search(r'(\d+_\d+)$', name)
    return m.group(1) if m else None

def infer_float_threshold(suffix: str):
    """Convert '0_5' -> 0.5 (float)."""
    try:
        whole, frac = suffix.split('_')
        return float(f"{whole}.{frac}")
    except Exception:
        return None

def parse_filename(file_path: str):
    """
    Expect filenames like '1_ABAHIW_init_pop.txt' or '12_RULJAM_init_pop.txt'.
    Returns (index:int, target:str, base_stem:str without extension), e.g. (1,'ABAHIW','1_ABAHIW').
    """
    stem = os.path.splitext(os.path.basename(file_path))[0]
    m = re.match(r'^(\d+)_([^_]+)_init_pop$', stem)
    if m:
        idx = int(m.group(1))
        tgt = m.group(2)
        return idx, tgt, f"{idx}_{tgt}"
    m2 = re.match(r'^(\d+)_([^_]+)', stem)
    if m2:
        idx = int(m2.group(1))
        tgt = m2.group(2)
        return idx, tgt, f"{idx}_{tgt}"
    return 9999, stem, stem

def load_similarities(file_path: str) -> pd.DataFrame:
    """
    Reads the last whitespace-separated token as similarity.
    Assumes lines like: <id> <smiles> <similarity>
    Returns DataFrame with a 'similarity' column.
    """
    try:
        df = pd.read_csv(file_path, sep=r"\s+", header=None, usecols=[2],
                         names=["similarity"], engine="python")
        return df
    except Exception:
        sims = []
        with open(file_path, "r", encoding="utf-8") as f:
            for line in f:
                parts = line.split()
                if not parts:
                    continue
                try:
                    sims.append(float(parts[-1]))
                except Exception:
                    pass
        return pd.DataFrame({"similarity": sims})

def style_axes(ax):
    """Apply the same visual style as your other plotting script."""
    ax.grid(axis='y', linestyle='--', linewidth=0.4, alpha=0.7)
    ax.set_axisbelow(True)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)

def process_folder(in_dir: str) -> dict:
    """
    Process one threshold folder:
      - write outputs in PARENT_OUT / (Similarity_distributions_<suffix>)
      - return a dict with metadata and combined similarities for panel plot
    """
    files = sorted(glob.glob(os.path.join(in_dir, "*.txt")))
    if not files:
        print(f"⚠️  No .txt files found in {in_dir}")
        return {}

    suffix = infer_suffix_from_dir(in_dir)  # e.g. '0_5'
    sub_out = PARENT_OUT / (f"Similarity_distributions_{suffix}" if suffix
                            else f"Similarity_distributions_{Path(in_dir).name}")
    sub_out.mkdir(parents=True, exist_ok=True)

    sim_cap = infer_float_threshold(suffix) if suffix else None

    all_rows = []
    for fp in files:
        idx, tgt, base = parse_filename(fp)
        df = load_similarities(fp)
        if df.empty:
            continue
        df["index"] = idx
        df["target"] = tgt
        all_rows.append(df)

        # Per-target histogram
        hist_range = (0.0, sim_cap) if sim_cap is not None else (
            max(0.0, float(df["similarity"].min())),
            float(df["similarity"].max())
        )
        plt.figure(figsize=(6, 4))
        plt.hist(df["similarity"], bins=PER_TARGET_BINS, range=hist_range, edgecolor="black")
        plt.title(f"Similarity distribution — {base}")
        plt.xlabel("HSR similarity")
        plt.ylabel("Count")
        style_axes(plt.gca())
        plt.tight_layout()
        plt.savefig(sub_out / f"{base}.svg", format="svg")
        plt.close()

    if not all_rows:
        print(f"⚠️  No similarity data parsed for {in_dir}.")
        return {}

    combined = pd.concat(all_rows, ignore_index=True)

    # Summary stats per target (sorted by 'index')
    summary = (combined
               .groupby(["index", "target"], as_index=False)["similarity"]
               .agg(count="count", mean="mean", median="median", std="std", min="min", max="max")
               .sort_values(["index", "target"]))
    summary.to_csv(sub_out / "summary_stats.csv", index=False)

    # Save *all* similarities as CSV for panel plotting
    combined.to_csv(sub_out / "all_targets.csv", index=False)

    # Combined histogram (for per-threshold figure)
    if sim_cap is not None:
        agg_range = (0.0, sim_cap)
        cap_label = f"{sim_cap:.2f}"
    else:
        agg_range = (float(combined["similarity"].min()), float(combined["similarity"].max()))
        cap_label = "auto"

    plt.figure(figsize=(6, 4))
    plt.hist(combined["similarity"], bins=COMBINED_BINS, range=agg_range, edgecolor="black")
    plt.title(f"All targets — similarity distribution (cap={cap_label})")
    plt.xlabel("HSR similarity")
    plt.ylabel("Count")
    style_axes(plt.gca())
    plt.tight_layout()
    plt.savefig(sub_out / "all_targets.svg", format="svg")
    plt.close()

    return {
        "suffix": suffix,
        "sim_cap": sim_cap,
        "out_dir": sub_out,
        "combined_csv": sub_out / "all_targets.csv",
        "combined_values": combined["similarity"].to_numpy(),
    }

def panel_plot(results: list):
    """
    Build a single horizontal 3-panel histogram figure using each threshold's all-targets CSV.
    Saves to PARENT_OUT / 'similarity_distributions_panel.svg'
    """
    # Load from CSVs to honor “use the all targets csv”
    datasets, labels, ranges = [], [], []
    for res in results:
        if not res:
            continue
        df = pd.read_csv(res["combined_csv"])
        vals = df["similarity"].to_numpy()
        datasets.append(vals)
        labels.append(res["suffix"])  # e.g., '0_3'
        cap = res["sim_cap"]
        ranges.append((0.0, cap) if cap is not None else (float(vals.min()), float(vals.max())))

    if not datasets:
        print("⚠️  No datasets for panel plot.")
        return

    # Common y-limit based on *raw* counts (no scaling of data)
    ymax = 0
    for vals, rng in zip(datasets, ranges):
        counts, _ = np.histogram(vals, bins=COMBINED_BINS, range=rng)
        ymax = max(ymax, counts.max())
    ymax *= 1.10  # headroom

    fig, axes = plt.subplots(
        nrows=1, ncols=len(datasets),
        figsize=(6.93, 2.50),
        sharey=True, sharex=True
    )
    if len(datasets) == 1:
        axes = [axes]

    # Tick formatter: display counts in ×10^5 units, without changing the underlying data
    to_1e5 = FuncFormatter(lambda y, _:
                           f"{(y/1e5):.1f}".rstrip('0').rstrip('.'))

    for idx, (ax, vals, rng, lab) in enumerate(zip(axes, datasets, ranges, labels)):
        ax.hist(
            vals,
            bins=COMBINED_BINS,
            range=rng,
            color=COLORS[idx % len(COLORS)],
            edgecolor="black"
        )
        ax.set_ylim(0, ymax)

        # Style (grid, spines, etc.)
        style_axes(ax)

        # X-axis fixed: 0.0 to 0.5 in 0.1 steps
        ax.set_xlim(0.0, 0.5)
        ax.set_xticks(np.arange(0.0, 0.51, 0.1))
        ax.xaxis.set_major_formatter(FuncFormatter(lambda x, _: f"{x:.1f}"))

        # Y ticks shown in 1e5 units; hide any scientific-notation offset text
        ax.yaxis.set_major_formatter(to_1e5)
        ax.yaxis.get_offset_text().set_visible(False)

        # No per-axes x-labels or titles
        ax.set_xlabel("")
        ax.set_title("")

        # Threshold badge (slightly right & down)
        cap = rng[1]
        badge = f"threshold = {cap:.1f}" if cap is not None else "threshold = auto"
        ax.text(
            0.12, 0.88, badge,
            transform=ax.transAxes, ha='left', va='top',
            fontsize=8, fontname='Arial',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor='black', linewidth=0.4)
        )

    # Shared labels
    axes[0].set_ylabel(r"Count ($\times 10^5$)", fontsize=10, fontname='Arial')
    fig.supxlabel("HSR similarity", fontsize=10, fontname='Arial')

    fig.tight_layout(pad=0.5)
    out_path = PARENT_OUT / "similarity_distributions_panel.svg"
    fig.savefig(out_path, format="svg")
    plt.close(fig)
    print(f"✅ Panel figure saved to {out_path.resolve()}")

def main():
    PARENT_OUT.mkdir(parents=True, exist_ok=True)

    results = []
    for folder in INPUT_FOLDERS:
        print(f"Processing: {folder}")
        res = process_folder(folder)
        if res:
            results.append(res)

    # Build the 3-panel figure from saved all_targets.csv files
    panel_plot(results)

if __name__ == "__main__":
    main()
