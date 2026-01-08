import os
import argparse
import pandas as pd

def find_benchmark_folders(top_folder):
    """Finds folders starting with '3DOpt' exactly two levels below top_folder."""
    result = []
    for first in os.listdir(top_folder):
        level1 = os.path.join(top_folder, first)
        if not os.path.isdir(level1):
            continue
        for second in os.listdir(level1):
            level2 = os.path.join(level1, second)
            if os.path.isdir(level2) and os.path.basename(level2).startswith("3DOpt"):
                result.append(level2)
    return result

def analyze_similarity_scores(top_folder, print_per_task=False):
    run_results = {}
    overall = {"zero": 0, "expected": 0}
    per_run_success_rates = []

    benchmark_dirs = find_benchmark_folders(top_folder)
    if not benchmark_dirs:
        print("⚠️ No '3DOpt*' folders found in the structure.")
        return

    for bench_path in benchmark_dirs:
        run_name = os.path.basename(bench_path)
        run_stats = {}
        task_success_rates = []
        per_task_output = []

        for task_folder in os.listdir(bench_path):
            task_path = os.path.join(bench_path, task_folder)
            if not os.path.isdir(task_path):
                continue

            csv_file = os.path.join(task_path, "scores.csv")
            if not os.path.exists(csv_file):
                continue

            try:
                df = pd.read_csv(csv_file)
            except Exception:
                continue

            hsr_cols = [c for c in df.columns if c.endswith("_HSR_score")]
            if not hsr_cols:
                continue

            col = hsr_cols[0]
            zero_count = int((df[col] == 0).sum())
            expected = len(df)
            if expected > 0:
                success_rate = (expected - zero_count) / expected
                task_success_rates.append(success_rate)
            else:
                success_rate = float('nan')

            run_stats[task_folder] = {
                "zero": zero_count,
                "expected": expected,
                "success_rate": success_rate
            }

            overall["zero"] += zero_count
            overall["expected"] += expected

            if print_per_task:
                per_task_output.append(
                    f"Task: {task_folder}\n"
                    f"  • {zero_count} zero‐score entries\n"
                    f"  • {expected} expected entries\n"
                    f"  • Success rate = {success_rate:.1%}"
                )

        run_results[run_name or "single_run"] = run_stats
        if task_success_rates:
            mean_run_success = sum(task_success_rates) / len(task_success_rates)
            per_run_success_rates.append(mean_run_success)
            if print_per_task and per_task_output:
                print(f"\n📦 Summary for run: {run_name}")
                print('\n'.join(per_task_output))
        else:
            per_run_success_rates.append(float('nan'))

    # Per-run stats
    print("\n=== Per-run average success rates ===")
    for i, rate in enumerate(per_run_success_rates):
        print(f"Run {i+1}: {rate:.2%}" if not pd.isnull(rate) else f"Run {i+1}: N/A")

    # Summary statistics
    rates_clean = [r for r in per_run_success_rates if not pd.isnull(r)]
    if rates_clean:
        mean_rate = sum(rates_clean) / len(rates_clean)
        min_rate = min(rates_clean)
        max_rate = max(rates_clean)
        print(f"\nMean of per-run success rates: {mean_rate:.2%}")
        print(f"Min of per-run success rates: {min_rate:.2%}")
        print(f"Max of per-run success rates: {max_rate:.2%}")
    else:
        print("\nNo valid runs to summarize.")

    # Overall summary
    print("\n=== ✅ Overall Summary ===")
    if overall["expected"] > 0:
        rate = (overall["expected"] - overall["zero"]) / overall["expected"]
        print(f"Total expected: {overall['expected']}")
        print(f"Total zero scores: {overall['zero']}")
        print(f"Overall success rate = {rate:.1%}")
    else:
        print("⚠️ No valid tasks found. Check folder structure or filenames.")

    return per_run_success_rates, (mean_rate if rates_clean else None), (min_rate if rates_clean else None), (max_rate if rates_clean else None)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyze ChemGE 3DOpt results. Computes per-run and overall success rates. Optionally print results per task."
    )
    parser.add_argument("folder", help="Folder containing ChemGE runs to analyze (e.g. Results_ChemGE_ccdc)")
    parser.add_argument("--per-task", action="store_true", default=False,
                        help="Print per-task results for each run")
    args = parser.parse_args()
    analyze_similarity_scores(args.folder, print_per_task=args.per_task)
