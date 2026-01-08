import os
import argparse
import pandas as pd

def analyze_results(top_folder, print_per_task=False):
    run_results = {}
    overall = {"sdf": 0, "zero": 0, "expected": 0}
    per_run_success_rates = []

    for first_level in os.listdir(top_folder):
        level1_path = os.path.join(top_folder, first_level)
        if not os.path.isdir(level1_path):
            continue

        for run_dir in os.listdir(level1_path):
            run_path = os.path.join(level1_path, run_dir)
            if not run_path.startswith(os.path.join(level1_path, "3DOpt")):
                continue
            if not os.path.isdir(run_path):
                continue

            run_stats = {}
            task_success_rates = []
            any_flag = False  # Did we flag any weird cases for this run?

            # Store per-task output if requested
            per_task_output = []

            for task_folder in os.listdir(run_path):
                task_path = os.path.join(run_path, task_folder)
                if not os.path.isdir(task_path):
                    continue

                csv_file = os.path.join(task_path, "scores.csv")
                if not os.path.exists(csv_file):
                    continue

                # Look for the subfolder ending in _HSR
                sdf_subdirs = [
                    d for d in os.listdir(task_path)
                    if os.path.isdir(os.path.join(task_path, d)) and d.endswith("_HSR")
                ]

                if len(sdf_subdirs) != 1:
                    continue

                sdf_dir = os.path.join(task_path, sdf_subdirs[0])

                try:
                    df = pd.read_csv(csv_file)
                except Exception:
                    continue

                hsr_cols = [c for c in df.columns if c.endswith("_HSR_score")]
                if not hsr_cols:
                    continue

                col = hsr_cols[0]
                sdf_count = 0
                for root, dirs, files in os.walk(sdf_dir):
                    sdf_count += sum(f.endswith(".sdf") for f in files)

                expected = len(df)
                zero_count = int((df[col] == 0).sum())

                # Consistency check
                remainder = expected - sdf_count - zero_count
                if remainder != 0:
                    print(f"⚠️ WARNING: In run '{run_dir}', task '{task_folder}': "
                          f"expected ({expected}) != sdf_count ({sdf_count}) + zero_count ({zero_count}) (remainder={remainder})")
                    any_flag = True

                # Success rate formula (updated)
                if expected > 0:
                    success_rate = sdf_count / expected
                    task_success_rates.append(success_rate)
                else:
                    success_rate = float('nan')

                run_stats[task_folder] = {
                    "sdf": sdf_count,
                    "zero": zero_count,
                    "expected": expected,
                    "success_rate": success_rate
                }

                overall["sdf"] += sdf_count
                overall["zero"] += zero_count
                overall["expected"] += expected

                if print_per_task:
                    per_task_output.append(
                        f"Task: {task_folder}\n"
                        f"  • {sdf_count} SDF files\n"
                        f"  • {zero_count} zero‐score entries\n"
                        f"  • {expected} expected entries\n"
                        f"  • Success rate = {success_rate:.1%}"
                    )

            run_results[run_dir] = run_stats

            if task_success_rates:
                mean_run_success = sum(task_success_rates) / len(task_success_rates)
                per_run_success_rates.append(mean_run_success)
                if print_per_task and per_task_output:
                    print(f"\n📦 Summary for run: {run_dir}")
                    print('\n'.join(per_task_output))
                if any_flag:
                    print(f"⚠️ Some inconsistencies were detected in run {run_dir} (see above warnings).")
            else:
                per_run_success_rates.append(float('nan'))
                print(f"⚠️ No valid tasks found in run {run_dir}.")

    # -- Print per-run stats --
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

    # -- Print overall --
    print("\n=== ✅ Overall Summary ===")
    if overall["expected"] > 0:
        overall_rate = overall["sdf"] / overall["expected"]
        print(f"Total expected: {overall['expected']}")
        print(f"Total SDF files: {overall['sdf']}")
        print(f"Total zero scores: {overall['zero']}")
        print(f"Overall success rate = {overall_rate:.1%}")
    else:
        print("⚠️ No valid tasks found. Check folder structure or filenames.")

    return per_run_success_rates, (mean_rate if rates_clean else None), (min_rate if rates_clean else None), (max_rate if rates_clean else None)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Analyze 3DOpt benchmark results. Computes per-run and overall success rates. "
                    "Optionally print results per task."
    )
    parser.add_argument("folder", help="Folder containing runs to analyze (e.g. Results_RndSam_RDKit_10Reps)")
    parser.add_argument("--per-task", action="store_true", default=False,
                        help="Print per-task results for each run")

    args = parser.parse_args()
    analyze_results(args.folder, print_per_task=args.per_task)
