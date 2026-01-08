import pandas as pd

# Load CSV file
df = pd.read_csv("./3D_Generators_analysis/average_scores_and_times.csv")

def compute_means(sim_col, time_col):
    # Similarity means
    mean_sim_success = df[sim_col].mean(skipna=True)         # only successes
    mean_sim_all = df[sim_col].fillna(0).mean()              # all 32 targets (0 for failures)
    
    # Time mean (only successes)
    mean_time_success = df[time_col].mean(skipna=True)
    
    return mean_sim_success, mean_sim_all, mean_time_success

# Compute for each generator
ccdc_succ_sim, ccdc_all_sim, ccdc_time = compute_means('ccdc_sim_mean', 'ccdc_time_mean')
obabel_succ_sim, obabel_all_sim, obabel_time = compute_means('obabel_sim_mean', 'obabel_time_mean')
rdkit_succ_sim, rdkit_all_sim, rdkit_time = compute_means('rdkit_sim_mean', 'rdkit_time_mean')

# Print results
print(f"CCDC:   Mean similarity (success only) = {ccdc_succ_sim:.3f}, Mean similarity (all targets) = {ccdc_all_sim:.3f}, Mean time = {ccdc_time:.3f} s")
print(f"OBabel: Mean similarity (success only) = {obabel_succ_sim:.3f}, Mean similarity (all targets) = {obabel_all_sim:.3f}, Mean time = {obabel_time:.3f} s")
print(f"RDKit:  Mean similarity (success only) = {rdkit_succ_sim:.3f}, Mean similarity (all targets) = {rdkit_all_sim:.3f}, Mean time = {rdkit_time:.3f} s")
import pandas as pd

# Load CSV file
df = pd.read_csv("./3D_Generators_analysis/average_scores_and_times.csv")

def compute_means(sim_col, time_col):
    # Similarity means
    mean_sim_success = df[sim_col].mean(skipna=True)         # only successes
    mean_sim_all = df[sim_col].fillna(0).mean()              # all 32 targets (0 for failures)
    
    # Time mean (only successes)
    mean_time_success = df[time_col].mean(skipna=True)
    
    # Success rate
    success_rate = df[sim_col].notna().mean()  # fraction of targets with non-NaN similarity
    
    return mean_sim_success, mean_sim_all, mean_time_success, success_rate

# Compute for each generator
ccdc_succ_sim, ccdc_all_sim, ccdc_time, ccdc_sr = compute_means('ccdc_sim_mean', 'ccdc_time_mean')
obabel_succ_sim, obabel_all_sim, obabel_time, obabel_sr = compute_means('obabel_sim_mean', 'obabel_time_mean')
rdkit_succ_sim, rdkit_all_sim, rdkit_time, rdkit_sr = compute_means('rdkit_sim_mean', 'rdkit_time_mean')

# Find common set of successes (intersection)
common_success = df.dropna(subset=['ccdc_sim_mean', 'obabel_sim_mean', 'rdkit_sim_mean'])
common_size = len(common_success)
common_rate = common_size / len(df)

# Mean similarities on the common set
ccdc_common_sim = common_success['ccdc_sim_mean'].mean()
obabel_common_sim = common_success['obabel_sim_mean'].mean()
rdkit_common_sim = common_success['rdkit_sim_mean'].mean()

# Print results
print(f"CCDC:   Mean similarity (success only) = {ccdc_succ_sim:.3f}, "
      f"Mean similarity (all targets) = {ccdc_all_sim:.3f}, "
      f"Mean time = {ccdc_time:.3f} s, "
      f"Success rate = {ccdc_sr*100:.1f}%")

print(f"OBabel: Mean similarity (success only) = {obabel_succ_sim:.3f}, "
      f"Mean similarity (all targets) = {obabel_all_sim:.3f}, "
      f"Mean time = {obabel_time:.3f} s, "
      f"Success rate = {obabel_sr*100:.1f}%")

print(f"RDKit:  Mean similarity (success only) = {rdkit_succ_sim:.3f}, "
      f"Mean similarity (all targets) = {rdkit_all_sim:.3f}, "
      f"Mean time = {rdkit_time:.3f} s, "
      f"Success rate = {rdkit_sr*100:.1f}%")

print(f"\nCommon set: {common_size} out of {len(df)} targets "
      f"({common_rate*100:.1f}% of the benchmark)")

print(f"  CCDC mean similarity (common set)   = {ccdc_common_sim:.3f}")
print(f"  OBabel mean similarity (common set) = {obabel_common_sim:.3f}")
print(f"  RDKit mean similarity (common set)  = {rdkit_common_sim:.3f}")
