import os
import glob
import re
from ccdc.io import MoleculeReader

def extract_cae_index(filename):
    match = re.search(r'_cae(\d+)_matches\.sdf$', filename)
    return int(match.group(1)) if match else -1

def extract_field(sdf_text, field):
    pattern = rf"<{field}>\s*(.*?)\s*(?:\n|\r|\$\$\$\$)"
    match = re.search(pattern, sdf_text)
    return match.group(1).strip() if match else None

# --- Target list ---
targets = ['ABAHIW', 'ABAKIZ', 'ABADOX', 'ABABIP', 'GASQOK', 'ABEKIE',
           'NIWPUE01', 'ABEKIF', 'APUFEX', 'ABEHAU', 'TITTUO', 'EGEYOG',
           'ABOBUP', 'XIDTOW', 'ACNCOB10', 'TACXUQ', 'ACAZFE', 'NIVHEJ',
           'ADUPAS', 'DAJLAC', 'OFOWIS', 'CATSUL', 'HESMUQ01', 'GUDQOL',
           'ABEVAG', 'AKOQOH', 'ADARUT', 'AFECIA', 'ACOVUL', 'AFIXEV', 'ABAYAF', 'RULJAM']

folder_path = "../CAEs/cae_comparison_results_0_5" # <--- CHANGE THIS IF NEEDED

total_all = 0
matched_all = 0

print()

for idx, target in enumerate(targets, start=1):
    base_name = f"{idx}_{target}"

    print(f"🔍 Processing: {base_name}")

    cae_pattern = os.path.join(folder_path, f"{base_name}_cae*_matches.sdf")
    cae_files = sorted(glob.glob(cae_pattern), key=extract_cae_index)

    total_files = len(cae_files)
    matched_caes = 0
    skipped_caes = []
    distorted_caes = []

    for cae_file in cae_files:
        with open(cae_file, 'r') as f:
            sdf_text = f.read()

        mols = list(MoleculeReader(cae_file))
        if len(mols) < 2:
            skipped_caes.append(os.path.basename(cae_file))
            continue

        # If "DistortedMatch" is present and set to True, count as distorted (not a match)
        is_distorted = False
        # The DistortedMatch field is present for matches below the threshold (see your script)
        distorted_value = extract_field(sdf_text, "DistortedMatch")
        if distorted_value and distorted_value.strip().lower() == "true":
            sim = extract_field(sdf_text, "Similarity")
            dist = extract_field(sdf_text, "DistanceDifference")
            distorted_caes.append((os.path.basename(cae_file), sim or dist))
            is_distorted = True

        if not is_distorted:
            matched_caes += 1

    num_skipped = len(skipped_caes)
    num_distorted = len(distorted_caes)
    num_matches = matched_caes

    # CAE-level reporting
    if skipped_caes:
        print(f"  🚫 Missing formula (1 mol only):")
        for f in skipped_caes:
            print(f"     - {f}")

    if distorted_caes:
        print(f"  ⚠️ Distorted match (below threshold):")
        for fname, val in distorted_caes:
            print(f"     - {fname} (value: {val})")

    print(f"  → Number of CAEs: {total_files}")
    print(f"  → Number of matches: {num_matches}")
    print(f"  → Number of distorted matches: {num_distorted}")
    print(f"  → Number of missing formulas: {num_skipped}")
    match_rate = 100 * num_matches / total_files if total_files else 0
    print(f"  → Match rate: {match_rate:.2f}%\n")

    # Global counter
    total_all += total_files
    matched_all += num_matches

# --- Final summary ---
print("======= FINAL SUMMARY =======")
print(f"Total CAEs: {total_all}")
print(f"Matched CAEs: {matched_all}")
overall_match_rate = (100 * matched_all / total_all) if total_all else 0
print(f"Overall match rate: {overall_match_rate:.2f}%")
