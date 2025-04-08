import os
import sys
import time
import numpy as np
import pandas as pd
import multiprocessing as mp
from collections import defaultdict, Counter
from ccdc.molecule import Molecule
from ccdc.io import MoleculeWriter
from ccdc import io
from ccdc.entry import Entry
from hsr import fingerprint as fp
from hsr import similarity as sim

os.environ["QT_QPA_PLATFORM"] = "xcb"  # For headless systems

### ---------- Utility Functions ---------- ###

def get_array_from_ccdcmol(ccdcmol):
    coords = np.array([
        [atom.coordinates[0], atom.coordinates[1], atom.coordinates[2], np.sqrt(atom.atomic_number)]
        for atom in ccdcmol.atoms
    ])
    return coords - coords.mean(axis=0)

def fingerprint_key(fp_obj):
    return tuple(fp_obj.tolist() if isinstance(fp_obj, np.ndarray) else fp_obj)

def formula_signature(fragment):
    atoms = fragment.atoms
    central = atoms[0].atomic_symbol
    n_atoms = len(atoms)
    counts = Counter(atom.atomic_symbol for atom in atoms)
    formula = ''.join(f"{el}{counts[el]}" for el in sorted(counts))
    return (central, n_atoms, formula)


def generate_fp_data(fragment):
    return {
        "sdf": fragment.to_string("sdf"),
        "fp": fp.generate_fingerprint_from_data(get_array_from_ccdcmol(fragment)),
        "formula": formula_signature(fragment),
        "n_atoms": len(fragment.atoms),
        "central_atom": fragment.atoms[0].atomic_symbol,
    }

def interatomic_distance(sdf_string):
    mol = Molecule.from_string(sdf_string, format="sdf")
    coords = [atom.coordinates for atom in mol.atoms]
    return np.linalg.norm(np.array(coords[0]) - np.array(coords[1]))

### ---------- Fragmentation Functions ---------- ###

def component_of_interest(molecule):
    components = molecule.components
    if not components:
        return None
    props = [{
        "component": c,
        "is_organometallic": c.is_organometallic,
        "mw": sum(atom.atomic_weight for atom in c.atoms),
        "atom_count": len(c.atoms)
    } for c in components]
    heaviest = max(props, key=lambda x: x["mw"])
    most_atoms = max(props, key=lambda x: x["atom_count"])
    for prop in props:
        if sum([
            prop["is_organometallic"],
            prop["component"] == heaviest["component"],
            prop["component"] == most_atoms["component"]
        ]) >= 2 and prop["atom_count"] >= 5:
            return prop["component"]
    return None

def create_fragment(central_atom):
    frag = Molecule(identifier=f"{central_atom.label}_frag")
    atom_map = {central_atom: frag.add_atom(central_atom)}
    for neighbor in central_atom.neighbours:
        atom_map[neighbor] = frag.add_atom(neighbor)
    for bond in central_atom.bonds:
        a1, a2 = bond.atoms
        if a1 in atom_map and a2 in atom_map:
            try:
                frag.add_bond(bond.bond_type, atom_map[a1], atom_map[a2])
            except:
                pass
    return frag

def get_fragments(mol):
    return [create_fragment(atom) for atom in mol.atoms]

### ---------- Target Processing ---------- ###

def process_target(entry_id, threshold=0.999):
    reader = io.EntryReader("CSD")
    mol = component_of_interest(reader.entry(entry_id).molecule)
    fragments = get_fragments(mol)

    grouped = defaultdict(list)
    for frag in fragments:
        data = generate_fp_data(frag)
        grouped[data["formula"]].append(data)

    unique_by_formula = defaultdict(list)
    for formula, frags in grouped.items():
        for frag in frags:
            if all(sim.compute_similarity_score(frag["fp"], other["fp"]) < threshold for other in unique_by_formula[formula]):
                unique_by_formula[formula].append(frag)

    return unique_by_formula

### ---------- Chunked Loader ---------- ###

import csv
import ast
import pandas as pd
from collections import defaultdict

def load_population_subset_from_index(pop_ids, formulas, csv_path, index_path="fragment_index.csv", chunk_size=100_000):
    """
    Returns a dict:
       {
          (central_atom, n_atoms, formula_str): [
              {
                 'fp': [float, float, ...],
                 'sdf': str,
                 'n_atoms': int,
                 'formula': (central_atom, n_atoms, formula_str)
              },
              ...
          ],
          ...
       }

    pop_ids:   an iterable of entry IDs we want (case-insensitive).
    formulas:  an iterable of (central_atom, n_atoms, formula_str) from the target.
    csv_path:  path to the big CSV file containing the actual fragment data (with columns fp, sdf, formula, etc.)
    index_path: path to 'fragment_index.csv' which has columns:
                chunk_id, row_in_chunk, entry_id, formula
    chunk_size: how many rows at a time to read from csv_path.
    """

    print(f"🔎 Filtering index for {len(pop_ids)} entry_ids and {len(formulas)} formulas...")

    # Normalize entry IDs and formulas for consistent comparison:
    pop_ids = set(pid.strip().upper() for pid in pop_ids)
    normalized_formulas = set((str(f[0]), int(f[1]), str(f[2])) for f in formulas)

    # For debugging, show a sample formula and some sample IDs
    print("🔬 Sample normalized formula:", next(iter(normalized_formulas)))
    print("📋 Sample normalized pop_ids:", list(pop_ids)[:5])

    # We'll track rows to read from the big CSV in rows_by_chunk:
    #   rows_by_chunk[chunk_id] = [row_index1, row_index2, ...]
    rows_by_chunk = defaultdict(list)
    matched_ids = set()

    # 1) Safely read fragment_index.csv using csv.DictReader
    with open(index_path, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                chunk_id = int(row['chunk_id'])
                row_in_chunk = int(row['row_in_chunk'])
                entry_id = row['entry_id'].strip().upper()
                formula_str = row['formula']  # e.g. "('C', 5, 'C4H1')"

                # Convert the string in "formula" column into a tuple
                formula_raw = ast.literal_eval(formula_str)
                # Force it into the shape (str, int, str)
                formula_parsed = (str(formula_raw[0]), int(formula_raw[1]), str(formula_raw[2]))

                # If this entry is one of our population IDs and the formula is one of our target formulas:
                if entry_id in pop_ids and formula_parsed in normalized_formulas:
                    rows_by_chunk[chunk_id].append(row_in_chunk)
                    matched_ids.add(entry_id)

            except Exception:
                # If something is malformed, skip it
                continue

    # If we found no matching lines in the index, just return empty.
    if not rows_by_chunk:
        print("⚠️ No matching fragments found in index!")
        return {}

    # Print some stats about how many rows and how many unique entry_ids we matched
    print(f"✅ Found {sum(len(v) for v in rows_by_chunk.values())} rows across {len(rows_by_chunk)} chunks.")
    print(f"🧩 Total unique matched entry_ids: {len(matched_ids)}")

    # 2) Now we actually load the big CSV (csv_path), chunk by chunk,
    #    but only the rows we need from each chunk.
    pop_group = defaultdict(list)
    reader = pd.read_csv(csv_path, chunksize=chunk_size)
    for chunk_id, chunk in enumerate(reader):
        if chunk_id not in rows_by_chunk:
            continue

        # row indexes from index file
        row_idxs = rows_by_chunk[chunk_id]

        try:
            subset = chunk.loc[row_idxs]
        except KeyError:
            print(f"⚠️ Skipping chunk {chunk_id} due to missing row indexes.")
            continue

        # Build the pop_group dictionary from relevant rows
        for _, row in subset.iterrows():
            try:
                # row["formula"] is like "('C', 5, 'C4H1')"
                parsed = ast.literal_eval(row["formula"])
                formula_tuple = (str(parsed[0]), int(parsed[1]), str(parsed[2]))

                # Convert the fingerprint string (like "[0.1, 0.2, 0.05]") to float list
                fp_floats = [float(x) for x in row["fp"].strip("[] ").split(",")]

                frag = {
                    "fp": fp_floats,
                    "sdf": row["sdf"],
                    "n_atoms": row["n_atoms"],
                    "formula": formula_tuple
                }
                pop_group[formula_tuple].append(frag)

            except Exception:
                # skip malformed lines in the big CSV
                continue

    return pop_group



# def load_population_subset_from_index(pop_ids, formulas, csv_path, index_path="fragment_index.csv", chunk_size=100_000):
#     import ast

#     print(f"🔎 Filtering index for {len(pop_ids)} entry_ids and {len(formulas)} formulas...")

#     # Convert to sets for fast lookup
#     pop_ids = set(pop_ids)
#     formulas = set(formulas)

#     pop_group = defaultdict(list)
#     rows_by_chunk = defaultdict(list)

#     # Process index file line-by-line to save memory
#     with open(index_path) as f:
#         next(f)  # skip header
#         for line in f:
#             try:
#                 chunk_id, row_idx, entry_id, formula_str = line.strip().split(",", 3)
#                 entry_id = entry_id.strip()
#                 formula = ast.literal_eval(formula_str.strip())
#                 if entry_id in pop_ids and formula in formulas:
#                     rows_by_chunk[int(chunk_id)].append(int(row_idx))
#             except Exception:
#                 continue

#     if not rows_by_chunk:
#         print("⚠️ No matching fragments found in index!")
#         return {}

#     print(f"✅ Found {sum(len(v) for v in rows_by_chunk.values())} rows across {len(rows_by_chunk)} chunks.")

#     # Load only relevant rows from CSV
#     reader = pd.read_csv(csv_path, chunksize=chunk_size)

#     for chunk_id, chunk in enumerate(reader):
#         if chunk_id not in rows_by_chunk:
#             continue

#         row_idxs = rows_by_chunk[chunk_id]

#         try:
#             subset = chunk.loc[row_idxs]
#         except KeyError:
#             print(f"⚠️ Skipping chunk {chunk_id} due to missing row indexes.")
#             continue

#         for _, row in subset.iterrows():
#             try:
#                 formula = ast.literal_eval(row["formula"])
#                 frag = {
#                     "fp": [float(x) for x in row["fp"].strip("[] ").split(",")],
#                     "sdf": row["sdf"],
#                     "n_atoms": row["n_atoms"],
#                     "formula": formula
#                 }
#                 pop_group[formula].append(frag)
#             except Exception:
#                 continue

#     return pop_group

### ---------- Similarity Comparison ---------- ###

def compare_group(args):
    target_group, pop_group, threshold = args
    target_status = {
        fingerprint_key(t["fp"]): {
            "data": t,
            "top_matches": [],
            "matched": False
        } for t in target_group
    }
    unmatched_keys = set(target_status.keys())

    for pop in pop_group:
        pop_fp = pop["fp"]
        pop_sdf = pop["sdf"]
        pop_formula = pop["formula"]
        pop_atoms = pop["n_atoms"]

        to_remove = []

        for key in unmatched_keys:
            t = target_status[key]["data"]
            is_both_biatomic = t["n_atoms"] == 2 and pop_atoms == 2

            if is_both_biatomic and t["formula"] == pop_formula:
                try:
                    t_dist = interatomic_distance(t["sdf"])
                    p_dist = interatomic_distance(pop_sdf)
                    diff = abs(t_dist - p_dist)
                    if diff <= 0.01:
                        score = 1.0 - diff
                        target_status[key]["top_matches"].append((score, pop_sdf))
                        target_status[key]["matched"] = True
                        to_remove.append(key)
                except:
                    continue
            elif not is_both_biatomic:
                score = sim.compute_similarity_score(t["fp"], pop_fp)
                if score >= threshold:
                    target_status[key]["top_matches"].append((score, pop_sdf))
                    target_status[key]["matched"] = True
                    to_remove.append(key)

            target_status[key]["top_matches"].sort(key=lambda x: -x[0])
            target_status[key]["top_matches"] = target_status[key]["top_matches"][:3]

        unmatched_keys -= set(to_remove)
        if not unmatched_keys:
            break

    return [
        (key, {
            "target_sdf": info["data"]["sdf"],
            "top_matches": info["top_matches"],
            "matched": info["matched"]
        }) for key, info in target_status.items()
    ]

def compare_fragments_parallel(target_frags, pop_fragments_gen, threshold=0.999, n_processes=8):
    from time import perf_counter

    all_results = defaultdict(lambda: {
        "target_sdf": None,
        "top_matches": [],
        "matched": False
    })

    print(f"🧩 Preparing group comparison tasks...")
    prep_start = perf_counter()

    tasks = []
    chunk_idx = 0
    total_pop_groups = 0

    for pop_group in pop_fragments_gen:
        chunk_idx += 1
        chunk_start = perf_counter()

        matches_this_chunk = 0
        for key, pop_list in pop_group.items():
            total_pop_groups += 1
            if key in target_frags:
                tasks.append((target_frags[key], pop_list, threshold))
                matches_this_chunk += 1

        chunk_end = perf_counter()
        print(f"  🔹 Processed chunk {chunk_idx} in {chunk_end - chunk_start:.2f}s "
              f"(matched {matches_this_chunk} target formulas)", flush=True)

    prep_end = perf_counter()
    print(f"🛠️ Prepared {len(tasks)} group tasks from {chunk_idx} chunks in {prep_end - prep_start:.2f}s")

    print(f"🚀 Starting multiprocessing with {n_processes} processes...")
    mp_start = perf_counter()

    with mp.Pool(n_processes) as pool:
        for batch in pool.imap_unordered(compare_group, tasks):
            for key, result in batch:
                all_results[key]["target_sdf"] = result["target_sdf"]
                all_results[key]["top_matches"].extend(result["top_matches"])
                all_results[key]["matched"] |= result["matched"]
                all_results[key]["top_matches"].sort(key=lambda x: -x[0])
                all_results[key]["top_matches"] = all_results[key]["top_matches"][:3]

    mp_end = perf_counter()
    print(f"✅ Finished multiprocessing in {mp_end - mp_start:.2f}s")

    return dict(all_results)

### ---------- Main ---------- ###

def run_analysis(entry_id, population_file, idx, csv_path="all_fragments_data_all.csv", threshold=0.999):
    start = time.time()
    print(f"\n🔍 Target: {entry_id}")

    output_dir = "frag_comparison_results_chunked"
    os.makedirs(output_dir, exist_ok=True)

    t0 = time.time()
    target_frags = process_target(entry_id, threshold)
    t1 = time.time()
    total = sum(len(v) for v in target_frags.values())
    print(f"✅ Unique target fragments: {total} (retrieved in {t1 - t0:.2f}s)")
    print(f'✅Unique formulas: {target_frags.keys()}')

    with open(population_file) as f:
        pop_ids = [line.split()[0] for line in f if line.strip()]
    print(f"⚙️ Population size: {len(pop_ids)}")
    
    t2 = time.time()
    pop_group = load_population_subset_from_index(
                pop_ids=pop_ids,
                formulas=target_frags.keys(),
                csv_path="all_fragments_data_all.csv",
                index_path="fragment_index.csv")
    pop_fragments_gen = [pop_group]  # single-chunk generator format

    print(f"📦 Started loading population fragments...")

    t3 = time.time()
    comparisons = compare_fragments_parallel(target_frags, pop_fragments_gen, threshold=0.99, n_processes=8)
    t4 = time.time()
    print(f"🔗 Similarity comparison done in {t4 - t3:.2f}s")

    matched = sum(1 for v in comparisons.values() if v["matched"])
    print(f"📊 Matched {matched}/{total} fragments.")

    t5 = time.time()
    with MoleculeWriter(os.path.join(output_dir, f"{idx}_{entry_id}_target_unique_fragments.sdf")) as w:
        for group in target_frags.values():
            for frag in group:
                w.write(Molecule.from_string(frag["sdf"], format="sdf"))
    print(f"🧪 Wrote target fragments SDF in {time.time() - t5:.2f}s")

    t6 = time.time()
    for i, (fp_key, comp) in enumerate(comparisons.items(), start=1):
        output_path = os.path.join(output_dir, f"{idx}_{entry_id}_frag{i}_matches.sdf")
        with MoleculeWriter(output_path) as writer:
            target_mol = Molecule.from_string(comp["target_sdf"], format="sdf")
            writer.write(target_mol)

            for sim_score, sdf in comp["top_matches"]:
                match_mol = Molecule.from_string(sdf, format="sdf")
                match_entry = Entry.from_molecule(match_mol)
                if len(match_mol.atoms) == 2 and len(target_mol.atoms) == 2:
                    try:
                        d1 = interatomic_distance(target_mol.to_string("sdf"))
                        d2 = interatomic_distance(match_mol.to_string("sdf"))
                        match_entry.attributes["DistanceDifference"] = f"{abs(d1 - d2):.4f}"
                    except:
                        match_entry.attributes["DistanceDifference"] = "ERROR"
                else:
                    match_entry.attributes["Similarity"] = f"{sim_score:.4f}"
                writer.write_entry(match_entry)
    print(f"📁 Wrote match SDFs in {time.time() - t6:.2f}s")

    print(f"✅ Done with target {entry_id} in {time.time() - start:.2f}s")


### ---------- Entry ---------- ###

if __name__ == "__main__":
    overall_start = time.time()
    targets = ['ABAHIW', 'ABAKIZ', 'ABADOX', 'ABABIP', 'GASQOK', 'ABEKIE',
               'NIWPUE01', 'ABEKIF', 'APUFEX', 'ABEHAU', 'TITTUO', 'EGEYOG',
               'ABOBUP', 'XIDTOW', 'ACNCOB10', 'TACXUQ', 'ACAZFE', 'NIVHEJ',
               'ADUPAS', 'DAJLAC', 'OFOWIS', 'CATSUL', 'HESMUQ01', 'GUDQOL',
               'ABEVAG', 'AKOQOH', 'ADARUT', 'AFECIA', 'ACOVUL', 'AFIXEV']

    for i, target in enumerate(targets, start=1):
        run_analysis(target, f"targets/init_populations_protons_2/{i}_{target}_init_pop.txt", i)

    print(f"\n⏱️ Total time: {time.time() - overall_start:.2f} seconds")
