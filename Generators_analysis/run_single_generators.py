#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Single-run benchmarking for CSD to 3D structure generation, for aggregation.
Each run writes results to its own folder, for post-processing by multi-run scripts.
"""

import os, sys
import time
import numpy as np
import pandas as pd
import random
import multiprocessing
from ccdc import io
from ccdc.conformer import ConformerGenerator
from ccdc.molecule import Molecule
from openbabel import pybel as pb
from rdkit import Chem
from rdkit.Chem import AllChem

from hsr import pre_processing as pp
from hsr import fingerprint as fp
from hsr import similarity as sim
from hsr.utils import PROTON_FEATURES

# =========================== Parameters ============================
TARGETS = [
    'ABAHIW', 'ABAKIZ', 'ABADOX', 'ABABIP', 'GASQOK', 'ABEKIE', 'NIWPUE01',
    'ABEKIF', 'APUFEX', 'ABEHAU', 'TITTUO', 'EGEYOG', 'ABOBUP', 'XIDTOW',
    'ACNCOB10', 'TACXUQ', 'ACAZFE', 'NIVHEJ', 'ADUPAS', 'DAJLAC', 'OFOWIS',
    'CATSUL', 'HESMUQ01', 'GUDQOL', 'ABEVAG', 'AKOQOH', 'ADARUT', 'AFECIA',
    'ACOVUL', 'AFIXEV', 'ABAYAF', 'RULJAM'
]

if len(sys.argv) >= 3:
    BASE_OUT = sys.argv[1]
    run_seed = int(sys.argv[2])
    RUN_ID = os.path.basename(BASE_OUT)
else:
    RUN_ID = os.environ.get('RUN_ID', 'run1')
    BASE_OUT = os.path.abspath(f"3D_Generators_single_runs/{RUN_ID}")
    run_seed = np.random.randint(0, 2**32-1)

os.makedirs(BASE_OUT, exist_ok=True)
CCDC_TIMEOUT = 10

NA = "N.A."

def fmt(v):
    return NA if v is None or (isinstance(v, float) and (np.isnan(v) or np.isinf(v))) else f"{v:.2f}"

def component_of_interest(molecule):
    comps = molecule.components
    if not comps: return None
    props = [{
        "c": c,
        "is_organomet": c.is_organometallic,
        "mw": sum(a.atomic_weight for a in c.atoms),
        "nat": len(c.atoms)
    } for c in comps]
    heaviest = max(props, key=lambda x: x["mw"])
    most_atoms = max(props, key=lambda x: x["nat"])
    for p in props:
        score = sum([p["is_organomet"], p["c"] is heaviest["c"], p["c"] is most_atoms["c"]])
        if score >= 2 and p["nat"] >= 5: return p["c"]
    return None

def get_array_from_ccdcmol(ccdcmol):
    arr = np.array([[a.coordinates[0], a.coordinates[1], a.coordinates[2], np.sqrt(a.atomic_number)] for a in ccdcmol.atoms])
    arr -= arr.mean(axis=0)
    return arr

def get_array_from_pybelmol(pybelmol):
    arr = np.array([[a.coords[0], a.coords[1], a.coords[2], np.sqrt(a.atomicnum)] for a in pybelmol.atoms])
    arr -= arr.mean(axis=0)
    return arr

# ============= Timeout-protected CCDC, with process isolation =============
def _ccdc_worker(smiles, q, run_seed):
    try:
        random.seed(run_seed)
        cg = ConformerGenerator()
        nconf = random.randint(2, 100)
        cg.settings.max_conformers = nconf
        mol0 = Molecule.from_string(smiles)
        hits = cg.generate(mol0).hits
        if not hits:
            q.put(None)
        else:
            q.put(hits[0].molecule.to_string("sdf"))
    except Exception:
        q.put(None)

def generate_ccdc_sdf(smiles, timeout=CCDC_TIMEOUT, run_seed=None):
    ctx = multiprocessing.get_context("fork")
    q = ctx.Queue()
    p = ctx.Process(target=_ccdc_worker, args=(smiles, q, run_seed))
    p.start()
    p.join(timeout)
    if p.is_alive():
        p.kill()
        p.join()
    return q.get() if not q.empty() else None

def gen_with_ccdc(smiles, run_seed):
    sdf = generate_ccdc_sdf(smiles, timeout=CCDC_TIMEOUT, run_seed=run_seed)
    if sdf is None:
        raise RuntimeError("CCDC generation failed or timed-out")
    mol = Molecule.from_string(sdf)
    return mol, sdf

def gen_with_rdkit(smiles, run_seed):
    # --- PRINT SMILES for debug ---
    print(f"  [RDKit] SMILES: {smiles}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"  [RDKit] MolFromSmiles failed for SMILES: {smiles}")
        raise ValueError("SMILES cannot be parsed by RDKit")
    mol = Chem.AddHs(mol)
    params = AllChem.ETKDG()
    params.randomSeed = run_seed
    AllChem.EmbedMolecule(mol, params)
    AllChem.MMFFOptimizeMolecule(mol)
    return mol, Chem.MolToMolBlock(mol)

def gen_with_obabel(smiles, run_seed):
    os.environ["OB_RANDOM_SEED"] = str(run_seed)
    mol = pb.readstring("smi", smiles)
    mol.addh()
    mol.make3D()
    mol.localopt()
    return mol, mol.write("sdf")

# =========================== Main Run Loop ================================
def main():
    print(f"==== RUN SEED: {run_seed} ====")
    methods = ["ccdc", "rdkit", "obabel"]
    sdf_dirs = {m: os.path.join(BASE_OUT, m + "_sdfs") for m in methods}
    for d in sdf_dirs.values():
        os.makedirs(d, exist_ok=True)
    run_rows = []

    reader = io.EntryReader("CSD")
    for i, refcode in enumerate(TARGETS, 1):
        print(f"\nTarget {i}/{len(TARGETS)}  {refcode}")
        try:
            mol_cryst = component_of_interest(reader.entry(refcode).molecule)
            smiles    = mol_cryst.smiles
            print(f"  [DEBUG] Using SMILES: {smiles}")  # Print for all targets
            if smiles is None:
                print(f"Could not obtain SMILES from target {refcode}")
            arr_orig  = get_array_from_ccdcmol(mol_cryst)
            fp_orig   = fp.generate_fingerprint_from_data(arr_orig)
        except Exception as e:
            print(f"  ❌  Could not load crystal structure: {e}")
            row = {"Target": refcode}
            for m in methods:
                row[f"{m}_sim"] = NA
                row[f"{m}_time"] = NA
                row[f"{m}_ok"] = "–"
            run_rows.append(row)
            continue

        # Try each generator
        res = {}
        for m, func in zip(methods, [gen_with_ccdc, gen_with_rdkit, gen_with_obabel]):
            t0 = time.time()
            try:
                mol3d, sdf = func(smiles, run_seed)
                elapsed = time.time() - t0
                if m == "obabel":
                    arr = get_array_from_pybelmol(mol3d)
                elif m == "rdkit":
                    arr = pp.molecule_to_ndarray(mol3d, features=PROTON_FEATURES, removeHs=False)
                else:  # ccdc
                    arr = get_array_from_ccdcmol(mol3d)
                fp_gen  = fp.generate_fingerprint_from_data(arr)
                sim_val = sim.compute_similarity_score(fp_orig, fp_gen)
                res[m] = {"time": elapsed, "sim": sim_val, "sdf": sdf}
                # Save SDF in dedicated subfolder for this method
                sdf_filename = os.path.join(sdf_dirs[m], f"{i}_{refcode}_{m}.sdf")
                with open(sdf_filename, "w") as fh:
                    fh.write(sdf)
            except Exception as e:
                print(f"  ⚠️  {m.upper()} failed: {e}")
                # Add exception message to result for further analysis
                res[m] = None

        # Report
        row = {"Target": refcode}
        for m in methods:
            row[f"{m}_sim"] = fmt(res[m]["sim"]) if res[m] else NA
            row[f"{m}_time"] = fmt(res[m]["time"]) if res[m] else NA
            row[f"{m}_ok"] = "✓" if res[m] else "–"
            print(f"  {m.upper():6} time: {row[f'{m}_time']} s   similarity: {row[f'{m}_sim']}")
        run_rows.append(row)

    # Save CSV with all per-target scores and times for this run
    df = pd.DataFrame(run_rows)
    run_csv = os.path.join(BASE_OUT, f"scores_and_times_{RUN_ID}.csv")
    df.to_csv(run_csv, index=False)
    print(f"\n[Run {RUN_ID}] Results saved to {run_csv}")

    # Save as tab-separated TXT
    summary_table = df[["Target"] + [f"{m}_ok" for m in methods]]
    summary_txt = os.path.join(BASE_OUT, "success_summary.txt")
    summary_table.to_csv(summary_txt, sep="\t", index=False)
    print(f"[Run {RUN_ID}] Success table saved to {summary_txt}")

    # Save as CSV
    summary_csv = os.path.join(BASE_OUT, "success_summary.csv")
    summary_table.to_csv(summary_csv, sep=",", index=False)
    print(f"[Run {RUN_ID}] CSV tick file saved to {summary_csv}")

    print("\n✅ Single run finished.")

if __name__ == "__main__":
    main()
