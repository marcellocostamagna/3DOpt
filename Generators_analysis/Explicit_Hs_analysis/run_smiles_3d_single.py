
"""
Run a single benchmark of SMILES→3D generation (implicit vs explicit-H SMILES).

- Inputs: CSD_SMILES (implicit H) and RDKit_SMILES_EXPL (explicit H).
- Methods: CCDC, RDKit, Open Babel (with optional H addition flags).
- For each (target × source × method): generate 3D structure, compute HSR
  similarity vs CSD reference.
- Outputs (run root):
    * similarity_comparison.csv   (pivot table; "N.A." = failure)
    * atom_counts_comparison.csv  (reference, SMILES, generated atom counts)

Usage:
    python run_smiles_3d_single.py --outdir OUT --seed 12345
"""

import os
import time
import argparse
import random
import multiprocessing
from typing import Optional, Dict, List, Tuple

import numpy as np
import pandas as pd

# CCDC
from ccdc import io
from ccdc.conformer import ConformerGenerator
from ccdc.molecule import Molecule

# RDKit
from rdkit import Chem
from rdkit.Chem import AllChem

# Open Babel (pybel)
from openbabel import pybel as pb

# HSR stack
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

CCDC_TIMEOUT = 10
NA = "N.A."


# =========================== Utilities ============================
def fmt_num(v: Optional[float]) -> str:
    if v is None:
        return NA
    if isinstance(v, float) and (np.isnan(v) or np.isinf(v)):
        return NA
    return f"{v:.2f}"


def fmt_count(v) -> str:
    return NA if v is None or (isinstance(v, float) and (np.isnan(v) or np.isinf(v))) else str(int(v))


def ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def component_of_interest(molecule: Molecule) -> Optional[Molecule]:
    """Heuristic: organometallic OR heaviest OR most atoms, with >=5 atoms."""
    comps = molecule.components
    if not comps:
        return None

    props = []
    for c in comps:
        try:
            mw = sum(a.atomic_weight for a in c.atoms)
            nat = len(c.atoms)
            is_org = c.is_organometallic
            props.append({"c": c, "mw": mw, "nat": nat, "is_org": is_org})
        except Exception:
            continue

    if not props:
        return None

    heaviest = max(props, key=lambda x: x["mw"])
    most_atoms = max(props, key=lambda x: x["nat"])

    for p in props:
        score = int(p["is_org"]) + int(p["c"] is heaviest["c"]) + int(p["c"] is most_atoms["c"])
        if score >= 2 and p["nat"] >= 5:
            return p["c"]
    return None


# -------------------- SMILES extractors --------------------
def smiles_from_csd(comp: Molecule) -> Optional[str]:
    """Canonical SMILES directly from CCDC component."""
    try:
        s = (comp.smiles or "").strip()
        return s or None
    except Exception:
        return None


def smiles_from_rdkit_explicit_from_sdf(sdf_block: str) -> Optional[str]:
    """
    Build SMILES from the SDF using RDKit without sanitization and without removing H,
    then write with all hydrogens explicit.
    """
    try:
        m = Chem.MolFromMolBlock(sdf_block, removeHs=False, sanitize=False)
        if m is None:
            return None
        smi = Chem.MolToSmiles(m, isomericSmiles=True, allHsExplicit=True)
        smi = (smi or "").strip()
        return smi or None
    except Exception:
        return None


# -------------------- Arrays for HSR --------------------
def arr_from_ccdc_mol(ccdcmol: Molecule) -> np.ndarray:
    arr = np.array([[a.coordinates[0], a.coordinates[1], a.coordinates[2], np.sqrt(a.atomic_number)]
                    for a in ccdcmol.atoms])
    arr -= arr.mean(axis=0)
    return arr


def arr_from_pybel_mol(pybelmol) -> np.ndarray:
    arr = np.array([[a.coords[0], a.coords[1], a.coords[2], np.sqrt(a.atomicnum)]
                    for a in pybelmol.atoms])
    arr -= arr.mean(axis=0)
    return arr


# -------------------- Timeout-protected CCDC generator --------------------
def _ccdc_worker(smiles: str, q, run_seed: int):
    try:
        random.seed(run_seed)
        cg = ConformerGenerator()
        cg.settings.max_conformers = random.randint(2, 100)
        mol0 = Molecule.from_string(smiles)
        hits = cg.generate(mol0).hits
        q.put(hits[0].molecule.to_string("sdf") if hits else None)
    except Exception:
        q.put(None)


def generate_ccdc_sdf(smiles: str, timeout: int, run_seed: int) -> Optional[str]:
    ctx = multiprocessing.get_context("fork")
    q = ctx.Queue()
    p = ctx.Process(target=_ccdc_worker, args=(smiles, q, run_seed))
    p.start()
    p.join(timeout)
    if p.is_alive():
        p.kill()
        p.join()
    return q.get() if not q.empty() else None


# -------------------- 3D builders (with configurable H addition) --------------------
def gen_with_ccdc(smiles: str, run_seed: int):
    sdf = generate_ccdc_sdf(smiles, timeout=CCDC_TIMEOUT, run_seed=run_seed)
    if sdf is None:
        raise RuntimeError("CCDC generation failed or timed out")
    mol = Molecule.from_string(sdf)
    return mol, sdf


def gen_with_rdkit(smiles: str, run_seed: int, add_hs: bool = True):
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        raise ValueError("RDKit MolFromSmiles failed")
    if add_hs:
        m = Chem.AddHs(m)
    params = AllChem.ETKDG()
    params.randomSeed = int(run_seed) & 0xffffffff
    AllChem.EmbedMolecule(m, params)
    AllChem.MMFFOptimizeMolecule(m)
    return m, Chem.MolToMolBlock(m)


def gen_with_obabel(smiles: str, run_seed: int, add_hs: bool = True):
    os.environ["OB_RANDOM_SEED"] = str(run_seed)
    pm = pb.readstring("smi", smiles)
    if add_hs:
        pm.addh()
    pm.make3D()
    pm.localopt()
    return pm, pm.write("sdf")


# -------------------- Robust SMILES atom counting --------------------
def _manual_smiles_atom_count(smi: str) -> Optional[int]:
    if not smi:
        return None
    i, n = 0, 0
    L = len(smi)

    two_letter = {
        "Cl","Br","Si","As","Se","Li","Na","Al","Ca","Ir","Pt","Pd","Fe","Co","Ni","Cu","Zn","Ag",
        "Sn","Sb","Xe","Hg","Pb","Bi","Mg","Mn","Mo","Cr","Re","Os","Ru","Rh","Ga","Ge","Hf","Ta",
        "Te","Zr","Sc","La","Ce","Pr","Nd","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Th","Pa",
        "se","as"
    }
    single_upper = set("BCNOPSFIKVYWHU")
    single_aromatic = set("bcnops")

    while i < L:
        ch = smi[i]
        if ch == '[':
            j = i + 1
            depth = 1
            while j < L and depth > 0:
                if smi[j] == '[':
                    depth += 1
                elif smi[j] == ']':
                    depth -= 1
                j += 1
            n += 1
            i = j
            continue

        if ch.isdigit() or ch in '()=#$-+@\\/.:*%':
            i += 1
            continue

        if i + 1 < L:
            tok2 = smi[i:i+2]
            if tok2 in two_letter:
                n += 1
                i += 2
                continue

        if ch in single_upper or ch in single_aromatic:
            n += 1
        i += 1

    return n


def count_atoms_from_smiles(smi: str) -> Optional[int]:
    if not smi:
        return None
    try:
        m = Chem.MolFromSmiles(smi, sanitize=False)
        if m is not None:
            return int(m.GetNumAtoms())
    except Exception:
        pass
    try:
        pm = pb.readstring("smi", smi)
        if pm is not None:
            return int(len(list(pm.atoms)))
    except Exception:
        pass
    try:
        return _manual_smiles_atom_count(smi)
    except Exception:
        return None


# =========================== Main ============================
def main():
    ap = argparse.ArgumentParser(
        description="Extract SMILES (CSD + RDKit explicit-H), generate 3D (ccdc/rdkit/obabel), and evaluate vs CSD."
    )
    ap.add_argument("--outdir", default="SMILES_and_3D_out", help="Base output directory")
    ap.add_argument("--seed", type=int, default=None, help="Random seed for 3D generation")

    # Control H addition (defaults: add H for RDKit and OBabel)
    ap.add_argument("--no-rdkit-addhs", dest="rdkit_addhs", action="store_false",
                    help="Disable AddHs in RDKit generation")
    ap.add_argument("--no-obabel-addhs", dest="obabel_addhs", action="store_false",
                    help="Disable addh in Open Babel generation")
    ap.set_defaults(rdkit_addhs=True, obabel_addhs=True)

    args = ap.parse_args()

    base_out = os.path.abspath(args.outdir)
    ensure_dir(base_out)

    run_seed = args.seed if args.seed is not None else np.random.randint(0, 2**32 - 1)
    print(f"==== RUN SEED: {run_seed} ====")
    print(f"RDKit AddHs: {args.rdkit_addhs} | OBabel addh: {args.obabel_addhs}")

    reader = io.EntryReader("CSD")

    results_rows: List[Dict] = []
    atom_rows_all: List[Dict] = []

    for i, refcode in enumerate(TARGETS, 1):
        target_label = f"{i}_{refcode}"
        print(f"\n[{i}/{len(TARGETS)}] Target: {target_label}")

        target_dir = os.path.join(base_out, target_label)
        ensure_dir(target_dir)

        results_rows_target: List[Dict] = []

        try:
            entry = reader.entry(refcode)
            comp = component_of_interest(entry.molecule)
            if comp is None:
                raise RuntimeError("No suitable component found")
        except Exception as e:
            print(f"  ❌ Could not load main component: {e}")
            fail_row = {
                "Target": target_label, "SMILES_Source": "-", "SMILES": "-",
                "Method": "-", "Time_s": None, "Similarity": None, "OK": "–", "SDF_Path": ""
            }
            results_rows.append(fail_row)
            results_rows_target.append(fail_row)
            atom_rows_all.append({
                "Target": target_label, "Category": "reference",
                "SMILES_Source": "", "Method": "-", "Atom_Count": None
            })
            with open(os.path.join(target_dir, "smiles_sources.txt"), "w") as fh:
                fh.write("Could not load main component\n")
            pd.DataFrame([{"SMILES_Source":"-", "SMILES":""}]).to_csv(
                os.path.join(target_dir, "smiles_sources.csv"), index=False
            )
            continue

        ref_sdf = comp.to_string("sdf")

        arr_ref = arr_from_ccdc_mol(comp)
        fp_ref = fp.generate_fingerprint_from_data(arr_ref)

        n_ref = len(comp.atoms)
        atom_rows_all.append({
            "Target": target_label, "Category": "reference",
            "SMILES_Source": "", "Method": "-", "Atom_Count": int(n_ref)
        })

        csd_smiles = smiles_from_csd(comp)
        rdk_smiles_expl = smiles_from_rdkit_explicit_from_sdf(ref_sdf)

        smiles_sources: List[Tuple[str, Optional[str]]] = [
            ("CSD_SMILES", csd_smiles),
            ("RDKit_SMILES_EXPL", rdk_smiles_expl),
        ]

        # Save SMILES files in the target folder (CSV)
        pd.DataFrame(
            [{"SMILES_Source": name, "SMILES": smi or ""} for name, smi in smiles_sources]
        ).to_csv(os.path.join(target_dir, "smiles_sources.csv"), index=False)
        print(f"  💾 Saved SMILES list in {target_dir}")

        # Record SMILES atom counts
        for src_name, smi in smiles_sources:
            n_smi = count_atoms_from_smiles(smi) if smi else None
            atom_rows_all.append({
                "Target": target_label, "Category": "smiles_parsed",
                "SMILES_Source": src_name, "Method": "-", "Atom_Count": n_smi
            })

        # Generate 3D structures and SAVE SDFs under each target/src folder
        for src_name, smi in smiles_sources:
            src_dir = os.path.join(target_dir, src_name)
            ensure_dir(src_dir)

            if not smi:
                print(f"  ⚠️  {src_name}: unavailable; skipping 3D generation")
                row = {
                    "Target": target_label, "SMILES_Source": src_name, "SMILES": "",
                    "Method": "-", "Time_s": None, "Similarity": None, "OK": "–", "SDF_Path": ""
                }
                results_rows.append(row)
                results_rows_target.append(row)
                continue

            print(f"  → {src_name}: {smi}")

            gens = {
                "ccdc":   lambda s, seed: gen_with_ccdc(s, seed),
                "rdkit":  lambda s, seed: gen_with_rdkit(s, seed, add_hs=args.rdkit_addhs),
                "obabel": lambda s, seed: gen_with_obabel(s, seed, add_hs=args.obabel_addhs),
            }

            for mname, gen_func in gens.items():
                t0 = time.time()
                sdf_path = os.path.join(src_dir, f"{refcode}__{src_name}__{mname}.sdf")
                try:
                    mol3d, sdf = gen_func(smi, run_seed)
                    elapsed = time.time() - t0

                    if mname == "obabel":
                        arr = arr_from_pybel_mol(mol3d)
                        n_gen = int(len(mol3d.atoms))
                    elif mname == "rdkit":
                        arr = pp.molecule_to_ndarray(mol3d, features=PROTON_FEATURES, removeHs=False)
                        n_gen = int(mol3d.GetNumAtoms())
                    else:
                        arr = arr_from_ccdc_mol(mol3d)
                        n_gen = int(len(mol3d.atoms))

                    fp_gen = fp.generate_fingerprint_from_data(arr)
                    sim_val = sim.compute_similarity_score(fp_ref, fp_gen)

                    # SAVE the generated SDF
                    with open(sdf_path, "w") as fh:
                        fh.write(sdf)

                    row = {
                        "Target": target_label,
                        "SMILES_Source": src_name,
                        "SMILES": smi,
                        "Method": mname,
                        "Time_s": elapsed,
                        "Similarity": float(sim_val),
                        "OK": "✓",
                        "SDF_Path": sdf_path
                    }
                    results_rows.append(row)
                    results_rows_target.append(row)

                    atom_rows_all.append({
                        "Target": target_label, "Category": "generated",
                        "SMILES_Source": src_name, "Method": mname,
                        "Atom_Count": n_gen
                    })

                    print(f"    {mname.upper():6}  time={elapsed:5.2f}s  sim={sim_val:5.2f}  → {os.path.relpath(sdf_path, base_out)}")

                except Exception as e:
                    elapsed = time.time() - t0
                    row = {
                        "Target": target_label,
                        "SMILES_Source": src_name,
                        "SMILES": smi,
                        "Method": mname,
                        "Time_s": None,
                        "Similarity": None,
                        "OK": "–",
                        "SDF_Path": "",
                    }
                    results_rows.append(row)
                    results_rows_target.append(row)
                    print(f"    {mname.upper():6}  FAILED after {elapsed:4.1f}s: {e}")

        # (Optional)
        # pd.DataFrame(results_rows_target)[["SMILES_Source","Method","Similarity","Time_s","OK","SDF_Path"]] \
        #   .to_csv(os.path.join(target_dir, "summary.csv"), index=False)


    # 1) similarity comparison 
    results_df = pd.DataFrame(results_rows)
    try:
        pivot = results_df.pivot_table(
            index="Target",
            columns=["SMILES_Source", "Method"],
            values="Similarity",
            aggfunc="first"
        ).sort_index(axis=1)
        order_labels = [f"{i}_{rc}" for i, rc in enumerate(TARGETS, 1)]
        pivot = pivot.reindex(order_labels)
        pivot_fmt = pivot.applymap(fmt_num)

        sim_cmp_path = os.path.join(base_out, "similarity_comparison.csv")
        pivot_fmt.to_csv(sim_cmp_path)
        print(f"\n💾 Saved similarity comparison → {sim_cmp_path}")
    except Exception as e:
        print(f"\n⚠️ Could not build/save similarity comparison: {e}")

    # 2) atom-count comparison 
    try:
        atom_all_df = pd.DataFrame(atom_rows_all)

        idx = [f"{i}_{rc}" for i, rc in enumerate(TARGETS, 1)]
        cols = [
            "Orig_Atoms",
            "SMILES_CSD",
            "SMILES_RDKit_explicit",
            "CCDC_from_CSD",
            "RDKit_from_CSD",
            "OBabel_from_CSD",
            "CCDC_from_Explicit",
            "RDKit_from_Explicit",
            "OBabel_from_Explicit",
        ]
        comp = pd.DataFrame(index=idx, columns=cols, dtype=float)
        ref_df = atom_all_df[atom_all_df["Category"] == "reference"]
        smi_df = atom_all_df[atom_all_df["Category"] == "smiles_parsed"]
        gen_df = atom_all_df[atom_all_df["Category"] == "generated"]

        comp.loc[ref_df["Target"], "Orig_Atoms"] = ref_df.set_index("Target")["Atom_Count"]
        comp.loc[smi_df[smi_df["SMILES_Source"] == "CSD_SMILES"]["Target"], "SMILES_CSD"] = \
            smi_df[smi_df["SMILES_Source"] == "CSD_SMILES"].set_index("Target")["Atom_Count"]
        comp.loc[smi_df[smi_df["SMILES_Source"] == "RDKit_SMILES_EXPL"]["Target"], "SMILES_RDKit_explicit"] = \
            smi_df[smi_df["SMILES_Source"] == "RDKit_SMILES_EXPL"].set_index("Target")["Atom_Count"]

        for method, colname in [("ccdc", "CCDC_from_CSD"), ("rdkit", "RDKit_from_CSD"), ("obabel", "OBabel_from_CSD")]:
            sel = gen_df[(gen_df["SMILES_Source"] == "CSD_SMILES") & (gen_df["Method"] == method)]
            comp.loc[sel["Target"], colname] = sel.set_index("Target")["Atom_Count"]

        for method, colname in [("ccdc", "CCDC_from_Explicit"), ("rdkit", "RDKit_from_Explicit"), ("obabel", "OBabel_from_Explicit")]:
            sel = gen_df[(gen_df["SMILES_Source"] == "RDKit_SMILES_EXPL") & (gen_df["Method"] == method)]
            comp.loc[sel["Target"], colname] = sel.set_index("Target")["Atom_Count"]

        comp_pretty = comp.applymap(fmt_count)
        comp_pretty_csv = os.path.join(base_out, "atom_counts_comparison.csv")
        comp_pretty.to_csv(comp_pretty_csv)

        print(f"💾 Saved atom counts comparison → {comp_pretty_csv}")

    except Exception as e:
        print(f"\n⚠️ Could not build/save atom-count comparison table: {e}")

    print("\n✅ Finished.")


if __name__ == "__main__":
    main()
