#!/usr/bin/env python3
import os
from hsr.pre_processing import load_molecules_from_sdf
from hsr.fingerprint import generate_fingerprint_from_molecule
from hsr.similarity import compute_similarity_score
from hsr.utils import PROTON_FEATURES

cwd = os.getcwd()

def first_mol(sdf_path, removeHs=False):
    mols = load_molecules_from_sdf(sdf_path, removeHs=removeHs, sanitize=False)
    if not mols:
        raise ValueError(f"No molecules found in {sdf_path}")
    return mols[0]

def get_fp(mol, chirality: bool):
    """HSR API returns (fp, extra) when chirality=True; just fp otherwise."""
    out = generate_fingerprint_from_molecule(mol, features=PROTON_FEATURES, chirality=chirality)
    return out[0] if chirality else out

# Load molecules
mols = [
    ("Target (CSD)", first_mol(f"{cwd}/OFOWIS.sdf")),
    ("CCDC",          first_mol(f"{cwd}/OFOWIS_ccdc.sdf")),
    ("OBabel",        first_mol(f"{cwd}/OFOWIS_obabel.sdf")),
    ("RDKit",         first_mol(f"{cwd}/OFOWIS_rdkit.sdf")),
]

for chirality in (True, False):
    label = "ON " if chirality else "OFF"
    print(f"\nChirality: {label}")
    print(f"{'Method':<14} {'HSR similarity vs Target':>26}")
    print("-" * 40)

    # Build fingerprints in the same order; first entry is the reference
    fps = [get_fp(m, chirality) for _, m in mols]
    ref_fp = fps[0]

    for (name, _), fp in zip(mols, fps):
        sim = compute_similarity_score(ref_fp, fp)
        print(f"{name:<14} {sim:>26.4f}")
