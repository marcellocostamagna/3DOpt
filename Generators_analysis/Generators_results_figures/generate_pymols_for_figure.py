import os
from rdkit import Chem
from rdkit.Geometry import Point3D
import pymol
from pymol import cmd
from pymol.cgo import *

# ========= EDIT HERE ==========
SDF_FILE = "OFOWIS_rdkit.sdf"  # <-- Set this to your filename: OFOWIS.sdf, OFOWIS_ccdc.sdf, OFOWIS_obabel.sdf, OFOWIS_rdkit.sdf
# ==============================

# Auto set carbon exclusion index
CARBON_IDX = 0 if SDF_FILE == "OFOWIS.sdf" else 12

def create_dotted_line(start, end, color=[0.4, 0.4, 0.4], segments=18, thickness=0.045):
    obj = []
    for i in range(segments):
        fraction = i / float(segments)
        next_fraction = (i + 1) / float(segments)
        intermediate_point = [start[j] + (end[j] - start[j]) * fraction for j in range(3)]
        next_intermediate_point = [start[j] + (end[j] - start[j]) * next_fraction for j in range(3)]
        if i % 2 == 0:
            obj.extend([
                CYLINDER,
                intermediate_point[0], intermediate_point[1], intermediate_point[2],
                next_intermediate_point[0], next_intermediate_point[1], next_intermediate_point[2],
                thickness,
                color[0], color[1], color[2],
                color[0], color[1], color[2]
            ])
    return obj

def get_iridium_bond_indices(mol, carbon_idx):
    ir_bonds = []
    for bond in mol.GetBonds():
        a1 = bond.GetBeginAtom()
        a2 = bond.GetEndAtom()
        # Exclude if either atom is Cl (atomic number 17)
        if (a1.GetAtomicNum() == 77 and a2.GetAtomicNum() == 17) or (a2.GetAtomicNum() == 77 and a1.GetAtomicNum() == 17):
            continue
        # Exclude if either atom is C and atom ID == carbon_idx (ID 1 in SDF/PyMOL is 0 in RDKit)
        if (a1.GetAtomicNum() == 77 and a2.GetIdx() == carbon_idx and a2.GetAtomicNum() == 6) or \
           (a2.GetAtomicNum() == 77 and a1.GetIdx() == carbon_idx and a1.GetAtomicNum() == 6):
            continue
        if a1.GetAtomicNum() == 77 or a2.GetAtomicNum() == 77:
            ir_bonds.append((bond.GetIdx(), a1.GetIdx(), a2.GetIdx()))
    return ir_bonds

def get_atom_coords(mol, idx):
    pos = mol.GetConformer().GetAtomPosition(idx)
    return [pos.x, pos.y, pos.z]

# ===== Main Script =====

cwd = os.getcwd()
sdf_path = os.path.join(cwd, SDF_FILE)

# Load with RDKit
mol = Chem.MolFromMolFile(sdf_path, removeHs=False, sanitize=False)

# Initialize PyMOL
pymol.finish_launching(['pymol','-q'])  # "-q" is quiet

# Load molecule in PyMOL
cmd.load(sdf_path, "molecule1")

# 1. Find all Iridium bonds (with exclusions)
ir_bonds = get_iridium_bond_indices(mol, carbon_idx=CARBON_IDX)

# 2. Replace each with a gray dotted line
for bond_idx, a1_idx, a2_idx in ir_bonds:
    a1_pymol = a1_idx + 1
    a2_pymol = a2_idx + 1
    cmd.unbond(f"molecule1 and id {a1_pymol}", f"molecule1 and id {a2_pymol}")
    start = get_atom_coords(mol, a1_idx)
    end = get_atom_coords(mol, a2_idx)
    gray = [0.2, 0.2, 0.2]
    cgo_obj = create_dotted_line(start, end, color=gray, segments=18, thickness=0.045)
    cmd.load_cgo(cgo_obj, f"ir_bond_{a1_pymol}_{a2_pymol}")

cmd.zoom("molecule1")

print(f"Substituted {len(ir_bonds)} Ir–X bonds (except Ir–Cl and Ir–C{CARBON_IDX+1}) with gray dotted lines!")
