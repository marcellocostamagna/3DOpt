# Explicit Hydrogens Analysis

This folder contains scripts for analyzing the impact of implicit vs. explicit hydrogens in SMILES on 3D structure generation.
We compare canonical CSD SMILES (implicit H) with explicit-H SMILES (from RDKit) as inputs to three generators: CCDC, RDKit, and Open Babel.

## Scripts

### `run_smiles_3d_single.py`    
Runs a single analysis run over all targets using both SMILES types and the three generators.
**Usage:**  
```bash
python run_smiles_3d_single.py --outdir OUT --seed 12345
```  
**Outputs:** 
- `similarity_comparison.csv`: per-target similarities to CSD reference (with “N.A.” for failures)
- `atom_counts_comparison.csv`: atom counts for references, SMILES, and generated structures

### `run_smiles_3d_multi.py`  
Repeats the single-run analysis multiple times with different random seeds and aggregates the results.
**Usage:**  
```bash
python run_smiles_3d_multi.py
```
**Outputs:**  
- Mean and standard deviation of similarities
- Success-rate summaries (per target, per run, overall)
- Scatter plots comparing implicit vs. explicit SMILES