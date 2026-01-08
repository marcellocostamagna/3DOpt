# Filtering

This folder contains the script `filter_csd.py` for filtering the Cambridge Structural Database (CSD) to retrieve all structures that meet strict suitability criteria for downstream benchmarking.

## Filtering criteria

Structures are retained only if they:

- Have 3D coordinates
- Are not polymeric
- Have explicit hydrogens
- Have a "component of interest" with more than 5 atoms
- The component of interest has no disordered or overlapping atoms
- The component of interest has a valid SMILES string

## Usage

Run the script with:

```bash
python filter_csd.py