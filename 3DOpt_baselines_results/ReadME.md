# 3DOpt Baselines Results

This folder provides scripts for aggregating, analyzing, and visualizing the results of 3DOpt benchmarks using two classes of baseline methods: Random Sampler and ChemGE (evolutionary algorithm). The workflow assumes you have already collected results from multiple independent runs (typically 10) for each method/generator combination, organized in the folders below. The scripts produce easy-to-plot per-task averages and publication-quality SVG figures.

## Obtaining Baselines Results

The complete set of baseline results for the 3DOpt benchmark (all runs, methods, and generators) is available as compressed archives on [Zenodo](https://zenodo.org/records/16602642), DOI: 10.5281/zenodo.16602641.

To download all results files directly to this folder, you can use `wget` as follows:

```bash
wget -P . -i zenodo_links.txt
```

Once the archives are downloaded, extract them with:

```bash
tar -xzf Results_ChemGE_CCDC.tar.gz
tar -xzf Results_ChemGE_OBabel.tar.gz
tar -xzf Results_ChemGE_RDKit.tar.gz
tar -xzf Results_RndSam_CCDC_Entries.tar.gz
tar -xzf Results_RndSam_CCDC_Smiles.tar.gz
tar -xzf Results_RndSam_OBabel.tar.gz
tar -xzf Results_RndSam_RDKit.tar.gz
```

## Scripts

### `aggregate_3DOpt_results.py`  
Aggregates all `results.csv` files for each method, averages the 3DOpt scores by task, and outputs easy-to-plot CSVs for each method.  
**Usage:**  
```bash
python aggregate_3DOpt_results.py
```
**Inputs:**

Results folders:

- `Results_RndSam_CCDC_Entries`, `Results_RndSam_CCDC_Smiles`, `Results_RndSam_OBabel`, `Results_RndSam_RDKit`

- `Results_ChemGE_CCDC`, `Results_ChemGE_OBabel`, `Results_ChemGE_RDKit`

**Outputs:**

- `Rnd_ccdc_entries.csv`
- `Rnd_ccdc.csv`
- `Rnd_obabel.csv`
- `Rnd_rdkit.csv`
- `ChemGE_ccdc.csv`
- `ChemGE_obabel.csv`
- `ChemGE_rdkit.csv`

### `plot_RndSam.py`
Generates a grouped bar plot of the averaged 3DOpt scores for the Random Sampler runs using different 3D generators.  
**Usage:**  
```bash
python plot_RndSam.py
```
**Inputs:**
- `Rnd_ccdc_entries.csv`, `Rnd_ccdc.csv`, `Rnd_obabel.csv`, `Rnd_rdkit.csv`

**Outputs:**
- `3DOpt_RndSam_Scores.svg` (grouped barplot comparing all methods across all tasks)

### `plot_ChemGE.py`
Generates a grouped bar plot of the averaged 3DOpt scores for ChemGE-based runs with different conformer generators.  
**Usage:**  
```bash
python plot_ChemGE.py
```
**Inputs:**
- `ChemGE_ccdc.csv`, `ChemGE_obabel.csv`, `ChemGE_rdkit.csv`

**Outputs:**
- `3DOpt_ChemGE_Scores.svg` (grouped barplot comparing ChemGE with all 3D generators across tasks)

### `get_success_rate_RndSam.py` and `get_success_rate_ChemGE.py`
Scripts for calculating and summarizing success rates of baseline runs.
They report, for each run and task, the number of structures scored, number of zero-score entries, total expected molecules, and the overall and per-run success rates.
**Usage:**  
```bash
python get_success_rate_RndSam.py Results_RndSam_RDKit/ [--per-task]
python get_success_rate_ChemGE.py Results_ChemGE_OBabel/ [--per-task]
```
The `--per-task` option (optional) prints detailed per-task statistics.

### `plot_delta_scores.py`
Scripts to compute and visualize the differences in scores between the random sampler and the ChemGE baselines.
**Usage:**  
```bash
python plot_delta_scores.py
```

## Notes

- Edit the list of results folders in `aggregate_3DOpt_results.py` if your folder names differ.

- All scripts output publication-ready SVG figures by default.

- Make sure all input results.csv files contain the columns `task` and `3DOpt_Score`.

- You can re-run the aggregation and plotting scripts as needed when new runs are added.