# 3DOpt

**3DOpt** is the first benchmark for molecular de novo design to test generative methods for their ability to handle and optimize both organic and inorganic molecules in 3D. This benchmark is now part of the benchmarking framework [MolScore](https://github.com/MorganCThomas/MolScore).

This repository contains the scripts, workflows, and analyses developed for 3DOpt, as well as supporting work presented in the associated article ([TODO: add reference]).  
It includes workflows for filtering and preparing chemical data (`Filtering`, `Starting_populations`), scripts to test and compare such data (`CAEs`, `CAEs_analysis`, `Metal_targets_statistics`), tests on tools used for baseline methods (`Generators_analysis`), scripts to generate configuration files for MolScore (`3DOpt_configuration`), and scripts to aggregate, analyze and visualize baseline results (`3DOpt_baselines_results`).

## Set Up
To run any script in the repo, first create and activate the provided conda environment:

```bash
conda env create -f environment.yml
conda activate 3DOpt_devel
```
*Note: Some scripts require a valid CCDC license*

## Overview
Each folder contains its own README describing the aim and usage of the scripts. Below is a brief overview:

- `Filtering`:  
Scripts for selecting structures of interest from the Cambridge Structural Database (CSD).

- `Starting_populations`:  
Scripts to generate and analyze starting populations after defining a set of targets.

- `CAEs`:  
Scripts to generate and compare Connected Atom Environments (CAEs) for targets and candidate starting populations (to check if the starting set contains the building blocks needed to reconstruct targets). 

- `CAEs_analysis`:  
Analysis scripts for CAE results to inform starting population selection.

- `Metal_targets_statistics`:  
Scripts for statistical analysis and visualization of metal types and geometry diversity among inorganic targets.

- `Generators_analysis`: 
Benchmarking of the three 3D structure generators used  (`CCDC`, `OBabel`, `RDKit`) for the chosen 3DOpt targets.

- `3DOpt_configuration`:  
Scripts to generate MolScore configuration files for 3DOpt, as well as 2D diagrams for the selected targets.

- `3DOpt_baselines_results`:  
Scripts to aggregate, analyze, and plot results from 3DOpt runs using two baseline methods: `RandomSampler` and `ChemGE`. For further information, see [Molscore_baselines](https://github.com/MorganCThomas/MolScore_baselines).


## References 
For more details and relevant citations, see the article (TODO: add reference), [MolScore](https://github.com/MorganCThomas/MolScore), its [publication](https://doi.org/10.1186/s13321-024-00861-w), and [MolScore_baselines](https://github.com/MorganCThomas/MolScore_baselines).

For questions or issues, open an issue in this repo or contact the corresponding author.