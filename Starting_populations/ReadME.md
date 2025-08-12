# Starting populations
This folder contains scripts to generate and analyze starting populations of molecules for each target in the 3DOpt benchmark. The populations are selected from the `viable_structures.csv` file based on a similarity threshold: only molecules below the specified similarity threshold to the target are included.

## Scripts

The scripts should be run in the following order:  

- `get_starting_populations.py`:  
For each target, collects all molecules (CSD entry code and SMILES) from `viable_structures.csv` that are below the similarity threshold with respect to the target structure. Outputs starting population files (one per target) in a dedicated folder for each threshold.

- `parse_populations_logs.py`:  
Parse the logs generated from the previous step for three similarity thresholds (0.3, 0.4, 0.5), and summarizes the number of selected molecules per target in `Starting_populations.csv`.

- `plot_starting_populations_numbers.py`:  
Visualizes the summary statistics from `Starting_populations.csv` to compare the effect of different thresholds on starting population size.

- `get_similarity_distributions.py`:
Builds HSR-similarity distributions for the starting populations for each threshold folder (0.3, 0.4, and 0.5). The outputs are are organized under the `Similarity_distributions`folder.

## Usage 
All scripts in this folder can be run directly with:  
```bash
python script.py
```
Replace `script.py` with the name of the script you want to run. Make sure that output from previous steps (e.g., `viable_structures.csv` and log files) is available in the expected locations.