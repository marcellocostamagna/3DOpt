#!/usr/bin/env python3
"""
count_init_populations.py
Count the number of non-empty lines in init_pop.txt files across different similarity thresholds.

Reads files from Starting_populations_0_3, Starting_populations_0_4, Starting_populations_0_5
and creates a CSV similar to Starting_populations.csv showing counts for each target.
"""

import csv
import re
from pathlib import Path
from collections import defaultdict

# Define directories and their corresponding column names
directories = {
    "Starting_populations_0_3": "0.3",
    "Starting_populations_0_4": "0.4",
    "Starting_populations_0_5": "0.5",
}

# Pattern to extract index and refcode from filename: <index>_<refcode>_init_pop.txt
filename_pattern = re.compile(r"^(\d+)_([A-Z0-9]+)_init_pop\.txt$", re.IGNORECASE)

table = defaultdict(dict)
order = []  # Will hold targets in first-seen order

for dir_name, column in directories.items():
    dir_path = Path(dir_name)
    
    if not dir_path.is_dir():
        print(f"⚠️  Skipping {dir_name}: directory not found")
        continue
    
    # Process all files in the directory
    for file_path in sorted(dir_path.iterdir()):
        if not file_path.is_file():
            continue
        
        match = filename_pattern.match(file_path.name)
        if not match:
            continue
        
        index, refcode = match.groups()
        index = int(index)
        
        # Count non-empty lines in the file
        try:
            with file_path.open('r') as f:
                count = sum(1 for line in f if line.strip())
        except Exception as e:
            print(f"⚠️  Error reading {file_path.name}: {e}")
            count = 0
        
        # Store count and maintain order
        target = refcode
        if target not in order:
            order.append(target)
        
        table[target][column] = count
        table[target]['index'] = index

# Output CSV
columns = ["0.3", "0.4", "0.5"]
out_file = Path("Starting_populations.csv")

with out_file.open("w", newline="") as fh:
    writer = csv.writer(fh)
    writer.writerow(["target", *columns])
    
    for target in order:
        # Use the index from the first directory processed
        idx = table[target].get('index', '')
        display_name = f"{idx}_{target}"
        row = [display_name] + [table[target].get(col, "") for col in columns]
        writer.writerow(row)

print(f"✅  Wrote {out_file.resolve()}")
print(f"    Processed {len(order)} targets across {len(directories)} directories")

