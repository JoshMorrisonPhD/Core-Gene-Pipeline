# the inputs for this file are the $species_presence_matrix.xlsx outputs from generate_all_presence_matrices_V2.py as well as the "no_results_all_species.xlsx" file which was made manually by just copy pasting all the 
# blast query ID's with no hits from each species into a spreadsheet with the same columns/column order,and leaving the rest of the row empty for each query
# the ouput is a single deduplicated_full_matrix.xlsx file
# this script creates the essential gene presence/absence spreadsheet used as input for generate_upset_input.py and for essential_all_extractor.py

import pandas as pd
import os

species_files = {
    "equi": "equi_presence_matrix.xlsx",
    "iniae": "iniae_presence_matrix.xlsx",
    "uberis": "uberis_presence_matrix.xlsx",
    "pneumo": "pneumo_presence_matrix.xlsx",
    "suis": "suis_presence_matrix.xlsx",
    "agal": "agal_presence_matrix.xlsx"
}

# Get all possible columns in order
all_cols = ['equi', 'iniae', 'uberis', 'pneumo', 'suis', 'agal']

# Collect all rows from all files
rows = []
for file in species_files.values():
    df = pd.read_excel(file)
    # Ensure all columns exist in correct order
    for col in all_cols:
        if col not in df.columns:
            df[col] = pd.NA
    df = df[all_cols]
    # Store every row as a tuple (for perfect uniqueness)
    rows.extend([tuple(row) for row in df.values])

# Also add from no_results
no_results_file = "no_results_all_species.xlsx"
if os.path.exists(no_results_file):
    df = pd.read_excel(no_results_file)
    for col in all_cols:
        if col not in df.columns:
            df[col] = pd.NA
    df = df[all_cols]
    rows.extend([tuple(row) for row in df.values])

# Convert to DataFrame, dropping only exact duplicate rows (not per-gene)
merged_df = pd.DataFrame(rows, columns=all_cols).drop_duplicates()

# Save final result
merged_df.to_excel("deduplicated_full_matrix.xlsx", index=False)
print("Merged file with all mapping contexts saved to: deduplicated_full_matrix.xlsx")