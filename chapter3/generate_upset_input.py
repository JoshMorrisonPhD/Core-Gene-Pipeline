# the input for this script is the deduplicated_full_matrix.xlsx created by merge2.py
# the output is a .txt file showing the number of essential genes for each intersection of species
# this script was used to generate the inputs for essential_upset.R

import pandas as pd
import itertools

df = pd.read_excel("deduplicated_full_matrix.xlsx")

species_columns = {
    "S.Pneumoniae": "pneumo",
    "S.Suis": "suis",
    "S.Agalactiae": "agal",
    "S.Equi_sb_equi": "equi",
    "S.Iniae": "iniae",
    "S.Uberis": "uberis"
}

column_to_species = {v: k for k, v in species_columns.items()}

available_columns = [col for col in df.columns if col in column_to_species]
species_names = [column_to_species[col] for col in available_columns]

input_dict = {}

for i in range(1, len(available_columns) + 1):
    for combo in itertools.combinations(available_columns, i):
        combo_species = [column_to_species[c] for c in combo]
        key = "&".join(combo_species)
        # Rows where all combo columns have non-empty, non-whitespace values
        valid_rows = df[list(combo)].dropna()
        valid_rows = valid_rows[valid_rows.apply(lambda row: all(str(x).strip() != "" for x in row), axis=1)]
        # For 1 column: count unique entries
        if len(combo) == 1:
            col = combo[0]
            unique_genes = set(valid_rows[col].astype(str).str.strip())
            input_dict[key] = len(unique_genes)
        else:
            # For multi-column: count unique tuples across those columns
            unique_pairs = set(tuple(str(row[c]).strip() for c in combo) for _, row in valid_rows.iterrows())
            input_dict[key] = len(unique_pairs)

# Write to .txt file
with open("upset_input_dataset.txt", "w") as f:
    f.write("# Dataset\n")
    f.write("input <- c(\n")
    for idx, (k, v) in enumerate(input_dict.items()):
        comma = "," if idx < len(input_dict) - 1 else ""
        f.write(f'  "{k}" = {v}{comma}\n')
    f.write(")\n")

print("upset_input_dataset.txt generated.")
