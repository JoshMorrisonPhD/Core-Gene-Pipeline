#this scripts inputs are the essential gene spreadsheet "deduplicated_full_matrix.xlsx"
# the output is a single spreadsheet showing the 62 superpangenome core and essential genes' tags for each species
# this script was used to get a list of the 62 superpangenome core and essential genes, which could then be extracted from any of the species by using it as a keyfile for fastafinder_v2.py

import pandas as pd

# Load the combined matrix
df = pd.read_excel("deduplicated_full_matrix.xlsx")

# Specify the columns in the order you want (each as a species)
species_columns = ['equi', 'iniae', 'uberis', 'pneumo', 'suis', 'agal']

# Filter to rows where every species column is filled (not NaN and not blank)
core_rows = df.dropna(subset=species_columns)
core_rows = core_rows[core_rows.apply(lambda row: all(str(row[c]).strip() != "" for c in species_columns), axis=1)]

# Drop duplicate ortholog groups if any (shouldn't happen, but just in case)
core_rows = core_rows.drop_duplicates(subset=species_columns)

# Save as Excel file, each column = species, each row = gene
core_rows[species_columns].to_excel("core_orthologous_groups_all_6_species.xlsx", index=False)

print(f"Saved {len(core_rows)} core groups to all_essential.xlsx")
