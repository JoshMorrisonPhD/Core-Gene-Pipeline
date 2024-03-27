import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from matplotlib_venn import venn2
import argparse
import sys

def get_args():
    try:
        parser = argparse.ArgumentParser(
            description="Generates histogram and venn diagrams for a given species."
        )
        parser.add_argument(
            "-s",
            "--species",
            action="store",
            required=True,
            help="The name of the species, e.g equi or agal",
        )

        if len(sys.argv) == 1:
            parser.print_help(sys.stderr)
            sys.exit(1)

    except NameError:
        sys.stderr.write(
            "An exception occured with argument parsing. Check your provided options."
        )
        sys.exit(1)

    return parser.parse_args()

args = get_args()
species_name = args.species

# Create a subdirectory for the output
output_dir = 'Species_Output'
os.makedirs(output_dir, exist_ok=True)

# Define the paths to the uploaded recap files for your species
recap_files = {
    'agalcoredb': f'{species_name}_essential_vs_agalcoredb_recap.txt',
    'allcoredb': f'{species_name}_essential_vs_allcoredb_recap.txt',
    'equicoredb': f'{species_name}_essential_vs_equicoredb_recap.txt',
    'iniaecoredb': f'{species_name}_essential_vs_iniaecoredb_recap.txt',
    'pneumocoredb': f'{species_name}_essential_vs_pneumocoredb_recap.txt',
    'suiscoredb': f'{species_name}_essential_vs_suiscoredb_recap.txt',
    'uberiscoredb': f'{species_name}_essential_vs_uberiscoredb_recap.txt'
}

# Initialize a dictionary to store the data
data_dict1 = {'Database': [], 'Genes with Hits': [], 'Genes with Zero Hits': []}

# Read each uploaded recap file and extract the required information
for db, path in recap_files.items():
    with open(path, 'r') as file:
        lines = file.readlines()
        hits = int(lines[0].strip().split(": ")[1])
        zero_hits = int(lines[1].strip().split(": ")[1])
        data_dict1['Database'].append(db)
        data_dict1['Genes with Hits'].append(hits)
        data_dict1['Genes with Zero Hits'].append(zero_hits)

# Convert the data to a DataFrame for easier plotting
df_data_dict1 = pd.DataFrame(data_dict1)

# Plotting the histogram
fig, ax = plt.subplots(figsize=(12, 8))
width = 0.35
x = np.arange(len(df_data_dict1))

bars1 = ax.bar(x - width/2, df_data_dict1['Genes with Hits'], width, label='Essential and Core', color='lightblue')
bars2 = ax.bar(x + width/2, df_data_dict1['Genes with Zero Hits'], width, label='Essential but not Core', color='darkblue')

ax.set_xlabel('Database')
ax.set_ylabel('Number of Genes')
ax.set_title(f'Essential vs Core gene comparison in Streptococcus {species_name.capitalize()}')
ax.set_xticks(x)
ax.set_xticklabels(df_data_dict1['Database'], rotation=45, ha="right")
ax.legend()

# Adding the numbers of genes on the bars for readability
for bars in [bars1, bars2]:
    for bar in bars:
        height = bar.get_height()
        ax.annotate('{}'.format(height),
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 3),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom')

plt.tight_layout(rect=[0,0,0.85,1])
plt.savefig(f'{output_dir}/{species_name}_histogram.png')
plt.close()

# Total genes in each database
total_genes = {
    'agalcoredb': 1556,
    'allcoredb': 716,
    'equicoredb': 1781,
    'iniaecoredb': 1561,
    'pneumocoredb': 1459,
    'suiscoredb': 1279,
    'uberiscoredb': 1422
}

# Function to plot and save Venn diagrams
def plot_and_save_venn(db_name, essential_genes, core_genes, essential_and_core, output_dir):
    plt.figure(figsize=(8, 6))
    venn2(subsets=(essential_genes - essential_and_core, core_genes, essential_and_core),
          set_labels=('Essential Genes', f'{db_name.split("_")[0].capitalize()} Core Genes'))
    plt.title(f'{db_name} Essential vs Core Gene Comparison')
    plt.savefig(f'{output_dir}/{db_name}_venn.png')
    plt.close()

# Generate and save Venn diagrams
for db in df_data_dict1['Database']:
    essential_genes = df_data_dict1.loc[df_data_dict1['Database'] == db, ['Genes with Hits', 'Genes with Zero Hits']].sum(axis=1).values[0]
    core_genes = total_genes[db]
    essential_and_core = df_data_dict1.loc[df_data_dict1['Database'] == db, 'Genes with Hits'].values[0]
    plot_and_save_venn(db, essential_genes, core_genes, essential_and_core, output_dir)

print(f"Histogram and Venn diagrams saved in '{output_dir}' directory.")
