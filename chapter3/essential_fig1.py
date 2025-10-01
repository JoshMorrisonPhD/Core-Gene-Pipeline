# input is hardcoded essential gene numbers compared to total number of genes in reference genome
# output is a stacked bar chart of essential and non essential genes in each species
# script was used to generate Figure 3.1 chapter 3

import matplotlib.pyplot as plt
import numpy as np

# Example data (replace with your actual numbers and species names if needed)
species = ['S. agalactiae', 'S. equi sbsp equi', 'S. iniae', 'S. pneumoniae', 'S. suis', 'S. uberis']

essential_genes = np.array([160, 401, 481, 536, 342, 384])
non_essential_genes = np.array([1970, 1764, 1535, 1664, 1615, 1489])

x = np.arange(len(species))
width = 0.6

fig, ax = plt.subplots(figsize=(10, 6))

# Plot stacked bars
p1 = ax.bar(x, essential_genes, width, label='Essential genes', color='#D33F6A')  # Red
p2 = ax.bar(x, non_essential_genes, width, bottom=essential_genes, label='Non-essential genes', color='#7F7F7F')  # Grey

# Annotate inside bars
for idx in range(len(species)):
    # Essential genes
    ax.text(x[idx], essential_genes[idx]/2, str(essential_genes[idx]),
            ha='center', va='center', color='white', fontweight='bold', fontsize=9)
    # Non-essential genes
    total_height = essential_genes[idx] + non_essential_genes[idx]
    ax.text(x[idx], essential_genes[idx] + non_essential_genes[idx]/2, str(non_essential_genes[idx]),
            ha='center', va='center', color='white', fontweight='bold', fontsize=9)
    # Total genes above bar
    ax.text(x[idx], total_height + 30, f'Total: {total_height}',
            ha='center', va='bottom', fontsize=9, fontweight='bold')

# Labels and aesthetics
ax.set_ylabel('Number of Genes')
ax.set_xlabel('Species')
ax.set_title('Essential and Non-Essential Genes in Different Streptococcus Species')
ax.set_xticks(x)
ax.set_xticklabels(species, rotation=15)
ax.legend(loc='upper left', bbox_to_anchor=(1.02, 1))

plt.tight_layout()

# Save high-resolution for thesis
plt.savefig("streptococcus_essential_nonessential_thesis_palette.png", dpi=600)
plt.savefig("streptococcus_essential_nonessential_thesis_palette.svg")

plt.show()
