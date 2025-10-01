# this script has hardcoded inputs: they are the CEGs, core genes without CEGs, essential genes and total genome sizes for each species.
# the outputs of this script are all the statistics for CEGs (chi-square test, Fisher test per-species, and a spearman correlation)
# this script runs all the stats for the CEGs in chapter 4

import numpy  as np
import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact, spearmanr
from statsmodels.stats.multitest import multipletests
from tabulate import tabulate

# Number of CEGs
core_essential = {          
    'agal'  : 160,
    'equi'  : 398,
    'iniae' : 419,
    'pneumo': 370,
    'suis'  : 308,
    'uberis': 324
}
# Number of Core genes that are not essential
core_non_essential = {      
    'agal'  : 1396,
    'equi'  : 1383,
    'iniae' : 1142,
    'pneumo': 1090,
    'suis'  :  971,
    'uberis': 1098
}
# Number of essential genes
total_essential = {         
    'agal'  : 160,
    'equi'  : 401,
    'iniae' : 476,
    'pneumo': 536,
    'suis'  : 342,
    'uberis': 384
}

# Genome sizes
genome_size = {
    'agal'  : 5300,
    'equi'  : 2347,
    'iniae' : 3047,
    'pneumo': 4341,
    'suis'  : 6645,
    'uberis': 6334
}

species = list(core_essential.keys())

# 2 x 6 x² test 
tbl = np.vstack([[core_essential[s]  for s in species],
                 [core_non_essential[s] for s in species]])
chi2, p_chi, dof, _ = chi2_contingency(tbl)

# Per-species Fisher exact tests (enrichment in core)
records = []
for sp in species:
    a = core_essential[sp]
    b = core_non_essential[sp]
    c = total_essential[sp] - a
    d = genome_size[sp] - (a + b + c)          
    OR, p = fisher_exact([[a, b], [c, d]])
    records.append(dict(Species=sp, a=a, b=b, c=c, d=d, OR=OR, p_raw=p))

df_fisher = pd.DataFrame(records)
df_fisher["q_BH"] = multipletests(df_fisher["p_raw"], method="fdr_bh")[1]
df_fisher.to_csv("species_fisher.tsv", sep="\t", index=False)

# Spearman correlation (core-essential vs total essential)
ce_vec  = np.array([core_essential[s] for s in species])
tot_vec = np.array([total_essential[s] for s in species])
rho, p_rho = spearmanr(ce_vec, tot_vec)

# Console & file report

report.append("=== 2 x 6 x² test (core-essential vs species) ===")
report.append(f"x² = {chi2:.2f},  df = {dof},  p = {p_chi:.3e}\n")

report.append("=== Fisher tests (core enrichment per species) ===")
report.append(tabulate(
    df_fisher[["Species","OR","p_raw","q_BH"]], headers="keys",
    tablefmt="github", floatfmt=".3g"))
report.append("")

report.append("=== Spearman correlation (core-essential vs total essential) ===")
report.append(f"p = {rho:.3f},  p = {p_rho:.3e}")

txt = "\n".join(report)
print(txt)

with open("summary.txt", "w") as fh:
    fh.write(txt)

print("\nFiles written:  summary.txt   species_fisher.tsv")
