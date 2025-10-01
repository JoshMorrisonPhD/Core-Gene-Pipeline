# the inputs for this are the number of genes found in the described overlap (hardcoded) obtained by running then desired overlap within the pangenome.R script
# the outputs are a stacked bar plot of the genes in the described overlap and statistics for the described overlap in a .csv file (chi-squared test, fishers test)
# this script was used to generate the plots and stats for the zoonosis risk genes (pneumoniae vs all non-pneumoniae genes)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


# INPUT 

categories = ["Strict core", "Collapsed core", "Shell", "Cloud"]
total   = np.array([1462, 186, 754, 1939])
matched = np.array([356,   9,   155, 412])


# Prep

matched   = np.minimum(matched, total)
unmatched = np.maximum(total - matched, 0)

# Normalized proportions (%)
matched_pct   = matched / total * 100.0
unmatched_pct = unmatched / total * 100.0

# Colors
col_strict_core    = "#532C6B"  # purple
col_collapsed_core = "#377EB8"  # blue
col_shell          = "#2CA25F"  # green
col_cloud          = "#F2CC45"  # yellow
unmatched_color    = "#7F7F7F"  # grey

cat_color = {
    "Strict core": col_strict_core,
    "Collapsed core": col_collapsed_core,
    "Shell": col_shell,
    "Cloud": col_cloud,
}
matched_colors = [cat_color[c] for c in categories]


# Plot (normalized 0–100%)

x = np.arange(len(categories))
width = 0.6

fig, ax = plt.subplots(figsize=(9, 6))

# Bars (no labels; we build a custom legend)
bars_matched = ax.bar(
    x, matched_pct, width,
    color=matched_colors, edgecolor="black", linewidth=0.5
)
bars_unmatched = ax.bar(
    x, unmatched_pct, width, bottom=matched_pct,
    color=unmatched_color, edgecolor="black", linewidth=0.5
)

# Inside annotations (counts + %)
for i, bar in enumerate(bars_matched):
    h = bar.get_height()
    if h > 0:
        ax.annotate(f"{matched[i]} ({matched_pct[i]:.1f}%)",
                    xy=(bar.get_x() + bar.get_width()/2, bar.get_y() + h/2),
                    ha="center", va="center", fontsize=9, color="white", fontweight="bold")
for i, bar in enumerate(bars_unmatched):
    h = bar.get_height()
    if h > 0:
        ax.annotate(f"{unmatched[i]} ({unmatched_pct[i]:.1f}%)",
                    xy=(bar.get_x() + bar.get_width()/2, matched_pct[i] + h/2),
                    ha="center", va="center", fontsize=9, color="white", fontweight="bold")

# Totals above bars
for i in range(len(x)):
    ax.annotate(f"Total: {int(total[i])}",
                xy=(x[i], 102), ha="center", va="bottom",
                fontsize=10, fontweight="bold")

# Axes/labels
ax.set_ylabel("Percentage of genes (%)")
ax.set_xlabel("Category")
ax.set_title("Genes in S. pneumoniae with orthologs in zoonotic species")
ax.set_xticks(x)
ax.set_xticklabels(categories, rotation=10)

# Custom legend
legend_handles = [
    Patch(facecolor=col_strict_core,    edgecolor="black", label="Matched – Strict core"),
    Patch(facecolor=col_collapsed_core, edgecolor="black", label="Matched – Collapsed core"),
    Patch(facecolor=col_shell,          edgecolor="black", label="Matched – Shell"),
    Patch(facecolor=col_cloud,          edgecolor="black", label="Matched – Cloud"),
    Patch(facecolor=unmatched_color,    edgecolor="black", label="Unmatched"),
]
ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False)

ax.set_ylim(0, 115)
ax.grid(axis="y", linestyle=":", alpha=0.4)
ax.set_axisbelow(True)
plt.tight_layout()


# Stats (with safe fallbacks)
p_global = None
pairwise_df = pd.DataFrame()
try:
    from scipy.stats import chi2_contingency, fisher_exact
    from itertools import combinations

    contingency = np.vstack([matched, unmatched])  # 2 x 4
    chi2, p_global, dof, expected = chi2_contingency(contingency)

    rows = []
    for (i, j) in combinations(range(len(categories)), 2):
        table_2x2 = np.array([[matched[i],   unmatched[i]],
                              [matched[j],   unmatched[j]]])
        odds, p = fisher_exact(table_2x2, alternative="two-sided")
        rows.append({
            "cat_A": categories[i], "cat_B": categories[j],
            "matched_A": int(matched[i]), "unmatched_A": int(unmatched[i]),
            "matched_B": int(matched[j]), "unmatched_B": int(unmatched[j]),
            "odds_ratio": odds, "p_value": p
        })
    pairwise_df = pd.DataFrame(rows)

    # FDR correction if statsmodels is available
    try:
        from statsmodels.stats.multitest import multipletests
        rej, p_adj, _, _ = multipletests(pairwise_df["p_value"].values, method="fdr_bh")
        pairwise_df["p_adj_fdr"] = p_adj
        pairwise_df["reject_fdr_0.05"] = rej
    except Exception:
        pairwise_df["p_adj_fdr"] = np.nan
        pairwise_df["reject_fdr_0.05"] = np.nan


except Exception as e:
    # SciPy not present or version mismatch
    print("[warn] Stats skipped:", e)
    ax.text(0.01, 0.98,
            "Chi-squared: not computed (SciPy unavailable)",
            ha="left", va="top", transform=ax.transAxes,
            fontsize=10, bbox=dict(boxstyle="round,pad=0.3", fc="white", ec="gray", alpha=0.85))


# Save outputs

plt.savefig("pneumo_genes_stacked_normalized.png", dpi=600)
plt.savefig("pneumo_genes_stacked_normalized.svg")
plt.show()

# Export counts + percents
pd.DataFrame({
    "category": categories,
    "total": total.astype(int),
    "matched": matched.astype(int),
    "unmatched": unmatched.astype(int),
    "matched_pct": matched_pct,
    "unmatched_pct": unmatched_pct
}).to_csv("pneumo_genes_counts_and_pcts.csv", index=False)

if not pairwise_df.empty:
    pairwise_df.to_csv("pneumo_genes_pairwise_fishers.csv", index=False)
    print(f"Global chi-squared p-value: {p_global:.3e}")
    print("Pairwise Fisher’s (FDR where available):")
    print(pairwise_df)
else:
    print("Pairwise tests not computed.")
