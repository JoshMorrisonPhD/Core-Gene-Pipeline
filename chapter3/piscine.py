# the inputs for this are the number of genes found in the described overlap (hardcoded) obtained by running then desired overlap within the pangenome.R script
# the outputs are a stacked bar plot of the genes in the described overlap and statistics for the described overlap in a .csv file (chi-squared test, fishers test)
# this script was used to generate the plots and stats for the piscine genes (agalactiae + iniae)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Patch


# TOGGLE: True → normalize bars to 100%; False → absolute counts

NORMALIZE_TO_PERCENT = True


# INPUT

categories = ["Strict core", "Collapsed core", "Shell", "Cloud"]

# PLACEHOLDER numbers — put yours here
total   = np.array([1137,  0,  50, 393], dtype=float)
matched = np.array([ 719,   0,  0,  7], dtype=float)

# Safety
matched   = np.minimum(matched, total)
unmatched = np.maximum(total - matched, 0) 

# Colors (category colors now applied to UNMATCHED)
col_strict_core    = "#532C6B"  # purple
col_collapsed_core = "#377EB8"  # blue
col_shell          = "#2CA25F"  # green
col_cloud          = "#F2CC45"  # yellow
matched_grey       = "#7F7F7F"  # grey for the superpangenome-overlap part

cat_color = {
    "Strict core": col_strict_core,
    "Collapsed core": col_collapsed_core,
    "Shell": col_shell,
    "Cloud": col_cloud,
}
unmatched_colors = [cat_color[c] for c in categories]

# Data to plot (counts or percents)
if NORMALIZE_TO_PERCENT:
    matched_plot   = matched / total * 100.0
    unmatched_plot = unmatched / total * 100.0
    y_label = "Percentage of genes (%)"
    y_top   = 115  # headroom
else:
    matched_plot   = matched
    unmatched_plot = unmatched
    y_label = "Number of genes"
    y_top   = (matched + unmatched).max() * 1.18


# Plot

x = np.arange(len(categories))
width = 0.6
fig, ax = plt.subplots(figsize=(9, 6))

# Bottom (grey) = Matched (also in superpangenome core)
bars_matched = ax.bar(
    x, matched_plot, width,
    color=matched_grey, edgecolor="black", linewidth=0.5
)
# Top (colored) = Unmatched (Piscine-specific core)
bars_unmatched = ax.bar(
    x, unmatched_plot, width, bottom=matched_plot,
    color=unmatched_colors, edgecolor="black", linewidth=0.5
)

# Annotations
# Inside grey (matched)
for i, bar in enumerate(bars_matched):
    h = bar.get_height()
    if h > 0:
        if NORMALIZE_TO_PERCENT:
            label = f"{int(matched[i])} ({matched_plot[i]:.1f}%)"
            y_pos = bar.get_y() + h/2
        else:
            pct = matched[i] / total[i] * 100 if total[i] else 0
            label = f"{int(matched[i])} ({pct:.1f}%)"
            y_pos = bar.get_y() + h/2
        ax.annotate(label,
                    xy=(bar.get_x() + bar.get_width()/2, y_pos),
                    ha="center", va="center", fontsize=9, color="white", fontweight="bold")

# Inside colored (unmatched = Piscine-specific)
for i, bar in enumerate(bars_unmatched):
    h = bar.get_height()
    if h > 0:
        if NORMALIZE_TO_PERCENT:
            label = f"{int(unmatched[i])} ({unmatched_plot[i]:.1f}%)"
            y_pos = matched_plot[i] + h/2
        else:
            pct = unmatched[i] / total[i] * 100 if total[i] else 0
            label = f"{int(unmatched[i])} ({pct:.1f}%)"
            y_pos = matched_plot[i] + h/2
        ax.annotate(label,
                    xy=(bar.get_x() + bar.get_width()/2, y_pos),
                    ha="center", va="center", fontsize=9, color="black", fontweight="bold")

# Totals above bars
for i in range(len(x)):
    if NORMALIZE_TO_PERCENT:
        ax.annotate(f"Total: {int(total[i])}",
                    xy=(x[i], 102), ha="center", va="bottom",
                    fontsize=10, fontweight="bold")
    else:
        ax.annotate(f"Total: {int(total[i])}",
                    xy=(x[i], matched_plot[i] + unmatched_plot[i] + y_top*0.02),
                    ha="center", va="bottom", fontsize=10, fontweight="bold")

# Labels/axes
ax.set_ylabel(y_label)
ax.set_xlabel("Category")
ax.set_title("Piscine-specific portion of the super-pangnome")
ax.set_xticks(x)
ax.set_xticklabels(categories, rotation=10)

# Custom legend (colored = UNMATCHED/Piscine-specific; grey = matched/super-core)
legend_handles = [
    Patch(facecolor=col_strict_core,    edgecolor="black", label="Piscine Strict core"),
    Patch(facecolor=col_collapsed_core, edgecolor="black", label="Piscine Collapsed core"),
    Patch(facecolor=col_shell,          edgecolor="black", label="Piscine Shell"),
    Patch(facecolor=col_cloud,          edgecolor="black", label="Piscine Cloud"),
    Patch(facecolor=matched_grey,       edgecolor="black", label="Super-pangenome"),
]
ax.legend(handles=legend_handles, loc="upper left", bbox_to_anchor=(1.02, 1), frameon=False)

ax.set_ylim(0, y_top)
ax.grid(axis="y", linestyle=":", alpha=0.4)
ax.set_axisbelow(True)
plt.tight_layout()


plt.tight_layout()
plt.savefig("piscine_by_category.png", dpi=600)
plt.savefig("piscine_by_category.svg")
plt.show()

# (Optional) also save a CSV of the counts/percents used
pd.DataFrame({
    "category": categories,
    "total": total.astype(int),
    "matched": matched.astype(int),
    "unmatched": unmatched.astype(int),
    "matched_pct": (matched/total*100).round(2),
    "unmatched_pct": (unmatched/total*100).round(2),
}).to_csv("piscine_by_category.csv", index=False)
