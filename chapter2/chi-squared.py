#This script runs a chi-square test on the gene counts per category, for each species

import numpy as np
from scipy.stats import chi2_contingency

# Define the observed data (gene counts for each species)
# Rows represent species, columns represent gene categories
data = np.array([
    [1556, 96, 691, 2957],  # S. agalactiae
    [1786, 20, 315, 226],    # S. equi
    [1561, 77, 258, 1151],   # S. iniae
    [1462, 186, 754, 1939],  # S. pneumoniae
    [1279, 132, 1070, 4164], # S. suis
    [1423, 57, 599, 4255]    # S. uberis
])

# Perform Chi-Square test testing H0 = all gene category proportions are independant within each species. If P<=0.05 then H0 is false.
chi2_stat, p_value, dof, expected = chi2_contingency(data)

# Save results to a text file
with open("chi_square_results.txt", "w") as file:
    file.write("Chi-Square Test for Gene Category Proportions Across Species\n")
    file.write("=" * 60 + "\n")
    file.write(f"Chi-Square Statistic: {chi2_stat:.4f}\n")
    file.write(f"P-value: {p_value:.4f}\n")
    file.write(f"Degrees of Freedom: {dof}\n")
    file.write("\nExpected Frequencies (If Proportions Were the Same):\n")
    file.write(str(expected) + "\n")
    file.write("=" * 60 + "\n")

print("Chi-Square test results have been saved to 'chi_square_results.txt'.")
