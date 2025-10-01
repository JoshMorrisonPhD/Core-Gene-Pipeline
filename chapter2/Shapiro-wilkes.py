# this script runs a Shapiro-wilk test on the gene category data

from scipy.stats import shapiro

# Define the data
data = {
    "Strict Core Genes": [1556, 1786, 1561, 1462, 1279, 1423],
    "Soft Core Genes": [96, 20, 77, 186, 132, 57],
    "Shell Genes": [691, 315, 258, 754, 1070, 599],
    "Cloud Genes": [2957, 226, 1151, 1939, 4164, 4255],
    "Total Genes": [5300, 2347, 3047, 4341, 6645, 6334]
}

# Open a text file to save the results
with open("normality_results.txt", "w") as file:
    file.write("Shapiro-Wilk Normality Test Results\n")
    file.write("=" * 40 + "\n")

    # Perform Shapiro-Wilk test for each category
    for column, values in data.items():
        stat, p_value = shapiro(values)
        conclusion = "Normally Distributed" if p_value > 0.05 else "Not Normally Distributed"

        # Write results to the file
        file.write(f"Category: {column}\n")
        file.write(f"  Shapiro-Wilk Test Statistic: {stat:.4f}\n")
        file.write(f"  P-value: {p_value:.4f}\n")
        file.write(f"  Conclusion: {conclusion}\n")
        file.write("-" * 40 + "\n")

print("Shapiro-Wilk test results have been saved to 'normality_results.txt'.")
