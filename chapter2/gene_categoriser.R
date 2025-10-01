# Load necessary packages
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)

# Load the data
data <- read.table("PEPPAN.PEPPAN.gene_content.Rtab", header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

# Count the number of genomes (excluding the first column, which is gene name)
num_genomes <- ncol(data) - 1

# Convert all values greater than 1 to 1, excluding the first column
data[,-1] <- lapply(data[,-1], function(x) ifelse(x > 1, 1, x))

# Calculate the presence percentage for each gene
data <- data %>%
  rowwise() %>%
  mutate(presence_percentage = sum(c_across(-Gene)) / num_genomes * 100)

# Categorize genes based on their presence percentage
data <- data %>%
  mutate(category = case_when(
    presence_percentage == 100 ~ "Strict core gene",
    presence_percentage >= 99 & presence_percentage < 100 ~ "Core gene",
    presence_percentage >= 95 & presence_percentage < 99 ~ "Soft core gene",
    presence_percentage >= 15 & presence_percentage < 95 ~ "Shell gene",
    presence_percentage < 15 ~ "Cloud gene"
  ))

# Select only relevant columns: Gene, presence_percentage, and category
final_data <- data %>% select(Gene, presence_percentage, category)

# Save the categorized data to a single file
write.table(final_data, "categorized_genes_with_names.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Create separate files for each category
categories <- c("Strict core gene", "Core gene", "Soft core gene", "Shell gene", "Cloud gene")
file_names <- c("strict_core_genes.txt", "core_genes.txt", "soft_core_genes.txt", "shell_genes.txt", "cloud_genes.txt")

for (i in seq_along(categories)) {
  # Extract only genes that belong to the current category
  genes <- final_data %>% filter(category == categories[i]) %>% select(Gene)
  
  # Save genes as a line-separated file
  write.table(genes, file_names[i], sep = "\n", row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# Print message indicating completion
print("Categorization completed. Files saved:")
print(file_names)
