#This uses the same .tsv files as COGSR_ungrouped.R as inputs, as well as the .txt results from COGSR_ungrouped.R. If the ,tsv files dont work, run tsv_to_xlsx.R
#The ouputs are color coded stacked bar charts with KW Pvalues for each species in the PWD


# Load necessary libraries
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("readr")) install.packages("readr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("openxlsx")) install.packages("openxlsx")
if (!require("stringr")) install.packages("stringr")

library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)
library(openxlsx)
library(stringr)

# Define color codes
colors <- c(
  "Core"  = "#532C6B",   # deep purple
  "Shell" = "#6BBF59",   # green
  "Cloud" = "#F1D22B"    # yellow
)

# extract p-values from the results file
extract_p_values <- function(file_path) {
  if (!file.exists(file_path)) {
    return(NA)
  }
  
  lines <- readLines(file_path)
  
  # Extract Kruskal-Wallis p-value
  kruskal_line <- grep("P-value", lines, ignore.case = TRUE, value = TRUE)
  
  # Check if a line was matched
  if (length(kruskal_line) == 0) {
    return(NA)
  }
  
  # interate over Pvalue patterns (reads the .txt files and finds the first match)
  p_value_patterns <- c("P-value\\s+([0-9.]+)", "P.value\\s+([0-9.]+)", "p.value\\s+([0-9.]+)",
                        "P-value\\s*[:=]\\s*([0-9.]+)", "p.value\\s*[:=]\\s*([0-9.]+)")
  
  for (pattern in p_value_patterns) {
    match <- regmatches(kruskal_line, regexpr(pattern, kruskal_line, ignore.case = TRUE))
    if (length(match) > 0 && !is.na(match)) {
      extracted_value <- str_extract(kruskal_line, pattern)
      if (!is.na(extracted_value)) {
        # Extract the numeric part
        p_value <- as.numeric(str_extract(extracted_value, "[0-9.]+"))
        return(p_value)
      }
    }
  }
  
  return(NA)
}

# Function to plot the data for a species
plot_species_data <- function(species_data, species_name, kruskal_p) {
  # Define the order of gene groups
  species_data$GENECAT <- factor(species_data$GENECAT, levels = c("Cloud", "Shell", "Core"))
  
  p <- ggplot(species_data, aes(x = DESCRIPTION, y = COUNT, fill = GENECAT)) +
    geom_col(position = "stack") +
    scale_fill_manual(values = colors) +
    labs(title = bquote(paste("Stacked Bar Chart of ", italic(paste("S.", .(species_name))), " Gene Group Distribution Within COG Categories")),
         x = "COG Categories",
         y = "Count",
         fill = "Gene Group") +
    theme_minimal() +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(hjust = 0.5))
  
  # Add p-value annotation only if we have a valid value
  if (!is.na(kruskal_p)) {
    p <- p + annotate("text", x = Inf, y = Inf, label = paste("Kruskal-Wallis p-value =", format(kruskal_p, digits = 4)),
                      hjust = 1.1, vjust = 1.1, size = 3.5)
  }
  
  return(p)
}

# Get all .xlsx files in the working directory
xlsx_files <- list.files(pattern = "\\.xlsx$")

# Get all .txt files in the working directory
txt_files <- list.files(pattern = "\\.txt$")

# Loop through each file
for (file in xlsx_files) {
  # Extract species name from file name (remove _total_classifier_reorg.xlsx)
  species_name <- gsub("_total_classifier_reorg.xlsx", "", file)
  
  # Define result file path by replacing .xlsx with _results.txt
  results_file <- gsub("\\.xlsx$", "_results.txt", file)
  
  # Read data
  species_data <- read.xlsx(file)
  
  # Remove rows where DESCRIPTION is related to "no COG function" or similar
  species_data <- species_data %>%
    filter(!grepl("Function unknown|not in COGs|no COG", DESCRIPTION, ignore.case = TRUE))
  
  # Check if the corresponding results file exists
  if (results_file %in% txt_files) {
    kruskal_p <- extract_p_values(results_file)
  } else {
    kruskal_p <- NA
    print(paste("Warning: Results file not found for", file))
  }
  
  # Print file paths for debugging
  print(paste("Excel file:", file))
  print(paste("Looking for results file:", results_file))
  
  # Plot
  plot <- plot_species_data(species_data, species_name, kruskal_p)
  
  # Print the plot
  print(plot)
}
