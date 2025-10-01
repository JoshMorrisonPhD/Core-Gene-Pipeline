# This code's input is all the .tsv files in the current working directory. These .TSV files were created by appending each category of genes' COG results  file classifier_count.tsv to eachother, 
# resulting in a  single .tsv file. A new column was added to this in excel called GENECAT to assign gene category to that result row. The NEWCAT columns are used in another script.
# The output of this code is a.txt file containing the results of the statistical analysis ($species_total_classifier_reorg_results.txt)
# The purpose of this script is to run statistical tests on COGS data, first a Shapiro_wilk normality test, then a Kruskal-Wallis test and finally a Dunn's post-hoc test.

# Load necessary libraries
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("FSA")) install.packages("FSA")
if (!require("ggpubr")) install.packages("ggpubr")
if (!require("readr")) install.packages("readr")
if (!require("knitr")) install.packages("knitr")

library(tidyr)
library(dplyr)
library(FSA)      # For Dunn's test
library(ggpubr)   # For Kruskal-Wallis test
library(readr)    # For properly reading CSV/TSV files
library(knitr)    # For pretty table output

# List all .tsv files in the current directory
file_list <- list.files(path=".", pattern="\\.tsv$", full.names=TRUE, recursive=FALSE)

# Process each file
for (file in file_list) {
  # Extract species name from filename (removing path and extension)
  species_name <- tools::file_path_sans_ext(basename(file))
  
  # Load dataset
  data <- read.delim(file, sep="\t", header=TRUE) %>%
    select(c(1, 2, 3)) %>%  # Remove columns 4 and 5
    drop_na()                # Remove NA values
  
  # Ensure data has required columns
  if (!all(c("GENECAT", "LETTER", "COUNT") %in% colnames(data))) {
    message(paste("Skipping", file, "- Missing required columns"))
    next
  }
  
  # Ensure categorical variables are factors
  data$GENECAT <- as.factor(data$GENECAT)
  data$LETTER<- as.factor(data$LETTER)
  
  # Check if dataset is empty
  if (nrow(data) == 0) {
    message(paste("Skipping", file, "- No valid data after filtering"))
    next
  }
  
  # Output file paths
  output_text <- paste0(species_name, "_results.txt")
  
  # Open a file connection to save results
  sink(output_text)  # Redirect output to a text file
  
  # Shapiro-Wilk normality test (per GENECAT)
  test_results <- data %>%
    group_by(GENECAT) %>%
    summarise(shapiro_p = ifelse(n() > 2, shapiro.test(COUNT)$p.value, NA_real_)) %>%
    ungroup()
  
  # Print Shapiro-Wilk test results
  cat("\n===== Shapiro-Wilk Test Results for", species_name, "=====\n")
  print(test_results)
  
  # Kruskal-Wallis Test
  kw_test <- kruskal.test(COUNT ~ GENECAT, data = data)
  
  # Print Kruskal-Wallis result table
  cat("\n===== Kruskal-Wallis Test Results for", species_name, "=====\n")
  print(kw_test)
  
  # Run Dunn's test if Kruskal-Wallis is significant
   posthoc_results <- NULL
  if (!is.null(kw_test$p.value) && kw_test$p.value < 0.05) {
   dunn_test <- dunnTest(COUNT ~ GENECAT, data = data, method = "bh")
   posthoc_results <- dunn_test$res
  
   # Print Dunn's post-hoc test results
   cat("\n===== Dunn's Post-Hoc Test Results for", file, "=====\n")
    print(kable(posthoc_results, format = "simple"))
  } else {
     cat("\nKruskal-Wallis test was not significant. Dunn's test was not performed.\n")
}

  # Close the text output file
  sink()
}