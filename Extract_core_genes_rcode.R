# Load necessary libraries
if (!requireNamespace("dplyr", quietly = TRUE)) {
  install.packages("dplyr")
}
library(dplyr)

# Function to read in the .RTAB file and remove any columns without a name
read_rtab <- function(file_path) {
  data <- read.table(file_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
  # Remove unnamed columns
  data <- data[, colnames(data) != ""]
  return(data)
}

# Function to filter genes present in all samples
filter_genes_present_in_all <- function(data) {
  data %>%
    filter(apply(.[, -1], 1, function(row) all(row != 0)))
}

# Function to save the first column of a data frame as a comma-separated list
save_genes_to_txt <- function(data, output_file) {
  write.table(data[, 1], file = output_file, row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
}

# Iterate over every .Rtab file in all subdirectories
process_rtab_files <- function() {
  # Find all .Rtab files in subdirectories
  file_paths <- list.files(path = ".", pattern = "\\.Rtab$", recursive = TRUE, full.names = TRUE)
  
  for (file_path in file_paths) {
    # Reading the file
    rtab_data <- read_rtab(file_path)
    
    # Filtering the data
    filtered_data <- filter_genes_present_in_all(rtab_data)
    
    # Define output file name based on input file
    output_file <- gsub("\\.Rtab$", "_core_peppan_gene_locuses.txt", file_path)
    
    # Save the first column of filtered_data
    save_genes_to_txt(filtered_data, output_file)
  }
}

# Call the function to process all .Rtab files
process_rtab_files()