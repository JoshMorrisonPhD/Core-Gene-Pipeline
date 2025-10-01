Install and load the openxlsx package if not already installed
if (!require("openxlsx")) install.packages("openxlsx")
library(openxlsx)

Get all .tsv files in the working directory
tsv_files <- list.files(pattern = "\.tsv$")

Loop through each file
for (file in tsv_files) {
  # Read the .tsv file
  data <- read.delim(file, sep = "\t")
  
  # Create the new file name with .xlsx extension
  xlsx_file <- sub("\.tsv$", ".xlsx", file)
  
  # Write the data to an .xlsx file
  write.xlsx(data, file = xlsx_file)
}