#The input for this is the same .xlsx files as for UG_plotting, it needs to be .xlsx or it won't work.
# The output of this file is plots based on the stats it runs.
# this script runs a chi-square test then plots correlation matrices of the chi-square residuals, then plots a bar chart of the COG category distribution.

# Load necessary libraries
if (!require("ggplot2")) install.packages("ggplot2")
if (!require("dplyr")) install.packages("dplyr")
if (!require("tidyr")) install.packages("tidyr")
if (!require("openxlsx")) install.packages("openxlsx")
if (!require("stats")) install.packages("stats")
if (!require("corrplot")) install.packages("corrplot")

library(ggplot2)
library(dplyr)
library(tidyr)
library(openxlsx)
library(stats)
library(corrplot)

# Define color codes
colors <- c(
  "Core"  = "#532C6B",   # deep purple
  "Shell" = "#6BBF59",   # green
  "Cloud" = "#F1D22B"    # yellow
)

# Function to plot the data and show residuals
plot_and_show_residuals <- function(species_data, species_name) {
  # Filter out NA values
  species_data <- species_data %>%
    filter(!is.na(NEWCAT_COUNT))
  
  # Define the order of gene groups
  species_data$GENECAT <- factor(species_data$GENECAT, levels = c("Cloud", "Shell", "Core"))
  
  # Create a contingency table using xtabs
  table_data <- xtabs(NEWCAT_COUNT ~ GENECAT + New_Category, data = species_data)
  
  # Check if we have enough data for chi-square test
  if (dim(table_data)[1] > 1 && dim(table_data)[2] > 1) {
    # Perform chi-square test
    chi_sq_test <- chisq.test(table_data)
    chi_sq_p <- chi_sq_test$p.value
    
    # Get residuals and handle special values
    chisq_residuals <- chi_sq_test$residuals
    chisq_residuals[is.na(chisq_residuals)] <- 0
    chisq_residuals[is.infinite(chisq_residuals)] <- 0
    
    # Rename specific category
    rownames(chisq_residuals) <- gsub("DNA/RNA Storage, Translation, and Transcription",
                                      "DNA/RNA Storage & Use",
                                      rownames(chisq_residuals))
    colnames(chisq_residuals) <- gsub("DNA/RNA Storage, Translation, and Transcription",
                                      "DNA/RNA Storage & Use",
                                      colnames(chisq_residuals))
    
    # Save chi-square test results to text file
    sink(paste0(species_name, "_ChiSquare_Grouped.txt"))
    cat("===== Chi-Square Test Results for", species_name, "=====\n")
    print(chi_sq_test)
    cat("\n===== Chi-Square Residuals for", species_name, "=====\n")
    print(round(chisq_residuals, 3))
    sink()
    
    # Print Chi-Square residuals to console
    cat("\n===== Chi-Square Residuals for", species_name, "=====\n")
    print(round(chisq_residuals, 3))
    
    # Plot residuals heatmap
    print("Chi-Square Residuals Heatmap:")
    corrplot(chisq_residuals,
             is.cor = FALSE,
             tl.srt = 45,       # Rotate labels
             tl.cex = 1,        # Font size
             tl.col = "black",
             cl.cex = 1,
             title = "")  # Removed title
    
    # Save residuals plot
    png(paste0(species_name, "_Group_Corr.png"), width = 1000, height = 800)
    corrplot(chisq_residuals,
             is.cor = FALSE,
             tl.srt = 45,
             tl.cex = 1,
             tl.col = "black",
             cl.cex = 1,
             title = "")
    dev.off()
  } else {
    chi_sq_p <- NA
    print(paste("Not enough categories for chi-square test in", species_name))
    
    # Still create a text file with the error message
    sink(paste0(species_name, "_ChiSquare_Grouped.txt"))
    cat("Not enough categories for chi-square test in", species_name, "\n")
    sink()
  }
  
  # Create unstacked bar chart
  p <- ggplot(species_data, aes(x = New_Category, y = NEWCAT_COUNT, fill = GENECAT)) +
    geom_col(position = position_dodge(width = 0.7), width = 0.7) +
    scale_fill_manual(values = colors) +
    labs(title = bquote(paste("Gene category distribution proportions by COG supercategory in ", italic(paste("S.", .(species_name))))),
         x = "COG Supercategory",
         y = "Count",
         fill = "Gene Category") +
    theme_minimal() +
    theme(legend.position = "right",
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
          plot.title = element_text(hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  
  # Add chi-square p-value annotation
  if (!is.na(chi_sq_p)) {
    p <- p + annotate("text", x = Inf, y = Inf,
                      label = paste("Chi-square p-value =", format(chi_sq_p, digits = 4)),
                      hjust = 1.1, vjust = 1.1, size = 3.5)
  }
  
  # Save bar chart plot
  ggsave(paste0(species_name, "_Grouped_Barchart.png"), plot = p, width = 10, height = 8)
  
  return(p)
}

# Get all .xlsx files in the working directory
xlsx_files <- list.files(pattern = "\\.xlsx$")

# Loop through each file
for (file in xlsx_files) {
  # Extract species name from file name (remove _total_classifier_reorg.xlsx)
  species_name <- gsub("_total_classifier_reorg.xlsx", "", file)
  
  # Read data
  species_data <- read.xlsx(file)
  
  print(paste("\n\nProcessing file:", file))
  
  # Plot and show residuals
  plot <- plot_and_show_residuals(species_data, species_name)
  
  # Print the plot
  print(plot)
}
