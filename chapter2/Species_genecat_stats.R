# Load necessary libraries
if (!requireNamespace("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")

library(ggplot2)
library(dplyr)

# Define gene counts for each category across species
categories <- list(
  "Strict Core Genes" = c(1556, 1786, 1561, 1462, 1279, 1423),
  "Soft Core Genes" = c(44, 20, 77, 92, 132, 57) + c(52, 0, 0, 94, 0, 0),  # Merging Core Genes into Soft Core
  "Shell Genes" = c(691, 315, 258, 754, 1070, 599),
  "Cloud Genes" = c(2957, 226, 1151, 1939, 4164, 4255)
)

# Convert data into a data frame for visualization
data_list <- lapply(names(categories), function(category) {
  data.frame(
    Category = category,
    Species = factor(1:length(categories[[category]])),
    Value = categories[[category]]
  )
})

df <- bind_rows(data_list)

# Initialize an empty results dataframe
anova_results <- data.frame(
  Category = character(),
  ANOVA_Statistic = numeric(),
  P_Value = numeric(),
  Conclusion = character(),
  stringsAsFactors = FALSE
)

# Perform ANOVA for each gene category
for (category in names(categories)) {
  values <- categories[[category]]
  species <- factor(rep(1:length(values), each = 1))  # Factor representing species groups
  
  # Check if there is variance in the data
  if (length(unique(values)) > 1) {
    anova_model <- aov(values ~ species)  # Perform one-way ANOVA
    anova_summary <- summary(anova_model)
    
    anova_stat <- anova_summary[[1]]$`F value`[1]
    anova_p_value <- anova_summary[[1]]$`Pr(>F)`[1]
    
    conclusion <- ifelse(anova_p_value < 0.05, "Significantly Different", "Not Significantly Different")
  } else {
    anova_stat <- NA
    anova_p_value <- NA
    conclusion <- "No variation in data"
  }
  
  # Append results to dataframe
  anova_results <- rbind(
    anova_results, 
    data.frame(
      Category = category,
      ANOVA_Statistic = anova_stat,
      P_Value = anova_p_value,
      Conclusion = conclusion,
      stringsAsFactors = FALSE
    )
  )
}

# Save ANOVA results to a text file
write.table(anova_results, "anova_all_categories.txt", row.names = FALSE, sep = "\t", quote = FALSE)

# Boxplot of gene counts across species
p1 <- ggplot(df, aes(x = Species, y = Value, fill = Category)) +
  geom_boxplot() +
  labs(title = "Gene Count Distribution Across Species", x = "Species", y = "Gene Count") +
  theme_minimal()

ggsave("boxplot_gene_counts.png", plot = p1, width = 8, height = 6)

# Bar plot of ANOVA statistics
p2 <- ggplot(anova_results, aes(x = Category, y = ANOVA_Statistic, fill = Conclusion)) +
  geom_bar(stat = "identity", na.rm = TRUE) +  # Ignore NA values
  labs(title = "ANOVA Statistics for Each Gene Category", x = "Gene Category", y = "ANOVA Statistic") +
  theme_minimal() +
  coord_flip()

ggsave("anova_statistics_barplot.png", plot = p2, width = 8, height = 6)

# P-value plot
p3 <- ggplot(anova_results, aes(x = Category, y = P_Value, fill = Conclusion)) +
  geom_bar(stat = "identity", na.rm = TRUE) +  # Ignore NA values
  labs(title = "ANOVA P-values for Each Gene Category", x = "Gene Category", y = "P-Value") +
  theme_minimal() +
  coord_flip()

ggsave("anova_pvalues_barplot.png", plot = p3, width = 8, height = 6)

print("ANOVA test results and plots have been saved!")
