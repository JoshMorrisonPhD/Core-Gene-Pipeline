# The input for this is the .Rtab file for all species combined.
# The ouput of this is the supertest results, broken down into categories. SUpplementary outputs can be created based on user choices within this code (refer to lines 89 to 95 and 250 to 268 for an example)
# This code also creates plots for all of the chosen intersections (both linear and circular)


#Load packages

if (!require("dplyr")) install.packages("dplyr")
if (!require("SuperExactTest")) install.packages("SuperExactTest")

library(dplyr)
library(SuperExactTest)

# Load data
data <- read.table("agal_pneumo_iniae_uberis_equi_suis.PEPPAN.gene_content.Rtab",
                   header = TRUE, sep = "\t", check.names = FALSE)

# Convert counts > 1 into presence/absence (1/0)
data[,-1] <- lapply(data[,-1], function(x) ifelse(x > 1, 1, x))

# Function to categorize genes for one species
categorize_species_genes <- function(data, species_prefix) {
  # Extract prefixes from column names (everything before first "_")
  toks <- tolower(sub("_.*$", "", colnames(data)[-1]))
  
  # Get columns belonging to this species
  species_cols <- colnames(data)[-1][toks == species_prefix]
  
  if (length(species_cols) == 0) {
    message("No columns found for species: ", species_prefix)
    return(NULL)
  }
  
  # Subset data for this species
  sub <- data[, c("Gene", species_cols), drop = FALSE]
  num_genomes <- ncol(sub) - 1
  
  # Calculate presence percentage across genomes of this species
  sub <- sub %>%
    rowwise() %>%
    mutate(presence_percentage = sum(c_across(-Gene)) / num_genomes * 100) %>%
    ungroup()
  
  # Categorize into strict core / collapsed_core / shell / cloud
  sub <- sub %>%
    mutate(category = case_when(
      presence_percentage == 100 ~ "Strict core gene",
      presence_percentage >= 95 & presence_percentage < 100 ~ "collapsed_core",
      presence_percentage >= 15 & presence_percentage < 95 ~ "Shell gene",
      presence_percentage > 0 & presence_percentage < 15 ~ "Cloud gene"
    ))
  
  # Return a named list of vectors (category → gene names)
  split(sub$Gene, sub$category)
}

# Apply to all species 
species_order <- c("agal", "pneumo", "iniae", "uberis", "equi", "suis")

species_gene_categories <- lapply(species_order, function(sp) {
  categorize_species_genes(data, sp)
})
names(species_gene_categories) <- species_order

#  Reorganize: categories → species 
all_categories <- c("Strict core gene", "collapsed_core", "Shell gene", "Cloud gene")

category_sets <- lapply(all_categories, function(cat) {
  sapply(species_gene_categories, function(sp_list) {
    if (!is.null(sp_list[[cat]])) sp_list[[cat]] else character(0)
  }, simplify = FALSE)
})
names(category_sets) <- all_categories

# Map short species codes to readable names 
species_labels <- c(
  agal   = "S. agalactiae",
  pneumo = "S. pneumoniae",
  iniae  = "S. iniae",
  uberis = "S. uberis",
  equi   = "S. equi",
  suis   = "S. suis"
)

#  Rename category_sets with full names
category_sets <- lapply(category_sets, function(cat_list) {
  setNames(cat_list, species_labels[names(cat_list)])
})


# Helper: intersection between a subset (group) and a single species
# group_species: vector of short codes or readable names
# single_species: short code or readable name
# mode: "all" (intersection across the group) or "any" (union across the group)
# returns a list with:
#   $per_category: named list of character vectors (genes in the intersection)
#   $supertests:   named list of SuperExactTest results (Group vs Single) or NULL if empty
subset_vs_single_intersections <- function(category_sets,
                                           group_species,
                                           single_species,
                                           mode = c("all", "any"),
                                           species_labels = NULL,
                                           out_prefix = "group_vs_single") {
  mode <- match.arg(mode)
  
  # Map short codes -> readable names if needed
  map_to_readable <- function(x) {
    if (!is.null(species_labels) && !is.null(species_labels[[x]])) {
      species_labels[[x]]
    } else {
      x
    }
  }
  
  group_readable  <- unique(vapply(group_species,  map_to_readable, character(1)))
  single_readable <- map_to_readable(single_species)
  
  # Ensure species exist in the category sets
  species_in_sets <- names(category_sets[[1]])
  group_readable  <- intersect(group_readable,  species_in_sets)
  single_readable <- intersect(single_readable, species_in_sets)
  
  if (length(group_readable) == 0 || length(single_readable) == 0) {
    stop("No valid species after resolving names. Group: ",
         paste(group_readable, collapse=", "),
         " | Single: ", paste(single_readable, collapse=", "))
  }
  
  # Build intersections per category
  per_category <- setNames(vector("list", length(category_sets)), names(category_sets))
  supertests   <- setNames(vector("list", length(category_sets)), names(category_sets))
  
  # Output dirs
  sum_dir  <- file.path("supertest_summaries", paste0(out_prefix, "_summaries"))
  plot_dir <- file.path("supertest_plots",     paste0(out_prefix, "_plots"))
  dir.create(sum_dir,  showWarnings = FALSE, recursive = TRUE)
  dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)
  
  for (cat in names(category_sets)) {
    cat_sets <- category_sets[[cat]]
    
    # Collect sets for the group and the single species
    group_sets  <- cat_sets[group_readable]
    single_set  <- cat_sets[[single_readable]]
    
    # Collapse the group
    group_collapsed <- if (mode == "all") {
      # genes present in ALL species in the group
      if (length(group_sets) == 1) {
        group_sets[[1]]
      } else {
        Reduce(intersect, group_sets)
      }
    } else {
      # genes present in ANY species in the group
      unique(unlist(group_sets, use.names = FALSE))
    }
    
    # Intersection with the single species
    inter_vec <- intersect(group_collapsed, single_set)
    per_category[[cat]] <- inter_vec
    
    # Optionally, run SuperExactTest on just two sets: Group vs Single
    # (Only if at least one of the sets is non-empty)
    if (length(group_collapsed) + length(single_set) > 0) {
      two_sets <- list(Group = group_collapsed, Single = single_set)
      universe_size <- length(unique(c(group_collapsed, single_set)))
      if (universe_size > 0) {
        res <- SuperExactTest::supertest(two_sets, n = universe_size)
        supertests[[cat]] <- res
        
        # Write summary CSV
        s   <- summary(res)
        tab <- s$Table
        fn  <- file.path(sum_dir,
                         paste0(gsub(" ", "_", cat), "_",
                                mode, "_",
                                gsub("[^A-Za-z0-9]+", "_", paste(group_readable, collapse="_")),
                                "_vs_",
                                gsub("[^A-Za-z0-9]+", "_", single_readable),
                                "_summary.csv"))
        write.csv(tab, fn, row.names = FALSE)
        
   
      } else {
        supertests[[cat]] <- NULL
      }
    } else {
      supertests[[cat]] <- NULL
    }
  }
  
  list(per_category = per_category, supertests = supertests)
}


# Run SuperExactTest for each category 
supertest_results <- lapply(names(category_sets), function(cat) {
  sets <- category_sets[[cat]]   # species → gene vectors
  
  # Use union of all genes in this category as universe
  universe_size <- length(unique(unlist(sets)))
  
  Result <- supertest(sets, n = universe_size)
  
  list(
    name = cat,
    result = Result
  )
})
names(supertest_results) <- names(category_sets)

# Print summaries
  # make an output folder
dir.create("supertest_summaries", showWarnings = FALSE)

  # one CSV per category
for (cat in names(supertest_results)) {
  res <- supertest_results[[cat]]$result
  s   <- summary(res)              # summary.supertest
  tab <- s$Table                   # data.frame of stats
  
  fn <- file.path(
    "supertest_summaries",
    paste0(gsub(" ", "_", cat), "_summary.csv")
  )
  write.csv(tab, fn, row.names = FALSE)
}

# Export plots
dir.create("supertest_plots", showWarnings = FALSE)

for (cat in names(supertest_results)) {
  res <- supertest_results[[cat]]$result
  safe_name <- gsub(" ", "_", cat)
  
  # 1. Linear plot
  pdf(file.path("supertest_plots", paste0(safe_name, "_linear.pdf")),
      width = 12, height = 10)
  plot(res, Layout="landscape", degree = 2:6, sort.by="size", margin=c(0.5,5,1,2), legend.pos=c(0.95,0.05), cex.label = 0.8, overlap.size.cex = 0.6)
  dev.off()
  
  # 2. Circular plot (legend bottom right)
  pdf(file.path("supertest_plots", paste0(safe_name, "_circular.pdf")),
      width = 12, height = 10)
  plot(res, sort.by="size", degree = 2:6, margin=c(2,2,2,2), color.scale.pos=c(0.85,1), legend.pos=c(0.90,0.02), cex.label = 0.9)
  dev.off()
}



# Example A: genes common to (S. agalactiae, S. uberis) AND present in S. pneumoniae
res_all <- subset_vs_single_intersections(
  category_sets  = category_sets,
  group_species  = c("agal", "uberis", "suis", "iniae", "equi"),   # or c("S. agalactiae","S. uberis")
  single_species = "pneumo",              # or "S. pneumoniae"
  mode           = "all",
  species_labels = species_labels,
  out_prefix     = "non_pneumo_vs_pneumo_all"
)

# Example B: genes present in ANY of (S. equi, S. suis) AND present in S. pneumoniae
res_any <- subset_vs_single_intersections(
  category_sets  = category_sets,
  group_species  = c("equi", "suis"),
  single_species = "pneumo",
  mode           = "any",
  species_labels = species_labels,
  out_prefix     = "non_pneumo_vs_pneumo_any"
)

# Access intersections (vectors of gene IDs) for each category:
# Directory to save results

non_pneumo_vs_pneumo_any <- res_any

dir.create("setdiff_lists_non_pneumo_vs_pneumo_any", showWarnings = FALSE)

lapply(names(category_sets), function(cat) {
  # intersection result from your non_pneumo_vs_pneumo_any run
  group_vs_single_inter <- non_pneumo_vs_pneumo_any$per_category[[cat]]
  
  # 6-species intersection for this category
  all6_inter <- Reduce(intersect, category_sets[[cat]])
  
  # set difference
  diff_genes <- setdiff(group_vs_single_inter, all6_inter)
  
  # write to file, one gene per line
  fn <- file.path("setdiff_lists_non_pneumo_vs_pneumo_any",
                  paste0(gsub(" ", "_", cat), "_diff.txt"))
  writeLines(as.character(diff_genes), fn)
  
  # return for inspection in R
  diff_genes
})


#get setdiff of any intersection for host comparisons

writeLines(
  as.character(
    setdiff(
      Reduce(intersect, category_sets[["Strict core gene"]][c("S. uberis","S. equi","S. suis")]),
      Reduce(intersect, category_sets[["Strict core gene"]])
    )
  ),
  "ruminant_Strict_core_gene_diff.txt"
)

# get gene list for any species/any category
write_species_category <- function(species_code, category, outfile) {
  writeLines(
    as.character(species_gene_categories[[species_code]][[category]]),
    outfile
  )
}

# Example usage:
write_species_category("pneumo", "Cloud gene", "pneumo_Cloud_genes.txt")
