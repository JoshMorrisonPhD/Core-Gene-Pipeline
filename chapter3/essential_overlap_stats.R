# the input for this is the deduplicated_full_matrix.xlsx file containing all of the essential genes for each species and their matches in other species that was created using
# generate_all_presence_matrices_V2.py and merge2.py.
# The output of this code is all relevant statistics about the overlap of essential genes between combinations of species

# Load necessary libraries
if (!require("readxl")) install.packages("readxl")
if (!require("dplyr")) install.packages("dplyr")
if (!require("SuperExactTest")) install.packages("SuperExactTest")
if (!require("stats")) install.packages("stats")
if (!require("stats4")) install.packages("stats4")


library(readxl)
library(dplyr)
library(SuperExactTest)
library(stats)
library(stats4)

# import data
file_xlsx    <- "deduplicated_full_matrix.xlsx"   # path to master sheet
sheet_name   <- 1                                 # sheet index or name
species_cols <- c("agal", "equi", "iniae",
                  "pneumo", "suis", "uberis")     # EXACT column headers
out_tsv      <- "overlap_stats.tsv"

# read + convert to binary presence matrix 
raw <- read_excel(file_xlsx, sheet = sheet_name)

# rows = ortholog clusters; TRUE if that species column is non-empty
binary <- raw |>
  select(all_of(species_cols)) |>
  mutate(across(everything(), ~ !is.na(.x) & .x != "")) |>
  mutate(cluster = row_number())

#  build species-specific sets of cluster IDs 
ess_sets <- lapply(species_cols, \(sp) binary$cluster[binary[[sp]]])
names(ess_sets) <- species_cols
N <- nrow(binary)                          # universe size

## SuperExactTest  +  BH correction 
st <- supertest(ess_sets, n = N)           # global & per-intersection
res_df <- as.data.frame(summary(st)$Table) # convert to ordinary DF

# add Benjamini–Hochberg q-values
res_df$q_BH <- p.adjust(res_df$P.value, method = "BH")

# build a readable key like "agal&equi&iniae"
res_df$key <- sapply(res_df$Intersections, \(x)
                     paste(sort(unlist(strsplit(x, "&"))), collapse = "&"))

# tidy columns, compute log2OR, add k, sort
res_df <- res_df |>
  rename(
    Observed = Observed.Overlap,
    Expected = Expected.Overlap,
    p_raw    = P.value
  )

#  log2 enrichment
if ("log2odds" %in% names(res_df)) {
  res_df <- res_df |> mutate(log2OR = log2odds)
} else if ("FE" %in% names(res_df)) {          # fold-enrichment column
  res_df <- res_df |> mutate(log2OR = log2(FE))
} else {
  res_df <- res_df |> mutate(
    log2OR = log2((Observed + 1e-9) / (Expected + 1e-9))
  )
}

#  k = number of species in the combination
k_col <- intersect(
  c("Set.size", "set_size", "Degree", "degree"), 
  names(res_df)
)[1]

if (is.na(k_col))
  stop("Cannot find a column that stores intersection size (k).")

res_df <- res_df |>
  mutate(k = .data[[k_col]]) |>
  arrange(k, p_raw)


## write results & print top lines
write.table(res_df, out_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nSignificant overlaps (BH < 0.05):\n")
print(res_df %>% filter(q_BH < 0.05) %>%
        select(key, k, Observed, Expected, log2OR, p_raw, q_BH),
      row.names = FALSE)

cat("\nGlobal multi-set P-value :", signif(st$P.value, 3),
    "\nFull table written to   :", out_tsv, "\n")

