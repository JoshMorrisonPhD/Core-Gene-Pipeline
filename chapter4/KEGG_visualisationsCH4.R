#This code takes KOFAMkoala output .txt files as inputs
#the output is a plot of the relevant KEGG pathways glowing on the global KEGG map
#The logic is super complex, as I had to get rid of a lot of labels to make the plots readable. I would discourage using this code, and advise finding a different way to visualise the KEGG data


# Had a really hard time installing packages with both CRAN and BioC
ensure_pkg <- function(pkgs, bioc = FALSE){
  need <- pkgs[!pkgs %in% rownames(installed.packages())]
  if (length(need)) {
    if (bioc) {
      if (!"BiocManager" %in% rownames(installed.packages())) install.packages("BiocManager")
      BiocManager::install(need, ask = FALSE, update = FALSE)
    } else {
      install.packages(need)
    }
  }
  invisible(lapply(pkgs, require, character.only = TRUE))
}

ensure_pkg(c("readr","dplyr","tidyr","stringr","ggplot2","tibble","purrr","ggrepel",
             "ggfx","ggraph","tidygraph"))
ensure_pkg(c("ggkegg","KEGGREST"), bioc = TRUE)

suppressPackageStartupMessages({
  library(readr); library(dplyr); library(tidyr); library(stringr)
  library(ggplot2); library(tibble); library(purrr); library(ggrepel)
  library(ggfx); library(ggraph); library(tidygraph)
  library(ggkegg); library(KEGGREST)
})

`%||%` <- function(a,b) if (!is.null(a)) a else b

# Species table
species <- tibble::tribble(
  ~label,              ~prefix,
  "agalactiae",        "agal",
  "uberis",            "uberis",
  "suis",              "suis",
  "iniae",             "iniae",
  "equi sbsp equi",    "equi",
  "pneumo",            "pneumo",
  "all",               "all"
)
dir.create("metabolic_maps", showWarnings = FALSE)

# Pulls all KEGG Orthology IDs (format K#####) from a string/vector of strings and returns them as a character vector.
extract_kos <- function(x) unlist(stringr::str_extract_all(x, "K\\d{5}"), use.names = FALSE)

# Reads a KofamKOALA mapper TSV file, extracts KO hits per gene, and returns a de-duplicated two-column tibble: gene_id, KO. Robust to one-column input lines.
parse_kofam_mapper <- function(path){
  if(!file.exists(path)) return(tibble(gene_id=character(), KO=character()))
  df <- suppressWarnings(readr::read_tsv(path, col_names = FALSE, show_col_types = FALSE, progress = FALSE))
  if(ncol(df) == 1) df$X2 <- NA_character_
  df |>
    transmute(gene_id = as.character(X1),
              raw = apply(across(-1), 1, \(r) paste(na.omit(r), collapse=" "))) |>
    rowwise() |> mutate(KO = list(extract_kos(raw))) |> ungroup() |>
    tidyr::unnest_longer(KO, keep_empty = TRUE) |>
    mutate(KO = ifelse(is.na(KO), NA_character_, KO)) |>
    select(gene_id, KO) |> distinct()
}

# Reads a KofamKOALA detail text report, collapses each row, extracts KO hits per gene, honors “starred” hits (if any are present, keeps only starred, non-NA KOs), and returns gene_id, KO unique pairs.
parse_kofam_detail <- function(path){
  if(!file.exists(path)) return(tibble(gene_id=character(), KO=character()))
  raw <- readr::read_lines(path)
  raw <- raw[!grepl("^\\s*$|^#", raw)]
  if(!length(raw)) return(tibble(gene_id=character(), KO=character()))
  parts <- stringr::str_split_fixed(raw, "\\t", n = 20); parts[parts==""] <- NA
  df <- tibble::as_tibble(parts, .name_repair = ~paste0("V", seq_along(.x)))
  gene <- apply(df, 1, function(r){ idx <- which(!is.na(r))[1]; if(length(idx)) r[idx] else NA_character_ })
  row_txt <- apply(df, 1, function(r) paste(na.omit(r), collapse=" "))
  kos <- lapply(row_txt, extract_kos)
  has_star <- grepl("\\*", row_txt)
  out <- tibble::tibble(gene_id = gene, KO = kos, star = has_star) |>
    tidyr::unnest_longer(KO, keep_empty = TRUE)
  if(any(out$star, na.rm=TRUE)) out <- dplyr::filter(out, star, !is.na(KO)) else out <- dplyr::filter(out, !is.na(KO))
  dplyr::select(out, gene_id, KO) |> dplyr::distinct()
}

# Combines a species’ KO assignments from either/both mapper and detail files, preferring mapper if present; appends a species column with label. Returns gene_id, KO, species.
parse_species <- function(label, prefix){
  mapper_file <- paste0(prefix, "_result_ko_mapper.txt")
  detail_file <- paste0(prefix, "_result_kofam.txt")
  if(!file.exists(mapper_file) && !file.exists(detail_file)){
    message("Skipping ", label, " (no files found).")
    return(tibble(gene_id=character(), KO=character(), species=character()))
  }
  gk_mapper <- parse_kofam_mapper(mapper_file)
  gk_detail <- parse_kofam_detail(detail_file)
  gene2ko <- if(nrow(gk_mapper)) gk_mapper else gk_detail
  mutate(gene2ko, species = label)
}

# Given KO IDs, fetches their KEGG pathway mappings via KEGGREST::keggGet, caches results in an RDS-backed environment, and returns a named list where each KO maps to a tibble (pathway, title).
.cache_path_file <- "metabolic_maps/ko_to_pathways_cache_v2.rds"
KO2PATH <- if (file.exists(.cache_path_file)) readRDS(.cache_path_file) else new.env(parent = emptyenv())

ko_to_pathways <- function(kos){
  kos <- unique(kos)
  have <- kos[kos %in% ls(KO2PATH)]
  need <- setdiff(kos, have)
  if (length(need)) {
    for (k in need) {
      rec <- tryCatch(KEGGREST::keggGet(paste0("ko:", k))[[1]], error = function(e) NULL)
      if (!is.null(rec) && "PATHWAY" %in% names(rec) && length(rec$PATHWAY)) {
        ids <- names(rec$PATHWAY)
        titles <- unname(rec$PATHWAY)
        KO2PATH[[k]] <- tibble(pathway = ids, title = titles)
      } else {
        KO2PATH[[k]] <- tibble(pathway = character(0), title = character(0))
      }
      Sys.sleep(0.15)
    }
    saveRDS(KO2PATH, .cache_path_file)
  }
  setNames(lapply(kos, function(k) KO2PATH[[k]]), kos)
}

# Global map
message("Preparing global map (ko01100)...")
g_global <- ggkegg::pathway("ko01100") |> ggkegg::process_line()

# Node labels (deduplicated)
node_pathway_labels <- function(g_base, kos){
  if (!length(kos)) return(list(graph = g_base, labels = tibble()))
  
  hits <- paste0("ko:", kos)
  g2 <- g_base |>
    activate(nodes) |> mutate(node_hit = ggkegg::highlight_set_nodes(hits)) |>
    activate(edges) |> mutate(edge_hit = ggkegg::highlight_set_edges(hits)) |>
    activate(edges)
  
  et <- as_tibble(g2)
  if (!nrow(et)) return(list(graph = g2, labels = tibble()))
  
  et$edge_kos <- lapply(et$name, function(nm){
    ks <- extract_kos(nm %||% "")
    intersect(ks, kos)
  })
  et$edge_has_ko <- (et$edge_hit %in% TRUE) & lengths(et$edge_kos) > 0
  
  ko_map <- ko_to_pathways(unique(unlist(et$edge_kos)))
  nt <- g2 |> activate(nodes) |> as_tibble()
  n <- nrow(nt)
  if (!n) return(list(graph = g2, labels = tibble()))
  
  titles_by_node <- vector("list", n)
  if (any(et$edge_has_ko)) {
    for (i in which(et$edge_has_ko)) {
      fr <- et$from[i]; to <- et$to[i]
      ks <- et$edge_kos[[i]]
      ptitles <- unique(unlist(lapply(ks, function(k){
        km <- ko_map[[k]]
        if (!is.null(km) && nrow(km)) km$title else character(0)
      })))
      if (length(ptitles)) {
        titles_by_node[[fr]] <- unique(c(titles_by_node[[fr]], ptitles))
        titles_by_node[[to]] <- unique(c(titles_by_node[[to]], ptitles))
      }
    }
  }
  
  lab_tbl <- tibble(
    x = nt$x,
    y = nt$y,
    label_text = vapply(titles_by_node, function(v){
      if (!length(v)) return("")
      v <- gsub(" pathway$", "", v, ignore.case = TRUE)
      v <- gsub(" - [A-Za-z].*$", "", v)
      utils::head(v, 1)
    }, character(1))
  ) |> 
    filter(label_text != "") |> 
    # deduplicate: only keep first instance of each label_text
    group_by(label_text) |> 
    slice_head(n = 1) |> 
    ungroup()
  
  list(graph = g2, labels = lab_tbl)
}

# Renders the global KEGG map with glowing highlights for KO-matched nodes/edges and adds non-overlapping pathway labels from node_pathway_labels(). Saves a PNG to outfile.
render_glow_with_pathway_labels <- function(g_base, kos, outfile,
                                            label_size = 3,
                                            width = 10, height = 7, dpi = 300){
  
  if (!length(kos)) {
    message("No KOs; skipping plot for ", outfile)
    return(invisible(NULL))
  }
  
  comp <- node_pathway_labels(g_base, kos)
  g2 <- comp$graph
  lab_tbl <- comp$labels
  
  p <- ggraph(g2, x = x, y = y) +
    geom_node_point(size = 1, aes(color = I(fgcolor),
                                  filter = fgcolor != "none" & type != "line"), na.rm = TRUE) +
    geom_edge_link(width = 0.1, aes(color = I(fgcolor),
                                    filter = type == "line" & fgcolor != "none"), na.rm = TRUE) +
    ggfx::with_outer_glow(
      geom_edge_link(width = 1, aes(color = I(fgcolor), filter = edge_hit),
                     lineend = "round", na.rm = TRUE),
      colour = "red", expand = 3
    ) +
    ggfx::with_outer_glow(
      geom_node_point(size = 2, aes(color = I(fgcolor),
                                    filter = node_hit), na.rm = TRUE),
      colour = "red", expand = 3
    ) +
    ggrepel::geom_text_repel(
      data = lab_tbl, aes(x = x, y = y, label = label_text),
      size = label_size, box.padding = 0.25, point.padding = 0.1,
      min.segment.length = 0, max.overlaps = Inf, seed = 42
    ) +
    theme_void()
  
  ggsave(outfile, p, width = width, height = height, dpi = dpi)
  message("Saved with pathway labels: ", outfile)
}

# ---------- 7) Run ----------
counts <- list()

for(i in seq_len(nrow(species))){
  lbl <- species$label[i]; pre <- species$prefix[i]
  mapper_file <- paste0(pre, "_result_ko_mapper.txt")
  detail_file <- paste0(pre, "_result_kofam.txt")
  
  if(!file.exists(mapper_file) && !file.exists(detail_file)){
    message("Skipping ", lbl, " (no files found)."); next
  }
  
  g2k <- parse_species(lbl, pre)
  kos <- g2k$KO |> stats::na.omit() |> unique()
  
  out_png <- file.path("metabolic_maps", paste0("ko01100_glow_pathway_labels_", gsub("\\s+","_", lbl), ".png"))
  render_glow_with_pathway_labels(g_global, kos, out_png, label_size = 3)
  
  counts[[length(counts)+1]] <- tibble::tibble(
    species = lbl,
    genes = dplyr::n_distinct(g2k$gene_id),
    hits = sum(!is.na(g2k$KO)),
    unique_KOs = length(kos)
  )
}

if (length(counts)) {
  counts <- dplyr::bind_rows(counts)
  readr::write_tsv(counts, file.path("metabolic_maps", "summary_counts.tsv"))
  print(counts)
}

message("Done.")
