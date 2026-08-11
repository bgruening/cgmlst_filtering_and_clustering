#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

get_arg <- function(flag, default = NULL) {
  idx <- match(flag, args)

  if (is.na(idx)) {
    return(default)
  }

  if (idx == length(args)) {
    stop(paste("Missing value for argument:", flag), call. = FALSE)
  }

  args[idx + 1]
}

input <- get_arg("--input")
max_missing <- get_arg("--max-missing")
clustering_method <- get_arg("--clustering-method")

dissimilarity_matrix_out <- get_arg(
  "--dissimilarity-matrix",
  "dissimilarity_matrix.tsv"
)
hamming_distances_out <- get_arg("--hamming-distances", "hamming_distances.tsv")
tree_out <- get_arg("--tree", "dendrogram.nwk")

if (is.null(input)) {
  stop("Missing required argument: --input", call. = FALSE)
}

if (is.null(clustering_method)) {
  stop("Missing required argument: --clustering-method", call. = FALSE)
}

# Load libraries
library(readr)
library(dplyr)
library(tidyr)
library(stringr)
library(tibble)
library(cluster)
library(purrr)
library(ape)

# Import data
data <- read_delim(
  input,
  delim = "\t",
  col_types = cols(.default = "c")
)

# Clean data
data_clean <- data %>%
  mutate_at(vars(-FILE), function(x) str_remove_all(x, "INF-")) %>%
  mutate_at(vars(-FILE), function(x) {
    str_replace_all(
      x,
      "(PLOT5)|(PLOT3)|(LNF)|(ASM)|(ALM)|(NIPH)|(NIPHEM)|(PAMA)|(PLNF)|(LOTSC)",
      NA_character_
    )
  }) %>%
  mutate_at(vars(-FILE), as.factor)

# Filter out samples with too many
# missing alleles
data_filtered <- data_clean %>%
  mutate(NA_count = rowSums(is.na(select(., -FILE)))) %>%
  filter(NA_count <= as.numeric(max_missing)) %>%
  select(-NA_count) %>%
  mutate_at(vars(-FILE), as.factor) %>%
  column_to_rownames("FILE")

if (nrow(data_filtered) < 2) {
  stop("Not enough samples after filtering (need at least 2). Adjust --max-missing to retain more samples.", call. = FALSE)
}

# Calculate dissimilarity matrix
dissimilarity <- daisy(
  data_filtered,
  metric = "gower"
)

dissimilarity_output <- as.data.frame(
  as.matrix(
    dissimilarity
  )
) %>%
  rownames_to_column("FILE")

write_delim(
  dissimilarity_output,
  dissimilarity_matrix_out,
  delim = "\t"
)

# Calculate hamming distance
transposed_data <- as.data.frame(as.matrix(t(data_filtered)))

hamming <- combn(colnames(transposed_data), 2, simplify = FALSE) %>%
  map_df(function(col) {
    data.frame(
      isolate1 = col[1],
      isolate2 = col[2],
      hamming = sum(
        transposed_data[, col[1]] != transposed_data[, col[2]],
        na.rm = TRUE
      ),
      compared_alleles_pair = sum(
        !is.na(transposed_data[, col[1]]) & !is.na(transposed_data[, col[2]])
      ),
      typed_alleles_pair = sum(
        !is.na(transposed_data[, col[1]]) | !is.na(transposed_data[, col[2]])
      ),
      missing_alleles_pair = sum(
        is.na(transposed_data[, col[1]]) | is.na(transposed_data[, col[2]])
      ),
      typed_alleles_isolate1 = sum(
        !is.na(transposed_data[, col[1]])
      ),
      typed_alleles_isolate2 = sum(
        !is.na(transposed_data[, col[2]])
      ),
      missing_alleles_isolate1 = sum(
        is.na(transposed_data[, col[1]])
      ),
      missing_alleles_isolate2 = sum(
        is.na(transposed_data[, col[2]])
      )
    )
  })

write_delim(
  hamming,
  hamming_distances_out,
  delim = "\t"
)

# Run clustering
if (clustering_method == "single") {
  tree <- as.phylo(hclust(dissimilarity, "single"))
} else if (clustering_method == "nj") {
  tree <- nj(dissimilarity)
} else {
  stop("Invalid clustering_method. Use 'single' or 'nj'.", call. = FALSE)
}

write.tree(
  tree,
  tree_out
)
