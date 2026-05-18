#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

input <- args[1]
max_missing <- args[2]
clustering_method <- args[3]

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
  "results_alleles.tsv",
  delim = "\t",
  col_types = cols(.default = "c")
)

# Clean data
data_clean <- data %>%
  mutate_at(vars(-FILE),
            function(x) str_remove_all(x, "INF-")) %>%
  mutate_at(vars(-FILE),
            function(x) str_replace_all(
              x,
              "(PLOT5)|(PLOT3)|(LNF)|(ASM)|(ALM)|(NIPH)|(NIPHEM)|(PAMA)|(PLNF)|(LOTSC)",
              NA_character_)) %>%
  mutate_at(vars(-FILE),
            as.factor)

# Filter out samples with too many 
# missing alleles
data_filtered <- data_clean %>%
  mutate(NA_count = apply(., 1, function(x) sum(is.na(x)))) %>%
  filter(NA_count <= as.numeric("38")) %>%
  select(-NA_count) %>%
  mutate_at(vars(-FILE),
            as.factor) %>%
  column_to_rownames("FILE")



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
  "dissimilarity_matrix.tsv",
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
        transposed_data[,col[1]] != transposed_data[,col[2]],
        na.rm = TRUE
      ),
      compared_alleles_pair = sum(
        !is.na(transposed_data[,col[1]]) & !is.na(transposed_data[,col[2]])
      ),
      typed_alleles_pair = sum(
        !is.na(transposed_data[,col[1]]) | !is.na(transposed_data[,col[2]])
      ),
      missing_alleles_pair = sum(
        is.na(transposed_data[,col[1]]) | is.na(transposed_data[,col[2]])
      ),
      typed_alleles_isolate1 = sum(
        !is.na(transposed_data[,col[1]])
      ),
      typed_alleles_isolate2 = sum(
        !is.na(transposed_data[,col[2]])
      ),
      missing_alleles_isolate1 = sum(
        is.na(transposed_data[,col[1]])
      ),
      missing_alleles_isolate2 = sum(
        is.na(transposed_data[,col[2]])
      )
    )
  }
  )

write_delim(
  hamming,
  "hamming_distances.tsv",
  delim = "\t"
)

# Run clustering
if (clustering_method == "single") {
  tree <- as.phylo(hclust(dissimilarity, "single"))
} else if (clustering_method == "nj") {
  tree <- nj(dissimilarity)
}

write.tree(
  tree,
  "dendrogram.nwk"
)