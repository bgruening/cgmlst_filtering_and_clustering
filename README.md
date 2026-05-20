# cgMLST filtering and clustering

Galaxy wrapper for filtering and clustering ChewBBACA cgMLST allele call results.

## Overview

This tool processes a ChewBBACA allele profile matrix, such as a `results_alleles.tsv` file. It removes samples with too many missing allele calls, calculates pairwise distances, and generates a clustering tree.

The tool produces:

- A Gower dissimilarity matrix
- A pairwise Hamming distance table
- A Newick tree

The clustering can be performed using either:

- Single-linkage hierarchical clustering
- Neighbor joining

## Input

The input file should be a tab-delimited ChewBBACA allele call result file.

The file is expected to contain a column named:

```text
FILE
