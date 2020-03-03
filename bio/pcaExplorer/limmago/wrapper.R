#!/usr/bin/R

# This script takes a deseq2 dataset object and performs
# a regularized-logarithmic transformation on it

base::library(package = "DESeq2");        # Differential Gene expression
base::library(package = "pcaExplorer");   # Handles PCA
base::library(package = "DelayedArray");  # Handle in-memory array-like datasets

# Load specified input files
dds_path <- base::as.character(
  x = snakemake@input[["dds"]]
);
dds <- base::readRDS(file = dds_path);

dst_path <- base::as.character(
  x = snakemake@input[["dst"]]
);
dst <- base::readRDS(file = dst_path);


# Building limmago
bg_ids <- IRanges::rownames(x = dds)[
  DelayedArray::rowSums(x = DESeq2::counts(dds)) > 0
];

limmago <- pcaExplorer::limmaquickpca2go(
  se = dst,
  organism = snakemake@params[["organism"]],
  background_genes = bg_ids,
  inputType = "ENSEMBL"
);
limmago_output <- base::as.character(
  x = snakemake@output[["limmago"]]
);
base::saveRDS(
  object = limmago,
  file = limmago_output
);
