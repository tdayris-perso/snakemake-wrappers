#!/usr/bin/R

# This script takes a deseq2 dataset object, a deseq2
# transformed data counts object, and an organism
# name, then returns a limma quick pca 2 gene onthology

base::library(package = "DESeq2");        # Differential Gene expression
base::library(package = "pcaExplorer");   # Handles PCA
base::library(package = "DelayedArray");  # Handle in-memory array-like datasets
base::library(package = "org.Hs.eg.db");  # Human genome annotation
base::library(package = "org.Mm.eg.db");  # Mouse genome annotation

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
print(head(bg_ids))

limmago <- pcaExplorer::limmaquickpca2go(
  se = dst,
  organism = snakemake@params[["organism"]],
  background_genes = bg_ids,
  inputType = "HUGO"
);
limmago_output <- base::as.character(
  x = snakemake@output[["limmago"]]
);
base::saveRDS(
  object = limmago,
  file = limmago_output
);
