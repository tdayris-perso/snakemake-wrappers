#!/usr/bin/R

# This script takes a tximport object and builds a deseq2 dataset
# for each formula given to snakemake.

base::library("tximport");   # Perform actual count importation in R
base::library("readr");      # Read faster!
base::library("jsonlite");   # Importing inferential replicates
base::library("DESeq2");     # Differential Gene expression

# Load txi object
txi_rds_path <- base::as.character(x = snakemake@input[["tximport"]]);
txi <- base::readRDS(
  file = txi_rds_path
);

# Load experimental design
design_path <- base::as.character(x = snakemake@input[["coldata"]]);
design <- utils::read.table(
  file = design_path,
  sep = "\t",
  header = TRUE
);

# Cast formula as formula instead of string
formula <- stats::as.formula(
  object = snakemake@params[["design"]]
);

# Create dds object
dds <- DESeq2::DESeqDataSetFromTximport(
  txi = txi,
  colData = design,
  desogn = formula
);

# Save as RDS
output_path <- base::as.character(x = snakemake@output[["dds"]]);
base::saveRDS(
  obj = dds,
  file = output_path
);
