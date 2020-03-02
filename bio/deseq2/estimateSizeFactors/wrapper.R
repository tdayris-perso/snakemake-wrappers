#!/usr/bin/R

# This script takes a deseq2 dataset object and estimates
# size factors for further normalization

base::library("DESeq2");     # Differential Gene expression

# Cast input path as character
dds_path <- base::as.character(x = snakemake@input[["dds"]]);
dds <- base::readRDS(dds_path);

# Cast locfunc as function name
extra <- base::as.character(x = snakemake@input[["extra"]]);

# Create object
dds <- base::eval(
  base::parse(
    text = base::paste0(
      "DESeq2::estimateSizeFactors(dds, ", extra, ");"
    )
  )
);
dds <- DESeq2::estimateSizeFactors(dds));
print(dds);

# Save as RDS
output_path <- base::as.character(snakemake@output[["esf"]]);
base::saveRDS(
  obj = dds,
  file = output_path
);
