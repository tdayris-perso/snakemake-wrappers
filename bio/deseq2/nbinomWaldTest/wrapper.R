#!/usr/bin/R

# This script takes a deseq2 dataset object and performs
# a negative binomial wald test on it

base::library("DESeq2");     # Differential Gene expression

# Load DESeq2 dataset
dds_path <- base::as.character(
  x = snakemake@input[["dds"]]
);
dds <- base::readRDS(file = dds_path);

# Build extra parameters for DESeq2 nbinomWaldTest
nbinom_extra <- "";
if ("nbinom_extra" %in% snakemake@params) {
  nbinom_extra <- base::paste0(
    ", ",
    base::as.character(x = snakemake@params[["nbinom_extra"]])
  );
}

# Create object
wald <- base::eval(
  base::parse(
    text = base::paste0(
      "DESeq2::nbinomWaldTest(object = dds", nbinom_extra, ");"
    )
  )
);

# Save results
output_rds <- base::as.character(
  x = snakemake@output[["rds"]]
);
base::saveRDS(
  obj = wald,
  file = output_rds
);


names <- DESeq2::resultsNames(
  object = wald
);

output_prefix <- snakemake@output[["tsv"]];
if (! base::file.exists(output_prefix)) {
  base::dir.create(
    path = output_prefix,
    recursive = TRUE
  );
}

alpha_threshold <- 0.5;
if ("alpha_threshold" %in% names(snakemake@params)) {
  alpha_threshold <- base::as.numeric(
    x = snakemake@params[["alpha_threshold"]]
  );
}

fc_threshold <- 0.001;
if ("fc_threshold" %in% names(snakemake@params)) {
  fc_threshold <- base::as.numeric(
    x = snakemake@params[["fc_threshold"]]
  );
}

for (resultname in names) {
  results_frame <- DESeq2::results(
    object = wald,
    name = resultname,
    independentFiltering = TRUE,
    alpha = alpha_threshold,
    lfcThreshold = fc_threshold,
    pAdjustMethod = "BH",
    cooksCutoff = TRUE
  );

  results_path <- base::file.path(
    output_prefix,
    base::paste0("Deseq2_", resultname)
  );

  utils::write.table(
    x = results_frame,
    file = results_path,
    quote = FALSE,
    sep = "\t",
    row.names = FALSE
  );
}
