#!/usr/bin/env Rscript

library(minfi)
library(ENmix)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("Usage: normalize_ENmix.R <TISSUE>")

tissue <- args[1]

input_file  <- paste0("/scratch/gent/501/vsc50116/AHTA_methylation/scripts/output_merge/", tissue, "_mset_filtered.RDS")
output_file <- paste0("/scratch/gent/501/vsc50116/AHTA_methylation/scripts/output_merge/", tissue, "_mset_ENmix_norm.RDS")

cb <- readRDS(input_file)

mset_cb <- MethylSet(
  Meth = getMeth(cb),
  Unmeth = getUnmeth(cb),
  colData = colData(cb),
  annotation = annotation(cb)
)

cb_norm <- norm.quantile(mset_cb, method = "quantile2")

saveRDS(cb_norm, output_file, compress = "xz")

cat("Finished ENmix normalization for", tissue, "\n")
