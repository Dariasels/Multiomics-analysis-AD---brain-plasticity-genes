#!/usr/bin/env Rscript

cat("=== split_tissues.R ===\n")

suppressMessages({
  library(minfi)
  library(S4Vectors)
})

# Input files located in output_merge/
merged_file <- "/scratch/gent/501/vsc50116/AHTA_methylation/scripts/output_merge/merged_mset_filtered.RDS"

# Output in same directory
out_dir <- "/scratch/gent/501/vsc50116/AHTA_methylation/scripts/output_merge"

dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)

merged <- readRDS(merged_file)
pd <- pData(merged)

stopifnot("Tissue" %in% colnames(pd))

cat("Available tissue types:\n")
print(table(pd$Tissue))

# Split
CBL <- merged[, pd$Tissue == "CBL"]
MTG <- merged[, pd$Tissue == "MTG"]

cat("CBL:", ncol(CBL), "samples\n")
cat("MTG:", ncol(MTG), "samples\n")

# Save in output_merge/
saveRDS(CBL, file.path(out_dir, "CBL_mset_filtered.RDS"), compress = "xz")
saveRDS(MTG, file.path(out_dir, "MTG_mset_filtered.RDS"), compress = "xz")

cat("=== split_tissues.R DONE ===\n")
