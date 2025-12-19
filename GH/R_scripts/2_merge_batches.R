#!/usr/bin/env Rscript

cat("=== merge_batches.R: Merging filtered batches ===\n")

suppressMessages({
  library(minfi)
  library(S4Vectors)
  library(BiocGenerics) 
})

# ------------------------
# DIRECTORIES
# ------------------------
input_dir <- "/scratch/gent/501/vsc50116/AHTA_methylation/scripts"
#input_dir  <- "/data/gent/501/vsc50116/Meth_analysis_seperate_steps/output_data/filtered_for_merge"
#output_dir <- "/data/gent/501/vsc50116/Meth_analysis_seperate_steps/output_data/merged"
output_dir <- "/scratch/gent/501/vsc50116/AHTA_methylation/scripts/output_merge"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ------------------------
# FIND FILES
# ------------------------
files <- Sys.glob(file.path(input_dir, "*_mset_filtered.RDS"))

cat("Found filtered batch files:\n")
print(files)

if (length(files) == 0) {
  stop("ERROR: No *_mset_filtered.RDS files found in filtered_for_merge directory!")
}

# ------------------------
# LOAD BATCHES
# ------------------------
msets <- lapply(files, function(f) {
  cat("Loading:", f, "\n")
  readRDS(f)
})

cat("Loaded", length(msets), "batches\n\n")

# ------------------------
# INTERSECT PROBES
# ------------------------
probe_lists <- lapply(msets, rownames)
common_probes <- Reduce(intersect, probe_lists)

cat("Common probes:", format(length(common_probes), big.mark=","), "\n")

if (length(common_probes) < 300000) {
  warning("WARNING: common probes unexpectedly low — check inputs!")
}

# Filter each batch to common probes
msets <- lapply(msets, function(m) m[common_probes, ])

# ------------------------
# MERGE SAMPLES SAFELY
# ------------------------
cat("Merging batches...\n")

merged <- msets[[1]]
if (length(msets) > 1) {
  for (i in 2:length(msets)) {
    cat(" -> Adding batch", i, "\n")
    merged <- BiocGenerics::combine(merged, msets[[i]])
#merged <- minfi::combine(merged, msets[[i]])
  }
}

cat("Merged dataset dim:", paste(dim(merged), collapse=" x "), "\n")

# ------------------------
# SAVE OUTPUT
# ------------------------
saveRDS(
  merged,
  file.path(output_dir, "merged_mset_filtered.RDS"),
  compress = "xz"
)

writeLines(
  common_probes,
  file.path(output_dir, "merged_common_probes.txt")
)

cat("=== merge_batches.R DONE ===\n")
