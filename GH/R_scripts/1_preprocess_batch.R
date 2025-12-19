#!/usr/bin/env Rscript

cat("=== preprocess_batch.R: Batch preprocessing ===\n")

suppressMessages({
  library(minfi)
  library(ENmix)
  library(S4Vectors)
  library(dplyr)
  library(IlluminaHumanMethylation450kmanifest)
  library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
})

# -----------------------------
# Parse arguments
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: preprocess_batch.R <rgSet_batchX.RDS(.gz)>\n", call. = FALSE)
}
batch_file <- args[1]
batch_name <- tools::file_path_sans_ext(basename(batch_file))

input_dir  <- "/data/gent/501/vsc50116/Meth_analysis_seperate_steps/input_data"
output_dir <- file.path(Sys.getenv("VSC_DATA"), "Meth_analysis_seperate_steps/output_data", batch_name)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Batch file:  ", batch_file, "\n")
cat("Batch name:  ", batch_name, "\n")
cat("Output dir:  ", output_dir, "\n")

# -----------------------------
# Load batch (supports .gz)
# -----------------------------
cat("Loading RGChannelSet...\n")
if (grepl("\\.gz$", batch_file)) {
  rg <- readRDS(gzfile(batch_file))
} else {
  rg <- readRDS(batch_file)
}

cat("Loaded RGChannelSet with dim:", paste(dim(rg), collapse = " x "), "\n")

# -----------------------------
# Attach updated phenotype
# -----------------------------
pheno_file <- file.path(input_dir, "pheno_filled.csv")
if (file.exists(pheno_file)) {
  cat("Loading pheno_filled.csv and matching by GSM_ID...\n")

  pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)
  ph_rg <- as.data.frame(colData(rg))

  if (!("GSM_ID" %in% colnames(ph_rg)))
    stop("GSM_ID missing in RGSet colData.", call. = FALSE)
  if (!("GSM_ID" %in% colnames(pheno)))
    stop("GSM_ID missing in pheno file.", call. = FALSE)

  idx <- match(ph_rg$GSM_ID, pheno$GSM_ID)
  if (any(is.na(idx)))
    stop("Some samples missing in pheno file.", call. = FALSE)

  pheno_batch <- pheno[idx, ]
  rownames(pheno_batch) <- colnames(rg)
  pData(rg) <- S4Vectors::DataFrame(pheno_batch)

  cat("Attached pheno info.\n")
}

# -----------------------------
# Detection P-value QC
# -----------------------------
cat("Computing detection P-values...\n")
detP <- detectionP(rg)

sample_fail_rate <- colMeans(detP > 0.01)
probe_fail_rate  <- rowMeans(detP > 0.01)

bad_samples <- names(sample_fail_rate)[sample_fail_rate > 0.05]
bad_probes  <- rownames(detP)[probe_fail_rate > 0.05]

cat("Bad samples (>5% failed probes):", length(bad_samples), "\n")
cat("Bad probes  (>5% failed samples):", length(bad_probes), "\n")

if (length(bad_samples) > 0) {
  rg <- rg[, !colnames(rg) %in% bad_samples]
  detP <- detP[, !colnames(detP) %in% bad_samples]
  sample_fail_rate <- sample_fail_rate[!names(sample_fail_rate) %in% bad_samples]
}

write.csv(
  data.frame(sample = names(sample_fail_rate),
             sample_fail_rate = sample_fail_rate),
  file.path(output_dir, paste0(batch_name, "_sample_detection_fail_rate.csv")),
  row.names = FALSE
)

writeLines(
  bad_probes,
  file.path(output_dir, paste0(batch_name, "_bad_probes_detectionP.txt"))
)

# -----------------------------
# Minfi QC report (safe)
# -----------------------------
cat("Generating minfi QC report...\n")
qcReport(rg, pdf = file.path(output_dir, paste0(batch_name, "_minfi_QC_report.pdf")))

# -----------------------------
# PreprocessRaw (needed for ENmix)
# -----------------------------
cat("Preprocessing raw signals...\n")
mset_raw <- preprocessRaw(rg)

# -----------------------------
# ENmix background + dye correction
# -----------------------------
cat("Running ENmix (QCinfo disabled)...\n")

NA_CpGs <- rownames(detP)[is.na(rowMeans(getB(mset_raw)))]

mset_bc <- preprocessENmix(
  rg,
  bgParaEst = "oob",
  dyeCorr   = "RELIC",
  QCinfo    = NULL,       # important fix
  exQCsample = TRUE,
  exQCcpg    = TRUE,
  exCpG      = NA_CpGs
)

saveRDS(
  mset_bc,
  file.path(output_dir, paste0(batch_name, "_mset_ENmix_BC.RDS"))
)

# -----------------------------
# Probe filtering (requires genomic)
# -----------------------------
cat("Filtering probes (mapping to genome first)...\n")

gmset_bc <- mapToGenome(mset_bc)   # REQUIRED FIX

mset_f <- dropLociWithSnps(gmset_bc)

keep <- !(grepl("^rs", rownames(mset_f)) | grepl("^ch", rownames(mset_f)))
mset_f <- mset_f[keep, ]

ann <- getAnnotation(mset_f)
autosomes <- !(ann$chr %in% c("chrX", "chrY"))
mset_f <- mset_f[autosomes, ]

saveRDS(
  mset_f,
  file.path(output_dir, paste0(batch_name, "_mset_filtered.RDS"))
)

writeLines(
  rownames(mset_f),
  file.path(output_dir, paste0(batch_name, "_kept_probes.txt"))
)

cat("Filtered probes:", nrow(mset_f), "\n")
cat("=== DONE ===\n")
