#!/usr/bin/env Rscript

cat("=== dmr_tissue.R: DMRcate with Effect Size Filtering ===\n")
cat("Start time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("R version:", R.version.string, "\n\n")

suppressMessages({
  library(minfi)
  library(DMRcate)
})

## -----------------------------
## 1. Parse tissue argument
## -----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: dmr_tissue.R <TISSUE>", call. = FALSE)
}

tissue <- args[1]
cat("Tissue:", tissue, "\n")

## -----------------------------
## 2. Paths
## -----------------------------
norm_dir <- "/scratch/gent/501/vsc50116/AHTA_methylation/scripts/output_merge"
res_dir  <- file.path("/scratch/gent/501/vsc50116/AHTA_methylation/scripts/results", tissue)
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

file_in <- file.path(norm_dir, paste0(tissue, "_mset_ENmix_norm.RDS"))
if (!file.exists(file_in)) {
  stop("Normalized MethylSet not found: ", file_in, call. = FALSE)
}

cat("Loading normalized MethylSet from:\n  ", file_in, "\n")
mset <- readRDS(file_in)
pd   <- pData(mset)

## -----------------------------
## 3. Check Condition column
## -----------------------------
if (!("Condition" %in% colnames(pd))) {
  stop("Condition column not found in pData(mset).", call. = FALSE)
}

cat("\nRaw Condition table:\n")
print(table(pd$Condition))

pd$Condition <- factor(pd$Condition, levels = c("ND", "AD"))
cat("\nReleveled Condition table (ND=reference, AD=case):\n")
print(table(pd$Condition))

if (sum(pd$Condition == "ND") < 2 || sum(pd$Condition == "AD") < 2) {
  stop("Need at least 2 samples per condition.", call. = FALSE)
}

## -----------------------------
## 4. Compute M-values and Beta values
## -----------------------------
cat("\nComputing methylation values...\n")
beta <- getBeta(mset)
M <- log2((getMeth(mset) + 1) / (getUnmeth(mset) + 1))
cat("M-value matrix dimensions:", nrow(M), "CpGs x", ncol(M), "samples\n")

## -----------------------------
## 5. EFFECT SIZE FILTERING (NEW!)
## -----------------------------
cat("\n>>> FILTERING BY EFFECT SIZE <<<\n")

# Calculate mean methylation per group
mean_nd <- rowMeans(beta[, pd$Condition == "ND"], na.rm = TRUE)
mean_ad <- rowMeans(beta[, pd$Condition == "AD"], na.rm = TRUE)
delta_beta <- mean_ad - mean_nd

# Report global statistics
cat("\nGlobal methylation statistics:\n")
cat("  Mean beta (ND):", mean(mean_nd, na.rm = TRUE), "\n")
cat("  Mean beta (AD):", mean(mean_ad, na.rm = TRUE), "\n")
cat("  Mean delta beta:", mean(delta_beta, na.rm = TRUE), "\n")
cat("  Median delta beta:", median(delta_beta, na.rm = TRUE), "\n")

# Count CpGs by effect size
cat("\nCpGs by absolute effect size:\n")
cat("  |delta beta| > 0.01:", sum(abs(delta_beta) > 0.01, na.rm = TRUE), "\n")
cat("  |delta beta| > 0.02:", sum(abs(delta_beta) > 0.02, na.rm = TRUE), "\n")
cat("  |delta beta| > 0.05:", sum(abs(delta_beta) > 0.05, na.rm = TRUE), "\n")
cat("  |delta beta| > 0.10:", sum(abs(delta_beta) > 0.10, na.rm = TRUE), "\n")
cat("  |delta beta| > 0.20:", sum(abs(delta_beta) > 0.20, na.rm = TRUE), "\n")

# Adaptive threshold selection based on data
max_abs_delta <- max(abs(delta_beta), na.rm = TRUE)
cat("\n  Maximum |delta beta|:", max_abs_delta, "\n")

# Choose threshold based on what's available in the data
if (sum(abs(delta_beta) > 0.05) >= 100) {
  DELTA_BETA_THRESHOLD <- 0.05
  cat("\n>>> Using threshold: 0.05 (moderate changes) <<<\n")
} else if (sum(abs(delta_beta) > 0.02) >= 1000) {
  DELTA_BETA_THRESHOLD <- 0.02
  cat("\n>>> Using threshold: 0.02 (small but meaningful changes) <<<\n")
} else if (sum(abs(delta_beta) > 0.01) >= 1000) {
  DELTA_BETA_THRESHOLD <- 0.01
  cat("\n>>> Using threshold: 0.01 (subtle changes) <<<\n")
} else {
  # If very few CpGs have even 0.01 difference, use top 5% by effect size
  DELTA_BETA_THRESHOLD <- quantile(abs(delta_beta), 0.95, na.rm = TRUE)
  cat("\n>>> Using adaptive threshold (top 5% by effect size):", DELTA_BETA_THRESHOLD, "<<<\n")
}

# Apply filter
keep_cpgs <- abs(delta_beta) > DELTA_BETA_THRESHOLD

cat("\n>>> Applying effect size filter: |delta beta| >", DELTA_BETA_THRESHOLD, "<<<\n")
cat("CpGs before filtering:", nrow(M), "\n")
cat("CpGs after filtering:", sum(keep_cpgs), "\n")
cat("Percent retained:", round(100 * sum(keep_cpgs) / nrow(M), 2), "%\n")

# Filter M-values
M <- M[keep_cpgs, ]

if (nrow(M) == 0) {
  cat("\n!!! No CpGs pass effect size threshold !!!\n")
  cat("This indicates very small methylation differences between groups.\n")
  
  summary_file <- file.path(res_dir, "DMR_summary_DMRcate.txt")
  writeLines(
    c(
      paste("Tissue:", tissue),
      paste("Analysis date:", Sys.Date()),
      "",
      "=== EFFECT SIZE ANALYSIS ===",
      paste("Maximum |delta beta|:", max_abs_delta),
      paste("CpGs with |delta beta| > 0.01:", sum(abs(delta_beta) > 0.01)),
      paste("CpGs with |delta beta| > 0.02:", sum(abs(delta_beta) > 0.02)),
      paste("CpGs with |delta beta| > 0.05:", sum(abs(delta_beta) > 0.05)),
      "",
      "CONCLUSION:",
      "No CpGs show large methylation differences (> 0.01).",
      "This suggests:",
      "- Very subtle methylation changes in this tissue",
      "- High within-group variability",
      "- Possible lack of strong epigenetic signal",
      "",
      "Recommendations:",
      "1. This tissue may not be suitable for DMR analysis",
      "2. Consider increasing sample size",
      "3. Check for batch effects or confounders",
      "4. Try analyzing other tissues with larger effects"
    ),
    summary_file
  )
  cat("Summary written to:", summary_file, "\n")
  cat("\n=== dmr_tissue.R DONE (no meaningful effect sizes) ===\n")
  quit(save = "no")
}

if (nrow(M) < 100) {
  cat("\n!!! WARNING: Very few CpGs remain (", nrow(M), ") !!!\n")
  cat("DMR analysis may not be meaningful with so few CpGs.\n\n")
}

## -----------------------------
## 6. Design matrix
## -----------------------------

# Ensure Sex is a factor
pd$Sex <- factor(pd$Sex)

# Build design matrix adjusting for Sex (mandatory) and Age (recommended)
design <- model.matrix(~ Condition + Age + Sex, data = pd)

cat("\nDesign matrix including covariates (Condition + Age + Sex):\n")
print(head(design))
#design <- model.matrix(~ Condition, data = pd)
#cat("\nDesign matrix (first 6 rows):\n")
#print(head(design))

## -----------------------------
## 7. DMRcate: cpg.annotate
## -----------------------------
cat("\n>>> Running cpg.annotate() <<<\n")
cat("Parameters: FDR=0.05, arraytype=450K, annotation=ilmn12.hg19\n")

myAnnotation <- suppressMessages(
  cpg.annotate(
    object        = M,
    datatype      = "array",
    what          = "M",
    analysis.type = "differential",
    design        = design,
    coef          = "ConditionAD",
    fdr           = 0.05,
    arraytype     = "450K",
    annotation    = "ilmn12.hg19"
  )
)

n_cpg <- length(myAnnotation@ranges)
cat("\nNumber of significant CpGs (FDR < 0.05):", n_cpg, "\n")

if (n_cpg == 0) {
  cat("\n!!! No CpGs significant at FDR 0.05 !!!\n")
  summary_file <- file.path(res_dir, "DMR_summary_DMRcate.txt")
  writeLines(
    c(
      paste("Tissue:", tissue),
      paste("Analysis date:", Sys.Date()),
      paste("Effect size threshold:", DELTA_BETA_THRESHOLD),
      paste("CpGs tested:", nrow(M)),
      "",
      "No CpGs significant at FDR 0.05 after effect size filtering.",
      "",
      "Try:",
      "- Lowering delta beta threshold to 0.05",
      "- Increasing FDR to 0.10"
    ),
    summary_file
  )
  cat("Summary written to:", summary_file, "\n")
  cat("\n=== dmr_tissue.R DONE (no significant CpGs) ===\n")
  quit(save = "no")
}

## -----------------------------
## 8. DMRcate: dmrcate
## -----------------------------
cat("\n>>> Running dmrcate() <<<\n")
cat("Parameters: lambda=1000, C=2\n")

dmr <- NULL
tryCatch({
  dmr <- dmrcate(myAnnotation, lambda = 1000, C = 2)
}, error = function(e) {
  err_msg <- conditionMessage(e)
  cat("\nError in dmrcate():", err_msg, "\n")
  
  if (grepl("no significant CpGs|no DMRs", err_msg, ignore.case = TRUE)) {
    summary_file <- file.path(res_dir, "DMR_summary_DMRcate.txt")
    writeLines(
      c(
        paste("Tissue:", tissue),
        paste("Analysis date:", Sys.Date()),
        paste("Effect size threshold:", DELTA_BETA_THRESHOLD),
        paste("CpGs tested:", nrow(M)),
        paste("Significant CpGs:", n_cpg),
        "",
        "dmrcate() could not identify DMRs.",
        "Error message:",
        err_msg
      ),
      summary_file
    )
    cat("Summary written to:", summary_file, "\n")
    cat("\n=== dmr_tissue.R DONE (dmrcate failed) ===\n")
    quit(save = "no")
  } else {
    stop(e)
  }
})

## -----------------------------
## 9. Check and process results
## -----------------------------
cat("\nChecking DMR results...\n")
cat("DMR object class:", class(dmr)[1], "\n")

# DMResults objects from newer DMRcate need special handling
dmr_df <- NULL

# Try to convert to data.frame
tryCatch({
  if (class(dmr)[1] == "DMResults") {
    # DMResults is a special S4 class that extends data.frame
    # Access it as a matrix first, then convert
    dmr_df <- as.data.frame(as.matrix(dmr))
    
    # Alternative: try direct slot access
    if (is.null(dmr_df) || nrow(dmr_df) == 0) {
      # Get the actual data from the object
      dmr_df <- data.frame(
        coord = dmr@coord,
        no.cpgs = dmr@no.cpgs,
        min_smoothed_fdr = dmr@min_smoothed_fdr,
        Stouffer = dmr@Stouffer,
        HMFDR = dmr@HMFDR,
        Fisher = dmr@Fisher,
        maxdiff = dmr@maxdiff,
        meandiff = dmr@meandiff,
        stringsAsFactors = FALSE
      )
    }
  } else if (is.data.frame(dmr)) {
    dmr_df <- dmr
  } else {
    # Fallback for other formats
    dmr_df <- as.data.frame(dmr)
  }
}, error = function(e) {
  cat("First conversion attempt failed:", conditionMessage(e), "\n")
  
  # Try alternative method: slot by slot
  tryCatch({
    dmr_df <<- data.frame(
      coord = slot(dmr, "coord"),
      no.cpgs = slot(dmr, "no.cpgs"),
      min_smoothed_fdr = slot(dmr, "min_smoothed_fdr"),
      Stouffer = slot(dmr, "Stouffer"),
      HMFDR = slot(dmr, "HMFDR"),
      Fisher = slot(dmr, "Fisher"),
      maxdiff = slot(dmr, "maxdiff"),
      meandiff = slot(dmr, "meandiff"),
      stringsAsFactors = FALSE
    )
    cat("Successfully extracted via slot access\n")
  }, error = function(e2) {
    cat("Slot access also failed:", conditionMessage(e2), "\n")
  })

})

# Check if we successfully got results
if (is.null(dmr_df) || nrow(dmr_df) == 0) {
  cat("\n!!! Could not extract DMR results or no DMRs found !!!\n")
  cat("Object structure:\n")
  print(str(dmr, max.level = 1))
  
  summary_file <- file.path(res_dir, "DMR_summary_DMRcate.txt")
  writeLines(
    c(
      paste("Tissue:", tissue),
      paste("Analysis date:", Sys.Date()),
      paste("Effect size threshold:", DELTA_BETA_THRESHOLD),
      paste("CpGs tested:", nrow(M)),
      paste("Significant CpGs:", n_cpg),
      "",
      "Could not extract DMR results or no DMRs found.",
      "",
      "Debug object saved for troubleshooting."
    ),
    summary_file
  )
  
  saveRDS(dmr, file.path(res_dir, "DMRcate_object_debug.RDS"), compress = "xz")
  cat("Summary and debug object saved\n")
  cat("\n=== dmr_tissue.R DONE (extraction failed) ===\n")
  quit(save = "no")
}

n_dmrs <- nrow(dmr_df)
cat("Successfully extracted", n_dmrs, "DMRs\n")

if (n_dmrs == 0) {
  cat("\n!!! No DMRs found (0 rows) !!!\n")
  summary_file <- file.path(res_dir, "DMR_summary_DMRcate.txt")
  writeLines(
    c(
      paste("Tissue:", tissue),
      paste("Analysis date:", Sys.Date()),
      paste("Effect size threshold:", DELTA_BETA_THRESHOLD),
      paste("CpGs tested:", nrow(M)),
      paste("Significant CpGs:", n_cpg),
      "",
      "dmrcate() returned 0 DMRs.",
      "",
      "Recommendations:",
      "- Decrease C parameter to 1",
      "- Decrease lambda to 500"
    ),
    summary_file
  )
  cat("Summary written to:", summary_file, "\n")
  cat("\n=== dmr_tissue.R DONE (no DMRs) ===\n")
  quit(save = "no")
} else {
  cat("\n>>> Number of DMRs identified:", n_dmrs, "<<<\n")

  # Summary statistics
  if (n_dmrs > 0) {
    cat("\nDMR Summary Statistics:\n")

    if ("Stouffer" %in% colnames(dmr_df)) {
      cat("  Min Stouffer p-value:", min(dmr_df$Stouffer, na.rm = TRUE), "\n")
      cat("  Max Stouffer p-value:", max(dmr_df$Stouffer, na.rm = TRUE), "\n")
    }

    if ("no.cpgs" %in% colnames(dmr_df)) {
      cat("  Mean CpGs per DMR:", mean(dmr_df$no.cpgs, na.rm = TRUE), "\n")
    }

    # Calculate width from coord if available
    if ("coord" %in% colnames(dmr_df)) {
      coords <- strsplit(as.character(dmr_df$coord), "[:-]")

      widths <- sapply(coords, function(x) {
        if (length(x) == 3) {
          as.numeric(x[3]) - as.numeric(x[2])
        } else {
          NA
        }
      })

      if (sum(!is.na(widths)) > 0) {
        cat("  Median DMR width (bp):", median(widths, na.rm = TRUE), "\n")
      }
    }
  }   # <-- closes: if (n_dmrs > 0)

}     # <-- closes: else {

## -----------------------------
## 10. Try to extract ranges
## -----------------------------
dmr_ranges <- NULL

tryCatch({
  cat("\nAttempting to extract GRanges object...\n")
  dmr_ranges <- extractRanges(dmr, genome = "hg19")
  cat("Successfully extracted", length(dmr_ranges), "ranges\n")
}, error = function(e) {
  cat("Note: Could not extract ranges:\n")
  cat("  ", conditionMessage(e), "\n")
})

## -----------------------------
## 11. Save outputs
## -----------------------------
cat("\n>>> Saving results <<<\n")

saveRDS(dmr, file.path(res_dir, "DMRcate_object.RDS"), compress = "xz")
write.csv(dmr_df, file.path(res_dir, "DMRcate_results_table.csv"), row.names = FALSE)

if (!is.null(dmr_ranges)) {
  saveRDS(dmr_ranges, file.path(res_dir, "DMRcate_ranges.RDS"), compress = "xz")
  write.csv(as.data.frame(dmr_ranges), 
            file.path(res_dir, "DMRcate_ranges.csv"), 
            row.names = FALSE)
}

summary_file <- file.path(res_dir, "DMR_summary_DMRcate.txt")
writeLines(
  c(
    paste("Tissue:", tissue),
    paste("Analysis date:", Sys.Date()),
    "",
    "=== Filtering Parameters ===",
    paste("Effect size threshold: |delta beta| >", DELTA_BETA_THRESHOLD),
    "",
    "=== Analysis Parameters ===",
    "Array type: 450K",
    "Annotation: ilmn12.hg19",
    "FDR threshold: 0.05",
    "Lambda (smoothing): 1000",
    "C (CpG threshold): 2",
    "",
    "=== Sample Information ===",
    paste("Total samples:", ncol(M)),
    paste("ND (control):", sum(pd$Condition == "ND")),
    paste("AD (case):", sum(pd$Condition == "AD")),
    "",
    "=== Results ===",
    paste("Total CpGs on array:", length(keep_cpgs)),
    paste("CpGs passing effect size filter:", nrow(M)),
    paste("Significant CpGs (FDR < 0.05):", n_cpg),
    paste("Number of DMRs:", n_dmrs),
    "",
    "=== Output Files ===",
    "- DMRcate_object.RDS (full R object)",
    "- DMRcate_results_table.csv (main results)",
    if (!is.null(dmr_ranges)) "- DMRcate_ranges.RDS (GRanges object)" else NULL,
    if (!is.null(dmr_ranges)) "- DMRcate_ranges.csv (GRanges table)" else NULL
  ),
  summary_file
)

cat("Saved all outputs to:", res_dir, "\n")
cat("\n=== dmr_tissue.R COMPLETED SUCCESSFULLY ===\n")
cat("End time:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

