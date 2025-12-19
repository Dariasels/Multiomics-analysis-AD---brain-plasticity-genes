#!/usr/bin/env Rscript

cat("=== analyze_tissue.R ===\n")

# Load required packages (on HPC must load Bioconductor bundle before running)
suppressMessages({
  library(minfi)
  library(limma)
  library(ggplot2)
  library(matrixStats)
})

# -----------------------------
# Parse argument
# -----------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: analyze_tissue.R <TISSUE>", call. = FALSE)
}

tissue <- args[1]
cat("Analyzing tissue:", tissue, "\n")

# -----------------------------
# Correct HPC paths
# -----------------------------
input_dir  <- "/scratch/gent/501/vsc50116/AHTA_methylation/scripts/output_merge"
output_dir <- file.path("/scratch/gent/501/vsc50116/AHTA_methylation/scripts/results", tissue)
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

file_in <- file.path(input_dir, paste0(tissue, "_mset_ENmix_norm.RDS"))
if (!file.exists(file_in)) stop("Normalized file not found:\n", file_in)

cat("Loading:", file_in, "\n")
mset <- readRDS(file_in)
pd <- pData(mset)

cat("Samples:", ncol(mset), "\n")

# -----------------------------
# Density plot of beta values
# -----------------------------
pdf(file.path(output_dir, "density_beta.pdf"))
densityPlot(getBeta(mset), 
            main = paste("Density after ENmix quantile normalization -", tissue))
dev.off()

# -----------------------------
# Compute Beta and M values
# -----------------------------
beta <- getBeta(mset)
M <- log2((getMeth(mset) + 1) / (getUnmeth(mset) + 1))

# -----------------------------
# PCA using 20k most variable probes
# -----------------------------
cat("Running PCA...\n")

vsel <- order(matrixStats::rowVars(M), decreasing = TRUE)[1:20000]
pc <- prcomp(t(M[vsel, ]), scale. = TRUE)

pca_df <- data.frame(
  PC1 = pc$x[, 1],
  PC2 = pc$x[, 2],
  Condition = factor(pd$Condition)
)

pdf(file.path(output_dir, "PCA.pdf"))
ggplot(pca_df, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 3) +
  theme_minimal() +
  ggtitle(paste("PCA -", tissue))
dev.off()

# -----------------------------
# LIMMA differential methylation
# -----------------------------
cat("Running LIMMA...\n")

# Your conditions are AD vs ND, so reorder levels
pd$Condition <- factor(pd$Condition, levels = c("ND", "AD"))

pd$Sex <- factor(pd$Sex)
pd$Age <- as.numeric(pd$Age)
design <- model.matrix(~ Condition + Sex + Age, data = pd)

#design <- model.matrix(~ Condition, data = pd)

fit <- lmFit(M, design)
fit <- eBayes(fit)

tt <- topTable(fit, coef = "ConditionAD", number = Inf, adjust = "BH")

# beta means & delta beta
tt$beta_ND <- rowMeans(beta[, pd$Condition == "ND"])
tt$beta_AD <- rowMeans(beta[, pd$Condition == "AD"])
tt$delta_beta <- tt$beta_AD - tt$beta_ND

write.csv(tt, file.path(output_dir, "DMP_limma_results.csv"))

# -----------------------------
# Summary file
# -----------------------------
writeLines(c(
  paste("Samples:", ncol(mset)),
  paste("Probes:", nrow(mset)),
  paste("Significant probes (FDR < 0.05):", sum(tt$adj.P.Val < 0.05))
), file.path(output_dir, "summary.txt"))

cat("=== analyze_tissue.R DONE ===\n")
