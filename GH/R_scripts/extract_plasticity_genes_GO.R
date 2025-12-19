###################################
# LOOK FOR PLASTICITY RELATED GENES:
###################################

install.packages("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")

suppressMessages({
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(dplyr)
})

plasticity_go_terms <- c(
  "GO:0060291",
  "GO:0048167",
  "GO:0048169",
  "GO:0007611",
  "GO:0007613",
  "GO:0048812",
  "GO:0048919",
  "GO:0060292",
  "GO:0048699",
  "GO:0043679"
)

genes <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = plasticity_go_terms,
  columns = c("SYMBOL"),
  keytype = "GO"
)

plasticity_genes <- sort(unique(genes$SYMBOL[!is.na(genes$SYMBOL)]))

writeLines(plasticity_genes, "plasticity_genes.txt")
write.csv(data.frame(Gene = plasticity_genes),
          "plasticity_genes.csv",
          row.names = FALSE)

plasticity_genes
