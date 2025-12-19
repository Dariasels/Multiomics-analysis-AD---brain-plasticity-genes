###### ANAYLYSE OBTAINED RESULT - ANSWER HYPOTHESIS
###### HYPOTHESIS: AD vs control, are there genes related to brainplasticity significantly higher or lower expressed?
###### both arrays: RNAseq GSE184942 and methylation - GSE134379
###### QUESTION OF THIS SCRIPT: overlap significant probes with plasticity related genes
###### MAIN ISSUE OF THIS STUDY: TISSUES ARE DIFFERENT - CBL VS HIPPOCAMPUS VS TEMPORAL GYRUS

# -----------------------------------------------------------------------------------------------------------------
### PACKAGES
### LIBRARIES
library(dplyr)
library(stringr)
library(tidyr)
library(minfi)
library(limma)
library(DMRcate)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# -----------------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------------
### RESULT DIRECTORY
setwd("/home/daria/Applied_high_throughput_analyses/final_scripts_used/results_data")
# -----------------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------------
### OUTPUT DIRECTORY VOOR GENORMALISEERDE, GEFILTERDE EN MERGED RDS FILES + pheno-metadata
directory_merge <- "/home/daria/Applied_high_throughput_analyses/final_scripts_used/results_data/output_merge"
setwd(directory_merge)

# -----------------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------------
## pheno metadata
pheno_data <- "pheno_filled.csv"
pheno <- read.csv(pheno_data)
head(pheno)
# -----------------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------------
### here general info on the dataset: 

# -----------------------------------------------------------------------------------------------------------------
### MERGING WAS NEEDED BECAUSE I HAVE 808 SAMPLES --> 8 batches of 100 and one of 8 to run through filter step

## common probes all batches
merged_common_probes_all <- "merged_common_probes.txt"
common_probes_all <- read.csv(merged_common_probes_all, header = TRUE)
colnames(common_probes_all) <- "probe"
head(common_probes_all)
count_probes <- (nrow(common_probes_all))
count_probes
#453.525 probes
# -----------------------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------------------
merged_mset_filtered <- "merged_mset_filtered.RDS"
# merged_mset_filtered <- readRDS(merged_mset_filtered) # if you want to check the RDS file
# -----------------------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------------------
### AFTER THE MERGE STEP I SPLIT INTO TWO TISSUES: CBL and MTG (cerrebellum and temporal gyrus)
# other stuff in this directory: filtered.RDS per tissue (from last step)
# the filtered is used to normalize:
# check normalized files:
#mtg_norm <- readRDS("CBL_mset_ENmix_norm.RDS")
#cbl_norm <- readRDS("MTG_mset_ENmix_norm.RDS")

# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
### ANALYSIS OF NORMALIZED DATA
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------


## 1. CBL
# -----------------------------------------------------------------------------------------------------------------
directory_limma_resultsCBL <- "/home/daria/Applied_high_throughput_analyses/final_scripts_used/results_data/results/CBL"
setwd(directory_limma_resultsCBL)

# file 1: 
density_plot_CBL <- "density_beta.pdf"
# density_plot_CBL <- pdftools::pdf_convert(density_plot_CBL, dpi = 200)
density_plot_CBL

# file 2:
# ensure last column is named 'probe'
DMP_limma_results_CBL <- "DMP_limma_results.csv"
DMP_limma_results_CBL <- read.csv(DMP_limma_results_CBL)
# head(DMP_limma_results_CBL)
colnames(DMP_limma_results_CBL)[1] <- "probe" 
head(DMP_limma_results_CBL)


# does not contain gene annotations, so first have to map 
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno_df <- data.frame(
  probe = rownames(anno),
  GeneSymbol = anno$UCSC_RefGene_Name,
  stringsAsFactors = FALSE
)

cbl_annot <- merge(DMP_limma_results_CBL, anno_df, by = "probe", all.x = TRUE)
cbl_annot

# file 3:
PCA_CBL <- "PCA.pdf"  


# file 4:
Summary_CBL <- "summary.txt"
read.csv(Summary_CBL)
# 0 significant probes with FDR < 0.05
# expected because this region is not associated with AD

# -------------------------------------------------------------------------------------------------------------
# DMR data:
# -------------------------------------------------------------------------------------------------------------

directory_DMR_resultsCBL <- "~/Applied_high_throughput_analyses/final_scripts_used/results_data/results_beforeage_sex_correction/CBL"

setwd(directory_DMR_resultsCBL)

# file 5:
DMR_summary_DMRcate_CBL <- "DMR_summary_DMRcate.txt" 
DMR_summary_DMRcate_CBL <- read.csv(DMR_summary_DMRcate_CBL)
DMR_summary_DMRcate_CBL

# file 5:
# DMR_cate_objectRDS_CBL <- "DMRcate_object.RDS"
# DMRcate_RDS_CBL <- readRDS(DMR_cate_objectRDS_CBL)
# DMRcate_RDS_CBL
# View(DMRcate)

# file 6: same results as RDS file
DMR_cate_results_table_CBL <- "DMRcate_results_table.csv" 
DMRcate_cbl <- read.csv(DMR_cate_results_table_CBL)
print(head(DMRcate_cbl))
# two significant regions!!

# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
## 2. MTG
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------

directory_limma_resultsMTG <- "/home/daria/Applied_high_throughput_analyses/final_scripts_used/results_data/results/MTG/"
setwd(directory_limma_resultsMTG)

# file 1: 
density_plot_MTG <- "density_beta.pdf"
# density_plot_CBL <- pdftools::pdf_convert(density_plot_CBL, dpi = 200)
density_plot_MTG

# file 2:
# ensure last column is named 'probe'
DMP_limma_results_MTG <- "DMP_limma_results.csv"
DMP_limma_results_MTG <- read.csv(DMP_limma_results_MTG)
# head(DMP_limma_results_CBL)
colnames(DMP_limma_results_MTG)[1] <- "probe" 
head(DMP_limma_results_MTG)


# does not contain gene annotations, so first have to map 
anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno_df <- data.frame(
  probe = rownames(anno),
  GeneSymbol = anno$UCSC_RefGene_Name,
  stringsAsFactors = FALSE
)

mtg_annot <- merge(DMP_limma_results_MTG, anno_df, by = "probe", all.x = TRUE)
mtg_annot

# file 3:
PCA_MTG <- "PCA.pdf"  


# file 4:
Summary_MTG <- "summary.txt"
read.csv(Summary_MTG)
###Probes: 453526,  Significant probes (FDR < 0.05): 40975

# -------------------------------------------------------------------------------------------------------------
#  DMR data
# -------------------------------------------------------------------------------------------------------------


directory_DMR_resultsMTG <- "~/Applied_high_throughput_analyses/final_scripts_used/results_data/results_beforeage_sex_correction/MTG"

setwd(directory_DMR_resultsMTG)

# file 5:
DMR_summary_DMRcate_MTG <- "DMR_summary_DMRcate.txt" 
DMR_summary_DMRcate_MTG <- read.csv(DMR_summary_DMRcate_MTG)
DMR_summary_DMRcate_MTG

# file 5:
# DMR_cate_objectRDS_MTG <- "DMRcate_object.RDS"
# DMRcate_RDS_MTG <- readRDS(DMR_cate_objectRDS_MTG)
# DMRcate_RDS_MTG
# View(DMRcate)

# file 6: same results as RDS file
DMR_cate_results_table_MTG <- "DMRcate_results_table.csv" 
DMRcate_mtg <- read.csv(DMR_cate_results_table_MTG)
print(head(DMRcate_mtg))
length(DMRcate_mtg)
# 8 significant DMR regions


# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
#### LOADING EXTRACTED PLASTICITY GENES GENERATED FROM GO --> see earlier R script
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------


plasticity_genes <- "/home/daria/Applied_high_throughput_analyses/final_scripts_used/results_data/plasticity_genes.csv"
plasticity_genes <- read.delim(plasticity_genes, header = FALSE)
colnames(plasticity_genes) <- "gene_name"
plasticity_genes$gene_name <- as.character(plasticity_genes$gene_name)
print(head(plasticity_genes))
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
#### START ANALYSIS OF PLASTICITY GENES

###1.FIX RELEVANT TABLES:
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# CLEAN CBL and MTG annotated TABLES:

# clean table --> handles multiple genes per probe, takes into account the NA and empty spots

# funtion:
clean_dmp_gene_table <- function(dmp_df, gene_col = "GeneSymbol") {
  dmp_df %>%
    filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "") %>%
    # split multi-gene annotations: "GENE1;GENE2"
    separate_rows(all_of(gene_col), sep = "[;,]") %>%
    mutate(!!gene_col := str_trim(.data[[gene_col]])) %>%
    filter(.data[[gene_col]] != "")
}


dmp_cbl_clean <- clean_dmp_gene_table(cbl_annot, gene_col = "GeneSymbol")
dmp_cbl_clean


dmp_mtg_clean <- clean_dmp_gene_table(mtg_annot, gene_col = "GeneSymbol")
dmp_mtg_clean

# -------------------------------------------------------------------------------------------------------------
# overlap genes from methylation analysis: CBL
# -------------------------------------------------------------------------------------------------------------
# add a filter to detect only those genes that have a p value lower than 5%

plasticity_hits_cbl <- dmp_cbl_clean %>%
  dplyr::filter(
    adj.P.Val < 0.05,                   # significance filter
    GeneSymbol %in% plasticity_genes$gene_name   # overlap with plasticity genes = 0?
  )
plasticity_hits_cbl
# 0 cbl hits

# count how many genes I have: 20.260
unique_genes_cbl <- dmp_cbl_clean %>%
  filter(GeneSymbol != "" & !is.na(GeneSymbol)) %>%
  distinct(GeneSymbol)


# Number of probes with no gene annotation:
empty_gene_probes <- dmp_cbl_clean %>%
  filter(GeneSymbol == "" | is.na(GeneSymbol)) %>%
  distinct(probe)

n_empty_gene_probes <- nrow(empty_gene_probes)
n_empty_gene_probes
#110158

# convert:

plasticity_gene_hits_cbl <- plasticity_hits_cbl %>%
  dplyr::distinct(GeneSymbol)

length(unique(plasticity_gene_hits_cbl$GeneSymbol))
# --> 0


# extract which genes are up and down regulated :

# hypermethylated in AD (logFC > 0), hypomethylated in AD (logFC < 0), mixed (some probes up, some down)
plasticity_direction <- plasticity_hits_cbl %>%
  group_by(GeneSymbol) %>%
  summarise(
    n_probes = n(),
    n_hyper  = sum(logFC > 0),
    n_hypo   = sum(logFC < 0),
    direction = case_when(
      n_hyper > 0 & n_hypo == 0 ~ "hyper_in_AD",
      n_hypo  > 0 & n_hyper == 0 ~ "hypo_in_AD",
      n_hyper > 0 & n_hypo  > 0 ~ "mixed",
      TRUE ~ "no_signal"
    ),
    min_adjP = min(adj.P.Val),
    median_logFC = median(logFC),
    median_delta_beta = median(delta_beta),
    .groups = "drop"
  ) %>%
  arrange(min_adjP)

# summary: 
plasticity_gene_summary_cbl <- plasticity_hits_cbl %>%
  group_by(GeneSymbol) %>%
  summarise(
    n_probes = n(),
    n_hyper = sum(logFC > 0),
    n_hypo  = sum(logFC < 0),
    min_adjP = min(adj.P.Val),
    median_delta_beta = median(delta_beta),
    .groups = "drop"
  ) %>%
  arrange(min_adjP)

# amount of genes that did not hit:
not_hit <- plasticity_genes$gene_name[
  !(plasticity_genes$gene_name %in% plasticity_gene_hits_cbl$GeneSymbol)
]
length(not_hit)
# 374 - all
# -----------------------------------------------------------------------------------------------------------------


# -----------------------------------------------------------------------------------------------------------------
### PLASTICITY HITS MTG

plasticity_hits_MTG <- dmp_mtg_clean %>%
  dplyr::filter(
    adj.P.Val < 0.05,                   # significance filter
    GeneSymbol %in% plasticity_genes$gene_name   # overlap with plasticity genes
  )
# overlapping genes nr of rows: 2.112 out of total rows: 643.739 ~ 20.260 genes


# count how many genes I have:
unique_genes_mtg <- dmp_mtg_clean %>%
  filter(GeneSymbol != "" & !is.na(GeneSymbol)) %>%
  distinct(GeneSymbol)

n_unique_genes_mtg <- nrow(unique_genes_mtg)
n_unique_genes_mtg
# 20260 genes in total

# Number of probes with no gene annotation:

empty_gene_probes <- mtg_annot%>%
  filter(GeneSymbol == "" | is.na(GeneSymbol)) %>%
  distinct(probe)

n_empty_gene_probes <- nrow(empty_gene_probes)
n_empty_gene_probes

#110158 in annot table before filtering


# probe level hits with plasticity genes: 2112

# convert:

plasticity_gene_hits_mtg <- plasticity_hits_MTG %>%
  dplyr::distinct(GeneSymbol)
length(unique(plasticity_gene_hits_mtg$GeneSymbol))

#  246 gene level hits!

### intermezzo: save table in latex format:
genes <- plasticity_genes %>%
  distinct(gene_name) %>%
  arrange(gene_name) %>%
  pull(gene_name)

n_cols <- 7

# Calculate required rows
n_rows <- ceiling(length(genes) / n_cols)

# Pad gene list with NA to fill matrix exactly
genes_padded <- c(
  genes,
  rep(NA, n_rows * n_cols - length(genes))
)

gene_table_wide <- as.data.frame(
  matrix(
    genes_padded,
    ncol = n_cols,
    byrow = FALSE
  )
)

colnames(gene_table_wide) <- paste0("Gene_", seq_len(n_cols))


library(knitr)
library(kableExtra)

latex_table <- kable(
  gene_table_wide,
  format = "latex",
  booktabs = TRUE,
  caption = "Plasticity-related genes extracted from GEO",
  label = "tab:plasticity_genes",
  align = "l"
) %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down"),
    font_size = 9
  )

writeLines(latex_table, "plasticity_genes.tex")


# let's do the same for plasticity genes list

### end intermezzo



# extract which genes are up and down regulated :

# hypermethylated in AD (logFC > 0), hypomethylated in AD (logFC < 0), mixed (some probes up, some down)
plasticity_direction_mtg <- plasticity_hits_MTG %>%
  group_by(GeneSymbol) %>%
  summarise(
    n_probes = n(),
    n_hyper  = sum(logFC > 0),
    n_hypo   = sum(logFC < 0),
    direction = case_when(
      n_hyper > 0 & n_hypo == 0 ~ "hyper_in_AD",
      n_hypo  > 0 & n_hyper == 0 ~ "hypo_in_AD",
      n_hyper > 0 & n_hypo  > 0 ~ "mixed",
      TRUE ~ "no_signal"
    ),
    min_adjP = min(adj.P.Val),
    median_logFC = median(logFC),
    median_delta_beta = median(delta_beta),
    .groups = "drop"
  ) %>%
  arrange(min_adjP)
plasticity_direction_mtg

# amount of genes that did not hit:
not_hit <- plasticity_genes$gene_name[
  !(plasticity_genes$gene_name %in% plasticity_gene_hits_mtg$GeneSymbol)
]
length(not_hit)
#128
# -----------------------------------------------------------------------------------------------------------------

# -----------------------------------------------------------------------------------------------------------------

### PLOTS
fig_dir <- "/home/daria/Applied_high_throughput_analyses/final_scripts_used/figures"
dir.create(fig_dir, showWarnings = FALSE)

           

library(ggplot2)
# -----------------------------------------------------------------------------------------------------------------
# --- Volcano plot: CBL ---

# use annotated CBL table
# useless plot:
# cbl_volcano <- cbl_annot %>%
#   mutate(
#     neg_log10_adjP = -log10(adj.P.Val),
#     is_sig = adj.P.Val < 0.05,
#     is_plasticity = GeneSymbol %in% plasticity_genes$gene_name,
#     group = dplyr::case_when(
#       !is_sig ~ "NS",
#       is_sig & !is_plasticity ~ "Sig_non_plasticity",
#       is_sig & is_plasticity ~ "Sig_plasticity"
#     )
#   )
# 
# p_cbl <- ggplot(cbl_volcano, aes(x = logFC, y = neg_log10_adjP)) +
#   geom_point(aes(color = group), alpha = 0.6, size = 0.5) +
#   scale_color_manual(
#     values = c(
#       "NS" = "grey70",
#       "Sig_non_plasticity" = "red",
#       "Sig_plasticity" = "blue"
#     )
#   ) +
#   theme_minimal() +
#   labs(
#     title = "Volcano plot – CBL DMPs (AD vs Control)",
#     x = "logFC (methylation difference)",
#     y = "-log10(adj.P.Val)",
#     color = "Category"
#   )
# 
# ggsave(
#   filename = file.path(fig_dir, "volcano_CBL_DMPs.png"),
#   plot = p_cbl,
#   width = 6, height = 5, dpi = 300
# )
# -----------------------------------------------------------------------------------------------------------------
#  in paper:
mtg_volcano <- mtg_annot %>%
  mutate(
    neg_log10_adjP = -log10(adj.P.Val),
    is_sig = adj.P.Val < 0.05,
    is_plasticity = GeneSymbol %in% plasticity_genes$gene_name,
    group = dplyr::case_when(
      !is_sig ~ "NS",
      is_sig & !is_plasticity ~ "Sig_non_plasticity",
      is_sig & is_plasticity ~ "Sig_plasticity"
    )
  )

p_mtg <- ggplot(mtg_volcano, aes(x = logFC, y = neg_log10_adjP)) +
  geom_point(aes(color = group), alpha = 0.6, size = 0.5) +
  scale_color_manual(
    values = c(
      "NS" = "grey70",
      "Sig_non_plasticity" = "red",
      "Sig_plasticity" = "blue"
    )
  ) +
  theme_minimal() +
  labs(
    title = "Volcano plot – MTG DMPs (AD vs Control)",
    x = "logFC (methylation difference)",
    y = "-log10(adj.P.Val)",
    color = "Category"
  )

ggsave(
  filename = file.path(fig_dir, "volcano_MTG_DMPs.png"),
  plot = p_mtg,
  width = 6, height = 5, dpi = 300
)

# -------------------------------------------------------------------------------------------------------------
plasticity_mtg_gene_directions <- plasticity_hits_MTG %>%
  group_by(GeneSymbol) %>%
  summarise(
    n_probes = n(),
    n_hyper = sum(logFC > 0),
    n_hypo  = sum(logFC < 0),
    meth_direction = case_when(
      n_hyper > 0 & n_hypo == 0 ~ "hyper_in_AD",
      n_hypo  > 0 & n_hyper == 0 ~ "hypo_in_AD",
      n_hyper > 0 & n_hypo  > 0 ~ "mixed",
      TRUE ~ "no_signal"
    ),
    .groups = "drop"
  )
plasticity_mtg_direction_counts <- dplyr::count(plasticity_mtg_gene_directions, meth_direction)
plasticity_mtg_direction_counts

### meth-direction:
# 1 hyper_in_AD       61
# 2 hypo_in_AD       122
# 3 mixed             63

###---------------------------------------------------------------------------------------------

## 1. barplot:
p_bar <- ggplot(plasticity_mtg_direction_counts, 
                aes(x = meth_direction, y = n, fill = meth_direction)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c(
      "hyper_in_AD" = "#fc8d62",   # orange
      "hypo_in_AD"  = "#66c2a5",   # green-teal
      "mixed"       = "#8da0cb"    # purple
    )
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Direction of DNA Methylation Changes in Plasticity Genes (MTG)",
    x = "Methylation Direction",
    y = "Number of Plasticity Genes",
    fill = "Direction"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 20, hjust = 1)
  )

ggsave(
  filename = file.path(fig_dir, "barplot_plasticity_direction_MTG.png"),
  plot = p_bar,
  width = 7, height = 5, dpi = 300
)

### PIEPLOT

p_pie <- ggplot(plasticity_mtg_direction_counts, 
                aes(x = "", y = n, fill = meth_direction)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(
    values = c(
      "hyper_in_AD" = "#fc8d62",
      "hypo_in_AD"  = "#66c2a5",
      "mixed"       = "#8da0cb"
    )
  ) +
  theme_void(base_size = 14) +
  labs(
    title = "Methylation Direction in Plasticity Genes (MTG)",
    fill = "Direction"
  )

ggsave(
  filename = file.path(fig_dir, "pie_plasticity_direction_MTG.png"),
  plot = p_pie,
  width = 6, height = 6, dpi = 300
)
# still have to add percentages to the pieplot


#### 3. HORIZONTL BARPLOT
p_bar_h <- ggplot(plasticity_mtg_direction_counts, 
                  aes(x = n, y = meth_direction, fill = meth_direction)) +
  geom_col(width = 0.7) +
  scale_fill_manual(
    values = c(
      "hyper_in_AD" = "#fc8d62",
      "hypo_in_AD"  = "#66c2a5",
      "mixed"       = "#8da0cb"
    )
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Direction of DNA Methylation Changes in Plasticity Genes (MTG)",
    x = "Number of Genes",
    y = "Direction",
    fill = "Direction"
  ) +
  theme(
    plot.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(fig_dir, "bar_h_plasticity_direction_MTG.png"),
  plot = p_bar_h,
  width = 7, height = 5, dpi = 300
)



# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------

#### NEXT STEP: overlap DMR 
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------


library(GenomicRanges)
library(tidyr)
library(dplyr)

# split the coordinates from the table --> this turns a string like "chr1:12345-12456" into:chr = "chr1", start = 12345 and end = 12456

dmr_fixed <- DMRcate_mtg %>%
  separate(coord, into = c("chr", "range"), sep = ":") %>%
  separate(range, into = c("start", "end"), sep = "-") %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end)
  )


dmr_gr <- GRanges(
  seqnames = dmr_fixed$chr,
  ranges = IRanges(start = dmr_fixed$start, end = dmr_fixed$end),
  DMR_id = seq_len(nrow(dmr_fixed))
)

anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

anno_df <- data.frame(
  probe = anno$Name,
  chr   = anno$chr,
  pos   = anno$pos,
  GeneSymbol = anno$UCSC_RefGene_Name,
  stringsAsFactors = FALSE
)

mtg_full <- left_join(mtg_annot, anno_df, by = "probe")
mtg_full_clean <- mtg_full %>%
  dplyr::mutate(GeneSymbol = GeneSymbol.x) %>%
  dplyr::select(probe, chr, pos, logFC, adj.P.Val, GeneSymbol, delta_beta) 

probe_gr <- GRanges(
  seqnames = mtg_full_clean$chr,
  ranges = IRanges(start = mtg_full_clean$pos, width = 1),
  mcols = mtg_full_clean %>%
    dplyr::select(probe, logFC, adj.P.Val, GeneSymbol, delta_beta)
)

length(probe_gr)
length(mcols(probe_gr)$mcols.logFC)
mcols(probe_gr)

colnames(mcols(probe_gr)) <- c(
  "probe",
  "logFC",
  "adj.P.Val",
  "GeneSymbol",
  "delta_beta"
)


hits <- findOverlaps(probe_gr, dmr_gr)
length(hits)


overlap_tbl <- data.frame(
  probe = probe_gr$probe[queryHits(hits)],
  DMR_id = dmr_gr$DMR_id[subjectHits(hits)],
  logFC = probe_gr$logFC[queryHits(hits)],
  adj.P.Val = probe_gr$adj.P.Val[queryHits(hits)],
  GeneSymbol = probe_gr$GeneSymbol[queryHits(hits)],
  delta_beta = probe_gr$delta_beta[queryHits(hits)],
  stringsAsFactors = FALSE
)
head(overlap_tbl)

overlap_tbl$direction <- ifelse(overlap_tbl$logFC > 0, "hyper_in_AD", "hypo_in_AD")

dmr_summary <- overlap_tbl %>%
  group_by(DMR_id) %>%
  summarise(
    n_probes = n(),
    n_hyper = sum(direction == "hyper_in_AD"),
    n_hypo  = sum(direction == "hypo_in_AD"),
    mean_logFC = mean(logFC, na.rm = TRUE),
    mean_delta_beta = mean(delta_beta, na.rm = TRUE),
    direction = case_when(
      n_hyper > 0 & n_hypo == 0 ~ "hyper_in_AD",
      n_hypo  > 0 & n_hyper == 0 ~ "hypo_in_AD",
      n_hyper > 0 & n_hypo > 0  ~ "mixed"
    ),
    genes = paste(unique(GeneSymbol[GeneSymbol != ""]), collapse = ";"),
    .groups = "drop"
  )

head(dmr_summary)
# Many DMRs are hypo_in_AD

dmr_gene_summary <- overlap_tbl %>%
  filter(GeneSymbol != "") %>%
  group_by(GeneSymbol) %>%
  summarise(
    n_probes = n(),
    n_hyper = sum(direction == "hyper_in_AD"),
    n_hypo  = sum(direction == "hypo_in_AD"),
    mean_logFC = mean(logFC),
    mean_delta_beta = mean(delta_beta),
    direction = case_when(
      n_hyper > 0 & n_hypo == 0 ~ "hyper_in_AD",
      n_hypo  > 0 & n_hyper == 0 ~ "hypo_in_AD",
      n_hyper > 0 & n_hypo  > 0 ~ "mixed"
    ),
    .groups = "drop"
  ) %>%
  arrange(mean_logFC)

dmr_gene_summary
# 1861 genes that overlap with DMRs (total not filtered yet)

top_hypo <- dmr_gene_summary %>% 
  filter(direction == "hypo_in_AD") %>%
  arrange(mean_logFC) %>%
  head(20)

top_hyper <- dmr_gene_summary %>%
  filter(direction == "hyper_in_AD") %>%
  arrange(desc(mean_logFC)) %>%
  head(20)

top_mixed <- dmr_gene_summary %>%
  filter(direction == "mixed") %>%
  head(20)

head(top_hyper)
head(top_hypo)
head(top_mixed)


dmr_plasticity_hits <- dmr_gene_summary %>%
  filter(GeneSymbol %in% plasticity_genes$gene_name)

dmr_plasticity_hits

#### intermzzo extract table:
genes <- dmr_plasticity_hits %>%
  distinct(GeneSymbol) %>%
  arrange(GeneSymbol) %>%
  pull(GeneSymbol)

n_cols <- 3

# Calculate required rows
n_rows <- ceiling(length(genes) / n_cols)

# Pad gene list with NA to fill matrix exactly
genes_padded <- c(
  genes,
  rep(NA, n_rows * n_cols - length(genes))
)

gene_table_wide <- as.data.frame(
  matrix(
    genes_padded,
    ncol = n_cols,
    byrow = FALSE
  )
)

colnames(gene_table_wide) <- paste0("Gene_", seq_len(n_cols))


library(knitr)
library(kableExtra)

latex_table <- kable(
  gene_table_wide,
  format = "latex",
  booktabs = TRUE,
  caption = "Plasticity genes overlapping significant MTG DMR regions.",
  label = "tab:plasticity_genes",
  align = "l"
) %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down"),
    font_size = 9
  )

writeLines(latex_table, "dmr_plasticity_genes.tex")

#### intermzzo end extract table

plasticity_dmr_direction_counts <- dmr_plasticity_hits %>%
 dplyr::count(direction)

plasticity_dmr_direction_counts


# Data frame
df <- data.frame(
  direction = c("hyper_in_AD", "hypo_in_AD", "mixed"),
  count = c(10, 12, 2)
)

# Add percentages
df <- df %>%
  mutate(
    percent = round(100 * count / sum(count), 1),
    label = paste0(direction, "\n", percent, "%")
  )

# Pie plot
ggplot(df, aes(x = "", y = count, fill = direction)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  geom_label(aes(label = label),
             position = position_stack(vjust = 0.5),
             show.legend = FALSE) +
  labs(
    title = "Methylation Direction of Plasticity Genes in DMRs",
    fill = "Direction"
  ) +
  theme_void() +
  scale_fill_manual(values = c("#E69F00", "#56B4E9", "#009E73"))

### --------------------------------------------------------------------------------
### --------------------------------------------------------------------------------

### MANHATTAN PLOT

dmr_fixed$DMR_id <- seq_len(nrow(dmr_fixed))

## add this for name genes:
dmr_fixed <- dmr_fixed %>%
  left_join(
    dmr_summary %>% dplyr::select(DMR_id, genes),
    by = "DMR_id"
  )
# "NFATC4;NFATC4;NFATC4;NFATC4" --> thats how the genes look like rn
# -->  have to split into a vector and extract unique values and then collapse back into a clean label:

dmr_fixed$clean_genes <- sapply(dmr_fixed$genes, function(x) {
  unique_genes <- unique(unlist(strsplit(x, ";")))
  paste(unique_genes[unique_genes != ""], collapse = ";")
})

plasticity_points <- dmr_fixed %>%
  filter(DMR_id %in% plasticity_dmr_ids)

dmr_fixed$chr <- factor(
  dmr_fixed$chr,
  levels = paste0("chr", c(1:22, "X", "Y"))
)

ggplot(dmr_fixed, aes(x = start, y = chr)) +
  geom_point(
    aes(color = is_plasticity, size = length),
    alpha = 0.7
  ) +
  # Label RED dots only
  geom_text_repel(
    data = plasticity_points,
    aes(label = clean_genes),
    size = 3,
    color = "red",
    min.segment.length = 0,
    segment.color = "grey50",
    box.padding = 0.3,
    max.overlaps = Inf
  ) +
  scale_color_manual(
    values = c("Other DMR" = "grey70", "Plasticity DMR" = "red")
  ) +
  scale_size_continuous(range = c(0.5, 4)) +
  labs(
    title = "DMR Manhattan Plot with Plasticity Gene Labels",
    x = "Genomic Position",
    y = "Chromosome",
    color = "DMR Class",
    size = "DMR Length"
  ) +
  theme_minimal()

    ### --------------------------------------------------------------------------------
    # Genes that have a direct response to neuroplasticity:
    # NFATC4, SORCS2 / SORCS3, PRKCG (Protein kinase C gamma), SHANK2, CNTNAP1, EP300, CSMD1, TBK1, HRH1, UCN (Urocortin), KCNJ10,...
    ### --------------------------------------------------------------------------------

##### normal manhattan plot without the gene names:
plasticity_dmr_ids <- overlap_tbl %>%
  filter(GeneSymbol %in% plasticity_genes$gene_name) %>%
  pull(DMR_id) %>%
  unique()


dmr_fixed$direction <- dmr_summary$direction[match(dmr_fixed$DMR_id, dmr_summary$DMR_id)]

dmr_fixed$is_plasticity <- ifelse(
  dmr_fixed$DMR_id %in% plasticity_dmr_ids,
  "Plasticity DMR",
  "Other DMR"
)

dmr_fixed$chr <- factor(dmr_fixed$chr, levels = paste0("chr", c(1:22, "X", "Y")))
dmr_fixed$length <- dmr_fixed$end - dmr_fixed$start

ggplot(dmr_fixed, aes(x = start, y = chr)) +
  geom_point(
    aes(color = is_plasticity, size = length),
    alpha = 0.8
  ) +
  scale_color_manual(
    values = c("Other DMR" = "grey70", "Plasticity DMR" = "red")
  ) +
  scale_size_continuous(range = c(0.5, 4)) +
  labs(
    title = "Genomic Distribution of DMRs (Plasticity Genes Highlighted)",
    x = "Genomic Position",
    y = "Chromosome",
    color = "DMR Class",
    size = "DMR Length"
  ) +
  theme_minimal()



### without overlap of red dots:
dmr_fixed$chr <- factor(dmr_fixed$chr, levels = paste0("chr", c(1:22, "X", "Y")))

ggplot(dmr_fixed, aes(x = start, y = chr, color = end - start)) +
  geom_point(alpha = 0.7) +
  scale_color_viridis_c(name = "DMR length") +
  labs(
    title = "Genomic Distribution of DMRs",
    x = "Genomic Position",
    y = "Chromosome"
  ) +
  theme_minimal()



### DISTRIBUTION OF DMR LENGTHS
dmr_fixed$length <- dmr_fixed$end - dmr_fixed$start

ggplot(dmr_fixed, aes(x = length)) +
  geom_histogram(bins = 50, fill = "#56B4E9") +
  scale_x_log10() + 
  labs(
    title = "Distribution of DMR Lengths",
    x = "DMR Length (bp, log scale)",
    y = "Count"
  ) +
  theme_minimal()

##boxplots hyper vs hypo

ggplot(dmr_summary, aes(x = direction, y = mean_delta_beta, fill = direction)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Mean Δβ per DMR by Direction")




### --------------------------------------------------------------------------------
### --------------------------------------------------------------------------------

################IMPORTANT overlap DMR and DMP genes MTG
# 1. Genes significant at DMP level
dmp_plasticity_genes <- unique(plasticity_gene_hits_mtg$GeneSymbol)

# 2. Genes present in DMRs
dmr_plasticity_genes <- unique(dmr_plasticity_hits$GeneSymbol)

# 3. Intersection
overlap_dmr_dmp <- intersect(dmp_plasticity_genes, dmr_plasticity_genes)

overlap_dmr_dmp
length(overlap_dmr_dmp)

## 23 genes overlap DMR and DMP genes

### --------------------------------------------------------------------------------
### --------------------------------------------------------------------------------

# I STILL NEED TO LOOK AT THOSE 2 DMR FROM CEREBELLUM
library(dplyr)
library(tidyr)
library(GenomicRanges)

# Fix coordinate columns
dmr_cbl_fixed <- DMRcate_cbl %>%
  separate(coord, into = c("chr", "range"), sep = ":", remove = FALSE) %>%
  separate(range, into = c("start", "end"), sep = "-", remove = TRUE) %>%
  mutate(
    start = as.numeric(start),
    end = as.numeric(end)
  )
dmr_cbl_gr <- GRanges(
  seqnames = dmr_cbl_fixed$chr,
  ranges = IRanges(start = dmr_cbl_fixed$start, end = dmr_cbl_fixed$end),
  DMR_id = seq_len(nrow(dmr_cbl_fixed))
)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

anno <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)

anno_df <- data.frame(
  chr = anno$chr,
  pos = anno$pos,
  gene = anno$UCSC_RefGene_Name,
  stringsAsFactors = FALSE
)

probe_gr <- GRanges(
  seqnames = anno_df$chr,
  ranges = IRanges(start = anno_df$pos, width = 1),
  gene = anno_df$gene
)
hits <- findOverlaps(dmr_cbl_gr, probe_gr)
overlap_genes_raw <- anno_df$gene[subjectHits(hits)]

# Split, deduplicate, and remove empty strings
overlap_genes_cbl <- unique(unlist(strsplit(overlap_genes_raw, ";")))
overlap_genes_cbl <- overlap_genes_cbl[overlap_genes_cbl != ""]
overlap_genes_cbl #--> 4 genes overlap the rgion

cbl_plasticity_overlap <- intersect(overlap_genes_cbl, plasticity_genes)

cbl_plasticity_overlap
length(cbl_plasticity_overlap)


# --> ZERO
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------

# RNAseq analyse:

# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------

library(DESeq2)
library(limma)
library(ggplot2)
library(pheatmap)
library(tidyverse)

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

raw_counts <- read.csv("/home/daria/Applied_high_throughput_analyses/RNAseq_vs2/GSE184942_raw_counts_GRCh38_p13_NCBI.tsv", check.names = FALSE, sep = "\t")

head(raw_counts)
raw_counts
# only has 9 samples
# so one will be removed for downstream analysis

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

metadata <- read.csv("/home/daria/Applied_high_throughput_analyses/RNAseq_vs2/SraRunTable.csv", sep = ",")
head(metadata)
metadata
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

# 5 AD vs 5 control
# the authors flagged that age and sex metadata were missing for this dataset
# paired-end sequencing
# homo sapiens
# hippocampus
# extract columns I need: 'tissue', 'Run'
metadata_small <- metadata[, c("Run", "tissue", "GEO_Accession..exp.")]
head(metadata_small)
metadata_small
# there are no duplicated rows


#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------


# match metadata and count matrix
rownames(raw_counts) <- raw_counts$GeneID
raw_counts <- raw_counts[, -1]   # keep only GSM columns
raw_counts

metadata_small$GSM <- metadata_small$GEO_Accession..exp.

# align metadata to raw_counts
metadata_small <- metadata_small[match(colnames(raw_counts), metadata_small$GSM), ]

# use GSM as rownames, this matches raw_counts columns
rownames(metadata_small) <- metadata_small$GSM

# make tissue a factor
metadata_small$tissue <- factor(metadata_small$tissue)
metadata_small$tissue <- relevel(metadata_small$tissue, ref = "Health")  


# visual inspection
raw_counts_t <- as.data.frame(t(raw_counts))
raw_counts_t$GSM <- rownames(raw_counts_t)
merged_table <- merge(metadata_small, raw_counts_t, by = "GSM")
# View(merged_table)
# merged_table
#---------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------

# construct DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = raw_counts,
  colData = metadata_small,
  design = ~ tissue
)
dds

#---------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------------------------------------------------------------

# filter low-count genes

dds <- dds[rowSums(counts(dds)) >10, ]
dds
#---------------------------------------------------------------------------------------------------------------------------


#~8000 deleted genes
#---------------------------------------------------------------------------------------------------------------------------

# run DESeq2
dds <- DESeq(dds)
resultsNames(dds)

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

library_sizes <- colSums(counts(dds))
summary(library_sizes)
library_sizes
df_lib <- data.frame(
  sample = names(library_sizes),
  library_size = library_sizes,
  condition = metadata_small$tissue
)

ggplot(df_lib, aes(x = condition, y = library_size)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.1, size = 2) +
  scale_y_log10() +
  theme_minimal() +
  labs(
    title = "RNA-seq library size distribution",
    y = "Total reads per sample (log10)",
    x = "Condition"
  )

summary(colSums(counts(dds)))

#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------

# inspect results


res <- results(dds, alpha = 0.05)
res

# resultsNames(dds) 
# change Health vs AD to AD vs Health


# BiocManager::install("apeglm")
# library(apeglm)

# apply log2 fold change shrinkage for more accurate estimates
res <- lfcShrink(
  dds,
  coef = "tissue_AD_vs_Health",
  type = "apeglm"
)
resultsNames(dds)
levels(metadata_small$tissue)
#---------------------------------------------------------------------------------------------------------------------------



#apeglm provides Bayesian shrinkage estimators for effect sizes for a variety of GLM models, using approximation of the posterior for individual coefficients.

#---------------------------------------------------------------------------------------------------------------------------

rld <- rlog(dds)

plotPCA(rld, intgroup="tissue")
png(filename = "/home/daria/Applied_high_throughput_analyses/final_scripts_used/figures/PCA_ADvsH.png", width = 200, height = 900, res = 150)

#Save plot
outpath <- "/home/daria/Applied_high_throughput_analyses/final_scripts_used/figures/MA_plot_AD_vs_Health.png"
plotPCA(rld, intgroup="tissue")
dev.off()
#---------------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------------------------------------------------------------
plotMA(res, ylim=c(-6,6))

# Save plot
png(filename = outpath, width = 1200, height = 900, res = 150)
plotMA(res, ylim = c(-4, 4))
dev.off()

cat("Saved MA plot to:", outpath)
#---------------------------------------------------------------------------------------------------------------------------





#---------------------------------------------------------------------------------------------------------------------------

library(org.Hs.eg.db)
library(AnnotationDbi)
library(tidyverse)
library(ggrepel)

# 1. Convert DESeq2 results to data frame
res_df <- as.data.frame(res)

# 2. Add EntrezID as a column (needed for AnnotationDbi)
res_df$ENTREZID <- rownames(res_df)

# 3. Map Entrez IDs → Gene Symbols + Gene Names
anno <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = res_df$ENTREZID,
  keytype = "ENTREZID",
  columns = c("SYMBOL", "GENENAME")
)

# 4. Align annotation with res_df rows
anno <- anno[match(res_df$ENTREZID, anno$ENTREZID), ]

# 5. Combine annotation into results table
res_df$SYMBOL <- anno$SYMBOL
res_df$GENENAME <- anno$GENENAME

# 6. Add significance categories
res_df <- res_df %>%
  mutate(sig = case_when(
    padj < 0.05 & log2FoldChange > 1  ~ "Up",
    padj < 0.05 & log2FoldChange < -1 ~ "Down",
    TRUE                              ~ "NS"
  ))

# 7. Replace NA gene symbols with Entrez IDs (optional)
res_df$label <- ifelse(is.na(res_df$SYMBOL), res_df$ENTREZID, res_df$SYMBOL)

# 8. Volcano plot WITH GENE NAMES
ggplot(res_df, aes(log2FoldChange, -log10(padj), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(
    data = subset(res_df, sig != "NS"),
    aes(label = label),
    size = 3,
    max.overlaps = 40
  ) +
  scale_color_manual(values = c(
    "Up" = "red",
    "Down" = "blue",
    "NS" = "grey60"
  )) +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  theme_minimal() +
  labs(
    title = "Volcano Plot (AD vs Health)",
    x = "log2 Fold Change",
    y = "-log10(adjusted p-value)",
    color = "Regulation"
  )
#---------------------------------------------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
# head(res_df)
# 
# head(plasticity_gene_hits_mtg)
# # write a function:
# list_overlapping_plasticity_genes_meth_rnaseq = []
# for SYMBOL in res_df
#   for GeneSymbol in plasticity_gene_hits_mtg
#     if SYMBOL in GeneSymbol
#       add. SYMBOL to list_overlapping_plasticity_genes_meth_rnaseq
#       
# 
### ------------------------------------------------------------
### FIND OVERLAPPING GENES BETWEEN RNA-seq (res_df) AND 
### METHYLATION (plasticity_gene_hits_mtg)
### ------------------------------------------------------------

### ------------------------------------------------------------
### 1. Extract SIGNIFICANT RNA-seq genes
###    Criteria: padj < 0.05 and |log2FC| > 1
### ------------------------------------------------------------

rna_sig <- res_df %>%
  filter(!is.na(padj)) %>%                # remove NA
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(SYMBOL) %>% 
  unique()

length(rna_sig)   # how many significant DE genes in RNA-seq?
rna_sig

#### intermezzo: save table
library(dplyr)
library(knitr)
library(kableExtra)

genes <- sort(rna_sig)

n_cols <- 4

n_rows <- ceiling(length(genes) / n_cols)

genes_padded <- c(
  genes,
  rep(NA, n_rows * n_cols - length(genes))
)

gene_table_wide <- as.data.frame(
  matrix(
    genes_padded,
    ncol = n_cols,
    byrow = FALSE
  )
)

colnames(gene_table_wide) <- paste0("Gene_", seq_len(n_cols))


#------------------------------------------------------------------------------

## 1. Get significant DEGs 
deg_sig <- res_df %>%
  filter(!is.na(padj)) %>%
  filter(padj < 0.05, abs(log2FoldChange) > 1)

nrow(deg_sig)  # should be 20

### intermezzo save DEGs:

library(dplyr)
deg_sig_table <- deg_sig %>%
  tibble::rownames_to_column(var = "EntrezID") %>%   # optional
  dplyr::select(SYMBOL, log2FoldChange, padj, baseMean)

deg_sig_table <- deg_sig_table %>%
  dplyr::select(
    SYMBOL,
    log2FoldChange,
    padj,
    baseMean
  ) %>%
  arrange(padj)
library(knitr)
library(kableExtra)

latex_table <- kable(
  deg_sig_table,
  format = "latex",
  booktabs = TRUE,
  caption = "Significantly differentially expressed genes identified by RNA-seq in hippocampus (AD vs.Control).",
  label = "tab:rna_deg",
  align = c("l", "r", "r", "r")
) %>%
  kable_styling(
    latex_options = c("hold_position", "scale_down"),
    font_size = 9
  )

writeLines(latex_table, "rna_deg_table.tex")
### end intermezzo

## 2. Make a simple plasticity vector
plast_vec <- unique(plasticity_genes$gene_name)   # or just plasticity_genes if it's a character vector

## 3. Overlap DEGs with plasticity genes
deg_plasticity <- deg_sig %>%
  filter(SYMBOL %in% plast_vec)

deg_plasticity_genes <- unique(deg_plasticity$SYMBOL)

deg_plasticity_genes
length(deg_plasticity_genes)  # how many of the 20 DEGs are plasticity genes?
# zero


### ------------------------------------------------------------
### 2. Extract SIGNIFICANT methylation genes
###    (plasticity_gene_hits_mtg already contains methylation hits)
### ------------------------------------------------------------

meth_sig <- unique(plasticity_gene_hits_mtg$GeneSymbol)

length(meth_sig)  # how many significant methylation plasticity genes?


### ------------------------------------------------------------
### 3. Intersect: significant in RNA AND methylation
### ------------------------------------------------------------

sig_overlap <- intersect(rna_sig, meth_sig)

sig_overlap
length(sig_overlap)   # 0

# there are no overlapping significant genes related to brain pl
#-----------------------------------------------------------------------------------
#-----------------------------------------------------------------------------------

rna_top50 <- res_df %>%
  arrange(padj) %>%
  head(50) %>%
  pull(SYMBOL)

sig_overlap_top50 <- intersect(rna_top50, meth_sig)
sig_overlap_top50
# NTS only
#----------------------------------------------------------------------------------
#----------------------------------------------------------------------------------
library(dplyr)

### ------------------------------------------------------------
### 1. Get methylation-significant gene symbols
### ------------------------------------------------------------
meth_sig <- unique(plasticity_gene_hits_mtg$GeneSymbol)
meth_sig <- meth_sig[!is.na(meth_sig)]   # remove any NA

length(meth_sig)  
# ### intermezzo export table
# library(dplyr)
# 
# multiomics_clean <- multiomics_table %>%
#   dplyr::select(
#     GeneSymbol,
#     baseMean,
#     pvalue,
#     padj
#   ) %>%
#   mutate(
#     baseMean   = round(baseMean, 1),
#     pvalue     = signif(pvalue, 3),
#     padj       = signif(padj, 3)
#   ) %>%
#   arrange(padj)
# 
# 
# latex_table <- kable(
#   multiomics_clean,
#   format = "latex",
#   booktabs = TRUE,
#   longtable = TRUE,
#   caption = "Multi-omics integration table showing RNA-seq expression statistics for plasticity-related genes with significant DNA methylation changes in the middle temporal gyrus (MTG).",
#   label = "tab:multiomics_plasticity",
#   align = c("l", "r", "r", "r", "r")
# ) %>%
#   kable_styling(
#     latex_options = "repeat_header",
#     font_size = 8
#   )
# 
# writeLines(latex_table, "multiomics_plasticity_table.tex")
library(dplyr)
library(purrr)
library(knitr)
library(kableExtra)

# Split table
tables <- split(
  multiomics_clean,
  ceiling(seq_len(nrow(multiomics_clean)) / 50)
)

# Export each chunk
walk2(tables, seq_along(tables), ~{
  latex <- kable(
    .x,
    format = "latex",
    booktabs = TRUE,
    caption = paste(
      "Multi-omics integration of plasticity genes (part", .y, ")"
    ),
    label = paste0("tab:multiomics_", .y),
    align = c("l", "r", "r", "r", "r")
  ) %>%
    kable_styling(font_size = 8)
  
  writeLines(latex, paste0("multiomics_part_", .y, ".tex"))
})

### end intermzeeo

### ------------------------------------------------------------
### 2. Subset RNA-seq table to ONLY those methylation hits
###    Even if RNA-seq padj is NOT significant!
### ------------------------------------------------------------
rna_for_meth_hits <- res_df %>%
  filter(SYMBOL %in% meth_sig)

nrow(rna_for_meth_hits)   # how many were found in RNA-seq? 238
head(rna_for_meth_hits)


### ------------------------------------------------------------
### 3. Create a clean multi-omics table:
###    Methylation hits + RNA-seq stats (log2FC, padj, baseMean)
### ------------------------------------------------------------
multiomics_table <- plasticity_gene_hits_mtg %>%
  left_join(
    rna_for_meth_hits,
    by = c("GeneSymbol" = "SYMBOL")
  )

########!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# this table includes the significan genes from Methylation study but the statistics from the rnastudy
### Inspect the final table
head(multiomics_table)

## plots


ggplot(multiomics_table, aes(x = baseMean, y = log2FoldChange, label = GeneSymbol)) +
  geom_point(color = "steelblue", size = 3, alpha = 0.8) +
  geom_text_repel(size = 3.5, max.overlaps = 15) +
  scale_x_log10() +
  theme_minimal() +
  labs(
    title = "RNA-seq Expression vs Fold-Change for Methylation-Significant Plasticity Genes",
    x = "Mean RNA-seq Expression (baseMean, log10 scale)",
    y = "RNA-seq log2 Fold Change (AD vs Health)"
  )

########### GOOD PLOT
ggplot(multiomics_table, aes(x = log2FoldChange, y = -log10(padj), label = GeneSymbol)) +
  geom_point(color = "darkorange", size = 3, alpha = 0.8) +
  geom_text_repel(size = 3.5, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "RNA-seq Significance of Methylation-Significant Plasticity Genes",
    x = "RNA-seq log2 Fold Change",
    y = "-log10(adjusted p-value)"
  )
###########
# Each point = a methylation-significant plasticity gene, plotted with:
#   
#   x-axis: RNA-seq log2 fold change (AD vs healthy)
# 
# y-axis: RNA-seq significance (-log10 adjusted p-value)
# 
# So this tells you:
#   
#   For each plasticity gene that changed in methylation, how does its expression behave in RNA-seq?
# 1. Most methylation-significant genes show very small RNA-seq fold changes
# 
# Almost all points cluster tightly around:
#   
#   log2FC ≈ 0
# 
# –log10(padj) < 1 → corresponds to RNA adj. p > 0.1 (not significant)
# 
# → This means:
#   
#   Methylation plasticity changes do not translate into large expression changes in human hippocampus RNA-seq (for AD vs health).
# 
# This is expected because:
#   
#   Epigenetic marks often regulate long-term gene responsiveness
# 
# mRNA levels can remain stable even when regulatory methylation changes occur
# 
# AD is a heterogeneous condition with subtle transcriptional effects

# 🔸 NTS
# 
# log2FC slightly positive
# 
# very high –log10(padj)
# 
# trend toward RNA upregulation
# 
# known neuromodulatory gene
# 
# associated with synaptic plasticity, neuropeptide signaling
# 
# 🔸 F2R
# 
# strong positive log2FC
# 
# moderate RNA significance
# 
# linked to inflammatory and coagulation pathways
# 
# often involved in neurovascular function in disease
# 
# 🔸 ABHD6
# 
# downregulated in RNA
# 
# ABHD6 modulates endocannabinoid signaling
# 
# relevant in neuronal homeostasis
# 
# 🔸 REELIN (RELN), GRIN3A, ITGA8, MYO16
# 
# These sit slightly above the bulk cluster, suggesting:
#   
#   slight RNA response
# 
# possibly tied to their methylation state
# 
# These are classic synaptic or neuronal structural plasticity genes.
# ###########

multiomics_table %>%
  ggplot(aes(x = reorder(GeneSymbol, log2FoldChange), y = log2FoldChange)) +
  geom_col(fill = "darkred", alpha = 0.8) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "RNA-seq log2 Fold Change for Methylation-Significant Plasticity Genes",
    x = "Gene",
    y = "log2 Fold Change (AD vs Health)"
  )

#-------------------------------------------------------------------------------------------------------

multiomics_full <- plasticity_hits_MTG %>%
  left_join(
    res_df,
    by = c("GeneSymbol" = "SYMBOL")
  )

head(multiomics_full)
library(ggrepel)

ggplot(multiomics_full, aes(x = log2FoldChange, y = delta_beta, label = GeneSymbol)) +
  geom_point(color = "darkorange", alpha = 0.7, size = 3) +
  geom_text_repel(size = 3, max.overlaps = 20) +
  theme_minimal() +
  labs(
    title = "RNA-seq Expression vs Methylation Change (Δβ)",
    x = "RNA-seq log2 Fold Change (AD vs Health)",
    y = "Methylation Δβ (AD − ND)"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40")

#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------

### integrate  DMRs with RNA-seq:

#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------
# BiocManager::install("GenomicRanges")
# BiocManager::install("annotatr")

# I already have: dmr_ranges <- extractRanges(dmr, genome = "hg19")
library(GenomicRanges)
library(annotatr)
### 1. Extract genes overlapping DMRs
dmr_ranges <- GRanges(
  seqnames = dmr_fixed$chr,
  ranges   = IRanges(dmr_fixed$start, dmr_fixed$end)
)

annotations <- build_annotations(
  genome = "hg19",
  annotations = "hg19_basicgenes"
)

dmr_ann <- annotate_regions(
  regions = dmr_ranges,
  annotations = annotations,
  ignore.strand = TRUE
)

dmr_genes <- unique(dmr_ann$gene_symbol)
dmr_genes <- dmr_genes[dmr_genes != "" & !is.na(dmr_genes)]

### 2. Significant RNA-seq genes
rna_sig <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1) %>%
  pull(SYMBOL) %>% unique()

### 3. Intersect
dmr_rna_overlap <- intersect(dmr_genes, rna_sig)
dmr_rna_overlap

### ZERO results

#________________________________________
# CREATE final VENN --> used AI for this 5 venn-plot - i calculated the numbers myself
#________________________________________

library(ggplot2)
library(dplyr)


# Function to create circle data
make_circle <- function(x, y, r, n = 100) {
  theta <- seq(0, 2*pi, length.out = n)
  data.frame(
    x = x + r * cos(theta),
    y = y + r * sin(theta)
  )
}

# Create three main circles for MTG (large), HC, and CB
mtg_circle <- make_circle(0, 0, 4)
hc_circle <- make_circle(7, 2, 1.6)
cb_circle <- make_circle(7, -2.5, 2)

# Create nested circles for MTG (DMR vs DMP)
mtg_dmr <- make_circle(-1, 0, 1.7)  # DMR circle (contains 1 + 23)
mtg_dmp <- make_circle(1, 0, 2.8)   # DMP circle (contains 23 + 223)

# Plot
p1 <- ggplot() +
  # MTG outer boundary
  geom_polygon(data = mtg_circle, aes(x, y), 
               fill = NA, color = "purple", linewidth = 1) +
  
  # MTG DMR circle
  geom_polygon(data = mtg_dmr, aes(x, y), 
               fill = "lightblue", alpha = 0.3, color = "steelblue4", linewidth = 1) +
  
  # MTG DMP circle
  geom_polygon(data = mtg_dmp, aes(x, y), 
               fill = "pink", alpha = 0.3, color = "brown3", linewidth = 1) +
  
  # HC circle
  geom_polygon(data = hc_circle, aes(x, y), 
               fill = "darkolivegreen4", alpha = 0.3, color = "olivedrab", 
               linewidth = 1) +
  
  # CB circle
  geom_polygon(data = cb_circle, aes(x, y), 
               fill = "goldenrod2", alpha = 0.3, color = "orange", 
               linewidth = 1) +
  
  # Add numbers
  annotate("text", x = -2.5, y = 0, label = "1", size = 6, fontface = "bold") +
  annotate("text", x = 0, y = 0, label = "23", size = 6, fontface = "bold") +
  annotate("text", x = 2.5, y = 0, label = "223", size = 6, fontface = "bold") +
  annotate("text", x = 7, y = 2, label = "0", size = 6, fontface = "bold") +
  annotate("text", x = 7, y = -2.5, label = "4", size = 6, fontface = "bold") +
  
  # Add labels
  annotate("text", x = 0, y = 5, label = "Middle Temporal Gyrus", size = 5, fontface = "bold", color = "purple4") +
  annotate("text", x = -1, y = -3.2, label = "DMR", size = 5, color = "blue") +
  annotate("text", x = 1, y = -3.2, label = "DMP", size = 5, color = "red") +
  annotate("text", x = 7, y = 5, label = "Hippocampus (RNA-seq)", size = 5, fontface = "bold", color = "olivedrab") +
  annotate("text", x = 7, y = -5, label = "Cerebellum\nDMR", size = 5, fontface = "bold", color = "orange") +
  
  coord_equal() +
  theme_void() +
  labs(title = "Brain-Plasticity Genes: Significant in AD vs Control",
       subtitle = "\n DMP: Differentially Methylated Position | DMR: Differentially Methylated Region") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 10, margin = margin(b = 20))
  )

print(p1)


