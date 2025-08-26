#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# bulk_brca_tcga.R
# Process TCGA BRCA bulk RNA-seq (cBioPortal exports):
# - load expression matrices from three TCGA BRCA studies
# - align genes, union samples
# - TMM normalization (edgeR) -> logCPM
# - subset to 59-gene panel (PAM50 + CHD1, CLDN3/4/5/7, OCLN, ALDH1A1/2/3)
# - PAM clustering (k = 5)
# - write processed matrices + cluster assignments
#
# Inputs (expected relative paths; cBioPortal "RNA_Seq_v2_expression_median" txt):
#   data/raw/tcga/brca_tcga/data_RNA_Seq_v2_expression_median.txt
#   data/raw/tcga/brca_tcga_pan_can_atlas_2018/data_RNA_Seq_v2_expression_median.txt
#   data/raw/tcga/brca_tcga_pub2015/data_RNA_Seq_v2_expression_median.txt
#
# Outputs:
#   data/processed/tcga/logCPM_all.tsv
#   data/processed/tcga/logCPM_59genes.tsv
#   results/tables/tcga_pam_clusters.tsv
#
# Notes:
# - cBioPortal “RNA_Seq_v2_expression_median” tables are already normalized
#   expression matrices; we apply TMM+CPM here to match your published pipeline.
# - This script focuses on bulk RNA-seq only (no proteomics, RPPA, miRNA).
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(cluster)   # PAM
})

ensure_dirs <- function(paths) {
  for (p in paths) if (!dir.exists(p)) dir.create(p, recursive = TRUE)
}

# --- I/O paths ----------------------------------------------------------------

in_tcga_legacy <- "data/raw/tcga/data_RNA_Seq_v2_expression_median_1.txt.bz2" 
in_tcga_pancan <- "data/raw/tcga/data_RNA_Seq_v2_expression_median_2.txt.bz2" 
in_tcga_cell15 <- "data/raw/tcga/data_RNA_Seq_v2_expression_median_3.txt.bz2"

out_dir_processed <- "data/processed/tcga"
out_dir_results   <- "results/tables"

ensure_dirs(c(out_dir_processed, out_dir_results))

# --- load expression tables ---------------------------------------------------

read_cbio_expr <- function(path) {
  df <- read_tsv(path, show_col_types = FALSE, progress = FALSE)
  
  # basic sanity checks
  req <- c("Hugo_Symbol", "Entrez_Gene_Id")
  missing_cols <- setdiff(req, names(df))
  if (length(missing_cols)) {
    stop(sprintf("Missing required columns in %s: %s",
                 path, paste(missing_cols, collapse = ", ")))
  }
  
  # normalize a couple of historical quirks
  df <- df %>%
    mutate(
      Hugo_Symbol = ifelse(Hugo_Symbol == "ORC6L", "ORC6", Hugo_Symbol),
      Entrez_Gene_Id = as.character(Entrez_Gene_Id)
    )
  
  df
}

message("Reading TCGA BRCA expression matrices...")
legacy  <- read_cbio_expr(in_tcga_legacy)   # "..._1.txt.bz2"
pancan  <- read_cbio_expr(in_tcga_pancan)   # "..._2.txt.bz2"
cell15  <- read_cbio_expr(in_tcga_cell15)   # "..._3.txt.bz2"

# --- merge gene sets across studies (by Hugo_Symbol) --------------------------

# Intersect genes across the two largest tables (legacy & pancan) to be conservative
genes_common <- Reduce(intersect, list(legacy$Hugo_Symbol, pancan$Hugo_Symbol, cell15$Hugo_Symbol))

legacy  <- legacy  %>% filter(Hugo_Symbol %in% genes_common) %>% distinct(Hugo_Symbol, .keep_all = TRUE)
pancan  <- pancan  %>% filter(Hugo_Symbol %in% genes_common) %>% distinct(Hugo_Symbol, .keep_all = TRUE)
cell15  <- cell15  %>% filter(Hugo_Symbol %in% genes_common) %>% distinct(Hugo_Symbol, .keep_all = TRUE)

legacy  <- legacy %>% arrange(Hugo_Symbol)
pancan  <- pancan %>% arrange(Hugo_Symbol)
cell15  <- cell15 %>% arrange(Hugo_Symbol)

stopifnot(identical(legacy$Hugo_Symbol, pancan$Hugo_Symbol))
stopifnot(identical(legacy$Hugo_Symbol, cell15$Hugo_Symbol))

# union of samples (columns 3: end are samples)
samples_legacy <- setdiff(names(legacy), c("Hugo_Symbol","Entrez_Gene_Id"))
samples_pancan <- setdiff(names(pancan), c("Hugo_Symbol","Entrez_Gene_Id"))
samples_cell15 <- setdiff(names(cell15), c("Hugo_Symbol","Entrez_Gene_Id"))

# combine into a single matrix with all unique samples (prefer pancan values where duplicated)
expr_tcga <- pancan
extra_legacy <- setdiff(samples_legacy, names(expr_tcga))
extra_cell15 <- setdiff(samples_cell15, names(expr_tcga))

if (length(extra_legacy) > 0) {
  expr_tcga <- cbind(expr_tcga, legacy[, extra_legacy])
}
if (length(extra_cell15) > 0) {
  expr_tcga <- cbind(expr_tcga, cell15[, extra_cell15])
}

# --- TMM normalization -> logCPM ----------------------------------------------

mat <- expr_tcga %>%
  select(-Entrez_Gene_Id) %>%
  tibble::column_to_rownames("Hugo_Symbol") %>%
  as.matrix()

# numeric coerce (cBio often comes as numeric already, but be safe)
mode(mat) <- "numeric"

dge <- DGEList(counts = mat)
# We do not filterByExpr aggressively here to retain PAM50/ALDH genes downstream.
dge <- calcNormFactors(dge)
logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 1)

logCPM_df <- as.data.frame(logCPM) %>%
  tibble::rownames_to_column("Symbol")

# --- 59-gene panel ------------------------------------------------------------

pam50 <- c("ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1",
           "CDC20","CDC6","CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2",
           "ESR1","EXO1","FGFR4","FOXA1","FOXC1","GPR160","GRB7","KIF2C",
           "KRT14","KRT17","KRT5","MAPT","MDM2","MELK","MIA","MKI67","MLPH",
           "MMP11","MYBL2","MYC","NAT1","NDC80","NUF2","ORC6","PGR","PHGDH",
           "PTTG1","RRM2","SFRP1","SLC39A6","TMEM45B","TYMS","UBE2C","UBE2T")
claudin_low <- c("CHD1","CLDN3","CLDN4","CLDN5","CLDN7","OCLN")
aldh <- c("ALDH1A1","ALDH1A2","ALDH1A3")
panel59 <- unique(c(pam50, claudin_low, aldh))

logCPM_59 <- logCPM_df %>% filter(Symbol %in% panel59)

# --- PAM clustering (k = 5) ---------------------------------------------------

# Transpose to samples x genes
mat_59 <- logCPM_59 %>%
  column_to_rownames("Symbol") %>%
  t()

set.seed(5896)
pam_fit <- cluster::pam(mat_59, k = 5)

clusters <- tibble::tibble(
  sample = rownames(mat_59),
  cluster = as.integer(pam_fit$clustering)
)

# --- write outputs ------------------------------------------------------------

readr::write_tsv(logCPM_df, file.path(out_dir_processed, "logCPM_all.tsv"))
readr::write_tsv(logCPM_59, file.path(out_dir_processed, "logCPM_59genes.tsv"))
readr::write_tsv(clusters,  file.path(out_dir_results, "tcga_pam_clusters.tsv"))

message("Done. Files written to:")
message("  - ", file.path(out_dir_processed, "logCPM_all.tsv"))
message("  - ", file.path(out_dir_processed, "logCPM_59genes.tsv"))
message("  - ", file.path(out_dir_results, "tcga_pam_clusters.tsv"))