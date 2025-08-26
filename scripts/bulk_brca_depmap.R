#!/usr/bin/env Rscript
# ------------------------------------------------------------------------------
# bulk_brca_depmap.R
# Process DepMap/CCLE bulk RNA-seq for breast cancer cell lines:
# - load counts + metadata
# - TMM normalization (edgeR) -> logCPM
# - subset to 59-gene panel (PAM50 + CHD1, CLDN3/4/5/7, OCLN, ALDH1A1/2/3)
# - PAM clustering (k = 4)
# - write processed matrices + cluster assignments
#
# Inputs (expected relative paths):
#   data/raw/depmap/CCLE_RNAseq_genes_counts_20180929.gct.gz
#   data/raw/depmap/sample_info22Q2.csv
#   data/raw/annotations/gencode.v19/geneInfo.tab
#
# Outputs:
#   data/processed/depmap/logCPM_all.tsv
#   data/processed/depmap/logCPM_59genes.tsv
#   results/tables/depmap_pam_clusters.tsv
#
# Notes:
# - This script focuses on bulk RNA-seq only (no proteomics, RPPA, miRNA).
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(cluster)       # PAM
})

# --- helpers ------------------------------------------------------------------

ensure_dirs <- function(paths) {
  for (p in paths) if (!dir.exists(p)) dir.create(p, recursive = TRUE)
}

read_gct_counts <- function(path_gct_gz) {
  # Minimal GCT reader for counts: skip first 2 header lines
  # Expected first two columns are identifiers (use them), rest are samples
  df <- readr::read_tsv(path_gct_gz, skip = 2, show_col_types = FALSE)
  # Try to standardize column names if typical GCT fields are present
  # If GCT has columns: Name, Description, then counts start at col 3
  if (all(c("Name","Description") %in% names(df))) {
    df <- df %>% rename(Ensembl = Name)
  } else if ("Name" %in% names(df)) {
    df <- df %>% rename(Ensembl = Name)
  } else if (!("Ensembl" %in% names(df))) {
    # Fallback: assume first column is Ensembl-like identifier
    names(df)[1] <- "Ensembl"
  }
  as.data.frame(df)
}

# --- I/O paths ----------------------------------------------------------------

in_counts <- "data/raw/depmap/CCLE_RNAseq_genes_counts_20180929.gct.gz"
in_sampleinfo <- "data/raw/depmap/sample_info22Q2.csv"
in_geneinfo <- "data/raw/annotations/geneInfo.v19.tab"

out_dir_processed <- "data/processed/depmap"
out_dir_results <- "results/tables"

ensure_dirs(c(out_dir_processed, out_dir_results))

# --- load data ----------------------------------------------------------------

message("Reading DepMap/CCLE counts...")
counts_df <- read_gct_counts(in_counts)

message("Reading DepMap sample info...")
info <- readr::read_csv(in_sampleinfo, show_col_types = FALSE)

# Keep only breast lines; DepMap marks lineage in `lineage` or use name filter
info_breast <- info %>%
  filter(tolower(lineage) == "breast" | grepl("BREAST", CCLE_Name, ignore.case = TRUE))

# Sample columns in counts sometimes have suffices like "_RNAseq"; strip suffix after underscore
sample_cols <- setdiff(names(counts_df), c("Ensembl","Description"))
clean_names <- gsub("_.*", "", sample_cols)
colnames(counts_df)[match(sample_cols, names(counts_df))] <- clean_names

# Match cell lines present in both metadata and counts
keep_cells <- intersect(info_breast$stripped_cell_line_name, clean_names)
stopifnot(length(keep_cells) > 0)

counts_df <- counts_df[, c("Ensembl", intersect(names(counts_df), keep_cells))]

# Map Ensembl to HGNC symbols (Gencode v19 geneInfo.tab: Ensembl, Symbol, Biotype)
gene_info <- readr::read_tsv(in_geneinfo, col_names = c("Ensembl","Symbol","Biotype"),
                             skip = 1, show_col_types = FALSE)

counts_df <- counts_df %>%
  inner_join(gene_info, by = "Ensembl") %>%
  relocate(Symbol, Biotype, .after = Ensembl) %>%
  arrange(Symbol)

# Optional: keep protein_coding only (comment out if you prefer all)
# counts_df <- counts_df %>% filter(Biotype == "protein_coding")

# --- TMM normalization -> logCPM ----------------------------------------------

mat_counts <- counts_df %>%
  select(all_of(keep_cells)) %>%
  as.matrix()
rownames(mat_counts) <- counts_df$Symbol

dge <- DGEList(counts = mat_counts)
keep <- filterByExpr(dge, min.count = 10)
dge  <- dge[keep, , keep.lib.sizes = FALSE]
dge  <- calcNormFactors(dge)
logCPM <- edgeR::cpm(dge, log = TRUE, prior.count = 1)  # samples in columns

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
# Some symbols (e.g., ORC6) may be ORC6L in older refs; ensure we use ORC6 per your script
logCPM_59$Symbol[logCPM_59$Symbol == "ORC6L"] <- "ORC6"

# --- PAM clustering (k = 4) ---------------------------------------------------

# Transpose to samples x genes
mat_59 <- logCPM_59 %>%
  column_to_rownames("Symbol") %>%
  t()

set.seed(5896)
pam_fit <- cluster::pam(mat_59, k = 4)

clusters <- tibble::tibble(
  sample = rownames(mat_59),
  cluster = as.integer(pam_fit$clustering)
)

# --- write outputs ------------------------------------------------------------

readr::write_tsv(logCPM_df, file.path(out_dir_processed, "logCPM_all.tsv"))
readr::write_tsv(logCPM_59, file.path(out_dir_processed, "logCPM_59genes.tsv"))
readr::write_tsv(clusters,  file.path(out_dir_results, "depmap_pam_clusters.tsv"))

message("Done. Files written to:")
message("  - ", file.path(out_dir_processed, "logCPM_all.tsv"))
message("  - ", file.path(out_dir_processed, "logCPM_59genes.tsv"))
message("  - ", file.path(out_dir_results, "depmap_pam_clusters.tsv"))