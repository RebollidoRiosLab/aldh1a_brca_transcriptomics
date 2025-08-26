#!/usr/bin/env Rscript
# =============================================================================
# sc_brca_cell_lines.R
# Single-cell RNA-seq analysis of breast cancer cell lines (Figshare atlas)
# Paper scope: TNBC & HER2+ subsets, PAM50/ALDH panels, SCT v2 + glmGamPoi
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(readr)
  library(stringr)
  library(dplyr)
  library(tibble)
  library(Seurat)
  library(sctransform)
  library(glmGamPoi)
  library(ggplot2)
  library(cowplot)
  library(openxlsx)
})

theme_set(theme_cowplot())
set.seed(333)

# ------------------------------- Paths ---------------------------------------
# Place the three Figshare files here:
# matrix.mtx.gz, features.tsv.gz, barcodes.tsv.gz
fig_dir    <- here("data", "raw", "sc_cell.lines")
mtx_path   <- file.path(fig_dir, "matrix.mtx.gz")
feat_path  <- file.path(fig_dir, "features.tsv.gz")
barc_path  <- file.path(fig_dir, "barcodes.tsv.gz")

# GENCODE gene info (used v28); adjust if neccesary 
gencode_tab <- here("data", "annotations", "gencode_v28_geneInfo.tsv")

# Optional: CCLE / DepMap cell-line annotations
cellline_meta <- here("data", "metadata", "sc_cell.lines_metadata.csv")

# Outputs
out_dir   <- here("outputs", "sc_cell_lines")
fig_dir_o <- file.path(out_dir, "figures")
tab_dir_o <- file.path(out_dir, "tables")
dir.create(fig_dir_o, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir_o, recursive = TRUE, showWarnings = FALSE)

# ----------------------------- Parameters ------------------------------------
remove_gene <- "MALAT1"
npcs <- 10
umap_dims <- 1:5
resolution  <- 0.6  # you used 0.6–0.8; keep 0.6 as default
k_param <- 50

pam50 <- c("ACTR3B","ANLN","BAG1","BCL2","BIRC5","BLVRA","CCNB1","CCNE1",
           "CDC20","CDC6","CDH3","CENPF","CEP55","CXXC5","EGFR","ERBB2",
           "ESR1","EXO1","FGFR4","FOXA1","FOXC1","GPR160","GRB7","KIF2C",
           "KRT14","KRT17","KRT5","MAPT","MDM2","MELK","MIA","MKI67","MLPH",
           "MMP11","MYBL2","MYC","NAT1","NDC80","NUF2","ORC6","PGR","PHGDH",
           "PTTG1","RRM2","SFRP1","SLC39A6","TMEM45B","TYMS","UBE2C","UBE2T")
claudin_low <- c("CLDN3","CLDN4","CLDN5","CLDN7","OCLN","CHD1")
aldh        <- c("ALDH1A1","ALDH1A2","ALDH1A3")
gene_panel  <- unique(c(pam50, claudin_low, aldh))

# --------------------- Read + add gene symbols to features -------------------
# Figshare features are Ensembl IDs (often with version). Map to symbols.
gencode <- read.delim(gencode_tab, header = FALSE, skip = 1, sep = "\t",
                      col.names = c("Ensembl","Symbol","Biotype"))
gencode$Ensembl <- sub("\\..*$", "", gencode$Ensembl)

features_in <- read_tsv(feat_path, col_names = FALSE, show_col_types = FALSE) |>
  rename(Ensembl = X1) |>
  mutate(Ensembl = sub("\\..*$", "", Ensembl)) |>
  left_join(gencode, by = "Ensembl") |>
  distinct(Ensembl, .keep_all = TRUE)

# If missing symbols, keep Ensembl as fallback
features_in$Symbol[is.na(features_in$Symbol)] <- features_in$Ensembl

# ------------------------------- Load counts ---------------------------------
mat <- ReadMtx(mtx = mtx_path, features = feat_path, cells = barc_path)
# Replace feature names with mapped symbols (Seurat expects rownames = features)
rownames(mat) <- features_in$Symbol[match(rownames(mat), features_in$Ensembl)]

# --------------------------- Create Seurat object ----------------------------
obj <- CreateSeuratObject(counts = mat, project = "BRCA_cell_lines")
rm(mat); gc()

# ------------------------------ Cell-line meta -------------------------------
if (file.exists(cellline_meta)) {
  meta <- read.csv(cellline_meta, stringsAsFactors = FALSE)
  # Expect a column matching obj$orig.ident (e.g., "CCLE_ID" or "orig.ident")
  key_col <- intersect(colnames(meta), c("orig.ident","CCLE_ID","Cell.Line"))[1]
  stopifnot(!is.na(key_col))
  meta <- meta |>
    rename(orig.ident = !!key_col) |>
    mutate(across(everything(), ~ tidyr::replace_na(., "NA")))
  # Attach a few useful annotations if present
  add_if <- function(col, newname) {
    if (col %in% names(meta)) {
      obj[[newname]] <- meta[[col]][match(obj$orig.ident, meta$orig.ident)]
    }
  }
  add_if("mol_subtype", "mol.subtype")
  add_if("marker_expr", "marker.expr")
  add_if("Primary.Disease", "primary.site")
  add_if("Tumor.Type", "tumor.type")
}

# ------------------------------ Quick QC -------------------------------------
obj <- PercentageFeatureSet(obj, pattern = "^RP[SL]", col.name = "percent_ribo")

pdf(file.path(fig_dir_o, "qc_violin_by_cellline.pdf"), width = 12, height = 7)
print(VlnPlot(obj, group.by = "orig.ident",
              features = c("nFeature_RNA","nCount_RNA","percent_ribo"),
              pt.size = 0.05) + NoLegend())
dev.off()

# ------------------------ Remove MALAT1 (per manuscript) ---------------------
if (remove_gene %in% rownames(obj)) {
  obj <- obj[setdiff(rownames(obj), remove_gene), ]
}

# ---------------------- Normalize + reduce (SCT v2) --------------------------
obj <- SCTransform(obj, method = "glmGamPoi", vst.flavor = "v2", verbose = TRUE)
obj <- RunPCA(obj, npcs = npcs, verbose = TRUE)
obj <- RunUMAP(obj, dims = umap_dims, umap.method = "umap-learn")
obj <- FindNeighbors(obj, dims = umap_dims, k.param = k_param)
obj <- FindClusters(obj, algorithm = 2, resolution = resolution)

# ------------------------------- Plots ---------------------------------------
pdf(file.path(fig_dir_o, "umap_overview.pdf"), width = 12, height = 5)
p1 <- DimPlot(obj, reduction = "umap", label = TRUE, repel = TRUE,
              label.size = 4, pt.size = 0.3) + NoLegend() + ggtitle("Clusters")
p2 <- DimPlot(obj, reduction = "umap", group.by = "orig.ident", label = FALSE,
              pt.size = 0.3) + NoLegend() + ggtitle("Cell lines")
print(p1 + p2)
dev.off()

if ("mol.subtype" %in% colnames(obj@meta.data)) {
  ggsave(file.path(fig_dir_o, "umap_by_molecular_subtype.pdf"),
         DimPlot(obj, reduction="umap", group.by="mol.subtype", pt.size=0.3),
         width=6.5, height=5)
}
if ("marker.expr" %in% colnames(obj@meta.data)) {
  ggsave(file.path(fig_dir_o, "umap_by_marker_expr.pdf"),
         DimPlot(obj, reduction="umap", group.by="marker.expr", pt.size=0.3),
         width=6.5, height=5)
}

# -------------------------- PAM50 / ALDH panels ------------------------------
present <- intersect(gene_panel, rownames(obj))
if (length(present)) {
  pdf(file.path(fig_dir_o, "vln_pam50_claudin_aldh.pdf"), width = 10, height = 16)
  print(VlnPlot(obj, features = present, pt.size = 0, ncol = 4) + NoLegend())
  dev.off()
  
  pdf(file.path(fig_dir_o, "feature_ALDH1A_panel.pdf"), width = 9, height = 4)
  print(FeaturePlot(obj, features = intersect(aldh, rownames(obj)),
                    cols = c("lightgrey", "firebrick2"),
                    label = FALSE, ncol = 3, pt.size = 0.3))
  dev.off()
}

# --------------------------- Save objects & tables ---------------------------
saveRDS(obj, file = file.path(out_dir, "sc_brca_cell_lines.seurat.rds"))

sessionInfo_df <- utils::capture.output(sessionInfo())
writeLines(sessionInfo_df, con = file.path(out_dir, "sessionInfo.txt"))