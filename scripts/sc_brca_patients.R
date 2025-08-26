#!/usr/bin/env Rscript
# =============================================================================
# sc_brca_patients.R
# Single-cell RNA-seq analysis of breast cancer patient samples
# Integration by molecular subtype (Basal, HER2, LumA, LumB, LumB/HER2),
# as described in Methods
# =============================================================================

suppressPackageStartupMessages({
  library(here)
  library(Seurat)
  library(sctransform)
  library(glmGamPoi)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(cowplot)
  library(openxlsx)
})

theme_set(theme_cowplot())
set.seed(333)

# ------------------------------- Paths ---------------------------------------
raw_dir <- here("data", "raw", "sc_patients")
meta_csv <- here("data", "metadata", "sc_patients_metadata.csv")

out_dir <- here("outputs", "sc_patients")
fig_dir <- file.path(out_dir, "figures")
tab_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

# ----------------------------- Parameters ------------------------------------
nfeatures_integration <- 3000
dims_use <- 1:20
npcs <- 20
resolution <- 0.1

aldh <- c("ALDH1A1","ALDH1A2","ALDH1A3")

# --------------------------- Read per-sample data ----------------------------
sample_dirs <- list.dirs(raw_dir, full.names = TRUE, recursive = FALSE)
stopifnot(length(sample_dirs) > 0)

objs <- lapply(sample_dirs, function(sd){
  mtx  <- file.path(sd, "matrix.mtx.gz")
  feat <- file.path(sd, "features.tsv.gz")
  barc <- file.path(sd, "barcodes.tsv.gz")
  mat <- ReadMtx(mtx, features = feat, cells = barc)
  sobj <- CreateSeuratObject(mat, project = basename(sd))
  sobj$sample_id <- basename(sd)
  sobj
})

# ------------------------------- Metadata ------------------------------------
if (file.exists(meta_csv)) {
  meta <- read.csv(meta_csv, stringsAsFactors = FALSE)
  for (nm in c("mol_subtype")) {
    if (nm %in% names(meta)) {
      for (i in seq_along(objs)) {
        sid <- unique(objs[[i]]$sample_id)
        lab <- meta[[nm]][match(sid, meta$sample_id)]
        objs[[i]][[nm]] <- lab
      }
    }
  }
}

# ----------------------------- Per-group SCT ---------------------------------
objs <- lapply(objs, function(x){
  x <- PercentageFeatureSet(x, pattern = "^MT-", col.name = "percent_mito")
  SCTransform(x, method = "glmGamPoi", vst.flavor = "v2",
              vars.to.regress = "percent_mito", verbose = FALSE)
})

# -------------------------- Integration by subtype ---------------------------
# Split objects by molecular subtype
stopifnot("mol_subtype" %in% colnames(objs[[1]]@meta.data))
obj_list <- SplitObject(merge(x = objs[[1]], y = objs[-1]),
                        split.by = "mol_subtype")

features <- SelectIntegrationFeatures(object.list = obj_list, nfeatures = nfeatures_integration)
obj_list <- PrepSCTIntegration(object.list = obj_list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = obj_list,
                                  normalization.method = "SCT",
                                  anchor.features = features, dims = dims_use)
obj <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# -------------------------- Reduction & clustering ---------------------------
obj <- RunPCA(obj, npcs = npcs)
obj <- RunUMAP(obj, dims = dims_use, umap.method = "umap-learn")
obj <- FindNeighbors(obj, dims = dims_use, k.param = 50)
obj <- FindClusters(obj, resolution = resolution, algorithm = 2)

# ------------------------------- Plots ---------------------------------------
pdf(file.path(fig_dir, "umap_patients_by_subtype.pdf"), width = 7, height = 5)
print(DimPlot(obj, reduction = "umap", group.by = "mol_subtype", pt.size = 0.2))
dev.off()

if (length(intersect(aldh, rownames(obj))) > 0) {
  ggsave(file.path(fig_dir, "feature_ALDH1A_panel_patients.pdf"),
         FeaturePlot(obj, features = intersect(aldh, rownames(obj)),
                     cols = c("lightgrey","firebrick2"), pt.size = 0.25, ncol = 3),
         width = 9, height = 4)
}

# --------------------------- Save outputs ------------------------------------
saveRDS(obj, file = file.path(out_dir, "sc_brca_patients.integrated.seurat.rds"))

wb <- createWorkbook()
clust_by_subtype <- as.data.frame.matrix(table(obj$seurat_clusters, obj$mol_subtype)) %>%
  rownames_to_column("cluster") %>%
  mutate(total_cells = rowSums(across(-cluster))) %>%
  relocate(total_cells, .after = cluster)
addWorksheet(wb, "cells_by_mol_subtype")
writeData(wb, "cells_by_mol_subtype", clust_by_subtype)
saveWorkbook(wb, file.path(tab_dir, "sc_patients_summary.xlsx"), overwrite = TRUE)

# --------------------------- Reproducibility ---------------------------------
sessionInfo_df <- utils::capture.output(sessionInfo())
writeLines(sessionInfo_df, con = file.path(out_dir, "sessionInfo.txt"))