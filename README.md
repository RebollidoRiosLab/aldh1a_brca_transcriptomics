# BRCA Transcriptomics

This repository contains scripts and resources for the analysis of publicly available bulk and single-cell transcriptomic data from breast cancer cell lines and patient tumors. The analyses focus on investigating molecular subtypes and the expression patterns of key genes, including the **ALDH1A subfamily: ALDH1A1, ALDH1A2, and ALDH1A3 genes**.

---

## Data Sources

### 1. Breast Cancer Cell Lines (Bulk RNA-seq)
- Publicly available data retrieved from the **CCLE 2019** dataset via the [DepMap portal](https://depmap.org/portal/).
- Includes raw read counts for 50 human breast cancer cell lines.
- Molecular subtypes classified into basal A, basal B, HER2-enriched, and luminal.

### 2. Breast Tumor Patient Samples (Bulk RNA-seq)
- Publicly available mRNA gene expression data retrieved from [cBioPortal](https://www.cbioportal.org):
  - **Breast Invasive Carcinoma (TCGA, Cell 2015)**
  - **Breast Invasive Carcinoma (TCGA, PanCancer Atlas)**
  - **Breast Invasive Carcinoma (TCGA, Firehose Legacy)**
- A total of 1103 unique patient samples analyzed.

### 3. Single-cell RNA-seq (scRNA-seq)
- **Cell lines**: Publicly available data retrieved from [Figshare](https://figshare.com/articles/dataset/Single_Cell_Breast_Cancer_cell-line_Atlas/15022698), encompassing 32 breast cancer cell lines.
- **Patient samples**: Publicly available data retrieved from the [Broad Institute Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP103930), encompassing 26 individual patient samples.

---

This repository provides the necessary scripts to process, analyze, and visualize these datasets for further exploration of breast cancer molecular subtypes and gene expression profiles.
