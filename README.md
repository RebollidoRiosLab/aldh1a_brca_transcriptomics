# ALDH1A BRCA Transcriptomics

This repository contains scripts and resources for the analysis of publicly available bulk and single-cell transcriptomic data from breast cancer cell lines and patient tumors. The analyses focus on investigating molecular subtypes and the expression patterns of key genes, in particular the **ALDH1A subfamily: ALDH1A1, ALDH1A2, and ALDH1A3**. 

This work is associated with the publication:  

**Inter- and intra-tumoral ALDH1 heterogeneity in breast cancer identifies therapeutic opportunities for ALDH1A-specific inhibitors.**  
Raquel Pequerul, Andrada Constantinescu, Bassam Janji, Akinchan Kumar, Céline Baier, Iris Manosalva, Xavier Parés, Oscar Palacios, Salvatore Spicuglia, Delphine Colignon, Axelle Berrou, Guy Fournet, Paul Berchard, Guillaume Martin, Ismail Ceylan, Rocio Rebollido-Rios*, Jaume Farrés*, Mileidys Perez-Alea*  

---

## 📂 Repository Structure

```markdown
aldh1a_brca_transcriptomics/
├── scripts/                  # R scripts for bulk and single-cell RNA-seq analyses
│   ├── bulk_brca_depmap.R        # Bulk RNA-seq analysis of DepMap/CCLE breast cancer cell lines
│   ├── bulk_brca_tcga.R          # Bulk RNA-seq analysis of TCGA breast cancer cohorts
│   ├── scrna_brca_cell_lines.R   # scRNA-seq analysis of breast cancer cell lines
│   └── scrna_brca_patients.R     # scRNA-seq analysis of breast cancer patient tumors
├── data/
│   ├── raw/              # Raw input data (compressed when large)
│   │   ├── depmap/       # DepMap/CCLE bulk RNA-seq data
│   │   ├── tcga/         # TCGA BRCA bulk RNA-seq data
│   │   ├── sc_cell.lines/# Single-cell breast cancer cell line data
│   │   └── sc_patients/  # Single-cell breast cancer patient data
│   ├── metadata/         # Metadata tables (sample annotations, cluster info, etc.)
│   └── annotations/      # Gene annotations (e.g., GENCODE references)
└── README.md
```

---

## 📊 Data Sources

### 1. Bulk RNA-seq (Breast Cancer Cell Lines)
- **DepMap / CCLE 2019** dataset via the [DepMap portal](https://depmap.org/portal/).  
- Includes raw read counts for 50 breast cancer cell lines.  
- Molecular subtypes: basal A, basal B, HER2-enriched, and luminal.  

### 2. Bulk RNA-seq (Breast Tumor Patients)
- Publicly available mRNA expression data retrieved from [cBioPortal](https://www.cbioportal.org):  
  - **Breast Invasive Carcinoma (TCGA, Cell 2015)**  
  - **Breast Invasive Carcinoma (TCGA, PanCancer Atlas)**  
  - **Breast Invasive Carcinoma (TCGA, Firehose Legacy)**  
- 1103 unique patient samples analyzed.  

### 3. Single-cell RNA-seq (scRNA-seq)
- **Cell lines**: [Figshare dataset](https://figshare.com/articles/dataset/Single_Cell_Breast_Cancer_cell-line_Atlas/15022698) (32 breast cancer cell lines).  
- **Patient samples**: [Broad Institute Single Cell Portal](https://singlecell.broadinstitute.org/single_cell/study/SCP103930) (26 breast cancer patients).  

---

## ⚙️ Usage

### 1. Clone this repository
```bash
git clone https://github.com/RebollidoRiosLab/aldh1a_brca_transcriptomics.git
cd aldh1a_brca_transcriptomics
```  

### 2. Install required R packages
```r
install.packages(c("dplyr", "ggplot2", "readr"))
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("Seurat", "edgeR", "limma", "clusterProfiler", "org.Hs.eg.db"))
```   

### 3. Run analyses
```r
Rscript scripts/bulk_brca_tcga.R
Rscript scripts/bulk_brca_depmap.R
Rscript scripts/scrna_brca_cell_lines.R
Rscript scripts/scrna_brca_patients.R
``` 

---

## 📖 Citation

If you use this repository, please cite the associated publications:

	1.	Main paper (in proofs, corresponding author)
Pequerul R, Constantinescu A, Janji B, Kumar A, Baier C, Manosalva I, Parés X, Palacios O, Spicuglia S, Colignon D, Berrou A, Fournet G, Berchard P, Martin G, Ceylan I, Rebollido-Ríos R*, Farrés J*, Perez-Alea M*.
Inter- and intra-tumoral ALDH1 heterogeneity in breast cancer identifies therapeutic opportunities for ALDH1A-specific inhibitors.
	
 2.	Related work
Pequerul R, Constantinescu A, Janji B, Kumar A, Parés X, Palacios O, Colignon D, Berrou A, Fournet G, Berchard P, Martin G, Ceylan I, Rebollido-Ríos R, Farrés J, Perez-Alea M.
ALDH1 subtype-specific inhibitor targets key cancerous epithelial cell populations in aggressive subtypes of breast cancer.
bioRxiv 2024.10.18.619128. doi: https://doi.org/10.1101/2024.10.18.619128
	
 3.	Code & data DOI (Zenodo, to be added)
Will be updated once Zenodo DOI is minted.


⸻

## 📬 Contact
Rocio Rebollido-Rios, PhD (Corresponding Author) -> rocio.rebollido-rios@uni-koeln.de
