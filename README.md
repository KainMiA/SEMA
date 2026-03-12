# SEMA: A Spatial-Aware Framework for Multi-Pathway Enrichment Analysis

### 📖 Overview

SMEA (Spatial-Aware Framework for Multi-Pathway Enrichment Analysis) is an R package designed for comprehensive gene set enrichment analysis in spatial transcriptomics data. It integrates spatial information directly into gene set scoring, enabling the discovery of spatially-informed biological patterns that traditional methods might miss.

<div align="left"> <img src="https://github.com/KainMiA/SEMA/blob/main/SEMA_diagram.png" alt="SEMA" width="800"/> </div>


### 📦 Installation

From GitHub

    # Install from GitHub
    devtools::install_github("Kain/SMEA/SEMA")

From Source

    # Or install from local source
    devtools::install_local("path/to/SMEA")


### 🚀 Quick Start

Basic Usage with Seurat

    library(SMEA)
    library(Seurat)
    
    # Load your Seurat object
    seurat_obj <- readRDS("your_data.rds")
    
    # Define gene sets
    gene_sets <- list(
      "Immune_Response" = c("CD8A", "CD4", "FOXP3", "PDCD1", "CTLA4"),
      "Angiogenesis" = c("VEGFA", "KDR", "PECAM1", "CD34"),
      "Hypoxia" = c("HIF1A", "VEGFA", "SLC2A1", "PGK1")
    )
    
    # Run spatial enrichment analysis
    smea_results <- RunSMEA(
      object = seurat_obj,
      gene_sets = gene_sets,
      assay = "Spatial",
      method = "mean",
      k = 6,
      spatial.weight = 0.5
    )
    
    # View results
    print(smea_results)
    summary(smea_results)
    

