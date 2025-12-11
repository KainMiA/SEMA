# SEMA: Spatial Expression Multi-pathways Enrichment Analysis


### 📖 Overview

SMEA (Spatial Expression Multi-pathways Enrichment Analysis) is an R package designed for comprehensive gene set enrichment analysis in spatial transcriptomics data. It integrates spatial information directly into gene set scoring, enabling the discovery of spatially-informed biological patterns that traditional methods might miss.

<div align="center"> <img src="https://raw.githubusercontent.com/yourusername/SMEA/main/vignettes/figures/spatial_scores.png" alt="Spatial Scores Visualization" width="800"/> </div>


### 📦 Installation

From GitHub

    # Install from GitHub
    devtools::install_github("yourusername/SMEA")

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
      gene.sets = gene_sets,
      method = "aucell",
      spatial.method = "gaussian_kernel",
      k = 10,
      spatial.weight = 0.3,
      n.cores = 4,
      verbose = TRUE
    )
    
    # View results
    print(smea_results)
    summary(smea_results)
    
    # Visualize results
    plot(smea_results, gene.set = "Immune_Response", type = "spatial")
    plot(smea_results, gene.set = "Immune_Response", type = "histogram")


### 📄 Citation

If you use SMEA in your research, please cite:

    bibtex
    @software{smea_package,
      title = {SMEA: Spatial Expression Multi-pathways Enrichment Analysis},
      author = {Your Name and Contributors},
      year = {2024},
      version = {0.1.0},
      url = {https://github.com/yourusername/SMEA},
      note = {R package}
    }

