# SEMA: A Spatial-Aware Framework for Multi-Pathway Enrichment Analysis

SEMA is an R package designed for comprehensive gene set enrichment analysis in spatial transcriptomics data. It integrates spatial information directly into gene set scoring, enabling the discovery of spatially informed biological patterns that traditional methods might miss.

---

## Features

- **Gene-set enrichment** — score any list of gene sets on a spot/cell level.
- **Spatially aware** — enrichment scores are smoothed over `k` spatial
  neighbors, controlled by `spatial.weight`.

---

## Installation

```r
# install.packages("devtools")
devtools::install_github("your-org/SEMA")
```

```r
library(SEMA)
library(Seurat)
```

**Requirements:** R >= 4.1, Seurat >= 5.0.

---

## Quick Start

### 1. Load a Seurat object

Your object must contain a spatial assay (e.g. `Spatial` from Visium, or
`Xenium`/`slide-seq` equivalents) with coordinates available.

```r
seurat_obj <- readRDS("demo_seurat_data.rds")
```

### 2. Define gene sets

```r
gene_sets <- list(
  Immune_Response = c("CD8A", "CD4", "FOXP3", "PDCD1", "CTLA4"),
  Angiogenesis    = c("VEGFA", "KDR", "PECAM1", "CD34"),
  Hypoxia         = c("HIF1A", "VEGFA", "SLC2A1", "PGK1")
)
```

### 3. Run spatial enrichment analysis

```r
sema_results <- RunSEMA(
  object         = seurat_obj,
  gene_sets      = gene_sets,
  assay          = "Spatial",
  method         = "mean",
  k              = 6,
  spatial.weight = 0.5
)
```

### 4. Add enrichment scores to the Seurat object

```r
for (gs in names(sema_results)) {
  seurat_obj[[gs]] <- sema_results[[gs]]$spatial_scores
}
```

### 5. Visualize on the tissue

```r
SpatialFeaturePlot(
  seurat_obj,
  features = c("Immune_Response", "Angiogenesis", "Hypoxia"),
  ncol     = 3
)
```

---

## Parameters

| Argument | Description |
| --- | --- |
| `object` | A `Seurat` object containing spatial coordinates. |
| `gene_sets` | Named `list` of character vectors; names become metadata columns. |
| `assay` | Assay to pull expression from, e.g. `"Spatial"`. |
| `method` | Aggregation method used to collapse a gene set into one score (e.g. `"mean"`). |
| `k` | Number of spatial neighbors used for smoothing. |
| `spatial.weight` | Weight of the spatial smoothing term; usually between `0` (expression only) and `1` (spatial only). |

---

## Output

`RunSEMA()` returns a named `list` with one entry per gene set. Each entry
contains:

- `spatial_scores` — numeric vector of spatially smoothed enrichment scores,
  one value per spot, aligned to the columns of `object`.

Scores can be inspected directly:

```r
head(sema_results$Immune_Response$spatial_scores)
```

or plotted as a distribution:

```r
hist(sema_results$Hypoxia$spatial_scores,
     main = "Hypoxia spatial enrichment", xlab = "score")
```

---

## Tips

- **Choose `k` to match your resolution.** A larger `k` yields smoother maps
  but blurs fine-grained structure; for Visium, `k = 6` is a reasonable
  starting point.
- **Tune `spatial.weight`.** Increase it when the signal is noisy and you care
  about regional patterns; decrease it to preserve spot-level heterogeneity.
- **Low-quality spots** with few detected genes can produce unstable scores —
  consider filtering before running SEMA.

---
For a detailed tutorial, please visit: https://kainmia.github.io/SEMA/
