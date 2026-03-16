# Main function - compatible with Seurat objects

#' Run SEMA analysis
#'
#' Spatially aware gene set expression analysis, supports Seurat object input
#'
#' @param object Seurat object or expression matrix
#' @param gene_sets Gene set list or single gene set vector
#' @param assay Assay name to use (for Seurat objects only), default "RNA"
#' @param slot Slot to use (for Seurat objects only), default "data"
#' @param coords Coordinate matrix, if Seurat object and not provided, obtained from reductions$spatial
#' @param reduction Dimensionality reduction method to obtain coordinates (for Seurat objects only), default "spatial"
#' @param k Number of nearest neighbors, default 30
#' @param raw_method Raw score calculation method: "mean" (default) or "addmodulescore"
#' @param spatial_weight Spatial smoothing weight, between 0-1, default 0.5
#' @param bandwidth Gaussian kernel bandwidth, automatically estimated if NULL
#' @param addmodulescore_params AddModuleScore parameter list
#' @param verbose Whether to display detailed information, default TRUE
#' @param n_cores Number of parallel cores, default 1
#' @param ... Other parameters passed to CalculateSEMA or CalculateMultipleSEMAs
#'
#' @return SEMA result object
#'
#' @examples
#' \dontrun{
#' # Using Seurat object
#' library(Seurat)
#' data("pbmc_small")
#'
#' # Create simulated spatial coordinates
#' coords <- matrix(runif(ncol(pbmc_small)*2), ncol=2)
#' rownames(coords) <- colnames(pbmc_small)
#'
#' # Add coordinates to Seurat object
#' pbmc_small[["spatial"]] <- CreateDimReducObject(
#'   embeddings = coords,
#'   key = "spatial_",
#'   assay = "RNA"
#' )
#'
#' # Define gene sets
#' gene_sets <- list(
#'   cell_cycle = c("MCM5", "PCNA", "TYMS"),
#'   stress_response = c("HSPA1A", "HSPA1B", "DNAJB1")
#' )
#'
#' # Run SEMA
#' result <- RunSEMA(pbmc_small, gene_sets = gene_sets)
#' }
#'
#' @export
RunSEMA <- function(object, gene_sets,
                    assay = NULL,
                    slot = "data",
                    coords = NULL,
                    reduction = "spatial",
                    k = 30,
                    raw_method = "mean",
                    spatial_weight = 0.5,
                    bandwidth = NULL,
                    addmodulescore_params = list(),
                    verbose = TRUE,
                    n_cores = 1,
                    ...) {
  
  UseMethod("RunSEMA", object)
}

#' @method RunSEMA Seurat
#' @export
RunSEMA.Seurat <- function(object, gene_sets,
                          assay = NULL,
                          slot = "data",
                          coords = NULL,
                          reduction = "spatial",
                          k = 30,
                          raw_method = "mean",
                          spatial_weight = 0.5,
                          bandwidth = NULL,
                          addmodulescore_params = list(),
                          verbose = TRUE,
                          n_cores = 1,
                          ...) {
  
  if (verbose) message("Input object type: Seurat")
  
  # Get assay
  if (is.null(assay)) {
    assay <- Seurat::DefaultAssay(object)
    if (verbose) message(sprintf("Using default assay: %s", assay))
  }
  
  # Get expression matrix
  if (verbose) message("Extracting expression matrix...")
  expression_matrix <- Seurat::GetAssayData(object, assay = assay, slot = slot)
  
  # Get coordinates
  if (is.null(coords)) {
    if (reduction %in% names(object@reductions)) {
      if (verbose) message(sprintf("Extracting coordinates from reduction '%s'", reduction))
      coords <- Seurat::Embeddings(object, reduction = reduction)
    } else {
      stop(sprintf("Reduction '%s' not found", reduction))
    }
  }
  
  # Check coordinate dimensions
  if (ncol(coords) > 2) {
    warning(sprintf("Coordinates have %d dimensions, using only the first 2", ncol(coords)))
    coords <- coords[, 1:2]
  }
  
  # Ensure consistent cell order
  common_cells <- intersect(colnames(expression_matrix), rownames(coords))
  if (length(common_cells) == 0) {
    stop("Cell names do not match between expression matrix and coordinates")
  }
  
  if (length(common_cells) < ncol(expression_matrix)) {
    warning(sprintf("Using only %d common cells", length(common_cells)))
    expression_matrix <- expression_matrix[, common_cells]
    coords <- coords[common_cells, ]
  }
  
  # Run SEMA analysis
  if (verbose) message(sprintf("Number of cells: %d", ncol(expression_matrix)))
  if (verbose) message(sprintf("Number of genes: %d", nrow(expression_matrix)))
  
  if (is.list(gene_sets) && length(gene_sets) > 1) {
    # Multiple gene sets
    result <- CalculateMultipleSEMAs(
      expression_matrix = expression_matrix,
      gene_sets = gene_sets,
      coords = coords,
      k = k,
      raw_method = raw_method,
      spatial_weight = spatial_weight,
      bandwidth = bandwidth,
      addmodulescore_params = addmodulescore_params,
      verbose = verbose,
      n_cores = n_cores,
      ...
    )
  } else {
    # Single gene set
    result <- CalculateSEMA(
      expression_matrix = expression_matrix,
      gene_set = if (is.list(gene_sets)) gene_sets[[1]] else gene_sets,
      coords = coords,
      k = k,
      raw_method = raw_method,
      spatial_weight = spatial_weight,
      bandwidth = bandwidth,
      addmodulescore_params = addmodulescore_params,
      verbose = verbose,
      ...
    )
  }
  
  # Add Seurat object information
  attr(result, "seurat_info") <- list(
    assay = assay,
    slot = slot,
    reduction = reduction,
    n_cells = ncol(expression_matrix)
  )
  
  return(result)
}

#' @method RunSEMA default
#' @export
RunSEMA.default <- function(object, gene_sets,
                           assay = NULL,
                           slot = NULL,
                           coords = NULL,
                           reduction = NULL,
                           k = 30,
                           raw_method = "mean",
                           spatial_weight = 0.5,
                           bandwidth = NULL,
                           addmodulescore_params = list(),
                           verbose = TRUE,
                           n_cores = 1,
                           ...) {
  
  if (verbose) message("Input object type: expression matrix")
  
  # Directly use provided parameters
  expression_matrix <- object
  
  if (is.null(coords)) {
    stop("For non-Seurat object input, coords parameter must be provided")
  }
  
  # Check dimensions
  if (ncol(expression_matrix) != nrow(coords)) {
    stop(sprintf("Cell count mismatch: expression matrix has %d cells, coordinates have %d cells",
                 ncol(expression_matrix), nrow(coords)))
  }
  
  # Run SEMA analysis
  if (is.list(gene_sets) && length(gene_sets) > 1) {
    # Multiple gene sets
    result <- CalculateMultipleSEMAs(
      expression_matrix = expression_matrix,
      gene_sets = gene_sets,
      coords = coords,
      k = k,
      raw_method = raw_method,
      spatial_weight = spatial_weight,
      bandwidth = bandwidth,
      addmodulescore_params = addmodulescore_params,
      verbose = verbose,
      n_cores = n_cores,
      ...
    )
  } else {
    # Single gene set
    result <- CalculateSEMA(
      expression_matrix = expression_matrix,
      gene_set = if (is.list(gene_sets)) gene_sets[[1]] else gene_sets,
      coords = coords,
      k = k,
      raw_method = raw_method,
      spatial_weight = spatial_weight,
      bandwidth = bandwidth,
      addmodulescore_params = addmodulescore_params,
      verbose = verbose,
      ...
    )
  }
  
  return(result)
}

#' Add SEMA results to Seurat object
#'
#' Add SEMA-calculated spatial scores to Seurat object metadata
#'
#' @param seurat_obj Seurat object
#' @param sema_result SEMA result object
#' @param prefix Prefix to add to metadata column names, default "SEMA_"
#' @param use_spatial_scores Whether to use spatial scores (TRUE) or raw scores (FALSE), default TRUE
#'
#' @return Updated Seurat object
#' @export
AddSEMAtoSeurat <- function(seurat_obj, sema_result, prefix = "SEMA_", use_spatial_scores = TRUE) {
  
  if (!inherits(sema_result, c("SEMAresult", "SEMAresults"))) {
    stop("sema_result must be SEMAresult or SEMAresults object")
  }
  
  if (inherits(sema_result, "SEMAresult")) {
    # Single gene set
    if (use_spatial_scores) {
      scores <- sema_result$spatial_scores
    } else {
      scores <- sema_result$raw_scores
    }
    
    score_name <- if (!is.null(sema_result$raw_info$genes_used)) {
      paste("gene_set", length(sema_result$raw_info$genes_used), sep = "_")
    } else {
      "gene_set"
    }
    
    col_name <- paste0(prefix, score_name)
    seurat_obj[[col_name]] <- scores
    
  } else if (inherits(sema_result, "SEMAresults")) {
    # Multiple gene sets
    for (set_name in names(sema_result)) {
      if (use_spatial_scores) {
        scores <- sema_result[[set_name]]$spatial_scores
      } else {
        scores <- sema_result[[set_name]]$raw_scores
      }
      
      col_name <- paste0(prefix, set_name)
      seurat_obj[[col_name]] <- scores
    }
  }
  
  return(seurat_obj)
}