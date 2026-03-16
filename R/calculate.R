# Core calculation functions

#' Calculate spatially aware gene set scores
#'
#' Uses Gaussian kernel smoothing to compute spatially aware gene set expression scores
#'
#' @param expression_matrix Expression matrix, genes×cells
#' @param gene_set Gene set vector
#' @param coords Coordinate matrix, n×2
#' @param k Number of nearest neighbors, default 30
#' @param raw_method Raw score calculation method: "mean" (default) or "addmodulescore"
#' @param spatial_weight Spatial smoothing weight, between 0-1, default 0.5
#' @param bandwidth Gaussian kernel bandwidth, automatically estimated if NULL
#' @param addmodulescore_params AddModuleScore parameter list
#' @param verbose Whether to display detailed information, default TRUE
#' @return List containing raw scores, spatial scores, and related information
#' @export
CalculateSEMA <- function(expression_matrix, gene_set, coords,
                         k = 30,
                         raw_method = "mean",
                         spatial_weight = 0.5,
                         bandwidth = NULL,
                         addmodulescore_params = list(),
                         verbose = TRUE) {
  
  if (verbose) {
    message("Starting SEMA score calculation...")
    message(sprintf("Number of cells: %d", ncol(expression_matrix)))
    message(sprintf("Number of genes: %d", nrow(expression_matrix)))
  }
  
  # Step 1: Calculate raw gene set scores
  if (verbose) message("Step 1: Calculating raw gene set scores...")
  
  if (raw_method == "addmodulescore") {
    # Using AddModuleScore method
    addmodulescore_params$features <- list(module = gene_set)
    addmodulescore_params$expression_matrix <- expression_matrix
    
    raw_result <- do.call(AddModuleScoreSEMA, addmodulescore_params)
    raw_scores <- as.numeric(raw_result$scores[, 1])
    control_genes <- raw_result$control_genes[[1]]
    
    raw_info <- list(
      method = "addmodulescore",
      control_genes = control_genes,
      n_genes = length(gene_set),
      genes_used = gene_set
    )
    
  } else if (raw_method == "mean") {
    # Using simple mean method
    # Check if genes are in expression matrix
    missing_genes <- setdiff(gene_set, rownames(expression_matrix))
    if (length(missing_genes) > 0) {
      if (verbose) {
        warning(sprintf("%d genes are not in expression matrix", length(missing_genes)))
      }
      gene_set <- intersect(gene_set, rownames(expression_matrix))
    }
    
    if (length(gene_set) == 0) {
      stop("No available genes in gene set")
    }
    
    # Calculate mean expression
    expr_subset <- expression_matrix[gene_set, , drop = FALSE]
    raw_scores <- Matrix::colMeans(expr_subset, na.rm = TRUE)
    
    raw_info <- list(
      method = "mean",
      n_genes = length(gene_set),
      genes_used = gene_set
    )
    
  } else {
    stop("raw_method must be either 'mean' or 'addmodulescore'")
  }
  
  # Step 2: Create adjacency matrix
  if (verbose) message("Step 2: Creating spatial adjacency graph...")
  graph_result <- CreateSpatialGraph(coords = coords, k = k)
  adj_matrix <- graph_result$adjacency
  
  # Step 3: Estimate bandwidth (if needed)
  if (is.null(bandwidth)) {
    if (verbose) message("Step 3: Estimating bandwidth parameter...")
    bandwidth <- EstimateBandwidth(coords, method = "median_distance", k = min(10, k))
  }
  
  # Step 4: Apply Gaussian kernel smoothing
  if (verbose) message("Step 4: Applying Gaussian kernel smoothing...")
  
  # Convert adjacency matrix to sparse format
  adj_mat <- as(adj_matrix, "TsparseMatrix")
  
  # Calculate Gaussian kernel weights using Rcpp
  kernel_result <- sparse_gaussian_kernel_cpp(
    i = adj_mat@i + 1,
    j = adj_mat@j + 1,
    coords = as.matrix(coords),
    bandwidth = bandwidth
  )
  
  weights <- kernel_result$weights
  
  # Build weighted adjacency matrix
  weighted_adj <- Matrix::sparseMatrix(
    i = adj_mat@i + 1,
    j = adj_mat@j + 1,
    x = weights,
    dims = dim(adj_matrix),
    symmetric = TRUE
  )
  
  # Normalize weight matrix (row sums = 1)
  row_sums <- Matrix::rowSums(weighted_adj)
  row_sums[row_sums == 0] <- 1
  weighted_adj_norm <- weighted_adj / row_sums
  
  # Apply spatial smoothing
  neighbor_scores <- as.numeric(weighted_adj_norm %*% raw_scores)
  spatial_scores <- (1 - spatial_weight) * raw_scores + spatial_weight * neighbor_scores
  
  # Step 5: Prepare return results
  if (verbose) message("Step 5: Organizing results...")
  
  result <- list(
    raw_scores = as.numeric(raw_scores),
    spatial_scores = as.numeric(spatial_scores),
    neighbor_scores = as.numeric(neighbor_scores),
    parameters = list(
      k = k,
      raw_method = raw_method,
      spatial_weight = spatial_weight,
      bandwidth = bandwidth
    ),
    graph = graph_result$graph,
    adjacency_matrix = adj_matrix,
    weighted_adjacency = weighted_adj_norm,
    coords = coords,
    raw_info = raw_info
  )
  
  class(result) <- "SEMAresult"
  
  if (verbose) message("Calculation complete!")
  
  return(result)
}

#' Calculate multiple gene sets in batch
#'
#' Calculate spatial scores for multiple gene sets simultaneously
#'
#' @param expression_matrix Expression matrix
#' @param gene_sets Gene set list
#' @param coords Coordinate matrix
#' @param k Number of nearest neighbors, default 30
#' @param raw_method Raw score calculation method
#' @param spatial_weight Spatial smoothing weight
#' @param bandwidth Gaussian kernel bandwidth
#' @param addmodulescore_params AddModuleScore parameters
#' @param verbose Whether to display detailed information
#' @param n_cores Number of parallel cores, default 1 (serial)
#' @return List containing all gene set results
#' @export
CalculateMultipleSEMAs <- function(expression_matrix, gene_sets, coords,
                                  k = 30,
                                  raw_method = "mean",
                                  spatial_weight = 0.5,
                                  bandwidth = NULL,
                                  addmodulescore_params = list(),
                                  verbose = TRUE,
                                  n_cores = 1) {
  
  if (verbose) {
    message(sprintf("Calculating %d gene sets...", length(gene_sets)))
  }
  
  # Precompute adjacency graph and bandwidth (shared by all gene sets)
  if (verbose) message("Precomputing spatial adjacency graph...")
  graph_result <- CreateSpatialGraph(coords = coords, k = k)
  
  if (is.null(bandwidth)) {
    bandwidth <- EstimateBandwidth(coords, method = "median_distance", k = min(10, k))
  }
  
  results <- list()
  
  for (set_name in names(gene_sets)) {
    if (verbose) message(sprintf("Calculating gene set: %s", set_name))
    
    result <- CalculateSEMA(
      expression_matrix = expression_matrix,
      gene_set = gene_sets[[set_name]],
      coords = coords,
      k = k,
      raw_method = raw_method,
      spatial_weight = spatial_weight,
      bandwidth = bandwidth,
      addmodulescore_params = addmodulescore_params,
      verbose = FALSE
    )
    
    results[[set_name]] <- result
  }
  
  class(results) <- "SEMAresults"
  
  # Add summary information
  attr(results, "summary") <- list(
    n_sets = length(gene_sets),
    n_cells = ncol(expression_matrix),
    parameters = list(
      k = k,
      raw_method = raw_method,
      spatial_weight = spatial_weight,
      bandwidth = bandwidth
    )
  )
  
  return(results)
}