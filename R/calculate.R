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

  coords <- .CheckCoords(coords, ncol(expression_matrix))

  # Step 1: Calculate raw gene set scores
  if (verbose) message("Step 1: Calculating raw gene set scores...")
  raw <- .CalculateRawScores(expression_matrix, gene_set, raw_method,
                             addmodulescore_params, verbose)

  # Steps 2-4: Spatial graph, bandwidth and normalized Gaussian kernel
  kernel <- .BuildSpatialKernel(coords, k, bandwidth, verbose)

  .AssembleSEMAresult(raw, kernel, coords, k, raw_method, spatial_weight, verbose)
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
#' @param n_cores Number of parallel cores, default 1 (serial). Currently unused.
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

  if (!is.list(gene_sets)) {
    gene_sets <- list(gene_sets)
  }

  # Unnamed gene sets would otherwise be skipped by the loop below
  set_names <- names(gene_sets)
  if (is.null(set_names)) set_names <- rep("", length(gene_sets))
  missing_names <- is.na(set_names) | set_names == ""
  set_names[missing_names] <- paste0("GeneSet", which(missing_names))
  names(gene_sets) <- make.unique(set_names)

  if (verbose) {
    message(sprintf("Calculating %d gene sets...", length(gene_sets)))
  }

  coords <- .CheckCoords(coords, ncol(expression_matrix))

  # Precompute adjacency graph, bandwidth and kernel (shared by all gene sets)
  if (verbose) message("Precomputing spatial adjacency graph...")
  kernel <- .BuildSpatialKernel(coords, k, bandwidth, verbose = FALSE)

  results <- list()

  for (set_name in names(gene_sets)) {
    if (verbose) message(sprintf("Calculating gene set: %s", set_name))

    raw <- .CalculateRawScores(expression_matrix, gene_sets[[set_name]], raw_method,
                               addmodulescore_params, verbose = FALSE)

    results[[set_name]] <- .AssembleSEMAresult(raw, kernel, coords, k, raw_method,
                                               spatial_weight, verbose = FALSE)
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
      bandwidth = kernel$bandwidth
    )
  )

  return(results)
}

# Validate coordinates against the expression matrix; only the first two
# columns are used (the Gaussian kernel is 2D)
.CheckCoords <- function(coords, n_cells) {
  if (is.null(dim(coords)) || ncol(coords) < 2) {
    stop("coords must be a matrix or data.frame with at least 2 columns")
  }
  if (nrow(coords) != n_cells) {
    stop(sprintf("Cell count mismatch: expression matrix has %d cells, coordinates have %d cells",
                 n_cells, nrow(coords)))
  }
  if (ncol(coords) > 2) {
    warning(sprintf("Coordinates have %d dimensions, using only the first 2", ncol(coords)))
    coords <- coords[, 1:2, drop = FALSE]
  }
  coords
}

# Step 1: raw (non-spatial) gene set scores
.CalculateRawScores <- function(expression_matrix, gene_set, raw_method,
                                addmodulescore_params, verbose) {
  gene_set <- unique(gene_set)

  if (raw_method == "addmodulescore") {
    # Using AddModuleScore method
    addmodulescore_params$features <- list(module = gene_set)
    addmodulescore_params$expression_matrix <- expression_matrix

    raw_result <- do.call(AddModuleScoreSEMA, addmodulescore_params)
    raw_scores <- as.numeric(raw_result$scores[, 1])
    genes_used <- raw_result$features_used[[1]]

    raw_info <- list(
      method = "addmodulescore",
      control_genes = raw_result$control_genes[[1]],
      n_genes = length(genes_used),
      genes_used = genes_used
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

  list(raw_scores = raw_scores, raw_info = raw_info)
}

# Steps 2-4: KNN graph, bandwidth and row-normalized Gaussian kernel matrix
.BuildSpatialKernel <- function(coords, k, bandwidth, verbose) {
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
  adj_mat <- methods::as(adj_matrix, "TsparseMatrix")

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

  list(
    graph_result = graph_result,
    weighted_adj_norm = weighted_adj_norm,
    bandwidth = bandwidth
  )
}

# Step 5: combine raw scores with the spatial kernel into a SEMAresult
.AssembleSEMAresult <- function(raw, kernel, coords, k, raw_method,
                                spatial_weight, verbose) {
  raw_scores <- raw$raw_scores

  # Apply spatial smoothing
  neighbor_scores <- as.numeric(kernel$weighted_adj_norm %*% raw_scores)
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
      bandwidth = kernel$bandwidth
    ),
    graph = kernel$graph_result$graph,
    adjacency_matrix = kernel$graph_result$adjacency,
    weighted_adjacency = kernel$weighted_adj_norm,
    coords = coords,
    raw_info = raw$raw_info
  )

  class(result) <- "SEMAresult"

  if (verbose) message("Calculation complete!")

  return(result)
}
