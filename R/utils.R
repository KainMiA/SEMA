# Helper functions

#' Create spatial adjacency graph (KNN)
#'
#' Quickly create KNN adjacency graph using RANN
#'
#' @param coords Coordinate matrix, n×2
#' @param k Number of nearest neighbors, default 30
#' @param symmetric Whether to create symmetric adjacency matrix, default TRUE
#' @return List containing adjacency matrix and graph object
#' @export
CreateSpatialGraph <- function(coords, k = 30, symmetric = TRUE) {
  n_cells <- nrow(coords)
  
  # 调整k值
  if (k >= n_cells) {
    k <- n_cells - 1
    warning(sprintf("k value adjusted to %d (maximum available)", k))
  }
  
  nn <- RANN::nn2(coords, k = k + 1)  
  
  neighbor_indices <- nn$nn.idx[, 2:(k + 1), drop = FALSE]
  
  i_vec <- rep(1:n_cells, each = k)
  j_vec <- as.vector(t(neighbor_indices))
  
  valid_idx <- !is.na(j_vec) & j_vec >= 1 & j_vec <= n_cells
  i_vec <- i_vec[valid_idx]
  j_vec <- j_vec[valid_idx]
  
  if (symmetric) {
    adj <- Matrix::sparseMatrix(
      i = c(i_vec, j_vec),
      j = c(j_vec, i_vec),
      x = 1,
      dims = c(n_cells, n_cells),
      symmetric = FALSE
    )
    adj <- as(adj, "symmetricMatrix")
  } else {
    adj <- Matrix::sparseMatrix(
      i = i_vec,
      j = j_vec,
      x = 1,
      dims = c(n_cells, n_cells),
      symmetric = FALSE
    )
  }
  
  g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected", weighted = NULL)
  
  list(
    adjacency = adj,
    graph = g,
    knn_result = nn,
    parameters = list(k = k, symmetric = symmetric)
  )
}
#' Estimate Gaussian kernel bandwidth
#'
#' Automatically estimate bandwidth parameter for Gaussian kernel
#'
#' @param coords Coordinate matrix
#' @param method Estimation method: "median_distance" (default) or "mean_distance"
#' @param k Number of nearest neighbors for estimation, default 10
#' @return Estimated bandwidth value
#' @export
EstimateBandwidth <- function(coords, method = "median_distance", k = 10) {
  coords_mat <- as.matrix(coords)
  n_cells <- nrow(coords_mat)
  
  if (n_cells < 2) {
    stop("At least 2 cells are needed to estimate bandwidth")
  }
  
  # Adjust k value
  k_adj <- min(k, n_cells - 1)
  
  if (method == "median_distance") {
    # Use median of KNN distances
    nn <- RANN::nn2(coords_mat, k = k_adj + 1)
    distances <- nn$nn.dists[, 2]  # First column is self, distance 0
    bandwidth <- stats::median(distances)
    
  } else if (method == "mean_distance") {
    # Use mean of KNN distances
    nn <- RANN::nn2(coords_mat, k = k_adj + 1)
    distances <- nn$nn.dists[, 2]
    bandwidth <- mean(distances)
    
  } else {
    stop("Method must be either 'median_distance' or 'mean_distance'")
  }
  
  # Multiply by coefficient to get appropriate bandwidth
  bandwidth <- bandwidth * 2
  
  message(sprintf("Estimated bandwidth: %.4f (method: %s, k=%d)", bandwidth, method, k_adj))
  return(bandwidth)
}

#' Seurat-style AddModuleScore
#'
#' Implements Seurat-style gene set scoring, including background correction
#'
#' @param expression_matrix Expression matrix, genes×cells
#' @param features Gene set list
#' @param n_bins Number of expression bins, default 24
#' @param ctrl Number of control genes per bin, default 100
#' @param seed Random seed, default 123
#' @return List containing scores and related information
#' @export
AddModuleScoreSEMA <- function(expression_matrix, features,
                              n_bins = 24, ctrl = 100,
                              seed = 123) {
  
  set.seed(seed)
  
  # Ensure features is a list
  if (!is.list(features)) {
    features <- list(module = features)
  }
  
  # Filter genes not in expression matrix
  all_genes <- rownames(expression_matrix)
  features <- lapply(features, function(x) {
    x[x %in% all_genes]
  })
  
  # Remove empty gene sets
  features <- features[sapply(features, length) > 0]
  
  if (length(features) == 0) {
    stop("No genes available")
  }
  
  # Calculate gene expression levels (mean expression)
  gene_means <- Matrix::rowMeans(expression_matrix)
  
  # Only use genes with expression > 0
  non_zero_genes <- gene_means > 0
  
  if (sum(non_zero_genes) < 2) {
    stop("Too few genes with expression > 0")
  }
  
  # Adjust bin count
  if (sum(non_zero_genes) < n_bins) {
    n_bins <- sum(non_zero_genes)
    message(sprintf("Adjusting bin count to %d", n_bins))
  }
  
  # Bin genes by expression level
  breaks <- stats::quantile(gene_means[non_zero_genes],
                           probs = seq(0, 1, length.out = n_bins + 1))
  breaks[1] <- -Inf
  breaks[length(breaks)] <- Inf
  
  gene_bins <- cut(gene_means, breaks = breaks, labels = FALSE, include.lowest = TRUE)
  names(gene_bins) <- names(gene_means)
  
  # Calculate scores for each gene set
  scores <- matrix(0, nrow = ncol(expression_matrix), ncol = length(features))
  colnames(scores) <- names(features)
  rownames(scores) <- colnames(expression_matrix)
  
  control_genes_list <- list()
  
  for (i in seq_along(features)) {
    module_genes <- features[[i]]
    
    # Select control genes for each bin
    ctrl_genes <- character(0)
    
    for (bin in unique(gene_bins[module_genes])) {
      bin_genes <- names(gene_bins[gene_bins == bin])
      bin_module_genes <- intersect(module_genes, bin_genes)
      
      if (length(bin_module_genes) > 0) {
        bin_ctrl_genes <- setdiff(bin_genes, bin_module_genes)
        
        if (length(bin_ctrl_genes) > ctrl) {
          bin_ctrl_genes <- sample(bin_ctrl_genes, ctrl)
        }
        
        ctrl_genes <- c(ctrl_genes, bin_ctrl_genes)
      }
    }
    
    control_genes_list[[names(features)[i]]] <- ctrl_genes
    
    # Calculate module score: mean expression of module genes - mean expression of control genes
    if (length(module_genes) > 0) {
      module_expr <- expression_matrix[module_genes, , drop = FALSE]
      module_score <- Matrix::colMeans(module_expr)
      
      if (length(ctrl_genes) > 0) {
        ctrl_expr <- expression_matrix[ctrl_genes, , drop = FALSE]
        ctrl_score <- Matrix::colMeans(ctrl_expr)
        scores[, i] <- module_score - ctrl_score
      } else {
        scores[, i] <- module_score
      }
    }
  }
  
  # Return results
  list(
    scores = scores,
    features_used = features,
    control_genes = control_genes_list,
    gene_bins = gene_bins,
    parameters = list(n_bins = n_bins, ctrl = ctrl, seed = seed)
  )
}