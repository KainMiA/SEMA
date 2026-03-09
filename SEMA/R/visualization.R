# Visualization functions

#' Visualize SEMA results
#'
#' Plot spatial gene set expression scores
#'
#' @param sema_result SEMA result object
#' @param type Visualization type: "spatial" (default), "raw", or "comparison"
#' @param title Plot title
#' @param point_size Point size, default 2
#' @param alpha Transparency, default 0.8
#' @param color_palette Color scheme: "viridis" (default) or "rdbu"
#' @param ... Other parameters passed to ggplot
#'
#' @return ggplot object
#' @export
VisualizeSEMA <- function(sema_result, type = "spatial",
                         title = NULL,
                         point_size = 2,
                         alpha = 0.8,
                         color_palette = "viridis",
                         ...) {
  
  if (!inherits(sema_result, "SEMAresult")) {
    stop("sema_result must be SEMAresult object")
  }
  
  coords <- sema_result$coords
  
  if (type == "spatial") {
    scores <- sema_result$spatial_scores
    if (is.null(title)) {
      title <- "Spatial Gene Set Scores"
    }
    
  } else if (type == "raw") {
    scores <- sema_result$raw_scores
    if (is.null(title)) {
      title <- "Raw Gene Set Scores"
    }
    
  } else if (type == "comparison") {
    # Compare raw scores and spatial scores
    plot_data <- data.frame(
      x = rep(coords[, 1], 2),
      y = rep(coords[, 2], 2),
      score = c(sema_result$raw_scores, sema_result$spatial_scores),
      type = factor(rep(c("Raw Scores", "Spatial Scores"),
                       each = length(sema_result$raw_scores)),
                   levels = c("Raw Scores", "Spatial Scores"))
    )
    
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = score)) +
      ggplot2::geom_point(size = point_size, alpha = alpha) +
      ggplot2::facet_wrap(~ type, ncol = 2) +
      ggplot2::labs(x = "X Coordinate", y = "Y Coordinate") +
      ggplot2::theme_minimal() +
      ggplot2::coord_fixed() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        strip.text = ggplot2::element_text(face = "bold", size = 10)
      )
    
    if (is.null(title)) {
      p <- p + ggplot2::ggtitle("Comparison of Raw and Spatial Scores")
    } else {
      p <- p + ggplot2::ggtitle(title)
    }
    
    if (color_palette == "viridis") {
      p <- p + ggplot2::scale_color_viridis_c(option = "plasma", name = "Score")
    } else if (color_palette == "rdbu") {
      p <- p + ggplot2::scale_color_gradient2(
        low = "blue", mid = "white", high = "red",
        midpoint = stats::median(plot_data$score, na.rm = TRUE),
        name = "Score"
      )
    }
    
    return(p)
    
  } else {
    stop("type must be either 'spatial', 'raw', or 'comparison'")
  }
  
  # Create single score plot
  plot_data <- data.frame(
    x = coords[, 1],
    y = coords[, 2],
    score = scores
  )
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = score)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::labs(title = title, x = "X Coordinate", y = "Y Coordinate") +
    ggplot2::theme_minimal() +
    ggplot2::coord_fixed() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  if (color_palette == "viridis") {
    p <- p + ggplot2::scale_color_viridis_c(option = "plasma", name = "Score")
  } else if (color_palette == "rdbu") {
    p <- p + ggplot2::scale_color_gradient2(
      low = "blue", mid = "white", high = "red",
      midpoint = stats::median(scores, na.rm = TRUE),
      name = "Score"
    )
  }
  
  return(p)
}

#' Visualize spatial adjacency graph
#'
#' Plot spatial connection relationships between cells
#'
#' @param sema_result SEMA result object
#' @param n_edges Maximum number of edges to display (to avoid overcrowding), default 1000
#' @param point_size Point size, default 1.5
#' @param edge_alpha Edge transparency, default 0.1
#' @param color_by What to color by: "degree" (default, connectivity) or "score" (score)
#' @param title Plot title
#'
#' @return ggplot object
#' @export
VisualizeSpatialGraph <- function(sema_result, n_edges = 1000,
                                 point_size = 1.5,
                                 edge_alpha = 0.1,
                                 color_by = "degree",
                                 title = "Spatial Graph") {
  
  if (!inherits(sema_result, "SEMAresult")) {
    stop("sema_result must be SEMAresult object")
  }
  
  coords <- sema_result$coords
  graph <- sema_result$graph
  
  # Get edge coordinates
  edges <- igraph::as_edgelist(graph)
  
  # Randomly sample edges to avoid overcrowding
  if (nrow(edges) > n_edges) {
    set.seed(123)
    edges <- edges[sample(1:nrow(edges), n_edges), ]
  }
  
  # Prepare edge data
  edge_data <- data.frame(
    x = coords[edges[, 1], 1],
    y = coords[edges[, 1], 2],
    xend = coords[edges[, 2], 1],
    yend = coords[edges[, 2], 2]
  )
  
  # Prepare point data
  if (color_by == "degree") {
    degrees <- igraph::degree(graph)
    point_data <- data.frame(
      x = coords[, 1],
      y = coords[, 2],
      value = degrees
    )
    color_name <- "Degree"
  } else if (color_by == "score") {
    point_data <- data.frame(
      x = coords[, 1],
      y = coords[, 2],
      value = sema_result$spatial_scores
    )
    color_name <- "Score"
  } else {
    stop("color_by must be either 'degree' or 'score'")
  }
  
  # Create plot
  p <- ggplot2::ggplot() +
    # Draw edges
    ggplot2::geom_segment(
      data = edge_data,
      ggplot2::aes(x = x, y = y, xend = xend, yend = yend),
      alpha = edge_alpha,
      color = "gray50"
    ) +
    # Draw points
    ggplot2::geom_point(
      data = point_data,
      ggplot2::aes(x = x, y = y, color = value),
      size = point_size
    ) +
    ggplot2::scale_color_viridis_c(option = "plasma", name = color_name) +
    ggplot2::labs(title = title, x = "X Coordinate", y = "Y Coordinate") +
    ggplot2::theme_minimal() +
    ggplot2::coord_fixed() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(hjust = 0.5, face = "bold")
    )
  
  return(p)
}

#' Print SEMA result summary
#'
#' @param x SEMA result object
#' @param ... Other parameters
#' @export
print.SEMAresult <- function(x, ...) {
  cat("SEMA Result Object\n")
  cat("==================\n")
  cat(sprintf("Number of cells: %d\n", length(x$spatial_scores)))
  cat(sprintf("Raw score method: %s\n", x$parameters$raw_method))
  cat(sprintf("Spatial weight: %.2f\n", x$parameters$spatial_weight))
  cat(sprintf("Bandwidth: %.4f\n", x$parameters$bandwidth))
  cat(sprintf("KNN parameter k: %d\n", x$parameters$k))
  
  if (!is.null(x$raw_info$n_genes)) {
    cat(sprintf("Number of genes used: %d\n", x$raw_info$n_genes))
  }
  
  cat("\nScore statistics:\n")
  cat(sprintf("  Raw scores: mean=%.4f, range=[%.4f, %.4f]\n",
              mean(x$raw_scores), min(x$raw_scores), max(x$raw_scores)))
  cat(sprintf("  Spatial scores: mean=%.4f, range=[%.4f, %.4f]\n",
              mean(x$spatial_scores), min(x$spatial_scores), max(x$spatial_scores)))
}

#' @export
print.SEMAresults <- function(x, ...) {
  cat("Multiple SEMA Result Objects\n")
  cat("=============================\n")
  summary_info <- attr(x, "summary")
  
  if (!is.null(summary_info)) {
    cat(sprintf("Number of gene sets: %d\n", summary_info$n_sets))
    cat(sprintf("Number of cells: %d\n", summary_info$n_cells))
    cat(sprintf("Parameters: k=%d, weight=%.2f, bandwidth=%.4f\n",
                summary_info$parameters$k,
                summary_info$parameters$spatial_weight,
                summary_info$parameters$bandwidth))
  }
  
  cat("\nIncluded gene sets:\n")
  for (set_name in names(x)) {
    cat(sprintf("  - %s\n", set_name))
  }
}