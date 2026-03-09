#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>

using namespace Rcpp;

// Gaussian kernel calculation function
// bandwidth explanation:
// - Controls decay rate of neighbor cell influence
// - Larger bandwidth means more distant cells have greater influence
// - Smaller bandwidth means only very close cells have significant influence
// [[Rcpp::export]]
List sparse_gaussian_kernel_cpp(IntegerVector i, IntegerVector j,
                               NumericMatrix coords,
                               double bandwidth) {
  int n_edges = i.size();
  
  // Store weights and distances
  NumericVector weights(n_edges);
  NumericVector distances(n_edges);
  
  // Precompute constant
  double constant = -1.0 / (2.0 * bandwidth * bandwidth);
  
  for (int idx = 0; idx < n_edges; idx++) {
    int row_idx = i[idx] - 1; // Convert to 0-based index
    int col_idx = j[idx] - 1;
    
    // Calculate Euclidean distance
    double dx = coords(row_idx, 0) - coords(col_idx, 0);
    double dy = coords(row_idx, 1) - coords(col_idx, 1);
    double dist = std::sqrt(dx * dx + dy * dy);
    distances[idx] = dist;
    
    // Calculate Gaussian kernel weight: exp(-dist² / (2 * bandwidth²))
    weights[idx] = std::exp(constant * dist * dist);
  }
  
  return List::create(
    _["weights"] = weights,
    _["distances"] = distances,
    _["bandwidth"] = bandwidth
  );
}

// Fast sparse matrix normalization
// [[Rcpp::export]]
List normalize_sparse_weights(IntegerVector i, IntegerVector j,
                             NumericVector weights,
                             int nrow) {
  NumericVector row_sums(nrow, 0.0);
  int n_edges = i.size();
  
  // Calculate weight sum for each row
  for (int idx = 0; idx < n_edges; idx++) {
    int row_idx = i[idx] - 1;
    row_sums[row_idx] += weights[idx];
  }
  
  // Normalize weights
  NumericVector normalized_weights(n_edges);
  for (int idx = 0; idx < n_edges; idx++) {
    int row_idx = i[idx] - 1;
    if (row_sums[row_idx] > 0) {
      normalized_weights[idx] = weights[idx] / row_sums[row_idx];
    } else {
      normalized_weights[idx] = 0.0;
    }
  }
  
  return List::create(
    _["normalized_weights"] = normalized_weights,
    _["row_sums"] = row_sums
  );
}