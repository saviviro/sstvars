#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @name get_Bt_Cpp
//' @title Calculate the impact matrix \eqn{B_t} for all \eqn{t} for models with a non-Gaussian
//'  conditional distribution with mutually independent shocks.
//'
//' @description This internal function takes a cube of matrices (\code{all_Omegas}) and a matrix of weights (\code{alpha_mt}),
//' and calculates the weighted sums of the matrices in the cube. For each row in \code{alpha_mt}, it computes
//' a weighted sum of matrices, and returns the
//'
//' @inheritParams check_Bt_Cpp
//'
//' @return An arma::cube value (3D array in R) such that each slice contains the weighted sum of the matrices,
//'   i.e, the impact matrix \eqn{B_t} for all \eqn{t}.
//'
//' @keywords internal
// [[Rcpp::export(name = "get_Bt_Cpp")]]
arma::cube get_Bt_Cpp(const arma::cube& all_Omegas, const arma::mat& alpha_mt) {
  // Create a storage for the weighted sums of matrices
  arma::cube weightedSums = arma::zeros<arma::cube>(all_Omegas.n_rows, all_Omegas.n_cols, alpha_mt.n_rows);

  // Iterate through each row of alpha_mt
  for(unsigned int i = 0; i < alpha_mt.n_rows; ++i) {
    arma::mat weightedSum = arma::zeros<arma::mat>(all_Omegas.n_rows, all_Omegas.n_cols);

    // Compute the weighted sum of matrices in all_Omegas
    for(unsigned int j = 0; j < alpha_mt.n_cols; ++j) {
      weightedSum += all_Omegas.slice(j)*alpha_mt(i, j);
    }

    // Store the weighted sum of the matrices to weightedSums
    weightedSums.slice(i) = weightedSum;
  }

  return weightedSums;
}
