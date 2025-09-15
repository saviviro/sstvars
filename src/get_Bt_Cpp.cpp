#include "arma_current.h"
// [[Rcpp::depends(RcppArmadillo)]]

//' @name get_Bt_Cpp
//' @title Calculate the impact matrix \eqn{B_t} for all \eqn{t} for models with a non-Gaussian
//'  conditional distribution with mutually independent shocks.
//'
//' @description This internal function takes a cube of matrices (\code{all_Omegas}) and a matrix of weights (\code{alpha_mt}),
//' and calculates the weighted sums of the matrices in the cube. For each row in \code{alpha_mt}, it computes
//' a weighted sum of matrices, and returns the (3D array) such that each slice contains the weighted sum of the matrices.
//' Note that the argument \code{all_Omegas} should contain the impact matrices of the regimes (and not the covariance matrices).
//'
//' @inheritParams check_Bt_Cpp
//'
//' @return An arma::cube value (3D array in R) such that each slice contains the weighted sum of the matrices,
//'   i.e, the impact matrix \eqn{B_t} for all \eqn{t}.
//'
//' @keywords internal
// [[Rcpp::export(name = "get_Bt_Cpp")]]
arma::cube get_Bt_Cpp(const arma::cube& all_Omegas,
                      const arma::mat&  alpha_mt) {
  const auto n_rows   = all_Omegas.n_rows;
  const auto n_cols   = all_Omegas.n_cols;
  const auto n_slices = all_Omegas.n_slices;
  const auto n_times  = alpha_mt.n_rows;

  if(static_cast<unsigned>(alpha_mt.n_cols) != n_slices) {
    Rcpp::stop("get_Bt_Cpp(): alpha_mt.n_cols (%u) != all_Omegas.n_slices (%u)",
               alpha_mt.n_cols, n_slices);
  }
  arma::cube weightedSums(n_rows, n_cols, n_times, arma::fill::zeros);

  for(unsigned t = 0; t < n_times; ++t) {
    arma::mat sum(n_rows, n_cols, arma::fill::zeros);
    for(unsigned k = 0; k < n_slices; ++k) {
      sum += all_Omegas.slice(k)*alpha_mt(t, k);
    }
    weightedSums.slice(t) = sum;
  }
  return weightedSums;
}

