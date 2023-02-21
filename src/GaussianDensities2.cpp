#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @name GaussianDensities2
//' @title Calculate log multivariate Gaussian densities
//' @description Calculates logs of multivariate Gaussian densities with varying mean but
//'   constant covariance matrix.
//'
//' @inheritParams GaussianDensities
//' @param obs_minus_mean a \eqn{(T \times dp)} such that the i:th row of the matrix contains
//'  the vector \eqn{(y_{i-1}',...,y_{i-p}') - }\code{rep(all_mu[,m], times=p)}
//'  \eqn{((dp)x1)}, where \eqn{y_{i}=(y_{1i},...,y_{di})} \eqn{(dx1)}. That is, the initial values are
//'  included but the last observations not (used in relative dens transition weights).
//' @param covmat the \eqn{(dp \times dp)} covariance matrix that is the same for all
//'   observations.
//' @return a numeric vector containing the multivariate Gaussian densities
// [[Rcpp::export]]
arma::vec GaussianDensities2(double const_term, arma::mat obs_minus_mean, arma::mat covmat) {
  int T_obs = obs_minus_mean.n_rows; // The number of observations
  int dp = obs_minus_mean.n_cols; // The length of each individual observation vector
  arma::vec out(T_obs); // Contains the densities for each observation
  arma::mat tmp(dp, 1); // A column vector containing a row from obs_minus_mean
  // Determinantin voi laskea R:lläkin? Tai ehkä parempi täällä?

  for(int i1 = 0; i1 < T_obs; i1++) {
    tmp = obs_minus_mean.row(i1);
    arma::mat tmp2 = tmp.t();
    arma::mat testi = tmp*tmp2;
    out[i1] = testi(0, 0);
  }

  return out;
}


/*** R
GaussianDensities2(13, matrix(c(1:6), nrow=3), matrix(c(1, 0.1, 0.1, 1), nrow=2))
*/


