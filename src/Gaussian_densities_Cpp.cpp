#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @name Gaussian_densities_Cpp
//' @title Calculate log multivariate Gaussian densities
//' @description Calculates logs of multivariate Gaussian densities with varying mean
//'   and varying covariance matrix AND EXCLUDING the constant term of the density
//'   (the constant is calculated and added in R code).
//'
//' @param obs a \eqn{(T \times d)} matrix such that the i:th row contains the vector
//'  \eqn{y_{i}=(y_{1i},...,y_{di})} \eqn{(dx1)}. That is, the initial values are
//'  excluded but the last observations is included.
//' @param means a \eqn{(T \times d)} matrix such that the i:th row constraints the
//'   conditional mean of the process \eqn{\mu_{y,i}}.
//' @param covmats a \eqn{(d \times d \times T)} array such that the slice \code{[, , t]}
//'   contains the time t conditional covariance matrix.
//' @return a numeric vector containing the multivariate Gaussian densities, excluding the constant term.
//' @keywords internal
// [[Rcpp::export]]
arma::vec Gaussian_densities_Cpp(arma::mat obs, arma::mat means, arma::cube covmats) {
  int T_obs = obs.n_rows; // The number of observations
  int d = obs.n_cols; // The dimension d
  arma::vec vals(T_obs); // Contains the densities for each observation
  arma::mat tmp(d, 1);
  arma::mat tmp2(1, 1);
  arma::mat cholcovmat(d, d); //arma::chol(covmats);
  arma::mat inv_cholcovmat(d, d); //= arma::inv(trimatu(cholcovmat));
  // double det_term = -arma::accu(arma::log(cholcovmat.diag()));

  for(int i1 = 0; i1 < T_obs; i1++) {
    cholcovmat = arma::chol(covmats.slice(i1));
    inv_cholcovmat = arma::inv(trimatu(cholcovmat));
    tmp = (obs.row(i1) - means.row(i1))*inv_cholcovmat;
    tmp2 = dot(tmp, tmp);
    vals[i1] = -arma::accu(arma::log(cholcovmat.diag())) - 0.5*tmp2(0, 0); // The first index is zero in C++
  }

  return vals;
}


/*** R
Gaussian_densities_Cpp(obs=matrix(c(1:6), nrow=3), means=matrix(0, nrow=3, ncol=2), covmats=array(rep(c(1, 0.1, 0.1, 1), times=3), dim=c(2, 2, 3)))
*/
