#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


//' @name GaussianDensities2
//' @title Calculate log multivariate Gaussian densities
//' @description Calculates logs of multivariate Gaussian densities with varying mean but
//'   constant covariance matrix AND EXCLUDING the constant term and determinant.
//'   Ehkä laita laskemaan pelkät kvadraattiset muodot?
//'
//' @inheritParams GaussianDensities
//' @param obs a \eqn{(T \times dp)} such that the i:th row of the matrix contains
//'  the vector \eqn{(y_{i-1}',...,y_{i-p}') - }\code{rep(all_mu[,m], times=p)}
//'  \eqn{((dp)x1)}, where \eqn{y_{i}=(y_{1i},...,y_{di})} \eqn{(dx1)}. That is, the initial values are
//'  included but the last observations not (used in relative dens transition weights).
//' @param mean
//' @param covmat the \eqn{(dp \times dp)} covariance matrix that is the same for all
//'   observations.
//' @return a numeric vector containing the multivariate Gaussian densities, excluding the constant term.
// [[Rcpp::export]]
arma::vec GaussianDensities2(arma::mat obs, arma::mat mean, arma::mat covmat) {
  int T_obs = obs.n_rows; // The number of observations
  arma::vec vals(T_obs); // Contains the densities for each observation
  arma::mat tmp(obs.n_cols, 1);
  arma::mat tmp2(1, 1);
  arma::mat cholcovmat = arma::chol(covmat);
  arma::mat inv_cholcovmat = arma::inv(trimatu(cholcovmat));
  double det_term = -arma::accu(arma::log(cholcovmat.diag()));

  for(int i1 = 0; i1 < T_obs; i1++) {
    tmp = (obs.row(i1) - mean)*inv_cholcovmat;
    tmp2 = dot(tmp, tmp);
    vals[i1] = det_term - 0.5*tmp2(0, 0);
  }

  return vals;
}


/*** R
GaussianDensities2(obs=matrix(c(1:6), nrow=3), mean=matrix(c(0, 0), nrow=1), covmat=matrix(c(1, 0.1, 0.1, 1), nrow=2))
*/


