#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;

//' @name Student_densities_Cpp
//' @title Calculate log multivariate Student's t densities
//' @description Calculates logs of multivariate Student t densities with varying mean
//'   and varying covariance matrix AND EXCLUDING the constant term of the density
//'   (the constant is calculated and added in R code). The varying conditional covariance
//'   matrix is calculated within the function from the regime covariance matrices and
//'   transition weights.
//'
//' @inheritParams Gaussian_densities_Cpp
//' @param df the degrees of freedom parameter value (assumed larger than two).
//' @details Note that the parametrization is with the covariance matrix and not the scale matrix.
//' @return a numeric vector containing the multivariate Student's t densities, excluding the constant term.
//' @keywords internal
// [[Rcpp::export]]
arma::mat Student_densities_Cpp(arma::mat obs, arma::mat means, arma::cube covmats, arma::mat alpha_mt, double df) {
  int T_obs = obs.n_rows; // The number of observations
  int d = obs.n_cols; // The dimension d
  int M = alpha_mt.n_cols; // The number of regimes
  arma::vec vals(T_obs); // Contains the densities for each observation
  arma::mat tmp(d, 1);
  arma::mat tmp2(1, 1);
  arma::mat cholcovmat(d, d);
  arma::mat inv_cholcovmat(d, d);

  for(int i1 = 0; i1 < T_obs; i1++) {
    arma::mat condcovmat(d, d, arma::fill::zeros);
    for(int i2 = 0; i2 < M; i2++) {
      condcovmat += alpha_mt(i1, i2)*covmats.slice(i2);
    }
    // arma::symmatu or symmatl forces symmetricity but does not help with chol warning: this still produces warning sometimes it
    // not being symmetric:
    //cholcovmat = arma::chol(arma::symmatu(arma::trimatu(condcovmat).eval()));

    cholcovmat = arma::chol(condcovmat);
    inv_cholcovmat = arma::inv(trimatu(cholcovmat));
    tmp = (obs.row(i1) - means.row(i1))*inv_cholcovmat;
    tmp2 = dot(tmp, tmp);
    vals[i1] = -arma::accu(arma::log(cholcovmat.diag())) - 0.5*(d + df)*log(1 + tmp2(0, 0)/(df - 2)); // The first index is zero in C++
  }

  return vals;
}


