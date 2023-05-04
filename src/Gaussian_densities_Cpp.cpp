#include <RcppArmadillo.h>
#include <RcppEigen.h>
// [[Rcpp::depends(RcppArmadillo, RcppEigen)]]
using namespace Rcpp;

//' @name Gaussian_densities_Cpp
//' @title Calculate log multivariate Gaussian densities
//' @description Calculates logs of multivariate Gaussian densities with varying mean
//'   and varying covariance matrix AND EXCLUDING the constant term of the density
//'   (the constant is calculated and added in R code). The varying conditional covariance
//'   matrix is calculated within the function from the regime covariance matrices and
//'   transition weights.
//'
//' @param obs a \eqn{(T \times d)} matrix such that the i:th row contains the vector
//'  \eqn{y_{i}=(y_{1i},...,y_{di})} \eqn{(dx1)}. That is, the initial values are
//'  excluded but the last observations is included.
//' @param means a \eqn{(T \times d)} matrix such that the i:th row contains the
//'   conditional mean of the process \eqn{\mu_{y,i}}.
//' @param covmats a \eqn{(d \times d \times M)} array such that the slice \code{[, , m]}
//'   contains the conditional covariance matrix of regime m.
//' @param alpha_mt a \eqn{(T \times M)} matrix such that \code{[t, m]} contains the time t
//'   transition weights of the m:th regime.
//' @return a numeric vector containing the multivariate Gaussian densities, excluding the constant term.
//' @keywords internal
// [[Rcpp::export]]
arma::mat Gaussian_densities_Cpp(arma::mat obs, arma::mat means, arma::cube covmats, arma::mat alpha_mt) {
  int T_obs = obs.n_rows; // The number of observations
  int d = obs.n_cols; // The dimension d
  int M = alpha_mt.n_cols; // The number of regimes
  arma::vec vals(T_obs); // Contains the densities for each observation
  arma::mat tmp(d, 1);
  arma::mat tmp2(1, 1);
  //arma::mat tmp3(d, 1);
  arma::mat cholcovmat(d, d);
  arma::mat inv_cholcovmat(d, d);
  //arma::mat symcondcovmat(d, d);
  //arma::mat inv_condcovmat(d, d);

  for(int i1 = 0; i1 < T_obs; i1++) {
    arma::mat condcovmat(d, d, arma::fill::zeros);
    for(int i2 = 0; i2 < M; i2++) {
      condcovmat += alpha_mt(i1, i2)*covmats.slice(i2);
    }
    //symcondcovmat = arma::symmatl(arma::trimatl(condcovmat)); // Even this does not help with the chol warning

    //cholcovmat = arma::chol(condcovmat); // arma::symmatu or symmatl forces symmetricity but does not help with chol warning

    Eigen::MatrixXd condcovmat_eigen = Rcpp::as<Eigen::MatrixXd>(Rcpp::wrap(condcovmat));
    Eigen::LLT<Eigen::MatrixXd> llt(condcovmat_eigen);
    Eigen::MatrixXd L_eigen = llt.matrixU();
    cholcovmat = Rcpp::as<arma::mat>(Rcpp::wrap(L_eigen));

    inv_cholcovmat = arma::inv(trimatu(cholcovmat));
    tmp = (obs.row(i1) - means.row(i1))*inv_cholcovmat;
    tmp2 = dot(tmp, tmp);
    vals[i1] = -arma::accu(arma::log(cholcovmat.diag())) - 0.5*tmp2(0, 0); // The first index is zero in C++
    // inv_condcovmat = arma::inv_sympd(arma::symmatu(condcovmat)); // Slower but does not fix chol warning
    // tmp3 = (obs.row(i1) - means.row(i1));
    // tmp = tmp3*inv_condcovmat*arma::trans(tmp3);
    // vals[i1] = -0.5*arma::log_det_sympd(arma::symmatu(condcovmat)) - 0.5*tmp(0, 0);
  }

  return vals;
}


