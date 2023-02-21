#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

//' @name GaussianDensities
//' @title Calculate log multivariate Gaussian densities
//' @description Calculates log of multivariate Gaussian densities with varying mean and covariance matrix.
//'
//' @param const_term the constant term \code{-0.5*d*log(2*pi)} that does not vary.
//' @param obs_minus_cmean a \eqn{(T \times d)} matrix containing the observation minus
//'   the (conditional) mean in each row.
//' @return a numeric vector containing the multivariate Gaussian densities
// [[Rcpp::export]]
NumericVector GaussianDensities(double const_term, NumericMatrix obs_minus_cmean) {
  int T_obs = obs_minus_cmean.nrow(); // The number of observations
  NumericVector out(T_obs); // Contains the densities for each observation

  for(int i1 = 0; i1 < T_obs; i1++) {
    out[i1] = i1;
  }

  return out;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
// Compile with: sourceCpp("/Users/savi/Documents/sstvars/src/GaussianDensities.cpp") or with Source buttom

/*** R
GaussianDensities(1, matrix(c(1:6), nrow=3))
*/
