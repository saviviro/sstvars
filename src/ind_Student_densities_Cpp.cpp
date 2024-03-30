#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @name ind_Student_densities_Cpp
//' @title Calculate log independent multivariate Student's t densities
//' @description Calculates logs of independent multivariate Student t densities with varying mean
//'   and impact matrix AND EXCLUDING the constant term of the density
//'   (the constant is calculated and added in R code). The varying impact matrix is calculated within
//'   the function from the impact matrices of the regimes and transition weights.
//'
//' @inheritParams Student_densities_Cpp
//' @inheritParams loglikelihood
//' @param impact_matrices A a size \eqn{d\times d \times M} \code{arma::cube} (3D array in R), where each slice contains an
//'  invertible (d x d) impact matrix of each regime.
//' @param distpars A numeric vector of length \eqn{d}, containing the degrees of freedom parameters for each component.
//' @details Returns \code{minval} if the impact matrix \eqn{B_t} is not invertible for some t up to the numerical tolerance \code{posdef_tol}.
//' @return A numeric vector of length \eqn{T}, where each element represents the computed density component for the corresponding observation.
//' @keywords internal
// [[Rcpp::export]]
arma::vec ind_Student_densities_Cpp(const arma::mat& obs,
                                    const arma::mat& means,
                                    const arma::cube& impact_matrices,
                                    const arma::mat& alpha_mt,
                                    const arma::vec& distpars,
                                    const arma::vec& minval,
                                    const double posdef_tol) {
  int T_obs = obs.n_rows; // The number of observations
  int d = obs.n_cols; // Dimension
  arma::vec all_lt(T_obs, arma::fill::zeros); // Vector to store results

  for(int i1 = 0; i1 < T_obs; ++i1) {
    arma::vec tdens_i1(d, arma::fill::zeros);
    arma::mat Bt = arma::zeros<arma::mat>(d, d);

    // Compute Bt as the weighted sum of impact_matrices
    for(int j = 0; j < impact_matrices.n_slices; ++j) {
      Bt += impact_matrices.slice(j)*alpha_mt(i1, j);
    }

    //arma::vec obs_minus_cmean = (obs.row(i1) - means.row(i1)).t(); // Compute obs_minus_cmean for the current observation
    arma::vec invBt_obs_minus_cmean = arma::solve(Bt, (obs.row(i1) - means.row(i1)).t()); // Solve for invBt_obs_minus_cmean

    for(int i2 = 0; i2 < d; ++i2) {
      tdens_i1(i2) = 0.5*(1 + distpars(i2))*log(1 + std::pow(invBt_obs_minus_cmean(i2), 2)/(distpars(i2) - 2));
    }
    double absdetBt = std::abs(arma::det(Bt));
    if(absdetBt < posdef_tol) {
      return minval;
    }

    all_lt(i1) = -log(absdetBt) - arma::sum(tdens_i1); // Store the results, logCd is handled outside
  }

  return all_lt;
}
