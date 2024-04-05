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
                                    arma::mat& alpha_mt,
                                    const arma::vec& distpars,
                                    const arma::vec& minval,
                                    const double posdef_tol) {
  int T_obs = obs.n_rows; // The number of observations
  int d = obs.n_cols; // Dimension
  arma::vec all_lt(T_obs, arma::fill::zeros); // Vector to store results

  // Precompute parts of the formula that don't change with each observation
  arma::vec inverse_distpars_minus_two = 1/(distpars - 2);
  arma::vec one_plus_distpars_times_half = 0.5*(1 + distpars);

  // Adjust weights in alpha_mt: values close to 1 or 0 are set to 1 or 0 to improve computation speed
  for(arma::uword i = 0; i < alpha_mt.n_rows; ++i) {
    for(arma::uword j = 0; j < alpha_mt.n_cols; ++j) {
      if(alpha_mt(i, j) > 0.99) {
        alpha_mt(i, j) = 1;
      } else if(alpha_mt(i, j) < 0.01) {
        alpha_mt(i, j) = 0;
      }
    }
  }

  // Placeholder for potential Bt precalculations (consider using map or vector of matrices)
  std::vector<arma::mat> precalculatedBt(T_obs);
  std::vector<arma::mat> precalculatedInvBt(T_obs);
  std::vector<bool> precalculationAvailable(T_obs, false);
  std::vector<double> precalculatedlogAbsDetBt(T_obs, 0);

  // Precalculate Bt and invBt for simplified rows
  for(int i1 = 0; i1 < T_obs; ++i1) {
    for(int i2 = 0; i2 < alpha_mt.n_cols; ++i2) {
      if(alpha_mt(i1, i2) == 1) { // One weight is 1 and the rest are 0, so need to calculate for the one with 1 only.
        precalculatedBt[i1] = impact_matrices.slice(i2);
        precalculatedlogAbsDetBt[i1] = -log(std::abs(arma::det(precalculatedBt[i1])));
        precalculatedInvBt[i1] = arma::inv(precalculatedBt[i1]); // Bt is always invertible here since it is just B_m
        precalculationAvailable[i1] = true;
        break; // No need to check other weights for this row
      }
    }
  }

  arma::vec invBt_obs_minus_cmean; // Placeholder for invBt_obs_minus_cmean
  double absdetBt = 0.0; // Placeholder for abs det

  for(int i1 = 0; i1 < T_obs; ++i1) {
    arma::vec tdens_i1(d, arma::fill::zeros);
    arma::mat Bt = arma::zeros<arma::mat>(d, d);

    if(precalculationAvailable[i1]) { // B_t and its inverse already calculated
      Bt = precalculatedBt[i1];
      invBt_obs_minus_cmean = precalculatedInvBt[i1]*(obs.row(i1) - means.row(i1)).t();
    } else {
      // Compute Bt as the weighted sum of impact_matrices (weights are not approximately 0 or 1)
      for(int i2 = 0; i2 < impact_matrices.n_slices; ++i2) {
        Bt += impact_matrices.slice(i2)*alpha_mt(i1, i2);
      }

      //arma::vec obs_minus_cmean = (obs.row(i1) - means.row(i1)).t(); // Compute obs_minus_cmean for the current observation
      invBt_obs_minus_cmean = arma::solve(Bt, (obs.row(i1) - means.row(i1)).t()); // Solve for invBt_obs_minus_cmean
    }

    for(int i2 = 0; i2 < d; ++i2) {
      tdens_i1(i2) = one_plus_distpars_times_half(i2)*log(1 + std::pow(invBt_obs_minus_cmean(i2), 2)*inverse_distpars_minus_two(i2));
    }

    if(precalculationAvailable[i1]) { // The determinant of Bt already calculated, always abs det pos here
      all_lt(i1) = precalculatedlogAbsDetBt[i1] - arma::sum(tdens_i1); // Store the results, logCd is handled outside
    } else {
      absdetBt = std::abs(arma::det(Bt));
      if(absdetBt < posdef_tol) {
        return minval;
      }
      all_lt(i1) = -log(absdetBt) - arma::sum(tdens_i1); // Store the results, logCd is handled outside
    }
  }

  return all_lt;
}
