#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @name ind_skewed_t_densities_Cpp
//' @title Calculate log independent multivariate skewed t densities
//' @description Calculates logs of independent multivariate skewed t densities with varying mean
//'   and impact matrix (including the constant terms of the density). The varying impact matrix
//'   is calculated within the function from the impact matrices of the regimes and transition weights.
//'
//' @inheritParams ind_Student_densities_Cpp
//' @param all_nu a numeric vector of length \eqn{d}, containing the degrees of freedom parameters for each component.
//' @param all_lambda a numeric vector of length \eqn{d}, containing the skewness parameters for each component.
//' @details Returns \code{minval} if the impact matrix \eqn{B_t} is not invertible for some t up to the numerical tolerance
//'  \code{posdef_tol}.
//' @return A numeric vector of length \eqn{T}, where each element represents the computed density component for
//'  the corresponding observation.
//' @keywords internal
// [[Rcpp::export]]
arma::vec ind_skewed_t_densities_Cpp(const arma::mat& obs,
                                    const arma::mat& means,
                                    const arma::cube& impact_matrices,
                                    arma::mat& alpha_mt,
                                    const arma::vec& all_nu,
                                    const arma::vec& all_lambda,
                                    const arma::vec& minval,
                                    const double posdef_tol) {
  int T_obs = obs.n_rows; // The number of observations
  int d = obs.n_cols; // Dimension
  arma::vec all_lt(T_obs, arma::fill::zeros); // Vector to store results

  // Precompute parts of the formula that don't change with each observation
  arma::vec inverse_all_nu_minus_two = 1/(all_nu - 2); // (d x 1)
  arma::vec one_plus_all_nu_times_half = 0.5*(1 + all_nu); // (d x 1)
  arma::vec half_nu_plus1 = 0.5*(1 + all_nu);   // (d x 1)
  arma::vec half_nu = 0.5*all_nu;               // (d x 1)
  arma::vec half_log_nu_minus2 = 0.5*arma::log(all_nu - 2); // (d x 1)
  double half_log_pi = 0.5*std::log(arma::datum::pi); // Scalar
  arma::vec logc_i(d); // (d x 1)
  arma::vec logb_i(d); // (d x 1)
  arma::vec a_i(d); // (d x 1)
  arma::vec b_i(d); // (d x 1)
  for(int i = 0; i < d; ++i) {
    logc_i(i) = std::lgamma(half_nu_plus1(i)) - half_log_pi - half_log_nu_minus2(i) - std::lgamma(half_nu(i));
    a_i(i) = 4*all_lambda(i)*std::exp(logc_i(i))*(all_nu(i) - 2)/(all_nu(i) - 1);
    logb_i(i) = 0.5*std::log(1 + 3*std::pow(all_lambda(i), 2) - std::pow(a_i(i), 2));
    b_i(i) = std::exp(logb_i(i));
  }

  // Adjust weights in alpha_mt: values close to 1 or 0 are set to 1 or 0 to improve computation speed
  for(arma::uword i = 0; i < alpha_mt.n_rows; ++i) {
    for(arma::uword j = 0; j < alpha_mt.n_cols; ++j) {
      if(alpha_mt(i, j) > 0.999) {
        alpha_mt(i, j) = 1;
      } else if(alpha_mt(i, j) < 0.001) {
        alpha_mt(i, j) = 0;
      }
    }
  }

  // Placeholder for potential Bt precalculations
  std::vector<arma::mat> precalculatedInvBt(impact_matrices.n_slices);
  std::vector<double> precalculatedLogAbsDetBt(impact_matrices.n_slices);

  // Precalculate Bt and invBt for simplified rows
  for(arma::uword i2 = 0; i2 < impact_matrices.n_slices; ++i2) {
    precalculatedInvBt[i2] = arma::inv(impact_matrices.slice(i2));
    precalculatedLogAbsDetBt[i2] = -log(std::abs(arma::det(impact_matrices.slice(i2))));
  }

  arma::vec e_t; // Placeholder for e_t
  double absdetBt = 0.0; // Placeholder for abs det
  int which_weight_is_one = 0; // Placeholder for which weight is one

  for(int i1 = 0; i1 < T_obs; ++i1) {
    arma::vec tdens_i1(d, arma::fill::zeros);
    arma::mat Bt = arma::zeros<arma::mat>(d, d);

    bool precalc_used = false;
    for(arma::uword i2 = 0; i2 < alpha_mt.n_cols; ++i2) {
      if(alpha_mt(i1, i2) == 1) {
        e_t = precalculatedInvBt[i2]*(obs.row(i1) - means.row(i1)).t();
        which_weight_is_one = i2;
        precalc_used = true;
        break; // The rest of the weights are zero
      }
    }
    if(!precalc_used) {
      // Compute Bt as the weighted sum of impact_matrices (weights are not approximately 0 or 1)
      for(arma::uword i2 = 0; i2 < impact_matrices.n_slices; ++i2) {
        Bt += impact_matrices.slice(i2)*alpha_mt(i1, i2);
      }
      e_t = arma::solve(Bt, (obs.row(i1) - means.row(i1)).t()); // Solve for e_t
    }


    // Calculate the exp-part of the density
    for(int i2 = 0; i2 < d; ++i2) {
      tdens_i1(i2) = one_plus_all_nu_times_half(i2)*log(1 + inverse_all_nu_minus_two(i2)*
        std::pow((b_i(i2)*e_t(i2) + a_i(i2))/(1 + ((e_t(i2) < -a_i(i2)/b_i(i2)) ? -all_lambda(i2) : all_lambda(i2))), 2));
    }

    if(precalc_used) { // The determinant of Bt already calculated, always abs det pos here
      all_lt(i1) = precalculatedLogAbsDetBt[which_weight_is_one] + arma::sum(logb_i) + arma::sum(logc_i) - arma::sum(tdens_i1); // Store the results
    } else { // The determinant of Bt not yet calculated
      absdetBt = std::abs(arma::det(Bt));
      if(absdetBt < posdef_tol) {
        return minval;
      }

      all_lt(i1) = arma::sum(logb_i) + arma::sum(logc_i) -log(absdetBt) - arma::sum(tdens_i1); // Store the results, logCd is handled outside
    }
  }

  return all_lt;
}

