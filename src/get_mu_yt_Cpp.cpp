#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @name get_mu_yt_Cpp
//' @title Calculate the conditional means of the process
//' @description Calculates the conditional means \eqn{\mu_{y,t}} of the process
//'
//' @param obs a \eqn{(T \times dp)} matrix such that the i:th row contains the vector
//'   \eqn{(y_{i-1},...,y_{i-p})} \eqn{((dp)x1)}, where \eqn{y_{i}=(y_{1i},...,y_{di})}
//'   \eqn{(dx1)}. That is, the initial values are included but the last observations not.
//' @param all_phi0 a \eqn{(d \times M)} matrix such that the m:th column contains the
//'   intercept parameters of the m:th regime.
//' @param all_A a \eqn{(d \times dp \times M)} array such that the slice \code{[, , m]}
//'   contains the AR matrices of the m:th regime cbinded together: \eqn{[A_{m,1}:...:A_{m,p}]}.
//' @param alpha_mt a \eqn{(T \times M)} matrix such that \code{[t, m]} contains the time t
//'   transition weights of the m:th regime.
//' @return a \eqn{(T \times d)} matrix such that the i:th row contains the conditional
//'   mean of the process.
//' @keywords internal
// [[Rcpp::export]]
arma::mat get_mu_yt_Cpp(arma::mat obs, arma::mat all_phi0, arma::cube all_A, arma::mat alpha_mt) {
  int T_obs = obs.n_rows; // The number of observations
  int d = all_phi0.n_rows; // The dimension d
  int M = alpha_mt.n_cols; // The number of regimes
  arma::mat mu_yt(T_obs, d);
  arma::mat trans_obsrow(d, 1);

  for(int i1 = 0; i1 < T_obs; i1++) {
    arma::mat time_i1_mean(d, 1, arma::fill::zeros);
    trans_obsrow = arma::trans(obs.row(i1));

    for(int m = 0; m < M; m++) {
      time_i1_mean += alpha_mt(i1, m)*(all_phi0.col(m) + all_A.slice(m)*trans_obsrow);
    }
    mu_yt.row(i1) = time_i1_mean.t();
  }

  return mu_yt;
}

