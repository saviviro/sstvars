#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//' @name spectral_radius
//' @title Compute the spectral radius of a square matrix.
//'
//' @description This internal helper function computes the spectral radius of a square matrix,
//'   defined as the maximum modulus of its eigenvalues. The computation is carried
//'   out using an exact eigenvalue decomposition via \code{RcppArmadillo}, ensuring
//'   numerical accuracy and mathematical correctness.
//' @param A an \code{arma::mat} representing a square real-valued matrix.
//' @return A numeric scalar giving the spectral radius of \code{A}.
//'
//' @keywords internal
inline double spectral_radius(const arma::mat& A) {
  arma::cx_vec eigval = arma::eig_gen(A);
  return arma::abs(eigval).max();
}

//' @name bound_jsr_G_Cpp
//' @title Compute lower and upper bounds for the joint spectral radius
//'  using Gripenberg's branch-and-bound algorithm (C++ implementation).
//' @description This internal function computes lower and upper bounds for the joint
//'   spectral radius (JSR) of a finite set of square matrices using a branch-and-bound
//'   algorithm of the type proposed by Gripenberg (1996).
//' @param S an \code{arma::cube} of dimension \eqn{n \times n \times m} containing the
//'   set of \eqn{m} square matrices whose joint spectral radius is to be bounded.
//' @param epsilon a positive scalar specifying the desired tolerance for the
//'   difference between the upper and lower bounds. The algorithm terminates once
//'   the upper bound minus the lower bound is less than or equal to \code{epsilon}.
//' @param maxit An integer specifying the maximum number of iterations of the
//'   branch-and-bound algorithm.
//' @param print_progress Logical; if \code{TRUE}, progress information (iteration
//'   number, current bounds, and number of candidate products) is printed to the
//'   console during execution.
//'
//' @return A numeric vector of length two containing the certified lower and upper
//'   bounds for the joint spectral radius, in the order \code{c(lower_bound, upper_bound)}.
//'
//' @details The implementation mirrors the logic of the corresponding R function  \code{bound_jsr_G},
//'   but performs the entire iteration process in C++ using \code{Rcpp} and \code{RcppArmadillo}
//'   for improved performance.
//'
//'   The algorithm iteratively constructs candidate matrix products, prunes them using
//'   norm-based upper bounds, and tightens the lower and upper bounds on the joint
//'   spectral radius until the desired tolerance is reached or a maximum number of
//'   iterations is exceeded.
//' @keywords internal
// [[Rcpp::export(name = "bound_jsr_G_Cpp")]]
NumericVector bound_jsr_G_Cpp(const arma::cube& S, double epsilon = 0.01, int maxit = 1000, bool print_progress = true) {
  const int m = S.n_slices;

  // Candidate structure
  struct Candidate {
    arma::mat P;
    double mu;
  };

  std::vector<Candidate> T, T_next;

  // --- Initialization (k = 1) ---
  double alpha = 0.0;
  double beta  = 0.0;

  T.reserve(m);

  for(int i = 0; i < m; ++i) {
    arma::mat P = S.slice(i);
    double normP = arma::norm(P, 2); // spectral norm
    double rhoP  = spectral_radius(P);

    alpha = std::max(alpha, rhoP);
    beta  = std::max(beta,  normP);

    T.push_back({P, normP});
  }

  if(print_progress) {
    Rcout << "Iteration: 1, bounds: "
          << alpha << ", " << beta << "\n";
  }

  // --- Main loop ---
  for(int k = 2; k <= maxit; ++k) {
    T_next.clear();
    T_next.reserve(T.size() * m);

    double alpha_prev = alpha;
    double beta_k     = 0.0;

    // Branch & bound
    for(const auto& cand : T) {
      for (int i = 0; i < m; ++i) {
        arma::mat P_new = S.slice(i)*cand.P;

        double mu_new = std::min(cand.mu, std::pow(arma::norm(P_new, 2), 1.0/k));

        if (mu_new > alpha_prev + epsilon) {
          T_next.push_back({P_new, mu_new});
          beta_k = std::max(beta_k, mu_new);
        }
      }
    }

    if(T_next.empty()) {
      if(print_progress) {
        Rcout << "Finished (no candidates).\n";
      }
      break;
    }

    // Update alpha using exact eigenvalues
    for(const auto& cand : T_next) {
      double rho = spectral_radius(cand.P);
      alpha = std::max(alpha, std::pow(rho, 1.0/k));
    }

    beta = std::min(beta, std::max(alpha_prev + epsilon, beta_k));

    if(print_progress && k % 10 == 0) {
      Rcout << "Iteration: " << k
            << ", bounds: " << alpha << ", " << beta
            << ", candidates: " << T_next.size() << "\n";
    }

    if(beta - alpha <= epsilon) {
      if (print_progress) {
        Rcout << "Finished (epsilon reached).\n";
      }
      break;
    }

    T.swap(T_next);
  }

  return NumericVector::create(alpha, beta);
}
