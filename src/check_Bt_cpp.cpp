#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

//' @name check_Bt_Cpp
//' @title Check Matrix B Invertibility with C++ (Internal Function)
//'
//' @description This internal function takes a cube of matrices (\code{all_Omegas}),
//' a matrix of weights (\code{alpha_mt}), and a numerical tolerance (\code{posdef_tol})
//' to check the invertibility of weighted sums of the matrices in the cube. For each row
//' in \code{alpha_mt}, it computes a weighted sum of matrices, and checks if this sum is
//' invertible by verifying that its determinant is not within the specified tolerance of zero.
//'
//' @param all_Omegas A cube (3D array) of impact matrices, with each slice being an inveritble square matrix.
//' @param alpha_mt A matrix of weights, with as many columns as there are slices in \code{all_Omegas}.
//' @param posdef_tol A strictly positive small number used as a tolerance for checking
//'        the invertibility of the matrix. The matrix is considered non-invertible if
//'        its determinant is less than this tolerance.
//'
//' @return A boolean value: `TRUE` if all weighted sums are invertible up to the specified
//'         tolerance, `FALSE` otherwise.
//'
//' @keywords internal
// [[Rcpp::export(name = "check_Bt_Cpp")]]
bool check_Bt_Cpp(const arma::cube& all_Omegas, const arma::mat& alpha_mt, double posdef_tol) {

  // Iterate through each row of alpha_mt
  for(unsigned int i = 0; i < alpha_mt.n_rows; ++i) {
    arma::mat weightedSum = arma::zeros<arma::mat>(all_Omegas.n_rows, all_Omegas.n_cols);

    // Compute the weighted sum of matrices in all_Omegas
    for(unsigned int j = 0; j < alpha_mt.n_cols; ++j) {
      weightedSum += all_Omegas.slice(j)*alpha_mt(i, j);
    }

    // Check if the determinant of the resulting matrix is within the tolerance
    if(std::abs(arma::det(weightedSum)) < posdef_tol) {
      return false; // Not invertible within tolerance
    }
  }

  return true; // All matrices passed the invertibility check
}

