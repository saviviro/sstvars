// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <bitset>
using namespace Rcpp;

//' @name get_permanent_Cpp
//' @title Calculate the permanent of a matrix
//' @description Calculates the permanent of a square matrix using the Ryser's formula.
//'
//' @param A a square matrix
//' @details Calculates the permanent of a square matrix using Ryser's formula.
//' @return Returns the permanent of the matrix A (a real number)
//' @keywords internal
// [[Rcpp::export]]
double get_permanent_Cpp(const arma::mat& A) {
  int n = A.n_rows;
  double perm = 0;

  for (int subset = 1; subset < (1 << n); subset++) {
    double prod = 1;

    for (int row = 0; row < n; row++) {
      double sum = 0;
      for (int col = 0; col < n; col++) {
        if (subset & (1 << col)) {
          sum += A(row, col);
        }
      }
      prod *= sum;
    }

    int bits = std::bitset<32>(subset).count(); // count the number of bits set to 1 in subset
    perm += std::pow(-1, n - bits) * prod;
  }

  return perm;
}
