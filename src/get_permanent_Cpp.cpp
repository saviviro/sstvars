// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <cmath>
#include <vector>

using namespace Rcpp;
using namespace arma;

//' @name dec2binar
//' @title Transform a decimal to a binary array
//' @description Transforms a decimal to a binary array with an additional
//'   element saving the number of ones in the array. Helper function for
//'   \code{get_permanent_Cpp}.
//'
//' @param n the index determining the submatrix of the square matrix A the permanent is calculated for.
//' @param dim the number of rows (or columns) in the square matrix A the permanent is calculated for.
//' @details Helper function for \code{get_permanent_Cpp}, which calculates the
//'   permanent of a square matrix using the Ryser's formula.
//'   This function is an Rcpp implementation of the C++ code written by "hanzzoid"
//'   on the "Code Project" website
//'   \url{https://www.codeproject.com/Articles/21282/Compute-Permanent-of-a-Matrix-with-Ryser-s-Algorit}.
//'   Modified code reproduced under The Code Project Open licence
//'   (\url{https://www.codeproject.com/info/cpol10.aspx}).
//' @return Returns a characteric string that allows to select the submatrix.
//' @keywords internal
std::vector<int> dec2binarr(long n, int dim)
{
  // note: res[dim] will save the sum res[0]+...+res[dim-1]
  std::vector<int> res(dim + 1, 0);
  int pos = dim - 1;

  // note: this will crash if dim < log_2(n)...
  while (n > 0)
  {
    res[pos] = n % 2;
    res[dim] += res[pos];
    n = n / 2; // integer division
    pos--;
  }

  return res;
}

//' @name get_permanent_Cpp
//' @title Calculate the permanent of a matrix
//' @description Calculates the permanent of a square matrix using the Ryser's formula.
//'
//' @param A a square matrix
//' @details Calculates the permanent of a square matrix using Ryser's formula.
//'   This function is an Rcpp implementation of the C++ code written by "hanzzoid"
//'   on the "Code Project" website
//'   \url{https://www.codeproject.com/Articles/21282/Compute-Permanent-of-a-Matrix-with-Ryser-s-Algorit}.
//'   Modified code reproduced under The Code Project Open licence
//'   (\url{https://www.codeproject.com/info/cpol10.aspx}).
//' @return Returns the permanent of the matrix A (a real number)
//' @keywords internal
// [[Rcpp::export]]
double get_permanent_Cpp(const arma::mat& A) // expects n by n matrix
{
  int n = A.n_rows;
  double sum = 0;
  double rowsumprod, rowsum;
  std::vector<int> chi(n + 1);
  double C = std::pow(2.0, n);

  // loop all 2^n submatrices of A
  for (int k = 1; k < C; k++)
  {
    rowsumprod = 1;
    chi = dec2binarr(k, n); // characteristic vector

    // loop columns of submatrix #k
    for (int m = 0; m < n; m++)
    {
      rowsum = 0;

      // loop rows and compute rowsum
      for (int p = 0; p < n; p++)
        rowsum += chi[p] * A(m, p);

      // update product of rowsums
      rowsumprod *= rowsum;
    }

    sum += std::pow(-1, n - chi[n]) * rowsumprod;
  }

  return sum;
}



// Old implementation that sometimes returns zero for large matrices, and is slower
// // [[Rcpp::depends(RcppArmadillo)]]
// #include <RcppArmadillo.h>
// #include <bitset>
// using namespace Rcpp;
//
// //' @name get_permanent_Cpp
//  //' @title Calculate the permanent of a matrix
//  //' @description Calculates the permanent of a square matrix using the Ryser's formula.
//  //'
//  //' @param A a square matrix
//  //' @details Calculates the permanent of a square matrix using Ryser's formula.
//  //' @return Returns the permanent of the matrix A (a real number)
//  //' @keywords internal
//  // [[Rcpp::export]]
//  long double get_permanent_Cpp(const arma::mat& A) {
//    int n = A.n_rows;
//    long double perm = 0;
//
//    for (int subset = 1; subset < (1 << n); subset++) {
//      long double prod = 1;
//
//      for (int row = 0; row < n; row++) {
//        long double sum = 0;
//        for (int col = 0; col < n; col++) {
//          if (subset & (1 << col)) {
//            sum += A(row, col);
//          }
//        }
//        prod *= sum;
//      }
//
//      int bits = std::bitset<32>(subset).count(); // count the number of bits set to 1 in subset
//      perm += std::pow(-1, n - bits) * prod;
//    }
//
//    return perm;
//  }
