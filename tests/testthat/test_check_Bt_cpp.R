context("check_Bt_Cpp")
library(sstvars)

# Test with invertible matrices
test_that("check_Bt_Cpp works with invertible matrices", {
  all_Omegas <- array(c(34, 50, 12, 21, 40, 13, 22, 10, 23, 31, 46, 31, 34, 29,
                        16, 46, 36, 4, 48, 20, 47, 31, 44, 47, 27, 25, 20),
                      dim = c(3, 3, 3)) # Cube of 3 matrices, each 3x3
  alpha_mt <- matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), ncol=3) # Identity weights
  posdef_tol <- 1e-5
  expect_true(check_Bt_Cpp(all_Omegas=all_Omegas, alpha_mt=alpha_mt, posdef_tol=posdef_tol))
})

# Test with non-invertible matrices (determinant close to 0)
test_that("check_Bt_Cpp works with non-invertible matrices", {
  all_Omegas <- array(c(1, 2, 3, 2, 4, 6, 3, 6, 9), dim=c(3, 3, 1)) # Single 3x3 matrix, determinant = 0
  alpha_mt <- matrix(1, ncol=1) # Single weight
  posdef_tol <- 1e-5
  expect_false(check_Bt_Cpp(all_Omegas=all_Omegas, alpha_mt=alpha_mt, posdef_tol=posdef_tol))
})

# Test with a positive tolerance check
test_that("check_Bt_Cpp respects positive tolerance", {
  all_Omegas <- array(c(2, 0, 0, 0, 0.5, 0, 0, 0, 0.25), dim=c(3, 3, 1)) # Single 3x3 matrix, determinant = 0.25
  alpha_mt <- matrix(1, ncol=1) # Single weight
  posdef_tol <- 0.5 # Larger than determinant
  expect_false(check_Bt_Cpp(all_Omegas=all_Omegas, alpha_mt=alpha_mt, posdef_tol=posdef_tol))
})

test_that("Weighted sum of invertible matrices is non-invertible", {
  # Two invertible matrices
  all_Omegas <- array(0, dim=c(2, 2, 2))
  all_Omegas[, , 1] <- diag(2) # Identity matrix
  all_Omegas[, , 2] <- matrix(c(0, 2, 2, 0), nrow=2) # Another invertible matrix

  # Weights designed to produce a singular matrix when applied
  alpha_mt <- matrix(c(2, 2, 1, 1), ncol=2)

  # Tolerance
  posdef_tol <- 1e-5

  # Expect the function to return FALSE, indicating the weighted sum is not invertible
  expect_false(check_Bt_Cpp(all_Omegas=all_Omegas, alpha_mt=alpha_mt, posdef_tol=posdef_tol))
})
