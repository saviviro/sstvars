context("get_Bt_Cpp")
library(sstvars)

# Test for basic functionality with known output
test_that("get_Bt_Cpp calculates correct weighted sums", {
  all_Omegas <- array(1:18, dim=c(3, 3, 2)) # Two 3x3 matrices
  alpha_mt <- matrix(c(1, 2), nrow=1) # Single row of weights

  # Expected result is the sum of the first matrix plus two times the second matrix
  expected <- matrix(c(21, 24, 27, 30, 33, 36, 39, 42, 45), nrow = 3)
  result <- get_Bt_Cpp(all_Omegas, alpha_mt)
  expect_true(all.equal(result[, , 1], expected))

  alpha_mt <- matrix(c(1, 0.5, 2, 3), nrow=2) # Two rows of weights
  expected <- array(c(21, 24, 27, 30, 33, 36, 39, 42, 45, 30.5, 34.0, 37.5,
                      41.0, 44.5, 48.0, 51.5, 55.0, 58.5), dim=c(3, 3, 2))
  result <- get_Bt_Cpp(all_Omegas, alpha_mt)
  expect_true(all.equal(result, expected))

  all_Omegas <- array(1:12, dim=c(2, 2, 3)) # Three 2x2 matrices
  alpha_mt <- matrix(c(0.5, 0.3, 0.2), nrow=1) # Single row of weights

  expected <- array(c(3.8, 4.8, 5.8, 6.8), dim=c(2, 2, 1))
  result <- get_Bt_Cpp(all_Omegas, alpha_mt)
  expect_true(all.equal(result, expected))

  alpha_mt <-  rbind(c(0.5, 0.3, 0.2), c(0, 0, 1), c(0.1, 0.2, 0.7)) # Three rows of weights
  expected <- array(c(3.8, 4.8, 5.8, 6.8, 9, 10, 11, 12, 7.4, 8.4, 9.4, 10.4), dim=c(2, 2, 3))
  result <- get_Bt_Cpp(all_Omegas, alpha_mt)
  expect_true(all.equal(result, expected))

  alpha_mt <-  rbind(c(0.5, 0.3, 0.2), c(0, 0, 1), c(0.1, 0.2, 0.7), c(0.8, 0, 0.2)) # Four rows of weights
  expected <- array(c(3.8, 4.8, 5.8, 6.8, 9, 10, 11, 12, 7.4, 8.4, 9.4, 10.4, 2.6, 3.6, 4.6, 5.6), dim=c(2, 2, 4))
  result <- get_Bt_Cpp(all_Omegas, alpha_mt)
  expect_true(all.equal(result, expected))
})

# Test for handling of a single matrix and single weight
test_that("get_Bt_Cpp handles single matrix and weight correctly", {
  all_Omegas <- array(1:9, dim=c(3, 3, 1)) # Single 3x3 matrix
  alpha_mt <- matrix(1, nrow=1) # Single weight

  # Expected result is the same as the input matrix
  expected <- matrix(1:9, nrow=3)
  result <- get_Bt_Cpp(all_Omegas, alpha_mt)
  expect_true(all.equal(result[, , 1], expected))
})

# Test for output structure
test_that("get_Bt_Cpp output has correct dimensions", {
  all_Omegas <- array(rnorm(27), dim=c(3, 3, 3)) # Three 3x3 matrices
  alpha_mt <- matrix(rnorm(6), nrow=2) # Two sets of weights

  result <- get_Bt_Cpp(all_Omegas, alpha_mt)
  expect_equal(dim(result), c(3, 3, 2))
})
