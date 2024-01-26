context("get_permanent_Cpp")
library(sstvars)


calc_3x3_permanent <- function(A) {
  A[1,1]*A[2,2]*A[3,3] + A[1,2]*A[2,3]*A[3,1] + A[1,3]*A[2,1]*A[3,2] + A[1,3]*A[2,2]*A[3,1] + A[1,2]*A[2,1]*A[3,3] + A[1,1]*A[2,3]*A[3,2]
}

A1 <- matrix(1:4, nrow=2)
A2 <- matrix(c(1, -2, 0.3, 0.5), nrow=2)
A3 <- matrix(c(0.9, -2, 0.3, -0.5), nrow=2)
A4 <- matrix(c(-0.9, -2, 0.3, -0.5), nrow=2)
A5 <- matrix(1:9, nrow=3)
A6 <- matrix(c(0.4, -0.2, 1.2, 0.7, 0.1, -1.5, -0.7, 0.13, 0.35), nrow=3)
A7 <- matrix(c(-0.1, -0.13, 0.5, -0.6, 0.9, 1.1, 1.4, -0.12, -0.06), nrow=3)
A8 <- matrix(c(0, 0, 0.5, 0, 0.9, 0, -1, 0, 0), nrow=3)
A9 <- matrix(1:16, nrow=4, byrow=TRUE)
A10 <- matrix(1:16, nrow=4)
A11 <- matrix(c(-0.3, -0.12, 3:16), nrow=4)
A12 <- matrix(c(-0.3, -0.12, 3:14, 9999, 99999), nrow=4)
A13 <- matrix(c(0.1, -0.5, 2, 0.13, -0.12, -0.06, -0.7, -0.9, 0.11, 0.56, -0.23, -0.8, -2.0, 0.17, 0.82, 0.15), nrow=4)

test_that("get_permanent_Cpp work correctly", {
  expect_equal(get_permanent_Cpp(A1), 1*4 + 3*2, tolerance=1e-3)
  expect_equal(get_permanent_Cpp(-A1), 10, tolerance=1e-3)
  expect_equal(get_permanent_Cpp(A2), 1*0.5 - 2*0.3, tolerance=1e-3)
  expect_equal(get_permanent_Cpp(A3), 0.9*(-0.5) - 2*0.3, tolerance=1e-3)
  expect_equal(get_permanent_Cpp(A4), -0.9*(-0.5) - 2*0.3, tolerance=1e-3)
  expect_equal(get_permanent_Cpp(A5), 450, tolerance=1e-3)
  expect_equal(get_permanent_Cpp(-A5), -450, tolerance=1e-3)
  expect_equal(get_permanent_Cpp(A6), calc_3x3_permanent(A6), tolerance=1e-3)
  expect_equal(get_permanent_Cpp(-A6), calc_3x3_permanent(-A6), tolerance=1e-3)
  expect_equal(get_permanent_Cpp(A7), calc_3x3_permanent(A7), tolerance=1e-3)
  expect_equal(get_permanent_Cpp(-A7), calc_3x3_permanent(-A7), tolerance=1e-3)
  expect_equal(get_permanent_Cpp(A8), calc_3x3_permanent(A8), tolerance=1e-3)
  expect_equal(get_permanent_Cpp(A9), 55456, tolerance=1e-3)
  expect_equal(get_permanent_Cpp(A10), 55456, tolerance=1e-3)
  expect_equal(get_permanent_Cpp(A11), 33592.3, tolerance=1e-3)
  expect_equal(get_permanent_Cpp(A12), 29273150, tolerance=1e-3)
  expect_equal(get_permanent_Cpp(A13), 2.632122, tolerance=1e-3)
})


