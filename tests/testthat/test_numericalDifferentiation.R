context("numericalDifferentiation")
library(sstvars)

foo1 <- function(x) x^2
foo2 <- function(x, a=1, b=1) a*x[1]^2 - b*x[2]^2
foo3 <- function(x) x[1]^2 + log(x[2]) - x[3]^3

test_that("calc_gradient works correctly", {
  expect_equal(calc_gradient(x=0, fn=foo1), 0, tolerance=1e-5)
  expect_equal(calc_gradient(x=1, fn=foo1), 2, tolerance=1e-5)
  expect_equal(calc_gradient(x=-2, fn=foo1), -4, tolerance=1e-5)

  expect_equal(calc_gradient(x=c(0, 0), fn=foo2), c(0, 0), tolerance=1e-5)
  expect_equal(calc_gradient(x=c(0, 0), fn=foo2, a=2, b=3), c(0, 0), tolerance=1e-5)
  expect_equal(calc_gradient(x=c(1, 2), fn=foo2), c(2, -4), tolerance=1e-5)
  expect_equal(calc_gradient(x=c(1, 2), fn=foo2, a=2, b=3), c(4, -12), tolerance=1e-5)

  expect_equal(calc_gradient(x=c(1, 2, 3), fn=foo3), c(2.0, 0.5, -27.0), tolerance=1e-4)
})


test_that("calc_hessian works correctly", {
  expect_equal(calc_hessian(x=0, fn=foo1), as.matrix(2), tolerance=1e-4)
  expect_equal(calc_hessian(x=1, fn=foo1), as.matrix(2), tolerance=1e-4)
  expect_equal(calc_hessian(x=-2, fn=foo1), as.matrix(2), tolerance=1e-4)

  expect_equal(calc_hessian(x=c(0, 0), fn=foo2), diag(c(2, -2)), tolerance=1e-4)
  expect_equal(calc_hessian(x=c(0, 0), fn=foo2, a=2, b=3), diag(c(4, -6)), tolerance=1e-4)
  expect_equal(calc_hessian(x=c(1, 2), fn=foo2), diag(c(2, -2)), tolerance=1e-4)
  expect_equal(calc_hessian(x=c(1, 2), fn=foo2, a=2, b=3), diag(c(4, -6)), tolerance=1e-4)

  expect_equal(calc_hessian(x=c(1, 2, 3), fn=foo3), diag(c(1.99998, -0.2500222, -18.00002)), tolerance=1e-4)
})
