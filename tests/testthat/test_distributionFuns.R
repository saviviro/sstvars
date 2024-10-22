context("distributionFuns")
library(sstvars)

test_that("skewed_t_dens when lambda = 0", {
  expect_equal(skewed_t_dens(y=-5, nu=3, lambda=0), 0.0009417452, tolerance = 1e-4)
  expect_equal(skewed_t_dens(y=-3, nu=5, lambda=0), 0.007657346, tolerance = 1e-4)
  expect_equal(skewed_t_dens(y=0, nu=5, lambda=0), 0.4900701, tolerance = 1e-4)
  expect_equal(skewed_t_dens(y=1, nu=5, lambda=0), 0.2067483, tolerance = 1e-4)
})

test_that("skewed_t_dens gives correct results", {
  expect_equal(skewed_t_dens(y=-5.1, nu=2.1, lambda=0.1), 0.0002128741, tolerance = 1e-4)
  expect_equal(skewed_t_dens(y=-3.2, nu=7, lambda=-0.2), 0.008096013, tolerance = 1e-4)
  expect_equal(skewed_t_dens(y=0, nu=4.3, lambda=0.99), 0.4197044, tolerance = 1e-4)
  expect_equal(skewed_t_dens(y=1.4, nu=12, lambda=-0.3), 0.1404645, tolerance = 1e-4)
  expect_equal(skewed_t_dens(y=2.4, nu=30, lambda=0.33), 0.03230915, tolerance = 1e-4)
})

test_that("skewed_t_dens returns vector of same length as y", {
  expect_equal(length(skewed_t_dens(y=rnorm(100), nu=4, lambda=0.3)), 100)
  expect_equal(length(skewed_t_dens(y=rnorm(13), nu=4, lambda=0.3)), 13)
})

test_that("skewed_t_dens returns zero for infinite y values", {
  y <- c(-Inf, -1, 0, 1, Inf)
  result <- skewed_t_dens(y=y, nu=5, lambda=0.5)
  expect_equal(result[which(y == -Inf)], 0)
  expect_equal(result[which(y == Inf)], 0)
})

test_that("skewed_t_dens returns non-negative densities", {
  result <- skewed_t_dens(y=rnorm(100), nu=5, lambda=0.5)
  expect_true(all(result >= 0))
})


#
test_that("stand_t_dens gives correct results", {
  expect_equal(stand_t_dens(y=-5.1, nu=2.1), 0.0002923304, tolerance = 1e-4)
  expect_equal(stand_t_dens(y=-3.2, nu=7), 0.005277828, tolerance = 1e-4)
  expect_equal(stand_t_dens(y=0, nu=4.3), 0.5149267, tolerance = 1e-4)
  expect_equal(stand_t_dens(y=1.4, nu=12), 0.1337244, tolerance = 1e-4)
  expect_equal(stand_t_dens(y=2.4, nu=30), 0.02254161, tolerance = 1e-4)
})

test_that("stand_t_dens returns vector of same length as y", {
  expect_equal(length(stand_t_dens(y=rnorm(100), nu=4)), 100)
  expect_equal(length(stand_t_dens(y=rnorm(13), nu=4)), 13)
})

test_that("stand_t_dens returns zero for infinite y values", {
  y <- c(-Inf, -1, 0, 1, Inf)
  result <- stand_t_dens(y=y, nu=5)
  expect_equal(result[which(y == -Inf)], 0)
  expect_equal(result[which(y == Inf)], 0)
})

test_that("stand_t_dens returns non-negative densities", {
  result <- stand_t_dens(y=rnorm(100), nu=4)
  expect_true(all(result >= 0))
})

test_that("bounding_const_M returns a positive numeric scalar", {
  bM <- bounding_const_M(nu=10, lambda=-0.5)
  expect_true(is.numeric(bM))
  expect_length(bM, 1)
  expect_gt(bM, 0)
})

test_that("bounding_const_M returns expected value for known inputs", {
  # Since lambda = 0, the skewed t-distribution reduces to the standard t-distribution
  # Therefore, M should be close to 1.1 for nu < 3 (due to the safety margin)
  expect_true(abs(bounding_const_M(nu=2.8, lambda=0) - 1.1) < 0.01)
})

test_that("bounding_const_M works correctly", {
  expect_equal(bounding_const_M(nu=2.1, lambda=0.1), 1.411976, tolerance=1e-3)
  expect_equal(bounding_const_M(nu=2.000001, lambda=-0.999999999), 2.199342, tolerance=1e-3)
  expect_equal(bounding_const_M(nu=2.1, lambda=0.999999999), 5.389264, tolerance=1e-3)
  expect_equal(bounding_const_M(nu=10000000, lambda=0.999999999), 6.146267, tolerance=1e-3)
  expect_equal(bounding_const_M(nu=10000000, lambda=-0.999999999), 6.146267, tolerance=1e-3)
  expect_equal(bounding_const_M(nu=100, lambda=0.9999), 6.220674, tolerance=1e-3)
  expect_equal(bounding_const_M(nu=10, lambda=0.999), 5.798037, tolerance=1e-3)
  expect_equal(bounding_const_M(nu=7.1, lambda=-0.3), 2.252657, tolerance=1e-3)
  expect_equal(bounding_const_M(nu=9.2, lambda=0.8), 4.06109, tolerance=1e-3)
})


test_that("generate_skewed_t returns a numeric vector of length n", {
  set.seed(1); samples <- generate_skewed_t(n=113, nu=5, lambda=0.5)
  expect_true(is.numeric(samples))
  expect_length(samples,113)
})

test_that("generate_skewed_t works correctly", {
  set.seed(2); res1 <- generate_skewed_t(n=1, nu=2.1, lambda=0.1)
  set.seed(3); res2 <- generate_skewed_t(n=1, nu=10, lambda=0.9)
  set.seed(4); res3 <- generate_skewed_t(n=1, nu=5, lambda=-0.3, bc_M=10)
  expect_equal(res1, 0.28912, tolerance=1e-3)
  expect_equal(res2, -0.5275574, tolerance=1e-3)
  expect_equal(res3, 0.07797383, tolerance=1e-3)
})
