context("miscellaneous")
library(sstvars)

set.seed(1); data2 <- cbind(gdpdef, round(rnorm(nrow(gdpdef)), 3))

# p=1, M=1, d=2
theta_112relg <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103, 0.601786, -0.002945, 0.067224)

# p=1, M=2, d=2
theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863, -0.3124, 0.505251, 0.298483, 0.030096,
                   -0.176925, 0.838898, 0.310863, 0.007512, 0.018244, 0.949533, -0.016941, 0.121403, 0.573269)

# p=2, M=2, d=2
theta_222relg <- c(0.356914, 0.107436, 0.356386, 0.08633, 0.13996, 0.035172, -0.164575, 0.386816, 0.451675, 0.013086,
                   0.227882, 0.336084, 0.239257, 0.024173, -0.021209, 0.707502, 0.063322, 0.027287, 0.009182, 0.197066,
                   0.205831, 0.005157, 0.025877, 1.092094, -0.009327, 0.116449, 0.592446)



test_that("get_minval works correctly", {
  expect_equal(get_minval(gdpdef), -99999, tolerance=1e-6)
  expect_equal(get_minval(rbind(gdpdef, gdpdef)), -99999, tolerance=1e-6)
  expect_equal(get_minval(data2), -999999, tolerance=1e-6)
  expect_equal(get_minval(rbind(data2, data2)), -999999, tolerance=1e-6)
})

test_that("get_IC works correctly", {
  expect_equal(get_IC(loglik=-993, npars=12, T_obs=200)$AIC, 10.05, tolerance=1e-2)
  expect_equal(get_IC(loglik=-993, npars=12, T_obs=200)$HQIC, 10.13009, tolerance=1e-3)
  expect_equal(get_IC(loglik=-993, npars=12, T_obs=200)$BIC, 10.2479, tolerance=1e-3)
  expect_equal(get_IC(loglik=-1013, npars=23, T_obs=240)$AIC, 8.633333, tolerance=1e-3)
  expect_equal(get_IC(loglik=-1013, npars=23, T_obs=240)$HQIC, 8.767734, tolerance=1e-3)
  expect_equal(get_IC(loglik=-1013, npars=23, T_obs=240)$BIC, 8.966895, tolerance=1e-3)
})
