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


# p=1, M=2, d=2
phi10_122 <- c(0.55, 0.11)
A11_122 <- matrix(c(0.34, 0.05, -0.01, 0.72), nrow=2, byrow=FALSE)
Omega1_122 <- matrix(c(0.58, 0.01, 0.01, 0.06), nrow=2, byrow=FALSE)

phi20_122 <- c(0.17, 0.25)
A21_122 <- A11_122
Omega2_122 <- matrix(c(0.50, -0.01, -0.01, 0.20), nrow=2, byrow=FALSE)

# p=2, M=2, d=2
phi10_222 <- c(0.36, 0.12)
A11_222 <- matrix(c(0.22, 0.06, -0.15, 0.39), nrow=2, byrow=FALSE)
A12_222 <- matrix(c(0.41, -0.01, 0.08, 0.3), nrow=2, byrow=FALSE)
Omega1_222 <- matrix(c(0.21, 0.01, 0.01, 0.03), nrow=2, byrow=FALSE)

phi20_222 <- c(0.48, 0.07)
A21_222 <- matrix(c(0.22, 0.02, -0.12, 0.72), nrow=2, byrow=FALSE)
A22_222 <- matrix(c(0.09, 0.03, 0.04, 0.19), nrow=2, byrow=FALSE)
Omega2_222 <- matrix(c(1.10, 0.01, 0.01, 0.11), nrow=2, byrow=FALSE)

rpars122_1 <-c(phi20_122, vec(A21_122))
rpars122_2 <-c(phi10_122, vec(A11_122))

rpars122t_1 <- c(rpars122_1, 1001)
rpars122t_2 <- c(rpars122_2, 1010)

rpars332_1 <- c(phi10_122, vec(A11_122), vec(A11_122), vec(A11_122))
rpars332_2 <- c(phi20_222, vec(A21_222), vec(A22_222), vec(A22_222))

rpars332t_1 <- c(rpars332_1, 3)
rpars332t_2 <- c(rpars332_2, 13)

test_that("regime_distance works correctly", {
  expect_equal(regime_distance(regime_pars1=rpars122_1, regime_pars2=rpars122_2), 0.4049691, tol=1e-4)
  expect_equal(regime_distance(regime_pars1=rpars122t_1, regime_pars2=rpars122t_2), 0.4049701, tol=1e-4)
  expect_equal(regime_distance(regime_pars1=rpars332_1, regime_pars2=rpars332_2), 0.8523497, tol=1e-4)
  expect_equal(regime_distance(regime_pars1=rpars332t_1, regime_pars2=rpars332t_2), 0.8691375, tol=1e-4)
})


test_that("get_new_start works correctly", {
  expect_equal(get_new_start(y_start=c(1999, 1), y_freq=4, steps_forward=2), c(1999, 3))
  expect_equal(get_new_start(y_start=c(1999, 1), y_freq=4, steps_forward=5), c(2000, 2))
  expect_equal(get_new_start(y_start=c(1989, 1), y_freq=4, steps_forward=15), c(1992, 4))

  expect_equal(get_new_start(y_start=c(1999, 12), y_freq=12, steps_forward=1), c(2000, 1))
  expect_equal(get_new_start(y_start=c(1999, 12), y_freq=12, steps_forward=13), c(2001, 1))
  expect_equal(get_new_start(y_start=c(1999, 12), y_freq=12, steps_forward=14), c(2001, 2))

  expect_equal(get_new_start(y_start=c(1999, 50), y_freq=52, steps_forward=2), c(1999, 52))
  expect_equal(get_new_start(y_start=c(1999, 50), y_freq=52, steps_forward=4), c(2000, 2))
  expect_equal(get_new_start(y_start=c(2000, 51), y_freq=52, steps_forward=52), c(2001, 51))

  expect_equal(get_new_start(y_start=c(1999, 50), y_freq=250, steps_forward=200), c(1999, 250))
  expect_equal(get_new_start(y_start=c(1999, 200), y_freq=250, steps_forward=101), c(2000, 51))
  expect_equal(get_new_start(y_start=c(2000, 50), y_freq=250, steps_forward=500), c(2002, 50))
})

