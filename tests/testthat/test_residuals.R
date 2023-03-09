context("residuals")
library(sstvars)

# p=1, M=1, d=2
theta_112relg <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103, 0.601786, -0.002945, 0.067224)

# p=1, M=2, d=2
theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863, -0.3124, 0.505251, 0.298483, 0.030096,
                   -0.176925, 0.838898, 0.310863, 0.007512, 0.018244, 0.949533, -0.016941, 0.121403, 0.573269)

# p=2, M=2, d=2
theta_222relg <- c(0.356914, 0.107436, 0.356386, 0.08633, 0.13996, 0.035172, -0.164575, 0.386816, 0.451675, 0.013086,
                   0.227882, 0.336084, 0.239257, 0.024173, -0.021209, 0.707502, 0.063322, 0.027287, 0.009182, 0.197066,
                   0.205831, 0.005157, 0.025877, 1.092094, -0.009327, 0.116449, 0.592446)


test_that("get_residuals works correctly", {
  # Relative_dens Gausssian STVAR
  expect_equal(c(get_residuals(data=gdpdef, p=1, M=1, params=theta_112relg, weight_function="relative_dens",
                               standardize=TRUE)[c(1, 2, 113, 243),]),
               c(1.4266258, -1.5593422, 0.4003396, -0.4257862, 7.4896578, -0.7203759, 1.6478346, 0.3644909),
               tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=1, M=1, params=theta_112relg, weight_function="relative_dens",
                               standardize=FALSE)[c(1, 13, 200, 243),]),
               c(1.08538455, -0.18623807, -1.15119664, -0.33133796, 1.93771162, 0.33579095, -1.39753618, 0.09570955),
               tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=1, M=2, params=theta_122relg, weight_function="relative_dens",
                               standardize=TRUE)[c(1, 3, 213, 243),]),
               c(1.5559288, -0.5781259, -0.5039419, -0.5981072, 9.0972135, -0.8712142, -0.5371637, 0.5387503),
               tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=1, M=2, params=theta_122relg, weight_function="relative_dens",
                               standardize=FALSE)[c(1, 20, 242, 243),]),
               c(1.07689397, 1.40168826, -0.08105729, -0.34307860, 1.88952878, 1.34462994, 0.15957893, 0.07570420),
               tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222relg, weight_function="relative_dens",
                               standardize=TRUE)[c(1, 2, 210, 242),]),
               c(-1.31052809, -0.35054679, 0.79660833, -0.85379617, -0.97656565, -0.47672003, 2.83244577, 0.08922395),
               tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222relg, weight_function="relative_dens",
                               standardize=FALSE)[c(1, 54, 150, 242),]),
               c(-1.146330783, 0.484636084, -0.214896901, -0.409659185, -0.277609777, 0.837360313, 0.610229891, 0.008792963),
               tolerance=1e-3)
})
