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
               c(1.3968841, -1.5548925, 0.3948874, -0.4266263, -0.6188517, 0.4927395, 0.1613921, 0.1354454),
               tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=1, M=1, params=theta_112relg, weight_function="relative_dens",
                               standardize=FALSE)[c(1, 13, 200, 243),]),
               c(1.08538455, -0.18623807, -1.15119664, -0.33133796, -0.16441838, -0.40278905, -0.28125618, 0.03632955),
               tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=1, M=2, params=theta_122relg, weight_function="relative_dens",
                               standardize=TRUE)[c(1, 3, 213, 243),]),
               c(1.5833935, -0.5878140, -0.5001501, -0.5921047, -1.0392301, -0.1613458, -0.8607716, 0.1433740),
               tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=1, M=2, params=theta_122relg, weight_function="relative_dens",
                               standardize=FALSE)[c(1, 20, 242, 243),]),
               c(1.07689397, 1.40168826, -0.08105729, -0.34307860, -0.21260122, -0.44618006, -0.20047107, 0.01632420),
               tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222relg, weight_function="relative_dens",
                               standardize=TRUE)[c(1, 2, 210, 242),]),
               c(-1.3062112, -0.3485900, 0.8453746, -0.8484317, 0.1119452, -0.1325865, -2.5016219, -0.2633884),
               tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222relg, weight_function="relative_dens",
                               standardize=FALSE)[c(1, 54, 150, 242),]),
               c(-1.14633078, 0.48463608, -0.21489690, -0.40965919, 0.03689022, 0.43521031, 0.10575989, -0.05058704),
               tolerance=1e-3)
})
