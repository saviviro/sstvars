context("jointSpectralRadius")
library(sstvars)


test_that("bound_jsr_G works correctly", {
  set.seed(1); S1 <- array(rnorm(3*3*2), dim=c(3, 3, 2))
  set.seed(2); S2 <- array(rnorm(2*2*2), dim=c(2, 2, 2))
  set.seed(3); S3 <- array(rnorm(2*2*3), dim=c(2, 2, 3))
  expect_equal(bound_jsr_G(S1, epsilon=0.01, adaptive_eps=FALSE, print_progress=FALSE), c(1.670107, 1.696743), tol=1e-4)
  expect_equal(bound_jsr_G(S2, epsilon=0.05, adaptive_eps=TRUE, print_progress=FALSE), c(1.567845, 1.618631), tol=1e-4)
  expect_equal(bound_jsr_G(S3, epsilon=0.02, adaptive_eps=FALSE, print_progress=FALSE), c(1.524043, 1.549635), tol=1e-4)
})

# p=1, M=1, d=2
theta_112relg <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103, 0.601786, -0.002945, 0.067224)
mod112relg <- STVAR(data=gdpdef, p=1, M=1, params=theta_112relg, weight_function="relative_dens")

# p=1, M=2, d=2
theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863, -0.3124, 0.505251, 0.298483, 0.030096,
                   -0.176925, 0.838898, 0.310863, 0.007512, 0.018244, 0.949533, -0.016941, 0.121403, 0.573269)
mod122relg <- STVAR(p=1, M=2, d=2, params=theta_122relg, weight_function="relative_dens")

# p=2, M=2, d=2
theta_222relg <- c(0.356914, 0.107436, 0.356386, 0.08633, 0.13996, 0.035172, -0.164575, 0.386816, 0.451675, 0.013086,
                   0.227882, 0.336084, 0.239257, 0.024173, -0.021209, 0.707502, 0.063322, 0.027287, 0.009182, 0.197066,
                   0.205831, 0.005157, 0.025877, 1.092094, -0.009327, 0.116449, 0.592446)
mod222relg <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222relg, weight_function="relative_dens")

# p=1, M=2, d=3, usamone
theta_123relg <- c(0.10741, 0.13813, -0.12092, 3.48957, 0.60615, 0.45646, 0.87227, -0.01595, 0.14124,
                   -0.08611, 0.61865, 0.34311, -0.02047, 0.025, 0.97548, 0.74976, 0.02187, 0.29213,
                   -1.55165, 0.58245, -0.00696, -0.07261, 0.02021, 0.96883, 0.66149, 0.02279, 0.09207,
                   0.05544, 0.00212, 0.12708, 0.78618, 0.00922, 0.42627, 0.23765, 0.25386, 3.40834, 0.77357)
mod123relg <- STVAR(data=usamone, p=1, M=2, params=theta_123relg, weight_function="relative_dens")


test_that("bound_JSR works correctly", {
  # Lower and upper bound by the Gripenberg's (1996) branch-and-bound method
  expect_equal(bound_JSR(mod112relg, epsilon=0.01, adaptive_eps=TRUE, print_progress=FALSE), c(0.8919073, 0.9022341), tol=1e-4)
  expect_equal(bound_JSR(mod222relg, epsilon=0.02, adaptive_eps=FALSE, print_progress=FALSE), c(0.9205818, 0.9421805), tol=1e-4)
  expect_equal(bound_JSR(mod123relg, epsilon=0.3, adaptive_eps=TRUE, print_progress=FALSE), c(0.9609557, 1.3269630), tol=1e-4)
})

