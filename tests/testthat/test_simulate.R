context("simulate.stvar")
library(sstvars)

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

s112 <- simulate(mod112relg, nsim=1, seed=1, init_regime=1)
s112_2 <- simulate(mod112relg, nsim=2, seed=1, init_regime=1)
s122 <- simulate(mod122relg, nsim=5, seed=2, init_values=gdpdef)
s222 <- simulate(mod222relg, nsim=3, seed=3, init_regime=2)
s123 <- simulate(mod123relg, nsim=3, seed=4, init_regime=1)
s123_2 <- simulate(mod123relg, nsim=1, seed=5, init_values=usamone)


test_that("simulate.stvar works correctly", {
  # Relative_dens Gaussian STVAR
  expect_equal(s112$sample[1,], c(-0.07206511, 1.343205), tol=1e-4)
  expect_equal(s112$transition_weights, as.matrix(1), tol=1e-4)
  expect_equal(s112_2$sample[1,], c(-0.07206511, 1.343205), tol=1e-4)
  expect_equal(s112_2$sample[2,], c(0.6932277, 1.0562789), tol=1e-4)
  expect_equal(s112_2$sample[2,], c(0.6932277, 1.0562789), tol=1e-4)
  expect_equal(s112_2$transition_weights, as.matrix(c(1, 1)), tol=1e-4)
  expect_equal(s122$sample[5,], c(2.1049174, 0.4326244), tol=1e-4)
  expect_equal(s122$transition_weights[5,], c(0.94892269, 0.05107731), tol=1e-4)
  expect_equal(s222$sample[3,], c(0.1541768, 0.8947292), tol=1e-4)
  expect_equal(s222$transition_weights[1,], c(0.07983494, 0.92016506), tol=1e-4)
  expect_equal(s123$sample[3,], c(1.346929, 0.875446, 5.316535), tol=1e-4)
  expect_equal(s123$transition_weights[3,], c(0.98754515, 0.01245485), tol=1e-4)
  expect_equal(s123_2$sample[1,], c(1.684557, 2.126412, -2.721774), tol=1e-4)
  expect_equal(s123_2$transition_weights[1,], c(0.0001666309, 0.9998333691), tol=1e-4)
})

