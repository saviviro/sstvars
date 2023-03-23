context("predict.stvar")
library(sstvars)

# p=1, M=1, d=2
theta_112relg <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103, 0.601786, -0.002945, 0.067224)
mod112relg <- STVAR(data=gdpdef, p=1, M=1, params=theta_112relg, weight_function="relative_dens")

# p=1, M=2, d=2
theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863, -0.3124, 0.505251, 0.298483, 0.030096,
                   -0.176925, 0.838898, 0.310863, 0.007512, 0.018244, 0.949533, -0.016941, 0.121403, 0.573269)
mod122relg <- STVAR(data=gdpdef, p=1, M=2, d=2, params=theta_122relg, weight_function="relative_dens")

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

set.seed(1); p112 <- predict(mod112relg, nsteps=1, nsim=1, pred_type="mean")
set.seed(2); p122 <- predict(mod122relg, nsteps=3, nsim=3, pred_type="median")
set.seed(3); p222 <- predict(mod222relg, nsteps=2, nsim=10, pred_type="mean", pi=c(0.5))
set.seed(4); p123 <- predict(mod123relg, nsteps=4, nsim=7, pred_type="median", pi=c(0.3, 0.5, 0.7))


test_that("simulate.stvar works correctly", {
  # Relative_dens Gaussian STVAR
  expect_equal(unname(p112$pred[1,]), c(0.2393741, 0.4936729), tol=1e-4)
  expect_equal(unname(p112$pred_ints[1, 3,]), c(0.2393741, 0.4936729), tol=1e-4)
  expect_equal(unname(p112$trans_pred[1,]), c(1), tol=1e-4)
  expect_equal(unname(p112$trans_pred_ints[1, 1, 1]), c(1), tol=1e-4)

  expect_equal(unname(p122$pred[3,]), c(1.4321082, 0.4083084), tol=1e-4)
  expect_equal(unname(p122$pred_ints[3, 2,]), c(1.0924034, 0.2907614), tol=1e-4)

  expect_equal(unname(p222$trans_pred[2,]), c(0.95938141, 0.04061859), tol=1e-4)
  expect_equal(unname(p222$trans_pred_ints[2, 1,]), c(0.95597200, 0.03039205), tol=1e-4)

  expect_equal(unname(p123$pred[4,]), c(2.154655, 1.960712, 1.612034), tol=1e-4)
  expect_equal(unname(p123$pred_ints[4, 5,]), c(2.422710, 2.372279, 7.711421), tol=1e-4)
  expect_equal(unname(p123$trans_pred[4,]), c(0.0002149667, 0.9997850333), tol=1e-4)
  expect_equal(unname(p123$trans_pred_ints[3, 4,]), c(0.1428647, 1.0000000), tol=1e-4)
})

