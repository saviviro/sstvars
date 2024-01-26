context("diagTests")
library(sstvars)

# p=2, M=2, d=2, STVAR with relative_dens weight function:
theta_222relg <- c(0.357, 0.107, 0.356, 0.086, 0.14, 0.035, -0.165, 0.387, 0.452, 0.013, 0.228, 0.336, 0.239, 0.024,
                   -0.021, 0.708, 0.063, 0.027, 0.009, 0.197, 0.206, 0.005, 0.026, 1.092, -0.009, 0.116, 0.592)
mod222relg <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222relg, weight_function="relative_dens")


# p=3, M=2, d=2, Student's t Threhold VAR with the first lag of the second variable as the switching variable:
theta_322thres <- c(0.527, 0.039, 1.922, 0.154, 0.284, 0.053, 0.033, 0.453, 0.291, 0.024, -0.108, 0.153, -0.108,
                    0.003, -0.128, 0.219, 0.195, -0.03, -0.893, 0.686, 0.047, 0.016, 0.524, 0.068, -0.025, 0.044,
                    -0.435, 0.119, 0.359, 0.002, 0.038, 1.252, -0.041, 0.151, 1.196, 12.312)
mod322thres <- STVAR(data=gdpdef, p=3, M=2, d=2, params=theta_322thres, weight_function="threshold",
                     weightfun_pars=c(2, 1), cond_dist="Student")

test_that("Portmanteau_test works correctly", {
  # p=2, M=2, d=2, relative dens
  expect_error(Portmanteau_test(mod222relg, nlags=2))

  actest222relg <- Portmanteau_test(mod222relg, nlags=3, which_test="autocorr")
  expect_equal(actest222relg$test_stat, 13.88495, tolerance=1e-3)
  expect_equal(actest222relg$p_value, 0.00767145, tolerance=1e-3)
  expect_equal(actest222relg$df, 4, tolerance=1e-3)

  chtest222relg <- Portmanteau_test(mod222relg, nlags=31, which_test="het.sked")
  expect_equal(chtest222relg$test_stat, 140.5294, tolerance=1e-3)
  expect_equal(chtest222relg$p_value, 0.06030631, tolerance=1e-3)
  expect_equal(chtest222relg$df, 116, tolerance=1e-3)

  # p=3, M=2, d=2, Student's t threshold
  expect_error(Portmanteau_test(mod322thres, nlags=3, which_test="het.sked"))
  actest322thres <- Portmanteau_test(mod322thres, nlags=20, which_test="autocorr")
  expect_equal(actest322thres$test_stat, 80.655, tolerance=1e-3)
  expect_equal(actest322thres$p_value, 0.1399178, tolerance=1e-3)
  expect_equal(actest322thres$df, 68, tolerance=1e-3)

  chtest322thres <- Portmanteau_test(mod322thres, nlags=4, which_test="het.sked")
  expect_equal(chtest322thres$test_stat, 22.94353, tolerance=1e-3)
  expect_equal(chtest322thres$p_value, 0.0001299585, tolerance=1e-3)
  expect_equal(chtest322thres$df, 4, tolerance=1e-3)
})

