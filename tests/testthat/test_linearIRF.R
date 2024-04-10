context("linearIRF")
library(sstvars)

# p=1, M=1, d=2, linear Gaussian VAR model, shocks identified recursively.
theta_112relg <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103, 0.601786, -0.002945, 0.067224)
mod112 <- STVAR(data=gdpdef, p=1, M=1, params=theta_112relg, identification="recursive")


# p=1, M=2, d=2, Gaussian STVAR with relative dens weight function,
# shocks identified recursively.
theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863, -0.3124, 0.505251, 0.298483,
                   0.030096, -0.176925, 0.838898, 0.310863, 0.007512, 0.018244, 0.949533, -0.016941, 0.121403, 0.573269)
mod122 <- STVAR(data=gdpdef, p=1, M=2, params=theta_122relg, identification="recursive")

# p=3, M=2, d=3, Students't logistic STVAR model with the first lag of the second
# variable as the switching variable. Autoregressive dynamics restricted linear,
# but the volatility regime varies in time, allowing the shocks to be identified
# by conditional heteroskedasticity.
theta_322 <- c(0.7575, 0.6675, 0.2634, 0.031, -0.007, 0.5468, 0.2508, 0.0217, -0.0356, 0.171, -0.083, 0.0111,
               -0.1089, 0.1987, 0.2181, -0.1685, 0.5486, 0.0774, 5.9398, 3.6945, 1.2216, 8.0716, 8.9718)
mod322 <- STVAR(data=gdpdef, p=3, M=2, params=theta_322, weight_function="logistic",
  weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2),
  AR_constraints=rbind(diag(3*2^2), diag(3*2^2)), identification="heteroskedasticity",
  parametrization="mean")

# p=1, M=2, d=2, weight_function="exogenous", cond_dist="ind_Student", weightfun_pars=twmat, AR_constraints=C_122,
tw1 <- c(0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
         1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1,
         1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
twmat <- cbind(tw1, 1 - tw1)
C_122 <- rbind(diag(1*2^2), diag(1*2^2))
params122exocit <- c(-0.25329555, 0.10340501, 0.73585978, 0.05335907, 0.18282378, 0.03478017, -0.01629061, 0.87424509, 0.83939713,
                     0.20970725, 0.47256152, -0.23555538, 0.59277160, 0.09901218, -0.29605914, 0.23895480, 6.72161948, 4.10403450)

mod122exocit <- STVAR(data=gdpdef, p=1, M=2, d=2, params=params122exocit, weight_function="exogenous", cond_dist="ind_Student",
                      weightfun_pars=twmat, AR_constraints=C_122)

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="ind_Student", identification="non-Gaussianity"
params222logit <- c(0.428924, 0.153884, 1.455217, 0.264533, 0.314288, 0.124903, -2.633578, -0.313929, 0.293621, 0.036355, -0.273429,
                    0.053141, 0.210305, -0.023328, -0.516816, 0.501481, 0.190278, 0.035209, -0.054981, 0.340776, 1.384283, -0.159169,
                    0.449031, 0.522074, -1.213738, -0.086115, 0.353828, -0.442951, 0.313206, 2.807629, 4.924918, 7.722209)
mod222logit <- STVAR(data=gdpdef, p=2, M=2, params=params222logit, weight_function="logistic", weightfun_pars=c(2, 1),
                     cond_dist="ind_Student", identification="non-Gaussianity")


test_that("linear_IRF works correctly", {
  irf1 <- linear_IRF(stvar=mod112, N=4, regime=1, robust_method="none", seed=1)
  irf2 <- linear_IRF(stvar=mod122, N=7, regime=1, scale=cbind(c(1, 1, 0.3), c(2, 2, 0.5)))
  irf3 <- linear_IRF(stvar=mod122, N=1, regime=2)
  irf4 <- linear_IRF(stvar=mod322, N=2, regime=1, ncores=1, which_cumulative=2, scale=c(2, 1, 0.4), seed=1)
  irf5 <- linear_IRF(stvar=mod122exocit, N=3, regime=2, ncores=1, which_cumulative=1, scale=c(2, 1, 0.4), seed=1)
  irf6 <- linear_IRF(stvar=mod222logit, N=2, regime=1, ncores=1, which_cumulative=1:2, scale=c(1, 1, 0.4), seed=1)

  expect_equal(c(irf1$point_est), c(0.775748671, -0.003796333, 0.000000000, 0.259248120, 0.224370424, 0.013480020, -0.037337951,
                                    0.232572266, 0.062795255, 0.016976837, -0.044268958, 0.207828543, 0.015672992, 0.016596836,
                                    -0.042705043, 0.185480007, 0.002131723, 0.015230225, -0.039035088, 0.165465110), tolerance=1e-3)
  expect_equal(c(irf2$point_est)[c(1, 3, 20)], c(0.3000000, 0.0000000, 0.0327313), tolerance=1e-3)
  expect_equal(c(irf3$point_est), c(0.97443984, -0.01738537, 0.00000000, 0.34799533, 0.29392963, 0.01474219, -0.06156907,
                                    0.29193259), tolerance=1e-3)
  expect_equal(c(irf4$point_est), c(0.21810000, -0.16850000, 0.40000000, 0.05643456, 0.05862704, -0.25387470, 0.10496496, 0.09969298,
                                    0.07673807, -0.32282088, 0.12565589, 0.14493090), tolerance=1e-3)
  expect_equal(c(irf5$point_est), c(0.59277160, 0.09901218, 0.40000000, -0.32284739, 0.69953138, 0.10717761, 0.47838889, -0.26833568,
                                    0.71730361, 0.09741262, 0.49709160, -0.23186477, 0.71896589, 0.08578063, 0.50428812, -0.20205615),
               tolerance=1e-3)
  expect_equal(c(irf6$point_est), c(0.4000000, -0.0459932, 0.4490310, 0.5220740, 0.6468419, 0.0184066, -0.7847665, 0.4142652, 0.6848437,
                                    0.0411188, -0.8995165, 0.3380725), tolerance=1e-3)
})

