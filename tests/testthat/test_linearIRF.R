context("linearIRF")
library(sstvars)

# p=1, M=1, d=2, linear Gaussian VAR model, shocks identified recursively.
theta_112relg <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103, 0.601786, -0.002945, 0.067224)
mod112 <- STVAR(data=gdpdef, p=1, M=1, params=theta_112relg)


# p=1, M=2, d=2, Gaussian STVAR with relative dens weight function,
# shocks identified recursively.
theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863, -0.3124, 0.505251, 0.298483,
                   0.030096, -0.176925, 0.838898, 0.310863, 0.007512, 0.018244, 0.949533, -0.016941, 0.121403, 0.573269)
mod122 <- STVAR(data=gdpdef, p=1, M=2, params=theta_122relg)

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

test_that("linear_IRF works correctly", {
  irf1 <- linear_IRF(stvar=mod112, N=4, regime=1, ci=0.90, bootstrap_reps=3, robust_method="none", seed=1)
  irf2 <- linear_IRF(stvar=mod122, N=7, regime=1, scale=cbind(c(1, 1, 0.3), c(2, 2, 0.5)))
  irf3 <- linear_IRF(stvar=mod122, N=1, regime=2)
  irf4 <- linear_IRF(stvar=mod322, N=2, regime=1, ci=0.60, ncores=1, bootstrap_reps=2, which_cumulative=2,
                     scale=c(2, 1, 0.4), seed=1)

})

