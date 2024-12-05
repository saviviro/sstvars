context("moreEst")
library(sstvars)

# The estimation process computationally demanding. Therefore, we only employ minimal
# tests checking that the functions work correctly without producing any errors.

# p=1, M=2, d=2 relative_dens Gaussian STVAR with the means and AR matrices constrained to be identical in both regimes
mod12cm <- STVAR(gdpdef, p=1, M=2, params=c(0.831093, 0.474375, 0.317932, 0.016297, -0.125438, 0.886249, 0.972206,
                                            -0.034546, 0.143727, 0.324064, 0.013734, 0.010002, 0.417756),
                 AR_constraints=rbind(diag(1*2^2), diag(1*2^2)), mean_constraints=list(1:2), parametrization="mean",
                 calc_std_errors=FALSE)
fit12cm <- iterate_more(mod12cm, maxit=30, calc_std_errors=FALSE, print_trace=FALSE) # Not found later if defined inside test_that only

test_that("iterate_more works correctly", {
  expect_equal(fit12cm$params[1], 0.8340364, tolerance=1e-4)
})


test_that("fitSSTVAR works correctly", {
  fit12cm_rec <- fitSSTVAR(fit12cm, identification="recursive", calc_std_errors=FALSE)
  expect_equal(fit12cm_rec$params[1], 0.8340364, tolerance=1e-4)

  fit12cm_hsked <- fitSSTVAR(fit12cm, identification="heteroskedasticity", calc_std_errors=FALSE)
  expect_equal(fit12cm_hsked$params[1], 0.8340364, tolerance=1e-4)

  fit12cm_hsked2 <- fitSSTVAR(fit12cm_hsked, identification="heteroskedasticity", maxit=100, maxit_robust=100, print_res=FALSE,
                              B_constraints=matrix(c(1, 0, NA, 1), nrow=2), robust_method="Nelder-Mead", calc_std_errors=FALSE)
  expect_equal(fit12cm_hsked2$params[1], 0.8379392, tolerance=1e-4)

})

