context("MAINest")
library(sstvars)

# The estimation process computationally very demanding and sometimes numerical errors due to machine (in)accuracy
# yield varying results even with a fixed seed. Therefore, we only employ minimal tests that check that the function
# works without producing any errors.

test_that("fitSTVAR works without errors", {
  # p=1, M=2, d=2 relative_dens Gaussian STVAR with the means and AR matrices constrained to be identical in both regimes
  fit12cm <- fitSTVAR(gdpdef, p=1, M=2, AR_constraints=rbind(diag(1*2^2), diag(1*2^2)), mean_constraints=list(1:2),
                      parametrization="mean", nrounds=1, seeds=1, use_parallel=FALSE, print_res=FALSE, ngen=2, maxit=2)
  expect_equal(class(fit12cm), "stvar")

  # Three-phase procedure
  fit12thres <- fitSTVAR(gdpdef, p=1, M=2, weight_function="threshold", weightfun_pars=c(1, 1),
                         AR_constraints=rbind(diag(1*2^2), diag(1*2^2)), parametrization="intercept",
                         weight_constraints=list(R=0, r=1),
                         nrounds=1, seeds=1, use_parallel=FALSE, print_res=FALSE, ngen=2, maxit=2)
  expect_equal(class(fit12thres), "stvar")
})

