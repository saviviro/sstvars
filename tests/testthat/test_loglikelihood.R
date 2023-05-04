context("loglikelihood")
library(sstvars)

set.seed(1); data2 <- cbind(gdpdef, round(rnorm(nrow(gdpdef)), 3))

# p=1, M=1, d=2
phi10_112 <- c(0.65, 0.7)
A11_112 <- matrix(c(0.29, 0.02, -0.14, 0.9), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.60, 0.01, 0.01, 0.07), nrow=2, byrow=FALSE)

theta_112relg <- c(phi10_112, vec(A11_112), vech(Omega1_112))

# p=2, M=1, d=2
phi10_212 <- c(0.53, 0.03)
A11_212 <- matrix(c(0.23, 0.02, -0.17, 0.66), nrow=2, byrow=FALSE)
A12_212 <- matrix(c(0.18, 0.02, 0.04, 0.26), nrow=2, byrow=FALSE)
Omega1_212 <- matrix(c(0.58, 0.01, 0.01, 0.06), nrow=2, byrow=FALSE)

theta_212relg <- c(phi10_212, vec(A11_212), vec(A12_212), vech(Omega1_212))

# p=1, M=2, d=2
phi10_122 <- c(0.55, 0.11)
A11_122 <- matrix(c(0.34, 0.05, -0.01, 0.72), nrow=2, byrow=FALSE)
Omega1_122 <- matrix(c(0.58, 0.01, 0.01, 0.06), nrow=2, byrow=FALSE)

phi20_122 <- c(0.17, 0.25)
A21_122 <- A11_122
Omega2_122 <- matrix(c(0.50, -0.01, -0.01, 0.20), nrow=2, byrow=FALSE)

alpha1_122 <- 0.60
theta_122relg <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)

# p=2, M=2, d=2
phi10_222 <- c(0.36, 0.12)
A11_222 <- matrix(c(0.22, 0.06, -0.15, 0.39), nrow=2, byrow=FALSE)
A12_222 <- matrix(c(0.41, -0.01, 0.08, 0.3), nrow=2, byrow=FALSE)
Omega1_222 <- matrix(c(0.21, 0.01, 0.01, 0.03), nrow=2, byrow=FALSE)

phi20_222 <- c(0.48, 0.07)
A21_222 <- matrix(c(0.22, 0.02, -0.12, 0.72), nrow=2, byrow=FALSE)
A22_222 <- matrix(c(0.09, 0.03, 0.04, 0.19), nrow=2, byrow=FALSE)
Omega2_222 <- matrix(c(1.10, 0.01, 0.01, 0.11), nrow=2, byrow=FALSE)

alpha1_222 <- 0.37
theta_222relg <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                   vech(Omega1_222), vech(Omega2_222), alpha1_222)

# p=8, M=1, d=2
theta_812relg <- c(0.5427, -0.044, 0.265, 0.0261, -0.1494, 0.5883, 0.1619, 0.0056, 0.0387, 0.1264, -0.0225,
                   0.0216, -0.3024, 0.1455, 0.0752, 0.0346, -0.0899, 0.1902, -0.0944, -0.0124, 0.2447,
                   -0.0966, 0.0684, 0.0095, 0.0399, 0.0423, -0.0234, 0.0403, 0.1939, -0.0795, -0.0161,
                   -0.0034, -0.0901, 0.0252, 0.5273, -0.0014, 0.0574)

# p=4, M=1, d=3, data2
theta_413relg <- c(0.5334, -0.036, 0.0065, 0.2421, 0.0198, -0.1067, -0.1501, 0.573, 0.0079, -0.0118,
                   0.0039, -0.0329, 0.1896, 0.0092, -0.0428, 0.1309, 0.1249, -0.4213, 0.0265, 0.0179,
                   -0.0475, -0.0396, 0.028, -0.0352, -0.2699, 0.1325, 0.3774, -0.0336, 0.0229, 0.015,
                   0.0405, 0.0421, 0.1274, 0.1589, 0.1201, 0.1134, -0.0358, 0.0028, 0.0502, 0.5676,
                   -0.0024, -0.0356, 0.0582, -7e-04, 0.8796)

# p=3, M=1, d=3, usamone
theta_313relg <- c(0.14652, 0.07905, -0.06877, 0.85178, -0.0212, 0.15671, -0.05778, 0.49156, 0.12683,
                   0.04203, 0.06951, 1.19613, 0.08181, 0.02087, -0.03947, 0.1369, 0.31143, 0.52774,
                   -0.29614, -0.05912, -0.48768, -0.13135, -0.0019, 0.04491, -0.1174, 0.10878, -0.27825,
                   0.22554, -0.0108, 0.24196, 0.78022, 0.03763, 0.11717, 0.07114, 0.01935, 0.57962)

# p=1, M=2, d=3, usamone
theta_123relg <- c(0.10741, 0.13813, -0.12092, 3.48957, 0.60615, 0.45646, 0.87227, -0.01595, 0.14124,
                   -0.08611, 0.61865, 0.34311, -0.02047, 0.025, 0.97548, 0.74976, 0.02187, 0.29213,
                   -1.55165, 0.58245, -0.00696, -0.07261, 0.02021, 0.96883, 0.66149, 0.02279, 0.09207,
                   0.05544, 0.00212, 0.12708, 0.78618, 0.00922, 0.42627, 0.23765, 0.25386, 3.40834, 0.77357)

# p=3, M=2, d=3, usamone
theta_323relg <- c(0.98249, 0.66144, -1.17552, 0.50289, 0.17399, -0.01771, 0.96105, -0.11406, 0.41223,
                   -0.31217, 0.49067, 0.3958, 0.04185, 0.08454, 1.0977, -0.03208, 0.06398, -0.12298,
                   0.13382, 0.20166, 0.87613, -0.34591, -0.06254, -0.47386, -0.09049, 0.03109, 0.0347,
                   -0.16531, 0.0427, -0.31646, 0.25299, -0.04865, 0.33893, 0.69963, -0.02912, 0.03398,
                   -0.24344, 0.20815, 0.22566, 0.20582, 0.14774, 1.69008, 0.04375, -0.01018, -0.00947,
                   -0.19371, 0.26341, 0.22082, -0.08841, -0.18303, -0.86488, -0.06031, 0.00634, 0.00181,
                   -0.5559, 0.10249, -0.25146, -0.11875, 0.05153, 0.15267, 0.58151, -0.01903, 0.12236, 0.09327,
                   0.10245, 1.81845, 0.72719, 0.03235, 0.09857, 0.04826, 0.00908, 0.09761, 0.72127)

################### Speed test ##
#microbenchmark::microbenchmark(loglikelihood(data=usamone, p=3, M=2, d=3, params=theta_323relg), times=1000L)
###################


# p=1, M=3, d=2 (not tested atm)
phi10_132 <- phi10_122
phi20_132 <- phi20_122
phi30_132 <- c(12, 13)

A11_132 <- A11_122
A21_132 <- A21_122
A31_132 <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow=2)
Omega1_132 <- Omega1_122
Omega2_132 <- Omega2_122
Omega3_132 <- matrix(c(1, 0.5, 0.5, 1), nrow=2)
alpha1_132 <- 0.5
alpha2_132 <- 0.3

theta_132relg <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                   vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), alpha1_132, alpha2_132)

## weight_function = "logit"

# p=1, M=2, d=2, weightfun_pars=list(vars=1, lags=1)
gamma1_122_1_1 <- c(0.1, 0.2)
theta_122log_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), gamma1_122_1_1)

# p=1, M=2, d=2, weightfun_pars=list(vars=1:2, lags=1)
gamma1_122_12_1 <- c(0.1, 0.2, 0.3)
theta_122log_12_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), gamma1_122_12_1)

# p=2, M=2, d=2, weightfun_pars=list(vars=2, lags=1)
gamma1_222_2_1 <- c(0.1, 0.2)
theta_222log_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                      vech(Omega1_222), vech(Omega2_222), gamma1_222_2_1)

# p=2, M=2, d=2, weightfun_pars=list(vars=1:2, lags=2)
gamma1_222_12_2 <- c(0.1, 0.2, 0.11, 0.22, 0.33)
theta_222log_12_2 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                       vech(Omega1_222), vech(Omega2_222),gamma1_222_12_2)

# p=2, M=3, d=2, weightfun_pars=list(vars=1, lags=1)
phi10_232 <- phi10_132; phi20_232 <- phi20_132; phi30_232 <- phi30_132
A11_232 <- A11_132; A12_232 <- -A11_132
A21_232 <- A21_132; A22_232 <- -A21_132
A31_232 <- A31_132; A32_232 <- -A31_132
Omega1_232 <- Omega1_132; Omega2_232 <- Omega1_132; Omega3_232 <- Omega3_132
gamma1_232_1_1 <- c(0.1, 0.2); gamma2_232_1_1 <- c(0.11, 0.22)
theta_232log_1_1 <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                      vec(A31_232), vec(A32_232), vech(Omega1_232), vech(Omega2_232), vech(Omega3_232),
                      gamma1_232_1_1, gamma2_232_1_1)

# p=2, M=3, d=2, weightfun_pars=list(vars=1:2, lags=1)
gamma1_232_12_1 <- c(0.1, 0.2, 0.202); gamma2_232_12_1 <- c(0.11, 0.111, 0.222)
theta_232log_12_1 <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                       vec(A31_232), vec(A32_232), vech(Omega1_232), vech(Omega2_232), vech(Omega3_232),
                       gamma1_232_12_1, gamma2_232_12_1)

# p=2, M=3, d=2, weightfun_pars=list(vars=2, lags=2)
gamma1_232_2_2 <- c(0.1, 0.2, 0.3); gamma2_232_2_2 <- c(0.11, 0.22, 0.33)
theta_232log_2_2 <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                      vec(A31_232), vec(A32_232), vech(Omega1_232), vech(Omega2_232), vech(Omega3_232),
                      gamma1_232_2_2, gamma2_232_2_2)

# p=2, M=3, d=2, weightfun_pars=list(vars=1:2, lags=2)
gamma1_232_12_2 <- c(0.1, 0.2, 0.101, 0.202, 0.303); gamma2_232_12_2 <- c(0.11, 0.22, 0.111, 0.222, 0.333)
theta_232log_12_2 <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                       vec(A31_232), vec(A32_232), vech(Omega1_232), vech(Omega2_232), vech(Omega3_232),
                       gamma1_232_12_2, gamma2_232_12_2)


# p=1, M=2, d=3, weightfun_pars=list(vars=1, lags=1), usamone
theta_123log_noweightpars <- theta_123relg[-length(theta_123relg)]
gamma1_123_1_1 <- c(0.1, 0.2)
theta_123log_1_1 <- c(theta_123log_noweightpars, gamma1_123_1_1)

# p=1, M=2, d=3, weightfun_pars=list(vars=2:3, lags=1), usamone
gamma1_123_23_1 <- c(0.1, 0.2, 0.3)
theta_123log_23_1 <- c(theta_123log_noweightpars, gamma1_123_23_1)

# p=1, M=2, d=3, weightfun_pars=list(vars=1:3, lags=1), usamone
gamma1_123_123_1 <- c(0.1, 0.2, 0.3, 0.4)
theta_123log_123_1 <- c(theta_123log_noweightpars, gamma1_123_123_1)


## Constrained models

## A(M)(p)_(p)(M)(d)
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}

# p=1, M=1, d=2 AR_constraints
C_112 <- matrix(c(1, 0, -0.5, 0, 0, 0.35, 0, -1), nrow=4, ncol=2)
psi_112 <- c(0.3, 0.15)
theta_112relgc <- c(phi10_112, psi_112, vech(Omega1_112))
theta_112relgc_expanded <- c(phi10_112, C_112%*%psi_112, vech(Omega1_112))

# p=1, M=2, d=2 AR matrices identical in both regimes
C_122 <- rbind_diags(p=1, M=2, d=2)
theta_122relgc <- c(phi10_122, phi20_122, vec(A11_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)
theta_122relgc_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A11_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)

# p=2, M=2, d=2 AR matrices identical in both regimes
C_222 <- rbind_diags(p=2, M=2, d=2)
theta_222relgc <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222)
theta_222relgc_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                             vech(Omega1_222), vech(Omega2_222), alpha1_222)

# p=2, M=2, d=2, constrain AR-parameters to be the same for all regimes
# and constrain the of-diagonal elements of AR-matrices to be zero.
mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE)
C_222_2 <- rbind(mat0, mat0)
A21_222_c2 <- A11_222_c2 <- matrix(c(1.26, 0, 0, 1.34), nrow=2, byrow=FALSE)
A22_222_c2 <- A12_222_c2 <- matrix(c(-0.29, 0, 0, -0.36), nrow=2, byrow=FALSE)
phi10_222_c2 <- c(-0.11, 2.83)
phi20_222_c2 <- c(0.36, 3.19)
Omega1_222_c2 <- matrix(c(0.98, -0.33, -0.33, 5.24), nrow=2, byrow=FALSE)
Omega2_222_c2 <- matrix(c(5.60, 3.46, 3.46, 9.62), nrow=2, byrow=FALSE)
alpha1_222_c2 <- 0.35
theta_222relgc2 <- c(phi10_222_c2, phi20_222_c2, 1.26, 1.34, -0.29, -0.36, vech(Omega1_222_c2),
                     vech(Omega2_222_c2), alpha1_222_c2)
theta_222relgc2_expanded <- c(phi10_222_c2, phi20_222_c2, vec(A11_222_c2), vec(A12_222_c2), vec(A21_222_c2), vec(A22_222_c2),
                              vech(Omega1_222_c2), vech(Omega2_222_c2), alpha1_222_c2)

# p=1, M=2, p=3, AR_constraints
C_123 <- rbind_diags(p=1, M=2, d=3)
phi10_123 <- c(1, 2, 3); phi20_123 <- c(0.1, 0.2, 0.3)
A11_123 <- matrix(c(0.1, 0.02, 0.12, 0.3, 0.21, 0.11, 0.05, 0.03, 0.09), nrow=3)
Omega1_123 <- matrix(c(c(1, 0.2, 0.3, 0.2, 2, 0.4, 0.3, 0.4, 3)), nrow=3)
Omega2_123 <- matrix(c(c(1.1, -0.2, -0.3, -0.2, 2.2, -0.4, -0.3, -0.4, 3.3)), nrow=3)
alpha1_123 <- 0.6
theta_123relgc <- c(phi10_123, phi20_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), alpha1_123)
theta_123relgc_expanded <- c(phi10_123, phi20_123, vec(A11_123), vec(A11_123), vech(Omega1_123), vech(Omega2_123), alpha1_123)


# weight_function == "logit"

# p=1, M=2, d=2, weightfun_pars=list(vars=1:2, lags=1), C_122
gamma1_122_12_1 <- c(0.1, 0.2, 0.3)
theta_122logc_12_1 <- c(phi10_122, phi20_122, vec(A11_122), vech(Omega1_122), vech(Omega2_122), gamma1_122_12_1)
theta_122logc_12_1_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A11_122), vech(Omega1_122), vech(Omega2_122), gamma1_122_12_1)

# p=2, M=2, d=2, weightfun_pars=list(vars=2, lags=1), C_222
gamma1_222_2_1 <- c(0.1, 0.2)
theta_222logc_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), gamma1_222_2_1)
theta_222logc_2_1_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                vech(Omega1_222), vech(Omega2_222), gamma1_222_2_1)


## Models with mean_constraints

# p=1, M=1, p=2, mean_constraints=list(1)
theta_112relgm <- theta_112relg

# p=1, M=2, d=2, mean_constraints=list(1:2)
mu_122relgm <- c(0.5342688, 0.7382623)
theta_122relgm <- c(mu_122relgm, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)
theta_122relgm_expanded <- c(mu_122relgm, mu_122relgm, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)

# p=2, M=2, d=2, mean_constraints=list(1:2), C_222
mu_222relgcm <- c(0.7209658, 0.8108580)
theta_222relgcm <- c(mu_222relgcm, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222)
theta_222relgcm_expanded <- c(mu_222relgcm, mu_222relgcm, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                              vech(Omega1_222), vech(Omega2_222), alpha1_222)

# weightfunction == "logit"

# p=1, M=2, d=2, weigthfun_pars=list(vars=1, lags=1), mean_constraints=list(1:2)
theta_122logm_1_1 <- c(phi10_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), gamma1_122_1_1)
theta_122logm_1_1_expanded <- c(phi10_122, phi10_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), gamma1_122_1_1)

# p=2, M=2, d=2, weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), C_222
theta_222logcm_12_2 <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), gamma1_222_12_2)
theta_222logcm_12_2_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                  vech(Omega1_222), vech(Omega2_222), gamma1_222_12_2)


test_that("loglikelihood works correctly", {
  # Relative_dens Gaussian STVAR
  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112relg, weight_function="relative_dens"), -1000.653, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=1, params=theta_212relg, weight_function="relative_dens"), -286.5474, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=8, M=1, params=theta_812relg, weight_function="relative_dens"), -257.0505, tolerance=1e-3)
  expect_equal(loglikelihood(data=data2, p=4, M=1, params=theta_413relg, weight_function="relative_dens"), -596.6938, tolerance=1e-3)

  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122relg, weight_function="relative_dens"), -314.6693, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222relg, weight_function="relative_dens"), -239.3485, tolerance=1e-3)

  expect_equal(loglikelihood(data=usamone, p=3, M=1, params=theta_313relg, weight_function="relative_dens"), -669.5716, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123relg, weight_function="relative_dens"), -570.019, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=3, M=2, params=theta_323relg, weight_function="relative_dens"), -490.9401, tolerance=1e-3)

  # logit STVAR
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122log_1_1, weight_function="logit",
                             weightfun_pars=list(vars=1, lags=1)), -334.0843, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122log_12_1, weight_function="logit",
                             weightfun_pars=list(vars=1:2, lags=1)), -335.5687, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222log_2_1, weight_function="logit",
                             weightfun_pars=list(vars=2, lags=1)), -315.326, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222log_12_2, weight_function="logit",
                             weightfun_pars=list(vars=1:2, lags=2)), -344.0656, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=3, params=theta_232log_1_1, weight_function="logit",
                             weightfun_pars=list(vars=1, lags=1)), -4260.053, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=3, params=theta_232log_12_1, weight_function="logit",
                             weightfun_pars=list(vars=1:2, lags=1)), -3756.747, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=3, params=theta_232log_2_2, weight_function="logit",
                             weightfun_pars=list(vars=2, lags=2)), -3450.811, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=3, params=theta_232log_12_2, weight_function="logit",
                             weightfun_pars=list(vars=1:2, lags=2)), -2695.943, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123log_1_1, weight_function="logit",
                             weightfun_pars=list(vars=1, lags=1)), -998.6099, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123log_23_1, weight_function="logit",
                             weightfun_pars=list(vars=2:3, lags=1)), -1100.063, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123log_123_1, weight_function="logit",
                             weightfun_pars=list(vars=1:3, lags=1)), -1155.278, tolerance=1e-3)


  # Constrained models
  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112relgc, weight_function="relative_dens",
                             AR_constraints=C_112), -1008.678, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122relgc, weight_function="relative_dens",
                             AR_constraints=C_122), -314.6693, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222relgc, weight_function="relative_dens",
                             AR_constraints=C_222), -299.1835, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222relgc2, weight_function="relative_dens",
                             AR_constraints=C_222_2), -1077.319, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123relgc, weight_function="relative_dens",
                             AR_constraints=C_123), -2171.694, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122logc_12_1, weight_function="logit",
                             weightfun_pars=list(vars=1:2, lags=1), AR_constraints=C_122), -335.5687, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logc_2_1, weight_function="logit",
                             weightfun_pars=list(vars=2, lags=1), AR_constraints=C_222), -369.6914, tolerance=1e-3)

  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112relgm, weight_function="relative_dens", parametrization="mean",
                             mean_constraints=list(1)), -302.0188, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122relgm, weight_function="relative_dens", parametrization="mean",
                             mean_constraints=list(1:2)), -326.8462, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222relgcm, weight_function="relative_dens", parametrization="mean",
                             AR_constraints=C_222, mean_constraints=list(1:2)), -301.0143, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122logm_1_1, weight_function="logit", parametrization="mean",
                             weightfun_pars=list(vars=1, lags=1), mean_constraints=list(1:2)), -368.2033, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logcm_12_2, weight_function="logit", parametrization="mean",
                             weightfun_pars=list(vars=1:2, lags=2), AR_constraints=C_222, mean_constraints=list(1:2)),
               -453.5176, tolerance=1e-3)
})
