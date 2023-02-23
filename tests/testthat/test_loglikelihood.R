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


# ADD MORE M>1 TESTS WHEN SOME KIND OF ESTIMATES CAN BE OBTAINED!
# microbenchmark::microbenchmark(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222relg))
# p=2, M=2, d=2: alle 1 millisekuntti menee kaikkeen muuhun + reilu 4 millisekunttia menee multivariate normal arvoihin; vapply ei nopeuta.
# chol2inv on hieman nopeampi kuin backsolve

# mvfast: medion 5.0-5.3ms (usein noin 5.1)
# Oma: noin yhtä suuri.

test_that("loglikelihood works correctly", {
  # Relative_dens Gausssian STVAR
  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112relg, weight_function="relative_dens"), -1000.653, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=1, params=theta_212relg, weight_function="relative_dens"), -286.5474, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=8, M=1, params=theta_812relg, weight_function="relative_dens"), -257.0505, tolerance=1e-3)
  expect_equal(loglikelihood(data=data2, p=4, M=1, params=theta_413relg, weight_function="relative_dens"), -596.6938, tolerance=1e-3)

  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122relg, weight_function="relative_dens"), -314.6693, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222relg, weight_function="relative_dens"), -239.3485, tolerance=1e-3)
})
