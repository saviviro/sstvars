context("parameterReforms")
library(sstvars)

data0 <- t(t(gdpdef))
data2 <- unname(data0)
data3 <- cbind(1:244, data2)
n_obs <- nrow(data3) # 244
dat2_1 <- reform_data(data2, p=1)
dat2_2 <- reform_data(data2, p=2)
dat2_3 <- reform_data(data2, p=3)
dat3_1 <- reform_data(data3, p=1)
dat3_2 <- reform_data(data3, p=2)

test_that("reform_data works correctly", {
  expect_equal(dat2_1, data2)

  expect_equal(nrow(dat2_2), n_obs - 2 + 1)
  expect_equal(dat2_2[1,], c(data2[2,], data2[1,]))
  expect_equal(dat2_2[23,], c(data2[24,], data2[23,]))
  expect_equal(dat2_2[nrow(dat2_2),], c(data2[244,], data2[243,]))

  expect_equal(nrow(dat2_3), n_obs - 3 + 1)
  expect_equal(dat2_3[1,], c(data2[3,], data2[2,], data2[1,]))
  expect_equal(dat2_3[100,], c(data2[102,], data2[101,], data2[100,]))
  expect_equal(dat2_3[nrow(dat2_3),], c(data2[244,], data2[243,], data2[242,]))

  expect_equal(dat3_1, data3)

  expect_equal(nrow(dat3_2), n_obs - 2 + 1)
  expect_equal(dat3_2[1,], c(data3[2,], data3[1,]))
  expect_equal(dat3_2[13,], c(data3[14,], data3[13,]))
  expect_equal(dat3_2[nrow(dat3_2),], c(data3[244,], data3[243,]))
})


## A(M)(p)_(p)(M)(d)

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

# p=3, M=1, d=2
phi10_312 <- phi10_212; A11_312 <- A11_212; A12_312 <- A12_212; Omega1_312 <- Omega1_212
A13_312 <-  matrix(c(-0.15, -0.01, 0.06, -0.23), nrow=2, byrow=FALSE)
theta_312relg <- c(phi10_312, vec(A11_312), vec(A12_312), vec(A13_312), vech(Omega1_312))

# p=1, M=2, d=2
phi10_122 <- c(0.55, 0.11)
A11_122 <- matrix(c(0.34, 0.05, -0.01, 0.72), nrow=2, byrow=FALSE)
Omega1_122 <- matrix(c(0.58, 0.01, 0.01, 0.06), nrow=2, byrow=FALSE)

phi20_122 <- c(0.17, 0.25)
A21_122 <- A11_122
Omega2_122 <- matrix(c(0.50, -0.01, -0.01, 0.20), nrow=2, byrow=FALSE)

alpha1_122 <- 0.60
theta_122relg <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)
alpha1_122_2 <- 0.40
theta_122relg_2 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), alpha1_122_2)

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

alpha1_222_2 <- 0.1
theta_222relg_2 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                   vech(Omega1_222), vech(Omega2_222), alpha1_222_2)


# p=1, M=3, d=2
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

alpha1_132_2 <- 0.3; alpha2_132_2 <- 0.5
theta_132relg_2 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                   vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), alpha1_132_2, alpha2_132_2)

alpha1_132_3 <- 0.1; alpha2_132_3 <- 0.3
theta_132relg_3 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                     vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), alpha1_132_3, alpha2_132_3)

alpha1_132_4 <- 0.6; alpha2_132_4 <- 0.1
theta_132relg_4 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                     vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), alpha1_132_4, alpha2_132_4)


# p=1, M=1, d=3
phi10_113 <- c(1, 2, 3)
A11_113 <- matrix(c(0.1, 0.02, 0.12, 0.3, 0.21, 0.11, 0.05, 0.03, 0.09), nrow=3)
Omega1_113 <- matrix(c(c(1, 0.2, 0.3, 0.2, 2, 0.4, 0.3, 0.4, 3)), nrow=3)

theta_113relg <- c(phi10_113, vec(A11_113), vech(Omega1_113))

# p=2, M=1, d=3
phi10_213 <- phi10_113; A11_213 <- A11_113; Omega1_213 <- Omega1_113
A12_213 <- matrix(c(0.13, 0.03, 0.21, 0.03, 0.14, 0.15, 0.06, 0.07, 0.08), nrow=3)
theta_213relg <- c(phi10_213, vec(A11_213), vec(A12_213), vech(Omega1_213))

# p=1, M=2, d=3
phi10_123 <- phi10_113; A11_123 <- A11_113; A21_123 <- A12_213; Omega1_123 <- Omega1_113
phi20_123 <- c(0.1, 0.2, 0.3)
Omega2_123 <- matrix(c(c(1.1, -0.2, -0.3, -0.2, 2.2, -0.4, -0.3, -0.4, 3.3)), nrow=3)
alpha1_123 <- 0.6
theta_123relg <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vech(Omega1_123),
                   vech(Omega2_123), alpha1_123)
alpha1_123_2 <- 0.2
theta_123relg_2 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vech(Omega1_123),
                     vech(Omega2_123), alpha1_123_2)



allA_112 <- pick_allA(p=1, M=1, d=2, params=theta_112relg)
allA_212 <- pick_allA(p=2, M=1, d=2, params=theta_212relg)
allA_312 <- pick_allA(p=3, M=1, d=2, params=theta_312relg)
allA_122 <- pick_allA(p=1, M=2, d=2, params=theta_122relg)
allA_222 <- pick_allA(p=2, M=2, d=2, params=theta_222relg)
allA_132 <- pick_allA(p=1, M=3, d=2, params=theta_132relg)
allA_113 <- pick_allA(p=1, M=1, d=3, params=theta_113relg)
allA_213 <- pick_allA(p=2, M=1, d=3, params=theta_213relg)
allA_123 <- pick_allA(p=1, M=2, d=3, params=theta_123relg)

lower_part <- function(p, d) {
  cbind(diag(nrow=d*(p-1)), matrix(0, nrow=d*(p - 1), ncol=d))
}

test_that("form_boldA works correctly", {
  expect_equal(form_boldA(p=1, M=1, d=2, all_A=allA_112)[, , 1], A11_112)
  expect_equal(form_boldA(p=2, M=1, d=2, all_A=allA_212)[, , 1], rbind(cbind(A11_212, A12_212), lower_part(p=2, d=2)))
  expect_equal(form_boldA(p=3, M=1, d=2, all_A=allA_312)[, , 1], rbind(cbind(A11_312, A12_312, A13_312), lower_part(p=3, d=2)))

  expect_equal(form_boldA(p=1, M=2, d=2, all_A=allA_122)[, , 1], A11_122)
  expect_equal(form_boldA(p=1, M=2, d=2, all_A=allA_122)[, , 2], A21_122)
  expect_equal(form_boldA(p=2, M=2, d=2, all_A=allA_222)[, , 1], rbind(cbind(A11_222, A12_222), lower_part(p=2, d=2)))
  expect_equal(form_boldA(p=2, M=2, d=2, all_A=allA_222)[, , 2], rbind(cbind(A21_222, A22_222), lower_part(p=2, d=2)))
  expect_equal(form_boldA(p=1, M=3, d=2, all_A=allA_132)[, , 1], A11_132)
  expect_equal(form_boldA(p=1, M=3, d=2, all_A=allA_132)[, , 2], A21_132)
  expect_equal(form_boldA(p=1, M=3, d=2, all_A=allA_132)[, , 3], A31_132)

  expect_equal(form_boldA(p=1, M=1, d=3, all_A=allA_113)[, , 1], A11_113)
  expect_equal(form_boldA(p=2, M=1, d=3, all_A=allA_213)[, , 1], rbind(cbind(A11_213, A12_213), lower_part(p=2, d=3)))
  expect_equal(form_boldA(p=1, M=2, d=3, all_A=allA_123)[, , 1], A11_123)
  expect_equal(form_boldA(p=1, M=2, d=3, all_A=allA_123)[, , 2], A21_123)
})


calc_mu <- function(p, M, d, params, AR_constraints=NULL) {
  stopifnot(is.null(AR_constraints))
  #params <- reform_constrained_pars(p, M, d, params, model=model, constraints=constraints, same_means=NULL,
  #                                  structural_pars=structural_pars)
  all_A <- pick_allA(p=p, M=M, d=d, params=params)
  all_phi0 <- pick_phi0(M=M, d=d, params=params)
  vapply(1:M, function(m) solve(diag(d) - rowSums(all_A[, , , m, drop=FALSE], dims=2), all_phi0[,m]), numeric(d))
}

theta_112relg_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112relg, change_to="mean")
theta_212relg_mu <- change_parametrization(p=2, M=1, d=2, params=theta_212relg, change_to="mean")
theta_312relg_mu <- change_parametrization(p=3, M=1, d=2, params=theta_312relg, change_to="mean")
theta_122relg_mu <- change_parametrization(p=1, M=2, d=2, params=theta_122relg, change_to="mean")
theta_222relg_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222relg, change_to="mean")
theta_132relg_mu <- change_parametrization(p=1, M=3, d=2, params=theta_132relg, change_to="mean")
theta_113relg_mu <- change_parametrization(p=1, M=1, d=3, params=theta_113relg, change_to="mean")
theta_213relg_mu <- change_parametrization(p=2, M=1, d=3, params=theta_213relg, change_to="mean")
theta_123relg_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123relg, change_to="mean")

test_that("change_parametrization works correctly", {
  expect_equal(pick_phi0(M=1, d=2, params=theta_112relg_mu), calc_mu(p=1, M=1, d=2, params=theta_112relg))
  expect_equal(change_parametrization(p=1, M=1, d=2, params=theta_112relg_mu, change_to="intercept"), theta_112relg)
  expect_equal(pick_phi0(M=1, d=2, params=theta_212relg_mu), calc_mu(p=2, M=1, d=2, params=theta_212relg))
  expect_equal(change_parametrization(p=2, M=1, d=2, params=theta_212relg_mu, change_to="intercept"), theta_212relg)
  expect_equal(pick_phi0(M=1, d=2, params=theta_312relg_mu), calc_mu(p=3, M=1, d=2, params=theta_312relg))
  expect_equal(change_parametrization(p=3, M=1, d=2, params=theta_312relg_mu, change_to="intercept"), theta_312relg)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122relg_mu), calc_mu(p=1, M=2, d=2, params=theta_122relg))
  expect_equal(change_parametrization(p=1, M=2, d=2, params=theta_122relg_mu, change_to="intercept"), theta_122relg)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222relg_mu), calc_mu(p=2, M=2, d=2, params=theta_222relg))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222relg_mu, change_to="intercept"), theta_222relg)
  expect_equal(pick_phi0(M=3, d=2, params=theta_132relg_mu), calc_mu(p=1, M=3, d=2, params=theta_132relg))
  expect_equal(change_parametrization(p=1, M=3, d=2, params=theta_132relg_mu, change_to="intercept"), theta_132relg)
  expect_equal(pick_phi0(M=1, d=3, params=theta_113relg_mu), calc_mu(p=1, M=1, d=3, params=theta_113relg))
  expect_equal(change_parametrization(p=1, M=1, d=3, params=theta_113relg_mu, change_to="intercept"), theta_113relg)
  expect_equal(pick_phi0(M=1, d=3, params=theta_213relg_mu), calc_mu(p=2, M=1, d=3, params=theta_213relg))
  expect_equal(change_parametrization(p=2, M=1, d=3, params=theta_213relg_mu, change_to="intercept"), theta_213relg)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123relg_mu), calc_mu(p=1, M=2, d=3, params=theta_123relg))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123relg_mu, change_to="intercept"), theta_123relg)
})


test_that("sort_regimes works correctly", {
  expect_equal(sort_regimes(p=1, M=1, d=2, params=theta_112relg, weight_function="relative_dens"), theta_112relg)
  expect_equal(sort_regimes(p=2, M=1, d=2, params=theta_212relg, weight_function="relative_dens"), theta_212relg)
  expect_equal(sort_regimes(p=3, M=1, d=2, params=theta_312relg, weight_function="relative_dens"), theta_312relg)
  expect_equal(sort_regimes(p=1, M=2, d=2, params=theta_122relg, weight_function="relative_dens"), theta_122relg)
  expect_equal(sort_regimes(p=1, M=2, d=2, params=theta_122relg_2, weight_function="relative_dens"),
               c(phi20_122, phi10_122, vec(A21_122), vec(A11_122), vech(Omega2_122), vech(Omega1_122), 1-alpha1_122_2))
  expect_equal(sort_regimes(p=2, M=2, d=2, params=theta_222relg_2, weight_function="relative_dens"),
               c(phi20_222, phi10_222, vec(A21_222), vec(A22_222), vec(A11_222), vec(A12_222),
                 vech(Omega2_222), vech(Omega1_222), 1-alpha1_222_2))
  expect_equal(sort_regimes(p=1, M=3, d=2, params=theta_132relg, weight_function="relative_dens"), theta_132relg)
  expect_equal(sort_regimes(p=1, M=3, d=2, params=theta_132relg_2, weight_function="relative_dens"),
              c(phi20_132, phi10_132, phi30_132, vec(A21_132), vec(A11_132), vec(A31_132),
                vech(Omega2_132), vech(Omega1_132), vech(Omega3_132), alpha2_132_2, alpha1_132_2))
  expect_equal(sort_regimes(p=1, M=3, d=2, params=theta_132relg_3, weight_function="relative_dens"),
               c(phi30_132, phi20_132, phi10_132, vec(A31_132), vec(A21_132), vec(A11_132),
                 vech(Omega3_132), vech(Omega2_132), vech(Omega1_132), 1-alpha1_132_3-alpha2_132_3, alpha2_132_3))
  expect_equal(sort_regimes(p=1, M=3, d=2, params=theta_132relg_4, weight_function="relative_dens"),
               c(phi10_132, phi30_132, phi20_132, vec(A11_132), vec(A31_132), vec(A21_132),
                 vech(Omega1_132), vech(Omega3_132), vech(Omega2_132), alpha1_132_4, 1-alpha1_132_4-alpha2_132_4))

  expect_equal(sort_regimes(p=1, M=1, d=3, params=theta_113relg, weight_function="relative_dens"), theta_113relg)
  expect_equal(sort_regimes(p=2, M=1, d=3, params=theta_213relg, weight_function="relative_dens"), theta_213relg)
  expect_equal(sort_regimes(p=1, M=2, d=3, params=theta_123relg, weight_function="relative_dens"), theta_123relg)
  expect_equal(sort_regimes(p=1, M=2, d=3, params=theta_123relg_2, weight_function="relative_dens"),
               c(phi20_123, phi10_123, vec(A21_123), vec(A11_123), vech(Omega2_123),
                 vech(Omega1_123), 1-alpha1_123_2))
})


rpars122_1 <-c(phi20_122, vec(A21_122))
rpars122_2 <-c(phi10_122, vec(A11_122))

rpars122t_1 <- c(rpars122_1, 1001)
rpars122t_2 <- c(rpars122_2, 1010)

rpars332_1 <- c(phi10_122, vec(A11_122), vec(A11_122), vec(A11_122))
rpars332_2 <- c(phi20_222, vec(A21_222), vec(A22_222), vec(A22_222))

rpars332t_1 <- c(rpars332_1, 3)
rpars332t_2 <- c(rpars332_2, 13)

test_that("regime_distance works correctly", {
  expect_equal(regime_distance(regime_pars1=rpars122_1, regime_pars2=rpars122_2), 0.4049691, tol=1e-4)
  expect_equal(regime_distance(regime_pars1=rpars122t_1, regime_pars2=rpars122t_2), 0.4049701, tol=1e-4)
  expect_equal(regime_distance(regime_pars1=rpars332_1, regime_pars2=rpars332_2), 0.8523497, tol=1e-4)
  expect_equal(regime_distance(regime_pars1=rpars332t_1, regime_pars2=rpars332t_2), 0.8691375, tol=1e-4)
})
