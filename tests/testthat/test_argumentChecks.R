context("argumentChecks")
library(sstvars)

## A(M)(p)_(p)(M)(d)


## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2
phi10_112 <- c(0.65, 0.7)
A11_112 <- matrix(c(0.29, 0.02, -0.14, 0.9), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.60, 0.01, 0.01, 0.07), nrow=2, byrow=FALSE)
Omega1_112_notpd <-  matrix(c(0.60, -0.01, -0.01, 0.0000001), nrow=2, byrow=FALSE)

theta_112relg <- c(phi10_112, vec(A11_112), vech(Omega1_112))
theta_112relg_notpd <- c(phi10_112, vec(A11_112), vech(Omega1_112_notpd))

# p=2, M=1, d=2
phi10_212 <- c(0.53, 0.03)
A11_212 <- matrix(c(0.23, 0.02, -0.17, 0.66), nrow=2, byrow=FALSE)
A12_212 <- matrix(c(1.18, 0.02, 1.04, 0.26), nrow=2, byrow=FALSE)
Omega1_212 <- matrix(c(0.58, 0.01, 0.01, 0.06), nrow=2, byrow=FALSE)
A12_212_stab <- matrix(c(0.18, 0.02, 0.04, 0.26), nrow=2, byrow=FALSE)

theta_212relg_notstab <- c(phi10_212, vec(A11_212), vec(A12_212), vech(Omega1_212))
theta_212relg <- c(phi10_212, vec(A11_212), vec(A12_212_stab), vech(Omega1_212))


# p=3, M=1, d=2
phi10_312 <- phi10_212; A11_312 <- A11_212; A12_312 <- A12_212; Omega1_312 <- Omega1_212; A12_312_stab <- A12_212_stab
A13_312 <- matrix(c(-0.15, -0.01, 0.06, -0.23), nrow=2, byrow=FALSE)
A13_312_stab <- matrix(c(-0.015, -0.01, 0.06, -0.023), nrow=2, byrow=FALSE)

theta_312relg_notstab <- c(phi10_312, vec(A11_312), vec(A12_312), vec(A13_312), vech(Omega1_312))
theta_312relg <- c(phi10_312, vec(A11_312), vec(A12_312_stab), vec(A13_312_stab), vech(Omega1_312))

# p=1, M=2, d=2
phi10_122 <- c(0.55, 0.11)
A11_122 <- matrix(c(0.34, 0.05, -0.01, 0.72), nrow=2, byrow=FALSE)
Omega1_122 <- matrix(c(0.58, 0.01, 0.01, 0.06), nrow=2, byrow=FALSE)

phi20_122 <- c(0.17, 0.25)
A21_122 <- A11_122
Omega2_122 <- matrix(c(0.50, -0.01, -0.01, 0.20), nrow=2, byrow=FALSE)
Omega2_122_notpd <- matrix(c(0.010000001, 0.01, 0.01, 0.01), nrow=2, byrow=FALSE)

alpha1_122 <- 0.60
theta_122relg <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)
theta_122relg_notpd <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122_notpd), alpha1_122)

theta_122relg_badalphas <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122_notpd), 1)

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

# p=1, M=3, d=2
phi10_132 <- phi10_122
phi20_132 <- phi20_122
phi30_132 <- c(12, 13)

A11_132 <- A11_122
A21_132 <- A21_122
A31_132 <- matrix(c(0.1, 2.2, 0.3, 0.4), nrow=2)
A31_132_stab <-  matrix(c(0.1, 0.2, 0.3, 0.1), nrow=2)
Omega1_132 <- Omega1_122
Omega2_132 <- Omega2_122
Omega3_132 <- matrix(c(1, 0.5, 0.5, 1), nrow=2)
Omega3_132_notpd <- matrix(c(0.005, 0.5, 0.5, 0.5), nrow=2)
alpha1_132 <- 0.5
alpha2_132 <- 0.3

theta_132relg_notstab <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                           vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), alpha1_132, alpha2_132)

theta_132relg_notpd <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132_stab),
                           vech(Omega1_132), vech(Omega2_132), vech(Omega3_132_notpd), alpha1_132, alpha2_132)
theta_132relg_badalphas <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132_stab),
                             vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), 0.6, 0.4)
theta_132relg_badalphas2 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132_stab),
                              vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), 0.6, -0.2)

# p=1, M=1, d=3
phi10_113 <- c(1, 2, 3)
A11_113 <- matrix(c(0.1, 0.02, 0.12, 0.3, 0.21, 0.11, 0.05, 0.03, 0.09), nrow=3)
Omega1_113 <- matrix(c(c(1, 0.2, 0.3, 0.2, 2, 0.4, 0.3, 0.4, 3)), nrow=3)
Omega1_113_notpd <- matrix(c(c(1, 0.2, 0.3, 0.02, 0.2, 0.4, 0.3, 0.4, 0.3)), nrow=3)

theta_113relg <- c(phi10_113, vec(A11_113), vech(Omega1_113))
theta_113relg_notpd <- c(phi10_113, vec(A11_113), vech(Omega1_113_notpd))

# p=2, M=1, d=3
phi10_213 <- phi10_113; A11_213 <- A11_113; Omega1_213 <- Omega1_113
A12_213 <- matrix(c(0.13, 0.03, 0.21, 0.03, 0.14, 2.15, 0.06, -1.07, 0.08), nrow=3)
A12_213_stab <- matrix(c(0.13, 0.03, 0.21, 0.03, 0.14, 0.15, 0.06, -0.07, 0.08), nrow=3)

theta_213relg_notstab <- c(phi10_213, vec(A11_213), vec(A12_213), vech(Omega1_213))
theta_213relg <- c(phi10_213, vec(A11_213), vec(A12_213_stab), vech(Omega1_213))


# p=1, M=2, d=3
phi10_123 <- phi10_113; A11_123 <- A11_113; A21_123 <- A12_213; Omega1_123 <- Omega1_113; A21_123_stab <- A12_213_stab
phi20_123 <- c(0.1, 0.2, 0.3)
Omega2_123 <- matrix(c(c(1.1, -0.2, -0.3, -0.2, 2.2, -0.4, -0.3, -0.4, 3.3)), nrow=3)
Omega2_123_notpd <- matrix(c(c(0.1, -0.2, -0.3, -0.2, 2.2, -0.4, -0.3, -0.4, 0.3)), nrow=3)
alpha1_123 <- 0.6
theta_123relg_notstab <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vech(Omega1_123),
                           vech(Omega2_123), alpha1_123)
theta_123relg_notpd <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123_stab), vech(Omega1_123),
                           vech(Omega2_123_notpd), alpha1_123)

Omegas_112 <- pick_Omegas(p=1, M=1, d=2, params=theta_112relg)
Omegas_112_notpd <- pick_Omegas(p=1, M=1, d=2, params=theta_112relg_notpd)
Omegas_212 <- pick_Omegas(p=2, M=1, d=2, params=theta_212relg_notstab)
Omegas_312 <- pick_Omegas(p=3, M=1, d=2, params=theta_312relg_notstab)
Omegas_122 <- pick_Omegas(p=1, M=2, d=2, params=theta_122relg)
Omegas_122_notpd <- pick_Omegas(p=1, M=2, d=2, params=theta_122relg_notpd)
Omegas_222 <- pick_Omegas(p=2, M=2, d=2, params=theta_222relg)
Omegas_132 <- pick_Omegas(p=1, M=3, d=2, params=theta_132relg_notstab)
Omegas_132_notpd <- pick_Omegas(p=1, M=3, d=2, params=theta_132relg_notpd)
Omegas_113 <- pick_Omegas(p=1, M=1, d=3, params=theta_113relg)
Omegas_113_notpd <- pick_Omegas(p=1, M=1, d=3, params=theta_113relg_notpd)
Omegas_213 <- pick_Omegas(p=2, M=1, d=3, params=theta_213relg_notstab)
Omegas_123 <- pick_Omegas(p=1, M=2, d=3, params=theta_123relg_notstab)
Omegas_123_notpd <- pick_Omegas(p=1, M=2, d=3, params=theta_123relg_notpd)

boldA_112 <- form_boldA(p=1, M=1, d=2, all_A=pick_allA(p=1, M=1, d=2, params=theta_112relg))
boldA_212_notstab <- form_boldA(p=2, M=1, d=2, all_A=pick_allA(p=2, M=1, d=2, params=theta_212relg_notstab))
boldA_212 <- form_boldA(p=2, M=1, d=2, all_A=pick_allA(p=2, M=1, d=2, params=theta_212relg))
boldA_312_notstab <- form_boldA(p=3, M=1, d=2, all_A=pick_allA(p=3, M=1, d=2, params=theta_312relg_notstab))
boldA_312 <- form_boldA(p=3, M=1, d=2, all_A=pick_allA(p=3, M=1, d=2, params=theta_312relg))
boldA_122 <- form_boldA(p=1, M=2, d=2, all_A=pick_allA(p=1, M=2, d=2, params=theta_122relg))
boldA_222 <- form_boldA(p=2, M=2, d=2, all_A=pick_allA(p=2, M=2, d=2, params=theta_222relg))
boldA_132_notstab <- form_boldA(p=1, M=3, d=2, all_A=pick_allA(p=1, M=3, d=2, params=theta_132relg_notstab))
boldA_132 <- form_boldA(p=1, M=3, d=2, all_A=pick_allA(p=1, M=3, d=2, params=theta_132relg_notpd))
boldA_113 <- form_boldA(p=1, M=1, d=3, all_A=pick_allA(p=1, M=1, d=3, params=theta_113relg))
boldA_213_notstab <- form_boldA(p=2, M=1, d=3, all_A=pick_allA(p=2, M=1, d=3, params=theta_213relg_notstab))
boldA_213 <- form_boldA(p=2, M=1, d=3, all_A=pick_allA(p=2, M=1, d=3, params=theta_213relg))
boldA_123_notstab <- form_boldA(p=1, M=2, d=3, all_A=pick_allA(p=1, M=2, d=3, params=theta_123relg_notstab))
boldA_123 <- form_boldA(p=1, M=2, d=3, all_A=pick_allA(p=1, M=2, d=3, params=theta_123relg_notpd))

weightpars_112rel <- pick_weightpars(p=1, M=1, d=2, params=theta_112relg, weight_function="relative_dens", cond_dist="Gaussian")
weightpars_212rel <- pick_weightpars(p=2, M=1, d=2, params=theta_212relg, weight_function="relative_dens", cond_dist="Gaussian")
weightpars_312rel <- pick_weightpars(p=3, M=1, d=2, params=theta_312relg, weight_function="relative_dens", cond_dist="Gaussian")
weightpars_122rel <- pick_weightpars(p=1, M=2, d=2, params=theta_122relg, weight_function="relative_dens", cond_dist="Gaussian")
weightpars_222rel <- pick_weightpars(p=2, M=2, d=2, params=theta_222relg, weight_function="relative_dens", cond_dist="Gaussian")
weightpars_132rel <- pick_weightpars(p=1, M=3, d=2, params=theta_132relg_notpd, weight_function="relative_dens", cond_dist="Gaussian")
weightpars_113rel <- pick_weightpars(p=1, M=1, d=3, params=theta_113relg, weight_function="relative_dens", cond_dist="Gaussian")
weightpars_213rel <- pick_weightpars(p=2, M=1, d=3, params=theta_213relg, weight_function="relative_dens", cond_dist="Gaussian")
weightpars_123rel <- pick_weightpars(p=1, M=2, d=3, params=theta_123relg_notpd, weight_function="relative_dens", cond_dist="Gaussian")


test_that("stab_conds_satisfied works correctly", {
  expect_true(stab_conds_satisfied(p=1, M=1, d=2, all_boldA=boldA_112))
  expect_false(stab_conds_satisfied(p=2, M=1, d=2, all_boldA=boldA_212_notstab))
  expect_false(stab_conds_satisfied(p=3, M=1, d=2, all_boldA=boldA_312_notstab))
  expect_true(stab_conds_satisfied(p=1, M=2, d=2, all_boldA=boldA_122))
  expect_true(stab_conds_satisfied(p=2, M=2, d=2, all_boldA=boldA_222))
  expect_false(stab_conds_satisfied(p=1, M=3, d=2, all_boldA=boldA_132_notstab))
  expect_true(stab_conds_satisfied(p=1, M=1, d=3, all_boldA=boldA_113))
  expect_false(stab_conds_satisfied(p=2, M=1, d=3, all_boldA=boldA_213_notstab))
  expect_false(stab_conds_satisfied(p=1, M=2, d=3, all_boldA=boldA_123_notstab))

  # Without boldA
  expect_true(stab_conds_satisfied(p=1, M=1, d=2, params=theta_112relg))
  expect_false(stab_conds_satisfied(p=3, M=1, d=2, params=theta_312relg_notstab))
  expect_true(stab_conds_satisfied(p=1, M=2, d=2, params=theta_122relg))
  expect_true(stab_conds_satisfied(p=2, M=2, d=2, params=theta_222relg))
  expect_false(stab_conds_satisfied(p=1, M=2, d=3, params=theta_123relg_notstab))
})


test_that("in_paramspace work correctly", {
  # Checks stability conditions
  expect_true(in_paramspace(p=1, M=1, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                            all_boldA=boldA_112, all_Omegas=Omegas_112, weightpars=weightpars_112rel))
  expect_false(in_paramspace(p=2, M=1, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                             all_boldA=boldA_212_notstab, all_Omegas=Omegas_212, weightpars=weightpars_312rel))
  expect_false(in_paramspace(p=3, M=1, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                             all_boldA=boldA_312_notstab, all_Omegas=Omegas_312, weightpars=weightpars_312rel))
  expect_true(in_paramspace(p=1, M=2, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                            all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=weightpars_122rel))
  expect_true(in_paramspace(p=2, M=2, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                            all_boldA=boldA_222, all_Omegas=Omegas_222, weightpars=weightpars_222rel))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                             all_boldA=boldA_132_notstab, all_Omegas=Omegas_132, weightpars=weightpars_132rel))
  expect_true(in_paramspace(p=1, M=1, d=3, weight_function="relative_dens", cond_dist="Gaussian",
                            all_boldA=boldA_113, all_Omegas=Omegas_113, weightpars=weightpars_113rel))
  expect_false(in_paramspace(p=2, M=1, d=3, weight_function="relative_dens", cond_dist="Gaussian",
                             all_boldA=boldA_213_notstab, all_Omegas=Omegas_213, weightpars=weightpars_213rel))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="relative_dens", cond_dist="Gaussian",
                             all_boldA=boldA_123_notstab, all_Omegas=Omegas_123, weightpars=weightpars_123rel))

  # Check Omegas
  expect_true(in_paramspace(p=1, M=1, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                            all_boldA=boldA_112, all_Omegas=Omegas_112, weightpars=weightpars_112rel))
  expect_false(in_paramspace(p=1, M=1, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                            all_boldA=boldA_112, all_Omegas=Omegas_112_notpd, weightpars=weightpars_112rel))
  expect_true(in_paramspace(p=2, M=1, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                            all_boldA=boldA_212, all_Omegas=Omegas_212, weightpars=weightpars_212rel))
  expect_true(in_paramspace(p=1, M=2, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                            all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=weightpars_122rel))
  expect_false(in_paramspace(p=1, M=2, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                             all_boldA=boldA_122, all_Omegas=Omegas_122_notpd, weightpars=weightpars_122rel))
  expect_true(in_paramspace(p=2, M=2, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                            all_boldA=boldA_222, all_Omegas=Omegas_222, weightpars=weightpars_222rel))
  expect_true(in_paramspace(p=1, M=3, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                             all_boldA=boldA_132, all_Omegas=Omegas_132, weightpars=weightpars_132rel))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                             all_boldA=boldA_132, all_Omegas=Omegas_132_notpd, weightpars=weightpars_132rel))
  expect_false(in_paramspace(p=1, M=1, d=3, weight_function="relative_dens", cond_dist="Gaussian",
                             all_boldA=boldA_113, all_Omegas=Omegas_113_notpd, weightpars=weightpars_113rel))
  expect_true(in_paramspace(p=1, M=2, d=3, weight_function="relative_dens", cond_dist="Gaussian",
                             all_boldA=boldA_123, all_Omegas=Omegas_123, weightpars=weightpars_123rel))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="relative_dens", cond_dist="Gaussian",
                             all_boldA=boldA_123, all_Omegas=Omegas_123_notpd, weightpars=weightpars_123rel))

  # Check weightpars
  expect_false(in_paramspace(p=1, M=1, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                            all_boldA=boldA_112, all_Omegas=Omegas_112, weightpars=-0.001))
  expect_false(in_paramspace(p=1, M=2, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                            all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=c(1.0, 0.0001)))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                            all_boldA=boldA_132, all_Omegas=Omegas_132, weightpars=c(0.99, 0.01001, 0.00001)))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                            all_boldA=boldA_132, all_Omegas=Omegas_132, weightpars=c(0.99, -0.01, 0.01)))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="relative_dens", cond_dist="Gaussian",
                            all_boldA=boldA_123, all_Omegas=Omegas_123, weightpars=c(1.01, 0.01)))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="relative_dens", cond_dist="Gaussian",
                             all_boldA=boldA_123, all_Omegas=Omegas_123, weightpars=c(-0.01, 0.99)))

  # Checks df
  # TO BE FILLED IN
})

test_that("check_params work correctly", {
  # Checks stability conditions
  check_params(p=1, M=1, d=2, params=theta_112relg, weight_function="relative_dens", cond_dist="Gaussian")
  check_params(p=1, M=2, d=2, params=theta_122relg, weight_function="relative_dens", cond_dist="Gaussian")
  check_params(p=2, M=2, d=2, params=theta_222relg, weight_function="relative_dens", cond_dist="Gaussian")
  check_params(p=1, M=1, d=3, params=theta_113relg, weight_function="relative_dens", cond_dist="Gaussian")
  expect_error(check_params(p=2, M=1, d=2, params=theta_212relg_notstab, weight_function="relative_dens", cond_dist="Gaussian"))
  expect_error(check_params(p=3, M=1, d=2, params=theta_312relg_notstab, weight_function="relative_dens", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=3, d=2, params=theta_132relg_notstab, weight_function="relative_dens", cond_dist="Gaussian"))
  expect_error(check_params(p=2, M=1, d=3, params=theta_213relg_notstab, weight_function="relative_dens", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=2, d=3, params=theta_123relg_notstab, weight_function="relative_dens", cond_dist="Gaussian"))

  # Check Omegas
  check_params(p=1, M=1, d=2, params=theta_112relg, weight_function="relative_dens", cond_dist="Gaussian")
  check_params(p=2, M=1, d=2, params=theta_212relg, weight_function="relative_dens", cond_dist="Gaussian")
  check_params(p=1, M=2, d=2, params=theta_122relg, weight_function="relative_dens", cond_dist="Gaussian")
  check_params(p=2, M=2, d=2, params=theta_222relg, weight_function="relative_dens", cond_dist="Gaussian")

  expect_error(check_params(p=1, M=1, d=2, params=theta_112relg_notpd, weight_function="relative_dens", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=2, d=2, params=theta_122relg_notpd, weight_function="relative_dens", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=3, d=2, params=theta_132relg_notpd, weight_function="relative_dens", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=1, d=3, params=theta_113relg_notpd, weight_function="relative_dens", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=2, d=3, params=theta_123relg_notpd, weight_function="relative_dens", cond_dist="Gaussian"))

  # Check weightpars
  expect_error(check_params(p=1, M=2, d=2, params=theta_122relg_badalphas, weight_function="relative_dens", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=3, d=2, params=theta_132relg_badalphas, weight_function="relative_dens", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=3, d=2, params=theta_132relg_badalphas2, weight_function="relative_dens", cond_dist="Gaussian"))

  # Checks df
  # TO BE FILLED IN

  # Check weight pars with other weight functions
  # TO BE FILLED IN

  # Check various constraints
  # TO BE FILLED
})


test_that("check_pMd works correctly", {
  expect_error(check_pMd(p=1, M=1, d=1))
  expect_error(check_pMd(p=1, M=1.2, d=2))
  expect_error(check_pMd(p=0, M=1, d=2))
  expect_error(check_pMd(p=1, M=-1, d=2))
  expect_error(check_pMd(p=1.1, M=1, d=2))
  expect_error(check_pMd(p=1, M=1, d=2.2))
  expect_error(check_pMd(p=2, M=c(1, 1), d=2))
  expect_error(check_pMd(p=-1, M=1, d=2))
  expect_error(check_pMd(p=c(1, 1), M=1, d=2))
})

test_that("all_pos_ints works correctly", {
  expect_true(all_pos_ints(c(1, 2, 3)))
  expect_true(all_pos_ints(1))
  expect_true(all_pos_ints(list(1, 3, 100)))
  expect_false(all_pos_ints(c(1, 2, 0)))
  expect_false(all_pos_ints(-1))
  expect_false(all_pos_ints(0.1))
  expect_false(all_pos_ints(1.1))
  expect_false(all_pos_ints(list(1, 2, 3, 0.1)))
})

test_that("n_params works correctly", {
  expect_equal(n_params(p=1, M=1, d=2, weight_function="relative_dens", cond_dist="Gaussian"), 9)
  expect_equal(n_params(p=2, M=1, d=2, weight_function="relative_dens", cond_dist="Gaussian"), 13)
  expect_equal(n_params(p=3, M=1, d=2, weight_function="relative_dens", cond_dist="Gaussian"), 17)
  expect_equal(n_params(p=1, M=2, d=2, weight_function="relative_dens", cond_dist="Gaussian"), 19)
  expect_equal(n_params(p=1, M=3, d=2, weight_function="relative_dens", cond_dist="Gaussian"), 29)
  expect_equal(n_params(p=1, M=1, d=3, weight_function="relative_dens", cond_dist="Gaussian"), 18)
  expect_equal(n_params(p=2, M=1, d=3, weight_function="relative_dens", cond_dist="Gaussian"), 27)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="relative_dens", cond_dist="Gaussian"), 37)
})

