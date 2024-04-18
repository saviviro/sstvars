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

## weight_function = "logistic"

# p=1, M=2, d=2, weightfun_pars=c(1, 1)
c_and_gamma_122_1_1 <- c(0.1, 0.2)
theta_122logistic_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), c_and_gamma_122_1_1)

# p=1, M=2, d=2, weightfun_pars=c(2, 1)
c_and_gamma_122_2_1 <- c(0.11, 0.22)
theta_122logistic_2_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), c_and_gamma_122_2_1)

# p=2, M=2, d=2, weightfun_pars=c(2, 1)
c_and_gamma_222_2_1 <- c(0.1, 0.2)
theta_222logistic_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                           vech(Omega1_222), vech(Omega2_222), c_and_gamma_222_2_1)

# p=2, M=2, d=2, weightfun_pars=c(1, 2)
c_and_gamma_222_1_2 <- c(0.11, 0.22)
theta_222logistic_1_2 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                           vech(Omega1_222), vech(Omega2_222), c_and_gamma_222_1_2)

# p=1, M=2, d=3, weightfun_pars=c(1, 1)
c_and_gamma_123_1_1 <- c(0.1, 0.5)
theta_123logistic_1_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vech(Omega1_123),
                           vech(Omega2_123), c_and_gamma_123_1_1)

# p=1, M=2, d=3, weightfun_pars=c(3, 1)
c_and_gamma_123_3_1 <- c(0.1, 0.4)
theta_123logistic_3_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vech(Omega1_123),
                           vech(Omega2_123), c_and_gamma_123_3_1)


## weight_function = "mlogit"

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

## weight_function == "exponential"

# p=1, M=2, d=2, weight_function="exponential", weightfun_pars=c(1, 1)
theta_122exp_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), c_and_gamma_122_1_1)

# p=1, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1)
theta_122exp_2_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), c_and_gamma_122_2_1)

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1)
theta_222exp_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                      vech(Omega1_222), vech(Omega2_222), c_and_gamma_222_2_1)

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(1, 2)
theta_222exp_1_2 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                      vech(Omega1_222), vech(Omega2_222), c_and_gamma_222_1_2)

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1)
theta_123exp_1_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vech(Omega1_123),
                      vech(Omega2_123), c_and_gamma_123_1_1)

## weight_function == "threshold"

# p=1, M=2, d=2, weight_function="threshold", weightfun_pars=c(1, 1)
r1_122_1_1 <- c(0.5)
theta_122thres_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), r1_122_1_1)

# p=2, M=2, d=2, weight_function="threshold", weightfun_pars=c(2, 2)
r1_222_2_2 <- c(0.7)
theta_222thres_2_2 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                        vech(Omega1_222), vech(Omega2_222), r1_222_2_2)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1)
r1_132_1_1 <- 0.5
r2_132_1_1 <- 1.2
theta_132thres_1_1 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                        vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), r1_132_1_1, r2_132_1_1)

# p=2, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1)
phi10_232 <- phi10_132; phi20_232 <- phi20_132; phi30_232 <- phi30_132
A11_232 <- A11_132; A12_232 <- -A11_132
A21_232 <- A21_132; A22_232 <- -A21_132
A31_232 <- A31_132; A32_232 <- -A31_132
Omega1_232 <- Omega1_132; Omega2_232 <- Omega1_132; Omega3_232 <- Omega3_132
r1_232_1_1 <- r1_132_1_1; r2_232_1_1 <- r2_132_1_1
theta_232thres_1_1 <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                        vec(A31_232), vec(A32_232), vech(Omega1_232), vech(Omega2_232), vech(Omega3_232),
                        r1_232_1_1, r2_232_1_1)


# p=1, M=2, d=3, weight_function="threshold", weightfun_pars=c(2, 1)
r1_123_2_1 <- 1
theta_123thres_2_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vech(Omega1_123),
                        vech(Omega2_123), r1_123_2_1)


## weight_function == "exogenous"

# p=1, M=2, d=3, weight_function="exogenous", weightfun_pars=cbind(c(0.4, 0.2, 0.9), c(0.6, 0.8, 0.1))
theta_123exo <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vech(Omega1_123), vech(Omega2_123))

# p=2, M=2, d=2,  weight_function="exogenous", weightfun_pars=cbind(c(0.4, 0, 0.9), c(0.6, 1, 0.1))
theta_222exo <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222), vech(Omega1_222), vech(Omega2_222))

# p=1, M=3, d=2, weight_function="exogenous", weightfun_pars=cbind(c(0.4, 0, 0.5), c(0.3, 1, 0.2), c(0.3, 0, 0.3))
theta_132exo <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                  vech(Omega2_132), vech(Omega3_132))

# p=2, M=3, d=2, weight_function="exogenous", weightfun_pars=cbind(c(0.4, 0, 0.5), c(0.3, 1, 0.2), c(0.3, 0, 0.3))
theta_232exo <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                  vec(A31_232), vec(A32_232), vech(Omega1_232), vech(Omega2_232), vech(Omega3_232))


## cond_dist == "Student"

# p=2, M=2, d=2, cond_dist="Student", weight_function="logistic", weightfun_pars=c(2, 1)
df_222_2_1 <- 3
theta_222logistict_2_1 <- c(theta_222logistic_2_1, df_222_2_1)

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="Student"
df_122_1_1 <- 13
theta_122logt_1_1 <- c(theta_122log_1_1, df_122_1_1)

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="Student"
df_123_1_1 <- 10
theta_123expt_1_1 <- c(theta_123exp_1_1, df_123_1_1)

# p=2, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="Student"
df_232_1_1 <- 30
theta_232threst_1_1 <- c(theta_232thres_1_1, df_232_1_1)

### Constrained models

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

# p=1, M=3, d=2, AR matrices identical across the regimes
C_132 <- rbind_diags(p=1, M=3, d=2)
theta_132relgc <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vech(Omega1_132), vech(Omega2_132),
                    vech(Omega3_132), alpha1_132, alpha2_132)
theta_132relgc_expanded <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A11_132), vec(A11_132),
                             vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), alpha1_132, alpha2_132)

# p=1, M=1, d=3, AR_constraints
C_113 <- matrix(c(1, 0, -0.3, 0, 0.1, 0.35, 0, -1, 0.14, 0, -0.2, 0, 0.17, 0.12, -0.13, 0, 1, 0.2,
                  0.2, 0.1, 0, 0, -0.11, 0.32, 0.04, 0, 0.1), nrow=9, ncol=3)
psi_113 <- c(0.11, 0.17, -0.3)
theta_113relgc <- c(phi10_113, psi_113, vech(Omega1_113))
theta_113relgc_expanded <- c(phi10_113, C_113%*%psi_113, vech(Omega1_113))

# p=2, M=1, p=3, AR_constraints
C_213 <- matrix(c(0.2, 0.1, 1, 0, -0.1, -0.12, 0, 0, 3, -0.17, 0, 1.1, 0, -0.12, 0, 0.16, -0.3, 1), nrow=18, ncol=1)
psi_213 <- 0.7
theta_213relgc <- c(phi10_213, psi_213, vech(Omega1_213))
theta_213relgc_expanded <- c(phi10_213, C_213%*%psi_213, vech(Omega1_213))

# p=1, M=2, p=3, AR_constraints
C_123 <- rbind_diags(p=1, M=2, d=3)
theta_123relgc <- c(phi10_123, phi20_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), alpha1_123)
theta_123relgc_expanded <- c(phi10_123, phi20_123, vec(A11_123), vec(A11_123), vech(Omega1_123), vech(Omega2_123), alpha1_123)


## weight_function = "logistic"

# p=1, M=2, d=2, weightfun_pars=c(1, 1), AR_constraints=C_122
theta_122logisticc_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vech(Omega1_122), vech(Omega2_122), c_and_gamma_122_1_1)
theta_122logisticc_1_1_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A11_122), vech(Omega1_122), vech(Omega2_122), c_and_gamma_122_1_1)

# p=2, M=2, d=2, weightfun_pars=c(2, 1), AR_constraints=C_222
theta_222logisticc_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), c_and_gamma_222_2_1)
theta_222logisticc_2_1_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                vech(Omega1_222), vech(Omega2_222), c_and_gamma_222_2_1)

# p=1, M=2, d=3, weightfun_pars=c(3, 1), AR_constraints=C_123
theta_123logisticc_3_1 <- c(phi10_123, phi20_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), c_and_gamma_123_3_1)
theta_123logisticc_3_1_expanded <- c(phi10_123, phi20_123, vec(A11_123), vec(A11_123), vech(Omega1_123), vech(Omega2_123), c_and_gamma_123_3_1)


## weight_function = "mlogit"

# p=1, M=2, d=2, weightfun_pars=c(1, 1), C_122
gamma1_122_1_1 <- c(0.1, 0.2)
theta_122logc_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vech(Omega1_122), vech(Omega2_122), c_and_gamma_122_1_1)
theta_122logc_1_1_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A11_122), vech(Omega1_122), vech(Omega2_122), gamma1_122_1_1)

# p=1, M=2, d=2, weightfun_pars=list(vars=1:2, lags=1), C_122
gamma1_122_12_1 <- c(0.1, 0.2, 0.3)
theta_122logc_12_1 <- c(phi10_122, phi20_122, vec(A11_122), vech(Omega1_122), vech(Omega2_122), gamma1_122_12_1)
theta_122logc_12_1_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A11_122), vech(Omega1_122), vech(Omega2_122), gamma1_122_12_1)

# p=2, M=2, d=2, weightfun_pars=list(vars=2, lags=1), C_222
gamma1_222_2_1 <- c(0.1, 0.2)
theta_222logc_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), gamma1_222_2_1)
theta_222logc_2_1_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                vech(Omega1_222), vech(Omega2_222), gamma1_222_2_1)

# p=2, M=2, d=2, weightfun_pars=list(vars=1:2, lags=2), C_222
gamma1_222_12_2 <- c(0.1, 0.2, 0.11, 0.22, 0.33)
theta_222logc_12_2 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), gamma1_222_12_2)
theta_222logc_12_2_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                 vech(Omega1_222), vech(Omega2_222), gamma1_222_12_2)

# p=1, M=2, d=3, weightfun_pars=list(vars=1, lags=1), C_123
gamma1_123_1_1 <- c(0.1, 0.2)
theta_123logc_1_1 <- c(phi10_123, phi20_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), gamma1_123_1_1)
theta_123logc_1_1_expanded <- c(phi10_123, phi20_123, vec(A11_123), vec(A11_123), vech(Omega1_123), vech(Omega2_123), gamma1_123_1_1)

# p=1, M=2, d=3, weightfun_pars=list(vars=2:3, lags=1), C_123
gamma1_123_23_1 <- c(0.1, 0.2, 0.3)
theta_123logc_23_1 <- c(phi10_123, phi20_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), gamma1_123_23_1)
theta_123logc_23_1_expanded <- c(phi10_123, phi20_123, vec(A11_123), vec(A11_123), vech(Omega1_123), vech(Omega2_123), gamma1_123_23_1)

# p=1, M=2, d=3, weightfun_pars=list(vars=1:3, lags=1), AR matrices identical in both regimes
gamma1_123_123_1 <- c(0.1, 0.2, 0.3, 0.4)
theta_123logc_123_1 <- c(phi10_123, phi20_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), gamma1_123_123_1)
theta_123logc_123_1_expanded <- c(phi10_123, phi20_123, vec(A11_123), vec(A11_123), vech(Omega1_123), vech(Omega2_123), gamma1_123_123_1)

## weight_function == "exponential"

# p=2, M=2, d=2, weightfun_pars=c(2, 1), AR_constraints=C_222
theta_222expc_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), c_and_gamma_222_2_1)
theta_222expc_2_1_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                     vech(Omega1_222), vech(Omega2_222), c_and_gamma_222_2_1)

# p=1, M=2, d=3, weightfun_pars=c(3, 1), AR_constraints=C_123
theta_123expc_3_1 <- c(phi10_123, phi20_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), c_and_gamma_123_3_1)
theta_123expc_3_1_expanded <- c(phi10_123, phi20_123, vec(A11_123), vec(A11_123), vech(Omega1_123), vech(Omega2_123), c_and_gamma_123_3_1)

## weight_function == "threshold"

# p=2, M=2, d=2, weight_function="threshold", weightfun_pars=c(1, 1), AR_constraints=C_222
theta_222thresc_1_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), r1_222_2_2 )
theta_222thresc_1_1_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                  vech(Omega1_222), vech(Omega2_222), r1_222_2_2 )

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1,1), AR_constraints=C_132
theta_132thresc_1_1 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vech(Omega1_132), vech(Omega2_132),
                    vech(Omega3_132), r1_132_1_1, r2_132_1_1)
theta_132thresc_1_1_expanded <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A11_132), vec(A11_132),
                              vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), r1_132_1_1, r2_132_1_1)


## weight_function == "exogenous"

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.6, 0.3), c(0, 0.4, 07)), AR_constraints=C_222
theta_222exoc <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222))
theta_222exoc_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                            vech(Omega1_222), vech(Omega2_222))

# p=1, M=3, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.6, 0.0), c(0, 0.3, 07), c(0, 0.1, 0.3)), AR_constraints=C_132
theta_132exoc <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vech(Omega1_132), vech(Omega2_132), vech(Omega3_132))
theta_132exoc_expanded <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A11_132), vec(A11_132),
                            vech(Omega1_132), vech(Omega2_132), vech(Omega3_132))



## Models with mean_constraints

# p=1, M=1, p=2, mean_constraints=list(1)
theta_112relgm <- theta_112relg
theta_112relgm_expanded <- theta_112relgm

# p=1, M=1, p=2, mean_constraints=list(1), C_112
theta_112relgcm <- theta_112relgc
theta_112relgcm_expanded <- theta_112relgc_expanded

# p=1, M=2, d=2, mean_constraints=list(1:2)
theta_122relgm <- c(phi10_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)
theta_122relgm_expanded <- c(phi10_122, phi10_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)

# p=1, M=2, d=2, mean_constraints=list(1, 2)
theta_122relgm2 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)
theta_122relgm2_expanded <- theta_122relgm2

# p=1, M=2, d=2, mean_constraints=list(1:2), C_122
theta_122relgcm <- c(phi10_122, vec(A11_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)
theta_122relgcm_expanded <- c(phi10_122, phi10_122, vec(A11_122), vec(A11_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)

# p=2, M=2, d=2, mean_constraints=list(1:2), C_222
theta_222relgcm <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222)
theta_222relgcm_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                             vech(Omega1_222), vech(Omega2_222), alpha1_222)

# p=1, M=3, d=2, mean_constraints=list(1, 2:3)
theta_132relgm1 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                     vech(Omega2_132), vech(Omega3_132), alpha1_132, alpha2_132)
theta_132relgm1_expanded <- c(phi10_132, phi20_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                              vech(Omega2_132), vech(Omega3_132), alpha1_132, alpha2_132)

# p=1, M=3, d=2, mean_constraints=list(1:2, 3)
theta_132relgm2 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                     vech(Omega2_132), vech(Omega3_132), alpha1_132, alpha2_132)
theta_132relgm2_expanded <- c(phi10_132, phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                              vech(Omega2_132), vech(Omega3_132), alpha1_132, alpha2_132)

# p=1, M=3, d=2, mean_constraints=list(c(1, 3), 2)
theta_132relgm3 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                     vech(Omega2_132), vech(Omega3_132), alpha1_132, alpha2_132)
theta_132relgm3_expanded <- c(phi10_132, phi20_132, phi10_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                              vech(Omega2_132), vech(Omega3_132), alpha1_132, alpha2_132)

# p=1, M=3, d=2, mean_constraints=list(2, c(1, 3))
theta_132relgm4 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                     vech(Omega2_132), vech(Omega3_132), alpha1_132, alpha2_132)
theta_132relgm4_expanded <- c(phi20_132, phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                              vech(Omega2_132), vech(Omega3_132), alpha1_132, alpha2_132)

# p=1, M=3, d=2, mean_constraints=list(1:3), C_132
theta_132relgcm <- c(phi10_132, vec(A11_132), vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), alpha1_132, alpha2_132)
theta_132relgcm_expanded <- c(phi10_132, phi10_132, phi10_132, vec(A11_132), vec(A11_132), vec(A11_132),
                              vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), alpha1_132, alpha2_132)

# p=1, M=1, p=3, mean_constraints=list(1)
theta_113relgm <- theta_113relg
theta_113relgm_expanded <- theta_113relgm

# p=1, M=2, p=3, mean_constraints=list(1:2)
theta_123relgm <- c(phi10_123, vec(A11_123), vec(A21_123), vech(Omega1_123), vech(Omega2_123), alpha1_123)
theta_123relgm_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A21_123), vech(Omega1_123), vech(Omega2_123), alpha1_123)

# p=1, M=2, p=3, mean_constraints=list(1:2), C_123
theta_123relgcm <- c(phi10_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), alpha1_123)
theta_123relgcm_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vech(Omega1_123), vech(Omega2_123), alpha1_123)


## weight_function == "logistic"

# p=2, M=2, d=2, weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222
theta_222logisticcm_2_1 <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), c_and_gamma_222_2_1)
theta_222logisticcm_2_1_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                       vech(Omega1_222), vech(Omega2_222), c_and_gamma_222_2_1)

# p=1, M=2, p=3, weightfun_pars=c(3, 1), mean_constraints=list(1:2), AR_constraints=C_123
theta_123logisticcm_3_1 <- c(phi10_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), c_and_gamma_123_3_1)
theta_123logisticcm_3_1_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vech(Omega1_123), vech(Omega2_123), c_and_gamma_123_3_1)


## weight_function == "mlogit"

# p=1, M=2, d=2, weigthfun_pars=list(vars=1, lags=1), mean_constraints=list(1:2)
theta_122logm_1_1 <- c(phi10_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), gamma1_122_1_1)
theta_122logm_1_1_expanded <- c(phi10_122, phi10_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), gamma1_122_1_1)

# p=2, M=2, d=2, weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), C_222
theta_222logcm_12_2 <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), gamma1_222_12_2)
theta_222logcm_12_2_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                 vech(Omega1_222), vech(Omega2_222), gamma1_222_12_2)

# p=1, M=2, p=3, weightfun_pars=list(vars=1:3, lags=1), mean_constraints=list(1:2), C_123
theta_123logcm_123_1 <- c(phi10_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), gamma1_123_123_1)
theta_123logcm_123_1_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vech(Omega1_123), vech(Omega2_123), gamma1_123_123_1)

## weight_function = "exponential"

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222
theta_222expcm_2_1 <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), c_and_gamma_222_2_1)
theta_222expcm_2_1_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                      vech(Omega1_222), vech(Omega2_222), c_and_gamma_222_2_1)

# p=1, M=2, p=3, weight_function="exponential", weightfun_pars=c(3, 1), mean_constraints=list(1:2), AR_constraints=C_123
theta_123expcm_3_1 <- c(phi10_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), c_and_gamma_123_3_1)
theta_123expcm_3_1_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vech(Omega1_123), vech(Omega2_123), c_and_gamma_123_3_1)


## weight_function = "threshold"

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), mean_constraints=list(1, 2:3)
theta_132thresm_1_1 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                     vech(Omega2_132), vech(Omega3_132), r1_132_1_1, r2_132_1_1)
theta_132thresm_1_1_expanded <- c(phi10_132, phi20_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                                  vech(Omega2_132), vech(Omega3_132), r1_132_1_1, r2_132_1_1)

# p=1, M=2, p=3, weight_function="threshold", weightfun_pars=c(2, 1), AR_constraints=C_123, mean_connstraints=list(1:2)
theta_123threscm_2_1 <- c(phi10_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123), r1_123_2_1)
theta_123threscm_2_1_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vech(Omega1_123), vech(Omega2_123), r1_123_2_1)

## weight_function == "exogenous"
# p=1, M=2, p=3, weight_function="exogenous", weightfun_pars=cbind(c(0, 0.4, 0.9), c(1, 0.6, 0.1)), mean_constraints=list(1:2),
# AR_constraints=C_123
theta_123exocm <- c(phi10_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123))
theta_123exocm_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vech(Omega1_123), vech(Omega2_123))

## Models with weight_constraints

# p=1, M=3, d=2, weight_function="relative_dens", weight_constraints=list(R=matrix(c(0.9, 0.5), nrow=2), r=c(0.13, 0.13))
xi_132relgw <- 0.4
theta_132relgw <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                   vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), xi_132relgw)
theta_132relgw_expanded <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                             vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), matrix(c(0.9, 0.5), nrow=2)%*%xi_132relgw + c(0.13, 0.13))

# p=2, M=2, d=2, weight_function="relative_dens", mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=0, r=0.6)
theta_222relgcmw <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222)) # No weight param since replaced with r
theta_222relgcmw_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                              vech(Omega1_222), vech(Omega2_222), 0.6)

## mlogit

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), weight_constraints=list(R=0, r=c(0.12, 0.13))
theta_122logw_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122))
theta_122logw_1_1_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), c(0.12, 0.13))

# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5), r=c(0, 0.11, 0.12, 0.13, 0))
xi_222logcmw_12_2 <- c(0.22, 0.33)
theta_222logcmw_12_2 <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), xi_222logcmw_12_2 )
theta_222logcmw_12_2_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                  vech(Omega1_222), vech(Omega2_222), c(0.22, 0.11, 0.12, 0.13, 0.33))

## logistic

# p=1, M=2, d=2, weight_function="logistic", weightfun_pars=c(1, 1), weight_constraints=list(R=0, r=c(0.02, 0.13))
theta_122logisticw_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122))
theta_122logisticw_1_1_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), c(0.02, 0.13))

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))
xi_222logisticcmw_2_1 <- c(0.33)
theta_222logisticcmw_2_1 <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), xi_222logisticcmw_2_1)
theta_222logisticcmw_2_1_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                   vech(Omega1_222), vech(Omega2_222), c(0.01, 0.33))

## exponential

# p=1, M=2, d=2, weight_function="exponential", weightfun_pars=c(1, 1), weight_constraints=list(R=0, r=c(0.02, 0.13))
theta_122expw_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122))
theta_122expw_1_1_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), c(0.02, 0.13))

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))
xi_222expcmw_2_1 <- c(0.33)
theta_222expcmw_2_1 <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), xi_222expcmw_2_1)
theta_222expcmw_2_1_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                       vech(Omega1_222), vech(Omega2_222), c(0.01, 0.33))

## threshold

# p=2, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.1, 0))
xi_232thres_1_1 <- 0.5
theta_232thresw_1_1 <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                        vec(A31_232), vec(A32_232), vech(Omega1_232), vech(Omega2_232), vech(Omega3_232), xi_232thres_1_1)
theta_232thresw_1_1_expanded <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                                  vec(A31_232), vec(A32_232), vech(Omega1_232), vech(Omega2_232), vech(Omega3_232), 0.1, 0.5)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))
theta_132thresmw_1_1 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                         vech(Omega2_132), vech(Omega3_132))
theta_132thresmw_1_1_expanded <- c(phi10_132, phi20_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                                  vech(Omega2_132), vech(Omega3_132), 0, 1.2)


## All kinds of constraints, cond_dist == "Student"

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="Student",
# mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))
df_122_2_1 <- 4
theta_222expcmwt_2_1 <- c(theta_222expcmw_2_1, df_122_2_1)
theta_222expcmwt_2_1_expanded <- c(theta_222expcmw_2_1_expanded, df_122_2_1)

# p=2, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="Student",
# weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.1, 0))
theta_232threswt_1_1 <- c(theta_232thresw_1_1, df_232_1_1)
theta_232threswt_1_1_expanded <- c(theta_232thresw_1_1_expanded, df_232_1_1)

# p=1, M=2, p=3, weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="Student",
# mean_constraints=list(1:2), AR_constraints=C_123
df_123_3_1 <- 11
theta_123logisticcmt_3_1 <- c(theta_123logisticcm_3_1, df_123_3_1)
theta_123logisticcmt_3_1_expanded <- c(theta_123logisticcm_3_1_expanded, df_123_3_1)

# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="Student", AR_constraints=C_222
theta_222logct_2_1 <- c(theta_222logc_2_1, df_222_2_1)
theta_222logct_expanded <- c(theta_222logc_2_1_expanded , df_222_2_1)

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="Student", AR_constraints=C_222
theta_222exoct <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), df_222_2_1)
theta_222exoct_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                            vech(Omega1_222), vech(Omega2_222), df_222_2_1)



###############
### ind_Student

# p=1, M=1, d=2, cond_dist="ind_Student", weight_function="threshold", weightfun_pars=c(1, 1)
dfs_112 <- c(3, 7)
B1_112 <- matrix(c(0.5, 0.2, -0.7, 0.3), nrow=2)
theta_112it <- c(phi10_112, vec(A11_112), vec(B1_112), dfs_112)

# p=2, M=2, d=2, cond_dist="ind_Student", weight_function="logistic", weightfun_pars=c(2, 1)
dfs_222_2_1 <- c(3, 7)
B1_222 <- matrix(c(0.5, 0.2, -0.1, 0.3), nrow=2)
B2_222 <- matrix(c(0.4, -0.1, -0.2, 0.3), nrow=2)
theta_222logistit_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                           vec(B1_222), vec(B2_222), c_and_gamma_222_2_1, dfs_222_2_1)

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="ind_Student"
dfs_122_1_1 <- c(4, 13)
B1_122 <- matrix(c(1.2, -0.3, 0.7, 0.1), nrow=2)
B2_122 <- matrix(c(0.5, 0.2, -0.1, 3.1), nrow=2)
theta_122logit_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(B1_122), vec(B2_122), gamma1_122_1_1, dfs_122_1_1)

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="ind_Student"
dfs_123_1_1 <- c(10, 12, 3)
B1_123 <- matrix(c(1.0, 0.3, 0.1, -0.8, 1.1, -0.5, -0.1, -0.2, 0.4), nrow=3)
B2_123 <- matrix(c(0.3, -0.2, -0.7, -0.8, 1.2, 0.5, 0.1, -0.2, 1.1), nrow=3)
theta_123expit_1_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(B1_123),
                        vec(B2_123), c_and_gamma_123_1_1, dfs_123_1_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student"
dfs_132_1_1 <- c(30, 6)
B1_132 <- matrix(c(0.6, 0.2, -0.1, 0.7), nrow=2)
B2_132 <- matrix(c(0.4, -0.1, -0.2, 0.5), nrow=2)
B3_132 <- matrix(c(0.9, -0.5, 0.2, 0.4), nrow=2)
theta_132thresit_1_1 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                          vec(B1_132), vec(B2_132), vec(B3_132), r1_132_1_1, r2_132_1_1, dfs_132_1_1)

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="ind_Student",
# mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))
dfs_122_2_1 <- c(4, 13)
theta_222expcmwit_2_1 <- c(phi10_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222), xi_222expcmw_2_1, dfs_122_2_1)
theta_222expcmwit_2_1_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                    vec(B1_222), vec(B2_222), c(0.01, 0.33), dfs_122_2_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student",
# mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))
theta_132thresmwit_1_1 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(B1_132),
                            vec(B2_132), vec(B3_132), dfs_132_1_1)
theta_132thresmwit_1_1_expanded <- c(phi10_132, phi20_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(B1_132),
                                     vec(B2_132), vec(B3_132), 0, 1.2, dfs_132_1_1)

# p=1, M=2, p=3, weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="ind_Student",
# mean_constraints=list(1:2), AR_constraints=C_123
dfs_123_3_1 <- c(11, 3, 20)
theta_123logisticcmit_3_1 <- c(phi10_123, vec(A11_123), vec(B1_123), vec(B2_123), c_and_gamma_123_3_1, dfs_123_3_1)
theta_123logisticcmit_3_1_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vec(B1_123), vec(B2_123),
                                        c_and_gamma_123_3_1, dfs_123_3_1)

# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_Student", AR_constraints=C_222
theta_222logcit_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222), gamma1_222_2_1, dfs_222_2_1)
theta_222logcit_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                              vec(B1_222), vec(B2_222), gamma1_222_2_1, dfs_222_2_1)

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_Student", AR_constraints=C_222
theta_222exoit <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222), dfs_222_2_1)
theta_222exoit_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                             vec(B1_222), vec(B2_222), dfs_222_2_1)

##########


### Structural models
# (recursively identified models use the same parametrization as reduced form models)

# p=1, M=2, d=2, weight_function="relative_dens", identification="heteroskedasticity"
W_122 <- matrix(c(-0.03, 0.24, -0.76, -0.02), nrow=2, ncol=2, byrow=FALSE)
lambdas_122 <- c(3.36, 0.86)
theta_122relgsh <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), lambdas_122, alpha1_122)

# p=1, M=3, d=2, weight_function="relative_dens", identification="heteroskedasticity"
W_132 <- W_122; lambdas2_132 <- lambdas_122; lambdas3_132 <- c(2.1, 0.62)
theta_132relgsh <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                     vec(W_132), lambdas2_132, lambdas3_132, alpha1_132, alpha2_132)

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student",
# identification="heteroskedasticity"
W_222 <- W_122; lambdas_222 <- lambdas_122
theta_222logistictsh_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                              vec(W_222), lambdas_222, c_and_gamma_222_2_1, df_222_2_1)

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=1), identification="heteroskedasticity"
theta_122logsh_12_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), lambdas_122, gamma1_122_12_1)

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), identification="heteroskedasticity"
W_123 <- matrix(c(-0.47, -0.40, 1.25, 0.58, -1.01, 0.18, -0.66, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE)
lambdas_123 <- c(1.56, 1.44, 0.59)
theta_123expsh_1_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(W_123), lambdas_123, c_and_gamma_123_1_1)

# p=2, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="Student",
# identification="heteroskedasticity"
W_232 <- W_132; lambdas2_232 <- lambdas2_132; lambdas3_232 <- lambdas3_132
theta_232threstsh_1_1 <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                           vec(A31_232), vec(A32_232), vec(W_232), lambdas2_232, lambdas3_232, r1_232_1_1, r2_232_1_1,
                           df_232_1_1)

## Structural models imposing constraints

# p=1, M=2, d=2, weight_function="relative_dens", identification="heteroskedasticity", AR_constraints=C_122
W_122 <- matrix(c(-0.03, 0.24, -0.76, -0.02), nrow=2, ncol=2, byrow=FALSE)
lambdas_122 <- c(3.36, 0.86)
theta_122relgshc <- c(phi10_122, phi20_122, vec(A11_122), vec(W_122), lambdas_122, alpha1_122)
theta_122relgshc_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A11_122), vec(W_122), lambdas_122, alpha1_122)

# p=1, M=2, d=2, weight_function="relative_dens", identification="heteroskedasticity", other_constraints=list(fixed_lambdas=c(3, 1))
theta_122relgshl <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), alpha1_122)
theta_122relgshl_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), 3, 1, alpha1_122)

# p=1, M=2, d=2, weight_function="relative_dens", identification="heteroskedasticity", AR_constraints=C_122,
# other_constraints=list(fixed_lambdas=c(3, 1))
theta_122relgshcl <- c(phi10_122, phi20_122, vec(A11_122), vec(W_122), alpha1_122)
theta_122relgshcl_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A11_122), vec(W_122), 3, 1, alpha1_122)

# p=1, M=3, d=2, weight_function="relative_dens", identification="heteroskedasticity",
# B_constraints=matrix(c(0.1, NA, 0.3, 0), nrow=2)
W_132b <- matrix(c(0.11, 0.22, 0.33, 0), nrow=2);
theta_132relgshb <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                      Wvec(W_132b), lambdas2_132, lambdas3_132, alpha1_132, alpha2_132)
theta_132relgshb_expanded <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                               vec(W_132b), lambdas2_132, lambdas3_132, alpha1_132, alpha2_132)

# p=1, M=3, d=2, weight_function="relative_dens", identification="heteroskedasticity",
# B_constraints=matrix(c(0.1, NA, 0.3, 0), nrow=2), other_constraints=list(fixed_lambdas=c(3, 1, 2, 0.5))
theta_132relgshbl <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                       Wvec(W_132b), alpha1_132, alpha2_132)
theta_132relgshbl_expanded <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                                vec(W_132b), 3, 1, 2, 0.5, alpha1_132, alpha2_132)

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student",
# identification="heteroskedasticity", mean_constraints=list(1:2), B_constraints=matrix(c(0.1, 0, 0, 0.3), nrow=2)
W_222b <- matrix(c(0.12, 0, 0, 0.31), nrow=2)
theta_222logistictshmb_2_1 <- c(phi10_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                                Wvec(W_222b), lambdas_222, c_and_gamma_222_2_1, df_222_2_1)
theta_222logistictshmb_2_1_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                                         vec(W_222b), lambdas_222, c_and_gamma_222_2_1, df_222_2_1)

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student",
# identification="heteroskedasticity", mean_constraints=list(1:2), B_constraints=matrix(c(0.1, 0, 0, 0.3), nrow=2),
# other_constraints=list(fixed_lambdas=c(3, 1))
theta_222logistictshmbl <- c(phi10_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                             Wvec(W_222b), c_and_gamma_222_2_1, df_222_2_1)
theta_222logistictshmbl_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                                      vec(W_222b), 3, 1, c_and_gamma_222_2_1, df_222_2_1)


# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=1), identification="heteroskedasticity",
# weight_constraints=list(R=0, r=c(0.1, 0.2, 0.3)), B_constraints=matrix(c(0.1, 0.2, 0.3, 0), nrow=2)
W_122b <- matrix(c(0.11, 0.22, 0.33, 0), nrow=2)
theta_122logshwb_12_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), Wvec(W_122b), lambdas_122)
theta_122logshwb_12_1_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122b),
                                    lambdas_122, c(0.1, 0.2, 0.3))

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), identification="heteroskedasticity",
# AR_constraints=C_123, weight_constraints=list(R=matrix(c(1, 0.5), nrow=2), r=c(0, 0))
# B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE)
W_123b <- matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE)
lambdas_123 <- c(1.56, 1.44, 0.59)
theta_123expshcwb_1_1 <- c(phi10_123, phi20_123, vec(A11_123), Wvec(W_123b), lambdas_123, 0.6)
theta_123expshcwb_1_1_expanded <- c(phi10_123, phi20_123, vec(A11_123), vec(A11_123), vec(W_123b),
                                    lambdas_123, c(0.6, 0.3))

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), identification="heteroskedasticity",
# AR_constraints=C_123, weight_constraints=list(R=matrix(c(1, 0.5), nrow=2), r=c(0, 0))
# B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE)
# other_constraints=list(fixed_lambdas=c(3, 1, 2))
theta_123expshcwl <- c(phi10_123, phi20_123, vec(A11_123), Wvec(W_123b), 0.6)
theta_123expshcwl_expanded <- c(phi10_123, phi20_123, vec(A11_123), vec(A11_123), vec(W_123b), 3, 1, 2, 0.6, 0.3)

# p=2, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="Student",
# identification="heteroskedasticity", B_constraints=matrix(c(0.1, 0.2, -0.3, 0), nrow=2)
W_232b <- matrix(c(0.1, 0.2, -0.1, 0), nrow=2)
theta_232threstshb_1_1 <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                            vec(A31_232), vec(A32_232), Wvec(W_232b), lambdas2_232, lambdas3_232, r1_232_1_1, r2_232_1_1,
                            df_232_1_1)
theta_232threstshb_1_1_expanded <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232),
                                     vec(A22_232), vec(A31_232), vec(A32_232), vec(W_232b), lambdas2_232,
                                     lambdas3_232, r1_232_1_1, r2_232_1_1, df_232_1_1)

# p=1, M=2, d=3, weight_function="exogenous", weightfun_pars=cbind(c(0.2, 0.8, 0.4), c(0.8, 0.2, 0.6)), identification="heteroskedasticity",
# AR_constraints=C_123, weight_constraints=list(R=matrix(c(1, 0.5), nrow=2), r=c(0, 0))
# B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE)
W_123b <- matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE)
lambdas_123 <- c(1.56, 1.44, 0.59)
theta_123exoshcwb_1_1 <- c(phi10_123, phi20_123, vec(A11_123), Wvec(W_123b), lambdas_123)
theta_123exoshcwb_1_1_expanded <- c(phi10_123, phi20_123, vec(A11_123), vec(A11_123), vec(W_123b), lambdas_123)



###############
### ind_Student

# p=2, M=2, d=2, cond_dist="ind_Student", weight_function="logistic", weightfun_pars=c(2, 1), identification="non-Gaussianity",
# B_constraints=matrix(c(NA, -1, 0, 1), nrow=2)
B1_222c <- matrix(c(0.5, -0.2, 0, 0.1), nrow=2)
B2_222c <- matrix(c(-0.4, -0.1, 0, 0.2), nrow=2)
theta_222logistitngb_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                              Wvec(B1_222c), Wvec(B2_222c), c_and_gamma_222_2_1, dfs_222_2_1)
theta_222logistitngb_2_1_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                                       vec(B1_222c), vec(B2_222c), c_and_gamma_222_2_1, dfs_222_2_1)

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="ind_Student", identification="non-Gaussianity",
# B_constraints=matrix(c(1, NA, -1, 1) , nrow=2)
dfs_122_1_1 <- c(4, 13)
B1_122c <- matrix(c(1.2, 0.3, -0.7, 0.1), nrow=2)
B2_122c <- matrix(c(0.5, -0.9, -0.1, 3.1), nrow=2)
theta_122logitngb_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), Wvec(B1_122c), Wvec(B2_122c), gamma1_122_1_1, dfs_122_1_1)
theta_122logitngb_1_1_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(B1_122c), vec(B2_122c), gamma1_122_1_1, dfs_122_1_1)


# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="ind_Student", identification="non-Gaussianity",
# B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)
dfs_123_1_1 <- c(10, 12, 3)
B1_123c <- matrix(c(1.0, 0.3, 0.1, 0, 1.1, -0.5, 0, -0.2, 0.4), nrow=3)
B2_123c <- matrix(c(0.3, -0.2, -0.7, 0, 1.2, 0.5, 0, -0.2, 1.1), nrow=3)
theta_123expitngb_1_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), Wvec(B1_123c),
                           Wvec(B2_123c), c_and_gamma_123_1_1, dfs_123_1_1)
theta_123expitngb_1_1_expanded <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(B1_123c),
                                    vec(B2_123c), c_and_gamma_123_1_1, dfs_123_1_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student", identification="non-Gaussianity",
# B_constraints=matrix(c(1, 0, NA, 1), nrow=2)
dfs_132_1_1 <- c(30, 6)
B1_132c <- matrix(c(0.6, 0, -0.1, 0.7), nrow=2)
B2_132c <- matrix(c(0.4, 0, 0.2, 0.5), nrow=2)
B3_132c <- matrix(c(0.9, 0, -0.2, 0.4), nrow=2)
theta_132thresitngb_1_1 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                             Wvec(B1_132c), Wvec(B2_132c), Wvec(B3_132c), r1_132_1_1, r2_132_1_1, dfs_132_1_1)
theta_132thresitngb_1_1_expanded <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                                      vec(B1_132c), vec(B2_132c), vec(B3_132c), r1_132_1_1, r2_132_1_1, dfs_132_1_1)

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="ind_Student",
# mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
# identification="non-Gaussianity", B_constraints=matrix(c(NA, -1, 0, 1), nrow=2)
theta_222expcmwitngb_2_1 <- c(phi10_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c), xi_222expcmw_2_1, dfs_222_2_1)
theta_222expcmwitngb_2_1_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                       vec(B1_222c), vec(B2_222c), c(0.01, 0.33), dfs_222_2_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student", mean_constraints=list(1, 2:3),
# weight_constraints=list(R=0, r=c(0, 1.2)), identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2)
theta_132thresmwitngb_1_1 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), Wvec(B1_132c),
                               Wvec(B2_132c), Wvec(B3_132c), dfs_132_1_1)
theta_132thresmwitngb_1_1_expanded <- c(phi10_132, phi20_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(B1_132c),
                                        vec(B2_132c), vec(B3_132c), 0, 1.2, dfs_132_1_1)

# p=1, M=2, d=3, weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="ind_Student", mean_constraints=list(1:2),
# AR_constraints=C_123, identification="non-Gaussianity", B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)
dfs_123_3_1 <- c(11, 3, 20)
theta_123logisticcmitngb_3_1 <- c(phi10_123, vec(A11_123), Wvec(B1_123c), Wvec(B2_123c), c_and_gamma_123_3_1, dfs_123_3_1)
theta_123logisticcmitngb_3_1_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vec(B1_123c), vec(B2_123c),
                                           c_and_gamma_123_3_1, dfs_123_3_1)

# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_Student", AR_constraints=C_222,
# identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)
theta_222logcitngb_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c), gamma1_222_2_1, dfs_222_2_1)
theta_222logcitngb_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                 vec(B1_222c), vec(B2_222c), gamma1_222_2_1, dfs_222_2_1)

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_Student",
# AR_constraints=C_222, identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)
theta_222exoitngb_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c), dfs_222_2_1)
theta_222exoitngb_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                vec(B1_222c), vec(B2_222c), dfs_222_2_1)

###############


test_that("reform_constrained_pars works correctly", {
  # Models with AR_constraints
  expect_equal(reform_constrained_pars(p=1, M=1, d=2, params=theta_112relgc, weight_function="relative_dens",
                                       AR_constraints=C_112), theta_112relgc_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122relgc, weight_function="relative_dens",
                                       AR_constraints=C_122), theta_122relgc_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222relgc, weight_function="relative_dens",
                                       AR_constraints=C_222), theta_222relgc_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222relgc2, weight_function="relative_dens",
                                       AR_constraints=C_222_2), theta_222relgc2_expanded)
  expect_equal(reform_constrained_pars(p=1, M=3, d=2, params=theta_132relgc, weight_function="relative_dens",
                                       AR_constraints=C_132), theta_132relgc_expanded)
  expect_equal(reform_constrained_pars(p=1, M=1, d=3, params=theta_113relgc, weight_function="relative_dens",
                                       AR_constraints=C_113), theta_113relgc_expanded)
  expect_equal(reform_constrained_pars(p=2, M=1, d=3, params=theta_213relgc, weight_function="relative_dens",
                                       AR_constraints=C_213), theta_213relgc_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123relgc, weight_function="relative_dens",
                                       AR_constraints=C_123), theta_123relgc_expanded)

  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122logc_1_1, weight_function="mlogit",
                                       weightfun_pars=list(vars=1, lags=1), AR_constraints=C_122), theta_122logc_1_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122logc_12_1, weight_function="mlogit",
                                       weightfun_pars=list(vars=1:2, lags=1), AR_constraints=C_122), theta_122logc_12_1_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222logc_2_1, weight_function="mlogit",
                                       weightfun_pars=list(vars=2, lags=1), AR_constraints=C_222), theta_222logc_2_1_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222logc_12_2, weight_function="mlogit",
                                       weightfun_pars=list(vars=1:2, lags=2), AR_constraints=C_222), theta_222logc_12_2_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123logc_1_1, weight_function="mlogit",
                                       weightfun_pars=list(vars=1, lags=1), AR_constraints=C_123), theta_123logc_1_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123logc_23_1, weight_function="mlogit",
                                       weightfun_pars=list(vars=2:3, lags=1), AR_constraints=C_123), theta_123logc_23_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123logc_123_1, weight_function="mlogit",
                                       weightfun_pars=list(vars=1:3, lags=1), AR_constraints=C_123), theta_123logc_123_1_expanded)

  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122logisticc_1_1, weight_function="logistic",
                                       weightfun_pars=c(1, 1), AR_constraints=C_122), theta_122logisticc_1_1_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222logisticc_2_1, weight_function="logistic",
                                       weightfun_pars=c(2, 1), AR_constraints=C_222), theta_222logisticc_2_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123logisticc_3_1, weight_function="logistic",
                                       weightfun_pars=c(3, 1), AR_constraints=C_123), theta_123logisticc_3_1_expanded)

  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222expc_2_1, weight_function="exponential",
                                       weightfun_pars=c(2, 1), AR_constraints=C_222), theta_222expc_2_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123expc_3_1, weight_function="exponential",
                                       weightfun_pars=c(3, 1), AR_constraints=C_123), theta_123expc_3_1_expanded)

  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222thresc_1_1, weight_function="threshold",
                                       weightfun_pars=c(1, 1), AR_constraints=C_222), theta_222thresc_1_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=3, d=2, params=theta_132thresc_1_1, weight_function="threshold",
                                       weightfun_pars=c(1, 1), AR_constraints=C_132), theta_132thresc_1_1_expanded)

  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222exoc, weight_function="exogenous",
                                       weightfun_pars=cbind(c(1, 0.6, 0.3), c(0, 0.4, 07)), AR_constraints=C_222),
               theta_222exoc_expanded)
  expect_equal(reform_constrained_pars(p=1, M=3, d=2, params=theta_132exoc, weight_function="exogenous",
                                       weightfun_pars=cbind(c(1, 0.6, 0.0), c(0, 0.3, 07), c(0, 0.1, 0.3)), AR_constraints=C_132),
               theta_132exoc_expanded)

  # Models with mean_constraints
  expect_equal(reform_constrained_pars(p=1, M=1, d=2, params=theta_112relgm, weight_function="relative_dens",
                                       mean_constraints=list(1)), theta_112relgm_expanded)
  expect_equal(reform_constrained_pars(p=1, M=1, d=2, params=theta_112relgcm, weight_function="relative_dens",
                                       AR_constraints=C_112, mean_constraints=list(1)), theta_112relgc_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122relgm, weight_function="relative_dens",
                                       mean_constraints=list(1:2)), theta_122relgm_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122relgm2, weight_function="relative_dens",
                                       mean_constraints=list(1, 2)), theta_122relgm2_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122relgcm, weight_function="relative_dens",
                                       AR_constraints=C_122, mean_constraints=list(1:2)), theta_122relgcm_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222relgcm, weight_function="relative_dens",
                                       AR_constraints=C_222, mean_constraints=list(1:2)), theta_222relgcm_expanded)
  expect_equal(reform_constrained_pars(p=1, M=3, d=2, params=theta_132relgm1, weight_function="relative_dens",
                                       mean_constraints=list(1, 2:3)), theta_132relgm1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=3, d=2, params=theta_132relgm2, weight_function="relative_dens",
                                       mean_constraints=list(1:2, 3)), theta_132relgm2_expanded)
  expect_equal(reform_constrained_pars(p=1, M=3, d=2, params=theta_132relgm3, weight_function="relative_dens",
                                       mean_constraints=list(c(1, 3), 2)), theta_132relgm3_expanded)
  expect_equal(reform_constrained_pars(p=1, M=3, d=2, params=theta_132relgm4, weight_function="relative_dens",
                                       mean_constraints=list(2, c(1, 3))), theta_132relgm4_expanded)
  expect_equal(reform_constrained_pars(p=1, M=3, d=2, params=theta_132relgcm, weight_function="relative_dens",
                                       AR_constraints=C_132, mean_constraints=list(1:3)), theta_132relgcm_expanded)
  expect_equal(reform_constrained_pars(p=1, M=1, d=3, params=theta_113relgm, weight_function="relative_dens",
                                       mean_constraints=list(1)), theta_113relgm_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123relgm, weight_function="relative_dens",
                                       mean_constraints=list(1:2)), theta_123relgm_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123relgcm, weight_function="relative_dens",
                                       AR_constraints=C_123, mean_constraints=list(1:2)), theta_123relgcm_expanded)

  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122logm_1_1, weight_function="mlogit",
                                       weightfun_pars=list(vars=1, lags=1), mean_constraints=list(1:2)), theta_122logm_1_1_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222logcm_12_2, weight_function="mlogit",
                                       weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222),
               theta_222logcm_12_2_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123logcm_123_1, weight_function="mlogit",
                                       weightfun_pars=list(vars=1:3, lags=1), mean_constraints=list(1:2), AR_constraints=C_123),
               theta_123logcm_123_1_expanded)

  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222logisticcm_2_1, weight_function="logistic",
                                       weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222),
               theta_222logisticcm_2_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123logisticcm_3_1, weight_function="logistic",
                                       weightfun_pars=c(3, 1), mean_constraints=list(1:2), AR_constraints=C_123),
               theta_123logisticcm_3_1_expanded)

  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222expcm_2_1, weight_function="exponential",
                                       weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222),
               theta_222expcm_2_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123expcm_3_1, weight_function="exponential",
                                       weightfun_pars=c(3, 1), mean_constraints=list(1:2), AR_constraints=C_123),
               theta_123expcm_3_1_expanded)

  expect_equal(reform_constrained_pars(p=1, M=3, d=2, params=theta_132thresm_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                       mean_constraints=list(1, 2:3)), theta_132thresm_1_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123threscm_2_1, weight_function="threshold", weightfun_pars=c(2, 1),
                                       mean_constraints=list(1:2), AR_constraints=C_123), theta_123threscm_2_1_expanded)

  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123exocm, weight_function="exogenous",
                                       weightfun_pars=cbind(c(0, 0.4, 0.9), c(1, 0.6, 0.1)), mean_constraints=list(1:2),
                                       AR_constraints=C_123), theta_123exocm_expanded)


  # Models with weight_constraints
  expect_equal(reform_constrained_pars(p=1, M=3, d=2, params=theta_132relgw, weight_function="relative_dens",
                                       weight_constraints=list(R=matrix(c(0.9, 0.5), nrow=2), r=c(0.13, 0.13))), theta_132relgw_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222relgcmw, weight_function="relative_dens",
                                       mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=0, r=0.6)),
               theta_222relgcmw_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122logw_1_1, weight_function="mlogit",
                                       weightfun_pars=list(vars=1, lags=1), weight_constraints=list(R=0, r=c(0.12, 0.13))),
               theta_122logw_1_1_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222logcmw_12_2, weight_function="mlogit",
                                       weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222,
                                       weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5),
                                                               r=c(0, 0.11, 0.12, 0.13, 0))),
               theta_222logcmw_12_2_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122logw_1_1, weight_function="logistic", weightfun_pars=c(1, 1),
                                       weight_constraints=list(R=0, r=c(0.02, 0.13))), theta_122logisticw_1_1_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222logisticcmw_2_1, weight_function="logistic", weightfun_pars=c(2, 1),
                                       mean_constraints=list(1:2), AR_constraints=C_222,
                                       weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))), theta_222logisticcmw_2_1_expanded)

  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122expw_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                                       weight_constraints=list(R=0, r=c(0.02, 0.13))), theta_122expw_1_1_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222expcmw_2_1, weight_function="exponential", weightfun_pars=c(2, 1),
                                       mean_constraints=list(1:2), AR_constraints=C_222,
                                       weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))), theta_222expcmw_2_1_expanded)

  expect_equal(reform_constrained_pars(p=2, M=3, d=2, params=theta_232thresw_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                       weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.1, 0))), theta_232thresw_1_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=3, d=2, params=theta_132thresmw_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                       mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))), theta_132thresmw_1_1_expanded)

  expect_equal(reform_constrained_pars(p=1, M=3, d=2, params=theta_132exo, weight_function="exogenous",
                                       weightfun_pars=cbind(c(1, 0.6, 0.0), c(0, 0.3, 07), c(0, 0.1, 0.3)),
                                       weight_constraints=list(R=0, r=c(0, 1.2))), theta_132exo)

  # Student
  expect_equal(reform_constrained_pars(p=2, M=3, d=2, params=theta_232threswt_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                       weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.1, 0)), cond_dist="Student"),
               theta_232threswt_1_1_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222logct_2_1, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1),
                                       cond_dist="Student", AR_constraints=C_222), theta_222logct_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123logisticcmt_3_1, weight_function="logistic", weightfun_pars=c(3, 1),
                                       cond_dist="Student", mean_constraints=list(1:2), AR_constraints=C_123),
               theta_123logisticcmt_3_1_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222expcmwt_2_1, weight_function="exponential", weightfun_pars=c(2, 1),
                                       cond_dist="Student", mean_constraints=list(1:2), AR_constraints=C_222,
                                       weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))), theta_222expcmwt_2_1_expanded)

  # Ind student
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222expcmwit_2_1, weight_function="exponential", weightfun_pars=c(2, 1),
                                       cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_222,
                                       weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))), theta_222expcmwit_2_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=3, d=2, params=theta_132thresmwit_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                       cond_dist="ind_Student", mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))),
               theta_132thresmwit_1_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123logisticcmit_3_1, weight_function="logistic", weightfun_pars=c(3, 1),
                                       cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_123),
               theta_123logisticcmit_3_1_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222logcit_2_1, weight_function="mlogit",
                                       weightfun_pars=list(vars=2, lags=1), cond_dist="ind_Student", AR_constraints=C_222),
               theta_222logcit_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222exoit, weight_function="exogenous",
                                       weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_Student", AR_constraints=C_222),
               theta_222exoit_expanded)

  ## Structural models
  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122relgshc, weight_function="relative_dens",
                                       identification="heteroskedasticity", AR_constraints=C_122),
               theta_122relgshc_expanded)
  expect_equal(reform_constrained_pars(p=1, M=3, d=2, params=theta_132relgshb, weight_function="relative_dens",
                                       identification="heteroskedasticity", B_constraints=matrix(c(0.1, NA, 0.3, 0), nrow=2)),
               theta_132relgshb_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222logistictshmb_2_1, weight_function="logistic", weightfun_pars=c(2, 1),
                                       cond_dist="Student", identification="heteroskedasticity", mean_constraints=list(1:2),
                                       B_constraints=matrix(c(0.1, 0, 0, 0.3), nrow=2)),
               theta_222logistictshmb_2_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122logshwb_12_1, weight_function="mlogit",
                                       weightfun_pars=list(vars=1:2, lags=1), identification="heteroskedasticity",
                                       weight_constraints=list(R=0, r=c(0.1, 0.2, 0.3)),
                                       B_constraints=matrix(c(0.1, 0.2, 0.3, 0), nrow=2)),
               theta_122logshwb_12_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123expshcwb_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                                       identification="heteroskedasticity", AR_constraints=C_123,
                                       weight_constraints=list(R=matrix(c(1, 0.5), nrow=2), r=c(0, 0)),
                                       B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19),
                                                            nrow=3, ncol=3, byrow=FALSE)),
               theta_123expshcwb_1_1_expanded)
  expect_equal(reform_constrained_pars(p=2, M=3, d=2, params=theta_232threstshb_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                       cond_dist="Student", identification="heteroskedasticity",
                                       B_constraints=matrix(c(0.1, 0.2, -0.3, 0), nrow=2)),
               theta_232threstshb_1_1_expanded)

  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122relgshl, weight_function="relative_dens",
                                       identification="heteroskedasticity", other_constraints=list(fixed_lambdas=c(3, 1))),
               theta_122relgshl_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122relgshcl, weight_function="relative_dens",
                                       identification="heteroskedasticity", AR_constraints=C_122,
                                       other_constraints=list(fixed_lambdas=c(3, 1))),
               theta_122relgshcl_expanded)
  expect_equal(reform_constrained_pars(p=1, M=3, d=2, params=theta_132relgshbl, weight_function="relative_dens",
                                       identification="heteroskedasticity", B_constraints=matrix(c(0.1, NA, 0.3, 0), nrow=2),
                                       other_constraints=list(fixed_lambdas=c(3, 1, 2, 0.5))),
               theta_132relgshbl_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123expshcwl, weight_function="exponential", weightfun_pars=c(1, 1),
                                       identification="heteroskedasticity", AR_constraints=C_123,
                                       weight_constraints=list(R=matrix(c(1, 0.5), nrow=2), r=c(0, 0)),
                                       B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19),
                                                            nrow=3, ncol=3, byrow=FALSE),
                                       other_constraints=list(fixed_lambdas=c(3, 1, 2))),
               theta_123expshcwl_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123exoshcwb_1_1, weight_function="exogenous",
                                       weightfun_pars=cbind(c(0.2, 0.8, 0.4), c(0.8, 0.2, 0.6)), identification="heteroskedasticity",
                                       AR_constraints=C_123, weight_constraints=list(R=matrix(c(1, 0.5), nrow=2), r=c(0, 0)),
                                       B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19), nrow=3, ncol=3,
                                                            byrow=FALSE)),
               theta_123exoshcwb_1_1_expanded)

  # ind_Student with B_constraints
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222logistitngb_2_1, cond_dist="ind_Student", weight_function="logistic",
                                       weightfun_pars=c(2, 1), identification="non-Gaussianity", B_constraints=matrix(c(NA, -1, 0, 1), nrow=2)),
               theta_222logistitngb_2_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=2, params=theta_122logitngb_1_1, weight_function="mlogit",
                                       weightfun_pars=list(vars=1, lags=1), cond_dist="ind_Student", identification="non-Gaussianity",
                                       B_constraints=matrix(c(1, NA, -1, 1) , nrow=2)),
               theta_122logitngb_1_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123expitngb_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                                       cond_dist="ind_Student", identification="non-Gaussianity",
                                       B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)),
               theta_123expitngb_1_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=3, d=2, params=theta_132thresitngb_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                       cond_dist="ind_Student", identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2)),
               theta_132thresitngb_1_1_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222expcmwitngb_2_1, weight_function="exponential", weightfun_pars=c(2, 1),
                                       cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_222,
                                       weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
                                       identification="non-Gaussianity", B_constraints=matrix(c(NA, -1, 0, 1), nrow=2)),
               theta_222expcmwitngb_2_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=3, d=2, params=theta_132thresmwitngb_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                       cond_dist="ind_Student", mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2)),
                                       identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2)),
               theta_132thresmwitngb_1_1_expanded)
  expect_equal(reform_constrained_pars(p=1, M=2, d=3, params=theta_123logisticcmitngb_3_1, weight_function="logistic", weightfun_pars=c(3, 1),
                                       cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_123,
                                       identification="non-Gaussianity", B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)),
               theta_123logisticcmitngb_3_1_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222logcitngb_2_1, weight_function="mlogit",
                                       weightfun_pars=list(vars=2, lags=1), cond_dist="ind_Student", AR_constraints=C_222,
                                       identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)),
               theta_222logcitngb_expanded)
  expect_equal(reform_constrained_pars(p=2, M=2, d=2, params=theta_222exoitngb_2_1, weight_function="exogenous",
                                       weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_Student",
                                       AR_constraints=C_222, identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)),
               theta_222exoitngb_expanded)
})




test_that("change_regime works correctly", {
  expect_equal(change_regime(p=1, M=1, d=2, params=theta_112relg, m=1, regime_pars=1:length(theta_112relg)), 1:length(theta_112relg))
  expect_equal(change_regime(p=2, M=1, d=2, params=theta_212relg, m=1, regime_pars=1:length(theta_212relg)), 1:length(theta_212relg))
  expect_equal(change_regime(p=3, M=1, d=2, params=theta_312relg, m=1, regime_pars=1:length(theta_312relg)), 1:length(theta_312relg))
  expect_equal(change_regime(p=1, M=2, d=2, params=theta_122relg, m=1, regime_pars=c(phi10_112, A11_112, vech(Omega1_112))),
               c(phi10_112, phi20_122, A11_112, A21_122, vech(Omega1_112), vech(Omega2_122), alpha1_122))
  expect_equal(change_regime(p=1, M=2, d=2, params=theta_122relg, m=2, regime_pars=c(phi10_112, A11_112, vech(Omega1_112))),
               c(phi10_122, phi10_112, A11_122, A11_112, vech(Omega1_122), vech(Omega1_112), alpha1_122))
  expect_equal(change_regime(p=2, M=2, d=2, params=theta_222relg, m=1, regime_pars=c(phi10_112, A11_112, A11_122, vech(Omega1_112))),
               c(phi10_112, phi20_222, A11_112, A11_122, A21_222, A22_222, vech(Omega1_112), vech(Omega2_222), alpha1_222))
  expect_equal(change_regime(p=2, M=2, d=2, params=theta_222relg, m=2, regime_pars=c(phi10_112, A11_112, A11_122, vech(Omega1_112))),
               c(phi10_222, phi10_112, A11_222, A12_222, A11_112, A11_122, vech(Omega1_222), vech(Omega1_112), alpha1_222))
  expect_equal(change_regime(p=1, M=3, d=2, params=theta_132relg, m=1, regime_pars=c(phi10_112, A11_112, vech(Omega1_112))),
               c(phi10_112, phi20_132, phi30_132, A11_112, A21_132, A31_132, vech(Omega1_112), vech(Omega2_132), vech(Omega3_132),
                 alpha1_132, alpha2_132))
  expect_equal(change_regime(p=1, M=3, d=2, params=theta_132relg, m=2, regime_pars=c(phi10_112, A11_112, vech(Omega1_112))),
               c(phi10_132, phi10_112, phi30_132, A11_132, A11_112, A31_132, vech(Omega1_132), vech(Omega1_112), vech(Omega3_132),
                 alpha1_132, alpha2_132))
  expect_equal(change_regime(p=1, M=3, d=2, params=theta_132relg, m=3, regime_pars=c(phi10_112, A11_112, vech(Omega1_112))),
               c(phi10_132, phi20_132, phi10_112, A11_132, A21_132, A11_112, vech(Omega1_132), vech(Omega2_132), vech(Omega1_112),
                 alpha1_132, alpha2_132))
  expect_equal(change_regime(p=1, M=1, d=3, params=theta_113relg, m=1, regime_pars=1:length(theta_113relg)), 1:length(theta_113relg))
  expect_equal(change_regime(p=2, M=1, d=3, params=theta_113relg, m=1, regime_pars=1:length(theta_213relg)), 1:length(theta_213relg))
  expect_equal(change_regime(p=1, M=2, d=3, params=theta_123relg, m=1, regime_pars=c(1:3, 3:11, 12:17)),
               c(1:3, phi20_123, 3:11, A21_123, 12:17, vech(Omega2_123), alpha1_123))
  expect_equal(change_regime(p=1, M=2, d=3, params=theta_123relg, m=2, regime_pars=c(1:3, 3:11, 12:17)),
               c(phi10_123, 1:3, A11_123, 3:11, vech(Omega1_123), 12:17, alpha1_123))

  expect_equal(change_regime(p=1, M=2, d=2, params=theta_122log_1_1, m=1, regime_pars=c(phi10_112, A11_112, vech(Omega1_112))),
               c(phi10_112, phi20_122, A11_112, A21_122, vech(Omega1_112), vech(Omega2_122), gamma1_122_1_1))
  expect_equal(change_regime(p=2, M=2, d=2, params=theta_222log_12_2, m=2, regime_pars=c(phi10_112, A11_112, A11_122, vech(Omega1_112))),
               c(phi10_222, phi10_112, A11_222, A12_222, A11_112, A11_122, vech(Omega1_222), vech(Omega1_112), gamma1_222_12_2))

  expect_equal(change_regime(p=1, M=2, d=2, params=theta_122logistic_1_1, m=1, regime_pars=c(phi10_112, A11_112, vech(Omega1_112))),
               c(phi10_112, phi20_122, A11_112, A21_122, vech(Omega1_112), vech(Omega2_122), c_and_gamma_122_1_1))
  expect_equal(change_regime(p=2, M=2, d=2, params=theta_222logistic_2_1, m=2, regime_pars=c(phi10_112, A11_112, A11_122, vech(Omega1_112))),
               c(phi10_222, phi10_112, A11_222, A12_222, A11_112, A11_122, vech(Omega1_222), vech(Omega1_112), c_and_gamma_222_2_1))

  expect_equal(change_regime(p=1, M=2, d=2, params=theta_122exp_1_1, m=1, regime_pars=c(phi10_112, A11_112, vech(Omega1_112))),
               c(phi10_112, phi20_122, A11_112, A21_122, vech(Omega1_112), vech(Omega2_122), c_and_gamma_122_1_1))
  expect_equal(change_regime(p=2, M=2, d=2, params=theta_222exp_2_1, m=2, regime_pars=c(phi10_112, A11_112, A11_122, vech(Omega1_112))),
               c(phi10_222, phi10_112, A11_222, A12_222, A11_112, A11_122, vech(Omega1_222), vech(Omega1_112), c_and_gamma_222_2_1))

  expect_equal(change_regime(p=1, M=2, d=2, params=theta_122thres_1_1, m=1, regime_pars=c(phi10_112, A11_112, vech(Omega1_112))),
               c(phi10_112, phi20_122, A11_112, A21_122, vech(Omega1_112), vech(Omega2_122), r1_122_1_1))
  expect_equal(change_regime(p=2, M=3, d=2, params=theta_232thres_1_1, m=2, regime_pars=c(phi10_232, A11_232, A12_232, vech(Omega1_232))),
               c(phi10_232, phi10_232, phi30_232, vec(A11_232), vec(A12_232), vec(A11_232), vec(A12_232),
                 vec(A31_232), vec(A32_232), vech(Omega1_232), vech(Omega1_232), vech(Omega3_232),
                 r1_232_1_1, r2_232_1_1))

  expect_equal(change_regime(p=1, M=2, d=3, params=theta_123exo, m=1, regime_pars=c(phi20_123, A21_123, vech(Omega2_123))),
               c(phi20_123, phi20_123, vec(A21_123), vec(A21_123), vech(Omega2_123), vech(Omega2_123)))
  expect_equal(change_regime(p=2, M=3, d=2, params=theta_232exo, m=3, regime_pars=c(phi10_232, A11_232, A12_232, vech(Omega1_232))),
               c(phi10_232, phi20_232, phi10_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                 vec(A11_232), vec(A12_232), vech(Omega1_232), vech(Omega2_232), vech(Omega1_232)))

  # Student
  expect_equal(change_regime(p=2, M=2, d=2, params=theta_222logistict_2_1, m=1, regime_pars=c(phi10_112, A11_112, A11_122, vech(Omega1_112))),
               c(phi10_112, phi20_222, A11_112, A11_122, A21_222, A22_222, vech(Omega1_112), vech(Omega2_222), c_and_gamma_222_2_1, df_222_2_1))
  expect_equal(change_regime(p=1, M=2, d=3, params=theta_123expt_1_1, m=1, regime_pars=c(1:3, 3:11, 12:17)),
               c(1:3, phi20_123, 3:11, A21_123, 12:17, vech(Omega2_123), c_and_gamma_123_1_1, df_123_1_1))
  expect_equal(change_regime(p=1, M=2, d=3, params=theta_123expt_1_1, m=2, regime_pars=c(1:3, 3:11, 12:17)),
               c(phi10_123, 1:3, A11_123, 3:11, vech(Omega1_123), 12:17, c_and_gamma_123_1_1, df_123_1_1))
  expect_equal(change_regime(p=1, M=2, d=2, params=theta_122logt_1_1, m=1, regime_pars=c(phi10_112, A11_112, vech(Omega1_112))),
               c(phi10_112, phi20_122, A11_112, A21_122, vech(Omega1_112), vech(Omega2_122), gamma1_122_1_1, df_122_1_1))
  expect_equal(change_regime(p=2, M=3, d=2, params=theta_232threst_1_1, m=2, regime_pars=c(phi10_232, A11_232, A12_232, vech(Omega1_232))),
               c(phi10_232, phi10_232, phi30_232, vec(A11_232), vec(A12_232), vec(A11_232), vec(A12_232),
                 vec(A31_232), vec(A32_232), vech(Omega1_232), vech(Omega1_232), vech(Omega3_232),
                 r1_232_1_1, r2_232_1_1, df_232_1_1))

  # ind_Student
  expect_equal(change_regime(p=2, M=2, d=2, params=theta_222logistit_2_1, m=1, regime_pars=c(1:2, 3:10, 11:14), cond_dist="ind_Student"),
               c(1:2, phi20_222, 3:10, A21_222, A22_222, 11:14, vec(B2_222), c_and_gamma_222_2_1, dfs_222_2_1))
  expect_equal(change_regime(p=1, M=2, d=2, params=theta_122logit_1_1, m=1, regime_pars=c(phi10_112, A11_112, vec(B1_112)),
                             cond_dist="ind_Student"),
               c(phi10_112, phi20_122, A11_112, A21_122, vec(B1_112), vec(B2_122), gamma1_122_1_1, dfs_122_1_1))
  expect_equal(change_regime(p=1, M=2, d=3, params=theta_123expit_1_1, m=2, regime_pars=c(1:3, 3:11, 12:20), cond_dist="ind_Student"),
               c(phi10_123, 1:3, A11_123, 3:11, vec(B1_123), 12:20, c_and_gamma_123_1_1, dfs_123_1_1))
  expect_equal(change_regime(p=1, M=3, d=2, params=theta_132thresit_1_1, m=2, regime_pars=c(phi10_112, A11_112, vec(B1_112)),
                             cond_dist="ind_Student"),
               c(phi10_132, phi10_112, phi30_132, vec(A11_132), vec(A11_112), vec(A31_132),
                 vec(B1_132), vec(B1_112), vec(B3_132), r1_132_1_1, r2_132_1_1, dfs_132_1_1))
  expect_equal(change_regime(p=1, M=3, d=2, params=theta_132thresit_1_1, m=1, regime_pars=c(phi10_112, A11_112, vec(B1_112)),
                             cond_dist="ind_Student"),
               c(phi10_112, phi20_132, phi30_132, vec(A11_112), vec(A21_132), vec(A31_132),
                 vec(B1_112), vec(B2_132), vec(B3_132), r1_132_1_1, r2_132_1_1, dfs_132_1_1))
  expect_equal(change_regime(p=1, M=3, d=2, params=theta_132thresit_1_1, m=3, regime_pars=c(phi10_112, A11_112, vec(B1_112)),
                             cond_dist="ind_Student"),
               c(phi10_132, phi20_132, phi10_112, vec(A11_132), vec(A21_132), vec(A11_112),
                 vec(B1_132), vec(B2_132), vec(B1_112), r1_132_1_1, r2_132_1_1, dfs_132_1_1))
})


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


calc_mu <- function(p, M, d, params, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                    cond_dist = c("Gaussian", "Student", "ind_Student"),
                    identification = c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity"),
                    AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL,
                    weightfun_pars=NULL) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  identification <- match.arg(identification)
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, weight_function=weight_function, cond_dist=cond_dist,
                                    identification=identification, AR_constraints=AR_constraints,
                                    mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                    B_constraints=B_constraints, weightfun_pars=weightfun_pars)
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

theta_112relgc_mu <- change_parametrization(p=1, M=1, d=2, params=theta_112relgc, AR_constraints=C_112, change_to="mean")
theta_222relgc2_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222relgc2, AR_constraints=C_222_2, change_to="mean")
theta_132relgc_mu <- change_parametrization(p=1, M=3, d=2, params=theta_132relgc, AR_constraints=C_132, change_to="mean")
theta_123relgc_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123relgc, AR_constraints=C_123, change_to="mean")

theta_122log_1_1_mu <- change_parametrization(p=1, M=2, d=2, params=theta_122log_1_1, weight_function="mlogit",
                                              weightfun_pars=list(vars=1, lags=1), change_to="mean")
theta_222logc_12_2_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222logc_12_2, weight_function="mlogit",
                                                weightfun_pars=list(vars=1:2, lags=2), AR_constraints=C_222, change_to="mean")
theta_123logc_123_1_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123logc_123_1, weight_function="mlogit",
                                                 weightfun_pars=list(vars=1:3, lags=1), AR_constraints=C_123, change_to="mean")

theta_132relgw_mu <- change_parametrization(p=1, M=3, d=2, params=theta_132relgw, weight_function="relative_dens",
                                            weight_constraints=list(R=matrix(c(0.9, 0.5), nrow=2), r=c(0.13, 0.13)), change_to="mean")
theta_122logw_1_1_mu <- change_parametrization(p=1, M=2, d=2, params=theta_122logw_1_1, weight_function="mlogit",
                                               weightfun_pars=list(vars=1, lags=1), weight_constraints=list(R=0, r=c(0.12, 0.13)),
                                               change_to="mean")

theta_122logistic_1_1_mu <- change_parametrization(p=1, M=2, d=2, params=theta_122logistic_1_1, weight_function="logistic",
                                                   weightfun_pars=c(1, 1), change_to="mean")
theta_222logisticc_2_1_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222logisticc_2_1, weight_function="logistic",
                                                    weightfun_pars=c(2, 1), AR_constraints=C_222, change_to="mean")
theta_122logisticw_1_1_mu <- change_parametrization(p=1, M=2, d=2, params=theta_122logisticw_1_1, weight_function="logistic",
                                                    weightfun_pars=c(1, 1), weight_constraints=list(R=0, r=c(0.02, 0.13)), change_to="mean")


theta_122exp_1_1_mu <- change_parametrization(p=1, M=2, d=2, params=theta_122exp_1_1, weight_function="exponential",
                                              weightfun_pars=c(1, 1), change_to="mean")
theta_123expc_3_1_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123expc_3_1, weight_function="exponential",
                                               weightfun_pars=c(3, 1), AR_constraints=C_123, change_to="mean")

theta_122thres_1_1_mu <- change_parametrization(p=1, M=2, d=2, params=theta_122thres_1_1, weight_function="threshold",
                                                weightfun_pars=c(1, 1), change_to="mean")
theta_232thresw_1_1_mu <- change_parametrization(p=2, M=3, d=2, params=theta_232thresw_1_1, weight_function="threshold",
                                                 weightfun_pars=c(1, 1), weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.1, 0)),
                                                 change_to="mean")

theta_123exo_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123exo, weight_function="exogenous",
                                          weightfun_pars=cbind(c(0.4, 0.2, 0.9), c(0.6, 0.8, 0.1)), change_to="mean")
theta_222exoct_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222exoct, weight_function="exogenous",
                                            weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)),
                                            cond_dist="Student", AR_constraints=C_222, change_to="mean")


# Student
theta_122logt_1_1_mu <- change_parametrization(p=1, M=2, d=2, params=theta_122logt_1_1, weight_function="mlogit",
                                              weightfun_pars=list(vars=1, lags=1), cond_dist="Student", change_to="mean")
theta_232threswt_1_1_mu <- change_parametrization(p=2, M=3, d=2, params=theta_232threswt_1_1, weight_function="threshold",
                                                  weightfun_pars=c(1, 1), weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.1, 0)),
                                                  cond_dist="Student", change_to="mean")
theta_222logistict_2_1_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222logistict_2_1, weight_function="logistic",
                                                    weightfun_pars=c(2, 1), cond_dist="Student", change_to="mean")
theta_123expt_1_1_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123expt_1_1, weight_function="exponential",
                                               weightfun_pars=c(1, 1), cond_dist="Student", change_to="mean")

# Ind student
theta_132thresit_1_1_mu <- change_parametrization(p=1, M=3, d=2, params=theta_132thresit_1_1, weight_function="threshold",
                                                  weightfun_pars=c(1, 1), cond_dist="ind_Student", change_to="mean")
theta_222exoit_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222exoit, weight_function="exogenous",
                                                weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_Student",
                                                AR_constraints=C_222, change_to="mean")


## Structural
theta_122relgshc_mu <- change_parametrization(p=1, M=2, d=2, params=theta_122relgshc, weight_function="relative_dens",
                                              identification="heteroskedasticity", AR_constraints=C_122, change_to="mean")
theta_132relgshb_mu <- change_parametrization(p=1, M=3, d=2, params=theta_132relgshb,  weight_function="relative_dens",
                                              identification="heteroskedasticity", B_constraints=matrix(c(0.1, NA, 0.3, 0), nrow=2),
                                              change_to="mean")
theta_122logshwb_12_1_mu <- change_parametrization(p=1, M=2, d=2, params=theta_122logshwb_12_1, weight_function="mlogit",
                                                   weightfun_pars=list(vars=1:2, lags=1), identification="heteroskedasticity",
                                                   weight_constraints=list(R=0, r=c(0.1, 0.2, 0.3)),
                                                   B_constraints=matrix(c(0.1, 0.2, 0.3, 0), nrow=2), change_to="mean")
theta_123expshcwb_1_1_mu <- change_parametrization(p=1, M=2, d=3, params=theta_123expshcwb_1_1, weight_function="exponential",
                                                   weightfun_pars=c(1, 1), identification="heteroskedasticity",
                                                   AR_constraints=C_123, weight_constraints=list(R=matrix(c(1, 0.5), nrow=2), r=c(0, 0)),
                                                   B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19),
                                                                        nrow=3, ncol=3, byrow=FALSE), change_to="mean")
theta_232threstshb_1_1_mu <- change_parametrization(p=2, M=3, d=2, params=theta_232threstshb_1_1, weight_function="threshold",
                                                    weightfun_pars=c(1, 1), cond_dist="Student",
                                                    identification="heteroskedasticity",
                                                    B_constraints=matrix(c(0.1, 0.2, -0.3, 0), nrow=2),
                                                    change_to="mean")

theta_222logistitngb_2_1_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222logistitngb_2_1, cond_dist="ind_Student",
                                                      weight_function="logistic", weightfun_pars=c(2, 1), identification="non-Gaussianity",
                                                      B_constraints=matrix(c(NA, -1, 0, 1), nrow=2), change_to="mean")
theta_222exoitngb_2_1_mu <- change_parametrization(p=2, M=2, d=2, params=theta_222exoitngb_2_1, weight_function="exogenous",
                                                   weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_Student",
                                                   AR_constraints=C_222, identification="non-Gaussianity",
                                                   B_constraints=matrix(c(NA, NA, 0, 1), nrow=2), change_to="mean")


test_that("change_parametrization works correctly", {
  expect_equal(pick_phi0(M=3, d=2, params=theta_232thresw_1_1_mu),
               calc_mu(p=2, M=3, d=2, params=theta_232thresw_1_1, weight_function="threshold", weightfun_pars=c(1, 1)))
  expect_equal(change_parametrization(p=2, M=3, d=2, params=theta_232thresw_1_1_mu, weight_function="threshold",
                                      weightfun_pars=c(1, 1), weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.1, 0)),
                                      change_to="intercept"), theta_232thresw_1_1)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122thres_1_1_mu),
               calc_mu(p=1, M=2, d=2, params=theta_122thres_1_1, weight_function="threshold", weightfun_pars=c(1, 1)))
  expect_equal(change_parametrization(p=1, M=2, d=2, params=theta_122thres_1_1_mu, weight_function="exponential", weightfun_pars=c(1, 1),
                                      change_to="intercept"), theta_122thres_1_1)

  expect_equal(pick_phi0(M=2, d=3, params=theta_123expc_3_1_mu),
               calc_mu(p=1, M=2, d=3, params=theta_123expc_3_1, weight_function="exponential", weightfun_pars=c(3, 1), AR_constraints=C_123))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123expc_3_1_mu, weight_function="exponential", weightfun_pars=c(3, 1),
                                      AR_constraints=C_123, change_to="intercept"), theta_123expc_3_1)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122exp_1_1_mu),
               calc_mu(p=1, M=2, d=2, params=theta_122exp_1_1, weight_function="exponential", weightfun_pars=c(1, 1)))
  expect_equal(change_parametrization(p=1, M=2, d=2, params=theta_122exp_1_1_mu, weight_function="exponential", weightfun_pars=c(1, 1),
                                      change_to="intercept"), theta_122exp_1_1)

  expect_equal(pick_phi0(M=2, d=2, params=theta_122logistic_1_1_mu),
               calc_mu(p=1, M=2, d=2, params=theta_122logistic_1_1, weight_function="logistic", weightfun_pars=c(1, 1)))
  expect_equal(change_parametrization(p=1, M=2, d=2, params=theta_122logistic_1_1_mu, weight_function="logistic", weightfun_pars=c(1, 1),
                                      change_to="intercept"), theta_122logistic_1_1)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222logisticc_2_1_mu),
               calc_mu(p=2, M=2, d=2, params=theta_222logisticc_2_1, weight_function="logistic", weightfun_pars=c(2, 1), AR_constraints=C_222))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222logisticc_2_1_mu, weight_function="logistic", weightfun_pars=c(2, 1),
                                      AR_constraints=C_222, change_to="intercept"), theta_222logisticc_2_1)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122logisticw_1_1_mu),
               calc_mu(p=1, M=2, d=2, params=theta_122logisticw_1_1, weight_function="logistic", weightfun_pars=c(1, 1),
                       weight_constraints=list(R=0, r=c(0.02, 0.13))))
  expect_equal(change_parametrization(p=1, M=2, d=2, params=theta_122logisticw_1_1_mu, weight_function="logistic", weightfun_pars=c(1, 1),
                                      weight_constraints=list(R=0, r=c(0.02, 0.13)), change_to="intercept"), theta_122logisticw_1_1)


  expect_equal(pick_phi0(M=3, d=2, params=theta_132relgw_mu),
               calc_mu(p=1, M=3, d=2, params=theta_132relgw, weight_function="relative_dens",
                       weight_constraints=list(R=matrix(c(0.9, 0.5), nrow=2), r=c(0.13, 0.13))))
  expect_equal(change_parametrization(p=1, M=3, d=2, params=theta_132relgw_mu, weight_function="relative_dens",
                                      weight_constraints=list(R=matrix(c(0.9, 0.5), nrow=2), r=c(0.13, 0.13)),
                                      change_to="intercept"), theta_132relgw)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122logw_1_1_mu),
               calc_mu(p=1, M=2, d=2, params=theta_122logw_1_1, weight_function="mlogit", weight_constraints=list(R=0, r=c(0.12, 0.13))))
  expect_equal(change_parametrization(p=1, M=2, d=2, params=theta_122logw_1_1_mu, weight_function="mlogit",
                                      weight_constraints=list(R=0, r=c(0.12, 0.13)), change_to="intercept"), theta_122logw_1_1)

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

  expect_equal(pick_phi0(M=1, d=2, params=theta_112relgc_mu), calc_mu(p=1, M=1, d=2, params=theta_112relgc, AR_constraints=C_112))
  expect_equal(change_parametrization(p=1, M=1, d=2, params=theta_112relgc_mu, AR_constraints=C_112, change_to="intercept"),
               theta_112relgc)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222relgc2_mu), calc_mu(p=2, M=2, d=2, params=theta_222relgc2, AR_constraints=C_222_2))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222relgc2_mu, AR_constraints=C_222_2, change_to="intercept"),
               theta_222relgc2)
  expect_equal(pick_phi0(M=3, d=2, params=theta_132relgc_mu), calc_mu(p=1, M=3, d=2, params=theta_132relgc, AR_constraints=C_132))
  expect_equal(change_parametrization(p=1, M=3, d=2, params=theta_132relgc_mu, AR_constraints=C_132, change_to="intercept"),
               theta_132relgc)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123relgc_mu), calc_mu(p=1, M=2, d=3, params=theta_123relgc, AR_constraints=C_123))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123relgc_mu, AR_constraints=C_123, change_to="intercept"),
               theta_123relgc)

  expect_equal(pick_phi0(M=2, d=2, params=theta_122log_1_1_mu),
               calc_mu(p=1, M=2, d=2, params=theta_122log_1_1, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1)))
  expect_equal(change_parametrization(p=1, M=2, d=2, params=theta_122log_1_1_mu, weight_function="mlogit",
                                      weightfun_pars=list(vars=1, lags=1), change_to="intercept"), theta_122log_1_1)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222logc_12_2_mu),
               calc_mu(p=2, M=2, d=2, params=theta_222logc_12_2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2),
                       AR_constraints=C_222))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222logc_12_2_mu, weight_function="mlogit",
                                      weightfun_pars=list(vars=1:2, lags=2), AR_constraints=C_222,
                                      change_to="intercept"), theta_222logc_12_2)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123logc_123_1_mu),
               calc_mu(p=1, M=2, d=3, params=theta_123logc_123_1, weight_function="mlogit", weightfun_pars=list(vars=1:3, lags=1),
                       AR_constraints=C_123))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123logc_123_1_mu, weight_function="mlogit",
                                      weightfun_pars=list(vars=1:3, lags=1), AR_constraints=C_123,
                                      change_to="intercept"), theta_123logc_123_1)

  expect_equal(pick_phi0(M=2, d=3, params=theta_123exo_mu),
               calc_mu(p=1, M=2, d=3, params=theta_123exo, weight_function="exogenous",
                       weightfun_pars=cbind(c(0.4, 0.2, 0.9), c(0.6, 0.8, 0.1))))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123exo_mu, weight_function="exogenous",
                                      weightfun_pars=cbind(c(0.4, 0.2, 0.9), c(0.6, 0.8, 0.1)),
                                      change_to="intercept"), theta_123exo)


  # Student
  expect_equal(pick_phi0(M=3, d=2, params=theta_232threswt_1_1_mu),
               calc_mu(p=2, M=3, d=2, params=theta_232threswt_1_1, weight_function="threshold", weightfun_pars=c(1, 1)),
                       cond_dist="Student")
  expect_equal(change_parametrization(p=2, M=3, d=2, params=theta_232threswt_1_1_mu, weight_function="threshold",
                                      weightfun_pars=c(1, 1), weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.1, 0)),
                                      cond_dist="Student", change_to="intercept"), theta_232threswt_1_1)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122logt_1_1_mu),
               calc_mu(p=1, M=2, d=2, params=theta_122logt_1_1,  weight_function="mlogit",
                       weightfun_pars=list(vars=1, lags=1), cond_dist="Student"))
  expect_equal(change_parametrization(p=1, M=2, d=2, params=theta_122logt_1_1_mu,  weight_function="mlogit",
                                      weightfun_pars=list(vars=1, lags=1), cond_dist="Student",
                                      change_to="intercept"), theta_122logt_1_1)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222logistict_2_1_mu),
               calc_mu(p=2, M=2, d=2, params=theta_222logistict_2_1, weight_function="logistic",
                       weightfun_pars=c(2, 1), cond_dist="Student"))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222logistict_2_1_mu, weight_function="logistic",
                                      weightfun_pars=c(2, 1), cond_dist="Student",
                                      change_to="intercept"), theta_222logistict_2_1)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123expt_1_1_mu),
               calc_mu(p=1, M=2, d=3, params=theta_123expt_1_1, weight_function="exponential",
                       weightfun_pars=c(1, 1), cond_dist="Student"))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123expt_1_1_mu, weight_function="exponential",
                                      weightfun_pars=c(1, 1), cond_dist="Student",
                                      change_to="intercept"), theta_123expt_1_1)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222exoct_mu),
               calc_mu(p=2, M=2, d=2, params=theta_222exoct, weight_function="exogenous",
                       weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)),
                       cond_dist="Student", AR_constraints=C_222))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222exoct_mu,weight_function="exogenous",
                                      weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="Student",
                                      AR_constraints=C_222, change_to="intercept"), theta_222exoct)

  # Ind student
  expect_equal(pick_phi0(M=3, d=2, params=theta_132thresit_1_1_mu),
               calc_mu(p=1, M=3, d=2, params=theta_132thresit_1_1, weight_function="threshold",
                       weightfun_pars=c(1, 1), cond_dist="ind_Student"))
  expect_equal(change_parametrization(p=1, M=3, d=2, params=theta_132thresit_1_1_mu, weight_function="threshold",
                                      weightfun_pars=c(1, 1), cond_dist="ind_Student", change_to="intercept"), theta_132thresit_1_1)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222exoit_mu),
               calc_mu(p=2, M=2, d=2, params=theta_222exoit, weight_function="exogenous",
                       weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_Student", AR_constraints=C_222))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222exoit_mu, weight_function="exogenous",
                                      weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_Student",
                                      AR_constraints=C_222), theta_222exoit)

  # Structural
  expect_equal(pick_phi0(M=2, d=2, params=theta_122relgshc_mu),
               calc_mu(p=1, M=2, d=2, params=theta_122relgshc, weight_function="relative_dens",
                       identification="heteroskedasticity", AR_constraints=C_122))
  expect_equal(change_parametrization(p=1, M=2, d=2, params=theta_122relgshc_mu, weight_function="relative_dens",
                                      identification="heteroskedasticity", AR_constraints=C_122,
                                      change_to="intercept"), theta_122relgshc)
  expect_equal(pick_phi0(M=3, d=2, params=theta_132relgshb_mu),
               calc_mu(p=1, M=3, d=2, params=theta_132relgshb, weight_function="relative_dens",
                       identification="heteroskedasticity", B_constraints=matrix(c(0.1, NA, 0.3, 0), nrow=2)))
  expect_equal(change_parametrization(p=1, M=3, d=2, params=theta_132relgshb_mu,   weight_function="relative_dens",
                                      identification="heteroskedasticity", B_constraints=matrix(c(0.1, NA, 0.3, 0), nrow=2),
                                      change_to="intercept"), theta_132relgshb)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122logshwb_12_1_mu),
               calc_mu(p=1, M=2, d=2, params=theta_122logshwb_12_1, weight_function="mlogit",
                       weightfun_pars=list(vars=1:2, lags=1), identification="heteroskedasticity",
                       weight_constraints=list(R=0, r=c(0.1, 0.2, 0.3)),
                       B_constraints=matrix(c(0.1, 0.2, 0.3, 0), nrow=2)))
  expect_equal(change_parametrization(p=1, M=2, d=2, params=theta_122logshwb_12_1_mu, weight_function="mlogit",
                                      weightfun_pars=list(vars=1:2, lags=1), identification="heteroskedasticity",
                                      weight_constraints=list(R=0, r=c(0.1, 0.2, 0.3)),
                                      B_constraints=matrix(c(0.1, 0.2, 0.3, 0), nrow=2),
                                      change_to="intercept"), theta_122logshwb_12_1)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123expshcwb_1_1_mu),
               calc_mu(p=1, M=2, d=3, params=theta_123expshcwb_1_1, weight_function="exponential",
                       weightfun_pars=c(1, 1), identification="heteroskedasticity",
                       AR_constraints=C_123, weight_constraints=list(R=matrix(c(1, 0.5), nrow=2), r=c(0, 0)),
                       B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE)))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123expshcwb_1_1_mu, weight_function="exponential",
                                      weightfun_pars=c(1, 1), identification="heteroskedasticity",
                                      AR_constraints=C_123, weight_constraints=list(R=matrix(c(1, 0.5), nrow=2), r=c(0, 0)),
                                      B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19),
                                                           nrow=3, ncol=3, byrow=FALSE),
                                      change_to="intercept"), theta_123expshcwb_1_1)
  expect_equal(pick_phi0(M=3, d=2, params=theta_232threstshb_1_1_mu),
               calc_mu(p=2, M=3, d=2, params=theta_232threstshb_1_1, weight_function="threshold",
                       weightfun_pars=c(1, 1), cond_dist="Student", identification="heteroskedasticity",
                       B_constraints=matrix(c(0.1, 0.2, -0.3, 0), nrow=2)))
  expect_equal(change_parametrization(p=2, M=3, d=2, params=theta_232threstshb_1_1_mu, weight_function="threshold",
                                      weightfun_pars=c(1, 1), cond_dist="Student", identification="heteroskedasticity",
                                      B_constraints=matrix(c(0.1, 0.2, -0.3, 0), nrow=2), change_to="intercept"), theta_232threstshb_1_1)

  expect_equal(pick_phi0(M=2, d=2, params=theta_222logistitngb_2_1_mu),
               calc_mu(p=2, M=2, d=2, params=theta_222logistitngb_2_1, cond_dist="ind_Student", weight_function="logistic",
                       weightfun_pars=c(2, 1), identification="non-Gaussianity", B_constraints=matrix(c(NA, -1, 0, 1), nrow=2)))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222logistitngb_2_1_mu, cond_dist="ind_Student",
                                      weight_function="logistic", weightfun_pars=c(2, 1), identification="non-Gaussianity",
                                      B_constraints=matrix(c(NA, -1, 0, 1), nrow=2), change_to="intercept"), theta_222logistitngb_2_1)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222exoitngb_2_1_mu),
               calc_mu(p=2, M=2, d=2, params=theta_222exoitngb_2_1, weight_function="exogenous",
                       weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_Student",
                       AR_constraints=C_222, identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222exoitngb_2_1_mu, weight_function="exogenous",
                                      weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_Student",
                                      AR_constraints=C_222, identification="non-Gaussianity",
                                      B_constraints=matrix(c(NA, NA, 0, 1), nrow=2), change_to="intercept"), theta_222exoitngb_2_1)

  ##### Test change_parametrization to change between "orig" and "alt" (B_m and B_m* parametrization):
  expect_equal(change_parametrization(p=1, M=1, d=2, params=theta_112it, cond_dist="ind_Student", weight_function="threshold",
                                      weightfun_pars=c(1, 1), change_to="alt"), theta_112it)
  expect_equal(change_parametrization(p=1, M=1, d=2, params=theta_112it, cond_dist="ind_Student", weight_function="threshold",
                                      weightfun_pars=c(1, 1), change_to="orig"), theta_112it)
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222logistit_2_1, cond_dist="ind_Student", weight_function="logistic",
                                      weightfun_pars=c(2, 1), change_to="alt"),
               c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222), vec(B1_222), vec(B2_222)-vec(B1_222),
                 c_and_gamma_222_2_1, dfs_222_2_1))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222logistit_2_1, cond_dist="ind_Student", weight_function="logistic",
                                      weightfun_pars=c(2, 1), change_to="orig"),
               c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222), vec(B1_222), vec(B2_222)+vec(B1_222),
                 c_and_gamma_222_2_1, dfs_222_2_1))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123expit_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                                      cond_dist="ind_Student", change_to="alt"),
               c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(B1_123), vec(B2_123)-vec(B1_123), c_and_gamma_123_1_1, dfs_123_1_1))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123expit_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                                      cond_dist="ind_Student", change_to="orig"),
               c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(B1_123), vec(B2_123)+vec(B1_123), c_and_gamma_123_1_1, dfs_123_1_1))
  expect_equal(change_parametrization(p=1, M=3, d=2, params=theta_132thresit_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                      cond_dist="ind_Student", change_to="alt"),
               c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                 vec(B1_132), vec(B2_132)-vec(B1_132), vec(B3_132)-vec(B1_132), r1_132_1_1, r2_132_1_1, dfs_132_1_1))
  expect_equal(change_parametrization(p=1, M=3, d=2, params=theta_132thresit_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                      cond_dist="ind_Student", change_to="orig"),
               c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                 vec(B1_132), vec(B2_132)+vec(B1_132), vec(B3_132)+vec(B1_132), r1_132_1_1, r2_132_1_1, dfs_132_1_1))
  expect_equal(change_parametrization(p=1, M=3, d=2, params=theta_132thresmwit_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                      cond_dist="ind_Student", mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2)),
                                      change_to="alt"),
               c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(B1_132), vec(B2_132)- vec(B1_132),
                 vec(B3_132)-vec(B1_132), dfs_132_1_1))
  expect_equal(change_parametrization(p=1, M=3, d=2, params=theta_132thresmwit_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                      cond_dist="ind_Student", mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2)),
                                      change_to="orig"),
               c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(B1_132), vec(B2_132)+vec(B1_132),
                 vec(B3_132)+vec(B1_132), dfs_132_1_1))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123logisticcmit_3_1, weight_function="logistic", weightfun_pars=c(3, 1),
                                      cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_123, change_to="alt"),
               c(phi10_123, vec(A11_123), vec(B1_123), vec(B2_123)-vec(B1_123), c_and_gamma_123_3_1, dfs_123_3_1))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123logisticcmit_3_1, weight_function="logistic", weightfun_pars=c(3, 1),
                                      cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_123, change_to="orig"),
               c(phi10_123, vec(A11_123), vec(B1_123), vec(B2_123)+vec(B1_123), c_and_gamma_123_3_1, dfs_123_3_1))
  # B_constraints
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222logistitngb_2_1, cond_dist="ind_Student", weight_function="logistic",
                                      weightfun_pars=c(2, 1), identification="non-Gaussianity",
                                      B_constraints=matrix(c(NA, -1, 0, 1), nrow=2), change_to="alt"),
               c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                 Wvec(B1_222c), Wvec(B2_222c)-Wvec(B1_222c), c_and_gamma_222_2_1, dfs_222_2_1))
  expect_equal(change_parametrization(p=2, M=2, d=2, params=theta_222logistitngb_2_1, cond_dist="ind_Student", weight_function="logistic",
                                      weightfun_pars=c(2, 1), identification="non-Gaussianity",
                                      B_constraints=matrix(c(NA, -1, 0, 1), nrow=2), change_to="orig"),
               c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                 Wvec(B1_222c), Wvec(B2_222c)+Wvec(B1_222c), c_and_gamma_222_2_1, dfs_222_2_1))
  expect_equal(change_parametrization(p=1, M=2, d=2, params=theta_122logitngb_1_1, weight_function="mlogit",
                                      weightfun_pars=list(vars=1, lags=1), cond_dist="ind_Student", identification="non-Gaussianity",
                                      B_constraints=matrix(c(1, NA, -1, 1), nrow=2), change_to="alt"),
               c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), Wvec(B1_122c), Wvec(B2_122c)-Wvec(B1_122c), gamma1_122_1_1, dfs_122_1_1))
  expect_equal(change_parametrization(p=1, M=2, d=2, params=theta_122logitngb_1_1, weight_function="mlogit",
                                      weightfun_pars=list(vars=1, lags=1), cond_dist="ind_Student", identification="non-Gaussianity",
                                      B_constraints=matrix(c(1, NA, -1, 1), nrow=2), change_to="orig"),
               c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), Wvec(B1_122c), Wvec(B2_122c)+Wvec(B1_122c), gamma1_122_1_1, dfs_122_1_1))
  expect_equal(change_parametrization(p=1, M=3, d=2, params=theta_132thresitngb_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                      cond_dist="ind_Student", identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2),
                                      change_to="alt"),
               c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132), Wvec(B1_132c),
                 Wvec(B2_132c)-Wvec(B1_132c), Wvec(B3_132c)-Wvec(B1_132c), r1_132_1_1, r2_132_1_1, dfs_132_1_1))
  expect_equal(change_parametrization(p=1, M=3, d=2, params=theta_132thresitngb_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                      cond_dist="ind_Student", identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2),
                                      change_to="orig"),
               c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132), Wvec(B1_132c),
                 Wvec(B2_132c)+Wvec(B1_132c), Wvec(B3_132c)+Wvec(B1_132c), r1_132_1_1, r2_132_1_1, dfs_132_1_1))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123logisticcmitngb_3_1, weight_function="logistic", weightfun_pars=c(3, 1),
                                      cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_123,
                                      identification="non-Gaussianity",
                                      B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3), change_to="alt"),
               c(phi10_123, vec(A11_123), Wvec(B1_123c), Wvec(B2_123c)-Wvec(B1_123c), c_and_gamma_123_3_1, dfs_123_3_1))
  expect_equal(change_parametrization(p=1, M=2, d=3, params=theta_123logisticcmitngb_3_1, weight_function="logistic", weightfun_pars=c(3, 1),
                                      cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_123,
                                      identification="non-Gaussianity",
                                      B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3), change_to="orig"),
               c(phi10_123, vec(A11_123), Wvec(B1_123c), Wvec(B2_123c)+Wvec(B1_123c), c_and_gamma_123_3_1, dfs_123_3_1))
})



# p=1, M=3, d=2, weight_function="relative_dens", identification="heteroskedasticity", to sort
alpha1_132_22 <- 0.2; alpha2_132_22 <- 0.3
theta_132relgsh_2 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                       vec(W_132), lambdas2_132, lambdas3_132, alpha1_132_22, alpha2_132_22)
sortedWpars_132 <- redecompose_Omegas(M=3, d=2, W=W_132, lambdas=c(lambdas2_132, lambdas3_132), perm=c(3, 2, 1))
theta_132relgsh_2_sorted <- c(phi30_132, phi20_132, phi10_132, vec(A31_132), vec(A21_132), vec(A11_132),
                              sortedWpars_132, 1 - alpha1_132_22 - alpha2_132_22, alpha2_132_22)

# p=1, M=3, d=2, weight_function="relative_dens", identification="heteroskedasticity", to sort
alpha1_132_33 <- 0.5; alpha2_132_33 <- 0.2
theta_132relgsh_3 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                       vec(W_132), lambdas2_132, lambdas3_132, alpha1_132_33, alpha2_132_33)
sortedWpars_132_3 <- redecompose_Omegas(M=3, d=2, W=W_132, lambdas=c(lambdas2_132, lambdas3_132), perm=c(1, 3, 2))
theta_132relgsh_3_sorted <- c(phi10_132, phi30_132, phi20_132, vec(A11_132), vec(A31_132), vec(A21_132),
                              sortedWpars_132_3, alpha1_132_33, 1 - alpha1_132_33 - alpha2_132_33)

# p=1, M=3, d=2, weight_function="relative_dens", identification="heteroskedasticity",
# B_constraints=matrix(c(0.1, NA, 0.3, 0), nrow=2); to sort!
W_132b <- matrix(c(0.11, 0.22, 0.33, 0), nrow=2);
theta_132relgshb_2 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                        Wvec(W_132b), lambdas2_132, lambdas3_132, alpha1_132_22, alpha2_132_22)
sortedWpars_132b_2 <- redecompose_Omegas(M=3, d=2, W=W_132b, lambdas=c(lambdas2_132, lambdas3_132), perm=c(3, 2, 1))
theta_132relgshb_2_sorted <- c(phi30_132, phi20_132, phi10_132, vec(A31_132), vec(A21_132), vec(A11_132),
                               Wvec(sortedWpars_132b_2), 1 - alpha1_132_22 - alpha2_132_22, alpha2_132_22)

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

  # Struct models
  expect_equal(sort_regimes(p=1, M=2, d=2, params=theta_122relgsh, weight_function="relative_dens", identification="heteroskedasticity"),
               theta_122relgsh)
  expect_equal(sort_regimes(p=1, M=3, d=2, params=theta_132relgsh_2, weight_function="relative_dens", identification="heteroskedasticity"),
               theta_132relgsh_2_sorted)
  expect_equal(sort_regimes(p=1, M=3, d=2, params=theta_132relgsh_3, weight_function="relative_dens", identification="heteroskedasticity"),
               theta_132relgsh_3_sorted)
  expect_equal(sort_regimes(p=1, M=3, d=2, params=theta_132relgshb, weight_function="relative_dens", identification="heteroskedasticity",
                            B_constraints=matrix(c(0.1, NA, 0.3, 0), nrow=2)), theta_132relgshb)
  expect_equal(sort_regimes(p=1, M=3, d=2, params=theta_132relgshb_2, weight_function="relative_dens", identification="heteroskedasticity",
                            B_constraints=matrix(c(0.1, NA, 0.3, 0), nrow=2)), theta_132relgshb_2_sorted)

  # Does not sort:
  expect_equal(sort_regimes(p=1, M=2, d=2, params=theta_122log_1_1, weight_function="mlogit"), theta_122log_1_1)
  expect_equal(sort_regimes(p=1, M=2, d=2, params=theta_122logistic_1_1, weight_function="logistic"), theta_122logistic_1_1)
  expect_equal(sort_regimes(p=1, M=2, d=2, params=theta_122exp_1_1, weight_function="exponential"), theta_122exp_1_1)
  expect_equal(sort_regimes(p=1, M=2, d=2, params=theta_122thres_1_1, weight_function="threshold"), theta_122thres_1_1)
  expect_equal(sort_regimes(p=1, M=2, d=2, params=theta_122logt_1_1, weight_function="mlogit", cond_dist="Student"), theta_122logt_1_1)
  expect_equal(sort_regimes(p=2, M=2, d=2, params=theta_222logistict_2_1, weight_function="logistic", cond_dist="Student"),
               theta_222logistict_2_1)
  expect_equal(sort_regimes(p=1, M=2, d=3, params=theta_123expt_1_1, weight_function="exponential", cond_dist="Student"), theta_123expt_1_1)
  expect_equal(sort_regimes(p=2, M=3, d=2, params=theta_232threst_1_1, weight_function="threshold", cond_dist="Student"), theta_232threst_1_1)
  expect_equal(sort_regimes(p=2, M=2, d=2, params=theta_222logistictsh_2_1, weight_function="logistic", weightfun_pars=c(2, 1),
                            cond_dist="Student", identification="heteroskedasticity"), theta_222logistictsh_2_1)
  expect_equal(sort_regimes(p=1, M=2, d=2, params=theta_122logsh_12_1, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=1),
                            identification="heteroskedasticity"), theta_122logsh_12_1)
  expect_equal(sort_regimes(p=1, M=2, d=3, params=theta_123expsh_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                            identification="heteroskedasticity"), theta_123expsh_1_1)
  expect_equal(sort_regimes(p=2, M=3, d=2, params=theta_232threstsh_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                            cond_dist="Student", identification="heteroskedasticity"), theta_232threstsh_1_1)
  expect_equal(sort_regimes(p=1, M=2, d=3, params=theta_123exo, weight_function="exogenous",
                            weightfun_pars=cbind(c(0.4, 0.2, 0.9), c(0.6, 0.8, 0.1))), theta_123exo)

})

B1_112_sorted <- matrix(c(0.7, -0.3, 0.5, 0.2), nrow=2)
B1_222_sorted <- matrix(c(0.5, 0.2, 0.1, -0.3), nrow=2)
B2_222_sorted <- matrix(c(0.4, -0.1, 0.2, -0.3), nrow=2)
B1_123_sorted <- matrix(c(1.0, 0.3, 0.1, 0.8, -1.1, 0.5, 0.1, 0.2, -0.4), nrow=3)
B2_123_sorted <- matrix(c(0.3, -0.2, -0.7, 0.8, -1.2, -0.5, -0.1, 0.2, -1.1), nrow=3)
B1_123_2 <- matrix(c(1.0, 0.3, 0.1, -0.8, 1.1, -0.5, -1.1, -0.2, 0.4), nrow=3)
B2_123_2 <- matrix(c(0.3, -0.2, -0.7, -0.8, 1.2, 0.5, 0.1, -0.2, 1.1), nrow=3)
theta_123expit_1_1_2 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(B1_123_2),
                          vec(B2_123_2), c_and_gamma_123_1_1, dfs_123_1_1)
B1_123_2_sorted <- matrix(c(1.1, 0.2, -0.4, 1.0, 0.3, 0.1, 0.8, -1.1, 0.5), nrow=3)
B2_123_2_sorted <- matrix(c(-0.1, 0.2, -1.1, 0.3, -0.2, -0.7, 0.8, -1.2, -0.5), nrow=3)

B1_132_sorted <- matrix(c(0.6, 0.2, 0.1, -0.7), nrow=2)
B2_132_sorted <- matrix(c(0.4, -0.1, 0.2, -0.5), nrow=2)
B3_132_sorted <- matrix(c(0.9, -0.5, -0.2, -0.4), nrow=2)

B1_132_2 <- matrix(c(-0.6, 0.2, -0.8, -0.7), nrow=2)
B2_132_2 <- matrix(c(0.4, -0.1, -0.2, 0.5), nrow=2)
B3_132_2 <- matrix(c(0.9, -0.5, 0.2, 0.4), nrow=2)
theta_132thresit_1_1_2 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                            vec(B1_132_2), vec(B2_132_2), vec(B3_132_2), r1_132_1_1, r2_132_1_1, dfs_132_1_1)

B1_132_2_sorted <- matrix(c(0.8, 0.7, 0.6, -0.2), nrow=2)
B2_132_2_sorted <- matrix(c(0.2, -0.5, -0.4, 0.1), nrow=2)
B3_132_2_sorted <- matrix(c(-0.2, -0.4, -0.9, 0.5), nrow=2)

theta_132thresmwit_1_1_2 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(B1_132_2),
                              vec(B2_132_2), vec(B3_132_2), dfs_132_1_1)
theta_123logisticcmit_3_1_2 <- c(phi10_123, vec(A11_123), vec(B1_123_2), vec(B2_123_2), c_and_gamma_123_3_1, dfs_123_3_1)


test_that("sort_impactmats works correctly", {
  expect_equal(sort_impactmats(p=1, M=1, d=2, params=theta_112it, cond_dist="ind_Student", weight_function="threshold",
                               weightfun_pars=c(1, 1)), c(phi10_112, vec(A11_112), vec(B1_112_sorted), dfs_112[2:1]))
  expect_equal(sort_impactmats(p=2, M=2, d=2, params=theta_222logistit_2_1, cond_dist="ind_Student", weight_function="logistic",
                               weightfun_pars=c(2, 1)), c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                                                          vec(B1_222_sorted), vec(B2_222_sorted), c_and_gamma_222_2_1, dfs_222_2_1))
  expect_equal(sort_impactmats(p=1, M=2, d=2, params=theta_122logit_1_1, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1),
                               cond_dist="ind_Student"), theta_122logit_1_1)
  expect_equal(sort_impactmats(p=1, M=2, d=3, params=theta_123expit_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                               cond_dist="ind_Student"), c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(B1_123_sorted),
                                                           vec(B2_123_sorted), c_and_gamma_123_1_1, dfs_123_1_1))
  expect_equal(sort_impactmats(p=1, M=2, d=3, params=theta_123expit_1_1_2, weight_function="exponential", weightfun_pars=c(1, 1),
                               cond_dist="ind_Student"), c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(B1_123_2_sorted),
                                                           vec(B2_123_2_sorted), c_and_gamma_123_1_1, dfs_123_1_1[c(3, 1, 2)]))
  expect_equal(sort_impactmats(p=1, M=3, d=2, params=theta_132thresit_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                               cond_dist="ind_Student"), c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                                                           vec(B1_132_sorted), vec(B2_132_sorted), vec(B3_132_sorted), r1_132_1_1,
                                                           r2_132_1_1, dfs_132_1_1))
  expect_equal(sort_impactmats(p=1, M=3, d=2, params=theta_132thresit_1_1_2, weight_function="threshold", weightfun_pars=c(1, 1),
                               cond_dist="ind_Student"), c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                                                           vec(B1_132_2_sorted), vec(B2_132_2_sorted), vec(B3_132_2_sorted), r1_132_1_1,
                                                           r2_132_1_1, dfs_132_1_1[2:1]))
  expect_equal(sort_impactmats(p=2, M=2, d=2, params=theta_222expcmwit_2_1, weight_function="exponential", weightfun_pars=c(2, 1),
                               cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_222,
                               weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))),
               c(phi10_222, vec(A11_222), vec(A12_222), vec(B1_222_sorted), vec(B2_222_sorted), xi_222expcmw_2_1, dfs_122_2_1))
  expect_equal(sort_impactmats(p=1, M=3, d=2, params=theta_132thresmwit_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                               cond_dist="ind_Student", mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))),
               c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(B1_132_sorted), vec(B2_132_sorted),
                 vec(B3_132_sorted), dfs_132_1_1))
  expect_equal(sort_impactmats(p=1, M=3, d=2, params=theta_132thresmwit_1_1_2, weight_function="threshold", weightfun_pars=c(1, 1),
                               cond_dist="ind_Student", mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))),
               c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(B1_132_2_sorted), vec(B2_132_2_sorted),
                 vec(B3_132_2_sorted), dfs_132_1_1[2:1]))
  expect_equal(sort_impactmats(p=1, M=2, d=3, params=theta_123logisticcmit_3_1, weight_function="logistic", weightfun_pars=c(3, 1),
                               cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_123),
               c(phi10_123, vec(A11_123), vec(B1_123_sorted), vec(B2_123_sorted), c_and_gamma_123_3_1, dfs_123_3_1))
  expect_equal(sort_impactmats(p=1, M=2, d=3, params=theta_123logisticcmit_3_1_2, weight_function="logistic", weightfun_pars=c(3, 1),
                               cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_123),
               c(phi10_123, vec(A11_123), vec(B1_123_2_sorted), vec(B2_123_2_sorted), c_and_gamma_123_3_1, dfs_123_3_1[c(3, 1, 2)]))
  expect_equal(sort_impactmats(p=2, M=2, d=2, params=theta_222logcit_2_1, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1),
                               cond_dist="ind_Student", AR_constraints=C_222),
               c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222_sorted), vec(B2_222_sorted), gamma1_222_2_1, dfs_222_2_1))
  expect_equal(sort_impactmats(p=2, M=2, d=2, params=theta_222exoit, weight_function="exogenous",
                               weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_Student", AR_constraints=C_222),
               c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222_sorted), vec(B2_222_sorted), dfs_222_2_1))

  expect_equal(sort_impactmats(p=1, M=2, d=3, params=theta_123relg_2, weight_function="relative_dens", cond_dist="Gaussian"),
               theta_123relg_2)
  expect_equal(sort_impactmats(p=2, M=2, d=2, params=theta_222logistict_2_1, cond_dist="Student", weight_function="logistic",
                               weightfun_pars=c(2, 1)), theta_222logistict_2_1)

})


