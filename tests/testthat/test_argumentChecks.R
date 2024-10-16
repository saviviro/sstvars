context("argumentChecks")
library(sstvars)


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

# p=2, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="Student"
df_232_1_1 <- 30
theta_232threst_1_1 <- c(theta_232thres_1_1, df_232_1_1)

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
C_132 <- rbind_diags(p=1, M=3, d=2)
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



### Models with mean_constraints

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

# p=1, M=2, d=2, weightfun_pars=list(vars=1, lags=1), mean_constraints=list(1:2)
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


## weight_function == "threshold"

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), mean_constraints=list(1, 2:3)
theta_132thresm_1_1 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                         vech(Omega2_132), vech(Omega3_132), r1_132_1_1, r2_132_1_1)
theta_132thresm_1_1_expanded <- c(phi10_132, phi20_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                                  vech(Omega2_132), vech(Omega3_132), r1_132_1_1, r2_132_1_1)

# p=1, M=2, p=3, weight_function="threshold", weightfun_pars=c(2, 1), AR_constraints=C_123, mean_connstraints=list(1:2)
r1_123_2_1 <- 0.5
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
                             vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), matrix(c(0.9, 0.5), nrow=2)%*%xi_132relgw + 0.13)

theta_132relgw_notpd  <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132_stab),
                           vech(Omega1_132), vech(Omega2_132), vech(Omega3_132_notpd), xi_132relgw)


# p=2, M=2, d=2, weight_function="relative_dens", mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=0, r=0.6)
theta_222relgcmw <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222)) # No weight param since replaced with r
theta_222relgcmw_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                               vech(Omega1_222), vech(Omega2_222), 0.6)


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

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))
theta_132thresmw_1_1 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                          vech(Omega2_132), vech(Omega3_132))
theta_132thresmw_1_1_expanded <- c(phi10_132, phi20_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                                   vech(Omega2_132), vech(Omega3_132), 0, 1.2)


### Student

# p=2, M=2, d=2, cond_dist="Student", weight_function="logistic", weightfun_pars=c(2, 1)
df_222_2_1 <- 3
theta_222logistict_2_1 <- c(theta_222logistic_2_1, df_222_2_1)

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="Student"
df_122_1_1 <- 13
theta_122logt_1_1 <- c(theta_122log_1_1, df_122_1_1)

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="Student"
df_123_1_1 <- 10
theta_123expt_1_1 <- c(theta_123exp_1_1, df_123_1_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="Student"
df_132_1_1 <- 30
theta_132threst_1_1 <- c(theta_132thres_1_1, df_132_1_1)


# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="Student",
# mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))
df_122_2_1 <- 4
theta_222expcmwt_2_1 <- c(theta_222expcmw_2_1, df_122_2_1)
theta_222expcmwt_2_1_expanded <- c(theta_222expcmw_2_1_expanded, df_122_2_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="Student",
# mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))
theta_132thresmwt_1_1 <- c(theta_132thresmw_1_1, df_132_1_1)
theta_132thresmwt_1_1_expanded <- c(theta_132thresmw_1_1_expanded, df_132_1_1)

# p=1, M=2, p=3, weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="Student",
# mean_constraints=list(1:2), AR_constraints=C_123
df_123_3_1 <- 11
theta_123logisticcmt_3_1 <- c(theta_123logisticcm_3_1, df_123_3_1)
theta_123logisticcmt_3_1_expanded <- c(theta_123logisticcm_3_1_expanded, df_123_3_1)

# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="Student", AR_constraints=C_222
theta_222logct_2_1 <- c(theta_222logc_2_1, df_222_2_1)
theta_222logct_expanded <- c(theta_222logc_2_1_expanded, df_222_2_1)

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="Student", AR_constraints=C_222
theta_222exot_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), df_222_2_1)
theta_222exot_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
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
theta_222exoit_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222), dfs_222_2_1)
theta_222exoit_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                             vec(B1_222), vec(B2_222), dfs_222_2_1)


################
### ind_skewed_t

# p=1, M=1, d=2, cond_dist="ind_skewed_t", weight_function="threshold", weightfun_pars=c(1, 1)
dfls_112 <- c(3, 7, 0.1, -0.2)
B1_112 <- matrix(c(0.5, 0.2, -0.7, 0.3), nrow=2)
theta_112ikt <- c(phi10_112, vec(A11_112), vec(B1_112), dfls_112)

# p=2, M=2, d=2, cond_dist="ind_skewed_t", weight_function="logistic", weightfun_pars=c(2, 1)
dfls_222_2_1 <- c(3, 7, 0, 0.4)
B1_222 <- matrix(c(0.5, 0.2, -0.1, 0.3), nrow=2)
B2_222 <- matrix(c(0.4, -0.1, -0.2, 0.3), nrow=2)
theta_222logistikt_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                           vec(B1_222), vec(B2_222), c_and_gamma_222_2_1, dfls_222_2_1)

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="ind_skewed_t"
dfls_122_1_1 <- c(4, 13, -0.1, 0)
B1_122 <- matrix(c(1.2, -0.3, 0.7, 0.1), nrow=2)
B2_122 <- matrix(c(0.5, 0.2, -0.1, 3.1), nrow=2)
theta_122logikt_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(B1_122), vec(B2_122), gamma1_122_1_1, dfls_122_1_1)

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t"
dfls_123_1_1 <- c(10, 12, 3, 0.1, -0.2, 0.3)
B1_123 <- matrix(c(1.0, 0.3, 0.1, -0.8, 1.1, -0.5, -0.1, -0.2, 0.4), nrow=3)
B2_123 <- matrix(c(0.3, -0.2, -0.7, -0.8, 1.2, 0.5, 0.1, -0.2, 1.1), nrow=3)
theta_123expikt_1_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(B1_123),
                        vec(B2_123), c_and_gamma_123_1_1, dfls_123_1_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t"
dfls_132_1_1 <- c(30, 6, 0.6, -0.7)
B1_132 <- matrix(c(0.6, 0.2, -0.1, 0.7), nrow=2)
B2_132 <- matrix(c(0.4, -0.1, -0.2, 0.5), nrow=2)
B3_132 <- matrix(c(0.9, -0.5, 0.2, 0.4), nrow=2)
theta_132thresikt_1_1 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                          vec(B1_132), vec(B2_132), vec(B3_132), r1_132_1_1, r2_132_1_1, dfls_132_1_1)

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="ind_skewed_t",
# mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))
dfls_122_2_1 <- c(4, 13, -0.1, -0.2)
theta_222expcmwikt_2_1 <- c(phi10_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222), xi_222expcmw_2_1, dfls_122_2_1)
theta_222expcmwikt_2_1_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                    vec(B1_222), vec(B2_222), c(0.01, 0.33), dfls_122_2_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t",
# mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))
theta_132thresmwikt_1_1 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(B1_132),
                            vec(B2_132), vec(B3_132), dfls_132_1_1)
theta_132thresmwikt_1_1_expanded <- c(phi10_132, phi20_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(B1_132),
                                     vec(B2_132), vec(B3_132), 0, 1.2, dfls_132_1_1)

# p=1, M=2, p=3, weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="ind_skewed_t",
# mean_constraints=list(1:2), AR_constraints=C_123
dfls_123_3_1 <- c(11, 3, 20, 0.1, 0.2, 0.3)
theta_123logisticcmikt_3_1 <- c(phi10_123, vec(A11_123), vec(B1_123), vec(B2_123), c_and_gamma_123_3_1, dfls_123_3_1)
theta_123logisticcmikt_3_1_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vec(B1_123), vec(B2_123),
                                        c_and_gamma_123_3_1, dfls_123_3_1)

# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_skewed_t", AR_constraints=C_222
theta_222logcikt_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222), gamma1_222_2_1, dfls_222_2_1)
theta_222logcikt_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                              vec(B1_222), vec(B2_222), gamma1_222_2_1, dfls_222_2_1)

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_skewed_t", AR_constraints=C_222
theta_222exoikt_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222), dfls_222_2_1)
theta_222exoikt_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                             vec(B1_222), vec(B2_222), dfls_222_2_1)

#####################
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

# p=2, M=3, d=2, weight_function="exogenous", weightfun_pars=cbind(c(0.1, 0.4, 0.8), c(0.2, 0.5, 0.1), c(0.7, 0.1, 0.1)), cond_dist="Student",
# identification="heteroskedasticity"
dfs_232_1_1 <- c(0.3, 0.7)
W_232 <- W_132; lambdas2_232 <- lambdas2_132; lambdas3_232 <- lambdas3_132
theta_232thresitsh_1_1 <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                            vec(A31_232), vec(A32_232), vec(W_232), lambdas2_232, lambdas3_232, dfs_232_1_1)

##############
## ind_student

# p=2, M=2, d=2, cond_dist="ind_Student", weight_function="logistic", weightfun_pars=c(2, 1), identification="non-Gaussianity"
theta_222logistitng_2_1 <- theta_222logistit_2_1

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="ind_Student", identification="non-Gaussianity"
theta_122logitng_1_1 <- theta_122logit_1_1

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="ind_Student", identification="non-Gaussianity"
theta_123expitng_1_1 <- theta_123expit_1_1

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student", identification="non-Gaussianity"
theta_132thresitng_1_1 <- theta_132thresit_1_1

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="ind_Student", identification="non-Gaussianity",
# mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
theta_222expcmwitng_2_1 <- theta_222expcmwit_2_1
theta_222expcmwitng_2_1_expanded <- theta_222expcmwit_2_1_expanded

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student",
# mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2)), identification="non-Gaussianity"
theta_132thresmwitng_1_1 <- theta_132thresmwit_1_1
theta_132thresmwitng_1_1_expanded <- theta_132thresmwit_1_1_expanded

# p=1, M=2, p=3, weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="ind_Student",
# mean_constraints=list(1:2), AR_constraints=C_123, identification="non-Gaussianity"
theta_123logisticcmitng_3_1 <- theta_123logisticcmit_3_1
theta_123logisticcmitng_3_1_expanded <- theta_123logisticcmit_3_1_expanded

# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_Student", AR_constraints=C_222,
# identification="non-Gaussianity"
theta_222logcitng_2_1 <- theta_222logcit_2_1
theta_222logcitng_expanded <- theta_222logcit_expanded

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_Student",
# AR_constraints=C_222, identification="non-Gaussianity"
theta_222exoitng_2_1 <- theta_222exoit_2_1
theta_222exoitng_expanded <- theta_222exoit_expanded

###############
## ind_skewed_t

# p=2, M=2, d=2, cond_dist="ind_skewed_t", weight_function="logistic", weightfun_pars=c(2, 1), identification="non-Gaussianity"
theta_222logistiktng_2_1 <- theta_222logistikt_2_1

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="ind_skewed_t", identification="non-Gaussianity"
theta_122logiktng_1_1 <- theta_122logikt_1_1

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t", identification="non-Gaussianity"
theta_123expiktng_1_1 <- theta_123expikt_1_1

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t", identification="non-Gaussianity"
theta_132thresiktng_1_1 <- theta_132thresikt_1_1

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="ind_skewed_t", identification="non-Gaussianity",
# mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
theta_222expcmwiktng_2_1 <- theta_222expcmwikt_2_1
theta_222expcmwiktng_2_1_expanded <- theta_222expcmwikt_2_1_expanded

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t",
# mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2)), identification="non-Gaussianity"
theta_132thresmwiktng_1_1 <- theta_132thresmwikt_1_1
theta_132thresmwiktng_1_1_expanded <- theta_132thresmwikt_1_1_expanded

# p=1, M=2, p=3, weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="ind_skewed_t",
# mean_constraints=list(1:2), AR_constraints=C_123, identification="non-Gaussianity"
theta_123logisticcmiktng_3_1 <- theta_123logisticcmikt_3_1
theta_123logisticcmiktng_3_1_expanded <- theta_123logisticcmikt_3_1_expanded

# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_skewed_t", AR_constraints=C_222,
# identification="non-Gaussianity"
theta_222logciktng_2_1 <- theta_222logcikt_2_1
theta_222logciktng_expanded <- theta_222logcikt_expanded

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_skewed_t",
# AR_constraints=C_222, identification="non-Gaussianity"
theta_222exoiktng_2_1 <- theta_222exoikt_2_1
theta_222exoiktng_expanded <- theta_222exoikt_expanded


#########################################
## Structural models imposing constraints

# p=1, M=2, d=2, weight_function="relative_dens", identification="heteroskedasticity", AR_constraints=C_122
W_122 <- matrix(c(-0.03, 0.24, -0.76, -0.02), nrow=2, ncol=2, byrow=FALSE)
lambdas_122 <- c(3.36, 0.86)
theta_122relgshc <- c(phi10_122, phi20_122, vec(A11_122), vec(W_122), lambdas_122, alpha1_122)
theta_122relgshc_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A11_122), vec(W_122), lambdas_122, alpha1_122)

# p=1, M=3, d=2, weight_function="relative_dens", identification="heteroskedasticity",
# B_constraints=matrix(c(0.1, NA, 0.3, 0), nrow=2)
W_132b <- matrix(c(0.11, 0.22, 0.33, 0), nrow=2);
theta_132relgshb <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                      Wvec(W_132b), lambdas2_132, lambdas3_132, alpha1_132, alpha2_132)
theta_132relgshb_expanded <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                               vec(W_132b), lambdas2_132, lambdas3_132, alpha1_132, alpha2_132)

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student",
# identification="heteroskedasticity", mean_constraints=list(1:2), B_constraints=matrix(c(0.1, 0, 0, 0.3), nrow=2)
W_222b <- matrix(c(0.12, 0, 0, 0.31), nrow=2)
theta_222logistictshmb_2_1 <- c(phi10_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                                Wvec(W_222b), lambdas_222, c_and_gamma_222_2_1, df_222_2_1)
theta_222logistictshmb_2_1_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                                         vec(W_222b), lambdas_222, c_and_gamma_222_2_1, df_222_2_1)

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

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="ind_Student", identification="non-Gaussianity",
# B_constraints=matrix(c(1, NA, -1, 1) , nrow=2)
dfs_122_1_1 <- c(4, 13)
B1_122c <- matrix(c(1.2, 0.3, -0.7, 0.1), nrow=2)
B2_122c <- matrix(c(0.5, -0.9, -0.1, 3.1), nrow=2)
theta_122logitngb_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), Wvec(B1_122c), Wvec(B2_122c), gamma1_122_1_1, dfs_122_1_1)

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="ind_Student", identification="non-Gaussianity",
# B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)
dfs_123_1_1 <- c(10, 12, 3)
B1_123c <- matrix(c(1.0, 0.3, 0.1, 0, 1.1, -0.5, 0, -0.2, 0.4), nrow=3)
B2_123c <- matrix(c(0.3, -0.2, -0.7, 0, 1.2, 0.5, 0, -0.2, 1.1), nrow=3)
theta_123expitngb_1_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), Wvec(B1_123c),
                        Wvec(B2_123c), c_and_gamma_123_1_1, dfs_123_1_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student", identification="non-Gaussianity",
# B_constraints=matrix(c(1, 0, NA, 1), nrow=2)
dfs_132_1_1 <- c(30, 6)
B1_132c <- matrix(c(0.6, 0, -0.1, 0.7), nrow=2)
B2_132c <- matrix(c(0.4, 0, 0.2, 0.5), nrow=2)
B3_132c <- matrix(c(0.9, 0, -0.2, 0.4), nrow=2)
theta_132thresitngb_1_1 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                          Wvec(B1_132c), Wvec(B2_132c), Wvec(B3_132c), r1_132_1_1, r2_132_1_1, dfs_132_1_1)

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="ind_Student",
# mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
# identification="non-Gaussianity", matrix(c(NA, -1, 0, 1), nrow=2)
theta_222expcmwitngb_2_1 <- c(phi10_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c), xi_222expcmw_2_1, dfs_222_2_1)
theta_222expcmwitngb_2_1_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                       Wvec(B1_222c), Wvec(B2_222c), c(0.01, 0.33), dfs_222_2_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student", mean_constraints=list(1, 2:3),
# weight_constraints=list(R=0, r=c(0, 1.2)), identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2)
theta_132thresmwitngb_1_1 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), Wvec(B1_132c),
                               Wvec(B2_132c), Wvec(B3_132c), dfs_132_1_1)
theta_132thresmwitngb_1_1_expanded <- c(phi10_132, phi20_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), Wvec(B1_132c),
                                        Wvec(B2_132c), Wvec(B3_132c), 0, 1.2, dfs_132_1_1)

# p=1, M=2, d=3, weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="ind_Student", mean_constraints=list(1:2),
# AR_constraints=C_123, identification="non-Gaussianity", B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)
dfs_123_3_1 <- c(11, 3, 20)
theta_123logisticcmitngb_3_1 <- c(phi10_123, vec(A11_123), Wvec(B1_123c), Wvec(B2_123c), c_and_gamma_123_3_1, dfs_123_3_1)
theta_123logisticcmitngb_3_1_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), Wvec(B1_123c), Wvec(B2_123c),
                                           c_and_gamma_123_3_1, dfs_123_3_1)

# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_Student", AR_constraints=C_222,
# identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)
theta_222logcitngb_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c), gamma1_222_2_1, dfs_222_2_1)
theta_222logcitngb_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                 Wvec(B1_222c), Wvec(B2_222c), gamma1_222_2_1, dfs_222_2_1)

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_Student",
# AR_constraints=C_222, identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)
theta_222exoitngb_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c), dfs_222_2_1)
theta_222exoitngb_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                Wvec(B1_222c), Wvec(B2_222c), dfs_222_2_1)

################
### ind_skewed_t

# p=2, M=2, d=2, cond_dist="ind_skewed_t", weight_function="logistic", weightfun_pars=c(2, 1), identification="non-Gaussianity",
# B_constraints=matrix(c(NA, -1, 0, 1), nrow=2)
B1_222c <- matrix(c(0.5, -0.2, 0, 0.1), nrow=2)
B2_222c <- matrix(c(-0.4, -0.1, 0, 0.2), nrow=2)
theta_222logistiktngb_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                              Wvec(B1_222c), Wvec(B2_222c), c_and_gamma_222_2_1, dfls_222_2_1)

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="ind_skewed_t", identification="non-Gaussianity",
# B_constraints=matrix(c(1, NA, -1, 1) , nrow=2)
B1_122c <- matrix(c(1.2, 0.3, -0.7, 0.1), nrow=2)
B2_122c <- matrix(c(0.5, -0.9, -0.1, 3.1), nrow=2)
theta_122logiktngb_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), Wvec(B1_122c), Wvec(B2_122c), gamma1_122_1_1, dfls_122_1_1)

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t", identification="non-Gaussianity",
# B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)
B1_123c <- matrix(c(1.0, 0.3, 0.1, 0, 1.1, -0.5, 0, -0.2, 0.4), nrow=3)
B2_123c <- matrix(c(0.3, -0.2, -0.7, 0, 1.2, 0.5, 0, -0.2, 1.1), nrow=3)
theta_123expiktngb_1_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), Wvec(B1_123c),
                           Wvec(B2_123c), c_and_gamma_123_1_1, dfls_123_1_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t", identification="non-Gaussianity",
# B_constraints=matrix(c(1, 0, NA, 1), nrow=2)
B1_132c <- matrix(c(0.6, 0, -0.1, 0.7), nrow=2)
B2_132c <- matrix(c(0.4, 0, 0.2, 0.5), nrow=2)
B3_132c <- matrix(c(0.9, 0, -0.2, 0.4), nrow=2)
theta_132thresiktngb_1_1 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                             Wvec(B1_132c), Wvec(B2_132c), Wvec(B3_132c), r1_132_1_1, r2_132_1_1, dfls_132_1_1)

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="ind_skewed_t",
# mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
# identification="non-Gaussianity", matrix(c(NA, -1, 0, 1), nrow=2)
theta_222expcmwiktngb_2_1 <- c(phi10_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c), xi_222expcmw_2_1, dfls_222_2_1)
theta_222expcmwiktngb_2_1_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                       Wvec(B1_222c), Wvec(B2_222c), c(0.01, 0.33), dfls_222_2_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t", mean_constraints=list(1, 2:3),
# weight_constraints=list(R=0, r=c(0, 1.2)), identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2)
theta_132thresmwiktngb_1_1 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), Wvec(B1_132c),
                               Wvec(B2_132c), Wvec(B3_132c), dfls_132_1_1)
theta_132thresmwiktngb_1_1_expanded <- c(phi10_132, phi20_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), Wvec(B1_132c),
                                        Wvec(B2_132c), Wvec(B3_132c), 0, 1.2, dfls_132_1_1)

# p=1, M=2, d=3, weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="ind_skewed_t", mean_constraints=list(1:2),
# AR_constraints=C_123, identification="non-Gaussianity", B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)
theta_123logisticcmiktngb_3_1 <- c(phi10_123, vec(A11_123), Wvec(B1_123c), Wvec(B2_123c), c_and_gamma_123_3_1, dfls_123_3_1)
theta_123logisticcmiktngb_3_1_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), Wvec(B1_123c), Wvec(B2_123c),
                                           c_and_gamma_123_3_1, dfls_123_3_1)

# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_skewed_t", AR_constraints=C_222,
# identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)
theta_222logciktngb_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c), gamma1_222_2_1, dfls_222_2_1)
theta_222logciktngb_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                 Wvec(B1_222c), Wvec(B2_222c), gamma1_222_2_1, dfls_222_2_1)

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_skewed_t",
# AR_constraints=C_222, identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)
theta_222exoiktngb_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c), dfls_222_2_1)
theta_222exoiktngb_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                Wvec(B1_222c), Wvec(B2_222c), dfls_222_2_1)

##############################################
## Structural models with illegal W or lambdas

# p=1, M=2, d=2, weight_function="relative_dens", identification="heteroskedasticity"
W_122 <- matrix(c(-0.03, 0.24, -0.76, -0.02), nrow=2, ncol=2, byrow=FALSE)
lambdas_122_e <- c(3.36, 0.00)
theta_122relgsh_e <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), lambdas_122_e, alpha1_122)

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student",
# identification="heteroskedasticity"
W_222 <- W_122; lambdas_222_e <- c(-0.1, 1.2)
theta_222logistictsh_2_1_e <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                                vec(W_222), lambdas_222_e, c_and_gamma_222_2_1, df_222_2_1)

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student",
# identification="heteroskedasticity", mean_constraints=list(1:2), B_constraints=matrix(c(0.1, 0, 0, 0.3), nrow=2)
W_222b_e <- matrix(c(0.12, 0, 0, -0.001), nrow=2)
theta_222logistictshmb_2_1_e <- c(phi10_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                                  Wvec(W_222b_e), lambdas_222, c_and_gamma_222_2_1, df_222_2_1)
theta_222logistictshmb_2_1_e_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                                           vec(W_222b_e), lambdas_222, c_and_gamma_222_2_1, df_222_2_1)

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=1), identification="heteroskedasticity",
# weight_constraints=list(R=0, r=c(0.1, 0.2, 0.3)), B_constraints=matrix(c(0.1, 0.2, 0.3, 0), nrow=2)
W_122b_e <- matrix(c(0.11, -0.22, 0.33, 0), nrow=2)
theta_122logshwb_12_1_e <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), Wvec(W_122b_e), lambdas_122)

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), identification="heteroskedasticity",
# AR_constraints=C_123, weight_constraints=list(R=matrix(c(1, 0.5), nrow=2), r=c(0, 0))
# B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE)
W_123b_e <- matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, 0.19), nrow=3, ncol=3, byrow=FALSE)
theta_123expshcwb_1_1_e <- c(phi10_123, phi20_123, vec(A11_123), Wvec(W_123b_e), lambdas_123, 0.6)
theta_123expshcwb_1_1_e_expanded <- c(phi10_123, phi20_123, vec(A11_123), vec(A11_123), vec(W_123b_e), lambdas_123, 0.6, 0.3)


# p=1, M=2, d=2, weight_function="exogenous", weightfun_pars=list(vars=1:2, lags=1), identification="heteroskedasticity",
# weight_constraints=list(R=0, r=c(0.1, 0.2, 0.3)), B_constraints=matrix(c(0.1, 0.2, 0.3, 0), nrow=2)
theta_122exoshwb_12_1_e <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), Wvec(W_122b_e), lambdas_122)


### Non-Gaussian structural models with illegal impact matrices ###

##############
## ind_Student

# p=1, M=1, d=2, cond_dist="ind_Student", weight_function="threshold", weightfun_pars=c(1, 1)
dfs_112 <- c(3, 7)
B1_112e <- matrix(c(0.5, 0.2, -0.5, -0.2), nrow=2) # Singular
theta_112itnge <- c(phi10_112, vec(A11_112), vec(B1_112e), dfs_112)

# p=2, M=2, d=2, cond_dist="ind_Student", weight_function="logistic", weightfun_pars=c(2, 1), identification="non-Gaussianity",
# B_constraints=matrix(c(NA, NA, 0, -1), nrow=2)
B2_222c_e <- matrix(c(0.4, -0.1, 1, -0.2), nrow=2)
theta_222logistitngb_2_1e <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                               Wvec(B1_222c), Wvec(B2_222c_e), dfs_222_2_1)

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="ind_Student", identification="non-Gaussianity",
# B_constraints=matrix(c(1, NA, -1, 1) , nrow=2)
dfs_122_1_1 <- c(4, 13)
B1_122c_e <- matrix(c(1.2, 0.3, 0.7, 0.1), nrow=2)
theta_122logitngb_1_1e <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), Wvec(B1_122c_e), Wvec(B2_122c), gamma1_122_1_1, dfs_122_1_1)
theta_122logitngb_1_1e_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), Wvec(B1_122c_e), Wvec(B2_122c), gamma1_122_1_1, dfs_122_1_1)

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="ind_Student", identification="non-Gaussianity",
# B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)
dfs_123_1_1 <- c(10, 12, 3)
B2_123c_e <- matrix(c(0.3, -0.2, -0.7, 0, 1.2, 0.5, 0, -0.2, -1.1), nrow=3)
theta_123expitngb_1_1e <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), Wvec(B1_123c),
                            Wvec(B2_123c_e), c_and_gamma_123_1_1, dfs_123_1_1)
theta_123expitngb_1_1e_expanded <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(B1_123c), vec(B2_123c_e),
                                     c_and_gamma_123_1_1, dfs_123_1_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student", identification="non-Gaussianity",
# B_constraints=matrix(c(1, 0, NA, 1), nrow=2)
dfs_132_1_1 <- c(30, 6)
B3_132c_e <- matrix(c(0.9, -0.5, -0.2, 0.4), nrow=2)
theta_132thresitngb_1_1e <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                              Wvec(B1_132c), Wvec(B2_132c), Wvec(B3_132c_e), r1_132_1_1, r2_132_1_1, dfs_132_1_1)

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="ind_Student",
# mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
# identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)
dfs_122_2_1 <- c(4, 13)
B1_222c_e <- matrix(c(0, 1, 0, 0.1), nrow=2)
theta_222expcmwitngb_2_1e <- c(phi10_222, vec(A11_222), vec(A12_222), Wvec(B1_222c_e), Wvec(B2_222c), xi_222expcmw_2_1, dfs_122_2_1)
theta_222expcmwitngb_2_1e_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                        vec(B1_222c_e), vec(B2_222c), c(0.01, 0.33), dfs_122_2_1)

B1_222c_e2 <- matrix(c(-1, 1, 0, 0.1), nrow=2)
theta_222expcmwitngb_2_1e2 <- c(phi10_222, vec(A11_222), vec(A12_222), Wvec(B1_222c_e2), Wvec(B2_222c), xi_222expcmw_2_1, dfs_122_2_1)
theta_222expcmwitngb_2_1e2_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                        vec(B1_222c_e2), vec(B2_222c), c(0.01, 0.33), dfs_122_2_1)

B1_222c_e3 <- matrix(c(0.1, 1, 0, 0.2), nrow=2)
theta_222expcmwitngb_2_1e3 <- c(phi10_222, vec(A11_222), vec(A12_222), Wvec(B1_222c_e3), Wvec(B2_222c), xi_222expcmw_2_1, dfs_122_2_1)
theta_222expcmwitngb_2_1e3_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                         vec(B1_222c_e3), vec(B2_222c), c(0.01, 0.33), dfs_122_2_1)

B1_222c_ok <- matrix(c(1, 1, 0, 0.2), nrow=2)
theta_222expcmwitngb_2_1ok <- c(phi10_222, vec(A11_222), vec(A12_222), Wvec(B1_222c_ok), Wvec(B2_222c), xi_222expcmw_2_1, dfs_122_2_1)
theta_222expcmwitngb_2_1ok_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                         vec(B1_222c_ok), vec(B2_222c), c(0.01, 0.33), dfs_122_2_1)




# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student", mean_constraints=list(1, 2:3),
# weight_constraints=list(R=0, r=c(0, 1.2)), identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2)
B1_132c_e <- matrix(c(-0.9, 0, -0.2, 0.4), nrow=2)
theta_132thresmwitngb_1_1e <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), Wvec(B1_132c),
                                Wvec(B2_132c), Wvec(B3_132c), dfs_132_1_1)
theta_132thresmwitngb_1_1e_expanded <- c(phi10_132, phi20_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(B1_132c_e),
                                         vec(B2_132c), vec(B3_132c), 0, 1.2, dfs_132_1_1)

# p=1, M=2, p=3, weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="ind_Student", mean_constraints=list(1:2),
# AR_constraints=C_123, identification="non-Gaussianity", B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)
dfs_123_3_1 <- c(11, 3, 20)
B1_123c_e <- matrix(c(1, 0.1, 0.2, 0, 1, -0.1, 1, 0.7, 1), nrow=3)
theta_123logisticcmitngb_3_1e <- c(phi10_123, vec(A11_123), Wvec(B1_123c_e), Wvec(B2_123c), c_and_gamma_123_3_1, dfs_123_3_1)
theta_123logisticcmitngb_3_1e_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vec(B1_123c_e), vec(B2_123c),
                                            c_and_gamma_123_3_1, dfs_123_3_1)

B1_123c_e3 <- matrix(c(1, 0.1, 0.2, 0, 0.9, -0.1, 0, -0.7, 1), nrow=3)
theta_123logisticcmitngb_3_1e3 <- c(phi10_123, vec(A11_123), Wvec(B1_123c_e3), Wvec(B2_123c), c_and_gamma_123_3_1, dfs_123_3_1)
theta_123logisticcmitngb_3_1e3_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vec(B1_123c_e3), vec(B2_123c),
                                            c_and_gamma_123_3_1, dfs_123_3_1)

B1_123c_e4 <- matrix(c(1, 0.1, 0.2, 0, 0.7, -0.1, 0, 0.9, 1), nrow=3)
theta_123logisticcmitngb_3_1e4 <- c(phi10_123, vec(A11_123), Wvec(B1_123c_e4), Wvec(B2_123c), c_and_gamma_123_3_1, dfs_123_3_1)
theta_123logisticcmitngb_3_1e4_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vec(B1_123c_e4), vec(B2_123c),
                                             c_and_gamma_123_3_1, dfs_123_3_1)

B1_123c_ok <- matrix(c(1, 0.1, 0.2, 0, 0.7, -0.1, 0, 0.5, 1), nrow=3)
theta_123logisticcmitngb_3_1ok <- c(phi10_123, vec(A11_123), Wvec(B1_123c_ok), Wvec(B2_123c), c_and_gamma_123_3_1, dfs_123_3_1)
theta_123logisticcmitngb_3_1ok_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vec(B1_123c_ok), vec(B2_123c),
                                             c_and_gamma_123_3_1, dfs_123_3_1)


# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_Student", AR_constraints=C_222,
# identification="non-Gaussianity", B_constraints=matrix(c(1, NA, NA, 1), nrow=2)
B1_222c_e2 <- matrix(c(0.4, -0.1, 1, 0.2), nrow=2) # Ok but different to B1_222c
B2_222c_e2 <- matrix(c(0.4, 0.2, 0.2, 0.1), nrow=2) # Singular
theta_222logcitngb_2_1e <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222c_e2), vec(B2_222c_e2), gamma1_222_2_1, dfs_222_2_1)
theta_222logcitngbe_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                  vec(B1_222c_e2), vec(B2_222c_e2), gamma1_222_2_1, dfs_222_2_1)

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_Student",
# AR_constraints=C_222, identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, NA, NA), nrow=2)
B1_222c_e3 <- B2_222c_e2 # Singular
B2_222c_e3 <- B1_222c_e2 # Ok
theta_222exoitngb_2_1e <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222c_e3), vec(B2_222c_e3), dfs_222_2_1)
theta_222exoitngbe_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                 vec(B1_222c_e3), vec(B2_222c_e3), dfs_222_2_1)

# p=1, M=2, p=3, weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="ind_Student", mean_constraints=list(1:2),
# AR_constraints=C_123, identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA), nrow=3)
dfs_123_3_1 <- c(11, 3, 20)
B1_123c_e2 <- matrix(c(1, 0.1, 0.2, 0.1, 1, -0.1, 1, 0.7, 1), nrow=3) # Ok but different to B1_123c
B2_123c_e2 <- matrix(c(1, 0.1, 0.2, 0.1, 1, -0.1, 2, 0.2, 0.4), nrow=3) # Singular
theta_123logisticcmitngb_3_1e <- c(phi10_123, vec(A11_123), vec(B1_123c), vec(B2_123c), c_and_gamma_123_3_1, dfs_123_3_1)
theta_123logisticcmitngb_3_1e_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vec(B1_123c_e2), vec(B2_123c_e2),
                                            c_and_gamma_123_3_1, dfs_123_3_1)

###############
## ind_skewed_t

# p=1, M=1, d=2, cond_dist="ind_skewed_t", weight_function="threshold", weightfun_pars=c(1, 1)
B1_112e <- matrix(c(0.5, 0.2, -0.5, -0.2), nrow=2) # Singular
theta_112iktnge <- c(phi10_112, vec(A11_112), vec(B1_112e), dfls_112)

# p=2, M=2, d=2, cond_dist="ind_skewed_t", weight_function="logistic", weightfun_pars=c(2, 1), identification="non-Gaussianity",
# B_constraints=matrix(c(NA, NA, 0, -1), nrow=2)
B2_222c_e <- matrix(c(0.4, -0.1, 1, -0.2), nrow=2)
theta_222logistiktngb_2_1e <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                               Wvec(B1_222c), Wvec(B2_222c_e), dfls_222_2_1)

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="ind_skewed_t", identification="non-Gaussianity",
# B_constraints=matrix(c(1, NA, -1, 1) , nrow=2)
B1_122c_e <- matrix(c(1.2, 0.3, 0.7, 0.1), nrow=2)
theta_122logiktngb_1_1e <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), Wvec(B1_122c_e), Wvec(B2_122c), gamma1_122_1_1, dfls_122_1_1)
theta_122logiktngb_1_1e_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), Wvec(B1_122c_e), Wvec(B2_122c), gamma1_122_1_1, dfls_122_1_1)

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t", identification="non-Gaussianity",
# B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)
B2_123c_e <- matrix(c(0.3, -0.2, -0.7, 0, 1.2, 0.5, 0, -0.2, -1.1), nrow=3)
theta_123expiktngb_1_1e <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), Wvec(B1_123c),
                            Wvec(B2_123c_e), c_and_gamma_123_1_1, dfls_123_1_1)
theta_123expiktngb_1_1e_expanded <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(B1_123c), vec(B2_123c_e),
                                     c_and_gamma_123_1_1, dfls_123_1_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t", identification="non-Gaussianity",
# B_constraints=matrix(c(1, 0, NA, 1), nrow=2)
B3_132c_e <- matrix(c(0.9, -0.5, -0.2, 0.4), nrow=2)
theta_132thresiktngb_1_1e <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                              Wvec(B1_132c), Wvec(B2_132c), Wvec(B3_132c_e), r1_132_1_1, r2_132_1_1, dfls_132_1_1)

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="ind_skewed_t",
# mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
# identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)
B1_222c_e <- matrix(c(0, 1, 0, 0.1), nrow=2)
theta_222expcmwiktngb_2_1e <- c(phi10_222, vec(A11_222), vec(A12_222), Wvec(B1_222c_e), Wvec(B2_222c), xi_222expcmw_2_1, dfls_122_2_1)
theta_222expcmwiktngb_2_1e_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                        vec(B1_222c_e), vec(B2_222c), c(0.01, 0.33), dfls_122_2_1)

B1_222c_e2 <- matrix(c(-1, 1, 0, 0.1), nrow=2)
theta_222expcmwiktngb_2_1e2 <- c(phi10_222, vec(A11_222), vec(A12_222), Wvec(B1_222c_e2), Wvec(B2_222c), xi_222expcmw_2_1, dfls_122_2_1)
theta_222expcmwiktngb_2_1e2_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                         vec(B1_222c_e2), vec(B2_222c), c(0.01, 0.33), dfls_122_2_1)

B1_222c_e3 <- matrix(c(0.1, 1, 0, 0.2), nrow=2)
theta_222expcmwiktngb_2_1e3 <- c(phi10_222, vec(A11_222), vec(A12_222), Wvec(B1_222c_e3), Wvec(B2_222c), xi_222expcmw_2_1, dfls_122_2_1)
theta_222expcmwiktngb_2_1e3_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                         vec(B1_222c_e3), vec(B2_222c), c(0.01, 0.33), dfls_122_2_1)

B1_222c_ok <- matrix(c(1, 1, 0, 0.2), nrow=2)
theta_222expcmwiktngb_2_1ok <- c(phi10_222, vec(A11_222), vec(A12_222), Wvec(B1_222c_ok), Wvec(B2_222c), xi_222expcmw_2_1, dfls_122_2_1)
theta_222expcmwiktngb_2_1ok_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                         vec(B1_222c_ok), vec(B2_222c), c(0.01, 0.33), dfls_122_2_1)




# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t", mean_constraints=list(1, 2:3),
# weight_constraints=list(R=0, r=c(0, 1.2)), identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2)
B1_132c_e <- matrix(c(-0.9, 0, -0.2, 0.4), nrow=2)
theta_132thresmwiktngb_1_1e <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), Wvec(B1_132c),
                                Wvec(B2_132c), Wvec(B3_132c), dfls_132_1_1)
theta_132thresmwiktngb_1_1e_expanded <- c(phi10_132, phi20_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(B1_132c_e),
                                         vec(B2_132c), vec(B3_132c), 0, 1.2, dfls_132_1_1)

# p=1, M=2, p=3, weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="ind_skewed_t", mean_constraints=list(1:2),
# AR_constraints=C_123, identification="non-Gaussianity", B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)
B1_123c_e <- matrix(c(1, 0.1, 0.2, 0, 1, -0.1, 1, 0.7, 1), nrow=3)
theta_123logisticcmiktngb_3_1e <- c(phi10_123, vec(A11_123), Wvec(B1_123c_e), Wvec(B2_123c), c_and_gamma_123_3_1, dfls_123_3_1)
theta_123logisticcmiktngb_3_1e_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vec(B1_123c_e), vec(B2_123c),
                                            c_and_gamma_123_3_1, dfls_123_3_1)

B1_123c_e3 <- matrix(c(1, 0.1, 0.2, 0, 0.9, -0.1, 0, -0.7, 1), nrow=3)
theta_123logisticcmiktngb_3_1e3 <- c(phi10_123, vec(A11_123), Wvec(B1_123c_e3), Wvec(B2_123c), c_and_gamma_123_3_1, dfls_123_3_1)
theta_123logisticcmiktngb_3_1e3_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vec(B1_123c_e3), vec(B2_123c),
                                             c_and_gamma_123_3_1, dfls_123_3_1)

B1_123c_e4 <- matrix(c(1, 0.1, 0.2, 0, 0.7, -0.1, 0, 0.9, 1), nrow=3)
theta_123logisticcmiktngb_3_1e4 <- c(phi10_123, vec(A11_123), Wvec(B1_123c_e4), Wvec(B2_123c), c_and_gamma_123_3_1, dfls_123_3_1)
theta_123logisticcmiktngb_3_1e4_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vec(B1_123c_e4), vec(B2_123c),
                                             c_and_gamma_123_3_1, dfls_123_3_1)

B1_123c_ok <- matrix(c(1, 0.1, 0.2, 0, 0.7, -0.1, 0, 0.5, 1), nrow=3)
theta_123logisticcmiktngb_3_1ok <- c(phi10_123, vec(A11_123), Wvec(B1_123c_ok), Wvec(B2_123c), c_and_gamma_123_3_1, dfls_123_3_1)
theta_123logisticcmiktngb_3_1ok_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vec(B1_123c_ok), vec(B2_123c),
                                             c_and_gamma_123_3_1, dfls_123_3_1)


# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_skewed_t", AR_constraints=C_222,
# identification="non-Gaussianity", B_constraints=matrix(c(1, NA, NA, 1), nrow=2)
B1_222c_e2 <- matrix(c(0.4, -0.1, 1, 0.2), nrow=2) # Ok but different to B1_222c
B2_222c_e2 <- matrix(c(0.4, 0.2, 0.2, 0.1), nrow=2) # Singular
theta_222logciktngb_2_1e <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222c_e2), vec(B2_222c_e2), gamma1_222_2_1, dfls_222_2_1)
theta_222logciktngbe_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                  vec(B1_222c_e2), vec(B2_222c_e2), gamma1_222_2_1, dfls_222_2_1)

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_skewed_t",
# AR_constraints=C_222, identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, NA, NA), nrow=2)
B1_222c_e3 <- B2_222c_e2 # Singular
B2_222c_e3 <- B1_222c_e2 # Ok
theta_222exoiktngb_2_1e <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222c_e3), vec(B2_222c_e3), dfls_222_2_1)
theta_222exoiktngbe_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                 vec(B1_222c_e3), vec(B2_222c_e3), dfls_222_2_1)

# p=1, M=2, p=3, weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="ind_skewed_t", mean_constraints=list(1:2),
# AR_constraints=C_123, identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, NA, NA, NA, NA, NA, NA, NA), nrow=3)
B1_123c_e2 <- matrix(c(1, 0.1, 0.2, 0.1, 1, -0.1, 1, 0.7, 1), nrow=3) # Ok but different to B1_123c
B2_123c_e2 <- matrix(c(1, 0.1, 0.2, 0.1, 1, -0.1, 2, 0.2, 0.4), nrow=3) # Singular
theta_123logisticcmiktngb_3_1e <- c(phi10_123, vec(A11_123), vec(B1_123c), vec(B2_123c), c_and_gamma_123_3_1, dfls_123_3_1)
theta_123logisticcmiktngb_3_1e_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vec(B1_123c_e2), vec(B2_123c_e2),
                                            c_and_gamma_123_3_1, dfls_123_3_1)


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

Omegas_112it <- pick_Omegas(p=1, M=1, d=2, params=theta_112it, cond_dist="ind_Student")
Omegas_112it_singular <- pick_Omegas(p=1, M=1, d=2, params=theta_112itnge, cond_dist="ind_Student")
Omegas_122it <- pick_Omegas(p=1, M=2, d=2, params=theta_122logit_1_1, cond_dist="ind_Student")
Omegas_222it <- pick_Omegas(p=2, M=2, d=2, params=theta_222logistit_2_1, cond_dist="ind_Student")
Omegas_222it_singular <- pick_Omegas(p=2, M=2, d=2, params=theta_222logcitngbe_expanded, cond_dist="ind_Student")
Omegas_222it_singular2 <- pick_Omegas(p=2, M=2, d=2, params=theta_222exoitngbe_expanded, cond_dist="ind_Student")
Omegas_123it <- pick_Omegas(p=1, M=2, d=3, params=theta_123expit_1_1, cond_dist="ind_Student")
Omegas_123it_singular <- pick_Omegas(p=1, M=2, d=3, params=theta_123logisticcmitngb_3_1e_expanded, cond_dist="ind_Student")
Omegas_132it <- pick_Omegas(p=1, M=3, d=2, params=theta_132thresit_1_1, cond_dist="ind_Student")

Omegas_122ite <- pick_Omegas(p=1, M=2, d=2, params=theta_122logitngb_1_1e_expanded, cond_dist="ind_Student", identification="non-Gaussianity")
Omegas_123ite <- pick_Omegas(p=1, M=2, d=3, params=theta_123expitngb_1_1e_expanded, cond_dist="ind_Student", identification="non-Gaussianity")
Omegas_132ite <- pick_Omegas(p=1, M=3, d=2, params=theta_132thresmwitngb_1_1e_expanded, cond_dist="ind_Student", identification="non-Gaussianity")

Omegas_222ite2 <- pick_Omegas(p=2, M=2, d=2, params=theta_222expcmwitngb_2_1e2_expanded, cond_dist="ind_Student", identification="non-Gaussianity")
Omegas_222ite3 <- pick_Omegas(p=2, M=2, d=2, params=theta_222expcmwitngb_2_1e3_expanded, cond_dist="ind_Student", identification="non-Gaussianity")
Omegas_222itok <- pick_Omegas(p=2, M=2, d=2, params=theta_222expcmwitngb_2_1ok_expanded, cond_dist="ind_Student", identification="non-Gaussianity")
Omegas_123ite3 <- pick_Omegas(p=1, M=2, d=3, params=theta_123logisticcmitngb_3_1e3_expanded, cond_dist="ind_Student",
                              identification="non-Gaussianity")
Omegas_123ite4 <- pick_Omegas(p=1, M=2, d=3, params=theta_123logisticcmitngb_3_1e4_expanded, cond_dist="ind_Student",
                              identification="non-Gaussianity")
Omegas_123itok <- pick_Omegas(p=1, M=2, d=3, params=theta_123logisticcmitngb_3_1ok_expanded, cond_dist="ind_Student",
                              identification="non-Gaussianity")

Omegas_112ikt <- pick_Omegas(p=1, M=1, d=2, params=theta_112ikt, cond_dist="ind_skewed_t")
Omegas_112ikt_singular <- pick_Omegas(p=1, M=1, d=2, params=theta_112iktnge, cond_dist="ind_skewed_t")
Omegas_122ikt <- pick_Omegas(p=1, M=2, d=2, params=theta_122logikt_1_1, cond_dist="ind_skewed_t")
Omegas_222ikt <- pick_Omegas(p=2, M=2, d=2, params=theta_222logistikt_2_1, cond_dist="ind_skewed_t")
Omegas_222ikt_singular <- pick_Omegas(p=2, M=2, d=2, params=theta_222logciktngbe_expanded, cond_dist="ind_skewed_t")
Omegas_222ikt_singular2 <- pick_Omegas(p=2, M=2, d=2, params=theta_222exoiktngbe_expanded, cond_dist="ind_skewed_t")
Omegas_123ikt <- pick_Omegas(p=1, M=2, d=3, params=theta_123expikt_1_1, cond_dist="ind_skewed_t")
Omegas_123ikt_singular <- pick_Omegas(p=1, M=2, d=3, params=theta_123logisticcmiktngb_3_1e_expanded, cond_dist="ind_skewed_t")
Omegas_132ikt <- pick_Omegas(p=1, M=3, d=2, params=theta_132thresikt_1_1, cond_dist="ind_skewed_t")

Omegas_122ikte <- pick_Omegas(p=1, M=2, d=2, params=theta_122logiktngb_1_1e_expanded, cond_dist="ind_skewed_t", identification="non-Gaussianity")
Omegas_123ikte <- pick_Omegas(p=1, M=2, d=3, params=theta_123expiktngb_1_1e_expanded, cond_dist="ind_skewed_t", identification="non-Gaussianity")
Omegas_132ikte <- pick_Omegas(p=1, M=3, d=2, params=theta_132thresmwiktngb_1_1e_expanded, cond_dist="ind_skewed_t", identification="non-Gaussianity")

Omegas_222ikte2 <- pick_Omegas(p=2, M=2, d=2, params=theta_222expcmwiktngb_2_1e2_expanded, cond_dist="ind_skewed_t", identification="non-Gaussianity")
Omegas_222ikte3 <- pick_Omegas(p=2, M=2, d=2, params=theta_222expcmwiktngb_2_1e3_expanded, cond_dist="ind_skewed_t", identification="non-Gaussianity")
Omegas_222iktok <- pick_Omegas(p=2, M=2, d=2, params=theta_222expcmwiktngb_2_1ok_expanded, cond_dist="ind_skewed_t", identification="non-Gaussianity")
Omegas_123ikte3 <- pick_Omegas(p=1, M=2, d=3, params=theta_123logisticcmiktngb_3_1e3_expanded, cond_dist="ind_skewed_t",
                              identification="non-Gaussianity")
Omegas_123ikte4 <- pick_Omegas(p=1, M=2, d=3, params=theta_123logisticcmiktngb_3_1e4_expanded, cond_dist="ind_skewed_t",
                              identification="non-Gaussianity")
Omegas_123iktok <- pick_Omegas(p=1, M=2, d=3, params=theta_123logisticcmiktngb_3_1ok_expanded, cond_dist="ind_skewed_t",
                              identification="non-Gaussianity")


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

weightpars_122logit <- pick_weightpars(p=1, M=2, d=2, params=theta_122logit_1_1, weight_function="mlogit", weightfun_pars=c(1, 1),
                                       cond_dist="ind_Student")
weightpars_222logit <- pick_weightpars(p=2, M=2, d=2, params=theta_222logistit_2_1, weight_function="logistic", weightfun_pars=c(2, 1),
                                       cond_dist="ind_Student")
weightpars_123expit <- pick_weightpars(p=1, M=2, d=3, params=theta_123expit_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                                       cond_dist="ind_Student")
weightpars_132thresit <- pick_weightpars(p=1, M=3, d=2, params=theta_132thresit_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                        cond_dist="ind_Student")

weightpars_122logikt <- pick_weightpars(p=1, M=2, d=2, params=theta_122logikt_1_1, weight_function="mlogit", weightfun_pars=c(1, 1),
                                       cond_dist="ind_skewed_t")
weightpars_222logikt <- pick_weightpars(p=2, M=2, d=2, params=theta_222logistikt_2_1, weight_function="logistic", weightfun_pars=c(2, 1),
                                       cond_dist="ind_skewed_t")
weightpars_123expikt <- pick_weightpars(p=1, M=2, d=3, params=theta_123expikt_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                                       cond_dist="ind_skewed_t")
weightpars_132thresikt <- pick_weightpars(p=1, M=3, d=2, params=theta_132thresikt_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                         cond_dist="ind_skewed_t")

Omegas_222log_12_2 <- pick_Omegas(p=2, M=2, d=2, params=theta_222log_12_2)
boldA_222log_12_2 <- form_boldA(p=2, M=2, d=2, all_A=pick_allA(p=2, M=2, d=2, params=theta_222log_12_2))
weightpars_222log_12_2 <- pick_weightpars(p=2, M=2, d=2, params=theta_222log_12_2, weight_function="mlogit",
                                          weightfun_pars=list(vars=1:2, lags=2), cond_dist="Gaussian")


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

  # Check impact matrices
  expect_true(in_paramspace(p=1, M=1, d=2, weight_function="threshold", cond_dist="ind_Student",
                            all_boldA=boldA_112, all_Omegas=Omegas_112it, weightpars=numeric(0), distpars=dfs_112))
  expect_false(in_paramspace(p=1, M=1, d=2, weight_function="threshold", cond_dist="ind_Student",
                            all_boldA=boldA_112, all_Omegas=Omegas_112it_singular, weightpars=numeric(0), distpars=dfs_112))
  expect_true(in_paramspace(p=1, M=2, d=2, weight_function="mlogit", cond_dist="ind_Student",
                            all_boldA=boldA_122, all_Omegas=Omegas_122it, weightpars=weightpars_122logit, distpars=dfs_122_1_1))
  expect_true(in_paramspace(p=2, M=2, d=2, weight_function="logistic", cond_dist="ind_Student", all_boldA=boldA_222,
                            all_Omegas=Omegas_222it, weightpars=weightpars_222logit, distpars=dfs_222_2_1))
  expect_false(in_paramspace(p=2, M=2, d=2, weight_function="logistic", cond_dist="ind_Student", all_boldA=boldA_222,
                            all_Omegas=Omegas_222it_singular, weightpars=weightpars_222logit, distpars=dfs_222_2_1))
  expect_false(in_paramspace(p=2, M=2, d=2, weight_function="logistic", cond_dist="ind_Student", all_boldA=boldA_222,
                             all_Omegas=Omegas_222it_singular2, weightpars=weightpars_222logit, distpars=dfs_222_2_1))
  expect_true(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_Student", all_boldA=boldA_123,
                            all_Omegas=Omegas_123it, weightpars=weightpars_123expit, distpars=dfs_123_1_1))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_Student", all_boldA=boldA_123,
                            all_Omegas=Omegas_123it_singular, weightpars=weightpars_123expit, distpars=dfs_123_1_1))
  expect_true(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_Student", all_boldA=boldA_132,
                            all_Omegas=Omegas_132it, weightpars=weightpars_132thresit, distpars=dfs_132_1_1))

  expect_true(in_paramspace(p=1, M=1, d=2, weight_function="threshold", cond_dist="ind_skewed_t",
                            all_boldA=boldA_112, all_Omegas=Omegas_112ikt, weightpars=numeric(0), distpars=dfls_112))
  expect_false(in_paramspace(p=1, M=1, d=2, weight_function="threshold", cond_dist="ind_skewed_t",
                             all_boldA=boldA_112, all_Omegas=Omegas_112ikt_singular, weightpars=numeric(0), distpars=dfls_112))
  expect_true(in_paramspace(p=1, M=2, d=2, weight_function="mlogit", cond_dist="ind_skewed_t",
                            all_boldA=boldA_122, all_Omegas=Omegas_122ikt, weightpars=weightpars_122logikt, distpars=dfls_122_1_1))
  expect_true(in_paramspace(p=2, M=2, d=2, weight_function="logistic", cond_dist="ind_skewed_t", all_boldA=boldA_222,
                            all_Omegas=Omegas_222ikt, weightpars=weightpars_222logikt, distpars=dfls_222_2_1))
  expect_false(in_paramspace(p=2, M=2, d=2, weight_function="logistic", cond_dist="ind_skewed_t", all_boldA=boldA_222,
                             all_Omegas=Omegas_222ikt_singular, weightpars=weightpars_222logikt, distpars=dfls_222_2_1))
  expect_false(in_paramspace(p=2, M=2, d=2, weight_function="logistic", cond_dist="ind_skewed_t", all_boldA=boldA_222,
                             all_Omegas=Omegas_222ikt_singular2, weightpars=weightpars_222logikt, distpars=dfls_222_2_1))
  expect_true(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_skewed_t", all_boldA=boldA_123,
                            all_Omegas=Omegas_123ikt, weightpars=weightpars_123expikt, distpars=dfls_123_1_1))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_skewed_t", all_boldA=boldA_123,
                             all_Omegas=Omegas_123ikt_singular, weightpars=weightpars_123expikt, distpars=dfls_123_1_1))
  expect_true(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_skewed_t", all_boldA=boldA_132,
                            all_Omegas=Omegas_132ikt, weightpars=weightpars_132thresikt, distpars=dfls_132_1_1))


  # Check sign constraints in B_1,...,B_M
  expect_false(in_paramspace(p=1, M=2, d=2, params=theta_122logitngb_1_1e_expanded, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1),
                             cond_dist="ind_Student", identification="non-Gaussianity", B_constraints=matrix(c(1, NA, -1, 1), nrow=2),
                             all_boldA=boldA_122, all_Omegas=Omegas_122ite, weightpars=weightpars_122logit, distpars=dfs_122_1_1))
  expect_false(in_paramspace(p=1, M=2, d=3, params=theta_123expitngb_1_1e_expanded, weight_function="exponential", weightfun_pars=c(1, 1),
                             cond_dist="ind_Student", identification="non-Gaussianity",
                             B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3),
                             all_boldA=boldA_123, all_Omegas=Omegas_123ite, weightpars=weightpars_123expit, distpars=dfs_123_1_1))
  expect_false(in_paramspace(p=1, M=3, d=2, params=theta_132thresmwitngb_1_1e_expanded, weight_function="threshold", weightfun_pars=c(1, 1),
                             cond_dist="ind_Student", identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2),
                             all_boldA=boldA_132, all_Omegas=Omegas_132ite, weightpars=weightpars_132thresit, distpars=dfs_132_1_1))

  expect_false(in_paramspace(p=1, M=2, d=2, params=theta_122logiktngb_1_1e_expanded, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1),
                             cond_dist="ind_skewed_t", identification="non-Gaussianity", B_constraints=matrix(c(1, NA, -1, 1), nrow=2),
                             all_boldA=boldA_122, all_Omegas=Omegas_122ikte, weightpars=weightpars_122logikt, distpars=dfls_122_1_1))
  expect_false(in_paramspace(p=1, M=2, d=3, params=theta_123expiktngb_1_1e_expanded, weight_function="exponential", weightfun_pars=c(1, 1),
                             cond_dist="ind_skewed_t", identification="non-Gaussianity",
                             B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3),
                             all_boldA=boldA_123, all_Omegas=Omegas_123ikte, weightpars=weightpars_123expikt, distpars=dfls_123_1_1))
  expect_false(in_paramspace(p=1, M=3, d=2, params=theta_132thresmwiktngb_1_1e_expanded, weight_function="threshold", weightfun_pars=c(1, 1),
                             cond_dist="ind_skewed_t", identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2),
                             all_boldA=boldA_132, all_Omegas=Omegas_132ikte, weightpars=weightpars_132thresikt, distpars=dfls_132_1_1))

  # Check other_constraints: that the impact matrix of the first regime has the first non-zero element in each column
  # strictly positive and that they are in decreasing order.
  expect_false(in_paramspace(p=1, M=1, d=2, params=theta_112it, cond_dist="ind_Student", weight_function="threshold",
                             weightfun_pars=c(1, 1), identification="non-Gaussianity",
                             other_constraints=list(B1_constraints="fixed_sign_and_order"),
                             all_boldA=boldA_112, all_Omegas=Omegas_112it, weightpars=c(1, 2), distpars=c(7, 9)))
  Omegas_112it2 <- array(c(0.5, 0.2, 0.3, -0.1), dim=c(2, 2, 1))
  expect_true(in_paramspace(p=1, M=1, d=2, params=theta_112it, cond_dist="ind_Student", weight_function="threshold",
                            weightfun_pars=c(1, 1), identification="non-Gaussianity",
                            other_constraints=list(B1_constraints="fixed_sign_and_order"),
                            all_boldA=boldA_112, all_Omegas=Omegas_112it2, weightpars=c(1, 2), distpars=c(7, 9)))
  expect_false(in_paramspace(p=2, M=2, d=2, params=theta_222expcmwitngb_2_1e2_expanded, weight_function="exponential", weightfun_pars=c(2, 1),
                             cond_dist="ind_Student",  identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2),
                             other_constraints=list(B1_constraints="fixed_sign_and_order"),
                             all_boldA=boldA_222, all_Omegas=Omegas_222ite2, weightpars=c(1, 2), distpars=c(7, 9)))
  expect_false(in_paramspace(p=2, M=2, d=2, params=theta_222expcmwitngb_2_1e3_expanded, weight_function="exponential", weightfun_pars=c(2, 1),
                             cond_dist="ind_Student",  identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2),
                             other_constraints=list(B1_constraints="fixed_sign_and_order"),
                             all_boldA=boldA_222, all_Omegas=Omegas_222ite3, weightpars=c(1, 2), distpars=c(7, 9)))
  expect_true(in_paramspace(p=2, M=2, d=2, params=theta_222expcmwitngb_2_1ok_expanded, weight_function="exponential", weightfun_pars=c(2, 1),
                            cond_dist="ind_Student",  identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2),
                            other_constraints=list(B1_constraints="fixed_sign_and_order"),
                            all_boldA=boldA_222, all_Omegas=Omegas_222itok, weightpars=c(1, 2), distpars=c(7, 9)))
  expect_false(in_paramspace(p=1, M=2, d=3, params=theta_123logisticcmitngb_3_1e3_expanded, weight_function="logistic", weightfun_pars=c(3, 1),
                             cond_dist="ind_Student",  identification="non-Gaussianity",
                             B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3),
                             other_constraints=list(B1_constraints="fixed_sign_and_order"),
                             all_boldA=boldA_123, all_Omegas=Omegas_123ite3, weightpars=c(1, 2), distpars=c(7, 9, 10)))
  expect_false(in_paramspace(p=1, M=2, d=3, params=theta_123logisticcmitngb_3_1e4_expanded, weight_function="logistic", weightfun_pars=c(3, 1),
                             cond_dist="ind_Student",  identification="non-Gaussianity",
                             B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3),
                             other_constraints=list(B1_constraints="fixed_sign_and_order"),
                             all_boldA=boldA_123, all_Omegas=Omegas_123ite4, weightpars=c(1, 2), distpars=c(7, 9, 10)))
  expect_true(in_paramspace(p=1, M=2, d=3, params=theta_123logisticcmitngb_3_1ok_expanded, weight_function="logistic", weightfun_pars=c(3, 1),
                            cond_dist="ind_Student",  identification="non-Gaussianity",
                            B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3),
                            other_constraints=list(B1_constraints="fixed_sign_and_order"),
                            all_boldA=boldA_123, all_Omegas=Omegas_123itok, weightpars=c(1, 2), distpars=c(7, 9, 10)))

  expect_false(in_paramspace(p=1, M=1, d=2, params=theta_112ikt, cond_dist="ind_skewed_t", weight_function="threshold",
                             weightfun_pars=c(1, 1), identification="non-Gaussianity",
                             other_constraints=list(B1_constraints="fixed_sign_and_order"),
                             all_boldA=boldA_112, all_Omegas=Omegas_112ikt, weightpars=c(1, 2), distpars=c(7, 9, 0.1, 0.2)))
  Omegas_112ikt2 <- array(c(0.5, 0.2, 0.3, -0.1), dim=c(2, 2, 1))
  expect_true(in_paramspace(p=1, M=1, d=2, params=theta_112ikt, cond_dist="ind_skewed_t", weight_function="threshold",
                            weightfun_pars=c(1, 1), identification="non-Gaussianity",
                            other_constraints=list(B1_constraints="fixed_sign_and_order"),
                            all_boldA=boldA_112, all_Omegas=Omegas_112ikt2, weightpars=c(1, 2), distpars=c(7, 9, 0.1, 0.2)))
  expect_false(in_paramspace(p=2, M=2, d=2, params=theta_222expcmwiktngb_2_1e2_expanded, weight_function="exponential", weightfun_pars=c(2, 1),
                             cond_dist="ind_skewed_t",  identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2),
                             other_constraints=list(B1_constraints="fixed_sign_and_order"),
                             all_boldA=boldA_222, all_Omegas=Omegas_222ikte2, weightpars=c(1, 2), distpars=c(7, 9, 0.1, 0.2)))
  expect_false(in_paramspace(p=2, M=2, d=2, params=theta_222expcmwiktngb_2_1e3_expanded, weight_function="exponential", weightfun_pars=c(2, 1),
                             cond_dist="ind_skewed_t",  identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2),
                             other_constraints=list(B1_constraints="fixed_sign_and_order"),
                             all_boldA=boldA_222, all_Omegas=Omegas_222ikte3, weightpars=c(1, 2), distpars=c(7, 9, 0.1, 0.2)))
  expect_true(in_paramspace(p=2, M=2, d=2, params=theta_222expcmwiktngb_2_1ok_expanded, weight_function="exponential", weightfun_pars=c(2, 1),
                            cond_dist="ind_skewed_t",  identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2),
                            other_constraints=list(B1_constraints="fixed_sign_and_order"),
                            all_boldA=boldA_222, all_Omegas=Omegas_222iktok, weightpars=c(1, 2), distpars=c(7, 9, 0.1, 0.2)))
  expect_false(in_paramspace(p=1, M=2, d=3, params=theta_123logisticcmiktngb_3_1e3_expanded, weight_function="logistic", weightfun_pars=c(3, 1),
                             cond_dist="ind_skewed_t",  identification="non-Gaussianity",
                             B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3),
                             other_constraints=list(B1_constraints="fixed_sign_and_order"),
                             all_boldA=boldA_123, all_Omegas=Omegas_123ikte3, weightpars=c(1, 2), distpars=c(7, 9, 10, 0.1, 0.2, 0.3)))
  expect_false(in_paramspace(p=1, M=2, d=3, params=theta_123logisticcmiktngb_3_1e4_expanded, weight_function="logistic", weightfun_pars=c(3, 1),
                             cond_dist="ind_skewed_t",  identification="non-Gaussianity",
                             B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3),
                             other_constraints=list(B1_constraints="fixed_sign_and_order"),
                             all_boldA=boldA_123, all_Omegas=Omegas_123ikte4, weightpars=c(1, 2), distpars=c(7, 9, 10, 0.1, 0.2, 0.3)))
  expect_true(in_paramspace(p=1, M=2, d=3, params=theta_123logisticcmiktngb_3_1ok_expanded, weight_function="logistic", weightfun_pars=c(3, 1),
                            cond_dist="ind_skewed_t",  identification="non-Gaussianity",
                            B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3),
                            other_constraints=list(B1_constraints="fixed_sign_and_order"),
                            all_boldA=boldA_123, all_Omegas=Omegas_123iktok, weightpars=c(1, 2), distpars=c(7, 9, 10, 0.1, 0.2 ,0.3)))


  # Check signs constraints in W
  expect_false(in_paramspace(p=2, M=2, d=2, params=theta_222logistictshmb_2_1_e_expanded, weight_function="logistic",
                             weightfun_pars=c(2, 1), cond_dist="Student", identification="heteroskedasticity",
                             B_constraints=matrix(c(0.1, 0, 0, 0.3), nrow=2), all_boldA=boldA_222, all_Omegas=Omegas_222,
                             weightpars=weightpars_222logit, distpars=df_222_2_1))
  expect_false(in_paramspace(p=1, M=2, d=3, params=theta_123expshcwb_1_1_e_expanded, weight_function="exponential",
                             cond_dist="Gaussian", weightfun_pars=c(1, 1), identification="heteroskedasticity",
                             B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE),
                             all_boldA=boldA_123, all_Omegas=Omegas_123, weightpars=c(0.6, 0.3), distpars=numeric(0)))

  # Check invertibility of B_t for all t
  expect_true(in_paramspace(p=2, M=2, d=2, weight_function="logistic", cond_dist="ind_Student", all_boldA=boldA_222,
                            all_Omegas=Omegas_222it, weightpars=weightpars_222logit, distpars=dfs_222_2_1,
                            transition_weights=cbind(c(0.1, 0.2, 0.3), c(0.9, 0.8, 0.7))))
  Omegas_222it_singt <- array(c(1:4, -(1:4)), dim=c(2, 2, 2))
  expect_false(in_paramspace(p=2, M=2, d=2, weight_function="logistic", cond_dist="ind_Student", all_boldA=boldA_222,
                             all_Omegas=Omegas_222it_singt , weightpars=weightpars_222logit, distpars=dfs_222_2_1,
                             transition_weights=cbind(c(0.1, 0.2, 0.5), c(0.9, 0.8, 0.5))))
  expect_false(in_paramspace(p=2, M=2, d=2, weight_function="logistic", cond_dist="ind_Student", all_boldA=boldA_222,
                             all_Omegas=Omegas_222it_singt, weightpars=weightpars_222logit, distpars=dfs_222_2_1,
                             transition_weights=cbind(c(0.1, 0.5, 0.7), c(0.9, 0.5, 0.3))))
  Omegas_123it_singt <- array(c(1:9, -(1:9)), dim=c(3, 3, 2))
  expect_true(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_Student", all_boldA=boldA_123,
                            all_Omegas=Omegas_123it, weightpars=weightpars_123expit, distpars=dfs_123_1_1,
                            transition_weights=cbind(c(0.1, 0.2, 0.3), c(0.9, 0.8, 0.7))))
  expect_true(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_Student", all_boldA=boldA_123,
                            all_Omegas=Omegas_123it, weightpars=weightpars_123expit, distpars=dfs_123_1_1,
                            transition_weights=cbind(c(0.1, 0.5, 0), c(0.9, 0.5, 1))))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_Student", all_boldA=boldA_123,
                            all_Omegas=Omegas_123it_singt, weightpars=weightpars_123expit, distpars=dfs_123_1_1,
                            transition_weights=cbind(c(0.1, 0.5, 0), c(0.9, 0.5, 1))))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_Student", all_boldA=boldA_123,
                             all_Omegas=Omegas_123it_singt, weightpars=weightpars_123expit, distpars=dfs_123_1_1,
                             transition_weights=cbind(c(0.5, 0, 0), c(0.5, 1, 1))))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_Student", all_boldA=boldA_123,
                             all_Omegas=Omegas_123it_singt, weightpars=weightpars_123expit, distpars=dfs_123_1_1,
                             transition_weights=cbind(c(0.3, 0.3, 0.5), c(0.7, 0.7, 0.5))))
  Omegas_132it_singt <- array(c(1:4, 0.1, -0.2, 0.7, 0.5, -(1:4)), dim=c(2, 2, 3))
  Omegas_132it_singt2 <- array(c(1:4, 1:4, -2, -4, -6, -8), dim=c(2, 2, 3))
  expect_true(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_Student", all_boldA=boldA_132,
                            all_Omegas=Omegas_132it, weightpars=weightpars_132thresit, distpars=dfs_132_1_1,
                            transition_weights=cbind(c(0.1, 0.5, 0.5), c(0.9, 0.1, 0), c(0, 0.4, 0.5))))
  expect_true(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_Student", all_boldA=boldA_132,
                            all_Omegas=Omegas_132it, weightpars=weightpars_132thresit, distpars=dfs_132_1_1,
                            transition_weights=cbind(c(0.1, 0.5, 1, 0.5), c(0.9, 0.1, 0, 0), c(0, 0.4, 0, 0.5))))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_Student", all_boldA=boldA_132,
                            all_Omegas=Omegas_132it_singt, weightpars=weightpars_132thresit, distpars=dfs_132_1_1,
                            transition_weights=cbind(c(0.1, 0.5, 0.5), c(0.9, 0.1, 0), c(0, 0.4, 0.5))))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_Student", all_boldA=boldA_132,
                             all_Omegas=Omegas_132it_singt, weightpars=weightpars_132thresit, distpars=dfs_132_1_1,
                             transition_weights=cbind(c(0.1, 0.5, 1, 0.5), c(0.9, 0.1, 0, 0), c(0, 0.4, 0, 0.5))))
  expect_true(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_Student", all_boldA=boldA_132,
                             all_Omegas=Omegas_132it, weightpars=weightpars_132thresit, distpars=dfs_132_1_1,
                             transition_weights=cbind(c(0.1, 0.5, 0.25), c(0.9, 0.1, 0.25), c(0, 0.4, 0.25))))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_Student", all_boldA=boldA_132,
                             all_Omegas=Omegas_132it_singt2, weightpars=weightpars_132thresit, distpars=dfs_132_1_1,
                             transition_weights=cbind(c(0.1, 0.5, 0.25), c(0.9, 0.1, 0.25), c(0, 0.4, 0.25))))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_Student", all_boldA=boldA_132,
                             all_Omegas=Omegas_132it_singt2, weightpars=weightpars_132thresit, distpars=dfs_132_1_1,
                             transition_weights=cbind(c(0.1, 0.5, 0.25, 1), c(0.9, 0.1, 0.25, 0), c(0, 0.4, 0.25, 0))))

  expect_true(in_paramspace(p=2, M=2, d=2, weight_function="logistic", cond_dist="ind_skewed_t", all_boldA=boldA_222,
                            all_Omegas=Omegas_222ikt, weightpars=weightpars_222logikt, distpars=dfls_222_2_1,
                            transition_weights=cbind(c(0.1, 0.2, 0.3), c(0.9, 0.8, 0.7))))
  Omegas_222ikt_singt <- array(c(1:4, -(1:4)), dim=c(2, 2, 2))
  expect_false(in_paramspace(p=2, M=2, d=2, weight_function="logistic", cond_dist="ind_skewed_t", all_boldA=boldA_222,
                             all_Omegas=Omegas_222ikt_singt , weightpars=weightpars_222logikt, distpars=dfls_222_2_1,
                             transition_weights=cbind(c(0.1, 0.2, 0.5), c(0.9, 0.8, 0.5))))
  expect_false(in_paramspace(p=2, M=2, d=2, weight_function="logistic", cond_dist="ind_skewed_t", all_boldA=boldA_222,
                             all_Omegas=Omegas_222ikt_singt, weightpars=weightpars_222logikt, distpars=dfls_222_2_1,
                             transition_weights=cbind(c(0.1, 0.5, 0.7), c(0.9, 0.5, 0.3))))
  Omegas_123ikt_singt <- array(c(1:9, -(1:9)), dim=c(3, 3, 2))
  expect_true(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_skewed_t", all_boldA=boldA_123,
                            all_Omegas=Omegas_123ikt, weightpars=weightpars_123expikt, distpars=dfls_123_1_1,
                            transition_weights=cbind(c(0.1, 0.2, 0.3), c(0.9, 0.8, 0.7))))
  expect_true(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_skewed_t", all_boldA=boldA_123,
                            all_Omegas=Omegas_123ikt, weightpars=weightpars_123expikt, distpars=dfls_123_1_1,
                            transition_weights=cbind(c(0.1, 0.5, 0), c(0.9, 0.5, 1))))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_skewed_t", all_boldA=boldA_123,
                             all_Omegas=Omegas_123ikt_singt, weightpars=weightpars_123expikt, distpars=dfls_123_1_1,
                             transition_weights=cbind(c(0.1, 0.5, 0), c(0.9, 0.5, 1))))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_skewed_t", all_boldA=boldA_123,
                             all_Omegas=Omegas_123ikt_singt, weightpars=weightpars_123expikt, distpars=dfls_123_1_1,
                             transition_weights=cbind(c(0.5, 0, 0), c(0.5, 1, 1))))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_skewed_t", all_boldA=boldA_123,
                             all_Omegas=Omegas_123ikt_singt, weightpars=weightpars_123expikt, distpars=dfls_123_1_1,
                             transition_weights=cbind(c(0.3, 0.3, 0.5), c(0.7, 0.7, 0.5))))
  Omegas_132ikt_singt <- array(c(1:4, 0.1, -0.2, 0.7, 0.5, -(1:4)), dim=c(2, 2, 3))
  Omegas_132ikt_singt2 <- array(c(1:4, 1:4, -2, -4, -6, -8), dim=c(2, 2, 3))
  expect_true(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_skewed_t", all_boldA=boldA_132,
                            all_Omegas=Omegas_132ikt, weightpars=weightpars_132thresikt, distpars=dfls_132_1_1,
                            transition_weights=cbind(c(0.1, 0.5, 0.5), c(0.9, 0.1, 0), c(0, 0.4, 0.5))))
  expect_true(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_skewed_t", all_boldA=boldA_132,
                            all_Omegas=Omegas_132ikt, weightpars=weightpars_132thresikt, distpars=dfls_132_1_1,
                            transition_weights=cbind(c(0.1, 0.5, 1, 0.5), c(0.9, 0.1, 0, 0), c(0, 0.4, 0, 0.5))))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_skewed_t", all_boldA=boldA_132,
                             all_Omegas=Omegas_132ikt_singt, weightpars=weightpars_132thresikt, distpars=dfls_132_1_1,
                             transition_weights=cbind(c(0.1, 0.5, 0.5), c(0.9, 0.1, 0), c(0, 0.4, 0.5))))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_skewed_t", all_boldA=boldA_132,
                             all_Omegas=Omegas_132ikt_singt, weightpars=weightpars_132thresikt, distpars=dfls_132_1_1,
                             transition_weights=cbind(c(0.1, 0.5, 1, 0.5), c(0.9, 0.1, 0, 0), c(0, 0.4, 0, 0.5))))
  expect_true(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_skewed_t", all_boldA=boldA_132,
                            all_Omegas=Omegas_132ikt, weightpars=weightpars_132thresikt, distpars=dfls_132_1_1,
                            transition_weights=cbind(c(0.1, 0.5, 0.25), c(0.9, 0.1, 0.25), c(0, 0.4, 0.25))))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_skewed_t", all_boldA=boldA_132,
                             all_Omegas=Omegas_132ikt_singt2, weightpars=weightpars_132thresikt, distpars=dfls_132_1_1,
                             transition_weights=cbind(c(0.1, 0.5, 0.25), c(0.9, 0.1, 0.25), c(0, 0.4, 0.25))))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_skewed_t", all_boldA=boldA_132,
                             all_Omegas=Omegas_132ikt_singt2, weightpars=weightpars_132thresikt, distpars=dfls_132_1_1,
                             transition_weights=cbind(c(0.1, 0.5, 0.25, 1), c(0.9, 0.1, 0.25, 0), c(0, 0.4, 0.25, 0))))

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
  expect_false(in_paramspace(p=1, M=2, d=2, weight_function="logistic", cond_dist="Gaussian",
                             all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=c(0.1, -0.1)))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="logistic", cond_dist="Gaussian",
                             all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=c(0.1, 0)))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="logistic", cond_dist="Gaussian",
                             all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=c(0.1, 1e-10)))
  expect_true(in_paramspace(p=1, M=2, d=2, weight_function="logistic", cond_dist="Gaussian",
                             all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=c(0.0, 0.1)))
  expect_true(in_paramspace(p=1, M=2, d=3, weight_function="logistic", cond_dist="Gaussian",
                             all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=c(0.1, 0.01)))
  expect_false(in_paramspace(p=1, M=2, d=2, weight_function="exponential", cond_dist="Gaussian",
                             all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=c(0.1, -0.1)))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="Gaussian",
                             all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=c(0.1, 0)))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="Gaussian",
                             all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=c(0.1, 1e-10)))
  expect_true(in_paramspace(p=1, M=2, d=2, weight_function="exponential", cond_dist="Gaussian",
                            all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=c(0.0, 0.1)))
  expect_true(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="Gaussian",
                            all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=c(0.1, 0.01)))
  expect_true(in_paramspace(p=1, M=2, d=2, weight_function="threshold", cond_dist="Gaussian",
                            all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=c(-0.1)))
  expect_true(in_paramspace(p=1, M=2, d=2, weight_function="threshold", cond_dist="Gaussian",
                            all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=c(100)))
  expect_true(in_paramspace(p=1, M=2, d=3, weight_function="threshold", cond_dist="Gaussian",
                            all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=0))
  expect_true(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="Gaussian",
                            all_boldA=boldA_132, all_Omegas=Omegas_132, weightpars=c(1, 2)))
  expect_true(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="Gaussian",
                            all_boldA=boldA_132, all_Omegas=Omegas_132, weightpars=c(-0.1, 0.01)))
  expect_true(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="Gaussian",
                            all_boldA=boldA_132, all_Omegas=Omegas_132, weightpars=c(-0.1, 0.0)))
  expect_true(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="Gaussian",
                            all_boldA=boldA_132, all_Omegas=Omegas_132, weightpars=c(0, 0.01)))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="Gaussian",
                            all_boldA=boldA_132, all_Omegas=Omegas_132, weightpars=c(0.01, 0.0)))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="Gaussian",
                             all_boldA=boldA_132, all_Omegas=Omegas_132, weightpars=c(-0.01, -0.1)))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="Gaussian",
                             all_boldA=boldA_132, all_Omegas=Omegas_132, weightpars=c(1, -1)))

  # mlogit (nothing to check in weightpars, so just checks that the function runs)
  expect_true(in_paramspace(p=2, M=2, d=2, weight_function="mlogit", cond_dist="Gaussian",
                            all_boldA=boldA_222log_12_2, all_Omegas=Omegas_222log_12_2, weightpars=weightpars_222log_12_2))

  # Checks distribution params
  expect_true(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="Student",
                            all_boldA=boldA_132, all_Omegas=Omegas_132, weightpars=c(0, 0.01), distpars=3))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="Student",
                             all_boldA=boldA_132, all_Omegas=Omegas_132, weightpars=c(0, 0.01), distpars=2))
  expect_true(in_paramspace(p=1, M=2, d=2, weight_function="logistic", cond_dist="Student",
                            all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=c(0.0, 0.1), distpars=11))
  expect_false(in_paramspace(p=1, M=2, d=2, weight_function="logistic", cond_dist="Student",
                            all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=c(0.0, 0.1), distpars=0))
  expect_true(in_paramspace(p=1, M=2, d=2, weight_function="exponential", cond_dist="Student",
                            all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=c(0.0, 0.1), distpars=2.01))
  expect_false(in_paramspace(p=1, M=2, d=2, weight_function="exponential", cond_dist="Student",
                            all_boldA=boldA_122, all_Omegas=Omegas_122, weightpars=c(0.0, 0.1), distpars=1.99999))
  expect_true(in_paramspace(p=2, M=2, d=2, weight_function="mlogit", cond_dist="Student",
                            all_boldA=boldA_222log_12_2, all_Omegas=Omegas_222log_12_2, weightpars=weightpars_222log_12_2,
                            distpars=100))
  expect_false(in_paramspace(p=2, M=2, d=2, weight_function="mlogit", cond_dist="Student",
                            all_boldA=boldA_222log_12_2, all_Omegas=Omegas_222log_12_2, weightpars=weightpars_222log_12_2,
                            distpars=-3))

  expect_true(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_Student",
                            all_boldA=boldA_132, all_Omegas=Omegas_132, weightpars=c(0, 0.01), distpars=c(3, 7)))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_Student",
                             all_boldA=boldA_132, all_Omegas=Omegas_132, weightpars=c(0, 0.01), distpars=c(2, 8)))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_Student",
                             all_boldA=boldA_132, all_Omegas=Omegas_132, weightpars=c(0, 0.01), distpars=c(3, -8)))
  expect_true(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_Student",
                            all_boldA=boldA_123, all_Omegas=Omegas_123, weightpars=c(0.0, 0.1), distpars=c(2.01, 3, 99)))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_Student",
                             all_boldA=boldA_123, all_Omegas=Omegas_123, weightpars=c(0.0, 0.1), distpars=c(2.01, 3, 1.99)))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_Student",
                             all_boldA=boldA_123, all_Omegas=Omegas_123, weightpars=c(0.0, 0.1), distpars=c(2.01, 3, -2.99)))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_Student",
                             all_boldA=boldA_123, all_Omegas=Omegas_123, weightpars=c(0.0, 0.1), distpars=c(0, 3, 0)))
  expect_true(in_paramspace(p=2, M=2, d=2, weight_function="mlogit", cond_dist="ind_Student",
                            all_boldA=boldA_222log_12_2, all_Omegas=Omegas_222log_12_2, weightpars=weightpars_222log_12_2,
                            distpars=c(13, 13)))
  expect_false(in_paramspace(p=2, M=2, d=2, weight_function="mlogit", cond_dist="ind_Student",
                             all_boldA=boldA_222log_12_2, all_Omegas=Omegas_222log_12_2, weightpars=weightpars_222log_12_2,
                             distpars=c(13, -13)))

  expect_true(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_skewed_t",
                            all_boldA=boldA_132, all_Omegas=Omegas_132, weightpars=c(0, 0.01), distpars=c(3, 7, 0.1, 0.2)))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_skewed_t",
                             all_boldA=boldA_132, all_Omegas=Omegas_132, weightpars=c(0, 0.01), distpars=c(2, 8, 0.1, 0.2)))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_skewed_t",
                             all_boldA=boldA_132, all_Omegas=Omegas_132, weightpars=c(0, 0.01), distpars=c(3, -8, 0.1, 0.2)))
  expect_true(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_skewed_t",
                            all_boldA=boldA_123, all_Omegas=Omegas_123, weightpars=c(0.0, 0.1), distpars=c(2.01, 3, 99, 0.1, 0.2, 0.3)))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_skewed_t",
                             all_boldA=boldA_123, all_Omegas=Omegas_123, weightpars=c(0.0, 0.1), distpars=c(2.01, 3, 1.99, 0.1, 0.2, 0.3)))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_skewed_t",
                             all_boldA=boldA_123, all_Omegas=Omegas_123, weightpars=c(0.0, 0.1), distpars=c(2.01, 3, -2.99, 0.1, 0.2 ,0.3)))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_skewed_t",
                             all_boldA=boldA_123, all_Omegas=Omegas_123, weightpars=c(0.0, 0.1), distpars=c(0, 3, 0, 0.1 ,0.2 ,0.3)))
  expect_true(in_paramspace(p=2, M=2, d=2, weight_function="mlogit", cond_dist="ind_skewed_t",
                            all_boldA=boldA_222log_12_2, all_Omegas=Omegas_222log_12_2, weightpars=weightpars_222log_12_2,
                            distpars=c(13, 13, 0.1, 0.2)))
  expect_false(in_paramspace(p=2, M=2, d=2, weight_function="mlogit", cond_dist="ind_skewed_t",
                             all_boldA=boldA_222log_12_2, all_Omegas=Omegas_222log_12_2, weightpars=weightpars_222log_12_2,
                             distpars=c(13, -13, 0.1, 0.2)))

  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_skewed_t", all_boldA=boldA_132,
                             all_Omegas=Omegas_132ikt, weightpars=weightpars_132thresikt, distpars=c(11, 12, 0.1, 1)))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_skewed_t", all_boldA=boldA_123,
                             all_Omegas=Omegas_123ikt_singular, weightpars=weightpars_123expikt, distpars=c(11, 12, 13, 0.1, -1, 0.2)))
  expect_false(in_paramspace(p=1, M=3, d=2, weight_function="threshold", cond_dist="ind_skewed_t", all_boldA=boldA_132,
                             all_Omegas=Omegas_132ikt, weightpars=weightpars_132thresikt, distpars=c(11, 12, 0.1, -1.00001)))
  expect_false(in_paramspace(p=1, M=2, d=3, weight_function="exponential", cond_dist="ind_skewed_t", all_boldA=boldA_123,
                             all_Omegas=Omegas_123ikt_singular, weightpars=weightpars_123expikt, distpars=c(11, 12, 13, 0.1, 0.3, 1.00001)))
})

test_that("check_params work correctly", {
  # Check that no errors when all is fine
  check_params(p=2, M=2, d=2, params=theta_222relgcmw, weight_function="relative_dens", mean_constraints=list(1:2),
               AR_constraints=C_222, weight_constraints=list(R=0, r=0.6), cond_dist="Gaussian")
  check_params(p=1, M=2, d=2, params=theta_122logw_1_1, weight_function="mlogit",
               weightfun_pars=list(vars=1, lags=1), weight_constraints=list(R=0, r=c(0.12, 0.13)), cond_dist="Gaussian")
  check_params(p=2, M=2, d=2, params=theta_222logcmw_12_2, weight_function="mlogit",
               weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222,
               weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5), r=c(0, 0.11, 0.12, 0.13, 0)),
               cond_dist="Gaussian")
  check_params(p=2, M=2, d=2, params=theta_222logisticcmw_2_1, weight_function="logistic", weightfun_pars=c(2, 1),
               mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
               cond_dist="Gaussian")
  check_params(p=2, M=2, d=2, params=theta_222expcmw_2_1, weight_function="exponential", weightfun_pars=c(2, 1),
               mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
               cond_dist="Gaussian")
  check_params(p=2, M=2, d=2, params=theta_222thresc_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
               AR_constraints=C_222, cond_dist="Gaussian")

  check_params(p=2, M=2, d=2, params=c(theta_222logcmw_12_2, 3), weight_function="mlogit",
               weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222,
               weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5), r=c(0, 0.11, 0.12, 0.13, 0)),
               cond_dist="Student")
  check_params(p=2, M=2, d=2, params=c(theta_222logisticcmw_2_1, 2.01), weight_function="logistic", weightfun_pars=c(2, 1),
               mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
               cond_dist="Student")
  check_params(p=2, M=2, d=2, params=c(theta_222expcmw_2_1, 999), weight_function="exponential", weightfun_pars=c(2, 1),
               mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
               cond_dist="Student")
  check_params(p=2, M=2, d=2, params=c(theta_222thresc_1_1, 13), weight_function="threshold", weightfun_pars=c(1, 1),
               AR_constraints=C_222, cond_dist="Student")

  check_params(p=1, M=2, d=2, params=theta_122relgsh, weight_function="relative_dens", identification="heteroskedasticity")
  check_params(p=2, M=2, d=2, params=theta_222logistictsh_2_1, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student",
               identification="heteroskedasticity")
  check_params(p=2, M=2, d=2, params=theta_222logistictshmb_2_1, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student",
               identification="heteroskedasticity", mean_constraints=list(1:2), B_constraints=matrix(c(0.1, 0, 0, 0.3), nrow=2))
  check_params(p=1, M=2, d=2, params=theta_122logshwb_12_1, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=1),
               identification="heteroskedasticity", weight_constraints=list(R=0, r=c(0.1, 0.2, 0.3)),
               B_constraints=matrix(c(0.1, 0.2, 0.3, 0), nrow=2))
  check_params(p=1, M=2, d=3, params=theta_123expshcwb_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
               identification="heteroskedasticity", AR_constraints=C_123, weight_constraints=list(R=matrix(c(1, 0.5), nrow=2), r=c(0, 0)),
               B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE))

  check_params(p=1, M=1, d=2, params=theta_112it, cond_dist="ind_Student", weight_function="threshold", weightfun_pars=c(1, 1))
  check_params(p=2, M=2, d=2, params=theta_222logistit_2_1, cond_dist="ind_Student", weight_function="logistic", weightfun_pars=c(2, 1))
  check_params(p=2, M=2, d=2, params=theta_222logcit_2_1,  weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_Student",
               AR_constraints=C_222)
  check_params(data=matrix(0, nrow=5, ncol=2), p=2, M=2, d=2, params=theta_222exoit_2_1, weight_function="exogenous",
               weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_Student", AR_constraints=C_222)

  check_params(p=1, M=1, d=2, params=theta_112ikt, cond_dist="ind_skewed_t", weight_function="threshold", weightfun_pars=c(1, 1))
  check_params(p=2, M=2, d=2, params=theta_222logistikt_2_1, cond_dist="ind_skewed_t", weight_function="logistic", weightfun_pars=c(2, 1))
  check_params(p=2, M=2, d=2, params=theta_222logcikt_2_1,  weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_skewed_t",
               AR_constraints=C_222)
  check_params(data=matrix(0, nrow=5, ncol=2), p=2, M=2, d=2, params=theta_222exoikt_2_1, weight_function="exogenous",
               weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_skewed_t", AR_constraints=C_222)


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
  expect_error(check_params(p=1, M=3, d=2, params=theta_132relgw, weight_function="relative_dens",
                            weight_constraints=list(R=matrix(c(0.9, 0.5), nrow=2), r=c(0.13, 0.13)), cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=3, d=2, params=theta_132thresmw_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                            mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2)), cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=3, d=2, params=theta_132relgshb, weight_function="relative_dens", identification="heteroskedasticity",
                            B_constraints=matrix(c(0.1, NA, 0.3, 0), nrow=2)))
  expect_error(check_params(p=1, M=2, d=3, params=theta_123expsh_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                            identification="heteroskedasticity"))
  expect_error(check_params(p=2, M=3, d=2, params=theta_232threstshb_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                            cond_dist="Student", identification="heteroskedasticity", B_constraints=matrix(c(0.1, 0.2, -0.3, 0), nrow=2)))
  expect_error(check_params(p=1, M=2, d=3, params=theta_123expit_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                            cond_dist="ind_Student"))
  expect_error(check_params(p=1, M=3, d=2, params=theta_132thresmwit_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                            cond_dist="ind_Student", mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))))
  expect_error(check_params(p=1, M=3, d=2, params=theta_132thresitngb_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                            cond_dist="ind_Student", identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2)))

  expect_error(check_params(p=1, M=2, d=3, params=theta_123expikt_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                            cond_dist="ind_skewed_t"))
  expect_error(check_params(p=1, M=3, d=2, params=theta_132thresmwikt_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                            cond_dist="ind_skewed_t", mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))))
  expect_error(check_params(p=1, M=3, d=2, params=theta_132thresiktngb_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                            cond_dist="ind_skewed_t", identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2)))

  # Check with AR_constraints (reform params used inside so just checks that the function goes through)
  check_params(p=1, M=1, d=2, params=theta_112relgc, weight_function="relative_dens", cond_dist="Gaussian",
               AR_constraints=C_112)
  check_params(p=2, M=2, d=2, params=theta_222relgc, weight_function="relative_dens", cond_dist="Gaussian",
               AR_constraints=C_222)
  check_params(p=2, M=2, d=2, params=theta_222relgcm, weight_function="relative_dens", cond_dist="Gaussian",
               AR_constraints=C_222, mean_constraints=list(1:2))

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
  expect_error(check_params(p=1, M=3, d=2, params=theta_132relgw_notpd, weight_function="relative_dens",
                            weight_constraints=list(R=matrix(c(0.9, 0.5), nrow=2), r=c(0.13, 0.13)), cond_dist="Gaussian"))

  # Check the impact matrices of the regimes
  check_params(p=2, M=2, d=2, params=theta_222logistitngb_2_1, cond_dist="ind_Student", weight_function="logistic", weightfun_pars=c(2, 1),
               identification="non-Gaussianity", B_constraints=matrix(c(NA, -1, 0, 1), nrow=2))
  expect_error(check_params(p=1, M=1, d=2, params=theta_112itnge, cond_dist="ind_Student", weight_function="threshold", weightfun_pars=c(1, 1)))
  expect_error(check_params(p=1, M=2, d=2, params=theta_122logitngb_1_1e, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1),
                            cond_dist="ind_Student", identification="non-Gaussianity", B_constraints=matrix(c(1, NA, -1, 1) , nrow=2)))
  expect_error(check_params(p=2, M=2, d=2, params=theta_222logcitngb_2_1e, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1),
                            cond_dist="ind_Student", AR_constraints=C_222, identification="non-Gaussianity", matrix(c(1, NA, NA, 1), nrow=2)))
  expect_error(check_params(data=matrix(0, nrow=5, ncol=2), p=2, M=2, d=2, params=theta_222exoitngb_2_1e,  weight_function="exogenous",
                            weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_Student",
                            AR_constraints=C_222, identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, NA, NA), nrow=2)))

  check_params(p=2, M=2, d=2, params=theta_222logistiktngb_2_1, cond_dist="ind_skewed_t", weight_function="logistic", weightfun_pars=c(2, 1),
               identification="non-Gaussianity", B_constraints=matrix(c(NA, -1, 0, 1), nrow=2))
  expect_error(check_params(p=1, M=1, d=2, params=theta_112iktnge, cond_dist="ind_skewed_t", weight_function="threshold", weightfun_pars=c(1, 1)))
  expect_error(check_params(p=1, M=2, d=2, params=theta_122logiktngb_1_1e, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1),
                            cond_dist="ind_skewed_t", identification="non-Gaussianity", B_constraints=matrix(c(1, NA, -1, 1) , nrow=2)))
  expect_error(check_params(p=2, M=2, d=2, params=theta_222logciktngb_2_1e, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1),
                            cond_dist="ind_skewed_t", AR_constraints=C_222, identification="non-Gaussianity", matrix(c(1, NA, NA, 1), nrow=2)))
  expect_error(check_params(data=matrix(0, nrow=5, ncol=2), p=2, M=2, d=2, params=theta_222exoiktngb_2_1e,  weight_function="exogenous",
                            weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_skewed_t",
                            AR_constraints=C_222, identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, NA, NA), nrow=2)))

  # Check invertibility of B_t for all t
  check_params(p=2, M=2, d=2, params=theta_222logistit_2_1, weight_function="logistic", cond_dist="ind_Student", weightfun_pars=c(2, 1),
               transition_weights=cbind(c(0.1, 0.2, 0.3), c(0.9, 0.8, 0.7)))
  Omegas_222it_singt <- array(c(1:4, -(1:4)), dim=c(2, 2, 2))
  theta_222logistit_2_1_singular <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                                      vec(Omegas_222it_singt), c_and_gamma_222_2_1, dfs_222_2_1)
  expect_error(check_params(p=2, M=2, d=2, params=theta_222logistit_2_1_singular, weight_function="logistic", cond_dist="ind_Student",
                            weightfun_pars=c(2, 1), transition_weights=cbind(c(0.1, 0.2, 0.5), c(0.9, 0.8, 0.5))))
  expect_error(check_params(p=2, M=2, d=2, params=theta_222logistit_2_1_singular, weight_function="logistic", cond_dist="ind_Student",
                            weightfun_pars=c(2, 1), transition_weights=cbind(c(0.1, 0.5, 0.7), c(0.9, 0.5, 0.3))))

  check_params(p=2, M=2, d=2, params=theta_222logistikt_2_1, weight_function="logistic", cond_dist="ind_skewed_t", weightfun_pars=c(2, 1),
               transition_weights=cbind(c(0.1, 0.2, 0.3), c(0.9, 0.8, 0.7)))
  Omegas_222ikt_singt <- array(c(1:4, -(1:4)), dim=c(2, 2, 2))
  theta_222logistikt_2_1_singular <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                                      vec(Omegas_222ikt_singt), c_and_gamma_222_2_1, dfls_222_2_1)
  expect_error(check_params(p=2, M=2, d=2, params=theta_222logistikt_2_1_singular, weight_function="logistic", cond_dist="ind_skewed_t",
                            weightfun_pars=c(2, 1), transition_weights=cbind(c(0.1, 0.2, 0.5), c(0.9, 0.8, 0.5))))
  expect_error(check_params(p=2, M=2, d=2, params=theta_222logistikt_2_1_singular, weight_function="logistic", cond_dist="ind_skewed_t",
                            weightfun_pars=c(2, 1), transition_weights=cbind(c(0.1, 0.5, 0.7), c(0.9, 0.5, 0.3))))


  # Check W, lambdas, in B_constraints
  expect_error(check_params(p=1, M=2, d=2, params=theta_122relgsh_e, weight_function="relative_dens", identification="heteroskedasticity"))
  expect_error(check_params(p=2, M=2, d=2, params=theta_222logistictsh_2_1_e, weight_function="logistic", weightfun_pars=c(2, 1),
                            cond_dist="Student", identification="heteroskedasticity"))
  expect_error(check_params(p=2, M=2, d=2, params=theta_222logistictshmb_2_1_e, weight_function="logistic", weightfun_pars=c(2, 1),
                            cond_dist="Student", identification="heteroskedasticity", mean_constraints=list(1:2),
                            B_constraints=matrix(c(0.1, 0, 0, 0.3), nrow=2)))
  expect_error(check_params(p=1, M=2, d=2, params=theta_122logshwb_12_1_e, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=1),
                            identification="heteroskedasticity", weight_constraints=list(R=0, r=c(0.1, 0.2, 0.3)),
                            B_constraints=matrix(c(0.1, 0.2, 0.3, 0), nrow=2)))
  expect_error(check_params(p=1, M=2, d=3, params=theta_123expshcwb_1_1_e, weight_function="exponential", weightfun_pars=c(1, 1),
                            identification="heteroskedasticity", AR_constraints=C_123,
                            weight_constraints=list(R=matrix(c(1, 0.5), nrow=2), r=c(0, 0)),
                            B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE)))

  # Check weightpars
  expect_error(check_params(p=1, M=2, d=2, params=theta_122relg_badalphas, weight_function="relative_dens", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=3, d=2, params=theta_132relg_badalphas, weight_function="relative_dens", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=3, d=2, params=theta_132relg_badalphas2, weight_function="relative_dens", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=2, d=2, params=c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122),
                                                    0.1, 0), weight_function="logistic", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=2, d=2, params=c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122),
                                                    0.1), weight_function="logistic", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=2, d=2, params=c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122),
                                                    0.1, -0.2), weight_function="logistic", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=2, d=2, params=c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122),
                                                    0.1, 0), weight_function="exponential", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=2, d=2, params=c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122),
                                                    0.1, 0.1, 0.1), weight_function="exponential", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=2, d=2, params=c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122),
                                                    0.1, -0.2), weight_function="exponential", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=2, d=2, params=c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122),
                                                    0.1, 0, 13), weight_function="exponential", cond_dist="Student"))
  expect_error(check_params(p=1, M=3, d=2, params=c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132_stab),
                                                    vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), 1, 0.5),
                            weight_function="threshold", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=3, d=2, params=c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132_stab),
                                                    vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), 1, 0.9999),
                            weight_function="threshold", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=3, d=2, params=c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132_stab),
                                                    vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), -1, -1.9999),
                            weight_function="threshold", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=3, d=2, params=c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132_stab),
                                                    vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), 1, -1.9999),
                            weight_function="threshold", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=3, d=2, params=c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132_stab),
                                                    vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), 0, -1.9999),
                            weight_function="threshold", cond_dist="Gaussian"))

  # Checks distpars
  expect_error(check_params(p=2, M=2, d=2, params=c(theta_222logcmw_12_2, -3), weight_function="mlogit",
               weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222,
               weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5), r=c(0, 0.11, 0.12, 0.13, 0)),
               cond_dist="Student"))
  expect_error(check_params(p=2, M=2, d=2, params=c(theta_222logisticcmw_2_1, 1.999), weight_function="logistic", weightfun_pars=c(2, 1),
               mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
               cond_dist="Student"))
  expect_error(check_params(p=2, M=2, d=2, params=c(theta_222expcmw_2_1, 0), weight_function="exponential", weightfun_pars=c(2, 1),
               mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
               cond_dist="Student"))
  expect_error(check_params(p=2, M=2, d=2, params=c(theta_222thresc_1_1, 0.22), weight_function="threshold", weightfun_pars=c(1, 1),
               AR_constraints=C_222, cond_dist="Student"))

  # mlogit (nothing to check in weightpars, so just checks that the function runs)
  check_params(p=2, M=2, d=2, params=theta_222log_12_2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2),
               cond_dist="Gaussian")
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
  expect_error(check_pMd(p=1, M=3, d=2, weight_function = "logistic"))
  expect_error(check_pMd(p=1, M=3, d=2, weight_function = "exponential"))
  expect_error(check_pMd(p=1, M=1, d=2, identification="heteroskedasticity"))
  check_pMd(p=1, M=3, d=2, weight_function = "relative_dens")
  check_pMd(p=1, M=3, d=2, weight_function = "mlogit")
  check_pMd(p=1, M=3, d=2, weight_function = "threshold")
  check_pMd(p=1, M=1, d=2, identification="recursive")
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
  expect_equal(n_params(p=1, M=2, d=2, weight_function="threshold", cond_dist="Gaussian", weightfun_pars=c(1, 1)), 19)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="threshold", cond_dist="Gaussian", weightfun_pars=c(2, 2)), 27)
  expect_equal(n_params(p=1, M=3, d=2, weight_function="threshold", cond_dist="Gaussian", weightfun_pars=c(1, 1)), 29)
  expect_equal(n_params(p=1, M=3, d=2, weight_function="threshold", cond_dist="Gaussian", weightfun_pars=c(1, 1),
                        AR_constraints=C_132), 21)
  expect_equal(n_params(p=1, M=3, d=2, weight_function="threshold", cond_dist="Gaussian", weightfun_pars=c(1, 1),
                        mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))), 25)

  expect_equal(n_params(p=1, M=2, d=2, weight_function="exponential", cond_dist="Gaussian", weightfun_pars=c(1, 1)), 20)
  expect_equal(n_params(p=1, M=2, d=2, weight_function="exponential", cond_dist="Gaussian", weightfun_pars=c(2, 1)), 20)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="exponential", cond_dist="Gaussian", weightfun_pars=c(2, 1)), 28)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="exponential", cond_dist="Gaussian", weightfun_pars=c(1, 1)), 38)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="exponential", cond_dist="Gaussian", weightfun_pars=c(3, 1),
                        AR_constraints=C_123), 29)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="exponential", cond_dist="Gaussian", weightfun_pars=c(3, 1),
                        AR_constraints=C_123, mean_constraints=list(1:2)), 26)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="exponential", cond_dist="Gaussian", weightfun_pars=c(2, 1),
                        mean_constraints=list(1:2), AR_constraints=C_222,
                        weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))), 17)

  expect_equal(n_params(p=1, M=2, d=2, weight_function="logistic", cond_dist="Gaussian", weightfun_pars=c(1, 1)), 20)
  expect_equal(n_params(p=1, M=2, d=2, weight_function="logistic", cond_dist="Gaussian", weightfun_pars=c(2, 1)), 20)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="logistic", cond_dist="Gaussian", weightfun_pars=c(2, 1)), 28)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="logistic", cond_dist="Gaussian", weightfun_pars=c(1, 1)), 38)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="logistic", cond_dist="Gaussian", weightfun_pars=c(3, 1),
                        AR_constraints=C_123), 29)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="logistic", cond_dist="Gaussian", weightfun_pars=c(3, 1),
                        AR_constraints=C_123, mean_constraints=list(1:2)), 26)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="logistic", cond_dist="Gaussian", weightfun_pars=c(2, 1),
                        mean_constraints=list(1:2), AR_constraints=C_222,
                        weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))), 17)

  expect_equal(n_params(p=1, M=3, d=2, weight_function="relative_dens", cond_dist="Gaussian",
                        weight_constraints=list(R=matrix(c(0.9, 0.5), nrow=2), r=c(0.13, 0.13))), 28)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="relative_dens", cond_dist="Gaussian", mean_constraints=list(1:2),
                        AR_constraints=C_222, weight_constraints=list(R=0, r=0.6)), 16)
  expect_equal(n_params(p=1, M=2, d=2, weight_function="mlogit", cond_dist="Gaussian",
                        weightfun_pars=list(vars=1, lags=1), weight_constraints=list(R=0, r=c(0.12, 0.13))), 18)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2),
                        AR_constraints=C_222, weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5),
                                                                      r=c(0, 0.11, 0.12, 0.13, 0))), 18)

  expect_equal(n_params(p=1, M=1, d=2, weight_function="relative_dens", cond_dist="Gaussian"), 9)
  expect_equal(n_params(p=2, M=1, d=2, weight_function="relative_dens", cond_dist="Gaussian"), 13)
  expect_equal(n_params(p=3, M=1, d=2, weight_function="relative_dens", cond_dist="Gaussian"), 17)
  expect_equal(n_params(p=1, M=2, d=2, weight_function="relative_dens", cond_dist="Gaussian"), 19)
  expect_equal(n_params(p=1, M=3, d=2, weight_function="relative_dens", cond_dist="Gaussian"), 29)
  expect_equal(n_params(p=1, M=1, d=3, weight_function="relative_dens", cond_dist="Gaussian"), 18)
  expect_equal(n_params(p=2, M=1, d=3, weight_function="relative_dens", cond_dist="Gaussian"), 27)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="relative_dens", cond_dist="Gaussian"), 37)

  expect_equal(n_params(p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="Gaussian"), 20)
  expect_equal(n_params(p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=1), cond_dist="Gaussian"), 21)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="Gaussian"), 28)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2), cond_dist="Gaussian"), 31)

  expect_equal(n_params(p=1, M=1, d=2, weight_function="relative_dens", cond_dist="Gaussian", AR_constraints=C_112), 7)
  expect_equal(n_params(p=1, M=2, d=2, weight_function="relative_dens", cond_dist="Gaussian", AR_constraints=C_122), 15)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="relative_dens", cond_dist="Gaussian", AR_constraints=C_222_2), 15)
  expect_equal(n_params(p=1, M=3, d=2, weight_function="relative_dens", cond_dist="Gaussian", AR_constraints=C_132), 21)
  expect_equal(n_params(p=1, M=1, d=3, weight_function="relative_dens", cond_dist="Gaussian", AR_constraints=C_113), 12)
  expect_equal(n_params(p=2, M=1, d=3, weight_function="relative_dens", cond_dist="Gaussian", AR_constraints=C_213), 10)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="relative_dens", cond_dist="Gaussian", AR_constraints=C_123), 28)
  expect_equal(n_params(p=1, M=1, d=2, weight_function="relative_dens", cond_dist="Gaussian", AR_constraints=C_112,
                        mean_constraints=list(1)), 7)
  expect_equal(n_params(p=1, M=2, d=2, weight_function="relative_dens", cond_dist="Gaussian", mean_constraints=list(1:2)), 17)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="relative_dens", cond_dist="Gaussian", AR_constraints=C_222,
                        mean_constraints=list(1:2)), 17)
  expect_equal(n_params(p=1, M=3, d=2, weight_function="relative_dens", cond_dist="Gaussian", AR_constraints=C_132,
                        mean_constraints=list(1:3)), 17)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="relative_dens", cond_dist="Gaussian", mean_constraints=list(1:2)), 34)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="relative_dens", cond_dist="Gaussian", AR_constraints=C_123,
                        mean_constraints=list(1:2)), 25)

  expect_equal(n_params(p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="Gaussian",
                        AR_constraints=C_222), 20)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2), cond_dist="Gaussian",
                        AR_constraints=C_222), 23)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="mlogit", weightfun_pars=list(vars=2:3, lags=1), cond_dist="Gaussian",
                        AR_constraints=C_123), 30)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="mlogit", weightfun_pars=list(vars=1:3, lags=1), cond_dist="Gaussian",
                        AR_constraints=C_123, mean_constraints=list(1:2)), 28)

  expect_equal(n_params(p=1, M=2, d=3, weight_function="exogenous", weightfun_pars=cbind(c(0.4, 0.2, 0.9), c(0.6, 0.8, 0.1))), 36)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(0.4, 0.2, 0.9), c(0.6, 0.8, 0.1))), 26)
  expect_equal(n_params(p=2, M=3, d=2, weight_function="exogenous", weightfun_pars=cbind(c(0.4, 0, 0.5), c(0.3, 1, 0.2), c(0.3, 0, 0.3))),
               39)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.6, 0.3), c(0, 0.4, 07)), AR_constraints=C_222), 18)
  expect_equal(n_params(p=1, M=3, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.6, 0.0), c(0, 0.3, 07), c(0, 0.1, 0.3)),
                        AR_constraints=C_132), 19)

  # Student
  expect_equal(n_params(p=1, M=3, d=2, weight_function="threshold", cond_dist="Student", weightfun_pars=c(1, 1),
                        mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))), 26)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="Student",
                        AR_constraints=C_222), 21)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="logistic", cond_dist="Student", weightfun_pars=c(3, 1),
                        AR_constraints=C_123), 30)
  expect_equal(n_params(p=1, M=2, d=2, weight_function="exponential", cond_dist="Student", weightfun_pars=c(1, 1)), 21)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="Student",
                        AR_constraints=C_222), 19)

  # Ind student
  expect_equal(n_params(p=1, M=1, d=2, cond_dist="ind_Student", weight_function="threshold", weightfun_pars=c(1, 1)), 12)
  expect_equal(n_params(p=2, M=2, d=2, cond_dist="ind_Student", weight_function="logistic", weightfun_pars=c(2, 1)), 32)
  expect_equal(n_params(p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="ind_Student"), 24)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="ind_Student"), 47)
  expect_equal(n_params(p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student"), 34)
  expect_equal(n_params(p=2, M=2, d=2,  weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="ind_Student",
                        mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))), 21)
  expect_equal(n_params(p=1, M=3, d=2,  weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student",
                        mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))), 30)
  expect_equal(n_params(p=1, M=2, d=3,weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="ind_Student",
                        mean_constraints=list(1:2), AR_constraints=C_123), 35)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_Student",
                        AR_constraints=C_222), 24)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_Student",
                        AR_constraints=C_222), 22)

  # Ind skewed t
  expect_equal(n_params(p=1, M=1, d=2, cond_dist="ind_skewed_t", weight_function="threshold", weightfun_pars=c(1, 1)), 14)
  expect_equal(n_params(p=2, M=2, d=2, cond_dist="ind_skewed_t", weight_function="logistic", weightfun_pars=c(2, 1)), 34)
  expect_equal(n_params(p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="ind_skewed_t"), 26)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t"), 50)
  expect_equal(n_params(p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t"), 36)
  expect_equal(n_params(p=2, M=2, d=2,  weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="ind_skewed_t",
                        mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))), 23)
  expect_equal(n_params(p=1, M=3, d=2,  weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t",
                        mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))), 32)
  expect_equal(n_params(p=1, M=2, d=3,weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="ind_skewed_t",
                        mean_constraints=list(1:2), AR_constraints=C_123), 38)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_skewed_t",
                        AR_constraints=C_222), 26)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_skewed_t",
                        AR_constraints=C_222), 24)

  ## Structural
  expect_equal(n_params(p=1, M=3, d=2, weight_function="threshold", cond_dist="Student", identification="recursive", weightfun_pars=c(1, 1),
                        mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))), 26)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="logistic", cond_dist="Student", identification="recursive", weightfun_pars=c(3, 1),
                        AR_constraints=C_123), 30)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="mlogit", weightfun_pars=list(vars=1:3, lags=1), cond_dist="Gaussian",
                        identification="recursive", AR_constraints=C_123, mean_constraints=list(1:2)), 28)

  expect_equal(n_params(p=1, M=2, d=2, weight_function="relative_dens", identification="heteroskedasticity"), 19)
  expect_equal(n_params(p=1, M=3, d=2, weight_function="relative_dens", identification="heteroskedasticity"), 28)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student",
                        identification="heteroskedasticity"), 29)
  expect_equal(n_params(p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=1),
                        identification="heteroskedasticity"), 21)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), identification="heteroskedasticity"), 38)
  expect_equal(n_params(p=2, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="Student",
                        identification="heteroskedasticity"), 41)
  expect_equal(n_params(p=1, M=2, d=2, weight_function="relative_dens", identification="heteroskedasticity", AR_constraints=C_122), 15)
  expect_equal(n_params(p=1, M=3, d=2, weight_function="relative_dens", identification="heteroskedasticity",
                        B_constraints=matrix(c(0.1, NA, 0.3, 0), nrow=2)), 27)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student",
                        identification="heteroskedasticity", mean_constraints=list(1:2), B_constraints=matrix(c(0.1, 0, 0, 0.3), nrow=2)), 25)
  expect_equal(n_params(p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=1), identification="heteroskedasticity",
                        weight_constraints=list(R=0, r=c(0.1, 0.2, 0.3)), B_constraints=matrix(c(0.1, 0.2, 0.3, 0), nrow=2)), 17)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), identification="heteroskedasticity",
                        AR_constraints=C_123, weight_constraints=list(R=matrix(c(1, 0.5), nrow=2), r=c(0, 0)),
                        B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE)), 26)
  expect_equal(n_params(p=2, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="Student",
                        identification="heteroskedasticity", B_constraints=matrix(c(0.1, 0.2, -0.3, 0), nrow=2)), 40)

  # Ind Student struct
  expect_equal(n_params(p=2, M=2, d=2, cond_dist="ind_Student", weight_function="logistic", weightfun_pars=c(2, 1),
                        identification="non-Gaussianity"), 32)
  expect_equal(n_params(p=2, M=2, d=2, cond_dist="ind_Student", weight_function="logistic", weightfun_pars=c(2, 1),
                        identification="non-Gaussianity", B_constraints=matrix(c(NA, -1, 0, 1), nrow=2)), 30)
  expect_equal(n_params(p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="ind_Student",
                        identification="non-Gaussianity", B_constraints=matrix(c(1, NA, -1, 1) , nrow=2)), 24)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="ind_Student",
                        identification="non-Gaussianity", B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)), 43)
  expect_equal(n_params(p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student",
                        identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2)), 31)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="ind_Student",
                        mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
                        identification="non-Gaussianity", matrix(c(NA, -1, 0, 1), nrow=2)), 19)
  expect_equal(n_params(p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student",
                        mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2)), identification="non-Gaussianity",
                        B_constraints=matrix(c(1, 0, NA, 1), nrow=2)), 27)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_Student",
                        AR_constraints=C_222, identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)), 22)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_Student",
                        AR_constraints=C_222, identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)), 20)

  expect_equal(n_params(p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_Student",
                        AR_constraints=C_222, identification="non-Gaussianity", B_constraints=matrix(c(1, NA, NA, 1), nrow=2)), 24)

  # Ind skewed t struct
  expect_equal(n_params(p=2, M=2, d=2, cond_dist="ind_skewed_t", weight_function="logistic", weightfun_pars=c(2, 1),
                        identification="non-Gaussianity"), 34)
  expect_equal(n_params(p=2, M=2, d=2, cond_dist="ind_skewed_t", weight_function="logistic", weightfun_pars=c(2, 1),
                        identification="non-Gaussianity", B_constraints=matrix(c(NA, -1, 0, 1), nrow=2)), 32)
  expect_equal(n_params(p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="ind_skewed_t",
                        identification="non-Gaussianity", B_constraints=matrix(c(1, NA, -1, 1) , nrow=2)), 26)
  expect_equal(n_params(p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t",
                        identification="non-Gaussianity", B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)), 46)
  expect_equal(n_params(p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t",
                        identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2)), 33)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="ind_skewed_t",
                        mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
                        identification="non-Gaussianity", matrix(c(NA, -1, 0, 1), nrow=2)), 21)
  expect_equal(n_params(p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t",
                        mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2)), identification="non-Gaussianity",
                        B_constraints=matrix(c(1, 0, NA, 1), nrow=2)), 29)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_skewed_t",
                        AR_constraints=C_222, identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)), 24)
  expect_equal(n_params(p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(c(1, 0.9, 0.8), c(0, 0.1, 0.2)), cond_dist="ind_skewed_t",
                        AR_constraints=C_222, identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)), 22)

  expect_equal(n_params(p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_skewed_t",
                        AR_constraints=C_222, identification="non-Gaussianity", B_constraints=matrix(c(1, NA, NA, 1), nrow=2)), 26)
})


test_that("check_constraints works correctly", {
  expect_error(check_constraints(p=1, M=3, d=2, weight_function="relative_dens",
                                 weight_constraints=list(R=c(0.9, 0.5), r=c(0.13, 0.13))))
  expect_error(check_constraints(p=1, M=3, d=2, weight_function="relative_dens",
                                 weight_constraints=list(R=matrix(c(0.9, 0.5, 1), nrow=3), r=c(0.13, 0.13))))
  expect_error(check_constraints(p=1, M=3, d=2, weight_function="relative_dens",
                                 weight_constraints=list(R=matrix(c(0.9, 0.5, 1, 2, 3, 4), nrow=2), r=c(0.13, 0.13))))
  expect_error(check_constraints(p=1, M=3, d=2, weight_function="relative_dens",
                                 weight_constraints=list(R=matrix(c(0.9, 0.5), nrow=2), r=0.13)))
  expect_error(check_constraints(p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1),
                                 weight_constraints=list(R=0, r=c(0.12, 0.13, 0.14))))
  expect_error(check_constraints(p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1),
                                 weight_constraints=list(R=0)))
  expect_error(check_constraints(p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2),
                                 weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=2), r=c(0, 0.11, 0.12, 0.13, 0))))
  expect_error(check_constraints(p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2),
                                 weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5), r=c(1, 0, 0.11, 0.12, 0.13, 0))))
  expect_error(check_constraints(p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2),
                                 weight_constraints=list(r=c(0, 0.11, 0.12, 0.13, 0))))
  expect_error(check_constraints(p=1, M=2, d=2, weight_function="logistic", weightfun_pars=c(1, 1),
                                 weight_constraints=list(R=0, r=c(0.13))))
  expect_error(check_constraints(p=1, M=2, d=2, weight_function="logistic", weightfun_pars=c(1, 1),
                                 weight_constraints=list(R=0, r=c(0.13, 0.1, 0.1))))
  expect_error(check_constraints(p=1, M=2, d=2, weight_function="logistic", weightfun_pars=c(1, 1),
                                 weight_constraints=list(R=matrix(1:6, nrow=3), r=c(0.13, 0.1))))
  expect_error(check_constraints(p=1, M=2, d=2, weight_function="logistic", weightfun_pars=c(1, 1),
                                 weight_constraints=list(R=matrix(1:6, nrow=2), r=c(0.13, 0.1))))
  expect_error(check_constraints(p=1, M=2, d=2, weight_function="exponential", weightfun_pars=c(1, 1),
                                 weight_constraints=list(R=0, r=c(0.13))))
  expect_error(check_constraints(p=1, M=2, d=2, weight_function="exponential", weightfun_pars=c(1, 1),
                                 weight_constraints=list(R=0, r=c(0.13, 0.1, 0.1))))
  expect_error(check_constraints(p=1, M=2, d=2, weight_function="exponential", weightfun_pars=c(1, 1),
                                 weight_constraints=list(R=matrix(1:6, nrow=3), r=c(0.13, 0.1))))
  expect_error(check_constraints(p=1, M=2, d=2, weight_function="exponential", weightfun_pars=c(1, 1),
                                 weight_constraints=list(R=matrix(1:6, nrow=2), r=c(0.13, 0.1))))
  expect_error(check_constraints(p=1, M=2, d=2, weight_function="threshold", weightfun_pars=c(1, 1),
                                 weight_constraints=list(R=0, r=c(0.13, 0.1))))
  expect_error(check_constraints(p=2, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1),
                                 weight_constraints=list(R=matrix(1:6, nrow=3), r=c(0.1, 0.13))))
  expect_error(check_constraints(p=2, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1),
                                 weight_constraints=list(R=matrix(1:6, nrow=2), r=c(0.1, 0.13))))
  expect_error(check_constraints(p=2, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1),
                                 weight_constraints=list(R=matrix(1:4, nrow=2), r=c(0.13))))
  expect_warning(check_constraints(data=matrix(NA, nrow=3), p=2, M=2, d=2, weight_function="exogenous",
                                   weightfun_pars=cbind(0, 1), weight_constraints=list(R=0, r=c(0.13, 0.1))))
  expect_warning(check_constraints(p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(0, 1),
                                   weight_constraints=list(R=0, r=c(0.13, 0.1))))
  check_constraints(p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(0, 1))
  check_constraints(data=NULL, p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=cbind(0, 1))

  expect_error(check_constraints(p=1, M=1, d=2, AR_constraints=cbind(C_112, C_112)))
  expect_error(check_constraints(p=1, M=1, d=2, AR_constraints=rbind(C_112, C_112)))
  expect_error(check_constraints(p=1, M=2, d=2, AR_constraints=cbind(C_122[,1:3], rep(0, 8))))
  expect_error(check_constraints(p=1, M=2, d=2, AR_constraints=as.data.frame(C_122)))
  expect_error(check_constraints(p=2, M=2, d=2, AR_constraints=C_222[-1,]))
  expect_error(check_constraints(p=2, M=2, d=2, AR_constraints=C_222_2[-16,]))
  expect_error(check_constraints(p=1, M=2, d=3, AR_constraints=cbind(C_123[,1:8], C_123[,8])))
  expect_error(check_constraints(p=2, M=1, d=3, AR_constraints=as.vector(C_213)))

  expect_error(check_constraints(p=2, M=2, d=2, mean_constraints=1:2))
  expect_error(check_constraints(p=2, M=2, d=2, mean_constraints=1:2))
  expect_error(check_constraints(p=1, M=1, d=2, mean_constraints=1:2))
  expect_error(check_constraints(p=2, M=2, d=2, mean_constraints=list(1)))
  expect_error(check_constraints(p=2, M=2, d=2, mean_constraints=list(1, "2")))
  expect_error(check_constraints(p=2, M=3, d=2, mean_constraints=list(2:3)))
  expect_error(check_constraints(p=3, M=2, d=2, mean_constraints=list(1, 1:2)))
  expect_error(check_constraints(p=3, M=2, d=2, mean_constraints=list(1, 1:2)))
  expect_error(check_constraints(p=1, M=1, d=2, mean_constraints=list()))
  expect_error(check_constraints(p=1, M=1, d=2, mean_constraints=list(1:2)))

  # Check B_constraints
  expect_error(check_constraints(p=1, M=2, d=2, identification="heteroskedasticity", B_constraints=matrix(1:9, nrow=3)))
  expect_error(check_constraints(p=1, M=2, d=2, identification="heteroskedasticity", B_constraints=matrix(c(0, 0, 1, NA), nrow=2)))
  expect_error(check_constraints(p=2, M=2, d=3, identification="heteroskedasticity",
                                 B_constraints=matrix(c(1, 0, NA, NA, 0, 1, 1, 0, NA), nrow=3)))
  expect_error(check_constraints(p=1, M=2, d=3, identification="heteroskedasticity", B_constraints=matrix(1:4, nrow=2)))
  expect_error(check_constraints(p=2, M=2, d=2, identification="heteroskedasticity", B_constraints=1:4))

  expect_error(check_constraints(p=1, M=2, d=2, identification="non-Gaussianity", B_constraints=matrix(1:9, nrow=3)))
  expect_error(check_constraints(p=1, M=2, d=2, identification="non-Gaussianity", B_constraints=matrix(c(0, 0, 1, NA), nrow=2)))
  expect_error(check_constraints(p=2, M=2, d=3, identification="non-Gaussianity",
                                 B_constraints=matrix(c(1, 0, NA, NA, 0, 1, 1, 0, NA), nrow=3)))
  expect_error(check_constraints(p=1, M=2, d=3, identification="non-Gaussianity", B_constraints=matrix(1:4, nrow=2)))
  expect_error(check_constraints(p=2, M=2, d=2, identification="non-Gaussianity", B_constraints=1:4))
})


weightfun_pars1 <- list(3, 4)
weightfun_pars2 <- list(tmp1=4, tmp2=5)
weightfun_pars3 <- list(vars=2:3, lags=2:3)
weightfun_pars4 <- list(3:4, lags=2)
weightfun_pars5 <- list(vars=1:4, lags=3)
weightfun_pars6 <- list(vars=c(3, 1), 1)


test_that("check_weightfun_pars works correctly", {
  # relative_dens
  expect_equal(check_weightfun_pars(p=1, d=2, weight_function="relative_dens", weightfun_pars="testobj"), NULL)
  expect_error(check_weightfun_pars(p=2, d=2, weight_function="relative_dens", cond_dist="Student"))

  # logistic
  expect_error(check_weightfun_pars(p=5, d=4, weight_function="logistic", weightfun_pars=list(1, 1)))
  expect_error(check_weightfun_pars(p=2, d=2, weight_function="logistic", weightfun_pars=c(1, 1, 1)))
  expect_error(check_weightfun_pars(p=2, d=2, weight_function="logistic", weightfun_pars=c(3, 1)))
  expect_error(check_weightfun_pars(p=2, d=2, weight_function="logistic", weightfun_pars=c(1, 3)))

  # mlogit
  expect_error(check_weightfun_pars(p=5, d=2, weight_function="mlogit", weightfun_pars=weightfun_pars1))
  expect_error(check_weightfun_pars(p=3, d=3, weight_function="mlogit", weightfun_pars=weightfun_pars1))
  expect_equal(check_weightfun_pars(p=4, d=3, weight_function="mlogit", weightfun_pars=weightfun_pars1),
               list(vars=3, lags=4))
  expect_error(check_weightfun_pars(p=4, d=5, weight_function="mlogit", weightfun_pars=weightfun_pars2))
  expect_error(check_weightfun_pars(p=10, d=3, weight_function="mlogit", weightfun_pars=weightfun_pars2))
  expect_equal(check_weightfun_pars(p=5, d=4, weight_function="mlogit", weightfun_pars=weightfun_pars2),
               list(vars=4, lags=5))
  expect_error(check_weightfun_pars(p=3, d=2, weight_function="mlogit", weightfun_pars=weightfun_pars3))
  expect_error(check_weightfun_pars(p=2, d=3, weight_function="mlogit", weightfun_pars=weightfun_pars3))
  expect_error(check_weightfun_pars(p=5, d=3, weight_function="mlogit", weightfun_pars=weightfun_pars4))
  expect_error(check_weightfun_pars(p=1, d=5, weight_function="mlogit", weightfun_pars=weightfun_pars4))
  expect_equal(check_weightfun_pars(p=2, d=4, weight_function="mlogit", weightfun_pars=weightfun_pars4),
               list(vars=3:4, lags=2))
  expect_error(check_weightfun_pars(p=2, d=5, weight_function="mlogit", weightfun_pars=weightfun_pars5))
  expect_error(check_weightfun_pars(p=3, d=2, weight_function="mlogit", weightfun_pars=weightfun_pars5))
  expect_equal(check_weightfun_pars(p=3, d=4, weight_function="mlogit", weightfun_pars=weightfun_pars5),
               weightfun_pars5)
  expect_equal(check_weightfun_pars(p=2, d=3, weight_function="mlogit", weightfun_pars=weightfun_pars6),
               list(vars=c(1, 3), lags=1))

  # exponential
  expect_error(check_weightfun_pars(p=5, d=4, weight_function="exponential", weightfun_pars=list(1, 1)))
  expect_error(check_weightfun_pars(p=2, d=2, weight_function="exponential", weightfun_pars=c(1, 1, 1)))
  expect_error(check_weightfun_pars(p=2, d=2, weight_function="exponential", weightfun_pars=c(3, 1)))
  expect_error(check_weightfun_pars(p=2, d=2, weight_function="exponential", weightfun_pars=c(1, 3)))

  # threshold
  expect_error(check_weightfun_pars(p=5, d=4, weight_function="threshold", weightfun_pars=list(1, 1)))
  expect_error(check_weightfun_pars(p=2, d=2, weight_function="threshold", weightfun_pars=c(1, 1, 1)))
  expect_error(check_weightfun_pars(p=2, d=2, weight_function="threshold", weightfun_pars=c(3, 1)))
  expect_error(check_weightfun_pars(p=2, d=2, weight_function="threshold", weightfun_pars=c(1, 3)))

  # Exogenous
  expect_equal(check_weightfun_pars(data=matrix(NA, nrow=4, ncol=3), p=1, M=2, d=3, weight_function="exogenous",
                                    weightfun_pars=cbind(c(0.4, 0.2, 0.9), c(0.6, 0.8, 0.1))),
               cbind(c(0.4, 0.2, 0.9), c(0.6, 0.8, 0.1)))
  expect_error(check_weightfun_pars(data=matrix(NA, nrow=4, ncol=3), p=1, M=2, d=3, weight_function="exogenous",
                                    weightfun_pars=cbind(c(0.4, 0.2), c(0.6, 0.8))))
  expect_error(check_weightfun_pars(data=matrix(NA, nrow=4, ncol=3), p=1, M=2, d=3, weight_function="exogenous",
                                    weightfun_pars=cbind(c(0.4, 0.2, -0.9), c(0.6, 0.8, 0.1))))
  expect_equal(check_weightfun_pars(data=matrix(NA, nrow=4), p=2, M=2, d=2, weight_function="exogenous",
                                    weightfun_pars=cbind(c(0.4, 0.2), c(0.6, 0.8))), cbind(c(0.4, 0.2), c(0.6, 0.8)))
  expect_equal(check_weightfun_pars(data=matrix(NA, nrow=4), p=2, M=2, d=2, weight_function="exogenous",
                                    weightfun_pars=cbind(c(0, 0), c(1, 1))), cbind(c(0, 0), c(1, 1)))
  expect_error(check_weightfun_pars(data=matrix(NA, nrow=4), p=2, M=2, d=2, weight_function="exogenous",
                                    weightfun_pars=cbind(c("1", "a"), c(1, 1))))
  expect_warning(check_weightfun_pars(data=matrix(NA, nrow=4), p=2, M=2, d=2, weight_function="exogenous",
                                       weightfun_pars=cbind(c(0, 1), c(1, 1))))
  expect_equal(check_weightfun_pars(data=matrix(NA, nrow=2, ncol=3), p=1, M=3, d=2, weight_function="exogenous",
                                    weightfun_pars=cbind(c(0.4), c(0.6), c(0))), cbind(c(0.4), c(0.6), c(0)))
  expect_error(check_weightfun_pars(data=matrix(NA, nrow=2, ncol=3), p=1, M=3, d=2, weight_function="exogenous",
                                    weightfun_pars=cbind(c(0.4), c(0.6))))
  expect_error(check_weightfun_pars(data=matrix(NA, nrow=2, ncol=3), p=1, M=3, d=2, weight_function="exogenous",
                                    weightfun_pars=cbind(c(0.4), c(0.6), 0, 0)))
  expect_equal(check_weightfun_pars(data=matrix(NA, nrow=4, ncol=3), p=2, M=3, d=2, weight_function="exogenous",
                                    weightfun_pars=cbind(c(0.4, 0.5), c(0.6, 0.4), c(0, 0.1))),
               cbind(c(0.4, 0.5), c(0.6, 0.4), c(0, 0.1)))
  expect_warning(check_weightfun_pars(data=matrix(NA, nrow=4, ncol=3), p=2, M=3, d=2, weight_function="exogenous",
                                    weightfun_pars=cbind(c(0.4, 0.5), c(0.6, 0.4), c(0, 0))))
  expect_error(check_weightfun_pars(data=matrix(NA, nrow=4, ncol=3), p=2, M=3, d=2, weight_function="exogenous",
                                      weightfun_pars=cbind(c(0.4, 0.5, 0.1), c(0.6, 0.4, 0.9), c(0, 0.1, 0))))
  expect_error(check_weightfun_pars(data=matrix(NA, nrow=4, ncol=3), p=2, M=3, d=2, weight_function="exogenous",
                                    weightfun_pars=cbind(c(0.4, 0.5), c(-0.6, 0.4), c(0, 0.1))))

  # Exogenous without data
  expect_equal(check_weightfun_pars(data=NULL, p=2, M=3, d=2, weight_function="exogenous",
                                    weightfun_pars=cbind(c(0.4, 0.5), c(0.6, 0.4), c(0, 0.1))),
               cbind(c(0.4, 0.5), c(0.6, 0.4), c(0, 0.1)))
  expect_equal(check_weightfun_pars(p=2, M=3, d=2, weight_function="exogenous",
                                    weightfun_pars=cbind(c(0.4, 0.5), c(0.6, 0.4), c(0, 0.1))),
               cbind(c(0.4, 0.5), c(0.6, 0.4), c(0, 0.1)))
  expect_error(check_weightfun_pars(data=NULL, p=2, M=3, d=2, weight_function="exogenous",
                                    weightfun_pars=cbind(c(0.4, 0.5), c(0.6, 0.4))))
})


# p=1, M=1, d=2
theta_112relg <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103, 0.601786, -0.002945, 0.067224)
mod112relg <- STVAR(data=gdpdef, p=1, M=1, params=theta_112relg, weight_function="relative_dens")

test_that("check_stvar works correctly", {
   check_stvar(mod112relg)
   expect_error(check_stvar(theta_112relg))
})


test_that("check_exoweights works correctly", {
  expect_error(check_exoweights(M=2, exo_weights=NULL, how_many_rows=2, name_of_row_number="tmp"))
  expect_error(check_exoweights(M=2, exo_weights=cbind(c(0.5, 0.6), c(0.5, 0.4)), how_many_rows=3, name_of_row_number="tmp"))
  expect_error(check_exoweights(M=3, exo_weights=cbind(c(0.5, 0.6), c(0.5, 0.4)), how_many_rows=2, name_of_row_number="tmp"))
  expect_error(check_exoweights(M=2, exo_weights=cbind(c(0.5, 0.6), c(0.5, 0.401)), how_many_rows=2, name_of_row_number="tmp"))
  expect_error(check_exoweights(M=2, exo_weights=c(0.5, 0.6, 0.5, 0.4), how_many_rows=2, name_of_row_number="tmp"))
})
