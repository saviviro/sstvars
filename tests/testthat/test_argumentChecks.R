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

  # mlogit (nothing to check in weightpars, so just checks that the function runs)
  expect_true(in_paramspace(p=2, M=2, d=2, weight_function="mlogit", cond_dist="Gaussian",
                            all_boldA=boldA_222log_12_2, all_Omegas=Omegas_222log_12_2, weightpars=weightpars_222log_12_2))
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


  # Check weightpars
  expect_error(check_params(p=1, M=2, d=2, params=theta_122relg_badalphas, weight_function="relative_dens", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=3, d=2, params=theta_132relg_badalphas, weight_function="relative_dens", cond_dist="Gaussian"))
  expect_error(check_params(p=1, M=3, d=2, params=theta_132relg_badalphas2, weight_function="relative_dens", cond_dist="Gaussian"))

  # Checks df
  # TO BE FILLED IN

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
})
