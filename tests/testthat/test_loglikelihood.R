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

# p=3, M=1, d=3, allow_unstab=TRUE, samone
theta_313relg_unstab <- c(0.14652, 0.07905, -0.06877, 2.85178, -0.0212, 0.15671, -0.05778, 0.49156, 0.12683,
                          0.04203, 0.06951, 1.19613, 2.08181, 0.02087, -0.03947, 0.1369, 0.31143, 0.52774,
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

# # p=1, M=2, d=3, weightfun_pars=c(1, 1)
# c_and_gamma_123_1_1 <- c(0.1, 0.5)
# theta_123logistic_1_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vech(Omega1_123),
#                            vech(Omega2_123), c_and_gamma_123_1_1)
#
# # p=1, M=2, d=3, weightfun_pars=c(3, 1)
# c_and_gamma_123_3_1 <- c(0.1, 0.4)
# theta_123logistic_3_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vech(Omega1_123),
#                            vech(Omega2_123), c_and_gamma_123_3_1)


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
c_and_gamma_123_1_1 <- c(0.5, 1)
theta_123exp_1_1 <-  c(theta_123log_noweightpars, c_and_gamma_123_1_1)

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
r1_232_1_1 <- -0.01; r2_232_1_1 <- 1.01
theta_232thres_1_1 <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                        vec(A31_232), vec(A32_232), vech(Omega1_232), vech(Omega2_232), vech(Omega3_232),
                        r1_232_1_1, r2_232_1_1)


# p=1, M=2, d=3, weight_function="threshold", weightfun_pars=c(2, 1)
r1_123_2_1 <- 1
theta_123thres_2_1 <- c(theta_123log_noweightpars, r1_123_2_1)


## weight_function == "exogenous"

# p=2, M=2, d=2,  weight_function="exogenous", weightfun_pars=weightfun_pars222
set.seed(1); weightfun_pars222 <- abs(matrix(rnorm(2*(nrow(gdpdef)-2)), nrow=nrow(gdpdef)-2, ncol=2))
weightfun_pars222 <- weightfun_pars222/rowSums(weightfun_pars222)
theta_222exo <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222), vech(Omega1_222), vech(Omega2_222))

# p=1, M=3, d=2, weight_function="exogenous", weightfun_pars=weightfun_pars132
set.seed(2); weightfun_pars132 <- abs(matrix(rnorm(3*(nrow(gdpdef)-1)), nrow=nrow(gdpdef)-1, ncol=3))
weightfun_pars132 <- weightfun_pars132/rowSums(weightfun_pars132)
theta_132exo <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                  vech(Omega2_132), vech(Omega3_132))

# p=1, M=2, d=3, weight_function="exogenous", weightfun_pars=weightfun_pars123
phi10_123 <- c(0.5, -0.2, 1.5)
phi20_123 <- c(0.9, 2.0, -0.7)
A11_123 <- matrix(c(0.2, 0.1, -0.1, 0.3, 0.4, -0.2, 0.1, 0.2, 0.3), nrow=3)
A21_123 <- matrix(c(0.3, -0.1, 0.1, 0.4, 0.2, -0.3, -0.1, 0.2, 0.5), nrow=3)
Omega1_123 <- matrix(c(c(1, 0.2, 0.3, 0.2, 2, 0.4, 0.3, 0.4, 3)), nrow=3)
Omega2_123 <- matrix(c(c(1.1, -0.2, -0.3, -0.2, 2.2, -0.4, -0.3, -0.4, 3.3)), nrow=3)
set.seed(3); weightfun_pars123 <- abs(matrix(rnorm(2*(nrow(usamone)-1)), nrow=nrow(usamone)-1, ncol=2))
weightfun_pars123 <- weightfun_pars123/rowSums(weightfun_pars123)
theta_123exo <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vech(Omega1_123), vech(Omega2_123))

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


# weight_function == "mlogit"

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

# p=1, M=1, p=2, mean_constraints=list(1), parametrization="mean"
theta_112relgm <- theta_112relg

# p=1, M=2, d=2, mean_constraints=list(1:2), parametrization="mean"
mu_122relgm <- c(0.5342688, 0.7382623)
theta_122relgm <- c(mu_122relgm, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)
theta_122relgm_expanded <- c(mu_122relgm, mu_122relgm, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), alpha1_122)

# p=2, M=2, d=2, mean_constraints=list(1:2), C_222, parametrization="mean"
mu_222relgcm <- c(0.7209658, 0.8108580)
theta_222relgcm <- c(mu_222relgcm, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), alpha1_222)
theta_222relgcm_expanded <- c(mu_222relgcm, mu_222relgcm, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                              vech(Omega1_222), vech(Omega2_222), alpha1_222)

# weightfunction == "mlogit"

# p=1, M=2, d=2, weigthfun_pars=list(vars=1, lags=1), mean_constraints=list(1:2), parametrization="mean"
theta_122logm_1_1 <- c(phi10_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), gamma1_122_1_1)
theta_122logm_1_1_expanded <- c(phi10_122, phi10_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), gamma1_122_1_1)

# p=2, M=2, d=2, weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), C_222, parametrization="mean"
theta_222logcm_12_2 <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), gamma1_222_12_2)
theta_222logcm_12_2_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                  vech(Omega1_222), vech(Omega2_222), gamma1_222_12_2)


## Models with weight_constraints

# p=1, M=3, d=2, weight_function="relative_dens", weight_constraints=list(R=matrix(c(0.9, 0.5), nrow=2), r=c(0.13, 0.13))
# xi_132relgw <- 0.4
# theta_132relgw <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
#                     vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), xi_132relgw)
# theta_132relgw_expanded <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
#                              vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), matrix(c(0.9, 0.5), nrow=2)%*%xi_132relgw + 0.13)

# p=2, M=2, d=2, weight_function="relative_dens", mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=0, r=0.6), parametrization="mean"
theta_222relgcmw <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222)) # No weight param since replaced with r
theta_222relgcmw_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                               vech(Omega1_222), vech(Omega2_222), 0.6)


# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), weight_constraints=list(R=0, r=c(0.12, 0.13))
theta_122logw_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122))
theta_122logw_1_1_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), c(0.12, 0.13))

# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5), r=c(0, 0.11, 0.12, 0.13, 0)), parametrization="mean"
xi_222logcmw_12_2 <- c(0.22, 0.33)
theta_222logcmw_12_2 <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), xi_222logcmw_12_2 )
theta_222logcmw_12_2_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                   vech(Omega1_222), vech(Omega2_222), c(0.22, 0.11, 0.12, 0.13, 0.33))

# p=1, M=2, d=2, weight_function="logistic", weightfun_pars=c(1, 1), weight_constraints=list(R=0, r=c(0.02, 0.13))
theta_122logisticw_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122))
theta_122logisticw_1_1_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), c(0.02, 0.13))

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean"
xi_222logisticcmw_2_1 <- c(0.33)
theta_222logisticcmw_2_1 <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), xi_222logisticcmw_2_1)
theta_222logisticcmw_2_1_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                       vech(Omega1_222), vech(Omega2_222), c(0.01, 0.33))


## weight_function == "exponential"

# p=1, M=2, d=2, weight_function="exponential", weightfun_pars=c(1, 1), weight_constraints=list(R=0, r=c(0.02, 0.13))
theta_122expw_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122))
theta_122expw_1_1_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vech(Omega1_122), vech(Omega2_122), c(0.02, 0.13))

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean"
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

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), mean_constraints=list(1, 2:3),
# weight_constraints=list(R=0, r=c(0, 1.2)), parametrization="mean"
theta_132thresmw_1_1 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                          vech(Omega2_132), vech(Omega3_132))
theta_132thresmw_1_1_expanded <- c(phi10_132, phi20_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                                   vech(Omega2_132), vech(Omega3_132), 0, 1.2)

## Exogenous

# p=1, M=2, d=3, weight_function="exogenous", weightfun_pars=weightfun_pars123, mean_constraints=list(1:2), AR_constraints=C_123,
# parametrization="mean"
theta_123exoc <- c(phi10_123, vec(A11_123), vech(Omega1_123), vech(Omega2_123))
theta_123exoc_expanded <- c(phi10_123, phi20_123, vec(A11_123), vec(A11_123), vech(Omega1_123), vech(Omega2_123))

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=weightfun_pars222, mean_constraints=list(1:2), AR_constraints=C_222,
# parametrization="mean", cond_dist="Student"
theta_222exotcm <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), 13)
theta_222exotcm_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                       vech(Omega1_222), vech(Omega2_222), 13)


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
theta_222logistit_2_1_alt <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                               vec(B1_222), vec(B2_222)-vec(B1_222), c_and_gamma_222_2_1, dfs_222_2_1)

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="ind_Student"
dfs_122_1_1 <- c(4, 13)
B1_122 <- matrix(c(1.2, -0.3, 0.7, 0.1), nrow=2)
B2_122 <- matrix(c(0.5, 0.2, -0.1, 3.1), nrow=2)
theta_122logit_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(B1_122), vec(B2_122), gamma1_122_1_1, dfs_122_1_1)
theta_122logit_1_1_alt <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(B1_122), vec(B2_122)-vec(B1_122),
                            gamma1_122_1_1, dfs_122_1_1)


# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="ind_Student"
dfs_123_1_1 <- c(10, 12, 3)
B1_123 <- matrix(c(1.0, 0.3, 0.1, -0.8, 1.1, -0.5, -0.1, -0.2, 0.4), nrow=3)
B2_123 <- matrix(c(0.3, -0.2, -0.7, -0.8, 1.2, 0.5, 0.1, -0.2, 1.1), nrow=3)
theta_123expit_1_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(B1_123),
                        vec(B2_123), c_and_gamma_123_1_1, dfs_123_1_1)
theta_123expit_1_1_alt <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(B1_123),
                            vec(B2_123)-vec(B1_123), c_and_gamma_123_1_1, dfs_123_1_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student"
dfs_132_1_1 <- c(30, 6)
B1_132 <- matrix(c(0.6, 0.2, -0.1, 0.7), nrow=2)
B2_132 <- matrix(c(0.4, -0.1, -0.2, 0.5), nrow=2)
B3_132 <- matrix(c(0.9, -0.5, 0.2, 0.4), nrow=2)
theta_132thresit_1_1 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                          vec(B1_132), vec(B2_132), vec(B3_132), r1_132_1_1, r2_132_1_1, dfs_132_1_1)
theta_132thresit_1_1_alt <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                              vec(B1_132), vec(B2_132)-vec(B1_132), vec(B3_132)-vec(B1_132), r1_132_1_1, r2_132_1_1, dfs_132_1_1)

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="ind_Student",
# mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))
dfs_122_2_1 <- c(4, 13)
theta_222expcmwit_2_1 <- c(phi10_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222), xi_222expcmw_2_1, dfs_122_2_1)
theta_222expcmwit_2_1_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                    vec(B1_222), vec(B2_222), c(0.01, 0.33), dfs_122_2_1)
theta_222expcmwit_2_1_alt <- c(phi10_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222)-vec(B1_222),
                               xi_222expcmw_2_1, dfs_122_2_1)


# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student",
# mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))
theta_132thresmwit_1_1 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(B1_132),
                            vec(B2_132), vec(B3_132), dfs_132_1_1)
theta_132thresmwit_1_1_expanded <- c(phi10_132, phi20_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(B1_132),
                                     vec(B2_132), vec(B3_132), 0, 1.2, dfs_132_1_1)
theta_132thresmwit_1_1_alt <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(B1_132),
                                vec(B2_132)-vec(B1_132), vec(B3_132)-vec(B1_132), dfs_132_1_1)

# p=1, M=2, p=3, weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="ind_Student",
# mean_constraints=list(1:2), AR_constraints=C_123
dfs_123_3_1 <- c(11, 3, 20)
c_and_gamma_123_3_1 <- c(0.1, 0.4)
theta_123logisticcmit_3_1 <- c(phi10_123, vec(A11_123), vec(B1_123), vec(B2_123), c_and_gamma_123_3_1, dfs_123_3_1)
theta_123logisticcmit_3_1_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vec(B1_123), vec(B2_123),
                                        c_and_gamma_123_3_1, dfs_123_3_1)
theta_123logisticcmit_3_1_alt <- c(phi10_123, vec(A11_123), vec(B1_123), vec(B2_123)-vec(B1_123), c_and_gamma_123_3_1, dfs_123_3_1)


# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_Student", AR_constraints=C_222
theta_222logcit_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222), gamma1_222_2_1, dfs_222_2_1)
theta_222logcit_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                              vec(B1_222), vec(B2_222), gamma1_222_2_1, dfs_222_2_1)
theta_222logcit_2_1_alt <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222)-vec(B1_222),
                             gamma1_222_2_1, dfs_222_2_1)


# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=weightfun_pars222, cond_dist="ind_Student", AR_constraints=C_222
theta_222exoit_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222), dfs_222_2_1)
theta_222exoit_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                             vec(B1_222), vec(B2_222), dfs_222_2_1)
theta_222exoit_2_1_alt <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222)-vec(B1_222), dfs_222_2_1)


################
### ind_skewed_t

# p=1, M=1, d=2, cond_dist="ind_skewed_t", weight_function="threshold", weightfun_pars=c(1, 1)
dfls_112 <- c(3, 7, 0, 0)
B1_112 <- matrix(c(0.5, 0.2, -0.7, 0.3), nrow=2)
theta_112ikt <- c(phi10_112, vec(A11_112), vec(B1_112), dfls_112)

# p=2, M=2, d=2, cond_dist="ind_skewed_t", weight_function="logistic", weightfun_pars=c(2, 1)
dfls_222_2_1 <- c(3, 7, 0.3, 0.1)
B1_222 <- matrix(c(0.5, 0.2, -0.1, 0.3), nrow=2)
B2_222 <- matrix(c(0.4, -0.1, -0.2, 0.3), nrow=2)
theta_222logistikt_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                           vec(B1_222), vec(B2_222), c_and_gamma_222_2_1, dfls_222_2_1)
theta_222logistikt_2_1_alt <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                               vec(B1_222), vec(B2_222)-vec(B1_222), c_and_gamma_222_2_1, dfls_222_2_1)

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="ind_skewed_t"
dfls_122_1_1 <- c(4, 13, -0.1, 0)
B1_122 <- matrix(c(1.2, -0.3, 0.7, 0.1), nrow=2)
B2_122 <- matrix(c(0.5, 0.2, -0.1, 3.1), nrow=2)
theta_122logikt_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(B1_122), vec(B2_122), gamma1_122_1_1, dfls_122_1_1)
theta_122logikt_1_1_alt <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(B1_122), vec(B2_122)-vec(B1_122),
                            gamma1_122_1_1, dfls_122_1_1)


# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t"
dfls_123_1_1 <- c(10, 12, 3, 0.2, 0, -0.3)
B1_123 <- matrix(c(1.0, 0.3, 0.1, -0.8, 1.1, -0.5, -0.1, -0.2, 0.4), nrow=3)
B2_123 <- matrix(c(0.3, -0.2, -0.7, -0.8, 1.2, 0.5, 0.1, -0.2, 1.1), nrow=3)
theta_123expikt_1_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(B1_123),
                        vec(B2_123), c_and_gamma_123_1_1, dfls_123_1_1)
theta_123expikt_1_1_alt <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vec(B1_123),
                            vec(B2_123)-vec(B1_123), c_and_gamma_123_1_1, dfls_123_1_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t"
dfls_132_1_1 <- c(30, 6, -0.1, -0.2)
B1_132 <- matrix(c(0.6, 0.2, -0.1, 0.7), nrow=2)
B2_132 <- matrix(c(0.4, -0.1, -0.2, 0.5), nrow=2)
B3_132 <- matrix(c(0.9, -0.5, 0.2, 0.4), nrow=2)
theta_132thresikt_1_1 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                          vec(B1_132), vec(B2_132), vec(B3_132), r1_132_1_1, r2_132_1_1, dfls_132_1_1)
theta_132thresikt_1_1_alt <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                              vec(B1_132), vec(B2_132)-vec(B1_132), vec(B3_132)-vec(B1_132), r1_132_1_1, r2_132_1_1, dfls_132_1_1)

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="ind_skewed_t",
# mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))
dfls_122_2_1 <- c(4, 13, 0.1, 0.2)
theta_222expcmwikt_2_1 <- c(phi10_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222), xi_222expcmw_2_1, dfls_122_2_1)
theta_222expcmwikt_2_1_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                    vec(B1_222), vec(B2_222), c(0.01, 0.33), dfls_122_2_1)
theta_222expcmwikt_2_1_alt <- c(phi10_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222)-vec(B1_222),
                               xi_222expcmw_2_1, dfls_122_2_1)


# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t",
# mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))
theta_132thresmwikt_1_1 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(B1_132),
                            vec(B2_132), vec(B3_132), dfls_132_1_1)
theta_132thresmwikt_1_1_expanded <- c(phi10_132, phi20_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(B1_132),
                                     vec(B2_132), vec(B3_132), 0, 1.2, dfls_132_1_1)
theta_132thresmwikt_1_1_alt <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(B1_132),
                                vec(B2_132)-vec(B1_132), vec(B3_132)-vec(B1_132), dfls_132_1_1)

# p=1, M=2, p=3, weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="ind_skewed_t",
# mean_constraints=list(1:2), AR_constraints=C_123
dfls_123_3_1 <- c(11, 3, 20, 0, -0.1, 0.2)
c_and_gamma_123_3_1 <- c(0.1, 0.4)
theta_123logisticcmikt_3_1 <- c(phi10_123, vec(A11_123), vec(B1_123), vec(B2_123), c_and_gamma_123_3_1, dfls_123_3_1)
theta_123logisticcmikt_3_1_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), vec(B1_123), vec(B2_123),
                                        c_and_gamma_123_3_1, dfls_123_3_1)
theta_123logisticcmikt_3_1_alt <- c(phi10_123, vec(A11_123), vec(B1_123), vec(B2_123)-vec(B1_123), c_and_gamma_123_3_1, dfls_123_3_1)


# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_skewed_t", AR_constraints=C_222
theta_222logcikt_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222), gamma1_222_2_1, dfls_222_2_1)
theta_222logcikt_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                              vec(B1_222), vec(B2_222), gamma1_222_2_1, dfls_222_2_1)
theta_222logcikt_2_1_alt <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222)-vec(B1_222),
                             gamma1_222_2_1, dfls_222_2_1)


# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=weightfun_pars222, cond_dist="ind_skewed_t", AR_constraints=C_222
theta_222exoikt_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222), dfls_222_2_1)
theta_222exoikt_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                             vec(B1_222), vec(B2_222), dfls_222_2_1)
theta_222exoikt_2_1_alt <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(B1_222), vec(B2_222)-vec(B1_222), dfls_222_2_1)


#####################
### Structural models
# (recursively identified models use the same parametrization as reduced form models)

# p=1, M=2, d=2, weight_function="relative_dens", identification="heteroskedasticity"
W_122 <- matrix(c(-0.03, 0.24, -0.76, -0.02), nrow=2, ncol=2, byrow=FALSE)
lambdas_122 <- c(3.36, 0.86)
theta_122relgsh <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), lambdas_122, alpha1_122)

# p=1, M=3, d=2, weight_function="relative_dens", identification="heteroskedasticity"
df_222_2_1 <- 3
W_132 <- W_122; lambdas2_132 <- lambdas_122; lambdas3_132 <- c(2.1, 0.62)
theta_132relgsh <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                     vec(W_132), lambdas2_132, lambdas3_132, alpha1_132, alpha2_132)

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student",
# identification="heteroskedasticity"
W_222 <- W_122; lambdas_222 <- lambdas_122
theta_222logistictsh_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                              vec(W_222), lambdas_222, c_and_gamma_222_2_1, df_222_2_1)

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=1), identification="heteroskedasticity"
theta_122logsh_12_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), lambdas_122, gamma1_122_12_1)

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), identification="heteroskedasticity"
W_123 <- matrix(c(-0.47, -0.40, 1.25, 0.58, -1.01, 0.18, -0.66, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE)
lambdas_123 <- c(1.56, 1.44, 0.59)
theta_123expsh_1_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A11_123), vec(W_123), lambdas_123, c_and_gamma_123_1_1)

# p=2, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="Student",
# identification="heteroskedasticity"
df_232_1_1 <- 10
W_232 <- W_132; lambdas2_232 <- lambdas2_132; lambdas3_232 <- lambdas3_132
theta_232threstsh_1_1 <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                           vec(A31_232), vec(A32_232), vec(W_232), lambdas2_232, lambdas3_232, r1_232_1_1, r2_232_1_1,
                           df_232_1_1)

# p=1, M=3, d=2, weight_function="exogenous", weightfun_pars=weightfun_pars132, identification="heteroskedasticity"
theta_132exosh <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132), vec(W_132), lambdas2_132, lambdas3_132)

# p=1, M=1, d=2, cond_dist="ind_Student", weight_function="threshold", weightfun_pars=c(1, 1), identification="non-Gaussianity"
dfs_112 <- c(3, 7)
B1_112 <- matrix(c(0.5, 0.2, -0.7, 0.3), nrow=2)
theta_112it <- c(phi10_112, vec(A11_112), vec(B1_112), dfs_112)

# p=2, M=2, d=2, cond_dist="ind_Student", weight_function="logistic", weightfun_pars=c(2, 1), identification="non-Gaussianity"
dfs_222_2_1 <- c(3, 7)
B1_222 <- matrix(c(0.5, 0.2, -0.1, 0.3), nrow=2)
B2_222 <- matrix(c(0.4, -0.1, -0.2, 0.3), nrow=2)
theta_222logistit_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                           vec(B1_222), vec(B2_222), c_and_gamma_222_2_1, dfs_222_2_1)

## Structural models imposing constraints

# p=1, M=2, d=2, weight_function="relative_dens", identification="heteroskedasticity", AR_constraints=C_122
W_122 <- matrix(c(-0.03, 0.24, -0.76, -0.02), nrow=2, ncol=2, byrow=FALSE)
lambdas_122 <- c(3.36, 0.86)
theta_122relgshc <- c(phi10_122, phi20_122, vec(A11_122), vec(W_122), lambdas_122, alpha1_122)
theta_122relgshc_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A11_122), vec(W_122), lambdas_122, alpha1_122)

# p=1, M=2, d=2, weight_function="relative_dens", identification="heteroskedasticity", other_constraints=list(fixed_lambdas=c(3, 1))
theta_122relgshl <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), alpha1_122)
theta_122relgshl_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), 3, 1, alpha1_122)

# p=1, M=3, d=2, weight_function="relative_dens", identification="heteroskedasticity",
# B_constraints=matrix(c(0.1, NA, 0.3, 0), nrow=2)
W_132b <- matrix(c(0.11, 0.22, 0.33, 0), nrow=2);
theta_132relgshb <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                      Wvec(W_132b), lambdas2_132, lambdas3_132, alpha1_132, alpha2_132)
theta_132relgshb_expanded <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                               vec(W_132b), lambdas2_132, lambdas3_132, alpha1_132, alpha2_132)

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student", parametrization="mean"
# identification="heteroskedasticity", mean_constraints=list(1:2), B_constraints=matrix(c(0.1, 0, 0, 0.3), nrow=2)
W_222b <- matrix(c(0.12, 0, 0, 0.31), nrow=2)
theta_222logistictshmb_2_1 <- c(phi10_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                                Wvec(W_222b), lambdas_222, c_and_gamma_222_2_1, df_222_2_1)
theta_222logistictshmb_2_1_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                                         vec(W_222b), lambdas_222, c_and_gamma_222_2_1, df_222_2_1)

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student", parametrization="mean"
# identification="heteroskedasticity", mean_constraints=list(1:2), B_constraints=matrix(c(0.1, 0, 0, 0.3), nrow=2),
# other_constraints=list(fixed_lambdas=c(3, 1))
theta_222logistictshml_2_1 <- c(phi10_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                                Wvec(W_222b), c_and_gamma_222_2_1, df_222_2_1)
theta_222logistictshml_2_1_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
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

# p=2, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="Student",
# identification="heteroskedasticity", B_constraints=matrix(c(0.1, 0.2, -0.3, 0), nrow=2),
# other_constraints=list(fixed_lambdas=c(3, 1, 2, 0.5))
theta_232threstshl_1_1 <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                            vec(A31_232), vec(A32_232), Wvec(W_232b), r1_232_1_1, r2_232_1_1, df_232_1_1)
theta_232threstshl_1_1_expanded <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232),
                                     vec(A22_232), vec(A31_232), vec(A32_232), vec(W_232b), 3, 1, 2, 0.5, r1_232_1_1,
                                     r2_232_1_1, df_232_1_1)

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=weightfun_pars222, cond_dist="Student", parametrization="mean"
# identification="heteroskedasticity", mean_constraints=list(1:2), B_constraints=matrix(c(0.1, 0, 0, 0.3), nrow=2),
# other_constraints=list(fixed_lambdas=c(3, 1))
theta_222exotshml <- c(phi10_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222), Wvec(W_222b),  df_222_2_1)
theta_222exotshml_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222), vec(W_222b), 3, 1, df_222_2_1)

# p=1, M=2, d=3, weight_function="exogenous", weightfun_pars=weightfun_pars123, identification="heteroskedasticity",
# AR_constraints=C_123, B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE)
theta_123exoshcb <- c(phi10_123, phi20_123, vec(A11_123), Wvec(W_123b), lambdas_123)
theta_123exoshcb_expanded <- c(phi10_123, phi20_123, vec(A11_123), vec(A11_123), vec(W_123b), lambdas_123)

###############
### ind_Student

# p=2, M=2, d=2, cond_dist="ind_Student", weight_function="logistic", weightfun_pars=c(2, 1), identification="non-Gaussianity",
# B_constraints=matrix(c(NA, -1, 0, 1), nrow=2)
B1_222c <- matrix(c(0.5, -0.2, 0, 0.1), nrow=2)
B2_222c <- matrix(c(-0.4, -0.1, 0, 0.2), nrow=2)
theta_222logistitngb_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                              Wvec(B1_222c), Wvec(B2_222c), c_and_gamma_222_2_1, dfs_222_2_1)
theta_222logistitngb_2_1_alt <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                                  Wvec(B1_222c), Wvec(B2_222c)-Wvec(B1_222c), c_and_gamma_222_2_1, dfs_222_2_1)

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="ind_Student", identification="non-Gaussianity",
# B_constraints=matrix(c(1, NA, -1, 1) , nrow=2)
dfs_122_1_1 <- c(4, 13)
B1_122c <- matrix(c(1.2, 0.3, -0.7, 0.1), nrow=2)
B2_122c <- matrix(c(0.5, -0.9, -0.1, 3.1), nrow=2)
theta_122logitngb_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), Wvec(B1_122c), Wvec(B2_122c), gamma1_122_1_1, dfs_122_1_1)
theta_122logitngb_1_1_alt <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122),
                               Wvec(B1_122c), Wvec(B2_122c)-Wvec(B1_122c), gamma1_122_1_1, dfs_122_1_1)

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="ind_Student", identification="non-Gaussianity",
# B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)
dfs_123_1_1 <- c(10, 12, 3)
B1_123c <- matrix(c(1.0, 0.3, 0.1, 0, 1.1, -0.5, 0, -0.2, 0.4), nrow=3)
B2_123c <- matrix(c(0.3, -0.2, -0.7, 0, 1.2, 0.5, 0, -0.2, 1.1), nrow=3)
theta_123expitngb_1_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), Wvec(B1_123c),
                           Wvec(B2_123c), c_and_gamma_123_1_1, dfs_123_1_1)
theta_123expitngb_1_1_alt <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), Wvec(B1_123c),
                               Wvec(B2_123c)-Wvec(B1_123c), c_and_gamma_123_1_1, dfs_123_1_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student", identification="non-Gaussianity",
# B_constraints=matrix(c(1, 0, NA, 1), nrow=2)
dfs_132_1_1 <- c(30, 6)
B1_132c <- matrix(c(0.6, 0, -0.1, 0.7), nrow=2)
B2_132c <- matrix(c(0.4, 0, 0.2, 0.5), nrow=2)
B3_132c <- matrix(c(0.9, 0, -0.2, 0.4), nrow=2)
theta_132thresitngb_1_1 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                             Wvec(B1_132c), Wvec(B2_132c), Wvec(B3_132c), r1_132_1_1, r2_132_1_1, dfs_132_1_1)
theta_132thresitngb_1_1_alt <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                                 Wvec(B1_132c), Wvec(B2_132c)-Wvec(B1_132c), Wvec(B3_132c)-Wvec(B1_132c),
                                 r1_132_1_1, r2_132_1_1, dfs_132_1_1)

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="ind_Student",
# mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
# identification="non-Gaussianity", B_constraints=matrix(c(NA, -1, 0, 1), nrow=2)
dfs_222_2_1 <- c(4, 13)
theta_222expcmwitngb_2_1 <- c(phi10_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c), xi_222expcmw_2_1, dfs_222_2_1)
theta_222expcmwitngb_2_1_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                       Wvec(B1_222c), Wvec(B2_222c), c(0.01, 0.33), dfs_222_2_1)
theta_222expcmwitngb_2_1_alt <- c(phi10_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c)-Wvec(B1_222c),
                                  xi_222expcmw_2_1, dfs_222_2_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student", mean_constraints=list(1, 2:3),
# weight_constraints=list(R=0, r=c(0, 1.2)), identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2)
theta_132thresmwitngb_1_1 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), Wvec(B1_132c),
                               Wvec(B2_132c), Wvec(B3_132c), dfs_132_1_1)
theta_132thresmwitngb_1_1_expanded <- c(phi10_132, phi20_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), Wvec(B1_132c),
                                        Wvec(B2_132c), Wvec(B3_132c), 0, 1.2, dfs_132_1_1)
theta_132thresmwitngb_1_1_alt <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), Wvec(B1_132c),
                                   Wvec(B2_132c)-Wvec(B1_132c), Wvec(B3_132c)-Wvec(B1_132c), dfs_132_1_1)

# p=1, M=2, d=3, weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="ind_Student", mean_constraints=list(1:2),
# AR_constraints=C_123, identification="non-Gaussianity", B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)
dfs_123_3_1 <- c(11, 3, 20)
theta_123logisticcmitngb_3_1 <- c(phi10_123, vec(A11_123), Wvec(B1_123c), Wvec(B2_123c), c_and_gamma_123_3_1, dfs_123_3_1)
theta_123logisticcmitngb_3_1_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), Wvec(B1_123c), Wvec(B2_123c),
                                           c_and_gamma_123_3_1, dfs_123_3_1)
theta_123logisticcmitngb_3_1_alt <- c(phi10_123, vec(A11_123), Wvec(B1_123c), Wvec(B2_123c)-Wvec(B1_123c), c_and_gamma_123_3_1, dfs_123_3_1)


# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_Student", AR_constraints=C_222,
# identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)
theta_222logcitngb_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c), gamma1_222_2_1, dfs_222_2_1)
theta_222logcitngb_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                 Wvec(B1_222c), Wvec(B2_222c), gamma1_222_2_1, dfs_222_2_1)
theta_222logcitngb_2_1_alt <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c)-Wvec(B1_222c),
                                gamma1_222_2_1, dfs_222_2_1)

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=weightfun_pars222, cond_dist="ind_Student",
# AR_constraints=C_222, identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)
theta_222exoitngb_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c), dfs_222_2_1)
theta_222exoitngb_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                Wvec(B1_222c), Wvec(B2_222c), dfs_222_2_1)
theta_222exoitngb_2_1_alt <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c)-Wvec(B1_222c), dfs_222_2_1)


###############
### ind_skewd_t

# p=2, M=2, d=2, cond_dist="ind_skewed_t", weight_function="logistic", weightfun_pars=c(2, 1), identification="non-Gaussianity",
# B_constraints=matrix(c(NA, -1, 0, 1), nrow=2)
B1_222c <- matrix(c(0.5, -0.2, 0, 0.1), nrow=2)
B2_222c <- matrix(c(-0.4, -0.1, 0, 0.2), nrow=2)
theta_222logistiktngb_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                              Wvec(B1_222c), Wvec(B2_222c), c_and_gamma_222_2_1, dfls_222_2_1)
theta_222logistiktngb_2_1_alt <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222),
                                  Wvec(B1_222c), Wvec(B2_222c)-Wvec(B1_222c), c_and_gamma_222_2_1, dfls_222_2_1)

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="ind_skewed_t", identification="non-Gaussianity",
# B_constraints=matrix(c(1, NA, -1, 1) , nrow=2)
dfls_122_1_1 <- c(4, 13, -0.4, 0.5)
B1_122c <- matrix(c(1.2, 0.3, -0.7, 0.1), nrow=2)
B2_122c <- matrix(c(0.5, -0.9, -0.1, 3.1), nrow=2)
theta_122logiktngb_1_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), Wvec(B1_122c), Wvec(B2_122c), gamma1_122_1_1, dfls_122_1_1)
theta_122logiktngb_1_1_alt <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122),
                               Wvec(B1_122c), Wvec(B2_122c)-Wvec(B1_122c), gamma1_122_1_1, dfls_122_1_1)

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t", identification="non-Gaussianity",
# B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)
dfls_123_1_1 <- c(10, 12, 3, 0.3, -0.1, -0.2)
B1_123c <- matrix(c(1.0, 0.3, 0.1, 0, 1.1, -0.5, 0, -0.2, 0.4), nrow=3)
B2_123c <- matrix(c(0.3, -0.2, -0.7, 0, 1.2, 0.5, 0, -0.2, 1.1), nrow=3)
theta_123expiktngb_1_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), Wvec(B1_123c),
                           Wvec(B2_123c), c_and_gamma_123_1_1, dfls_123_1_1)
theta_123expiktngb_1_1_alt <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), Wvec(B1_123c),
                               Wvec(B2_123c)-Wvec(B1_123c), c_and_gamma_123_1_1, dfls_123_1_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t", identification="non-Gaussianity",
# B_constraints=matrix(c(1, 0, NA, 1), nrow=2)
dfls_132_1_1 <- c(30, 6, 0, 0)
B1_132c <- matrix(c(0.6, 0, -0.1, 0.7), nrow=2)
B2_132c <- matrix(c(0.4, 0, 0.2, 0.5), nrow=2)
B3_132c <- matrix(c(0.9, 0, -0.2, 0.4), nrow=2)
theta_132thresiktngb_1_1 <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                             Wvec(B1_132c), Wvec(B2_132c), Wvec(B3_132c), r1_132_1_1, r2_132_1_1, dfls_132_1_1)
theta_132thresiktngb_1_1_alt <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                                 Wvec(B1_132c), Wvec(B2_132c)-Wvec(B1_132c), Wvec(B3_132c)-Wvec(B1_132c),
                                 r1_132_1_1, r2_132_1_1, dfls_132_1_1)

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="ind_skewed_t",
# mean_constraints=list(1:2), AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
# identification="non-Gaussianity", B_constraints=matrix(c(NA, -1, 0, 1), nrow=2)
dfls_222_2_1 <- c(4, 13, -0.1, 0.3)
theta_222expcmwiktngb_2_1 <- c(phi10_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c), xi_222expcmw_2_1, dfls_222_2_1)
theta_222expcmwiktngb_2_1_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                       Wvec(B1_222c), Wvec(B2_222c), c(0.01, 0.33), dfls_222_2_1)
theta_222expcmwiktngb_2_1_alt <- c(phi10_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c)-Wvec(B1_222c),
                                  xi_222expcmw_2_1, dfls_222_2_1)

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_skewed_t", mean_constraints=list(1, 2:3),
# weight_constraints=list(R=0, r=c(0, 1.2)), identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2)
theta_132thresmwiktngb_1_1 <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), Wvec(B1_132c),
                               Wvec(B2_132c), Wvec(B3_132c), dfls_132_1_1)
theta_132thresmwiktngb_1_1_expanded <- c(phi10_132, phi20_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), Wvec(B1_132c),
                                        Wvec(B2_132c), Wvec(B3_132c), 0, 1.2, dfls_132_1_1)
theta_132thresmwiktngb_1_1_alt <- c(phi10_132, phi20_132, vec(A11_132), vec(A21_132), vec(A31_132), Wvec(B1_132c),
                                   Wvec(B2_132c)-Wvec(B1_132c), Wvec(B3_132c)-Wvec(B1_132c), dfls_132_1_1)

# p=1, M=2, d=3, weight_function="logistic", weightfun_pars=c(3, 1), cond_dist="ind_skewed_t", mean_constraints=list(1:2),
# AR_constraints=C_123, identification="non-Gaussianity", B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)
dfls_123_3_1 <- c(11, 3, 20, 0.1, 0.2, 0.3)
theta_123logisticcmiktngb_3_1 <- c(phi10_123, vec(A11_123), Wvec(B1_123c), Wvec(B2_123c), c_and_gamma_123_3_1, dfls_123_3_1)
theta_123logisticcmiktngb_3_1_expanded <- c(phi10_123, phi10_123, vec(A11_123), vec(A11_123), Wvec(B1_123c), Wvec(B2_123c),
                                           c_and_gamma_123_3_1, dfls_123_3_1)
theta_123logisticcmiktngb_3_1_alt <- c(phi10_123, vec(A11_123), Wvec(B1_123c), Wvec(B2_123c)-Wvec(B1_123c), c_and_gamma_123_3_1, dfls_123_3_1)


# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), cond_dist="ind_skewed_t", AR_constraints=C_222,
# identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)
theta_222logciktngb_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c), gamma1_222_2_1, dfls_222_2_1)
theta_222logciktngb_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                 Wvec(B1_222c), Wvec(B2_222c), gamma1_222_2_1, dfls_222_2_1)
theta_222logciktngb_2_1_alt <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c)-Wvec(B1_222c),
                                gamma1_222_2_1, dfls_222_2_1)

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=weightfun_pars222, cond_dist="ind_skewed_t",
# AR_constraints=C_222, identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)
theta_222exoiktngb_2_1 <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c), dfls_222_2_1)
theta_222exoiktngb_expanded <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                Wvec(B1_222c), Wvec(B2_222c), dfls_222_2_1)
theta_222exoiktngb_2_1_alt <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), Wvec(B1_222c), Wvec(B2_222c)-Wvec(B1_222c), dfls_222_2_1)


test_that("loglikelihood works correctly", {

  # Relative_dens Gaussian STVAR
  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112relg, weight_function="relative_dens"), -1000.653, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=1, params=theta_212relg, weight_function="relative_dens"), -286.5474, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=8, M=1, params=theta_812relg, weight_function="relative_dens"), -257.0505, tolerance=1e-3)
  expect_equal(loglikelihood(data=data2, p=4, M=1, params=theta_413relg, weight_function="relative_dens"), -596.6938, tolerance=1e-3)

  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122relg, weight_function="relative_dens"), -314.6693, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222relg, weight_function="relative_dens"), -239.3485, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222relg, weight_function="relative_dens",
                             penalized=TRUE, penalty_params=c(0.5, 0.1)), -254.771, tolerance=1e-3)

  expect_equal(loglikelihood(data=usamone, p=3, M=1, params=theta_313relg, weight_function="relative_dens"), -669.5716, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=3, M=1, params=theta_313relg_unstab, weight_function="threshold", weightfun_pars=c(1, 1),
                             allow_unstab=TRUE), -7491.638, tolerance=1e-2)

  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123relg, weight_function="relative_dens"), -570.019, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=3, M=2, params=theta_323relg, weight_function="relative_dens"), -490.9401, tolerance=1e-3)

  # logistic STVAR
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122logistic_1_1, weight_function="logistic",
                             weightfun_pars=c(1, 1)), -342.6918, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistic_2_1, weight_function="logistic",
                             weightfun_pars=c(2, 1)), -291.4979, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistic_1_2, weight_function="logistic",
                             weightfun_pars=c(1, 2)), -303.1628, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistic_1_2, weight_function="logistic",
                             weightfun_pars=c(1, 2), penalized=TRUE, penalty_params=c(0.3, 0.2)), -308.7045, tolerance=1e-3)
  tmp_pars <- theta_222logistic_1_2
  tmp_pars[8] <- 3
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=tmp_pars, weight_function="logistic", weightfun_pars=c(1, 2),
                             penalized=TRUE, penalty_params=c(0.3, 0.2), allow_unstab=TRUE), -3211.04, tolerance=1e-2)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=tmp_pars, weight_function="logistic", weightfun_pars=c(1, 2),
                             penalized=TRUE, penalty_params=c(0.01, 0.001), allow_unstab=TRUE), -2653.194, tolerance=1e-3)

  # mlogit STVAR
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122log_1_1, weight_function="mlogit",
                             weightfun_pars=list(vars=1, lags=1)), -334.0843, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122log_12_1, weight_function="mlogit",
                             weightfun_pars=list(vars=1:2, lags=1)), -335.5687, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222log_2_1, weight_function="mlogit",
                             weightfun_pars=list(vars=2, lags=1)), -315.326, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222log_12_2, weight_function="mlogit",
                             weightfun_pars=list(vars=1:2, lags=2)), -344.0656, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=3, params=theta_232log_1_1, weight_function="mlogit",
                             weightfun_pars=list(vars=1, lags=1)), -4260.053, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=3, params=theta_232log_12_1, weight_function="mlogit",
                             weightfun_pars=list(vars=1:2, lags=1)), -3756.747, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=3, params=theta_232log_2_2, weight_function="mlogit",
                             weightfun_pars=list(vars=2, lags=2)), -3450.811, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=3, params=theta_232log_12_2, weight_function="mlogit",
                             weightfun_pars=list(vars=1:2, lags=2)), -2695.943, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123log_1_1, weight_function="mlogit",
                             weightfun_pars=list(vars=1, lags=1)), -998.6099, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123log_23_1, weight_function="mlogit",
                             weightfun_pars=list(vars=2:3, lags=1)), -1100.063, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123log_123_1, weight_function="mlogit",
                             weightfun_pars=list(vars=1:3, lags=1)), -1155.278, tolerance=1e-3)

  # Exponential STVAR
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122exp_1_1, weight_function="exponential",
                             weightfun_pars=c(1, 1)), -336.8664, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222exp_2_1, weight_function="exponential",
                             weightfun_pars=c(2, 1)), -308.0051, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222exp_1_2, weight_function="exponential",
                             weightfun_pars=c(1, 2)), -407.7225, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123exp_1_1, weight_function="exponential",
                             weightfun_pars=c(1, 1)), -1115.766, tolerance=1e-3)

  # Threshold STVAR
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122thres_1_1, weight_function="threshold",
                             weightfun_pars=c(1, 1)), -391.7463, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122thres_1_1, weight_function="threshold",
                             weightfun_pars=c(1, 1), parametrization="mean"), -427.092, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222thres_2_2, weight_function="threshold",
                             weightfun_pars=c(2, 2)), -282.805, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=theta_132thres_1_1, weight_function="threshold",
                             weightfun_pars=c(1, 1)), -5944.922, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=3, params=theta_232thres_1_1, weight_function="threshold",
                             weightfun_pars=c(1, 1)), -8282.43, tolerance=1e-2)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123thres_2_1, weight_function="threshold",
                             weightfun_pars=c(2, 1)), -682.5003, tolerance=1e-3)

  # Exogenous STVAR
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222exo, weight_function="exogenous",
                             weightfun_pars=weightfun_pars222), -331.7587, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=theta_132exo, weight_function="exogenous",
                             weightfun_pars=weightfun_pars132), -6906.911, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123exo, weight_function="exogenous",
                             weightfun_pars=weightfun_pars123), -2045.486, tolerance=1e-3)

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
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122logc_12_1, weight_function="mlogit",
                             weightfun_pars=list(vars=1:2, lags=1), AR_constraints=C_122), -335.5687, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logc_2_1, weight_function="mlogit",
                             weightfun_pars=list(vars=2, lags=1), AR_constraints=C_222), -369.6914, tolerance=1e-3)

  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112relgm, weight_function="relative_dens", parametrization="mean",
                             mean_constraints=list(1)), -302.0188, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122relgm, weight_function="relative_dens", parametrization="mean",
                             mean_constraints=list(1:2)), -326.8462, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222relgcm, weight_function="relative_dens", parametrization="mean",
                             AR_constraints=C_222, mean_constraints=list(1:2)), -301.0143, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122logm_1_1, weight_function="mlogit", parametrization="mean",
                             weightfun_pars=list(vars=1, lags=1), mean_constraints=list(1:2)), -368.2033, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logcm_12_2, weight_function="mlogit", parametrization="mean",
                             weightfun_pars=list(vars=1:2, lags=2), AR_constraints=C_222, mean_constraints=list(1:2)),
               -453.5176, tolerance=1e-3)

  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222relgcmw, weight_function="relative_dens",
                             mean_constraints=list(1:2), AR_constraints=C_222, parametrization="mean",
                             weight_constraints=list(R=0, r=0.6)), -348.8925, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122logw_1_1, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1),
                             weight_constraints=list(R=0, r=c(0.12, 0.13))), -334.8481, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logcmw_12_2, weight_function="mlogit", parametrization="mean",
                             weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222,
                             weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5), r=c(0, 0.11, 0.12, 0.13, 0))),
               -458.0816, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122logisticw_1_1, weight_function="logistic", weightfun_pars=c(1, 1),
                             weight_constraints=list(R=0, r=c(0.02, 0.13))), -341.6723, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logisticcmw_2_1, weight_function="logistic", weightfun_pars=c(2, 1),
                             mean_constraints=list(1:2), AR_constraints=C_222, parametrization="mean",
                             weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))), -382.3062, tolerance=1e-3)

  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122expw_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                             weight_constraints=list(R=0, r=c(0.02, 0.13))), -336.0778, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222expcmw_2_1, weight_function="exponential", weightfun_pars=c(2, 1),
                             mean_constraints=list(1:2), AR_constraints=C_222, parametrization="mean",
                             weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))), -386.0939, tolerance=1e-3)

  expect_equal(loglikelihood(data=gdpdef, p=2, M=3, params=theta_232thresw_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                             weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.1, 0))), -16285.82, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=theta_132thresmw_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                             mean_constraints=list(1, 2:3), parametrization="mean", weight_constraints=list(R=0, r=c(0, 1.2))),
               -467.6997, tolerance=1e-3)

  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123exoc, weight_function="exogenous", weightfun_pars=weightfun_pars123,
                             mean_constraints=list(1:2), AR_constraints=C_123, parametrization="mean"), -2093.978, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222exotcm,  weight_function="exogenous", weightfun_pars=weightfun_pars222,
                             mean_constraints=list(1:2), AR_constraints=C_222, parametrization="mean", cond_dist="Student"),
               -377.2813, tolerance=1e-3)

  ## Student
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=c(theta_122logistic_1_1, 300), weight_function="logistic",
                             weightfun_pars=c(1, 1), cond_dist="Student"), -341.7576, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=c(theta_222logistic_2_1, 10), weight_function="logistic",
                             weightfun_pars=c(2, 1), cond_dist="Student"), -274.9028, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=c(theta_222logistic_1_2, 3), weight_function="logistic",
                             weightfun_pars=c(1, 2), cond_dist="Student"), -276.1146, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=3, params=c(theta_232log_2_2, 13), weight_function="mlogit",
                             weightfun_pars=list(vars=2, lags=2), cond_dist="Student"), -2300.686, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=3, params=c(theta_232log_12_2, 2.1), weight_function="mlogit",
                             weightfun_pars=list(vars=1:2, lags=2), cond_dist="Student"), -1987.197, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=c(theta_123log_1_1, 30), weight_function="mlogit",
                             weightfun_pars=list(vars=1, lags=1), cond_dist="Student"), -954.2666, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=c(theta_222exp_1_2, 7), weight_function="exponential",
                             weightfun_pars=c(1, 2), cond_dist="Student"), -312.2618, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=c(theta_123exp_1_1, 12), weight_function="exponential",
                             weightfun_pars=c(1, 1), cond_dist="Student"), -1001.058, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=3, params=c(theta_232thres_1_1, 10), weight_function="threshold",
                             weightfun_pars=c(1, 1), cond_dist="Student"), -2233.18, tolerance=1e-2)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=c(theta_123thres_2_1, 2.01), weight_function="threshold",
                             weightfun_pars=c(2, 1), cond_dist="Student"), -1448.457, tolerance=1e-2)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=c(theta_123thres_2_1, 2.01), weight_function="threshold",
                             weightfun_pars=c(2, 1), cond_dist="Student", penalized=TRUE, penalty_params=c(0.2, 0.1)),
               -1452.174, tolerance=1e-2)

  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=c(theta_122logc_12_1, 4.4), weight_function="mlogit",
                             weightfun_pars=list(vars=1:2, lags=1), AR_constraints=C_122, cond_dist="Student"), -307.7079, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=c(theta_222logisticcmw_2_1, 30), weight_function="logistic",
                             weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222, cond_dist="Student",
                             parametrization="mean", weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))),
               -368.2865, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=c(theta_122expw_1_1, 2.3), weight_function="exponential", weightfun_pars=c(1, 1),
                             weight_constraints=list(R=0, r=c(0.02, 0.13)), cond_dist="Student"), -380.4028, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=c(theta_132thresmw_1_1, 11), weight_function="threshold", weightfun_pars=c(1, 1),
                             mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2)), cond_dist="Student",
                             parametrization="mean"), -447.6652, tolerance=1e-3)
  expect_equal(c(loglikelihood(data=gdpdef, p=1, M=3, params=c(theta_132thresmw_1_1, 11), weight_function="threshold", weightfun_pars=c(1, 1),
                             mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2)), cond_dist="Student",
                             parametrization="mean", to_return="total_ccovs")[, , 1]), c(1.0, 0.5, 0.5, 1.0), tolerance=1e-3)

  # ind_Student
  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112it,  cond_dist="ind_Student", weight_function="threshold",
                             weightfun_pars=c(1, 1)), -851.9794, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistit_2_1, cond_dist="ind_Student", weight_function="logistic",
                             weightfun_pars=c(2, 1)), -336.6867, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122logit_1_1, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1),
                             cond_dist="ind_Student"), -580.561, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123expit_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                             cond_dist="ind_Student"), -3830.948, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222expcmwit_2_1, weight_function="exponential", weightfun_pars=c(2, 1),
                             cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_222,
                             weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))), -393.5391, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=theta_132thresmwit_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                             cond_dist="ind_Student", mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))),
               -540.3383, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123logisticcmit_3_1, weight_function="logistic", weightfun_pars=c(3, 1),
                             cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_123), -6872.711, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logcit_2_1, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1),
                             cond_dist="ind_Student", AR_constraints=C_222), -380.7479, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222exoit_2_1, weight_function="exogenous",
                             weightfun_pars=weightfun_pars222, cond_dist="ind_Student", AR_constraints=C_222),
               -397.4386, tolerance=1e-3)

  # ind_skewed_t
  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112ikt,  cond_dist="ind_skewed_t", weight_function="threshold",
                             weightfun_pars=c(1, 1)), -851.9794, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistikt_2_1, cond_dist="ind_skewed_t", weight_function="logistic",
                             weightfun_pars=c(2, 1)), -350.3691, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122logikt_1_1, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1),
                             cond_dist="ind_skewed_t"), -578.4793, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123expikt_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                             cond_dist="ind_skewed_t"), -4050.617, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222expcmwikt_2_1, weight_function="exponential", weightfun_pars=c(2, 1),
                             cond_dist="ind_skewed_t", mean_constraints=list(1:2), AR_constraints=C_222,
                             weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))), -391.6119, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=theta_132thresmwikt_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                             cond_dist="ind_skewed_t", mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2))),
               -556.0658, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123logisticcmikt_3_1, weight_function="logistic", weightfun_pars=c(3, 1),
                             cond_dist="ind_skewed_t", mean_constraints=list(1:2), AR_constraints=C_123), -7225.257, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logcikt_2_1, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1),
                             cond_dist="ind_skewed_t", AR_constraints=C_222), -403.1479, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222exoikt_2_1, weight_function="exogenous",
                             weightfun_pars=weightfun_pars222, cond_dist="ind_skewed_t", AR_constraints=C_222),
               -422.6878, tolerance=1e-3)


  # Structural models identified by heteroskedasticity
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122relgsh, weight_function="relative_dens",
                             identification="heteroskedasticity"), -313.0517, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122relgsh, weight_function="relative_dens",
                             identification="heteroskedasticity", parametrization="mean"), -357.8102, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=theta_132relgsh, weight_function="relative_dens",
                             identification="heteroskedasticity"), -311.6129, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistictsh_2_1, weight_function="logistic", weightfun_pars=c(2, 1),
                             cond_dist="Student", identification="heteroskedasticity"), -322.2773, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123expsh_1_1,  weight_function="exponential", weightfun_pars=c(1, 1),
                             identification="heteroskedasticity"), -2336.49, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=3, params=theta_232threstsh_1_1, weight_function="threshold",
                             weightfun_pars=c(1, 1), cond_dist="Student", identification="heteroskedasticity"),
               -2988.147, tolerance=1e-3)

  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122relgshc, weight_function="relative_dens",
                             identification="heteroskedasticity", AR_constraints=C_122), -313.0517, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=theta_132relgshb, weight_function="relative_dens",
                             identification="heteroskedasticity",
                             B_constraints=matrix(c(0.1, NA, 0.3, 0), nrow=2)), -755.7973, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistictshmb_2_1, weight_function="logistic", weightfun_pars=c(2, 1),
                             cond_dist="Student", identification="heteroskedasticity", mean_constraints=list(1:2), parametrization="mean",
                             B_constraints=matrix(c(0.1, 0, 0, 0.3), nrow=2)), -776.4632, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122logshwb_12_1, weight_function="mlogit",
                             weightfun_pars=list(vars=1:2, lags=1),
                             identification="heteroskedasticity", weight_constraints=list(R=0, r=c(0.1, 0.2, 0.3)),
                             B_constraints=matrix(c(0.1, 0.2, 0.3, 0), nrow=2)), -755.1974, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122logshwb_12_1, weight_function="mlogit",
                             weightfun_pars=list(vars=1:2, lags=1),
                             identification="heteroskedasticity", weight_constraints=list(R=0, r=c(0.1, 0.2, 0.3)), parametrization="mean",
                             B_constraints=matrix(c(0.1, 0.2, 0.3, 0), nrow=2)), -806.3316, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123expshcwb_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                             identification="heteroskedasticity", AR_constraints=C_123,
                             weight_constraints=list(R=matrix(c(1, 0.5), nrow=2), r=c(0, 0)),
                             B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE)),
               -9514.231, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123expshcwb_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                             identification="heteroskedasticity", AR_constraints=C_123,
                             weight_constraints=list(R=matrix(c(1, 0.5), nrow=2), r=c(0, 0)),
                             B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE),
                             parametrization="mean"), -9557.789, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=3, params=theta_232threstshb_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                             cond_dist="Student", identification="heteroskedasticity",
                             B_constraints=matrix(c(0.1, 0.2, -0.3, 0), nrow=2)), -4541.691, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=3, params=theta_232threstshl_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                             cond_dist="Student", identification="heteroskedasticity", B_constraints=matrix(c(0.1, 0.2, -0.3, 0), nrow=2),
                             other_constraints=list(fixed_lambdas=c(3, 1, 2, 0.5))), -4529.612, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=3, params=theta_232threstshl_1_1_expanded, weight_function="threshold",
                             weightfun_pars=c(1, 1), cond_dist="Student", identification="heteroskedasticity"), -4529.612, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122relgshl, weight_function="relative_dens",
                             identification="heteroskedasticity", other_constraints=list(fixed_lambdas=c(3, 1))), -305.7633, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122relgshl_expanded, weight_function="relative_dens",
                             identification="heteroskedasticity"), -305.7633, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistictshml_2_1, weight_function="logistic", weightfun_pars=c(2, 1),
                             cond_dist="Student", parametrization="mean", identification="heteroskedasticity", mean_constraints=list(1:2),
                             B_constraints=matrix(c(0.1, 0, 0, 0.3), nrow=2), other_constraints=list(fixed_lambdas=c(3, 1))),
               -806.8079, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistictshml_2_1_expanded, weight_function="logistic",
                             weightfun_pars=c(2, 1), cond_dist="Student", parametrization="mean", identification="heteroskedasticity"),
               -806.8079, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123expshcwl,weight_function="exponential", weightfun_pars=c(1, 1),
                             identification="heteroskedasticity", AR_constraints=C_123,
                             weight_constraints=list(R=matrix(c(1, 0.5), nrow=2), r=c(0, 0)),
                             B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE),
                             other_constraints=list(fixed_lambdas=c(3, 1, 2))),
               -6573.02, tolerance=1e-2)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123expshcwl_expanded, weight_function="exponential", weightfun_pars=c(1, 1),
                             identification="heteroskedasticity"),
               -6573.02, tolerance=1e-2)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=theta_132exosh, weight_function="exogenous",
                             weightfun_pars=weightfun_pars132, identification="heteroskedasticity"), -32505.69, tolerance=1e-3)
  expect_equal(c(loglikelihood(data=usamone, p=1, M=2, params=theta_123expshcwl_expanded, weight_function="exponential", weightfun_pars=c(1, 1),
                             identification="heteroskedasticity", to_return="total_ccovs")[, , 269]),
               c(0.6331981, -0.3332059, -0.3828000, -0.3332059, 2.2054354, 1.9355345, -0.3828000, 1.9355345, 2.0949759), tolerance=1e-3)

  ###  Structural models identified by non-Gaussianity

  # ind_Student
  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112it, cond_dist="ind_Student", weight_function="threshold",
                             weightfun_pars=c(1, 1), identification="non-Gaussianity"), -851.9794, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistit_2_1, cond_dist="ind_Student", weight_function="logistic",
                             weightfun_pars=c(2, 1), identification="non-Gaussianity"), -336.6867, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistitngb_2_1, cond_dist="ind_Student", weight_function="logistic",
                             weightfun_pars=c(2, 1), identification="non-Gaussianity", B_constraints=matrix(c(NA, -1, 0, 1), nrow=2)),
               -5345.796, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122logitngb_1_1, weight_function="mlogit",
                             weightfun_pars=list(vars=1, lags=1), cond_dist="ind_Student", identification="non-Gaussianity",
                            B_constraints=matrix(c(1, NA, -1, 1) , nrow=2)), -580.1779, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123expitngb_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                             cond_dist="ind_Student", identification="non-Gaussianity",
                             B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)), -2782.476, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=theta_132thresitngb_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                             cond_dist="ind_Student", identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2)),
               -3827.6, tolerance=1e-1)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222expcmwitngb_2_1, weight_function="exponential", weightfun_pars=c(2, 1),
                             cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_222,
                             weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), identification="non-Gaussianity",
                             B_constraints=matrix(c(NA, -1, 0, 1), nrow=2)),
               -1589.59, tolerance=1e-2)

  tmp_pars <- theta_222expcmwitngb_2_1
  tmp_pars[9] <- 3
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=tmp_pars, weight_function="exponential", weightfun_pars=c(2, 1),
                             cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_222,
                             weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), identification="non-Gaussianity",
                             B_constraints=matrix(c(NA, -1, 0, 1), nrow=2), allow_unstab=TRUE,
                             penalized=TRUE, penalty_params=c(0.1, 0.11)), -5393.25, tolerance=1e-2)


  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=theta_132thresmwitngb_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                             cond_dist="ind_Student", mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2)),
                             identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2)),
               -508.1507, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123logisticcmitngb_3_1, weight_function="logistic", weightfun_pars=c(3, 1),
                             cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_123, identification="non-Gaussianity",
                             B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)),
               -3988.082, tolerance=1e-2)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logcitngb_2_1, weight_function="mlogit",
                             weightfun_pars=list(vars=2, lags=1), cond_dist="ind_Student", AR_constraints=C_222,
                             identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)),
               -2820.663, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222exoitngb_2_1, weight_function="exogenous",
                             weightfun_pars=weightfun_pars222, cond_dist="ind_Student", AR_constraints=C_222,
                             identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)),
               -2657.869, tolerance=1e-3)
  expect_equal(c(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222exoitngb_2_1, weight_function="exogenous",
                             weightfun_pars=weightfun_pars222, cond_dist="ind_Student", AR_constraints=C_222,
                             identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2), to_return="total_ccovs")[, , 242]),
               c(0.07106374, -0.04640166, -0.04640166, 0.04615817), tolerance=1e-3)

  # ind_Student / non-Gaussianity alt parametrization
  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112it,  cond_dist="ind_Student", weight_function="threshold",
                             weightfun_pars=c(1, 1), alt_par=TRUE), -851.9794, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistit_2_1_alt,
                             cond_dist="ind_Student", weight_function="logistic",
                             weightfun_pars=c(2, 1), alt_par=TRUE), -336.6867, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122logit_1_1_alt, weight_function="mlogit",
                             weightfun_pars=list(vars=1, lags=1), cond_dist="ind_Student", alt_par=TRUE), -580.561, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123expit_1_1_alt, weight_function="exponential", weightfun_pars=c(1, 1),
                             cond_dist="ind_Student", alt_par=TRUE), -3830.948, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222expcmwit_2_1_alt, weight_function="exponential", weightfun_pars=c(2, 1),
                             cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_222,
                             weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), alt_par=TRUE), -393.5391, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=theta_132thresmwit_1_1_alt, weight_function="threshold", weightfun_pars=c(1, 1),
                             cond_dist="ind_Student", mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2)), alt_par=TRUE),
               -540.3383, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123logisticcmit_3_1_alt, weight_function="logistic", weightfun_pars=c(3, 1),
                             cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_123, alt_par=TRUE),
               -6872.711, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logcit_2_1_alt, weight_function="mlogit",
                             weightfun_pars=list(vars=2, lags=1), cond_dist="ind_Student", AR_constraints=C_222, alt_par=TRUE),
               -380.7479, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222exoit_2_1_alt, weight_function="exogenous",
                             weightfun_pars=weightfun_pars222, cond_dist="ind_Student", AR_constraints=C_222, alt_par=TRUE),
               -397.4386, tolerance=1e-3)

  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112it, cond_dist="ind_Student", weight_function="threshold",
                             weightfun_pars=c(1, 1), identification="non-Gaussianity", alt_par=TRUE), -851.9794, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistit_2_1_alt, cond_dist="ind_Student", weight_function="logistic",
                             weightfun_pars=c(2, 1), identification="non-Gaussianity", alt_par=TRUE), -336.6867, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistitngb_2_1_alt, cond_dist="ind_Student", weight_function="logistic",
                             weightfun_pars=c(2, 1), identification="non-Gaussianity", B_constraints=matrix(c(NA, -1, 0, 1), nrow=2),
                             alt_par=TRUE), -5345.796, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122logitngb_1_1_alt, weight_function="mlogit",
                             weightfun_pars=list(vars=1, lags=1), cond_dist="ind_Student", identification="non-Gaussianity",
                             B_constraints=matrix(c(1, NA, -1, 1) , nrow=2), alt_par=TRUE), -580.1779, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123expitngb_1_1_alt, weight_function="exponential", weightfun_pars=c(1, 1),
                             cond_dist="ind_Student", identification="non-Gaussianity",
                             B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3), alt_par=TRUE), -2782.476, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=theta_132thresitngb_1_1_alt, weight_function="threshold", weightfun_pars=c(1, 1),
                             cond_dist="ind_Student", identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2),
                             alt_par=TRUE), -3827.6, tolerance=1e-1)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222expcmwitngb_2_1_alt, weight_function="exponential", weightfun_pars=c(2, 1),
                             cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_222,
                             weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), identification="non-Gaussianity",
                             B_constraints=matrix(c(NA, -1, 0, 1), nrow=2), alt_par=TRUE), -1589.59, tolerance=1e-2)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=theta_132thresmwitngb_1_1_alt, weight_function="threshold", weightfun_pars=c(1, 1),
                             cond_dist="ind_Student", mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2)),
                             identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2), alt_par=TRUE),
               -508.1507, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123logisticcmitngb_3_1_alt, weight_function="logistic", weightfun_pars=c(3, 1),
                             cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_123, identification="non-Gaussianity",
                             B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3), alt_par=TRUE), -3988.082, tolerance=1e-2)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logcitngb_2_1_alt, weight_function="mlogit",
                             weightfun_pars=list(vars=2, lags=1), cond_dist="ind_Student", AR_constraints=C_222,
                             identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2), alt_par=TRUE),
               -2820.663, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222exoitngb_2_1_alt, weight_function="exogenous",
                             weightfun_pars=weightfun_pars222, cond_dist="ind_Student", AR_constraints=C_222,
                             identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2), alt_par=TRUE),
               -2657.869, tolerance=1e-3)
  expect_equal(c(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222exoitngb_2_1_alt, weight_function="exogenous",
                               weightfun_pars=weightfun_pars222, cond_dist="ind_Student", AR_constraints=C_222,
                               identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2), alt_par=TRUE,
                               to_return="total_ccovs")[, , 242]),
               c(0.07106374, -0.04640166, -0.04640166, 0.04615817), tolerance=1e-3)

  # ind_skewed_t
  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112ikt, cond_dist="ind_skewed_t", weight_function="threshold",
                             weightfun_pars=c(1, 1), identification="non-Gaussianity"), -851.9794, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistikt_2_1, cond_dist="ind_skewed_t", weight_function="logistic",
                             weightfun_pars=c(2, 1), identification="non-Gaussianity"), -350.3691, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistiktngb_2_1, cond_dist="ind_skewed_t", weight_function="logistic",
                             weightfun_pars=c(2, 1), identification="non-Gaussianity", B_constraints=matrix(c(NA, -1, 0, 1), nrow=2)),
               -5428.652, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122logiktngb_1_1, weight_function="mlogit",
                             weightfun_pars=list(vars=1, lags=1), cond_dist="ind_skewed_t", identification="non-Gaussianity",
                             B_constraints=matrix(c(1, NA, -1, 1) , nrow=2)), -612.9666, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123expiktngb_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                             cond_dist="ind_skewed_t", identification="non-Gaussianity",
                             B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)), -2939.885, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=theta_132thresiktngb_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                             cond_dist="ind_skewed_t", identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2)),
               -3827.6, tolerance=1e-1)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222expcmwiktngb_2_1, weight_function="exponential", weightfun_pars=c(2, 1),
                             cond_dist="ind_skewed_t", mean_constraints=list(1:2), AR_constraints=C_222,
                             weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), identification="non-Gaussianity",
                             B_constraints=matrix(c(NA, -1, 0, 1), nrow=2)),
               -1637.919, tolerance=1e-2)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=theta_132thresmwiktngb_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                             cond_dist="ind_skewed_t", mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2)),
                             identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2)),
               -508.1507, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123logisticcmiktngb_3_1, weight_function="logistic", weightfun_pars=c(3, 1),
                             cond_dist="ind_skewed_t", mean_constraints=list(1:2), AR_constraints=C_123, identification="non-Gaussianity",
                             B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3)),
               -4450.637, tolerance=1e-2)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logciktngb_2_1, weight_function="mlogit",
                             weightfun_pars=list(vars=2, lags=1), cond_dist="ind_skewed_t", AR_constraints=C_222,
                             identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)),
               -2978.483, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222exoiktngb_2_1, weight_function="exogenous",
                             weightfun_pars=weightfun_pars222, cond_dist="ind_skewed_t", AR_constraints=C_222,
                             identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2)),
               -2750.308, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222exoiktngb_2_1, weight_function="exogenous",
                             weightfun_pars=weightfun_pars222, cond_dist="ind_skewed_t", AR_constraints=C_222,
                             identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2),
                             penalized=TRUE, penalty_params=c(0.6, 0.2)),
               -2804.673, tolerance=1e-3)
  tmp_pars <- theta_222exoiktngb_2_1
  tmp_pars[9] <- 3
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=tmp_pars, weight_function="exogenous",
                             weightfun_pars=weightfun_pars222, cond_dist="ind_skewed_t", AR_constraints=C_222,
                             identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2),
                             penalized=TRUE, penalty_params=c(0.1, 0.2), allow_unstab=TRUE),
               -7895.042, tolerance=1e-3)

  expect_equal(c(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222exoiktngb_2_1, weight_function="exogenous",
                               weightfun_pars=weightfun_pars222, cond_dist="ind_skewed_t", AR_constraints=C_222,
                               identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2), to_return="total_ccovs")[, , 242]),
               c(0.07106374, -0.04640166, -0.04640166, 0.04615817), tolerance=1e-3)

  # ind_skewed_t / non-Gaussianity alt parametrization
  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112ikt,  cond_dist="ind_skewed_t", weight_function="threshold",
                             weightfun_pars=c(1, 1), alt_par=TRUE), -851.9794, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistikt_2_1_alt,
                             cond_dist="ind_skewed_t", weight_function="logistic",
                             weightfun_pars=c(2, 1), alt_par=TRUE), -350.3691, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122logikt_1_1_alt, weight_function="mlogit",
                             weightfun_pars=list(vars=1, lags=1), cond_dist="ind_skewed_t", alt_par=TRUE), -578.4793, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123expikt_1_1_alt, weight_function="exponential", weightfun_pars=c(1, 1),
                             cond_dist="ind_skewed_t", alt_par=TRUE), -4050.617, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222expcmwikt_2_1_alt, weight_function="exponential", weightfun_pars=c(2, 1),
                             cond_dist="ind_skewed_t", mean_constraints=list(1:2), AR_constraints=C_222,
                             weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), alt_par=TRUE), -391.6119, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=theta_132thresmwikt_1_1_alt, weight_function="threshold", weightfun_pars=c(1, 1),
                             cond_dist="ind_skewed_t", mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2)), alt_par=TRUE),
               -556.0658, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123logisticcmikt_3_1_alt, weight_function="logistic", weightfun_pars=c(3, 1),
                             cond_dist="ind_skewed_t", mean_constraints=list(1:2), AR_constraints=C_123, alt_par=TRUE),
               -7225.257, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logcikt_2_1_alt, weight_function="mlogit",
                             weightfun_pars=list(vars=2, lags=1), cond_dist="ind_skewed_t", AR_constraints=C_222, alt_par=TRUE),
               -403.1479, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222exoikt_2_1_alt, weight_function="exogenous",
                             weightfun_pars=weightfun_pars222, cond_dist="ind_skewed_t", AR_constraints=C_222, alt_par=TRUE),
               -422.6878, tolerance=1e-3)

  expect_equal(loglikelihood(data=gdpdef, p=1, M=1, params=theta_112ikt, cond_dist="ind_skewed_t", weight_function="threshold",
                             weightfun_pars=c(1, 1), identification="non-Gaussianity", alt_par=TRUE), -851.9794, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistikt_2_1_alt, cond_dist="ind_skewed_t", weight_function="logistic",
                             weightfun_pars=c(2, 1), identification="non-Gaussianity", alt_par=TRUE), -350.3691, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logistiktngb_2_1_alt, cond_dist="ind_skewed_t", weight_function="logistic",
                             weightfun_pars=c(2, 1), identification="non-Gaussianity", B_constraints=matrix(c(NA, -1, 0, 1), nrow=2),
                             alt_par=TRUE), -5428.652, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=2, params=theta_122logiktngb_1_1_alt, weight_function="mlogit",
                             weightfun_pars=list(vars=1, lags=1), cond_dist="ind_skewed_t", identification="non-Gaussianity",
                             B_constraints=matrix(c(1, NA, -1, 1) , nrow=2), alt_par=TRUE), -612.9666, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123expiktngb_1_1_alt, weight_function="exponential", weightfun_pars=c(1, 1),
                             cond_dist="ind_skewed_t", identification="non-Gaussianity",
                             B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3), alt_par=TRUE), -2939.885, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=theta_132thresiktngb_1_1_alt, weight_function="threshold", weightfun_pars=c(1, 1),
                             cond_dist="ind_skewed_t", identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2),
                             alt_par=TRUE), -3827.6, tolerance=1e-1)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222expcmwiktngb_2_1_alt, weight_function="exponential", weightfun_pars=c(2, 1),
                             cond_dist="ind_skewed_t", mean_constraints=list(1:2), AR_constraints=C_222,
                             weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), identification="non-Gaussianity",
                             B_constraints=matrix(c(NA, -1, 0, 1), nrow=2), alt_par=TRUE), -1637.919, tolerance=1e-2)
  expect_equal(loglikelihood(data=gdpdef, p=1, M=3, params=theta_132thresmwiktngb_1_1_alt, weight_function="threshold", weightfun_pars=c(1, 1),
                             cond_dist="ind_skewed_t", mean_constraints=list(1, 2:3), weight_constraints=list(R=0, r=c(0, 1.2)),
                             identification="non-Gaussianity", B_constraints=matrix(c(1, 0, NA, 1), nrow=2), alt_par=TRUE),
               -508.1507, tolerance=1e-3)
  expect_equal(loglikelihood(data=usamone, p=1, M=2, params=theta_123logisticcmiktngb_3_1_alt, weight_function="logistic", weightfun_pars=c(3, 1),
                             cond_dist="ind_skewed_t", mean_constraints=list(1:2), AR_constraints=C_123, identification="non-Gaussianity",
                             B_constraints=matrix(c(1, NA, NA, 0, 1, NA, 0, NA, 1), nrow=3), alt_par=TRUE), -4450.637, tolerance=1e-2)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222logciktngb_2_1_alt, weight_function="mlogit",
                             weightfun_pars=list(vars=2, lags=1), cond_dist="ind_skewed_t", AR_constraints=C_222,
                             identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2), alt_par=TRUE),
               -2978.483, tolerance=1e-3)
  expect_equal(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222exoiktngb_2_1_alt, weight_function="exogenous",
                             weightfun_pars=weightfun_pars222, cond_dist="ind_skewed_t", AR_constraints=C_222,
                             identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2), alt_par=TRUE),
               -2750.308, tolerance=1e-3)
  expect_equal(c(loglikelihood(data=gdpdef, p=2, M=2, params=theta_222exoiktngb_2_1_alt, weight_function="exogenous",
                               weightfun_pars=weightfun_pars222, cond_dist="ind_skewed_t", AR_constraints=C_222,
                               identification="non-Gaussianity", B_constraints=matrix(c(NA, NA, 0, 1), nrow=2), alt_par=TRUE,
                               to_return="total_ccovs")[, , 242]),
               c(0.07106374, -0.04640166, -0.04640166, 0.04615817), tolerance=1e-3)
})
