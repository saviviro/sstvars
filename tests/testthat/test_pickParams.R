context("pickParams")
library(sstvars)

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
A31_132 <- matrix(c(0.1, 0.2, 0.3, 0.4), nrow=2)
Omega1_132 <- Omega1_122
Omega2_132 <- Omega2_122
Omega3_132 <- matrix(c(1, 0.5, 0.5, 1), nrow=2)
alpha1_132 <- 0.5
alpha2_132 <- 0.3

theta_132relg <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                   vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), alpha1_132, alpha2_132)

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

# with cond_dist="Student"

# p=2, M=2, d=2, cond_dist="Student", weight_function="logistic", weightfun_pars=c(2, 1)
df_222_2_1 <- 3
theta_222logistict_2_1 <- c(theta_222logistic_2_1, df_222_2_1)


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

# p=1, M=2, d=3, weightfun_pars=list(vars=1, lags=1)
gamma1_123_1_1 <- c(0.1, 0.2)
theta_123log_1_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vech(Omega1_123),
                      vech(Omega2_123), gamma1_123_1_1)

# p=1, M=2, d=3, weightfun_pars=list(vars=2:3, lags=1)
gamma1_123_23_1 <- c(0.1, 0.2, 0.3)
theta_123log_23_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vech(Omega1_123),
                      vech(Omega2_123), gamma1_123_23_1)

# p=1, M=2, d=3, weightfun_pars=list(vars=1:3, lags=1)
gamma1_123_123_1 <- c(0.1, 0.2, 0.3, 0.4)
theta_123log_123_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vech(Omega1_123),
                       vech(Omega2_123), gamma1_123_123_1)

# with cond_dist="Student"

# p=1, M=2, d=3, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="Student"
df_123_1_1 <- 13
theta_123logt_1_1 <- c(theta_123log_1_1, df_123_1_1)


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

# with cond_dist="Student"

# p=1, M=2, d=2, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="Student"
df_122_1_1 <- 10
theta_122expt_1_1 <- c(theta_122exp_1_1, df_122_1_1)


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
r1_232_1_1 <- r1_132_1_1; r2_232_1_1 <- r2_132_1_1
theta_232thres_1_1 <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                        vec(A31_232), vec(A32_232), vech(Omega1_232), vech(Omega2_232), vech(Omega3_232),
                        r1_232_1_1, r2_232_1_1)


# p=1, M=2, d=3, weight_function="threshold", weightfun_pars=c(2, 1)
r1_123_2_1 <- 1
theta_123thres_2_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vech(Omega1_123),
                        vech(Omega2_123), r1_123_2_1)

# with cond_dist="Student"

# p=2, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="Student"
df_232_1_1 <- 30
theta_232threst_1_1 <- c(theta_232thres_1_1, df_232_1_1)


## weight_function == "exogenous"

# p=1, M=2, d=3, weight_function="exogenous", weightfun_pars=matrix(cbind(c(0.4, 0.2, 0.9), c(0.6, 0.8, 0.1)))
theta_123exo <- c(phi10_123, phi20_123, vec(A11_123), vec(A21_123), vech(Omega1_123), vech(Omega2_123))

# p=2, M=2, d=2,  weight_function="exogenous", weightfun_pars=matrix(cbind(c(0.4, 0, 0.9), c(0.6, 1, 0.1)))
theta_222exo <- c(phi10_222, phi20_222, vec(A11_222), vec(A12_222), vec(A21_222), vec(A22_222), vech(Omega1_222), vech(Omega2_222))

# p=1, M=3, d=2, weight_function="exogenous", weightfun_pars=matrix(cbind(c(0.4, 0, 0.5), c(0.3, 1, 0.2), c(0.3, 0, 0.3)))
theta_132exo <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132), vech(Omega1_132),
                  vech(Omega2_132), vech(Omega3_132))

# p=2, M=3, d=2, weight_function="exogenous", weightfun_pars=matrix(cbind(c(0.4, 0, 0.5), c(0.3, 1, 0.2), c(0.3, 0, 0.3)))
theta_232exo <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                  vec(A31_232), vec(A32_232), vech(Omega1_232), vech(Omega2_232), vech(Omega3_232))

# p=2, M=3, d=2, weight_function="exogenous", weightfun_pars=matrix(cbind(c(0.4, 0, 0.5), c(0.3, 1, 0.2), c(0.3, 0, 0.3))), cond_dist="Student"
theta_232exo <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                  vec(A31_232), vec(A32_232), vech(Omega1_232), vech(Omega2_232), vech(Omega3_232), df_232_1_1 )

###############
### ind_Student

# p=1, M=1, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student"
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

# p=2, M=3, d=2, weight_function="exogenous", weightfun_pars=matrix(cbind(c(0.4, 0, 0.5), c(0.3, 1, 0.2), c(0.3, 0, 0.3))),
# cond_dist="ind_Student"
dfs_232_1_1 <- c(3, 6)
B1_232 <- matrix(c(0.6, -0.2, 0.1, 0.7), nrow=2)
B2_232 <- matrix(c(0.4, 0.1, 0.2, 0.5), nrow=2)
B3_232 <- matrix(c(0.9, 0.3, -0.2, 0.4), nrow=2)
theta_232exoit <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                    vec(A31_232), vec(A32_232), vec(B1_232), vec(B2_232), vec(B3_232), dfs_232_1_1)

###############

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

### ind_Student

# p=2, M=2, d=2, cond_dist="ind_Student", weight_function="logistic", weightfun_pars=c(2, 1), identification="non-Gaussianity"
theta_222logistitng_2_1 <- theta_222logistit_2_1

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1, lags=1), cond_dist="ind_Student", identification="non-Gaussianity"
theta_122logitng_1_1 <- theta_122logit_1_1

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), cond_dist="ind_Student", identification="non-Gaussianity"
theta_123expitng_1_1 <- theta_123expit_1_1

# p=1, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="ind_Student", identification="non-Gaussianity"
theta_132thresitng_1_1 <- theta_132thresit_1_1

# p=2, M=3, d=2, weight_function="exogenous", weightfun_pars=matrix(cbind(c(0.4, 0, 0.5), c(0.3, 1, 0.2), c(0.3, 0, 0.3))),
# cond_dist="ind_Student", identification="non-Gaussianity"
theta_232exoitng <- theta_232exoit

test_that("pick_phi0 works correctly", {

  expect_equal(pick_phi0(M=2, d=3, params=theta_123exo)[,1], phi10_123)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123exo)[,2], phi20_123)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222exo)[,1], phi10_222)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222exo)[,2], phi20_222)
  expect_equal(pick_phi0(M=3, d=2, params=theta_132exo)[,1], phi10_132)
  expect_equal(pick_phi0(M=3, d=2, params=theta_132exo)[,2], phi20_132)
  expect_equal(pick_phi0(M=3, d=2, params=theta_132exo)[,3], phi30_132)
  expect_equal(pick_phi0(M=3, d=2, params=theta_232exo)[,1], phi10_232)
  expect_equal(pick_phi0(M=3, d=2, params=theta_232exo)[,2], phi20_232)
  expect_equal(pick_phi0(M=3, d=2, params=theta_232exo)[,3], phi30_232)

  expect_equal(pick_phi0(M=2, d=2, params=theta_122thres_1_1)[,1], phi10_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122thres_1_1)[,2], phi20_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222thres_2_2)[,1], phi10_222)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222thres_2_2)[,2], phi20_222)
  expect_equal(pick_phi0(M=3, d=2, params=theta_132thres_1_1)[,1], phi10_132)
  expect_equal(pick_phi0(M=3, d=2, params=theta_132thres_1_1)[,2], phi20_132)
  expect_equal(pick_phi0(M=3, d=2, params=theta_132thres_1_1)[,3], phi30_132)
  expect_equal(pick_phi0(M=3, d=2, params=theta_232thres_1_1)[,1], phi10_232)
  expect_equal(pick_phi0(M=3, d=2, params=theta_232thres_1_1)[,2], phi20_232)
  expect_equal(pick_phi0(M=3, d=2, params=theta_232thres_1_1)[,3], phi30_232)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123thres_2_1)[,1], phi10_123)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123thres_2_1)[,2], phi20_123)

  expect_equal(pick_phi0(M=2, d=2, params=theta_122exp_1_1)[,1], phi10_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122exp_1_1)[,2], phi20_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122exp_2_1)[,1], phi10_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122exp_2_1)[,2], phi20_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222exp_2_1)[,1], phi10_222)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222exp_2_1)[,2], phi20_222)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222exp_1_2)[,1], phi10_222)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222exp_1_2)[,2], phi20_222)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123exp_1_1)[,1], phi10_123)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123exp_1_1)[,2], phi20_123)

  expect_equal(pick_phi0(M=2, d=2, params=theta_122logistic_1_1)[,1], phi10_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122logistic_1_1)[,2], phi20_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122logistic_2_1)[,1], phi10_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122logistic_2_1)[,2], phi20_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222logistic_2_1)[,1], phi10_222)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222logistic_2_1)[,2], phi20_222)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222logistic_1_2)[,1], phi10_222)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222logistic_1_2)[,2], phi20_222)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123logistic_1_1)[,1], phi10_123)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123logistic_1_1)[,2], phi20_123)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123logistic_3_1)[,1], phi10_123)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123logistic_3_1)[,2], phi20_123)

  expect_equal(pick_phi0(M=1, d=2, params=theta_112relg), as.matrix(phi10_112))
  expect_equal(pick_phi0(M=1, d=2, params=theta_212relg), as.matrix(phi10_212))
  expect_equal(pick_phi0(M=1, d=2, params=theta_312relg), as.matrix(phi10_312))
  expect_equal(pick_phi0(M=2, d=2, params=theta_122relg)[,1], phi10_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122relg)[,2], phi20_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222relg)[,1], phi10_222)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222relg)[,2], phi20_222)
  expect_equal(pick_phi0(M=3, d=2, params=theta_132relg)[,1], phi10_132)
  expect_equal(pick_phi0(M=3, d=2, params=theta_132relg)[,2], phi20_132)
  expect_equal(pick_phi0(M=3, d=2, params=theta_132relg)[,3], phi30_132)
  expect_equal(pick_phi0(M=1, d=3, params=theta_113relg), as.matrix(phi10_113))
  expect_equal(pick_phi0(M=1, d=3, params=theta_213relg), as.matrix(phi10_213))
  expect_equal(pick_phi0(M=2, d=3, params=theta_123relg)[,1], phi10_123)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123relg)[,2], phi20_123)

  expect_equal(pick_phi0(M=2, d=2, params=theta_122log_12_1)[,1], phi10_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122log_12_1)[,2], phi20_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222log_2_1)[,1], phi10_222)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222log_2_1)[,2], phi20_222)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222log_12_2)[,1], phi10_222)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222log_12_2)[,2], phi20_222)
  expect_equal(pick_phi0(M=3, d=2, params=theta_232log_12_2)[,1], phi10_232)
  expect_equal(pick_phi0(M=3, d=2, params=theta_232log_12_2)[,2], phi20_232)
  expect_equal(pick_phi0(M=3, d=2, params=theta_232log_12_2)[,3], phi30_232)

  # Student
  expect_equal(pick_phi0(M=3, d=2, params=theta_232threst_1_1)[,2], phi20_232)
  expect_equal(pick_phi0(M=3, d=2, params=theta_232threst_1_1)[,3], phi30_232)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122expt_1_1)[,1], phi10_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122expt_1_1)[,2], phi20_122)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123logt_1_1)[,1], phi10_123)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123logt_1_1)[,2], phi20_123)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123logt_1_1)[,1], phi10_123)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123logt_1_1)[,2], phi20_123)

  # ind_Student
  expect_equal(pick_phi0(M=1, d=2, params=theta_112it)[,1], phi10_112)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122logit_1_1)[,1], phi10_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_122logit_1_1)[,2], phi20_122)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222logistit_2_1)[,1], phi10_222)
  expect_equal(pick_phi0(M=2, d=2, params=theta_222logistit_2_1)[,2], phi20_222)
  expect_equal(pick_phi0(M=3, d=2, params=theta_132thresit_1_1)[,1], phi10_132)
  expect_equal(pick_phi0(M=3, d=2, params=theta_132thresit_1_1)[,2], phi20_132)
  expect_equal(pick_phi0(M=3, d=2, params=theta_132thresit_1_1)[,3], phi30_132)
  expect_equal(pick_phi0(M=3, d=2, params=theta_232exoit)[,1], phi10_232)
  expect_equal(pick_phi0(M=3, d=2, params=theta_232exoit)[,2], phi20_232)
  expect_equal(pick_phi0(M=3, d=2, params=theta_232exoit)[,3], phi30_232)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123expit_1_1)[,1], phi10_123)
  expect_equal(pick_phi0(M=2, d=3, params=theta_123expit_1_1)[,2], phi20_123)
})

test_that("pick_Ami works correctly", {
  expect_equal(pick_Ami(p=1, M=1, d=2, m=1, i=1, params=theta_112relg), A11_112)
  expect_equal(pick_Ami(p=2, M=1, d=2, m=1, i=1, params=theta_212relg), A11_212)
  expect_equal(pick_Ami(p=2, M=1, d=2, m=1, i=2, params=theta_212relg), A12_212)
  expect_equal(pick_Ami(p=3, M=1, d=2, m=1, i=1, params=theta_312relg), A11_312)
  expect_equal(pick_Ami(p=3, M=1, d=2, m=1, i=2, params=theta_312relg), A12_312)
  expect_equal(pick_Ami(p=3, M=1, d=2, m=1, i=3, params=theta_312relg), A13_312)
  expect_equal(pick_Ami(p=1, M=2, d=2, m=1, i=1, params=theta_122relg), A11_122)
  expect_equal(pick_Ami(p=1, M=2, d=2, m=2, i=1, params=theta_122relg), A21_122)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=1, i=1, params=theta_222relg), A11_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=1, i=2, params=theta_222relg), A12_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=2, i=1, params=theta_222relg), A21_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=2, i=2, params=theta_222relg), A22_222)
  expect_equal(pick_Ami(p=1, M=3, d=2, m=1, i=1, params=theta_132relg), A11_132)
  expect_equal(pick_Ami(p=1, M=3, d=2, m=2, i=1, params=theta_132relg), A21_132)
  expect_equal(pick_Ami(p=1, M=3, d=2, m=3, i=1, params=theta_132relg), A31_132)
  expect_equal(pick_Ami(p=1, M=1, d=3, m=1, i=1, params=theta_113relg), A11_113)
  expect_equal(pick_Ami(p=2, M=1, d=3, m=1, i=1, params=theta_213relg), A11_213)
  expect_equal(pick_Ami(p=2, M=1, d=3, m=1, i=2, params=theta_213relg), A12_213)
  expect_equal(pick_Ami(p=1, M=2, d=3, m=1, i=1, params=theta_123relg), A11_123)
  expect_equal(pick_Ami(p=1, M=2, d=3, m=2, i=1, params=theta_123relg), A21_123)

  expect_equal(pick_Ami(p=1, M=2, d=2, m=1, i=1, params=theta_122log_1_1), A11_122)
  expect_equal(pick_Ami(p=1, M=2, d=2, m=2, i=1, params=theta_122log_1_1), A21_122)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=1, i=1, params=theta_222log_12_2), A11_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=1, i=2, params=theta_222log_12_2), A12_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=2, i=1, params=theta_222log_12_2), A21_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=2, i=2, params=theta_222log_12_2), A22_222)

  expect_equal(pick_Ami(p=1, M=2, d=2, m=1, i=1, params=theta_122logistic_1_1), A11_122)
  expect_equal(pick_Ami(p=1, M=2, d=2, m=2, i=1, params=theta_122logistic_1_1), A21_122)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=1, i=1, params=theta_222logistic_2_1), A11_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=1, i=2, params=theta_222logistic_2_1), A12_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=2, i=1, params=theta_222logistic_2_1), A21_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=2, i=2, params=theta_222logistic_2_1), A22_222)
  expect_equal(pick_Ami(p=1, M=2, d=3, m=1, i=1, params=theta_123logistic_3_1), A11_123)
  expect_equal(pick_Ami(p=1, M=2, d=3, m=2, i=1, params=theta_123logistic_3_1), A21_123)

  expect_equal(pick_Ami(p=1, M=2, d=2, m=1, i=1, params=theta_122exp_1_1), A11_122)
  expect_equal(pick_Ami(p=1, M=2, d=2, m=2, i=1, params=theta_122exp_1_1), A21_122)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=1, i=1, params=theta_222exp_2_1), A11_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=1, i=2, params=theta_222exp_2_1), A12_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=2, i=1, params=theta_222exp_2_1), A21_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=2, i=2, params=theta_222exp_2_1), A22_222)
  expect_equal(pick_Ami(p=1, M=2, d=3, m=1, i=1, params=theta_123exp_1_1), A11_123)
  expect_equal(pick_Ami(p=1, M=2, d=3, m=2, i=1, params=theta_123exp_1_1), A21_123)

  expect_equal(pick_Ami(p=2, M=2, d=2, m=1, i=1, params=theta_222thres_2_2), A11_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=1, i=2, params=theta_222thres_2_2), A12_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=2, i=1, params=theta_222thres_2_2), A21_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=2, i=2, params=theta_222thres_2_2), A22_222)
  expect_equal(pick_Ami(p=1, M=3, d=2, m=1, i=1, params=theta_132thres_1_1), A11_132)
  expect_equal(pick_Ami(p=1, M=3, d=2, m=2, i=1, params=theta_132thres_1_1), A21_132)
  expect_equal(pick_Ami(p=1, M=3, d=2, m=3, i=1, params=theta_132thres_1_1), A31_132)

  expect_equal(pick_Ami(p=2, M=2, d=2, m=1, i=1, params=theta_222exo), A11_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=1, i=2, params=theta_222exo), A12_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=2, i=1, params=theta_222exo), A21_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=2, i=2, params=theta_222exo), A22_222)
  expect_equal(pick_Ami(p=1, M=3, d=2, m=1, i=1, params=theta_132exo), A11_132)
  expect_equal(pick_Ami(p=1, M=3, d=2, m=2, i=1, params=theta_132exo), A21_132)
  expect_equal(pick_Ami(p=1, M=3, d=2, m=3, i=1, params=theta_132exo), A31_132)

  # Student
  expect_equal(pick_Ami(p=1, M=2, d=3, m=1, i=1, params=theta_123logt_1_1), A11_123)
  expect_equal(pick_Ami(p=1, M=2, d=3, m=2, i=1, params=theta_123logt_1_1), A21_123)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=1, i=1, params=theta_222logistict_2_1), A11_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=1, i=2, params=theta_222logistict_2_1), A12_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=2, i=1, params=theta_222logistict_2_1), A21_222)
  expect_equal(pick_Ami(p=2, M=2, d=2, m=2, i=2, params=theta_222logistict_2_1), A22_222)
  expect_equal(pick_Ami(p=1, M=2, d=2, m=1, i=1, params=theta_122expt_1_1), A11_122)
  expect_equal(pick_Ami(p=1, M=2, d=2, m=2, i=1, params=theta_122expt_1_1), A21_122)
  expect_equal(pick_Ami(p=2, M=3, d=2, m=1, i=1, params=theta_232threst_1_1), A11_232)
  expect_equal(pick_Ami(p=2, M=3, d=2, m=2, i=1, params=theta_232threst_1_1), A21_232)
  expect_equal(pick_Ami(p=2, M=3, d=2, m=3, i=1, params=theta_232threst_1_1), A31_232)
  expect_equal(pick_Ami(p=2, M=3, d=2, m=1, i=2, params=theta_232threst_1_1), A12_232)
  expect_equal(pick_Ami(p=2, M=3, d=2, m=2, i=2, params=theta_232threst_1_1), A22_232)
  expect_equal(pick_Ami(p=2, M=3, d=2, m=3, i=2, params=theta_232threst_1_1), A32_232)

  # ind_Student
  expect_equal(pick_Ami(p=1, M=1, d=2, m=1, i=1, params=theta_112it), A11_112)
  expect_equal(pick_Ami(p=1, M=2, d=3, m=1, i=1, params=theta_123expit_1_1), A11_123)
  expect_equal(pick_Ami(p=1, M=2, d=3, m=2, i=1, params=theta_123expit_1_1), A21_123)
  expect_equal(pick_Ami(p=2, M=3, d=2, m=1, i=1, params=theta_232exoit), A11_232)
  expect_equal(pick_Ami(p=2, M=3, d=2, m=2, i=1, params=theta_232exoit), A21_232)
  expect_equal(pick_Ami(p=2, M=3, d=2, m=3, i=1, params=theta_232exoit), A31_232)
  expect_equal(pick_Ami(p=2, M=3, d=2, m=1, i=2, params=theta_232exoit), A12_232)
  expect_equal(pick_Ami(p=2, M=3, d=2, m=2, i=2, params=theta_232exoit), A22_232)
  expect_equal(pick_Ami(p=2, M=3, d=2, m=3, i=2, params=theta_232exoit), A32_232)

  # unvec=FALSE
  expect_equal(pick_Ami(p=1, M=1, d=2, m=1, i=1, params=theta_112relg, unvec=FALSE), vec(A11_112))
  expect_equal(pick_Ami(p=1, M=2, d=2, m=1, i=1, params=theta_122relg, unvec=FALSE), vec(A11_122))
  expect_equal(pick_Ami(p=2, M=2, d=2, m=2, i=2, params=theta_222relg, unvec=FALSE), vec(A22_222))
  expect_equal(pick_Ami(p=2, M=1, d=3, m=1, i=1, params=theta_213relg, unvec=FALSE), vec(A11_213))
})

test_that("pick_Am works correctly", {
  expect_equal(pick_Am(p=1, M=1, d=2, m=1, params=theta_112relg)[, , 1], A11_112)
  expect_equal(pick_Am(p=2, M=1, d=2, m=1, params=theta_212relg)[, , 1], A11_212)
  expect_equal(pick_Am(p=2, M=1, d=2, m=1, params=theta_212relg)[, , 2], A12_212)
  expect_equal(pick_Am(p=3, M=1, d=2, m=1, params=theta_312relg)[, , 1], A11_312)
  expect_equal(pick_Am(p=3, M=1, d=2, m=1, params=theta_312relg)[, , 2], A12_312)
  expect_equal(pick_Am(p=3, M=1, d=2, m=1, params=theta_312relg)[, , 3], A13_312)
  expect_equal(pick_Am(p=1, M=2, d=2, m=1, params=theta_122relg)[, , 1], A11_122)
  expect_equal(pick_Am(p=1, M=2, d=2, m=2, params=theta_122relg)[, , 1], A21_122)
  expect_equal(pick_Am(p=2, M=2, d=2, m=1, params=theta_222relg)[, , 1], A11_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=1, params=theta_222relg)[, , 2], A12_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=2, params=theta_222relg)[, , 1], A21_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=2, params=theta_222relg)[, , 2], A22_222)
  expect_equal(pick_Am(p=1, M=3, d=2, m=1, params=theta_132relg)[, , 1], A11_132)
  expect_equal(pick_Am(p=1, M=3, d=2, m=2, params=theta_132relg)[, , 1], A21_132)
  expect_equal(pick_Am(p=1, M=3, d=2, m=3, params=theta_132relg)[, , 1], A31_132)
  expect_equal(pick_Am(p=1, M=1, d=3, m=1, params=theta_113relg)[, , 1], A11_113)
  expect_equal(pick_Am(p=2, M=1, d=3, m=1, params=theta_213relg)[, , 1], A11_213)
  expect_equal(pick_Am(p=2, M=1, d=3, m=1, params=theta_213relg)[, , 2], A12_213)
  expect_equal(pick_Am(p=1, M=2, d=3, m=1, params=theta_123relg)[, , 1], A11_123)
  expect_equal(pick_Am(p=1, M=2, d=3, m=2, params=theta_123relg)[, , 1], A21_123)

  expect_equal(pick_Am(p=1, M=2, d=2, m=1, params=theta_122log_12_1)[, , 1], A11_122)
  expect_equal(pick_Am(p=1, M=2, d=2, m=2, params=theta_122log_12_1)[, , 1], A21_122)
  expect_equal(pick_Am(p=2, M=2, d=2, m=1, params=theta_222log_2_1)[, , 1], A11_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=1, params=theta_222log_2_1)[, , 2], A12_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=2, params=theta_222log_2_1)[, , 1], A21_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=2, params=theta_222log_2_1)[, , 2], A22_222)

  expect_equal(pick_Am(p=1, M=2, d=2, m=1, params=theta_122logistic_1_1)[, , 1], A11_122)
  expect_equal(pick_Am(p=1, M=2, d=2, m=2, params=theta_122logistic_1_1)[, , 1], A21_122)
  expect_equal(pick_Am(p=2, M=2, d=2, m=1, params=theta_222logistic_2_1)[, , 1], A11_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=1, params=theta_222logistic_2_1)[, , 2], A12_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=2, params=theta_222logistic_2_1)[, , 1], A21_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=2, params=theta_222logistic_2_1)[, , 2], A22_222)
  expect_equal(pick_Am(p=1, M=2, d=3, m=1, params=theta_123logistic_1_1)[, , 1], A11_123)
  expect_equal(pick_Am(p=1, M=2, d=3, m=2, params=theta_123logistic_1_1)[, , 1], A21_123)

  expect_equal(pick_Am(p=1, M=2, d=2, m=1, params=theta_122exp_1_1)[, , 1], A11_122)
  expect_equal(pick_Am(p=1, M=2, d=2, m=2, params=theta_122exp_1_1)[, , 1], A21_122)
  expect_equal(pick_Am(p=2, M=2, d=2, m=1, params=theta_222exp_2_1)[, , 1], A11_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=1, params=theta_222exp_2_1)[, , 2], A12_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=2, params=theta_222exp_2_1)[, , 1], A21_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=2, params=theta_222exp_2_1)[, , 2], A22_222)
  expect_equal(pick_Am(p=1, M=2, d=3, m=1, params=theta_123exp_1_1)[, , 1], A11_123)
  expect_equal(pick_Am(p=1, M=2, d=3, m=2, params=theta_123exp_1_1)[, , 1], A21_123)

  expect_equal(pick_Am(p=1, M=2, d=2, m=1, params=theta_122thres_1_1)[, , 1], A11_122)
  expect_equal(pick_Am(p=1, M=2, d=2, m=2, params=theta_122thres_1_1)[, , 1], A21_122)
  expect_equal(pick_Am(p=2, M=2, d=2, m=1, params=theta_222thres_2_2)[, , 1], A11_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=1, params=theta_222thres_2_2)[, , 2], A12_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=2, params=theta_222thres_2_2)[, , 1], A21_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=2, params=theta_222thres_2_2)[, , 2], A22_222)
  expect_equal(pick_Am(p=1, M=3, d=2, m=1, params=theta_132thres_1_1)[, , 1], A11_132)
  expect_equal(pick_Am(p=1, M=3, d=2, m=2, params=theta_132thres_1_1)[, , 1], A21_132)
  expect_equal(pick_Am(p=1, M=3, d=2, m=3, params=theta_132thres_1_1)[, , 1], A31_132)

  expect_equal(pick_Am(p=1, M=2, d=3, m=1, params=theta_123exo)[, , 1], A11_123)
  expect_equal(pick_Am(p=1, M=2, d=3, m=2, params=theta_123exo)[, , 1], A21_123)
  expect_equal(pick_Am(p=2, M=2, d=2, m=1, params=theta_222exo)[, , 1], A11_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=1, params=theta_222exo)[, , 2], A12_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=2, params=theta_222exo)[, , 1], A21_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=2, params=theta_222exo)[, , 2], A22_222)
  expect_equal(pick_Am(p=1, M=3, d=2, m=1, params=theta_132exo)[, , 1], A11_132)
  expect_equal(pick_Am(p=1, M=3, d=2, m=2, params=theta_132exo)[, , 1], A21_132)
  expect_equal(pick_Am(p=1, M=3, d=2, m=3, params=theta_132exo)[, , 1], A31_132)

  # Student
  expect_equal(pick_Am(p=2, M=3, d=2, m=1, params=theta_232threst_1_1)[, , 1], A11_232)
  expect_equal(pick_Am(p=2, M=3, d=2, m=2, params=theta_232threst_1_1)[, , 1], A21_232)
  expect_equal(pick_Am(p=2, M=3, d=2, m=3, params=theta_232threst_1_1)[, , 1], A31_232)
  expect_equal(pick_Am(p=2, M=3, d=2, m=1, params=theta_232threst_1_1)[, , 2], A12_232)
  expect_equal(pick_Am(p=2, M=3, d=2, m=2, params=theta_232threst_1_1)[, , 2], A22_232)
  expect_equal(pick_Am(p=2, M=3, d=2, m=3, params=theta_232threst_1_1)[, , 2], A32_232)
  expect_equal(pick_Am(p=1, M=2, d=2, m=1, params=theta_122expt_1_1)[, , 1], A11_122)
  expect_equal(pick_Am(p=1, M=2, d=2, m=2, params=theta_122expt_1_1)[, , 1], A21_122)
  expect_equal(pick_Am(p=2, M=2, d=2, m=1, params=theta_222logistict_2_1)[, , 1], A11_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=1, params=theta_222logistict_2_1)[, , 2], A12_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=2, params=theta_222logistict_2_1)[, , 1], A21_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=2, params=theta_222logistict_2_1)[, , 2], A22_222)
  expect_equal(pick_Am(p=1, M=2, d=3, m=1, params=theta_123logt_1_1)[, , 1], A11_123)
  expect_equal(pick_Am(p=1, M=2, d=3, m=2, params=theta_123logt_1_1)[, , 1], A21_123)

  # ind_Student
  expect_equal(pick_Am(p=1, M=1, d=2, m=1, params=theta_112it)[, , 1], A11_112)
  expect_equal(pick_Am(p=2, M=3, d=2, m=1, params=theta_232exoit)[, , 1], A11_232)
  expect_equal(pick_Am(p=2, M=3, d=2, m=2, params=theta_232exoit)[, , 1], A21_232)
  expect_equal(pick_Am(p=2, M=3, d=2, m=3, params=theta_232exoit)[, , 1], A31_232)
  expect_equal(pick_Am(p=2, M=3, d=2, m=1, params=theta_232exoit)[, , 2], A12_232)
  expect_equal(pick_Am(p=2, M=3, d=2, m=2, params=theta_232exoit)[, , 2], A22_232)
  expect_equal(pick_Am(p=2, M=3, d=2, m=3, params=theta_232exoit)[, , 2], A32_232)
  expect_equal(pick_Am(p=1, M=2, d=2, m=1, params=theta_122logit_1_1)[, , 1], A11_122)
  expect_equal(pick_Am(p=1, M=2, d=2, m=2, params=theta_122logit_1_1)[, , 1], A21_122)
  expect_equal(pick_Am(p=2, M=2, d=2, m=1, params=theta_222logistit_2_1)[, , 1], A11_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=1, params=theta_222logistit_2_1)[, , 2], A12_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=2, params=theta_222logistit_2_1)[, , 1], A21_222)
  expect_equal(pick_Am(p=2, M=2, d=2, m=2, params=theta_222logistit_2_1)[, , 2], A22_222)
  expect_equal(pick_Am(p=1, M=2, d=3, m=1, params=theta_123expit_1_1)[, , 1], A11_123)
  expect_equal(pick_Am(p=1, M=2, d=3, m=2, params=theta_123expit_1_1)[, , 1], A21_123)
})

test_that("pick_allA works correctly", {
  expect_equal(pick_allA(p=1, M=1, d=2, params=theta_112relg)[, , 1, 1], A11_112)
  expect_equal(pick_allA(p=2, M=1, d=2, params=theta_212relg)[, , 1, 1], A11_212)
  expect_equal(pick_allA(p=2, M=1, d=2, params=theta_212relg)[, , 2, 1], A12_212)
  expect_equal(pick_allA(p=3, M=1, d=2, params=theta_312relg)[, , 1, 1], A11_312)
  expect_equal(pick_allA(p=3, M=1, d=2, params=theta_312relg)[, , 2, 1], A12_312)
  expect_equal(pick_allA(p=3, M=1, d=2, params=theta_312relg)[, , 3, 1], A13_312)
  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122relg)[, , 1, 1], A11_122)
  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122relg)[, , 1, 2], A21_122)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222relg)[, , 1, 1], A11_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222relg)[, , 2, 1], A12_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222relg)[, , 1, 2], A21_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222relg)[, , 2, 2], A22_222)
  expect_equal(pick_allA(p=1, M=3, d=2, params=theta_132relg)[, , 1, 2], A11_132)
  expect_equal(pick_allA(p=1, M=3, d=2, params=theta_132relg)[, , 1, 2], A21_132)
  expect_equal(pick_allA(p=1, M=3, d=2, params=theta_132relg)[, , 1, 3], A31_132)
  expect_equal(pick_allA(p=1, M=1, d=3, params=theta_113relg)[, , 1, 1], A11_113)
  expect_equal(pick_allA(p=2, M=1, d=3, params=theta_213relg)[, , 1, 1], A11_213)
  expect_equal(pick_allA(p=2, M=1, d=3, params=theta_213relg)[, , 2, 1], A12_213)
  expect_equal(pick_allA(p=1, M=2, d=3, params=theta_123relg)[, , 1, 1], A11_123)
  expect_equal(pick_allA(p=1, M=2, d=3, params=theta_123relg)[, , 1, 2], A21_123)

  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122log_1_1)[, , 1, 1], A11_122)
  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122log_1_1)[, , 1, 2], A21_122)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222log_12_2)[, , 1, 1], A11_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222log_12_2)[, , 2, 1], A12_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222log_12_2)[, , 1, 2], A21_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222log_12_2)[, , 2, 2], A22_222)

  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122logistic_1_1)[, , 1, 1], A11_122)
  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122logistic_1_1)[, , 1, 2], A21_122)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222logistic_1_2)[, , 1, 1], A11_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222logistic_1_2)[, , 2, 1], A12_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222logistic_1_2)[, , 1, 2], A21_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222logistic_1_2)[, , 2, 2], A22_222)
  expect_equal(pick_allA(p=1, M=2, d=3, params=theta_123logistic_1_1)[, , 1, 1], A11_123)
  expect_equal(pick_allA(p=1, M=2, d=3, params=theta_123logistic_1_1)[, , 1, 2], A21_123)

  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122exp_1_1)[, , 1, 1], A11_122)
  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122exp_1_1)[, , 1, 2], A21_122)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222exp_1_2)[, , 1, 1], A11_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222exp_1_2)[, , 2, 1], A12_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222exp_1_2)[, , 1, 2], A21_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222exp_1_2)[, , 2, 2], A22_222)
  expect_equal(pick_allA(p=1, M=2, d=3, params=theta_123exp_1_1)[, , 1, 1], A11_123)
  expect_equal(pick_allA(p=1, M=2, d=3, params=theta_123exp_1_1)[, , 1, 2], A21_123)

  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122thres_1_1)[, , 1, 1], A11_122)
  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122thres_1_1)[, , 1, 2], A21_122)
  expect_equal(pick_allA(p=1, M=3, d=2, params=theta_132thres_1_1)[, , 1, 2], A11_132)
  expect_equal(pick_allA(p=1, M=3, d=2, params=theta_132thres_1_1)[, , 1, 2], A21_132)
  expect_equal(pick_allA(p=1, M=3, d=2, params=theta_132thres_1_1)[, , 1, 3], A31_132)

  expect_equal(pick_allA(p=1, M=2, d=3, params=theta_123exo)[, , 1, 1], A11_123)
  expect_equal(pick_allA(p=1, M=2, d=3, params=theta_123exo)[, , 1, 2], A21_123)
  expect_equal(pick_allA(p=1, M=3, d=2, params=theta_132exo)[, , 1, 2], A11_132)
  expect_equal(pick_allA(p=1, M=3, d=2, params=theta_132exo)[, , 1, 2], A21_132)
  expect_equal(pick_allA(p=1, M=3, d=2, params=theta_132exo)[, , 1, 3], A31_132)

  # Student
  expect_equal(pick_allA(p=2, M=3, d=2, params=theta_232threst_1_1)[, , 1, 2], A11_232)
  expect_equal(pick_allA(p=2, M=3, d=2, params=theta_232threst_1_1)[, , 1, 2], A21_232)
  expect_equal(pick_allA(p=2, M=3, d=2, params=theta_232threst_1_1)[, , 1, 3], A31_232)
  expect_equal(pick_allA(p=2, M=3, d=2, params=theta_232threst_1_1)[, , 2, 2], A12_232)
  expect_equal(pick_allA(p=2, M=3, d=2, params=theta_232threst_1_1)[, , 2, 2], A22_232)
  expect_equal(pick_allA(p=2, M=3, d=2, params=theta_232threst_1_1)[, , 2, 3], A32_232)
  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122expt_1_1)[, , 1, 1], A11_122)
  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122expt_1_1)[, , 1, 2], A21_122)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222logistict_2_1)[, , 1, 1], A11_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222logistict_2_1)[, , 2, 1], A12_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222logistict_2_1)[, , 1, 2], A21_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222logistict_2_1)[, , 2, 2], A22_222)
  expect_equal(pick_allA(p=1, M=2, d=3, params=theta_123log_1_1)[, , 1, 1], A11_123)
  expect_equal(pick_allA(p=1, M=2, d=3, params=theta_123log_1_1)[, , 1, 2], A21_123)

  # ind_Student
  expect_equal(pick_allA(p=1, M=1, d=2, params=theta_112it)[, , 1, 1], A11_112)
  expect_equal(pick_allA(p=2, M=3, d=2, params=theta_232exoit)[, , 1, 2], A11_232)
  expect_equal(pick_allA(p=2, M=3, d=2, params=theta_232exoit)[, , 1, 2], A21_232)
  expect_equal(pick_allA(p=2, M=3, d=2, params=theta_232exoit)[, , 1, 3], A31_232)
  expect_equal(pick_allA(p=2, M=3, d=2, params=theta_232exoit)[, , 2, 2], A12_232)
  expect_equal(pick_allA(p=2, M=3, d=2, params=theta_232exoit)[, , 2, 2], A22_232)
  expect_equal(pick_allA(p=2, M=3, d=2, params=theta_232exoit)[, , 2, 3], A32_232)
  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122logit_1_1)[, , 1, 1], A11_122)
  expect_equal(pick_allA(p=1, M=2, d=2, params=theta_122logit_1_1)[, , 1, 2], A21_122)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222logistit_2_1)[, , 1, 1], A11_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222logistit_2_1)[, , 2, 1], A12_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222logistit_2_1)[, , 1, 2], A21_222)
  expect_equal(pick_allA(p=2, M=2, d=2, params=theta_222logistit_2_1)[, , 2, 2], A22_222)
  expect_equal(pick_allA(p=1, M=2, d=3, params=theta_123expit_1_1)[, , 1, 1], A11_123)
  expect_equal(pick_allA(p=1, M=2, d=3, params=theta_123expit_1_1)[, , 1, 2], A21_123)
})

test_that("pick_Omegas works correctly", {
  expect_equal(pick_Omegas(p=1, M=1, d=2, params=theta_112relg)[, , 1], Omega1_112)
  expect_equal(pick_Omegas(p=2, M=1, d=2, params=theta_212relg)[, , 1], Omega1_212)
  expect_equal(pick_Omegas(p=3, M=1, d=2, params=theta_312relg)[, , 1], Omega1_312)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122relg)[, , 1], Omega1_122)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122relg)[, , 2], Omega2_122)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222relg)[, , 1], Omega1_222)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222relg)[, , 2], Omega2_222)
  expect_equal(pick_Omegas(p=1, M=3, d=2, params=theta_132relg)[, , 1], Omega1_132)
  expect_equal(pick_Omegas(p=1, M=3, d=2, params=theta_132relg)[, , 2], Omega2_132)
  expect_equal(pick_Omegas(p=1, M=3, d=2, params=theta_132relg)[, , 3], Omega3_132)
  expect_equal(pick_Omegas(p=1, M=1, d=3, params=theta_113relg)[, , 1], Omega1_113)
  expect_equal(pick_Omegas(p=2, M=1, d=3, params=theta_213relg)[, , 1], Omega1_213)
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123relg)[, , 1], Omega1_123)
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123relg)[, , 2], Omega2_123)

  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122log_12_1)[, , 1], Omega1_122)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122log_12_1)[, , 2], Omega2_122)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222log_2_1)[, , 1], Omega1_222)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222log_2_1)[, , 2], Omega2_222)

  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122logistic_1_1)[, , 1], Omega1_122)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122logistic_1_1)[, , 2], Omega2_122)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222logistic_2_1)[, , 1], Omega1_222)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222logistic_2_1)[, , 2], Omega2_222)
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123logistic_3_1)[, , 1], Omega1_123)
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123logistic_3_1)[, , 2], Omega2_123)

  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122exp_1_1)[, , 1], Omega1_122)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122exp_1_1)[, , 2], Omega2_122)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222exp_2_1)[, , 1], Omega1_222)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222exp_2_1)[, , 2], Omega2_222)
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123exp_1_1)[, , 1], Omega1_123)
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123exp_1_1)[, , 2], Omega2_123)

  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122thres_1_1)[, , 1], Omega1_122)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122thres_1_1)[, , 2], Omega2_122)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222thres_2_2)[, , 1], Omega1_222)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222thres_2_2)[, , 2], Omega2_222)
  expect_equal(pick_Omegas(p=1, M=3, d=2, params=theta_132thres_1_1)[, , 1], Omega1_132)
  expect_equal(pick_Omegas(p=1, M=3, d=2, params=theta_132thres_1_1)[, , 2], Omega2_132)
  expect_equal(pick_Omegas(p=1, M=3, d=2, params=theta_132thres_1_1)[, , 3], Omega3_132)

  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123exo)[, , 1], Omega1_123)
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123exo)[, , 2], Omega2_123)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222exo)[, , 1], Omega1_222)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222exo)[, , 2], Omega2_222)
  expect_equal(pick_Omegas(p=1, M=3, d=2, params=theta_132exo, cond_dist="Student")[, , 1], Omega1_132)
  expect_equal(pick_Omegas(p=1, M=3, d=2, params=theta_132exo, cond_dist="Student")[, , 2], Omega2_132)
  expect_equal(pick_Omegas(p=1, M=3, d=2, params=theta_132exo, cond_dist="Student")[, , 3], Omega3_132)

  # Student
  expect_equal(pick_Omegas(p=2, M=3, d=2, params=theta_232threst_1_1, cond_dist="Student")[, , 1], Omega1_232)
  expect_equal(pick_Omegas(p=2, M=3, d=2, params=theta_232threst_1_1, cond_dist="Student")[, , 2], Omega2_232)
  expect_equal(pick_Omegas(p=2, M=3, d=2, params=theta_232threst_1_1, cond_dist="Student")[, , 3], Omega3_232)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122expt_1_1, cond_dist="Student")[, , 1], Omega1_122)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122expt_1_1, cond_dist="Student")[, , 2], Omega2_122)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222logistict_2_1, cond_dist="Student")[, , 1], Omega1_222)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222logistict_2_1, cond_dist="Student")[, , 2], Omega2_222)
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123logt_1_1, cond_dist="Student")[, , 1], Omega1_123)
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123logt_1_1, cond_dist="Student")[, , 2], Omega2_123)

  # ind_Student
  expect_equal(pick_Omegas(p=1, M=1, d=2, params=theta_112it, cond_dist="ind_Student")[, , 1], B1_112)
  expect_equal(pick_Omegas(p=2, M=3, d=2, params=theta_232exoit, cond_dist="ind_Student")[, , 1], B1_232)
  expect_equal(pick_Omegas(p=2, M=3, d=2, params=theta_232exoit, cond_dist="ind_Student")[, , 2], B2_232)
  expect_equal(pick_Omegas(p=2, M=3, d=2, params=theta_232exoit, cond_dist="ind_Student")[, , 3], B3_232)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122logit_1_1, cond_dist="ind_Student")[, , 1], B1_122)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122logit_1_1, cond_dist="ind_Student")[, , 2], B2_122)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222logistit_2_1, cond_dist="ind_Student")[, , 1], B1_222)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222logistit_2_1, cond_dist="ind_Student")[, , 2], B2_222)
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123expit_1_1, cond_dist="ind_Student")[, , 1], B1_123)
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123expit_1_1, cond_dist="ind_Student")[, , 2], B2_123)
  expect_equal(pick_Omegas(p=1, M=3, d=2, params=theta_132thresit_1_1, cond_dist="ind_Student")[, , 1], B1_132)
  expect_equal(pick_Omegas(p=1, M=3, d=2, params=theta_132thresit_1_1, cond_dist="ind_Student")[, , 2], B2_132)
  expect_equal(pick_Omegas(p=1, M=3, d=2, params=theta_132thresit_1_1, cond_dist="ind_Student")[, , 3], B3_132)

  # Structural
  expect_equal(pick_Omegas(p=2, M=3, d=2, params=theta_232threst_1_1, identification="recursive")[, , 1], Omega1_232)
  expect_equal(pick_Omegas(p=2, M=3, d=2, params=theta_232threst_1_1, identification="recursive")[, , 2], Omega2_232)
  expect_equal(pick_Omegas(p=2, M=3, d=2, params=theta_232threst_1_1, identification="recursive")[, , 3], Omega3_232)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122relgsh, identification="heteroskedasticity")[, , 1],
               tcrossprod(W_122), tolerance=1e-4)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122relgsh, identification="heteroskedasticity")[, , 2],
               W_122%*%tcrossprod(diag(lambdas_122), W_122), tolerance=1e-4)
  expect_equal(pick_Omegas(p=1, M=3, d=2, params=theta_132relgsh, identification="heteroskedasticity")[, , 1],
               tcrossprod(W_132), tolerance=1e-4)
  expect_equal(pick_Omegas(p=1, M=3, d=2, params=theta_132relgsh, identification="heteroskedasticity")[, , 2],
               W_132%*%tcrossprod(diag(lambdas2_132), W_132), tolerance=1e-4)
  expect_equal(pick_Omegas(p=1, M=3, d=2, params=theta_132relgsh, identification="heteroskedasticity")[, , 3],
               W_132%*%tcrossprod(diag(lambdas3_132), W_132), tolerance=1e-4)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222logistictsh_2_1, identification="heteroskedasticity")[, , 1],
               tcrossprod(W_222), tolerance=1e-4)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222logistictsh_2_1, identification="heteroskedasticity")[, , 2],
               W_222%*%tcrossprod(diag(lambdas_222), W_222), tolerance=1e-4)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122logsh_12_1, identification="heteroskedasticity")[, , 1],
               tcrossprod(W_122), tolerance=1e-4)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122logsh_12_1, identification="heteroskedasticity")[, , 2],
               W_122%*%tcrossprod(diag(lambdas_122), W_122), tolerance=1e-4)
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123expsh_1_1, identification="heteroskedasticity")[, , 1],
               tcrossprod(W_123), tolerance=1e-4)
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123expsh_1_1, identification="heteroskedasticity")[, , 2],
               W_123%*%tcrossprod(diag(lambdas_123), W_123), tolerance=1e-4)
  expect_equal(pick_Omegas(p=2, M=3, d=2, params=theta_232threstsh_1_1, identification="heteroskedasticity")[, , 1],
               tcrossprod(W_232), tolerance=1e-4)
  expect_equal(pick_Omegas(p=2, M=3, d=2, params=theta_232threstsh_1_1, identification="heteroskedasticity")[, , 2],
               W_232%*%tcrossprod(diag(lambdas2_232), W_232), tolerance=1e-4)
  expect_equal(pick_Omegas(p=2, M=3, d=2, params=theta_232threstsh_1_1, identification="heteroskedasticity")[, , 3],
               W_232%*%tcrossprod(diag(lambdas3_232), W_232), tolerance=1e-4)

  # Ind student structural
  expect_equal(pick_Omegas(p=1, M=1, d=2, params=theta_112it, cond_dist="ind_Student", identification="non-Gaussianity")[, , 1], B1_112)
  expect_equal(pick_Omegas(p=2, M=3, d=2, params=theta_232exoitng, cond_dist="ind_Student", identification="non-Gaussianity")[, , 1], B1_232)
  expect_equal(pick_Omegas(p=2, M=3, d=2, params=theta_232exoitng, cond_dist="ind_Student", identification="non-Gaussianity")[, , 2], B2_232)
  expect_equal(pick_Omegas(p=2, M=3, d=2, params=theta_232exoitng, cond_dist="ind_Student", identification="non-Gaussianity")[, , 3], B3_232)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122logitng_1_1, cond_dist="ind_Student", identification="non-Gaussianity")[, , 1], B1_122)
  expect_equal(pick_Omegas(p=1, M=2, d=2, params=theta_122logitng_1_1, cond_dist="ind_Student", identification="non-Gaussianity")[, , 2], B2_122)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222logistitng_2_1, cond_dist="ind_Student", identification="non-Gaussianity")[, , 1], B1_222)
  expect_equal(pick_Omegas(p=2, M=2, d=2, params=theta_222logistitng_2_1, cond_dist="ind_Student", identification="non-Gaussianity")[, , 2], B2_222)
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123expitng_1_1, cond_dist="ind_Student", identification="non-Gaussianity")[, , 1], B1_123)
  expect_equal(pick_Omegas(p=1, M=2, d=3, params=theta_123expitng_1_1, cond_dist="ind_Student", identification="non-Gaussianity")[, , 2], B2_123)
  expect_equal(pick_Omegas(p=1, M=3, d=2, params=theta_132thresitng_1_1, cond_dist="ind_Student", identification="non-Gaussianity")[, , 1], B1_132)
  expect_equal(pick_Omegas(p=1, M=3, d=2, params=theta_132thresitng_1_1, cond_dist="ind_Student", identification="non-Gaussianity")[, , 2], B2_132)
  expect_equal(pick_Omegas(p=1, M=3, d=2, params=theta_132thresitng_1_1, cond_dist="ind_Student", identification="non-Gaussianity")[, , 3], B3_132)
})

test_that("pick_weightpars works correctly", {
  # exogenous
  expect_equal(pick_weightpars(p=1, M=2, d=3, params=theta_123exo, weight_function="exogenous", cond_dist="Gaussian",
                               weightfun_pars=matrix(cbind(c(0.4, 0.2, 0.9), c(0.6, 0.8, 0.1)))), numeric(0))
  expect_equal(pick_weightpars(p=2, M=2, d=2, params=theta_222exo, weight_function="exogenous", cond_dist="Gaussian",
                               weightfun_pars=matrix(cbind(c(0.4, 0, 0.9), c(0.6, 1, 0.1)))), numeric(0))
  expect_equal(pick_weightpars(p=1, M=3, d=2, params=theta_132exo, weight_function="exogenous", cond_dist="Gaussian",
                               weightfun_pars=matrix(cbind(c(0.4, 0, 0.5), c(0.3, 1, 0.2), c(0.3, 0, 0.3)))), numeric(0))
  expect_equal(pick_weightpars(p=2, M=3, d=2, params=theta_232exo, weight_function="exogenous", cond_dist="Student",
                               weightfun_pars=matrix(cbind(c(0.4, 0, 0.5), c(0.3, 1, 0.2), c(0.3, 0, 0.3)))), numeric(0))
  expect_equal(pick_weightpars(p=2, M=3, d=2, params=theta_232exoit, weight_function="exogenous", cond_dist="ind_Student",
                               weightfun_pars=matrix(cbind(c(0.4, 0, 0.5), c(0.3, 1, 0.2), c(0.3, 0, 0.3)))), numeric(0))

  # threshold
  expect_equal(pick_weightpars(p=1, M=2, d=2, params=theta_122thres_1_1, weight_function="threshold", cond_dist="Gaussian",
                               weightfun_pars=c(1, 1)), r1_122_1_1)
  expect_equal(pick_weightpars(p=2, M=2, d=2, params=theta_222thres_2_2, weight_function="threshold", cond_dist="Gaussian",
                               weightfun_pars=c(2, 2)), r1_222_2_2)
  expect_equal(pick_weightpars(p=1, M=3, d=2, params=theta_132thres_1_1, weight_function="threshold", cond_dist="Gaussian",
                               weightfun_pars=c(1, 1)), c(r1_132_1_1, r2_132_1_1))
  expect_equal(pick_weightpars(p=2, M=3, d=2, params=theta_232thres_1_1, weight_function="threshold", cond_dist="Gaussian",
                               weightfun_pars=c(1, 1)), c(r1_232_1_1, r2_232_1_1))
  expect_equal(pick_weightpars(p=1, M=2, d=3, params=theta_123thres_2_1, weight_function="threshold", cond_dist="Gaussian",
                               weightfun_pars=c(2, 1)), r1_123_2_1)
  expect_equal(pick_weightpars(p=2, M=3, d=2, params=theta_232threst_1_1, weight_function="threshold", cond_dist="Student",
                               weightfun_pars=c(1, 1)), c(r1_232_1_1, r2_232_1_1))
  expect_equal(pick_weightpars(p=1, M=3, d=2, params=theta_132thresit_1_1, weight_function="threshold", cond_dist="ind_Student",
                               weightfun_pars=c(1, 1)), c(r1_132_1_1, r2_132_1_1))
  expect_equal(pick_weightpars(p=1, M=1, d=2, params=theta_122it, weight_function="threshold", cond_dist="ind_Student",
                               weightfun_pars=c(1, 1)), numeric(0))

  # exponential
  expect_equal(pick_weightpars(p=1, M=2, d=2, params=theta_122exp_1_1, weight_function="exponential", cond_dist="Gaussian",
                               weightfun_pars=c(1, 1)), c_and_gamma_122_1_1)
  expect_equal(pick_weightpars(p=1, M=2, d=2, params=theta_122exp_2_1, weight_function="exponential", cond_dist="Gaussian",
                               weightfun_pars=c(2, 1)), c_and_gamma_122_2_1)
  expect_equal(pick_weightpars(p=2, M=2, d=2, params=theta_222exp_2_1, weight_function="exponential", cond_dist="Gaussian",
                               weightfun_pars=c(2, 1)), c_and_gamma_222_2_1)
  expect_equal(pick_weightpars(p=2, M=2, d=2, params=theta_222exp_1_2, weight_function="exponential", cond_dist="Gaussian",
                               weightfun_pars=c(1, 2)), c_and_gamma_222_1_2)
  expect_equal(pick_weightpars(p=1, M=2, d=3, params=theta_123exp_1_1, weight_function="exponential", cond_dist="Gaussian",
                               weightfun_pars=c(1, 1)), c_and_gamma_123_1_1)
  expect_equal(pick_weightpars(p=1, M=2, d=2, params=theta_122expt_1_1, weight_function="exponential", cond_dist="Student",
                               weightfun_pars=c(1, 1)), c_and_gamma_122_1_1)
  expect_equal(pick_weightpars(p=1, M=2, d=3, params=theta_123expit_1_1, weight_function="exponential", cond_dist="ind_Student",
                               weightfun_pars=c(1, 1)), c_and_gamma_123_1_1)

  # logistic
  expect_equal(pick_weightpars(p=1, M=2, d=2, params=theta_122logistic_1_1, weight_function="logistic", cond_dist="Gaussian",
                               weightfun_pars=c(1, 1)), c_and_gamma_122_1_1)
  expect_equal(pick_weightpars(p=1, M=2, d=2, params=theta_122logistic_2_1, weight_function="logistic", cond_dist="Gaussian",
                               weightfun_pars=c(2, 1)), c_and_gamma_122_2_1)
  expect_equal(pick_weightpars(p=2, M=2, d=2, params=theta_222logistic_2_1, weight_function="logistic", cond_dist="Gaussian",
                               weightfun_pars=c(2, 1)), c_and_gamma_222_2_1)
  expect_equal(pick_weightpars(p=2, M=2, d=2, params=theta_222logistic_1_2, weight_function="logistic", cond_dist="Gaussian",
                               weightfun_pars=c(1, 2)), c_and_gamma_222_1_2)
  expect_equal(pick_weightpars(p=1, M=2, d=3, params=theta_123logistic_1_1, weight_function="logistic", cond_dist="Gaussian",
                               weightfun_pars=c(1, 1)), c_and_gamma_123_1_1)
  expect_equal(pick_weightpars(p=1, M=2, d=3, params=theta_123logistic_3_1, weight_function="logistic", cond_dist="Gaussian",
                               weightfun_pars=c(3, 1)), c_and_gamma_123_3_1)
  expect_equal(pick_weightpars(p=2, M=2, d=2, params=theta_222logistict_2_1, weight_function="logistic", cond_dist="Student",
                               weightfun_pars=c(2, 1)), c_and_gamma_222_2_1)
  expect_equal(pick_weightpars(p=2, M=2, d=2, params=theta_222logistit_2_1, weight_function="logistic", cond_dist="ind_Student",
                               weightfun_pars=c(2, 1)), c_and_gamma_222_2_1)

  # mlogit
  expect_equal(pick_weightpars(p=1, M=2, d=2, params=theta_122log_1_1, weight_function="mlogit", cond_dist="Gaussian",
                               weightfun_pars=list(vars=1, lags=1)), gamma1_122_1_1)
  expect_equal(pick_weightpars(p=1, M=2, d=2, params=theta_122log_12_1, weight_function="mlogit", cond_dist="Gaussian",
                               weightfun_pars=list(vars=1:2, lags=1)), gamma1_122_12_1)
  expect_equal(pick_weightpars(p=2, M=2, d=2, params=theta_222log_2_1, weight_function="mlogit", cond_dist="Gaussian",
                               weightfun_pars=list(vars=2, lags=1)), gamma1_222_2_1)
  expect_equal(pick_weightpars(p=2, M=2, d=2, params=theta_222log_12_2, weight_function="mlogit", cond_dist="Gaussian",
                               weightfun_pars=list(vars=1:2, lags=2)), gamma1_222_12_2)
  expect_equal(pick_weightpars(p=2, M=3, d=2, params=theta_232log_1_1, weight_function="mlogit", cond_dist="Gaussian",
                               weightfun_pars=list(vars=1, lags=1)), c(gamma1_232_1_1, gamma2_232_1_1))
  expect_equal(pick_weightpars(p=2, M=3, d=2, params=theta_232log_12_1, weight_function="mlogit", cond_dist="Gaussian",
                               weightfun_pars=list(vars=1:2, lags=1)), c(gamma1_232_12_1, gamma2_232_12_1))
  expect_equal(pick_weightpars(p=2, M=3, d=2, params=theta_232log_2_2, weight_function="mlogit", cond_dist="Gaussian",
                               weightfun_pars=list(vars=2, lags=2)), c(gamma1_232_2_2, gamma2_232_2_2))
  expect_equal(pick_weightpars(p=2, M=3, d=2, params=theta_232log_12_2, weight_function="mlogit", cond_dist="Gaussian",
                               weightfun_pars=list(vars=1:2, lags=2)), c(gamma1_232_12_2, gamma2_232_12_2))
  expect_equal(pick_weightpars(p=1, M=2, d=3, params=theta_123log_1_1, weight_function="mlogit", cond_dist="Gaussian",
                               weightfun_pars=list(vars=1, lags=1)), gamma1_123_1_1)
  expect_equal(pick_weightpars(p=1, M=2, d=3, params=theta_123log_23_1, weight_function="mlogit", cond_dist="Gaussian",
                               weightfun_pars=list(vars=2:3, lags=1)), gamma1_123_23_1)
  expect_equal(pick_weightpars(p=1, M=2, d=3, params=theta_123log_123_1, weight_function="mlogit", cond_dist="Gaussian",
                               weightfun_pars=list(vars=1:3, lags=1)), gamma1_123_123_1)
  expect_equal(pick_weightpars(p=1, M=2, d=3, params=theta_123logt_1_1, weight_function="mlogit", cond_dist="Student",
                               weightfun_pars=list(vars=1, lags=1)), gamma1_123_1_1)
  expect_equal(pick_weightpars(p=1, M=2, d=2, params=theta_122logit_1_1, weight_function="mlogit", cond_dist="ind_Student",
                               weightfun_pars=list(vars=1, lags=1)), gamma1_122_1_1)

  # relative_dens
  expect_equal(pick_weightpars(p=1, M=1, d=2, params=theta_112relg, weight_function="relative_dens", cond_dist="Gaussian"), 1)
  expect_equal(pick_weightpars(p=2, M=1, d=2, params=theta_212relg, weight_function="relative_dens", cond_dist="Gaussian"), 1)
  expect_equal(pick_weightpars(p=3, M=1, d=2, params=theta_312relg, weight_function="relative_dens", cond_dist="Gaussian"), 1)
  expect_equal(pick_weightpars(p=1, M=2, d=2, params=theta_122relg, weight_function="relative_dens", cond_dist="Gaussian"),
               c(alpha1_122, 1 - alpha1_122))
  expect_equal(pick_weightpars(p=2, M=2, d=2, params=theta_222relg, weight_function="relative_dens", cond_dist="Gaussian"),
               c(alpha1_222, 1 - alpha1_222))
  expect_equal(pick_weightpars(p=1, M=3, d=2, params=theta_132relg, weight_function="relative_dens", cond_dist="Gaussian"),
               c(alpha1_132, alpha2_132, 1 - alpha1_132 - alpha2_132))
  expect_equal(pick_weightpars(p=1, M=1, d=3, params=theta_113relg, weight_function="relative_dens", cond_dist="Gaussian"), 1)
  expect_equal(pick_weightpars(p=2, M=1, d=3, params=theta_213relg, weight_function="relative_dens", cond_dist="Gaussian"), 1)
  expect_equal(pick_weightpars(p=1, M=2, d=3, params=theta_123relg, weight_function="relative_dens", cond_dist="Gaussian"),
               c(alpha1_123, 1 - alpha1_123))
})


test_that("pick_regime works correctly", {
  expect_equal(pick_regime(p=1, M=1, d=2, m=1, params=theta_112relg), c(phi10_112, vec(A11_112), vech(Omega1_112)))
  expect_equal(pick_regime(p=2, M=1, d=2, m=1, params=theta_212relg), c(phi10_212, vec(A11_212), vec(A12_212), vech(Omega1_212)))
  expect_equal(pick_regime(p=3, M=1, d=2, m=1, params=theta_312relg), c(phi10_312, vec(A11_312), vec(A12_312),
                                                                        vec(A13_312), vech(Omega1_312)))
  expect_equal(pick_regime(p=1, M=2, d=2, m=1, params=theta_122relg), c(phi10_122, vec(A11_122), vech(Omega1_122)))
  expect_equal(pick_regime(p=1, M=2, d=2, m=2, params=theta_122relg), c(phi20_122, vec(A21_122), vech(Omega2_122)))
  expect_equal(pick_regime(p=2, M=2, d=2, m=1, params=theta_222relg), c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222)))
  expect_equal(pick_regime(p=2, M=2, d=2, m=2, params=theta_222relg), c(phi20_222, vec(A21_222), vec(A22_222), vech(Omega2_222)))
  expect_equal(pick_regime(p=1, M=3, d=2, m=1, params=theta_132relg), c(phi10_132, vec(A11_132), vech(Omega1_132)))
  expect_equal(pick_regime(p=1, M=3, d=2, m=2, params=theta_132relg), c(phi20_132, vec(A21_132), vech(Omega2_132)))
  expect_equal(pick_regime(p=1, M=3, d=2, m=3, params=theta_132relg), c(phi30_132, vec(A31_132), vech(Omega3_132)))
  expect_equal(pick_regime(p=1, M=1, d=3, m=1, params=theta_113relg), c(phi10_113, vec(A11_113), vech(Omega1_113)))
  expect_equal(pick_regime(p=2, M=1, d=3, m=1, params=theta_213relg), c(phi10_213, vec(A11_213), vec(A12_213), vech(Omega1_213)))
  expect_equal(pick_regime(p=1, M=2, d=3, m=1, params=theta_123relg), c(phi10_123, vec(A11_123), vech(Omega1_123)))
  expect_equal(pick_regime(p=1, M=2, d=3, m=2, params=theta_123relg), c(phi20_123, vec(A21_123), vech(Omega2_123)))

  expect_equal(pick_regime(p=1, M=2, d=2, m=1, params=theta_122log_1_1), c(phi10_122, vec(A11_122), vech(Omega1_122)))
  expect_equal(pick_regime(p=1, M=2, d=2, m=2, params=theta_122log_1_1), c(phi20_122, vec(A21_122), vech(Omega2_122)))
  expect_equal(pick_regime(p=2, M=2, d=2, m=1, params=theta_222log_12_2), c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222)))
  expect_equal(pick_regime(p=2, M=2, d=2, m=2, params=theta_222log_12_2), c(phi20_222, vec(A21_222), vec(A22_222), vech(Omega2_222)))

  expect_equal(pick_regime(p=1, M=2, d=2, m=1, params=theta_122logistic_1_1), c(phi10_122, vec(A11_122), vech(Omega1_122)))
  expect_equal(pick_regime(p=1, M=2, d=2, m=2, params=theta_122logistic_1_1), c(phi20_122, vec(A21_122), vech(Omega2_122)))
  expect_equal(pick_regime(p=2, M=2, d=2, m=1, params=theta_222logistic_1_2), c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222)))
  expect_equal(pick_regime(p=2, M=2, d=2, m=2, params=theta_222logistic_1_2), c(phi20_222, vec(A21_222), vec(A22_222), vech(Omega2_222)))
  expect_equal(pick_regime(p=1, M=2, d=3, m=1, params=theta_123logistic_1_1), c(phi10_123, vec(A11_123), vech(Omega1_123)))
  expect_equal(pick_regime(p=1, M=2, d=3, m=2, params=theta_123logistic_1_1), c(phi20_123, vec(A21_123), vech(Omega2_123)))

  expect_equal(pick_regime(p=1, M=2, d=2, m=1, params=theta_122exp_1_1), c(phi10_122, vec(A11_122), vech(Omega1_122)))
  expect_equal(pick_regime(p=1, M=2, d=2, m=2, params=theta_122exp_1_1), c(phi20_122, vec(A21_122), vech(Omega2_122)))
  expect_equal(pick_regime(p=2, M=2, d=2, m=1, params=theta_222exp_1_2), c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222)))
  expect_equal(pick_regime(p=2, M=2, d=2, m=2, params=theta_222exp_1_2), c(phi20_222, vec(A21_222), vec(A22_222), vech(Omega2_222)))
  expect_equal(pick_regime(p=1, M=2, d=3, m=1, params=theta_123exp_1_1), c(phi10_123, vec(A11_123), vech(Omega1_123)))
  expect_equal(pick_regime(p=1, M=2, d=3, m=2, params=theta_123exp_1_1), c(phi20_123, vec(A21_123), vech(Omega2_123)))

  expect_equal(pick_regime(p=1, M=2, d=2, m=1, params=theta_122thres_1_1), c(phi10_122, vec(A11_122), vech(Omega1_122)))
  expect_equal(pick_regime(p=1, M=2, d=2, m=2, params=theta_122thres_1_1), c(phi20_122, vec(A21_122), vech(Omega2_122)))
  expect_equal(pick_regime(p=2, M=2, d=2, m=1, params=theta_222thres_2_2), c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222)))
  expect_equal(pick_regime(p=2, M=2, d=2, m=2, params=theta_222thres_2_2), c(phi20_222, vec(A21_222), vec(A22_222), vech(Omega2_222)))
  expect_equal(pick_regime(p=1, M=3, d=2, m=1, params=theta_132thres_1_1), c(phi10_132, vec(A11_132), vech(Omega1_132)))
  expect_equal(pick_regime(p=1, M=3, d=2, m=2, params=theta_132thres_1_1), c(phi20_132, vec(A21_132), vech(Omega2_132)))
  expect_equal(pick_regime(p=1, M=3, d=2, m=3, params=theta_132thres_1_1), c(phi30_132, vec(A31_132), vech(Omega3_132)))
  expect_equal(pick_regime(p=1, M=2, d=3, m=1, params=theta_123thres_2_1), c(phi10_123, vec(A11_123), vech(Omega1_123)))
  expect_equal(pick_regime(p=1, M=2, d=3, m=2, params=theta_123thres_2_1), c(phi20_123, vec(A21_123), vech(Omega2_123)))

  expect_equal(pick_regime(p=2, M=2, d=2, m=1, params=theta_222exo), c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222)))
  expect_equal(pick_regime(p=2, M=2, d=2, m=2, params=theta_222exo), c(phi20_222, vec(A21_222), vec(A22_222), vech(Omega2_222)))
  expect_equal(pick_regime(p=1, M=3, d=2, m=1, params=theta_132exo), c(phi10_132, vec(A11_132), vech(Omega1_132)))
  expect_equal(pick_regime(p=1, M=3, d=2, m=2, params=theta_132exo), c(phi20_132, vec(A21_132), vech(Omega2_132)))
  expect_equal(pick_regime(p=1, M=3, d=2, m=3, params=theta_132exo), c(phi30_132, vec(A31_132), vech(Omega3_132)))
  expect_equal(pick_regime(p=1, M=2, d=3, m=1, params=theta_123exo), c(phi10_123, vec(A11_123), vech(Omega1_123)))
  expect_equal(pick_regime(p=1, M=2, d=3, m=2, params=theta_123exo), c(phi20_123, vec(A21_123), vech(Omega2_123)))

  # Student
  expect_equal(pick_regime(p=2, M=3, d=2, m=1, params=theta_232threst_1_1), c(phi10_232, vec(A11_232), vec(A12_232), vech(Omega1_232)))
  expect_equal(pick_regime(p=2, M=3, d=2, m=2, params=theta_232threst_1_1), c(phi20_232, vec(A21_232), vec(A22_232), vech(Omega2_232)))
  expect_equal(pick_regime(p=2, M=3, d=2, m=3, params=theta_232threst_1_1), c(phi30_232, vec(A31_232), vec(A32_232), vech(Omega3_232)))
  expect_equal(pick_regime(p=1, M=2, d=2, m=1, params=theta_122expt_1_1), c(phi10_122, vec(A11_122), vech(Omega1_122)))
  expect_equal(pick_regime(p=1, M=2, d=2, m=2, params=theta_122expt_1_1), c(phi20_122, vec(A21_122), vech(Omega2_122)))
  expect_equal(pick_regime(p=2, M=2, d=2, m=1, params=theta_222logistict_2_1), c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222)))
  expect_equal(pick_regime(p=2, M=2, d=2, m=2, params=theta_222logistict_2_1), c(phi20_222, vec(A21_222), vec(A22_222), vech(Omega2_222)))
  expect_equal(pick_regime(p=1, M=2, d=3, m=1, params=theta_123log_1_1), c(phi10_123, vec(A11_123), vech(Omega1_123)))
  expect_equal(pick_regime(p=1, M=2, d=3, m=2, params=theta_123log_1_1), c(phi20_123, vec(A21_123), vech(Omega2_123)))

  # ind_Student
  expect_equal(pick_regime(p=1, M=1, d=2, m=1, params=theta_112it, cond_dist="ind_Student"), c(phi10_112, vec(A11_112), vec(B1_112)))
  expect_equal(pick_regime(p=2, M=3, d=2, m=1, params=theta_232exoit, cond_dist="ind_Student"),
               c(phi10_232, vec(A11_232), vec(A12_232), vec(B1_232)))
  expect_equal(pick_regime(p=2, M=3, d=2, m=2, params=theta_232exoit, cond_dist="ind_Student"),
               c(phi20_232, vec(A21_232), vec(A22_232), vec(B2_232)))
  expect_equal(pick_regime(p=2, M=3, d=2, m=3, params=theta_232exoit, cond_dist="ind_Student"),
               c(phi30_232, vec(A31_232), vec(A32_232), vec(B3_232)))
  expect_equal(pick_regime(p=1, M=2, d=2, m=1, params=theta_122logit_1_1, cond_dist="ind_Student"), c(phi10_122, vec(A11_122), vec(B1_122)))
  expect_equal(pick_regime(p=1, M=2, d=2, m=2, params=theta_122logit_1_1, cond_dist="ind_Student"), c(phi20_122, vec(A21_122), vec(B2_122)))
  expect_equal(pick_regime(p=2, M=2, d=2, m=1, params=theta_222logistit_2_1, cond_dist="ind_Student"),
               c(phi10_222, vec(A11_222), vec(A12_222), vec(B1_222)))
  expect_equal(pick_regime(p=2, M=2, d=2, m=2, params=theta_222logistit_2_1, cond_dist="ind_Student"),
               c(phi20_222, vec(A21_222), vec(A22_222), vec(B2_222)))
  expect_equal(pick_regime(p=1, M=2, d=3, m=1, params=theta_123expit_1_1, cond_dist="ind_Student"), c(phi10_123, vec(A11_123), vec(B1_123)))
  expect_equal(pick_regime(p=1, M=2, d=3, m=2, params=theta_123expit_1_1, cond_dist="ind_Student"), c(phi20_123, vec(A21_123), vec(B2_123)))
})


test_that("pick_W works correctly", {
  expect_equal(pick_W(p=1, M=2, d=2, params=theta_122relgsh, identification="heteroskedasticity"), W_122)
  expect_equal(pick_W(p=1, M=3, d=2, params=theta_132relgsh, identification="heteroskedasticity"), W_132)
  expect_equal(pick_W(p=2, M=2, d=2, params=theta_222logistictsh_2_1, identification="heteroskedasticity"), W_222)
  expect_equal(pick_W(p=1, M=2, d=2, params=theta_122logsh_12_1, identification="heteroskedasticity"), W_122)
  expect_equal(pick_W(p=1, M=2, d=3, params=theta_123expsh_1_1, identification="heteroskedasticity"), W_123)
  expect_equal(pick_W(p=2, M=3, d=2, params=theta_232threstsh_1_1, identification="heteroskedasticity"), W_232)
})


test_that("pick_lambdas works correctly", {
  expect_equal(pick_lambdas(p=1, M=2, d=2, params=theta_122relgsh, identification="heteroskedasticity"), lambdas_122)
  expect_equal(pick_lambdas(p=1, M=3, d=2, params=theta_132relgsh, identification="heteroskedasticity"),
               c(lambdas2_132, lambdas3_132))
  expect_equal(pick_lambdas(p=2, M=2, d=2, params=theta_222logistictsh_2_1, identification="heteroskedasticity"), lambdas_222)
  expect_equal(pick_lambdas(p=1, M=2, d=2, params=theta_122logsh_12_1, identification="heteroskedasticity"), lambdas_122)
  expect_equal(pick_lambdas(p=1, M=2, d=3, params=theta_123expsh_1_1, identification="heteroskedasticity"), lambdas_123)
  expect_equal(pick_lambdas(p=2, M=3, d=2, params=theta_232threstsh_1_1, identification="heteroskedasticity"),
               c(lambdas2_232, lambdas3_232))
})

test_that("pick_distpars works correctly", {
  expect_equal(pick_distpars(d=2, params=theta_122relg, cond_dist="Gaussian"), numeric(0))
  expect_equal(pick_distpars(d=2, params=theta_222logistict_2_1, cond_dist="Student"), df_222_2_1)
  expect_equal(pick_distpars(d=2, params=theta_112it, cond_dist="ind_Student"), dfs_112)
  expect_equal(pick_distpars(d=2, params=theta_132thresit_1_1, cond_dist="ind_Student"), dfs_132_1_1)
  expect_equal(pick_distpars(d=2, params=theta_132thresit_1_1, cond_dist="ind_Student"), dfs_132_1_1)
  expect_equal(pick_distpars(d=2, params=theta_222logistit_2_1, cond_dist="ind_Student"), dfs_222_2_1)
})
