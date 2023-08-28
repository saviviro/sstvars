context("eigens")
library(sstvars)

## A(M)(p)_(p)(M)(d)

# p=1, M=1, d=2
phi10_112 <- c(0.65, 0.7)
A11_112 <- matrix(c(0.29, 0.02, -0.14, 0.9), nrow=2, byrow=FALSE)
Omega1_112 <- matrix(c(0.60, 0.01, 0.01, 0.07), nrow=2, byrow=FALSE)

theta_112relg <- c(phi10_112, vec(A11_112), vech(Omega1_112))

stvar112 <- STVAR(p=1, M=1, d=2, params=theta_112relg, weight_function="relative_dens")

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

stvar312 <- STVAR(p=3, M=1, d=2, params=theta_312relg, weight_function="relative_dens")

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

stvar222 <- STVAR(p=2, M=2, d=2, params=theta_222relg, weight_function="relative_dens")

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

stvar132 <- STVAR(p=1, M=3, d=2, params=theta_132relg, weight_function="relative_dens")


# p=1, M=1, d=3
phi10_113 <- c(1, 2, 3)
A11_113 <- matrix(c(0.1, 0.02, 0.12, 0.3, 0.21, 0.11, 0.05, 0.03, 0.09), nrow=3)
Omega1_113 <- matrix(c(c(1, 0.2, 0.3, 0.2, 2, 0.4, 0.3, 0.4, 3)), nrow=3)

theta_113relg <- c(phi10_113, vec(A11_113), vech(Omega1_113))

stvar113 <- STVAR(p=1, M=1, d=3, params=theta_113relg, weight_function="relative_dens")

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

stvar123 <- STVAR(p=1, M=2, d=3, params=theta_123relg, weight_function="relative_dens")


## weight_function = "mlogit"

# p=2, M=2, d=2, weightfun_pars=list(vars=1:2, lags=2)
gamma1_222_12_2 <- c(0.1, 0.2, 0.11, 0.22, 0.33)
theta_222log_12_2 <- c(theta_222relg[-length(theta_222relg)], gamma1_222_12_2)

stvar222mlogit <- STVAR(p=2, M=2, d=2, params=theta_222log_12_2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2))


## weight_function == "threshold"

# p=2, M=2, d=2, weight_function="threshold", weighfun_pars=c(2, 1)
theta_222thres_2_1 <- c(theta_222relg[-length(theta_222relg)], 1)
mod222thres_2_1 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222thres_2_1, weight_function="threshold",
                         weightfun_pars=c(2, 1))


## Constrained models
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}

# p=2, M=2, d=2, C_222
C_222 <- rbind_diags(p=2, M=2, d=2)
theta_222relgc <- c(0.36, 0.12, 0.48, 0.07, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, 0.21, 0.01, 0.03, 1.1, 0.01, 0.11, 0.37)
stvar222relgc <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222relgc, weight_function="relative_dens",
                       AR_constraints=C_222)

# p=2, M=2, d=2, mean_constraints=list(1:2), AR_constraints=C_222, parametrization="mean"
theta_222relgcm <- c(0.7209658, 0.810858, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, 0.21, 0.01,
                     0.03, 1.1, 0.01, 0.11, 0.37)
stvar222relgcm <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222relgcm, weight_function="relative_dens",
                        AR_constraints=C_222, mean_constraints=list(1:2), parametrization="mean")

## mlogit

# p=2, M=2, d=2, weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222, parametrization="mean"
theta_222logcm_12_2 <- c(phi10_222, vec(A11_222), vec(A12_222), vech(Omega1_222), vech(Omega2_222), gamma1_222_12_2)
theta_222logcm_12_2_expanded <- c(phi10_222, phi10_222, vec(A11_222), vec(A12_222), vec(A11_222), vec(A12_222),
                                  vech(Omega1_222), vech(Omega2_222), gamma1_222_12_2)

stvar222mlogitcm <- STVAR(p=2, M=2, d=2, params=theta_222logcm_12_2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2),
                         AR_constraints=C_222, mean_constraints=list(1:2), parametrization="mean")


## Models with weight_constraints

# p=1, M=3, d=2, weight_function="relative_dens", weight_constraints=list(R=matrix(c(0.9, 0.5), nrow=2), r=c(0.13, 0.13))
xi_132relgw <- 0.4
theta_132relgw <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                    vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), xi_132relgw)
theta_132relgw_expanded <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                             vech(Omega1_132), vech(Omega2_132), vech(Omega3_132), matrix(c(0.9, 0.5), nrow=2)%*%xi_132relgw + 0.13)

mod132relgw <- STVAR(data=gdpdef, p=1, M=3, d=2, params=theta_132relgw,  weight_function="relative_dens",
                     weight_constraints=list(R=matrix(c(0.9, 0.5), nrow=2), r=c(0.13, 0.13)))


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
mod222logcmw_12_2 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222logcmw_12_2, weight_function="mlogit", parametrization="mean",
                           weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222,
                           weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5), r=c(0, 0.11, 0.12, 0.13, 0)))


# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean"
xi_222expcmw_2_1 <- c(0.33)
theta_222expcmw_2_1 <- c(theta_222relgcm[-length(theta_222relgcm)], xi_222expcmw_2_1)
mod222expcmw_2_1 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222expcmw_2_1, weight_function="exponential",
                          weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222,
                          weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean")

## Student

# p=2, M=2, d=2, weight_function="threshold", weighfun_pars=c(2, 1), cond_dist="Student"
theta_222threst_2_1 <- c(theta_222thres_2_1, 13)
mod222threst_2_1 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222threst_2_1, weight_function="threshold",
                          weightfun_pars=c(2, 1), cond_dist="Student")

# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2), cond_dist="Student", mean_constraints=list(1:2),
# AR_constraints=C_222, parametrization="mean"
theta_222logcmt_12_2 <- c(theta_222logcm_12_2, 2.13)
mod222logcmt_12_2 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222logcmt_12_2, weight_function="mlogit", parametrization="mean",
                           weightfun_pars=list(vars=1:2, lags=2), cond_dist="Student", mean_constraints=list(1:2), AR_constraints=C_222)


# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean"
theta_222expcmwt_2_1 <- c(theta_222expcmw_2_1, 4)
mod222expcmwt_2_1 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222expcmwt_2_1, weight_function="exponential", parametrization="mean",
                           weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2), AR_constraints=C_222,
                           weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)))

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean"
theta_222logisticcmwt_2_1 <- theta_222expcmwt_2_1
mod222logisticcmwt_2_1 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222logisticcmwt_2_1, weight_function="logistic",
                                weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2), AR_constraints=C_222,
                                weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean")


## Structural models

# p=2, M=2, d=2, weight_function="threshold", weighfun_pars=c(2, 1), identification="recursive"
mod222thressr_2_1 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222thres_2_1, weight_function="threshold",
                           weightfun_pars=c(2, 1), identification="recursive")

# p=1, M=2, d=2, weight_function="relative_dens", identification="heteroskedasticity"
all_phi_122 <- c(0.734054, 0.225598, 0.705744, 0.187897)
all_A_122 <- c(0.259626, -0.000863, -0.3124, 0.505251, 0.298483, 0.030096, -0.176925, 0.838898)
W_122 <- matrix(c(-0.03, 0.24, -0.76, -0.02), nrow=2, ncol=2, byrow=FALSE)
lambdas_122 <- c(3.36, 0.86)
alpha1_122 <- 0.6
theta_122relgsh <- c(all_phi_122, all_A_122, vec(W_122), lambdas_122, alpha1_122)
mod122relgsh <- STVAR(data=gdpdef, p=1, M=2, d=2, params=theta_122relgsh, weight_function="relative_dens", identification="heteroskedasticity")

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student",
# identification="heteroskedasticity"
all_phi_222 <- c(0.356914, 0.107436, 0.356386, 0.08633)
all_A_222 <- c(0.13996, 0.035172, -0.164575, 0.386816, 0.451675, 0.013086, 0.227882, 0.336084, 0.239257, 0.024173,
               -0.021209, 0.707502, 0.063322, 0.027287, 0.009182, 0.197066)
W_222 <- W_122; lambdas_222 <- lambdas_122
c_and_gamma_222_2_1 <- c(0.1, 0.2)
df_222_2_1 <- 7
theta_222logistictsh_2_1 <- c(all_phi_222, all_A_222, vec(W_222), lambdas_222, c_and_gamma_222_2_1, df_222_2_1)
mod222logistictsh_2_1 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222logistictsh_2_1, weight_function="logistic", weightfun_pars=c(2, 1),
                               cond_dist="Student", identification="heteroskedasticity")


## Structural models imposing constraints

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean"
theta_222expcmwtsh_2_1 <- c(0.7209658, 0.810858, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, # mu + A
                            -0.03, 0.24, -0.76, -0.02, 3.36, 0.86, # W + lambdas
                            0.33, 4) # xi + nu
mod222expcmwtsh_2_1 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222expcmwtsh_2_1, weight_function="exponential",
                             weightfun_pars=c(2, 1), identification="heteroskedasticity", cond_dist="Student", mean_constraints=list(1:2),
                             AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean")



test_that("get_boldA_eigens work correctly", {
  expect_equal(get_boldA_eigens(stvar112), as.matrix(c(0.8953748, 0.2946252)), tol=1e-3)
  expect_equal(get_boldA_eigens(stvar312), as.matrix(c(0.6367684, 0.6367684, 0.5731850, 0.5466223, 0.5256310, 0.5256310)), tol=1e-3)
  expect_equal(get_boldA_eigens(stvar222)[,1], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)
  expect_equal(get_boldA_eigens(stvar222)[,2], c(0.9198650, 0.4317946, 0.2541513, 0.1575083), tol=1e-3)
  expect_equal(get_boldA_eigens(stvar132)[,1], c(0.7186796, 0.3413204), tol=1e-3)
  expect_equal(get_boldA_eigens(stvar132)[,2], c(0.7186796, 0.3413204), tol=1e-3)
  expect_equal(get_boldA_eigens(stvar132)[,3], c(0.53722813, 0.03722813), tol=1e-3)
  expect_equal(get_boldA_eigens(stvar113), as.matrix(c(0.29792217, 0.05646904, 0.05646904)), tol=1e-3)
  expect_equal(get_boldA_eigens(stvar123)[,1], c(0.29792217, 0.05646904, 0.05646904), tol=1e-3)
  expect_equal(get_boldA_eigens(stvar123)[,2], c(0.28054283, 0.10466989, 0.03521272), tol=1e-3)

  expect_equal(get_boldA_eigens(stvar222mlogit)[,1], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)
  expect_equal(get_boldA_eigens(stvar222mlogit)[,2], c(0.9198650, 0.4317946, 0.2541513, 0.1575083), tol=1e-3)

  expect_equal(get_boldA_eigens(stvar222relgc)[,1], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)
  expect_equal(get_boldA_eigens(stvar222relgc)[,2], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)
  expect_equal(get_boldA_eigens(stvar222relgcm)[,1], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)
  expect_equal(get_boldA_eigens(stvar222relgcm)[,2], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)

  expect_equal(get_boldA_eigens(stvar222mlogitcm)[,1], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)
  expect_equal(get_boldA_eigens(stvar222mlogitcm)[,2], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)

  expect_equal(get_boldA_eigens(mod132relgw)[1,], c(0.7186796, 0.7186796, 0.5372281), tol=1e-3)
  expect_equal(get_boldA_eigens(mod132relgw)[2,], c(0.34132038, 0.34132038, 0.03722813), tol=1e-3)

  expect_equal(get_boldA_eigens(mod222logcmw_12_2)[,1], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)
  expect_equal(get_boldA_eigens(mod222logcmw_12_2)[,2], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)

  expect_equal(get_boldA_eigens(mod222thres_2_1)[,1], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)
  expect_equal(get_boldA_eigens(mod222thres_2_1)[,2], c(0.9198650, 0.4317946, 0.2541513, 0.1575083), tol=1e-3)

  expect_equal(get_boldA_eigens(mod222expcmw_2_1)[,1], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)
  expect_equal(get_boldA_eigens(mod222expcmw_2_1)[,2], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)

  # Student
  expect_equal(get_boldA_eigens(mod222threst_2_1)[,1], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)
  expect_equal(get_boldA_eigens(mod222threst_2_1)[,2], c(0.9198650, 0.4317946, 0.2541513, 0.1575083), tol=1e-3)
  expect_equal(get_boldA_eigens(mod222logcmt_12_2)[,1], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)
  expect_equal(get_boldA_eigens(mod222logcmt_12_2)[,2], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)
  expect_equal(get_boldA_eigens(mod222expcmwt_2_1)[,1], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)
  expect_equal(get_boldA_eigens(mod222expcmwt_2_1)[,2], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)
  expect_equal(get_boldA_eigens(mod222logisticcmwt_2_1)[,1], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)
  expect_equal(get_boldA_eigens(mod222logisticcmwt_2_1)[,2], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)

  # Structural models
  expect_equal(get_boldA_eigens(mod222thressr_2_1)[,1], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)
  expect_equal(get_boldA_eigens(mod222thressr_2_1)[,2], c(0.9198650, 0.4317946, 0.2541513, 0.1575083), tol=1e-3)
  expect_equal(get_boldA_eigens(mod122relgsh)[,1], c(0.5063438, 0.2585332), tol=1e-3)
  expect_equal(get_boldA_eigens(mod122relgsh)[,2], c(0.8288584, 0.3085226), tol=1e-3)
  expect_equal(get_boldA_eigens(mod222logistictsh_2_1)[,1], c(0.8300715, 0.7143956, 0.5980901, 0.4196011), tol=1e-3)
  expect_equal(get_boldA_eigens(mod222logistictsh_2_1)[,2], c(0.9205818, 0.3981026, 0.2208400, 0.1510854), tol=1e-3)
  expect_equal(get_boldA_eigens(mod222expcmwtsh_2_1)[,1], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)
  expect_equal(get_boldA_eigens(mod222expcmwtsh_2_1)[,2], c(0.7667397, 0.7667397, 0.5076821, 0.4147943), tol=1e-3)
})


test_that("get_omega_eigens work correctly", {
  expect_equal(get_omega_eigens(stvar112), as.matrix(c(0.60018861, 0.06981139)), tol=1e-3)
  expect_equal(get_omega_eigens(stvar312), as.matrix(c(0.58019224, 0.05980776)), tol=1e-3)
  expect_equal(get_omega_eigens(stvar222)[,1], c(0.21055385, 0.02944615), tol=1e-3)
  expect_equal(get_omega_eigens(stvar222)[,2], c(1.100101, 0.109899), tol=1e-3)
  expect_equal(get_omega_eigens(stvar132)[,1], c(0.58019224, 0.05980776), tol=1e-3)
  expect_equal(get_omega_eigens(stvar132)[,2], c(0.500333, 0.199667), tol=1e-3)
  expect_equal(get_omega_eigens(stvar132)[,3], c(1.5, 0.5), tol=1e-3)
  expect_equal(get_omega_eigens(stvar113), as.matrix(c(3.1960208, 1.8679839, 0.9359953)), tol=1e-3)
  expect_equal(get_omega_eigens(stvar123)[,1], c(3.1960208, 1.8679839, 0.9359953), tol=1e-3)
  expect_equal(get_omega_eigens(stvar123)[,2], c(3.451839, 2.143882, 1.004279), tol=1e-3)

  expect_equal(get_omega_eigens(stvar222mlogit)[,1], c(0.21055385, 0.02944615), tol=1e-3)
  expect_equal(get_omega_eigens(stvar222mlogit)[,2], c(1.100101, 0.109899), tol=1e-3)

  expect_equal(get_omega_eigens(stvar222relgc)[,1], c(0.21055385, 0.02944615), tol=1e-3)
  expect_equal(get_omega_eigens(stvar222relgc)[,2], c(1.100101, 0.109899), tol=1e-3)
  expect_equal(get_omega_eigens(stvar222relgcm)[,1], c(0.21055385, 0.02944615), tol=1e-3)
  expect_equal(get_omega_eigens(stvar222relgcm)[,2], c(1.100101, 0.109899), tol=1e-3)

  expect_equal(get_omega_eigens(stvar222mlogitcm)[,1], c(0.21055385, 0.02944615), tol=1e-3)
  expect_equal(get_omega_eigens(stvar222mlogitcm)[,2], c(1.100101, 0.109899), tol=1e-3)

  expect_equal(get_omega_eigens(mod132relgw)[1,], c(0.5801922, 0.5003330, 1.5000000), tol=1e-3)
  expect_equal(get_omega_eigens(mod132relgw)[2,], c(0.05980776, 0.19966704, 0.50000000), tol=1e-3)

  expect_equal(get_omega_eigens(mod222logcmw_12_2)[,1], c(0.21055385, 0.02944615), tol=1e-3)
  expect_equal(get_omega_eigens(mod222logcmw_12_2)[,2], c(1.100101, 0.109899), tol=1e-3)

  expect_equal(get_omega_eigens(mod222thres_2_1)[,1], c(0.21055385, 0.02944615), tol=1e-3)
  expect_equal(get_omega_eigens(mod222thres_2_1)[,2], c(1.100101, 0.109899), tol=1e-3)

  expect_equal(get_omega_eigens(mod222expcmw_2_1)[,1], c(0.21055385, 0.02944615), tol=1e-3)
  expect_equal(get_omega_eigens(mod222expcmw_2_1)[,2], c(1.100101, 0.109899), tol=1e-3)

  # Student
  expect_equal(get_omega_eigens(mod222threst_2_1)[,1], c(0.21055385, 0.02944615), tol=1e-3)
  expect_equal(get_omega_eigens(mod222threst_2_1)[,2], c(1.100101, 0.109899), tol=1e-3)
  expect_equal(get_omega_eigens(mod222logcmt_12_2)[,1], c(0.21055385, 0.02944615), tol=1e-3)
  expect_equal(get_omega_eigens(mod222logcmt_12_2)[,2], c(1.100101, 0.109899), tol=1e-3)
  expect_equal(get_omega_eigens(mod222expcmwt_2_1)[,1], c(0.21055385, 0.02944615), tol=1e-3)
  expect_equal(get_omega_eigens(mod222expcmwt_2_1)[,2], c(1.100101, 0.109899), tol=1e-3)
  expect_equal(get_omega_eigens(mod222logisticcmwt_2_1)[,1], c(0.21055385, 0.02944615), tol=1e-3)
  expect_equal(get_omega_eigens(mod222logisticcmwt_2_1)[,2], c(1.100101, 0.109899), tol=1e-3)

  # Structural
  expect_equal(get_omega_eigens(mod222thressr_2_1)[,1], c(0.21055385, 0.02944615), tol=1e-3)
  expect_equal(get_omega_eigens(mod222thressr_2_1)[,2], c(1.100101, 0.109899), tol=1e-3)
  expect_equal(get_omega_eigens(mod122relgsh)[,1], c(0.57862293, 0.05787707), tol=1e-3)
  expect_equal(get_omega_eigens(mod122relgsh)[,2], c(0.5001637, 0.1934763), tol=1e-3)
  expect_equal(get_omega_eigens(mod222logistictsh_2_1)[,1], c(0.57862293, 0.05787707), tol=1e-3)
  expect_equal(get_omega_eigens(mod222logistictsh_2_1)[,2], c(0.5001637, 0.1934763), tol=1e-3)
  expect_equal(get_omega_eigens(mod222expcmwtsh_2_1)[,1], c(0.57862293, 0.05787707), tol=1e-3)
  expect_equal(get_omega_eigens(mod222expcmwtsh_2_1)[,2], c(0.5001637, 0.1934763), tol=1e-3)
})
