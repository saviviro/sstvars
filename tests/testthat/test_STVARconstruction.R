context("STVARconstruction")
library(sstvars)

# p=1, M=1, d=2
theta_112relg <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103, 0.601786, -0.002945, 0.067224)
mod112relg <- STVAR(data=gdpdef, p=1, M=1, params=theta_112relg, weight_function="relative_dens")

# p=1, M=2, d=2
theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863, -0.3124, 0.505251, 0.298483, 0.030096,
                   -0.176925, 0.838898, 0.310863, 0.007512, 0.018244, 0.949533, -0.016941, 0.121403, 0.573269)
mod122relg <- STVAR(p=1, M=2, d=2, params=theta_122relg, weight_function="relative_dens")

# p=2, M=2, d=2
theta_222relg <- c(0.356914, 0.107436, 0.356386, 0.08633, 0.13996, 0.035172, -0.164575, 0.386816, 0.451675, 0.013086,
                   0.227882, 0.336084, 0.239257, 0.024173, -0.021209, 0.707502, 0.063322, 0.027287, 0.009182, 0.197066,
                   0.205831, 0.005157, 0.025877, 1.092094, -0.009327, 0.116449, 0.592446)
mod222relg <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222relg, weight_function="relative_dens")

# p=1, M=2, d=3, usamone
theta_123relg <- c(0.10741, 0.13813, -0.12092, 3.48957, 0.60615, 0.45646, 0.87227, -0.01595, 0.14124,
                   -0.08611, 0.61865, 0.34311, -0.02047, 0.025, 0.97548, 0.74976, 0.02187, 0.29213,
                   -1.55165, 0.58245, -0.00696, -0.07261, 0.02021, 0.96883, 0.66149, 0.02279, 0.09207,
                   0.05544, 0.00212, 0.12708, 0.78618, 0.00922, 0.42627, 0.23765, 0.25386, 3.40834, 0.77357)
mod123relg <- STVAR(data=usamone, p=1, M=2, params=theta_123relg, weight_function="relative_dens")

# p=3, M=2, d=3, usamone
theta_323relg <- c(0.98249, 0.66144, -1.17552, 0.50289, 0.17399, -0.01771, 0.96105, -0.11406, 0.41223,
                   -0.31217, 0.49067, 0.3958, 0.04185, 0.08454, 1.0977, -0.03208, 0.06398, -0.12298,
                   0.13382, 0.20166, 0.87613, -0.34591, -0.06254, -0.47386, -0.09049, 0.03109, 0.0347,
                   -0.16531, 0.0427, -0.31646, 0.25299, -0.04865, 0.33893, 0.69963, -0.02912, 0.03398,
                   -0.24344, 0.20815, 0.22566, 0.20582, 0.14774, 1.69008, 0.04375, -0.01018, -0.00947,
                   -0.19371, 0.26341, 0.22082, -0.08841, -0.18303, -0.86488, -0.06031, 0.00634, 0.00181,
                   -0.5559, 0.10249, -0.25146, -0.11875, 0.05153, 0.15267, 0.58151, -0.01903, 0.12236, 0.09327,
                   0.10245, 1.81845, 0.72719, 0.03235, 0.09857, 0.04826, 0.00908, 0.09761, 0.72127)
mod323relg <- STVAR(data=usamone, p=3, M=2, params=theta_323relg, weight_function="relative_dens")


## weight_function == "threshold"

# p=2, M=2, d=2, weight_function="threshold", weighfun_pars=c(2, 1)
theta_222thres_2_1 <- c(theta_222relg[-length(theta_222relg)], 1)
mod222thres_2_1 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222thres_2_1, weight_function="threshold",
                         weightfun_pars=c(2, 1))

## weight_function == "exogenous"

# p=2, M=2, d=2, weight_function="exogenous", weighfun_pars=weightfun_pars222, cond_dist="Student"
set.seed(2); tw1 <- runif(nrow(gdpdef) - 2)
weightfun_pars222 <- cbind(tw1, 1-tw1)
theta_222exo <- c(theta_222relg[-length(theta_222relg)], 7)
mod222exo <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222exo, weight_function="exogenous",
                   weightfun_pars=weightfun_pars222, cond_dist="Student")

## Constrained models
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}

# p=2, M=2, d=2, mean_constraints=list(1:2), AR_constraints=C_222, parametrization="mean"
C_222 <- rbind_diags(p=2, M=2, d=2)
theta_222relgcm <- c(0.7209658, 0.810858, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, 0.21, 0.01,
                     0.03, 1.1, 0.01, 0.11, 0.37)
mod222relgcm <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222relgcm, weight_function="relative_dens",
                      AR_constraints=C_222, mean_constraints=list(1:2), parametrization="mean")


# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2),
# AR_constraints=C_222, parametrization="mean"
gamma1_222_12_2 <- c(0.1, 0.2, 0.11, 0.22, 0.33)
theta_222logcm_12_2 <- c(theta_222relgcm[-length(theta_222relgcm)], gamma1_222_12_2)
mod222logcm_12_2 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222logcm_12_2, weight_function="mlogit", parametrization="mean",
                          weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222)


# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5), r=c(0, 0.11, 0.12, 0.13, 0)), parametrization="mean"
xi_222logcmw_12_2 <- c(0.002, 1.33)
theta_222logcmw_12_2 <-  c(theta_222relgcm[-length(theta_222relgcm)], xi_222logcmw_12_2)
mod222logcmw_12_2 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222logcmw_12_2, weight_function="mlogit", parametrization="mean",
                          weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222,
                          weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5), r=c(0, 0.11, 0.12, 0.13, 0)))

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))
xi_222logisticcmw_2_1 <- c(0.33)
theta_222logisticcmw_2_1 <- c(theta_222relgcm[-length(theta_222relgcm)], xi_222logisticcmw_2_1)
mod222logisticcmw_2_1 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222logisticcmw_2_1, weight_function="logistic",
                               weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222, parametrization="mean",
                               weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)))

tmppars <- c(theta_222relg[-length(theta_222relg)], c(0.75, 100))
tmpmod <- STVAR(data=gdpdef, p=2, M=2, d=2, params=tmppars, weight_function="exponential", weightfun_pars=c(2, 1))
plot(tmpmod)

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0))
xi_222expcmw_2_1 <- c(0.33)
theta_222expcmw_2_1 <- c(theta_222relgcm[-length(theta_222relgcm)], xi_222expcmw_2_1)
mod222expcmw_2_1 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222expcmw_2_1, weight_function="exponential",
                          weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222, parametrization="mean",
                          weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)))


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

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2),
# AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean"
theta_222logisticcmwt_2_1 <- c(theta_222logisticcmw_2_1, 30)
mod222logisticcmwt_2_1 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222logisticcmwt_2_1, weight_function="logistic",
                               weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2), AR_constraints=C_222,
                               weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean")


# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2),
# AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean"
theta_222expcmwt_2_1 <- c(theta_222expcmw_2_1, 4)
mod222expcmwt_2_1 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222expcmwt_2_1, weight_function="exponential",
                           weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2), AR_constraints=C_222,
                           weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean")

## ind_Student

# p=1, M=2, d=3, weight_function="exogenous", weighfun_pars=weightfun_pars123, cond_dist="ind_Student"
set.seed(3); tw1 <- runif(nrow(usamone) - 1)
weightfun_pars123 <- cbind(tw1, 1-tw1)
set.seed(4); Bmatpars123 <- round(rnorm(18), 3)
theta_123exoit <- c(0.10741, 0.13813, -0.12092, 3.48957, 0.60615, 0.45646, 0.87227, -0.01595, 0.14124,
                  -0.08611, 0.61865, 0.34311, -0.02047, 0.025, 0.97548, 0.74976, 0.02187, 0.29213,
                  -1.55165, 0.58245, -0.00696, -0.07261, 0.02021, 0.96883, Bmatpars123, 7, 3, 13)
mod123exoit <- STVAR(data=usamone, p=1, M=2, d=3, params=theta_123exoit, weight_function="exogenous",
                   weightfun_pars=weightfun_pars123, cond_dist="ind_Student")

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="ind_Student", mean_constraints=list(1:2),
# AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean"
set.seed(5); Bmatpars222 <- round(rnorm(8), 3)
theta_222logistit <- c(0.7209658, 0.810858, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, Bmatpars222, 0.4, 7, 3)
mod222logistit <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222logistit, weight_function="logistic", weightfun_pars=c(2, 1),
                        cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_222,
                        weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean")

## ind_skewed_t

# p=1, M=2, d=3, weight_function="exogenous", weighfun_pars=weightfun_pars123, cond_dist="ind_skewed_t"
theta_123exoikt <- c(0.10741, 0.13813, -0.12092, 3.48957, 0.60615, 0.45646, 0.87227, -0.01595, 0.14124,
                    -0.08611, 0.61865, 0.34311, -0.02047, 0.025, 0.97548, 0.74976, 0.02187, 0.29213,
                    -1.55165, 0.58245, -0.00696, -0.07261, 0.02021, 0.96883, Bmatpars123, 7, 3, 13, -0.1, 0.2, 0.3)
mod123exoikt <- STVAR(data=usamone, p=1, M=2, d=3, params=theta_123exoikt, weight_function="exogenous",
                      weightfun_pars=weightfun_pars123, cond_dist="ind_skewed_t")

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="ind_skewed_t", mean_constraints=list(1:2),
# AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean"
theta_222logistikt <- c(0.7209658, 0.810858, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, Bmatpars222, 0.4, 7, 3, 0.1, 0)
mod222logistikt <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222logistikt, weight_function="logistic", weightfun_pars=c(2, 1),
                        cond_dist="ind_skewed_t", mean_constraints=list(1:2), AR_constraints=C_222,
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
mod122relgsh <- STVAR(data=gdpdef, p=1, M=2, d=2, params=theta_122relgsh, weight_function="relative_dens",
                      identification="heteroskedasticity")

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student",
# identification="heteroskedasticity"
all_phi_222 <- c(0.356914, 0.107436, 0.356386, 0.08633)
all_A_222 <- c(0.13996, 0.035172, -0.164575, 0.386816, 0.451675, 0.013086, 0.227882, 0.336084, 0.239257, 0.024173,
               -0.021209, 0.707502, 0.063322, 0.027287, 0.009182, 0.197066)
W_222 <- W_122; lambdas_222 <- lambdas_122
c_and_gamma_222_2_1 <- c(0.1, 0.2)
df_222_2_1 <- 7
theta_222logistictsh_2_1 <- c(all_phi_222, all_A_222, vec(W_222), lambdas_222, c_and_gamma_222_2_1, df_222_2_1)
mod222logistictsh_2_1 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222logistictsh_2_1,
                               weight_function="logistic", weightfun_pars=c(2, 1),
                               cond_dist="Student", identification="heteroskedasticity")


## Structural models imposing constraints

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2),
# AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean"
theta_222expcmwtsh_2_1 <- c(0.7209658, 0.810858, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, # mu + A
                            -0.03, 0.24, -0.76, -0.02, 3.36, 0.86, # W + lambdas
                            0.33, 4) # xi + nu
mod222expcmwtsh_2_1 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222expcmwtsh_2_1, weight_function="exponential",
                             weightfun_pars=c(2, 1), identification="heteroskedasticity", cond_dist="Student",
                             mean_constraints=list(1:2), AR_constraints=C_222,
                             weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean")

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2),
# AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean",
# B_constraints=matrix(c(-0.03, 0.24, 0, -0.02), nrow=2, ncol=2)
theta_222expcmwbtsh_2_1 <- c(0.7209658, 0.810858, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, # mu + A
                            -0.03, 0.24, -0.02, 3.36, 0.86, # W + lambdas (excludes zero constr element)
                            0.33, 4) # xi + nu
mod222expcmwbtsh_2_1 <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222expcmwbtsh_2_1, weight_function="exponential",
                              weightfun_pars=c(2, 1), identification="heteroskedasticity", cond_dist="Student",
                              mean_constraints=list(1:2), AR_constraints=C_222,
                              weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)),
                              B_constraints=matrix(c(-0.03, 0.24, 0, -0.02), nrow=2, ncol=2), parametrization="mean")

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="ind_Student", mean_constraints=list(1:2),
# AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean",
# B_constraints=matrix(c(1, NA, 0, 1), nrow=2, ncol=2)
theta_222logistitb <- c(0.7209658, 0.810858, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, # mu + A
                        0.1, 0.2, 0.7, 0.11, -0.22, 0.73, # B mats
                        0.4, 7, 3)
mod222logistitb <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222logistitb, weight_function="logistic", weightfun_pars=c(2, 1),
                        cond_dist="ind_Student", mean_constraints=list(1:2), AR_constraints=C_222,
                        weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), identification="non-Gaussianity",
                        parametrization="mean", B_constraints=matrix(c(1, NA, 0, 1), nrow=2, ncol=2))

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="ind_skewed_t", mean_constraints=list(1:2),
# AR_constraints=C_222, weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean",
# B_constraints=matrix(c(1, NA, 0, 1), nrow=2, ncol=2)
theta_222logistiktb <- c(0.7209658, 0.810858, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, # mu + A
                        0.1, 0.2, 0.7, 0.11, -0.22, 0.73, # B mats
                        0.4, 7, 3, 0.4, 0.7)
mod222logistiktb <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222logistiktb, weight_function="logistic", weightfun_pars=c(2, 1),
                         cond_dist="ind_skewed_t", mean_constraints=list(1:2), AR_constraints=C_222,
                         weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), identification="non-Gaussianity",
                         parametrization="mean", B_constraints=matrix(c(1, NA, 0, 1), nrow=2, ncol=2))


test_that("STVAR works correctly", {
  # Relative_dens Gaussian STVAR
  expect_equal(mod112relg$params, theta_112relg)
  expect_equal(mod122relg$params, theta_122relg)
  expect_equal(mod222relg$params, theta_222relg)
  expect_equal(mod123relg$params, theta_123relg)
  expect_equal(mod323relg$params, theta_323relg)
  expect_equal(mod222relgcm$params, theta_222relgcm)

  # Logistic
  expect_equal(mod222logisticcmw_2_1$params, theta_222logisticcmw_2_1)

  # Logit
  expect_equal(mod222logcm_12_2$params, theta_222logcm_12_2)
  expect_equal(mod222logcmw_12_2$params, theta_222logcmw_12_2)

  # Exponential
  expect_equal(mod222expcmw_2_1$params, theta_222expcmw_2_1)

  # Threshold
  expect_equal(mod222thres_2_1$params, theta_222thres_2_1)

  # Exogenous
  expect_equal(mod222exo$params, theta_222exo)

  # Student
  expect_equal(mod222logisticcmwt_2_1$params, theta_222logisticcmwt_2_1)
  expect_equal(mod222logcmt_12_2$params, theta_222logcmt_12_2)
  expect_equal(mod222expcmwt_2_1$params, theta_222expcmwt_2_1)
  expect_equal(mod222threst_2_1$params, theta_222threst_2_1)

  # ind_Student
  expect_equal(mod123exoit$params, theta_123exoit)
  expect_equal(mod222logistit$params, theta_222logistit)

  # ind_skewed_t
  expect_equal(mod123exoikt$params, theta_123exoikt)
  expect_equal(mod222logistikt$params, theta_222logistikt)

  # Structural
  expect_equal(mod222thressr_2_1$params, theta_222thres_2_1)
  expect_equal(mod122relgsh$params, theta_122relgsh)
  expect_equal(mod222logistictsh_2_1$params, theta_222logistictsh_2_1)
  expect_equal(mod222expcmwtsh_2_1$params, theta_222expcmwtsh_2_1)
  expect_equal(mod222expcmwbtsh_2_1$params, theta_222expcmwbtsh_2_1)
  expect_equal(mod222logistitb$params, theta_222logistitb)
  expect_equal(mod222logistiktb$params, theta_222logistiktb)
})

test_that("swap_parametrization works correctly", {
  # Linear Gaussian STVAR
  expect_equal(swap_parametrization(mod112relg, calc_std_errors=FALSE)$params,
               c(0.7499743, 0.8049962, 0.2885260, 0.0217670, -0.1440240, 0.8971030,
                 0.6017860, -0.0029450, 0.0672240), tolerance=1e-3)

  # Relative dens Gaussian STVAR
  expect_equal(swap_parametrization(mod222relg, calc_std_errors=FALSE)$params,
               c(0.960032, 0.554909, 0.490841, 1.1693, 0.13996, 0.035172, -0.164575,
                 0.386816, 0.451675, 0.013086, 0.227882, 0.336084, 0.239257, 0.024173,
                 -0.021209, 0.707502, 0.063322, 0.027287, 0.009182, 0.197066, 0.205831,
                 0.005157, 0.025877, 1.092094, -0.009327, 0.116449, 0.592446), tolerance=1e-3)

  # Student Threshold
  expect_equal(swap_parametrization(mod222threst_2_1, calc_std_errors=FALSE)$params,
               c(0.960032, 0.554909, 0.490841, 1.1693, 0.13996, 0.035172, -0.164575,
                 0.386816, 0.451675, 0.013086, 0.227882, 0.336084, 0.239257, 0.024173,
                 -0.021209, 0.707502, 0.063322, 0.027287, 0.009182, 0.197066, 0.205831,
                 0.005157, 0.025877, 1.092094, -0.009327, 0.116449, 1, 13), tolerance=1e-3)

  # Structural logistic
  expect_equal(swap_parametrization(mod222logistictsh_2_1, calc_std_errors=FALSE)$params,
               c(0.9600324, 0.5549088, 0.4908410, 1.1693004, 0.1399600, 0.0351720, -0.1645750,
                 0.3868160, 0.4516750, 0.0130860, 0.2278820, 0.3360840, 0.2392570, 0.0241730,
                 -0.0212090, 0.7075020, 0.0633220, 0.0272870, 0.0091820, 0.1970660, -0.0300000,
                 0.2400000, -0.7600000, -0.0200000, 3.3600000, 0.8600000, 0.1000000, 0.2000000, 7.0000000),
               tolerance=1e-3)

  # ind_Student exo
  expect_equal(swap_parametrization(mod123exoit, calc_std_errors=FALSE)$params,
               c(0.01872735, 0.54692881, 2.82960017, -0.50345541, 1.88535674, 9.50476379, 0.87227000, -0.01595000, 0.14124000,
                 -0.08611000, 0.61865000, 0.34311000, -0.02047000, 0.02500000, 0.97548000, 0.74976000, 0.02187000, 0.29213000,
                 -1.55165000, 0.58245000, -0.00696000, -0.07261000, 0.02021000, 0.96883000, 0.21700000, -0.54200000, 0.89100000,
                 0.59600000, 1.63600000, 0.68900000, -1.28100000, -0.21300000, 1.89700000, 1.77700000, 0.56700000, 0.01600000,
                 0.38300000, -0.04500000, 0.03400000, 0.16900000, 1.16500000, -0.04400000, 7.00000000, 3.00000000, 13.00000000),
               tolerance=1e-3)

  # ind_skewed_t exo
  expect_equal(swap_parametrization(mod123exoikt, calc_std_errors=FALSE)$params,
               c(0.01872735, 0.54692881, 2.82960017, -0.50345541, 1.88535674, 9.50476379, 0.87227000, -0.01595000, 0.14124000,
                 -0.08611000, 0.61865000, 0.34311000, -0.02047000, 0.02500000, 0.97548000, 0.74976000, 0.02187000, 0.29213000,
                 -1.55165000, 0.58245000, -0.00696000, -0.07261000, 0.02021000, 0.96883000, 0.21700000, -0.54200000, 0.89100000,
                 0.59600000, 1.63600000, 0.68900000, -1.28100000, -0.21300000, 1.89700000, 1.77700000, 0.56700000, 0.01600000,
                 0.38300000, -0.04500000, 0.03400000, 0.16900000, 1.16500000, -0.04400000, 7.00000000, 3.00000000, 13.00000000,
                 -0.1, 0.2, 0.3), tolerance=1e-3)
})



test_that("get_hetsked_sstvar works correctly", {
  # Relative_dens Gaussian STVAR
  expect_equal(get_hetsked_sstvar(mod123relg)$params,
               c(0.10741, 0.13813, -0.12092, 3.48957, 0.60615, 0.45646, 0.87227, -0.01595, 0.14124, -0.08611,
                 0.61865, 0.34311, -0.02047, 0.025, 0.97548, 0.74976, 0.02187, 0.29213, -1.55165, 0.58245,
                 -0.00696, -0.07261, 0.02021, 0.96883, 0.0325296, 0.0298782, 0.339136, -0.0670146, 0.2290013,
                 -0.0541424, 0.8099018, 0.0458877, 0.0955791, 29.4434945, 3.98537, 1.1237683, 0.77357), tolerance=1e-4)

  # Logistic
  expect_equal(get_hetsked_sstvar(mod222logisticcmw_2_1)$params,
               c(0.7209658, 0.810858, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, 0.435331,
                 -0.0329379, 0.1431325, 0.1700444, 5.4150627, 3.6010663, 0.33), tolerance=1e-4)

  # Logit
  expect_equal(get_hetsked_sstvar(mod222logcm_12_2)$params,
               c(0.7209658, 0.810858, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, 0.435331,
                 -0.0329379, 0.1431325, 0.1700444, 5.4150627, 3.6010663,0.10, 0.20, 0.11, 0.22, 0.33), tolerance=1e-4)

  # Exponential
  expect_equal(get_hetsked_sstvar(mod222expcmw_2_1)$params,
               c(0.7209658, 0.810858, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, 0.435331,
                 -0.0329379, 0.1431325, 0.1700444, 5.4150627, 3.6010663, 0.33), tolerance=1e-4)

  # Student (with Threshold)
  expect_equal(get_hetsked_sstvar(mod222threst_2_1)$params,
               c(0.356914, 0.107436, 0.356386, 0.08633, 0.13996, 0.035172, -0.164575, 0.386816, 0.451675, 0.013086,
                 0.227882, 0.336084, 0.239257, 0.024173, -0.021209, 0.707502, 0.063322, 0.027287, 0.009182, 0.197066,
                 0.4046369, -0.0624318, 0.2051828, 0.1482541, 5.5613353, 4.3119004, 1, 13), tolerance=4)
})

test_that("reorder_B_columns works correctly", {
  expect_equal(reorder_B_columns(mod222logistictsh_2_1, perm=1:2, calc_std_errors=FALSE)$params,
               mod222logistictsh_2_1$params, tolerance=1e-4)
  expect_equal(reorder_B_columns(mod222logistictsh_2_1, perm=2:1, calc_std_errors=FALSE)$params,
               c(0.356914, 0.107436, 0.356386, 0.086330, 0.139960, 0.035172, -0.164575, 0.386816, 0.451675, 0.013086,
                 0.227882, 0.336084, 0.239257, 0.024173, -0.021209, 0.707502, 0.063322, 0.027287, 0.009182, 0.197066,
                 -0.760000, -0.020000, -0.030000, 0.240000, 0.860000, 3.360000, 0.100000, 0.200000, 7.000000), tolerance=1e-4)

  expect_equal(reorder_B_columns(mod222expcmwbtsh_2_1, perm=2:1, calc_std_errors=FALSE)$params,
               c(0.7209658, 0.8108580, 0.2200000, 0.0600000, -0.1500000, 0.3900000, 0.4100000, -0.0100000, 0.0800000, 0.3000000,
                 -0.0200000, -0.0300000, 0.2400000, 0.8600000, 3.3600000, 0.3300000, 4.0000000), tolerance=1e-4)

  expect_equal(reorder_B_columns(mod123exoit, perm=c(3, 1, 2), calc_std_errors=FALSE)$params,
               c(0.10741, 0.13813, -0.12092, 3.48957, 0.60615, 0.45646, 0.87227, -0.01595, 0.14124, -0.08611, 0.61865, 0.34311,
                 -0.02047, 0.025, 0.97548, 0.74976, 0.02187, 0.29213, -1.55165, 0.58245, -0.00696, -0.07261, 0.02021, 0.96883,
                 -1.281, -0.213, 1.897, 0.217, -0.542, 0.891, 0.596, 1.636, 0.689, 0.169, 1.165, -0.044, 1.777, 0.567, 0.016,
                 0.383, -0.045, 0.034, 13, 7, 3), tolerance=1e-4)

  expect_equal(reorder_B_columns(mod222logistitb, perm=c(2, 1), calc_std_errors=FALSE)$params,
               c(0.720966, 0.810858, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, 0.7, 0.1,
                 0.2, 0.73, 0.11, -0.22, 0.4, 3, 7), tolerance=1e-4)

  expect_equal(reorder_B_columns(mod222logistitb, perm=c(2, 1), calc_std_errors=FALSE)$model$B_constraints,
               matrix(c(0, 1, 1, NA), nrow=2), tolerance=1e-4)

  # ind_skewed_t
  expect_equal(reorder_B_columns(mod123exoikt, perm=c(3, 1, 2), calc_std_errors=FALSE)$params,
               c(0.10741, 0.13813, -0.12092, 3.48957, 0.60615, 0.45646, 0.87227, -0.01595, 0.14124, -0.08611, 0.61865, 0.34311,
                 -0.02047, 0.025, 0.97548, 0.74976, 0.02187, 0.29213, -1.55165, 0.58245, -0.00696, -0.07261, 0.02021, 0.96883,
                 -1.281, -0.213, 1.897, 0.217, -0.542, 0.891, 0.596, 1.636, 0.689, 0.169, 1.165, -0.044, 1.777, 0.567, 0.016,
                 0.383, -0.045, 0.034, 13, 7, 3, 0.3, -0.1, 0.2), tolerance=1e-4)

  expect_equal(reorder_B_columns(mod222logistiktb, perm=c(2, 1), calc_std_errors=FALSE)$params,
               c(0.7209658, 0.810858, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, 0.7, 0.1, 0.2, 0.73, 0.11, -0.22,
                 0.4, 3, 7, 0.7, 0.4), tolerance=1e-4)

  expect_equal(reorder_B_columns(mod222logistiktb, perm=c(1, 2), calc_std_errors=FALSE)$params,
               mod222logistiktb$params, tolerance=1e-4)

  expect_equal(reorder_B_columns(mod222logistiktb, perm=c(2, 1), calc_std_errors=FALSE)$model$B_constraints,
               matrix(c(0, 1, 1, NA), nrow=2), tolerance=1e-4)
})


test_that("swap_B_signs works correctly", {
  expect_equal(swap_B_signs(mod122relgsh, which_to_swap=2, calc_std_errors=FALSE)$params,
               c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863, -0.312400, 0.505251, 0.298483, 0.030096, -0.176925,
                 0.838898, -0.030000, 0.240000, 0.760000, 0.020000, 3.360000, 0.860000, 0.600000), tolerance=1e-4)

  expect_equal(swap_B_signs(mod222expcmwbtsh_2_1, which_to_swap=c(2, 1), calc_std_errors=FALSE)$params,
               c(0.7209658, 0.8108580, 0.2200000, 0.0600000, -0.1500000, 0.3900000, 0.4100000, -0.0100000, 0.0800000, 0.3000000,
                 0.0300000, -0.2400000, 0.0200000, 3.3600000, 0.8600000, 0.3300000, 4.0000000), tolerance=1e-4)

  expect_equal(swap_B_signs(mod222expcmwbtsh_2_1, which_to_swap=1, calc_std_errors=FALSE)$params,
               c(0.7209658, 0.8108580, 0.2200000, 0.0600000, -0.1500000, 0.3900000, 0.4100000, -0.0100000, 0.0800000, 0.3000000,
                 0.0300000, -0.2400000, -0.0200000, 3.3600000, 0.8600000, 0.3300000, 4.0000000), tolerance=1e-4)

  expect_equal(swap_B_signs(mod123exoit, which_to_swap=1, calc_std_errors=FALSE)$params,
               c(0.10741, 0.13813, -0.12092, 3.48957, 0.60615, 0.45646, 0.87227, -0.01595, 0.14124, -0.08611, 0.61865, 0.34311,
                 -0.02047, 0.025, 0.97548, 0.74976, 0.02187, 0.29213, -1.55165, 0.58245, -0.00696, -0.07261, 0.02021, 0.96883,
                 -0.217, 0.542, -0.891, 0.596, 1.636, 0.689, -1.281, -0.213, 1.897, -1.777, -0.567, -0.016, 0.383, -0.045, 0.034,
                 0.169, 1.165, -0.044, 7, 3, 13), tolerance=1e-4)

  expect_equal(swap_B_signs(mod222logistitb, which_to_swap=2, calc_std_errors=FALSE)$params,
               c(0.720966, 0.810858, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, 0.1, 0.2, -0.7, 0.11, -0.22, -0.73, 0.4, 7, 3),
               tolerance=1e-4)

  expect_equal(swap_B_signs(mod222logistitb, which_to_swap=2, calc_std_errors=FALSE)$model$B_constraints,
               matrix(c(1, NA, 0, -1), nrow=2), tolerance=1e-4)
  expect_equal(swap_B_signs(mod222logistitb, which_to_swap=1:2, calc_std_errors=FALSE)$model$B_constraints,
               matrix(c(-1, NA, 0, -1), nrow=2), tolerance=1e-4)
  expect_equal(swap_B_signs(mod222logistitb, which_to_swap=1, calc_std_errors=FALSE)$model$B_constraints,
               matrix(c(-1, NA, 0, 1), nrow=2), tolerance=1e-4)

  # ind_skewed_t
  expect_equal(swap_B_signs(mod123exoikt, which_to_swap=1, calc_std_errors=FALSE)$params,
               c(0.10741, 0.13813, -0.12092, 3.48957, 0.60615, 0.45646, 0.87227, -0.01595, 0.14124, -0.08611, 0.61865, 0.34311,
                 -0.02047, 0.025, 0.97548, 0.74976, 0.02187, 0.29213, -1.55165, 0.58245, -0.00696, -0.07261, 0.02021, 0.96883,
                 -0.217, 0.542, -0.891, 0.596, 1.636, 0.689, -1.281, -0.213, 1.897, -1.777, -0.567, -0.016, 0.383, -0.045, 0.034,
                 0.169, 1.165, -0.044, 7, 3, 13, -0.10000, 0.2, 0.3), tolerance=1e-4)

  expect_equal(swap_B_signs(mod222logistiktb, which_to_swap=2, calc_std_errors=FALSE)$params,
               c(0.720966, 0.810858, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, 0.1, 0.2, -0.7, 0.11, -0.22, -0.73, 0.4, 7, 3,
                 0.4, 0.7), tolerance=1e-4)

  expect_equal(swap_B_signs(mod222logistiktb, which_to_swap=2, calc_std_errors=FALSE)$model$B_constraints,
               matrix(c(1, NA, 0, -1), nrow=2), tolerance=1e-4)
  expect_equal(swap_B_signs(mod222logistiktb, which_to_swap=1:2, calc_std_errors=FALSE)$model$B_constraints,
               matrix(c(-1, NA, 0, -1), nrow=2), tolerance=1e-4)
  expect_equal(swap_B_signs(mod222logistiktb, which_to_swap=1, calc_std_errors=FALSE)$model$B_constraints,
               matrix(c(-1, NA, 0, 1), nrow=2), tolerance=1e-4)
})
