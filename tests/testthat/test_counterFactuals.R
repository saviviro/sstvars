context("counterFactuals")
library(sstvars)

# p=1, M=1, d=2
theta_112relg <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103, 0.601786, -0.002945, 0.067224)
mod112relg <- STVAR(data=gdpdef[1:50,], p=1, M=1, params=theta_112relg, weight_function="relative_dens")

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
mod123relg <- STVAR(data=usamone[1:30,], p=1, M=2, params=theta_123relg, weight_function="relative_dens")

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
mod222thres_2_1 <- STVAR(data=gdpdef[1:30,], p=2, M=2, d=2, params=theta_222thres_2_1, weight_function="threshold",
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
mod222logcmt_12_2 <- STVAR(data=gdpdef[1:60,], p=2, M=2, d=2, params=theta_222logcmt_12_2, weight_function="mlogit", parametrization="mean",
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
mod222expcmwtsh_2_1 <- STVAR(data=gdpdef[1:25,], p=2, M=2, d=2, params=theta_222expcmwtsh_2_1, weight_function="exponential",
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
mod222logistitb <- STVAR(data=gdpdef[1:55,], p=2, M=2, d=2, params=theta_222logistitb, weight_function="logistic", weightfun_pars=c(2, 1),
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

test_that("get_phi_yt works correctly", {
  # Simple (2x2) matrix
  all_phi0 <- matrix(c(1, 2, 3, 4),  nrow=2, ncol=2, byrow=FALSE)
  alpha_mt <- c(0.3, 0.7)
  result <- get_phi_yt(all_phi0=all_phi0, alpha_mt=alpha_mt)
  expect_type(result, "double")
  expect_length(result, 2)
  expect_equal(result, c(1*0.3 + 3*0.7, 2*0.3 + 4*0.7), tolerance=1e-6)

  # Zero weights give zero vector
  set.seed(42); all_phi0 <- matrix(rnorm(6), nrow=3, ncol=2)
  expect_equal(get_phi_yt(all_phi0=all_phi0, alpha_mt=c(0, 0)), rep(0, 3), tolerance=1e-6)

  # Identity intercepts reproduces alpha: if all_phi0 is identity, output == alpha_mt
  alpha_mt <- c(0.1, 0.5, 0.4)
  expect_equal(get_phi_yt(all_phi0=diag(3), alpha_mt=alpha_mt), alpha_mt, tolerance=1e-6)

  # Test with random larger matrices:
  set.seed(123); d <- 4; M <- 5
  all_phi0 <- matrix(rnorm(d*M), nrow=d, ncol=M)
  alpha_mt <- runif(M)
  expect_equal(get_phi_yt(all_phi0, alpha_mt), as.numeric(all_phi0%*%alpha_mt), tolerance=1e-6)
})


test_that("get_allA_yti works correctly", {
  ## Simple (2x2x2x2) matrix
  all_A <- array(NA, dim=c(2, 2, 2, 2))

  # Regime 1
  all_A[, , 1, 1] <- matrix(c(1, 2, 3, 4), nrow=2, byrow=TRUE)
  all_A[, , 2, 1] <- matrix(c(5, 6, 7, 8), nrow=2, byrow=TRUE)

  # Regime 2
  all_A[, , 1, 2] <- matrix(c(2, 1, 0, 1), nrow=2, byrow=TRUE)
  all_A[, , 2, 2] <- matrix(c(0, 1, 1, 0), nrow=2, byrow=TRUE)
  alpha_mt <- c(0.6, 0.4)

  result <- get_allA_yti(all_A, alpha_mt)

  # lag 1: 0.6*regime1[, , 1] + 0.4*regime2[, , 1]
  expected1 <- 0.6*all_A[, , 1, 1] + 0.4*all_A[, , 1, 2]
  # lag 2: similarly
  expected2 <- 0.6*all_A[, , 2, 1] + 0.4*all_A[, , 2, 2]

  expect_equal(dim(result), c(2, 2, 2))
  expect_equal(result[, , 1], expected1, tolerance=1e-6)
  expect_equal(result[, , 2], expected2, tolerance=1e-6)

  ## (3x3x2x2) case
  # Regime 1: A_{1,i} = i1*I_3
  # Regime 2: A_{2,i} = 2*i1*I_3
  all_A <- array(0, dim = c(3, 3, 2, 2))
  for(i1 in 1:2) {
    all_A[, , i1, 1] <- i1*diag(3)
    all_A[, , i1, 2] <- 2*i1*diag(3)
  }
  alpha_mt <- c(0.2, 0.8)

  result <- get_allA_yti(all_A, alpha_mt)

  # expected lag i: 0.2*(i1*I) + 0.8*(2*i1*I) = (0.2 + 1.6)*i1 * I = 1.8*i1*I
  expected1 <- 1.8*1*diag(3)
  expected2 <- 1.8*2*diag(3)

  expect_equal(dim(result), c(3, 3, 2))
  expect_equal(result[, , 1], expected1, tolerance=1e-6)
  expect_equal(result[, , 2], expected2, tolerance=1e-6)

  ## (2x2x1x3) case
  all_A <- array(NA, dim=c(2, 2, 1, 3))
  # Regime 1
  all_A[, , 1, 1] <- matrix(c(1, 2, 3, 4), nrow=2, byrow=TRUE)
  # Regime 2
  all_A[, , 1, 2] <- matrix(c(2, 0, 0, 2), nrow=2, byrow=TRUE)
  # Regime 3
  all_A[, , 1, 3] <- matrix(c(0, 1, 1, 0), nrow=2, byrow=TRUE)

  alpha_mt <- c(0.2, 0.3, 0.5)
  result <- get_allA_yti(all_A, alpha_mt)

  expect_equal(dim(result), c(2, 2, 1))
  # expected = 0.2*regime1 + 0.3*regime2 + 0.5*regime3
  expected <- 0.2*all_A[, , 1, 1] + 0.3*all_A[, , 1, 2] + 0.5*all_A[, , 1, 3]
  expect_equal(result[ , , 1], expected, tolerance=1e-6)

  ## (3×3×2×1) case
  all_A <- array(NA, dim=c(3, 3, 2, 1))
  # lag 1 matrix
  all_A[, , 1, 1] <- matrix(c(1, 0, 2, 0, 3, 0, 4, 0, 5), nrow=3, byrow=TRUE)
  # lag 2 matrix
  all_A[, , 2, 1] <- matrix(c(5, 4, 3, 2, 1, 0, 0, 1, 2), nrow=3, byrow=TRUE)

  alpha_mt <- 1
  result <- get_allA_yti(all_A, alpha_mt)

  expect_equal(dim(result), c(3, 3, 2))
  expect_equal(result[, , 1], all_A[, , 1, 1], tolerance=1e-6)
  expect_equal(result[, , 2], all_A[, , 2, 1], tolerance=1e-6)
})


test_that("get_mu_yt works correctly", {
  # d=2, p=2 case
  phi_yt <- c(1, 2)
  all_A_yti <- array(0, dim=c(2, 2, 2))
  all_A_yti[, , 1] <- diag(2) # A_{y,t,1} = I
  all_A_yti[, , 2] <- matrix(c(0, 2, 3, 0), nrow=2, byrow=TRUE) # A_{y,t,2}
  bold_y_t_minus_1 <- c(1, 1, 2, 3) # [y_{t-1}; y_{t-2}]

  # Compute result
  result <- get_mu_yt(phi_yt=phi_yt, all_A_yti=all_A_yti, bold_y_t_minus_1=bold_y_t_minus_1)

  # Manual:
  # big_A = [I | A2] = [1 0 0 2; 0 1 3 0]
  # big_A %*% bold_y = c(1*1 + 0*1 + 0*2 + 2*3,
  #                      0*1 + 1*1 + 3*2 + 0*3) = c(7 ,7)
  # mu = phi0 + that = c(8, 9)
  expect_equal(result, c(8, 9), tolerance=1e-8)

  # d=3, p=1 case
  d <- 3; p <- 1
  phi_yt <- c(-1, 0, 2)
  all_A_yti <- array(0, dim=c(d, d, p))
  all_A_yti[, , 1] <- 2*diag(3)
  bold_y_t_minus_1 <- c(1, 2, 3)

  # mu = phi0 + 2*bold_y = c(-1+2, 0+4, 2+6) = c(1,4,8)
  result <- get_mu_yt(phi_yt=phi_yt, all_A_yti=all_A_yti, bold_y_t_minus_1=bold_y_t_minus_1)
  expect_equal(result, c(1, 4, 8), tolerance=1e-8)

  # d=4, p=3 case
  set.seed(42)
  d <- 4; p <- 3
  phi_yt <- runif(d)
  all_A_yti <- array(rnorm(d*d*p), dim=c(d, d, p))
  bold_y_t_minus_1 <- rnorm(d*p)

  # vectorized result
  vec_res <- get_mu_yt(phi_yt=phi_yt, all_A_yti=all_A_yti, bold_y_t_minus_1=bold_y_t_minus_1)

  # loop result
  loop_res <- phi_yt
  for(i1 in seq_len(p)) {
    ylag <- bold_y_t_minus_1[((i1 - 1)*d + 1):(i1*d)]
    loop_res <- loop_res + all_A_yti[, , i1]%*%ylag
  }
  expect_equal(vec_res, as.numeric(loop_res), tolerance=1e-6)

  # d = 2, p = 2
  phi_yt <- c(0.5, -0.5)
  all_A_yti <- array(0, dim=c(2, 2, 2))
  # A_{y,t,1} non-diagonal
  all_A_yti[, , 1] <- matrix(c(1, 2, 3, 4), nrow=2, byrow=TRUE)
  # A_{y,t,2} also non-diagonal
  all_A_yti[, , 2] <- matrix(c(-1, 0, 0, 1), nrow=2, byrow=TRUE)
  # y_{t-1} = (1,1), y_{t-2} = (-1,2)
  bold_y_t_minus_1 <- c(1, 1, -1, 2)

  result <- get_mu_yt(phi_yt=phi_yt, all_A_yti=all_A_yti, bold_y_t_minus_1=bold_y_t_minus_1)

  # manual:
  # A1 %*% c(1,1) = c(3,7)
  # A2 %*% c(-1,2) = c(1,2)
  # sum = c(4,9), + phi0 = c(4.5, 8.5)
  expect_equal(result, c(4.5, 8.5), tolerance=1e-6)
})

test_that("get_B_yt works correctly", {
  ## non-Gaussianity

  # d=2, M=2
  all_Bm <- array(NA, dim=c(2, 2, 2))
  all_Bm[, , 1] <- matrix(c(1, 2, 3, 4), nrow=2, byrow=TRUE)
  all_Bm[, , 2] <- matrix(c(5, 6, 7, 8), nrow=2, byrow=TRUE)
  alpha_mt <- c(0.3, 0.7)

  expect_equal(get_B_yt(all_Omegas=all_Bm, alpha_mt=alpha_mt, cond_dist="ind_Student"),
               0.3*all_Bm[, , 1] + 0.7*all_Bm[, , 2], tolerance=1e-6)

  # d=2, M=3
  all_Bm <- array(NA, dim = c(2, 2, 3))
  all_Bm[, , 1] <- matrix(c(1, 0, 0, 1), nrow=2, byrow=TRUE)
  all_Bm[, , 2] <- matrix(c(2, 1, 1, 2), nrow=2, byrow=TRUE)
  all_Bm[, , 3] <- matrix(c(0, 1, 2, 0), nrow=2, byrow=TRUE)
  alpha_mt <- c(0.2, 0.3, 0.5)

  expect_equal(get_B_yt(all_Omegas=all_Bm, alpha_mt=alpha_mt, cond_dist="ind_skewed_t"),
               0.2*all_Bm[, , 1] + 0.3*all_Bm[, , 2] + 0.5*all_Bm[, , 3], tolerance=1e-6)

  # d=3, M=2
  all_Bm <- array(NA, dim=c(3, 3, 2))
  all_Bm[, , 1] <- diag(3)
  all_Bm[, , 2] <- 2*diag(3)
  alpha_mt <- c(0.4, 0.6)

  # expected = 0.4*I + 0.6*(2I) = (0.4 + 1.2) * I = 1.6 * I
  expect_equal(get_B_yt(all_Omegas=all_Bm, alpha_mt=alpha_mt, cond_dist="ind_skewed_t"), 1.6*diag(3), tolerance=1e-6)

  ## heteroskedasticity
  # d=2, M=3
  d <- 2; M <- 3
  all_Omegas <- array(0, dim=c(d, d, M)) # not used in this test
  alpha_mt <- c(0.2, 0.3, 0.5)
  W <- matrix(c(1, 0, 0, 2), nrow=d)
  lambdas <- c(4, 6, 5, 7)

  result <- get_B_yt(all_Omegas=all_Omegas, alpha_mt=alpha_mt, W=W, lambdas=lambdas, cond_dist="Gaussian",
                     identification="heteroskedasticity")

  Lambda_sum <- alpha_mt[1]*diag(d) + alpha_mt[2]*diag(c(4, 6)) + alpha_mt[3]*diag(c(5, 7)) # manual calculation
  expected <- W%*%sqrt(Lambda_sum)
  expect_equal(result, expected, tolerance=1e-6)

  # d=3, M=2
  d <- 3; M <- 2
  all_Omegas <- array(0, dim=c(d, d, M))
  alpha_mt   <- c(0.6, 0.4)
  W <- matrix(c(2, 0, 0, 0, 1, 0, 1, 0, 3), nrow=d, byrow=TRUE)
  lambdas <- c(9, 4, 1)
  B_yt <- get_B_yt(all_Omegas=all_Omegas, alpha_mt=alpha_mt, W=W, lambdas=lambdas, cond_dist="Student",
                   identification="heteroskedasticity")

  # Manual calculation:
  Lambda_sum <- alpha_mt[1]*diag(d) + alpha_mt[2]*diag(lambdas)
  expected <- W%*%sqrt(Lambda_sum)
  expect_equal(B_yt, expected, tolerance=1e-6)

  ## recursive
  # d=3, M=2
  d <- 3; M <- 2
  all_Omegas <- array(NA, dim=c(d, d, M))

  # regime 1
  all_Omegas[, , 1] <- matrix(c(4, 1, 0, 1, 5, 2, 0, 2, 6), nrow=d, byrow=TRUE)

  # regime 2
  all_Omegas[, , 2] <- matrix(c(9, 0, 1, 0, 7, 0, 1, 0, 10), nrow=d, byrow=TRUE)
  alpha_mt <- c(0.7, 0.3)

  result <- get_B_yt(all_Omegas=all_Omegas, alpha_mt=alpha_mt, W=NULL, lambdas=NULL, cond_dist="Gaussian",
                     identification="recursive")

  # Manual calculation:
  Omega_t <- alpha_mt[1]*all_Omegas[, , 1] + alpha_mt[2]*all_Omegas[, , 2]
  expected <- t(chol(Omega_t)) # lower-triangular

  expect_equal(result, expected, tolerance=1e-6)

  # additional check:  B_t%*%t(B_t) == Ω_t
  expect_equal(tcrossprod(result, result), Omega_t, tolerance=1e-6)

  # d=2, M=3
  d <- 2; M <- 3
  all_Omegas <- array(NA, dim=c(d, d, M))
  all_Omegas[, , 1] <- matrix(c(4, 1, 1, 3), nrow=2, byrow=TRUE)
  all_Omegas[, , 2] <- matrix(c(10, 2, 2, 8), nrow=2, byrow=TRUE)
  all_Omegas[, , 3] <- matrix(c(6, 0, 0, 5), nrow=2, byrow=TRUE)
  alpha_mt <- c(0.2, 0.5, 0.3)

  B_yt <- get_B_yt(all_Omegas=all_Omegas, alpha_mt=alpha_mt, W=NULL, lambdas=NULL, cond_dist="Student",
                  identification="reduced_form")

  # expected
  Omega_t  <- alpha_mt[1]*all_Omegas[, , 1] + alpha_mt[2]*all_Omegas[, , 2] + alpha_mt[3]*all_Omegas[, , 3]
  expected <- t(chol(Omega_t)) # lower-triangular

  expect_equal(B_yt, expected, tolerance=1e-6)
  expect_equal(B_yt%*%t(B_yt), Omega_t, tolerance=1e-6)
})


test_that("cfact_hist works correctly", {

  # Linear SVAR
  mod <- mod112relg
  tmp <- cfact_hist(mod, cfact_type="fixed_path", policy_var=1, cfact_start=1, cfact_end=1, cfact_path=c(13))
  expect_equal(c(tmp$cfact_data[1:3,]), c(1.9195500, 13.0000000, 3.1782960, 0.2331600, 0.1004573, 0.5717706), tolerance=1e-3)
  expect_equal(tmp$cfact_alpha_mt[1:3], c(1, 1, 1), tolerance=1e-3)
  expect_equal(c(tmp$cfact_e_t[1:2,]), c(15.2500610, -1.5566894, -0.6137239, 0.4870329), tolerance=1e-3)

  # Relative_dens Gaussian STVAR
  mod <- mod123relg
  tmp <- cfact_hist(mod, cfact_type="fixed_path", policy_var=2, cfact_start=3, cfact_end=4, cfact_path=c(1, -1))
  expect_equal(c(tmp$cfact_data[4:6,]), c(0.58931119, 0.89152424, 1.15595397, 1.00000000, -1.00000000, 0.06629706,
                                          1.49926178, 1.88890671, 2.23420120), tolerance=1e-3)
  expect_equal(c(tmp$cfact_alpha_mt[4:6,]), c(0.907111044, 0.866407293, 0.996805256, 0.092888956, 0.133592707, 0.003194744), tolerance=1e-3)
  expect_equal(c(tmp$cfact_e_t[4:6,]), c(0.2715842, -0.4813848, -1.4705232, -6.8484605, 1.6189520, 0.8025761,
                                         0.9507210, 0.6733533, 0.1111197), tolerance=1e-3)

  # Logistic
  mod <- mod222logistitb
  tmp <- cfact_hist(mod, cfact_type="muted_response", policy_var=1, mute_var=2, cfact_start=2, cfact_end=3)
  expect_equal(c(tmp$cfact_data[2:5,]), c(2.2551700, 0.0705600, 0.5191495, 2.5195128, 0.1530400, 0.3850600, 0.5493942, 0.4989937), tolerance=1e-3)
  expect_equal(c(tmp$cfact_alpha_mt[2:5,]), c(0.4625642, 0.4462689, 0.4512559, 0.4446424, 0.5374358, 0.5537311, 0.5487441, 0.5553576), tolerance=1e-3)
  expect_equal(c(tmp$cfact_e_t[2:5,]), c(-6.63344471, 19.94040937, -11.92495765, -4.20250133, -0.02000102, 0.79938909, -0.72334832, 0.06099980), tolerance=1e-3)

  # Logit
  mod <- mod222logcmt_12_2
  tmp <- cfact_hist(mod, cfact_type="muted_response", policy_var=2, mute_var=1, cfact_start=2, cfact_end=2)
  expect_equal(c(tmp$cfact_data[2:4,]), c(2.2551700, 0.0705600, 0.4736337, 0.1530400, 0.3850600, 0.5935747), tolerance=1e-3)
  expect_equal(c(tmp$cfact_alpha_mt[2:4,]), c(0.6688768, 0.5993153, 0.7225385, 0.3311232, 0.4006847, 0.2774615), tolerance=1e-3)
  expect_equal(c(tmp$cfact_e_t[2:4,]), c(-1.0479836, 2.6930322, -1.7982568, 0.9072007, 0.1117878, -0.7601725), tolerance=1e-3)

  # Exponential
  mod <- mod222expcmwtsh_2_1
  tmp <- cfact_hist(mod, cfact_type="fixed_path", policy_var=2, cfact_start=1, cfact_end=2, cfact_path=c(-1, -2))
  expect_equal(c(tmp$cfact_data[2:5,]), c(2.255170, -58.633289, 43.205010, -12.043863, 0.153040, -1.000000, -2.000000, 2.380093), tolerance=1e-3)
  expect_equal(c(tmp$cfact_alpha_mt[1:3,]), c(0.993270797, 0.714170893, 0.263623587, 0.006729203, 0.285829107, 0.736376413), tolerance=1e-3)
  expect_equal(c(tmp$cfact_e_t[1:3,]), c(0.51127304, 0.73815826, 0.03739836, 79.27445159, -73.48881182, -2.70788018), tolerance=1e-3)

  # Threshold
  mod <- mod222thres_2_1
  tmp <- cfact_hist(mod, cfact_type="muted_response", policy_var=1, mute_var=2, cfact_start=10, cfact_end=10)
  expect_equal(c(tmp$cfact_data[11:13,]), c(1.9195900, 1.9501847, 1.7821673, 0.2562900, 0.3148065, 0.5151003), tolerance=1e-3)
  expect_equal(c(tmp$cfact_alpha_mt[9:11,]), c(1, 1, 1, 0, 0, 0), tolerance=1e-3)
  expect_equal(c(tmp$cfact_e_t[9:11,]), c(2.2251009, 1.2056860, 0.6142575, -0.6701476, -0.4580799, 0.6173660), tolerance=1e-3)

  # Exogenous
  mod <- mod123exoikt
  tmp <- cfact_hist(mod, cfact_type="fixed_path", policy_var=3, cfact_start=10, cfact_end=10, cfact_path=c(-5))
  expect_equal(c(tmp$cfact_data[10:12,]), c(-0.5522708, 4.3789027, 4.6584843, 0.4013990, -0.6250584, -0.6618279, 2.9266667,
                                            -5.0000000 ,-3.9782313), tolerance=1e-3)
  expect_equal(c(tmp$cfact_alpha_mt[8:10,]), c(0.2946009, 0.5776099, 0.6309793, 0.7053991, 0.4223901, 0.3690207), tolerance=1e-3)
  expect_equal(c(tmp$cfact_e_t[8:10,]), c(-1.8250055, 0.6423156, -1.2923466, 0.2980476, -0.8033733, 0.4184558, 0.8974124,
                                          0.1082131, -6.2418672), tolerance=1e-3)

})


test_that("cfact_fore works correctly", {
  # Linear SVAR
  mod <- mod112relg
  set.seed(1); tmp <- cfact_fore(mod, nsteps=2, nsim=1, pi=0.95, pred_type="mean", cfact_type="fixed_path", policy_var=1, cfact_start=1, cfact_end=1, cfact_path=c(13))
  expect_equal(c(tmp$cfact_pred[1:2,]), c(13.000000, 3.566238, 1.259185, 1.895067), tolerance=1e-3)
  expect_equal(tmp$cfact_pred_ints[c(1, 8)], c(13.000000, 1.895067), tolerance=1e-3)
  expect_equal(c(tmp$cfact_trans_pred[1:2,]), c(1, 1), tolerance=1e-3)
  expect_equal(c(tmp$cfact_trans_pred_ints[1:2, ,]), c(1, 1, 1, 1), tolerance=1e-3)
  expect_equal(c(tmp$pred[1:2,]), c(0.8743667, 1.1275558, 1.0433457, 1.2115593), tolerance=1e-3)
  expect_equal(tmp$pred_ints[c(1, 8)], c(0.8743667, 1.2115593), tolerance=1e-3)
  expect_equal(c(tmp$trans_pred[1:2,]), c(1, 1), tolerance=1e-3)
  expect_equal(c(tmp$trans_pred_ints[1:2, ,]), c(1, 1, 1, 1), tolerance=1e-3)

  # Relative_dens Gaussian STVAR
  mod <- mod123relg
  set.seed(2); tmp <- cfact_fore(mod, nsteps=3, nsim=3, pi=0.9, pred_type="median", cfact_type="muted_response", policy_var=3, mute_var=2, cfact_start=2, cfact_end=3)
  expect_equal(unname(c(tmp$cfact_pred[3,])), c(0.6479103, 0.4293333, 3.5803017), tolerance=1e-3)
  expect_equal(tmp$cfact_pred_ints[c(1, 11, 17)], c(1.010682, 0.826568, 4.279469), tolerance=1e-3)
  expect_equal(unname(c(tmp$cfact_trans_pred[3,])), c(0.994834639, 0.005165361), tolerance=1e-3)
  expect_equal(c(tmp$cfact_trans_pred_ints[c(1, 9)]), c(0.994135488, 0.004216029), tolerance=1e-3)
  expect_equal(c(tmp$pred[c(2, 5)]), c(1.6822459, 0.1661626), tolerance=1e-3)

  # Logistic
  mod <- mod222logistitb
  set.seed(3); tmp <- cfact_fore(mod, nsteps=3, nsim=2, pi=0.9, pred_type="median", cfact_type="fixed_path", policy_var=2, cfact_start=2,
                                 cfact_end=3, cfact_path=c(1, -1))
  expect_equal(unname(c(tmp$cfact_pred[2:3,])), c(0.8880295, 1.0383006, 1.0000000, -1.0000000), tolerance=1e-3)
  expect_equal(tmp$cfact_pred_ints[c(1, 7, 9)], c(1.3338063, 0.9898941, -1.0000000), tolerance=1e-3)
  expect_equal(unname(c(tmp$cfact_trans_pred[3,])), c(0.4022738, 0.5977262), tolerance=1e-3)
  expect_equal(c(tmp$cfact_trans_pred_ints[c(2, 8)]), c(0.3544091, 0.5966858), tolerance=1e-3)
  expect_equal(c(tmp$pred[c(2, 5)]), c(0.9215264, 1.0140519), tolerance=1e-3)

  # Logit
  mod <- mod222logcmt_12_2
  set.seed(4); tmp <- cfact_fore(mod, nsteps=3, nsim=2, pi=0.9, pred_type="median", cfact_type="muted_response", policy_var=2, cfact_start=1,
                                 cfact_end=1, mute_var=1)
  expect_equal(unname(c(tmp$cfact_pred[1:2,])), c(0.3849129, 0.8063971, 1.5854783, 1.4866707), tolerance=1e-3)
  expect_equal(tmp$cfact_pred_ints[c(1, 7, 9)], c(0.218715, 1.554466, 1.308472), tolerance=1e-3)

  # Exponential
  mod <- mod222expcmwtsh_2_1
  set.seed(5); tmp <- cfact_fore(mod, nsteps=3, nsim=3, pi=0.8, pred_type="mean", cfact_type="fixed_path", policy_var=2, cfact_start=1,
                                 cfact_end=1, cfact_path=c(-20))
  expect_equal(unname(c(tmp$cfact_pred[1:2,])), c(-788.7226, -169.2069, -20.0000, -54.1639), tolerance=1e-3)
  expect_equal(unname(c(tmp$cfact_trans_pred[1,])), c(0.92357051, 0.07642949), tolerance=1e-3)

  # Threshold
  mod <- mod222thres_2_1
  set.seed(3); tmp <- cfact_fore(mod, nsteps=3, nsim=1, pi=0.95, pred_type="mean", cfact_type="fixed_path", policy_var=2, cfact_start=1,
                                 cfact_end=3, cfact_path=c(-1, 2, -3))
  expect_equal(unname(c(tmp$cfact_pred[1:3,])), c(0.9936570, 1.2205136, 0.9502944, -1.0000000, 2.0000000, -3.0000000), tolerance=1e-3)
  expect_equal(tmp$cfact_pred_ints[c(1, 6, 9)], c(0.9936570, 0.9502944, -3.0000000), tolerance=1e-3)
  expect_equal(c(tmp$cfact_trans_pred_ints[c(2, 8)]), c(1, 0), tolerance=1e-3)

  # Exogenous
  mod <- mod123exoikt
  set.seed(6); tmp <- cfact_fore(mod, nsteps=3, nsim=2, pi=0.9, pred_type="mean", cfact_type="muted_response", policy_var=2, cfact_start=1,
                                 cfact_end=2, mute_var=3, exo_weights=cbind(c(0.8, 0, 0.3), c(0.2, 1, 0.7)))
  expect_equal(unname(c(tmp$cfact_pred[c(1, 2, 5, 7, 9)])), c(2.4966361, -7.3056026, 2.0205320, 0.5530254, -0.6760020), tolerance=1e-3)
  expect_equal(tmp$cfact_pred_ints[c(1, 7, 9, 16)], c(1.919356, 1.292601, 1.360463, 1.671384), tolerance=1e-3)
})


test_that("cfact_girf works correctly", {
  # Linear SVAR
  mod <- mod112relg
  mod$model$identification <- "recursive" # Reduced form model, flag to recursive to avoid messages
  tmp <- cfact_girf(mod, which_shocks=2, N=2, R1=2, R2=2, init_regime=1, scale=c(1, 1, 0.3), ci=0.9, seeds=1:2, use_parallel=FALSE,
                    cfact_type="muted_response", policy_var=1, mute_var=2, cfact_start=0, cfact_end=1)
  expect_equal(tmp$girf$girf_res$shock2$point_est[c(1, 2, 3, 4, 5)], c(0.000000, 0.000000, -0.067445, 0.522410, 0.468287), tolerance=1e-3)

  # Relative_dens Gaussian STVAR
  mod <- mod123relg
  mod$model$identification <- "recursive" # Reduced form model, flag to recursive to avoid messages
  tmp <- cfact_girf(mod, which_shocks=1, N=2, R1=2, R2=3, init_regime=1, scale=c(1, 1, 0.3), ci=0.9, seeds=1:3, use_parallel=FALSE,
                    cfact_type="muted_response", policy_var=3, mute_var=2, cfact_start=0, cfact_end=2)
  expect_equal(tmp$girf$girf_res$shock1$conf_ints[c(1, 6, 7, 14, 17)], c(0.30000000, 0.22600039, 0.01028607, 0.08353936, 0.08483720), tolerance=1e-3)

  # Logistic
  mod <- mod222logistitb
  tmp <- cfact_girf(mod, which_shocks=1, N=3, R1=3, R2=2, init_regime=1, scale=c(1, 2, 0.3), ci=0.9, seeds=4:5, use_parallel=FALSE,
                    cfact_type="fixed_path", policy_var=2, cfact_start=1, cfact_end=2, cfact_path=c(1, -1))
  expect_equal(tmp$girf$girf_res$shock1$point_est[c(1, 2, 5, 6, 7, 8)], c(0.939669, 0.161927, 0.300000, 0.000000, 0.000000, 0.025074), tolerance=1e-3)

  # Logit
  mod <- mod222logcmt_12_2
  mod$model$identification <- "recursive" # Reduced form model, flag to recursive to avoid messages
  tmp <- cfact_girf(mod, which_shocks=2, N=2, R1=3, R2=2, init_regime=1, scale=c(1, 2, 0.3), ci=0.9, seeds=1:2, use_parallel=FALSE,
                    cfact_type="fixed_path", policy_var=1, cfact_start=2, cfact_end=2, cfact_path=c(3))
  expect_equal(tmp$girf$girf_res$shock2$point_est[c(1, 2, 3, 5, 6, 8)], c(0.000000000, -0.052473993, 0.000000000, 0.137319453,
                                                                          0.156611490, 0.007277011), tolerance=1e-3)

  # Exponential
  mod <- mod222expcmwtsh_2_1
  tmp <- cfact_girf(mod, which_shocks=1, N=3, R1=3, R2=3, init_regime=1, scale=c(1, 1, 0.3), ci=0.9, seeds=1:3, use_parallel=FALSE,
                    cfact_type="muted_response", policy_var=2, mute_var=1, cfact_start=3, cfact_end=3)
  expect_equal(tmp$girf$girf_res$shock1$conf_ints[c(1, 6, 7, 14, 17)], c(0.3000000, 0.4284184, 0.1788517, -0.8757628, 0.0000000), tolerance=1e-3)

  # Threshold
  mod <- mod222thres_2_1
  mod$model$identification <- "recursive" # Reduced form model, flag to recursive to avoid messages
  tmp <- cfact_girf(mod, which_shocks=1, N=1, R1=1, R2=1, init_regime=1, scale=c(1, 1, 0.3), ci=0.9, seeds=1:3, use_parallel=FALSE,
                    cfact_type="muted_response", policy_var=2, mute_var=1, cfact_start=0, cfact_end=0)
  expect_equal(tmp$girf$girf_res$shock1$point_est[1:5], c(0.3000000, 0.0419880, 0.0000000, 0.0105516, 0.0000000), tolerance=1e-3)

  # Exogenous
  mod <- mod123exoikt
  tmp <- cfact_girf(mod, which_shocks=1:2, N=2, R1=3, R2=2, init_regime=1, scale=c(1, 2, 0.3), ci=0.9, seeds=1:2, use_parallel=FALSE,
                    exo_weights=cbind(c(0.8, 0, 0.3), c(0.2, 1, 0.7)), cfact_type="fixed_path", policy_var=3, cfact_start=0, cfact_end=0, cfact_path=c(3))
  expect_equal(tmp$girf$girf_res$shock1$point_est[c(1, 5, 7, 8, 9)], c(-0.856752, 0.155998, 0.000000, -0.252371, -0.503199), tolerance=1e-3)
  expect_equal(tmp$girf$girf_res$shock2$point_est[c(1, 5, 7, 8, 9)], c(0.844930, 0.701460, 0.000000, 0.238668, 0.007718), tolerance=1e-3)
})
