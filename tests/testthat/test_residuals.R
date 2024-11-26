context("residuals")
library(sstvars)

# p=1, M=1, d=2
theta_112relg <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103, 0.601786, -0.002945, 0.067224)

# p=1, M=2, d=2
theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863, -0.3124, 0.505251, 0.298483, 0.030096,
                   -0.176925, 0.838898, 0.310863, 0.007512, 0.018244, 0.949533, -0.016941, 0.121403, 0.573269)

# p=2, M=2, d=2
theta_222relg <- c(0.356914, 0.107436, 0.356386, 0.08633, 0.13996, 0.035172, -0.164575, 0.386816, 0.451675, 0.013086,
                   0.227882, 0.336084, 0.239257, 0.024173, -0.021209, 0.707502, 0.063322, 0.027287, 0.009182, 0.197066,
                   0.205831, 0.005157, 0.025877, 1.092094, -0.009327, 0.116449, 0.592446)

# p=3, M=2, d=3, usamone
theta_323relg <- c(0.98249, 0.66144, -1.17552, 0.50289, 0.17399, -0.01771, 0.96105, -0.11406, 0.41223,
                   -0.31217, 0.49067, 0.3958, 0.04185, 0.08454, 1.0977, -0.03208, 0.06398, -0.12298,
                   0.13382, 0.20166, 0.87613, -0.34591, -0.06254, -0.47386, -0.09049, 0.03109, 0.0347,
                   -0.16531, 0.0427, -0.31646, 0.25299, -0.04865, 0.33893, 0.69963, -0.02912, 0.03398,
                   -0.24344, 0.20815, 0.22566, 0.20582, 0.14774, 1.69008, 0.04375, -0.01018, -0.00947,
                   -0.19371, 0.26341, 0.22082, -0.08841, -0.18303, -0.86488, -0.06031, 0.00634, 0.00181,
                   -0.5559, 0.10249, -0.25146, -0.11875, 0.05153, 0.15267, 0.58151, -0.01903, 0.12236, 0.09327,
                   0.10245, 1.81845, 0.72719, 0.03235, 0.09857, 0.04826, 0.00908, 0.09761, 0.72127)

## weight_function == "mlogit"

# p=1, M=2, d=2, weightfun_pars=list(vars=1, lags=1)
gamma1_122_1_1 <- c(0.1, 0.2)
theta_122log_1_1 <- c(theta_122relg[-length(theta_122relg)], gamma1_122_1_1)

# p=2, M=2, d=2, weightfun_pars=list(vars=1:2, lags=2)
gamma1_222_12_2 <- c(0.1, 0.2, 0.11, 0.22, 0.33)
theta_222log_12_2 <- c(theta_222relg[-length(theta_222relg)], gamma1_222_12_2)


## weight_function == "threshold"

# p=2, M=2, d=2, weight_function="threshold", weightfun_pars=c(2, 1),
theta_222thres_2_1 <- c(theta_222relg[-length(theta_222relg)], 1)


## Constrained models
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}

# p=2, M=2, d=2, mean_constraints=list(1:2), AR_constraints=C_222, parametrization="mean"
C_222 <- rbind_diags(p=2, M=2, d=2)
theta_222relgcm <- c(0.7209658, 0.810858, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, 0.21, 0.01,
                     0.03, 1.1, 0.01, 0.11, 0.37)

# p=2, M=2, d=2, weigh_function="mlogit", weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), C_222
theta_222logcm_12_2 <- c(theta_222relgcm[-length(theta_222relgcm)], gamma1_222_12_2)

# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5), r=c(0, 0.11, 0.12, 0.13, 0)), parametrization="mean"
xi_222logcmw_12_2 <- c(0.002, 1.33)
theta_222logcmw_12_2 <-  c(theta_222relgcm[-length(theta_222relgcm)], xi_222logcmw_12_2)

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean"
xi_222logisticcmw_2_1 <- c(0.33)
theta_222logisticcmw_2_1 <- c(theta_222relgcm[-length(theta_222relgcm)], xi_222logisticcmw_2_1)

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean"
xi_222expcmw_2_1 <- c(0.33)
theta_222expcmw_2_1 <- c(theta_222relgcm[-length(theta_222relgcm)], xi_222logisticcmw_2_1)

## Structural models

# p=1, M=2, d=2, weight_function="relative_dens", identification="heteroskedasticity"
all_phi_122 <- c(0.734054, 0.225598, 0.705744, 0.187897)
all_A_122 <- c(0.259626, -0.000863, -0.3124, 0.505251, 0.298483, 0.030096, -0.176925, 0.838898)
W_122 <- matrix(c(-0.03, 0.24, -0.76, -0.02), nrow=2, ncol=2, byrow=FALSE)
lambdas_122 <- c(3.36, 0.86)
alpha1_122 <- 0.6
theta_122relgsh <- c(all_phi_122, all_A_122, vec(W_122), lambdas_122, alpha1_122)

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student",
# identification="heteroskedasticity"
all_phi_222 <- all_phi_122
all_A_222 <- c(0.13996, 0.035172, -0.164575, 0.386816, 0.451675, 0.013086, 0.227882, 0.336084, 0.239257, 0.024173,
               -0.021209, 0.707502, 0.063322, 0.027287, 0.009182, 0.197066)
W_222 <- W_122; lambdas_222 <- lambdas_122
c_and_gamma_222_2_1 <- c(0.1, 0.2)
df_222_2_1 <- 7
theta_222logistictsh_2_1 <- c(all_phi_222, all_A_222, vec(W_222), lambdas_222, c_and_gamma_222_2_1, df_222_2_1)


## Exogenous

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=weightfun_pars222
set.seed(1); weightfun_pars222 <- abs(matrix(rnorm(2*(nrow(gdpdef) - 2)), nrow=nrow(gdpdef) - 2, ncol=2))
weightfun_pars222 <- weightfun_pars222/rowSums(weightfun_pars222)
theta_222exo <- c(all_phi_222, all_A_222, 0.205831, 0.005157, 0.025877, 1.092094, -0.009327, 0.116449)

## ind_Student

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=weightfun_pars222, cond_dist="ind_Student"
B1_222 <- matrix(c(0.2, -0.3, 0.1, 0.4), nrow=2)
B2_222 <- matrix(c(-0.1, 0.2, 0.3, 0.4), nrow=2)
dfs_222 <- c(3, 7)
theta_222exoit <- c(all_phi_222, all_A_222, vec(B1_222), vec(B2_222), dfs_222)

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="ind_Student"
theta_222logistit <- c(all_phi_222, all_A_222, vec(B1_222), vec(B2_222), c_and_gamma_222_2_1, dfs_222)

## ind_skewed_t

# p=2, M=2, d=2, weight_function="exogenous", weightfun_pars=weightfun_pars222, cond_dist="ind_skewed_t"
B1_222 <- matrix(c(0.2, -0.3, 0.1, 0.4), nrow=2)
B2_222 <- matrix(c(-0.1, 0.2, 0.3, 0.4), nrow=2)
dfls_222 <- c(3, 7, 0.1, -0.2)
theta_222exoikt <- c(all_phi_222, all_A_222, vec(B1_222), vec(B2_222), dfls_222)

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="ind_skewed_t"
theta_222logistikt <- c(all_phi_222, all_A_222, vec(B1_222), vec(B2_222), c_and_gamma_222_2_1, dfls_222)


test_that("get_residuals works correctly", {
  # ind_Student
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222exoit, weight_function="exogenous",
                               weightfun_pars=weightfun_pars222, cond_dist="ind_Student", standardize=FALSE)[c(1, 4, 211, 242),]),
               c(-1.72187357, -1.77134204, -0.16040969, -0.70595177, -0.07764988, -0.35084729, 0.16804125, -0.16616452), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222exoit, weight_function="exogenous",
                               weightfun_pars=weightfun_pars222, cond_dist="ind_Student", standardize=TRUE)[c(1, 4, 211, 242),]),
               c(-9.1271938, -35.5238379, 2.2777582, -3.4403451, -4.0396570, -2.8173797, -0.2085455, -1.8803190), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222logistit, weight_function="logistic", weightfun_pars=c(2, 1),
                               cond_dist="ind_Student", standardize=FALSE)[c(1, 5, 212, 242),]),
               c(-1.58310517, -0.72319699, -0.65773933, -0.61862850, -0.07319559, -0.13840029, -0.21880463, -0.16385961), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222logistit, weight_function="logistic", weightfun_pars=c(2, 1),
                               cond_dist="ind_Student", standardize=TRUE)[c(1, 5, 212, 242),]),
               c(-21.009310, -9.210099, -8.864146, -7.759104, -2.739508, -1.408796,-1.374498, -1.270225), tolerance=1e-3)

  # Exogenous
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222exo, weight_function="exogenous",
                               weightfun_pars=weightfun_pars222, standardize=FALSE)[c(1, 2, 210, 242),]),
               c(-1.72187357, -0.72509735, 0.22855284, -0.70595177, -0.07764988, -0.14141620, -0.61638685, -0.16616452), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222exo, weight_function="exogenous",
                               weightfun_pars=weightfun_pars222, standardize=TRUE)[c(1, 2, 210, 242),]),
               c(-2.5984273, -0.7508414, 0.2692129, -1.0677359, -0.3305885, -0.4580447, -2.2484807, -0.7402457), tolerance=1e-3)
  # threshold
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222thres_2_1, weight_function="threshold",
                               weightfun_pars=c(2, 1), standardize=TRUE)[c(1, 2, 210, 242),]),
               c(-3.3073904, -2.3720290, 0.9687195, -0.9181623, 0.3946221, 0.4722489, -2.8858480, -0.2656000), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222thres_2_1, weight_function="threshold",
                               weightfun_pars=c(2, 1), standardize=FALSE)[c(1, 2, 210, 242),]),
               c(-1.49694675, -1.07200730, 0.41519088, -0.41871592, 0.03562626, 0.05594944, -0.45546140, -0.05037556), tolerance=1e-3)

  # exponential
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222expcmw_2_1, weight_function="exponential",
                               weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222, parametrization="mean",
                               weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), standardize=TRUE)[c(1, 2, 131, 242),]),
               c(-3.2931481, -1.8707419, 1.1592722, -0.6483687, -0.1403753, 0.1644838, -1.7407096, -0.7096609), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222expcmw_2_1, weight_function="exponential",
                               weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222, parametrization="mean",
                               weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), standardize=FALSE)[c(1, 2, 131, 242),]),
               c(-1.53180711, -0.93325451, 0.61629329, -0.32727751, -0.07600599, 0.00267501, -0.32524759, -0.13736299), tolerance=1e-3)

  # mlogit and logistic
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222logcmw_12_2, weight_function="mlogit",
                               weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222, parametrization="mean",
                               weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5), r=c(0, 0.11, 0.12, 0.13, 0)),
                               standardize=FALSE)[c(1, 2, 131, 242),]),
               c(-1.53180711, -0.93325451, 0.61629329, -0.32727751, -0.07600599, 0.00267501, -0.32524759, -0.13736299), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222logcmw_12_2, weight_function="mlogit", parametrization="mean",
                               weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222,
                               weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5), r=c(0, 0.11, 0.12, 0.13, 0)),
                               standardize=TRUE)[c(1, 2, 131, 242),]),
               c(-2.20635541, -1.27687738, 0.98627696, -0.47684704, -0.22411340, 0.06489152, -1.53428089, -0.58082203), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222logisticcmw_2_1, weight_function="logistic",
                               weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222, parametrization="mean",
                               weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), standardize=FALSE)[c(1, 2, 131, 242),]),
               c(-1.53180711, -0.93325451, 0.61629329, -0.32727751, -0.07600599, 0.00267501, -0.32524759, -0.13736299), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222logisticcmw_2_1, weight_function="logistic",
                               weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222, parametrization="mean",
                               weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), standardize=TRUE)[c(1, 2, 131, 242),]),
               c(-1.87533862, -1.13026252, 0.75108910, -0.39190481, -0.22041415, 0.04830197, -1.22258744, -0.49855698), tolerance=1e-3)


  # Relative_dens Gausssian STVAR
  expect_equal(c(get_residuals(data=gdpdef, p=1, M=1, params=theta_112relg, weight_function="relative_dens",
                               standardize=TRUE)[c(1, 2, 113, 243),]),
               c(1.3968841, -1.5548925, 0.3948874, -0.4266263, -0.6188517, 0.4927395, 0.1613921, 0.1354454), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=1, M=1, params=theta_112relg, weight_function="relative_dens",
                               standardize=FALSE)[c(1, 13, 200, 243),]),
               c(1.08538455, -0.18623807, -1.15119664, -0.33133796, -0.16441838, -0.40278905, -0.28125618, 0.03632955), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=1, M=2, params=theta_122relg, weight_function="relative_dens",
                               standardize=TRUE)[c(1, 3, 213, 243),]),
               c(1.5833935, -0.5878140, -0.5001501, -0.5921047, -1.0392301, -0.1613458, -0.8607716, 0.1433740), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=1, M=2, params=theta_122relg, weight_function="relative_dens",
                               standardize=FALSE)[c(1, 20, 242, 243),]),
               c(1.07689397, 1.40168826, -0.08105729, -0.34307860, -0.21260122, -0.44618006, -0.20047107, 0.01632420), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222relg, weight_function="relative_dens",
                               standardize=TRUE)[c(1, 2, 210, 242),]),
               c(-1.3062112, -0.3485900, 0.8453746, -0.8484317, 0.1119452, -0.1325865, -2.5016219, -0.2633884), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222relg, weight_function="relative_dens",
                               standardize=FALSE)[c(1, 54, 150, 242),]),
               c(-1.14633078, 0.48463608, -0.21489690, -0.40965919, 0.03689022, 0.43521031, 0.10575989, -0.05058704), tolerance=1e-3)

  expect_equal(c(get_residuals(data=gdpdef, p=1, M=2, params=theta_122log_1_1, weight_function="mlogit",
                               weightfun_pars=list(vars=1, lags=1), standardize=FALSE)[c(1, 22, 242, 243),]),
               c(1.06589778, 0.59845112, -0.11054445, -0.35996291, -0.22665544, 0.02447191, -0.26597485, -0.02033956), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=1, M=2, params=theta_122log_1_1, weight_function="mlogit",
                               weightfun_pars=list(vars=1, lags=1), standardize=TRUE)[c(2, 23, 242, 243),]),
               c(-1.66726467, -0.98636396, -0.14828201, -0.46799903, 0.21655693, -0.10698226, -1.05996051, -0.08650509), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222log_12_2, weight_function="mlogit",
                               weightfun_pars=list(vars=1:2, lags=2), standardize=FALSE)[c(1, 2, 142, 242),]),
               c(-1.35588849, -0.79135089, 0.23827929, -0.30835286, 0.03613477, 0.01874721, -0.02179699, -0.05295253), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222log_12_2, weight_function="mlogit",
                               weightfun_pars=list(vars=1:2, lags=2), standardize=TRUE)[c(1, 3, 143, 242),]),
               c(-2.06176907, 2.38657010, -0.58697015, -0.43064458, 0.17849521, -0.07701045, -0.15982137, -0.22108502), tolerance=1e-3)

  expect_equal(c(get_residuals(data=usamone, p=3, M=2, params=theta_323relg, weight_function="relative_dens",
                               standardize=TRUE)[c(1, 2, 120, 267),]),
               c(0.19617354, -0.10410342, 0.45792505, 0.06508034, -0.02659528, 1.28507226, 0.66641418, 0.20628137,
                 -0.71537249, 0.93914017, -0.32550858, -0.72564548), tolerance=1e-3)
  expect_equal(c(get_residuals(data=usamone, p=3, M=2, params=theta_323relg, weight_function="relative_dens",
                               standardize=FALSE)[c(1, 3, 141, 267),]),
               c(0.104542201, -0.342051868, -0.093341863, 0.002395015, -0.009236036, 0.457097339, 0.309139979, 0.014297527,
                 -0.201011975, 0.034911801, 0.195203026, -0.959691787), tolerance=1e-3)

  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222relgcm, weight_function="relative_dens", parametrization="mean",
                               AR_constraints=C_222, mean_constraints=list(1:2), standardize=TRUE)[c(1, 2, 111, 242),]),
               c(-1.45951694, -0.89019514, -0.06394894, -0.44363764, -0.19736637, 0.02752430, 0.27478373, -0.55008306), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222relgcm, weight_function="relative_dens", parametrization="mean",
                               AR_constraints=C_222, mean_constraints=list(1:2), standardize=FALSE)[c(1, 2, 113, 242),]),
               c(-1.53180711, -0.93325451, 0.05443799, -0.32727751, -0.07600599, 0.00267501, 0.02195891, -0.13736299), tolerance=1e-3)

  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222logcm_12_2, weight_function="mlogit", parametrization="mean",
                               weightfun_pars=list(vars=1:2, lags=2), AR_constraints=C_222, mean_constraints=list(1:2),
                               standardize=TRUE)[c(1, 2, 111, 242),]),
               c(-2.31162836, -1.31484188, -0.05949686, -0.44674751, -0.22244747, 0.06969202, 0.26055843, -0.55304014), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222logcm_12_2, weight_function="mlogit", parametrization="mean",
                               weightfun_pars=list(vars=1:2, lags=2), AR_constraints=C_222, mean_constraints=list(1:2),
                               standardize=FALSE)[c(1, 2, 131, 242),]),
               c(-1.53180711, -0.93325451, 0.61629329, -0.32727751, -0.07600599, 0.00267501, -0.32524759, -0.13736299), tolerance=1e-3)
  tmp_pars <- theta_222logcm_12_2
  tmp_pars[6:7] <- 3
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=tmp_pars, weight_function="mlogit", parametrization="mean",
                               weightfun_pars=list(vars=1:2, lags=2), AR_constraints=C_222, mean_constraints=list(1:2),
                               standardize=FALSE, allow_unstab=TRUE)[c(1, 2, 131, 242),]),
               c(-4.6361402, -4.9068434, 1.1747641, -0.5136389, 1.6408990, 1.1140078, 0.2422395, 1.1295519), tolerance=1e-3)

  # Student
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=c(theta_222thres_2_1, 100), weight_function="threshold",
                               weightfun_pars=c(2, 1), cond_dist="Student", standardize=TRUE)[c(1, 242),]),
               c(-3.3073904, -0.9181623, 0.3946221, -0.2656000), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=c(theta_222thres_2_1, 3), weight_function="threshold",
                               weightfun_pars=c(2, 1), cond_dist="Student", standardize=FALSE)[c(2, 242),]),
               c(-1.07200730, -0.41871592, 0.05594944, -0.05037556), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=c(theta_222expcmw_2_1, 20), weight_function="exponential", parametrization="mean",
                               weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2), AR_constraints=C_222,
                               weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), standardize=TRUE)[c(1, 242),]),
               c(-3.2931481, -0.6483687, -0.1403753, -0.7096609), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=c(theta_222expcmw_2_1, 3), weight_function="exponential", parametrization="mean",
                               weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2), AR_constraints=C_222,
                               weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), standardize=FALSE)[c(131, 242),]),
               c(0.6162933, -0.3272775, -0.3252476, -0.1373630), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=c(theta_222logcmw_12_2, 7), weight_function="mlogit", parametrization="mean",
                               weightfun_pars=list(vars=1:2, lags=2), cond_dist="Student", mean_constraints=list(1:2), AR_constraints=C_222,
                               weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5), r=c(0, 0.11, 0.12, 0.13, 0)),
                               standardize=FALSE)[c(1, 2),]),
               c(-1.53180711, -0.93325451, -0.07600599, 0.00267501), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=c(theta_222logcmw_12_2, 4), weight_function="mlogit", parametrization="mean",
                               weightfun_pars=list(vars=1:2, lags=2), cond_dist="Student", mean_constraints=list(1:2), AR_constraints=C_222,
                               weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5), r=c(0, 0.11, 0.12, 0.13, 0)),
                               standardize=TRUE)[c(1, 242),]),
               c(-2.2063554, -0.4768470, -0.2241134, -0.5808220), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=c(theta_222logisticcmw_2_1, 5), weight_function="logistic", parametrization="mean",
                               weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2), AR_constraints=C_222,
                               weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), standardize=FALSE)[c(2, 131),]),
               c(-0.93325451, 0.61629329, 0.00267501, -0.32524759), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=c(theta_222logisticcmw_2_1, 12), weight_function="logistic", parametrization="mean",
                               weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2), AR_constraints=C_222,
                               weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), standardize=TRUE)[c(1, 242),]),
               c(-1.8753386, -0.3919048, -0.2204141, -0.4985570), tolerance=1e-3)

  # Structural models
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=c(theta_222thres_2_1, 100), weight_function="threshold",
                               weightfun_pars=c(2, 1), cond_dist="Student", identification="recursive", standardize=TRUE)[c(1, 242),]),
               c(-3.3073904, -0.9181623, 0.3946221, -0.2656000), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=c(theta_222thres_2_1, 3), weight_function="threshold",
                               weightfun_pars=c(2, 1), cond_dist="Student", identification="recursive", standardize=FALSE)[c(2, 242),]),
               c(-1.07200730, -0.41871592, 0.05594944, -0.05037556), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=1, M=2, params=theta_122relgsh, weight_function="relative_dens",
                               identification="heteroskedasticity", standardize=TRUE)[c(1, 243),]),
               c(1.4454483, -0.4589756, -0.7624527, 0.0456093), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=1, M=2, params=theta_122relgsh, weight_function="relative_dens",
                               identification="heteroskedasticity", standardize=FALSE)[c(12, 243),]),
               c(0.628148740, -0.346101003, 0.117316575, 0.009761137), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222logistictsh_2_1, weight_function="logistic", weightfun_pars=c(2, 1),
                               cond_dist="Student", identification="heteroskedasticity", standardize=TRUE)[c(1, 200),]),
               c(-2.1569257, -0.1839126, -0.2149345, -0.9619645), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222logistictsh_2_1, weight_function="logistic", weightfun_pars=c(2, 1),
                               cond_dist="Student", identification="heteroskedasticity", standardize=FALSE)[c(1, 242),]),
               c(-1.58310517, -0.61862850, -0.07319559, -0.16385961), tolerance=1e-3)


  # Structural shocks
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=c(theta_222thres_2_1, 100), weight_function="threshold",
                               weightfun_pars=c(2, 1), cond_dist="Student", identification="recursive",
                               structural_shocks=TRUE)[c(1, 242),]),
               c(-3.2995215, -0.9229200, 0.4557586, -0.2485637), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=1, M=2, params=theta_122relgsh, weight_function="relative_dens",
                               identification="heteroskedasticity", structural_shocks=TRUE)[c(1, 243),]),
               c(-0.83735425, 0.06954961, -1.40338622, 0.45596239), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222logistictsh_2_1, weight_function="logistic", weightfun_pars=c(2, 1),
                               cond_dist="Student", identification="heteroskedasticity", structural_shocks=TRUE)[c(12, 142),]),
               c(-1.35036836, -0.42294593, 0.95211166, 0.06125506), tolerance=1e-3)

  # ind_Student
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222exoit, weight_function="exogenous",
                               weightfun_pars=weightfun_pars222, cond_dist="ind_Student", structural_shocks=TRUE)[c(1, 4, 211, 242),]),
               c(-9.1271938, -35.5238379, 2.2777582, -3.4403451, -4.0396570, -2.8173797, -0.2085455, -1.8803190), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222logistit, weight_function="logistic", weightfun_pars=c(2, 1),
                               cond_dist="ind_Student", structural_shocks=TRUE)[c(1, 5, 212, 242),]),
               c(-21.009310, -9.210099, -8.864146, -7.759104, -2.739508, -1.408796, -1.374498, -1.270225), tolerance=1e-3)

  # ind_skewed_t
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222exoikt, weight_function="exogenous",
                               weightfun_pars=weightfun_pars222, cond_dist="ind_skewed_t", structural_shocks=TRUE)[c(1, 4, 211, 242),]),
               c(-9.1271938, -35.5238379, 2.2777582, -3.4403451, -4.0396570, -2.8173797, -0.2085455, -1.8803190), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222logistikt, weight_function="logistic", weightfun_pars=c(2, 1),
                               cond_dist="ind_skewed_t", structural_shocks=TRUE)[c(1, 5, 212, 242),]),
               c(-21.009310, -9.210099, -8.864146, -7.759104, -2.739508, -1.408796, -1.374498, -1.270225), tolerance=1e-3)
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=theta_222logistikt, weight_function="logistic", weightfun_pars=c(2, 1),
                               cond_dist="ind_skewed_t", structural_shocks=FALSE)[c(1, 5, 212, 242),]),
               c(-21.009310, -9.210099, -8.864146, -7.759104, -2.739508, -1.408796, -1.374498, -1.270225), tolerance=1e-3)
  tmp_pars <- theta_222logistikt
  tmp_pars[7:10] <- 2
  expect_equal(c(get_residuals(data=gdpdef, p=2, M=2, params=tmp_pars, weight_function="logistic", weightfun_pars=c(2, 1),
                               cond_dist="ind_skewed_t", structural_shocks=FALSE, allow_unstab=TRUE)[c(1, 5, 212, 242),]),
               c(-29.573936, -20.124658, -19.706305, -13.917376, -8.830849, -8.668706, -6.193063, -4.519659), tolerance=1e-3)
})
