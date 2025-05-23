context("simulate.stvar")
library(sstvars)

# NOTE: predict uses simulate, so some of the setups are only tested in predict, which breaks if simulate breaks

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

# mlogit weights, gdpdef, weightfun_pars=list(vars=2, lags=1)
theta_322log_2_1 <- c(2.746765, 0.297951, 0.57546, 0.039418, 0.173881, -0.028861, -1.123912, 0.652867, -0.046741,
                      0.003972, 0.610594, 0.089587, -0.066095, 0.045594, -0.651105, 0.066679, 0.268968, 0.055636,
                      -0.284343, 0.368566, 0.321876, 0.026181, -0.134247, 0.128348, -0.063424, 0.012676, -0.043061,
                      0.296153, 1.158639, -0.039743, 0.153417, 0.356729, 0.005054, 0.031869, -8.970684, 7.518877)
mod322log <- STVAR(data=gdpdef, p=3, M=2, d=2, params=theta_322log_2_1, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1))

mod322logt <- STVAR(data=gdpdef, p=3, M=2, d=2, params=c(theta_322log_2_1, 10), weight_function="mlogit", weightfun_pars=list(vars=2, lags=1),
                    cond_dist="Student")

# p=3, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1), identification="recursive"
mod322logtr_2_1 <- STVAR(data=gdpdef, p=3, M=2, d=2, params=c(theta_322log_2_1, 10), weight_function="mlogit",
                          weightfun_pars=list(vars=2, lags=1), cond_dist="Student", identification="recursive")

# p=1, M=2, d=2, weight_function="relative_dens", identification="heteroskedasticity"
all_phi_122 <- c(0.734054, 0.225598, 0.705744, 0.187897)
all_A_122 <- c(0.259626, -0.000863, -0.3124, 0.505251, 0.298483, 0.030096, -0.176925, 0.838898)
W_122 <- matrix(c(-0.03, 0.24, -0.76, -0.02), nrow=2, ncol=2, byrow=FALSE)
lambdas_122 <- c(3.36, 0.86)
alpha1_122 <- 0.6
theta_122relgsh <- c(all_phi_122, all_A_122, vec(W_122), lambdas_122, alpha1_122)
mod122relgsh <- STVAR(data=gdpdef, p=1, M=2, d=2, params=theta_122relgsh, weight_function="relative_dens", identification="heteroskedasticity")


# p=1, M=2, d=2, weight_function="exogenous", cond_dist="ind_Student", weightfun_pars=twmat, AR_constraints=C_122,
tw1 <- c(0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1,
 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
twmat <- cbind(tw1, 1 - tw1)
C_122 <- rbind(diag(1*2^2), diag(1*2^2))
params122exocit <- c(-0.25329555, 0.10340501, 0.73585978, 0.05335907, 0.18282378, 0.03478017, -0.01629061, 0.87424509, 0.83939713,
                     0.20970725, 0.47256152, -0.23555538, 0.59277160, 0.09901218, -0.29605914, 0.23895480, 6.72161948, 4.10403450)

mod122exocit <- STVAR(data=gdpdef, p=1, M=2, d=2, params=params122exocit, weight_function="exogenous", cond_dist="ind_Student",
                      weightfun_pars=twmat, AR_constraints=C_122)

# p=3, M=2, weight_function="threshold", cond_dist="ind_Student", weightfun_pars=c(2, 1), identification="non-Gaussianity"
params32thresit <- c(0.53077, 0.003653, 1.299719, -0.004862, 0.276631, 0.045029, -0.029358, 0.427229, 0.296498, 0.027308,
                     -0.098684, 0.160642, -0.121114, 0.021843, -0.06632, 0.283899, 0.186097, -0.017753, -0.361943, 0.782694,
                     0.155553, 0.023192, 0.153155, 0.051293, -0.068252, 0.02092, -0.310782, 0.115128, 0.626432, 0.037041,
                     0.081263, -0.186364, -0.878845, -0.141764, 0.479942, -0.314881, 1, 4.7304, 9.933647)

mod32thresit <- STVAR(gdpdef, p=3, M=2, params=params32thresit, weight_function="threshold", weightfun_pars=c(2, 1),
                      cond_dist="ind_Student", identification="non-Gaussianity")

# p=1, M=2, d=2, weight_function="exogenous", cond_dist="ind_skewed_t", weightfun_pars=twmat, AR_constraints=C_122
params122exocikt <- c(params122exocit, 0.1, -0.11)
mod122exocikt <- STVAR(data=gdpdef, p=1, M=2, d=2, params=params122exocikt, weight_function="exogenous", cond_dist="ind_skewed_t",
                      weightfun_pars=twmat, AR_constraints=C_122)

# p=3, M=2, weight_function="threshold", cond_dist="ind_skewed_t", weightfun_pars=c(2, 1), identification="non-Gaussianity"
params32thresikt <- c(params32thresit, -0.11, 0.1)
mod32thresikt <- STVAR(gdpdef, p=3, M=2, params=params32thresikt, weight_function="threshold", weightfun_pars=c(2, 1),
                      cond_dist="ind_skewed_t", identification="non-Gaussianity")

test_that("simulate.stvar works correctly", {
  s112 <- simulate(mod112relg, nsim=1, seed=1, init_regime=1, use_stat_for_Gaus=TRUE)
  s112_2 <- simulate(mod112relg, nsim=2, seed=1, init_regime=1)
  s122 <- simulate(mod122relg, nsim=5, seed=2, init_values=gdpdef, use_stat_for_Gaus=TRUE)
  s222 <- simulate(mod222relg, nsim=3, seed=3, init_regime=2, use_stat_for_Gaus=TRUE)
  s123 <- simulate(mod123relg, nsim=3, seed=4, init_regime=1, use_stat_for_Gaus=TRUE)
  s123_2 <- simulate(mod123relg, nsim=1, seed=5, init_values=usamone, use_stat_for_Gaus=TRUE)
  s322 <- simulate(mod322log, nsim=3, seed=3, init_regime=1, use_stat_for_Gaus=TRUE)
  s322_2 <- simulate(mod322log, nsim=3, seed=3, init_values=gdpdef, use_stat_for_Gaus=TRUE)
  s322t <- simulate(mod322logt, nsim=3, seed=3, init_values=gdpdef)
  s322t_2 <- simulate(mod322logt, nsim=3, seed=3, init_regime=2)
  s322tr <- simulate(mod322logtr_2_1, nsim=4, seed=3, init_values=gdpdef)
  s122relgsh <- simulate(mod122relgsh, nsim=3, seed=4, init_regime=2, use_stat_for_Gaus=TRUE)
  s122exocit <- simulate(mod122exocit, nsim=3, seed=5, init_regime=2, exo_weights=cbind(c(0.9, 0.5, 0.2), c(0.1, 0.5, 0.8)))
  s322thresit <- simulate(mod32thresit, nsim=4, seed=6, init_regime=1)
  s122exocikt <- simulate(mod122exocikt, nsim=3, seed=5, init_regime=2, exo_weights=cbind(c(0.9, 0.5, 0.2), c(0.1, 0.5, 0.8)))
  s322thresikt <- simulate(mod32thresikt, nsim=4, seed=6, init_regime=1)

  # Relative_dens Gaussian STVAR
  expect_equal(s112$sample[1,], c(-0.07206511, 1.343205), tol=1e-4)
  expect_equal(s112$transition_weights, as.matrix(1), tol=1e-4)
  expect_equal(c(s112_2$sample[1:2,]), c(0.1300491, 0.2689557, -0.1769628, 0.2490280), tol=1e-4)
  expect_equal(s112_2$transition_weights, as.matrix(c(1, 1)), tol=1e-4)
  expect_equal(s122$sample[5,], c(2.1049174, 0.4326244), tol=1e-4)
  expect_equal(s122$transition_weights[5,], c(0.94892269, 0.05107731), tol=1e-4)
  expect_equal(s222$sample[3,], c(0.1541768, 0.8947292), tol=1e-4)
  expect_equal(s222$transition_weights[1,], c(0.07983494, 0.92016506), tol=1e-4)
  expect_equal(s123$sample[3,], c(1.346929, 0.875446, 5.316535), tol=1e-4)
  expect_equal(s123$transition_weights[3,], c(0.98754515, 0.01245485), tol=1e-4)
  expect_equal(s123_2$sample[1,], c(1.684557, 2.126412, -2.721774), tol=1e-4)
  expect_equal(s123_2$transition_weights[1,], c(0.0001666309, 0.9998333691), tol=1e-4)

  # Logit
  expect_equal(s322$sample[3,], c(-0.6565136, 1.7113249), tol=1e-4)
  expect_equal(s322$transition_weights[1,], c(0.8947074, 0.1052926), tol=1e-4)
  expect_equal(s322_2$sample[3,], c(0.7678891, 0.3236954), tol=1e-4)
  expect_equal(s322_2$transition_weights[2,], c(0.002522048, 0.997477952), tol=1e-4)

  # Student
  expect_equal(s322t$sample[3,], c(0.4749844, 0.6297811), tol=1e-4)
  expect_equal(s322t$transition_weights[1,], c(0.002759205, 0.997240795), tol=1e-4)
  expect_equal(s322t_2$sample[3,], c(1.30571561, 0.09767046), tol=1e-4)
  expect_equal(s322t_2$transition_weights[1,], c(0.0008428691, 0.9991571309), tol=1e-4)

  # ind_Student
  expect_equal(s122exocit$sample[3,], c(-0.02252565, 0.15108849), tol=1e-4)
  expect_equal(s122exocit$transition_weights, cbind(c(0.9, 0.5, 0.2), c(0.1, 0.5, 0.8)), tol=1e-4)
  expect_equal(s322thresit$sample[4,], c(1.5028673, 0.3365102), tol=1e-4)
  expect_equal(s322thresit$transition_weights[4,], c(1, 0), tol=1e-4)

  # ind_skewed_t
  expect_equal(s122exocikt$sample[3,], c(0.8788235, 0.7216581), tol=1e-4)
  expect_equal(s122exocikt$transition_weights, cbind(c(0.9, 0.5, 0.2), c(0.1, 0.5, 0.8)), tol=1e-4)
  expect_equal(s322thresikt$sample[4,], c(1.07865055, -0.05331638), tol=1e-4)
  expect_equal(s322thresikt$transition_weights[4,], c(1, 0), tol=1e-4)

  # Structural
  expect_equal(s322tr$sample[4,], c(-0.05019288, -0.07683787), tol=1e-4)
  expect_equal(s322tr$transition_weights[4,], c(0.01415435, 0.985845655), tol=1e-4)
  expect_equal(s122relgsh$sample[3,], c(0.64751483, 0.9216674), tol=1e-4)
  expect_equal(s122relgsh$transition_weights[1,], c(0.6535839, 0.3464161), tol=1e-4)
})



test_that("simulate_from_regime works correctly", {
  set.seed(1); sim_reg112_1 <- simulate_from_regime(mod112relg, regime=1, nsim=1, use_transweights=FALSE)
  set.seed(1); sim_reg112_2 <- simulate_from_regime(mod112relg, regime=1, nsim=1, init_values=gdpdef, use_transweights=FALSE)
  set.seed(2); sim_reg122_1 <- simulate_from_regime(mod122relg, regime=1, nsim=2, use_transweights=FALSE)
  set.seed(3); sim_reg122_2 <- simulate_from_regime(mod122relg, regime=2, nsim=3, use_transweights=FALSE)
  set.seed(4); sim_reg322_1 <- simulate_from_regime(mod322log, regime=1, nsim=4, use_transweights=FALSE)
  set.seed(5); sim_reg322t_1 <- simulate_from_regime(mod322logt, regime=2, nsim=5,  use_transweights=TRUE)
  set.seed(5); sim_regs322tr_1 <- simulate_from_regime(mod322logtr_2_1, regime=2, nsim=1, use_transweights=FALSE)
  set.seed(6); sim_regs122relgsh_1 <- simulate_from_regime(mod122relgsh, regime=1, nsim=6, use_transweights=TRUE)
  set.seed(1); sim_reg122exocit_1 <- simulate_from_regime(mod122exocit, regime=2, nsim=2, use_transweights=FALSE)
  set.seed(2); sim_reg322thresit_1 <- simulate_from_regime(mod32thresit, regime=1, nsim=5, use_transweights=TRUE)

  # Relative_dens Gaussian STVAR
  expect_equal(sim_reg112_1[1,], c(0.2634843, 0.8543901), tol=1e-4)
  expect_equal(sim_reg112_2[1,], c(0.2393741, 0.4936729), tol=1e-4)
  expect_equal(sim_reg122_1[2,], c(1.5384970, 0.3277234), tol=1e-4)
  expect_equal(sim_reg122_2[3,], c(0.9612695, 0.8755360), tol=1e-4)

  # Logit
  expect_equal(sim_reg322_1[4,], c(-0.3290652, 1.7006393), tol=1e-4)

  # Student
  expect_equal(sim_reg322t_1[3,], c(0.5666465, 0.2113145), tol=1e-4)

  # ind_Student
  expect_equal(sim_reg122exocit_1[2,], c(2.2260794, 0.5797427), tol=1e-4)
  expect_equal(sim_reg322thresit_1[3,], c(1.4671417, 0.9364542), tol=1e-4)

  # Structural
  expect_equal(sim_regs322tr_1[1,], c(0.05300774, 0.8396193), tol=1e-4)
  expect_equal(sim_regs122relgsh_1[1,], c(-1.0609420, 0.1987437), tol=1e-4)
})
