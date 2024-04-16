context("predict.stvar")
library(sstvars)

# p=1, M=1, d=2
theta_112relg <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103, 0.601786, -0.002945, 0.067224)
mod112relg <- STVAR(data=gdpdef, p=1, M=1, params=theta_112relg, weight_function="relative_dens")

# p=1, M=2, d=2
theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863, -0.3124, 0.505251, 0.298483, 0.030096,
                   -0.176925, 0.838898, 0.310863, 0.007512, 0.018244, 0.949533, -0.016941, 0.121403, 0.573269)
mod122relg <- STVAR(data=gdpdef, p=1, M=2, d=2, params=theta_122relg, weight_function="relative_dens")

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

# p=3, M=2, d=2, weight_function=logistic, weightfun_pars=c(2, 1)
theta_322logistic <- c(0.575516, 0.039404, 2.747048, 0.298252, 0.268959, 0.055646, -0.28457, 0.368543,
                       0.321873, 0.026186, -0.134235, 0.128367, -0.063395, 0.01267, -0.04303, 0.296143,
                       0.173863, -0.028923, -1.124229, 0.652773, -0.046647, 0.003917, 0.611047, 0.089456,
                       -0.066099, 0.045613, -0.651387, 0.066754, 0.356738, 0.005062, 0.031868, 1.15855,
                       -0.039808, 0.15342, 1.193073, 7.518277)
mod322logistic <- STVAR(data=gdpdef, p=3, M=2, params=theta_322logistic, weight_function="logistic", weightfun_pars=c(2, 1))

# p=3, M=2, d=2, weight_function=mlogit, weightfun_pars=list(vars=2, lags=1)
theta_322mlogit <- c(0.575457, 0.039417, 2.746742, 0.297945, 0.268969, 0.055635, -0.284341, 0.368566, 0.321876,
                     0.026181, -0.134244, 0.128348, -0.063425, 0.012676, -0.043061, 0.296153, 0.173878, -0.028861,
                     -1.123896, 0.652866, -0.046737, 0.003973, 0.61057, 0.089591, -0.066095, 0.045594, -0.651087,
                     0.066679, 0.356729, 0.005054, 0.031869, 1.158638, -0.039743, 0.153416, 8.970726, -7.518924)
mod322mlogit <- STVAR(data=gdpdef, p=3, M=2, params=theta_322mlogit, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1))

# p=3, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1)
theta_322exp <- c(0.39211, 0.01111, 2.36699, 0.46926, 0.25534, 0.03812, -0.03716, 0.37795, 0.35103, 0.02618, -0.09847,
                  0.16067, -0.03619, 0.02972, -0.08469, 0.28501, 0.20048, -0.01467, -0.96397, 0.66189, -0.10317, 0.00512,
                  0.72416, 0.05694, -0.06667, 0.02475, -0.70513, 0.01575, 0.34356, 0.00828, 0.02348, 1.29206, -0.05969,
                  0.19413, 0.48621, 1.01986)
mod322exp <- STVAR(data=gdpdef, p=3, M=2, params=theta_322exp, weight_function="exponential", weightfun_pars=c(2, 1))

# p=3, M=2, d=2, weight_function="threshold", weightfun_pars=c(2, 1)
theta_322thres <- c(0.50787, 0.04477, 1.9715, 0.11992, 0.27899, 0.0507, -0.0855, 0.4122, 0.29264, 0.02062, -0.12133, 0.11042,
                    -0.06851, 0.01426, -0.02835, 0.2834, 0.18424, -0.01164, -0.88906, 0.65569, -0.00658, 0.01358, 0.60105,
                    0.11503, -0.03238, 0.05692, -0.52757, 0.11841, 0.37367, 0.00339, 0.03619, 1.13522, -0.03236, 0.13564, 1.10785)
mod322thres <- STVAR(data=gdpdef, p=3, M=2, params=theta_322thres, weight_function="threshold", weightfun_pars=c(2, 1))

## Student
mod322logistict <- STVAR(data=gdpdef, p=3, M=2, params=c(theta_322logistic, 7), weight_function="logistic", weightfun_pars=c(2, 1),
                        cond_dist="Student")
mod322mlogitt <- STVAR(data=gdpdef, p=3, M=2, params=c(theta_322mlogit, 3), weight_function="mlogit", weightfun_pars=list(vars=2, lags=1),
                      cond_dist="Student")
mod322expt <- STVAR(data=gdpdef, p=3, M=2, params=c(theta_322exp, 20), weight_function="exponential", weightfun_pars=c(2, 1),
                    cond_dist="Student")
mod322threst <- STVAR(data=gdpdef, p=3, M=2, params=c(theta_322thres, 2.01), weight_function="threshold", weightfun_pars=c(2, 1),
                      cond_dist="Student")

## Structural

# p=3, M=2, d=2, weight_function="threshold", weightfun_pars=list(vars=2, lags=1), identification="recursive"
mod322thresr_2_1 <- STVAR(data=gdpdef, p=3, M=2, d=2, params=c(theta_322thres, 2.01), weight_function="threshold",
                          weightfun_pars=c(2, 1), cond_dist="Student", identification="recursive")

# p=1, M=2, d=2, weight_function="relative_dens", identification="heteroskedasticity"
all_phi_122 <- c(0.734054, 0.225598, 0.705744, 0.187897)
all_A_122 <- c(0.259626, -0.000863, -0.3124, 0.505251, 0.298483, 0.030096, -0.176925, 0.838898)
W_122 <- matrix(c(-0.03, 0.24, -0.76, -0.02), nrow=2, ncol=2, byrow=FALSE)
lambdas_122 <- c(3.36, 0.86)
alpha1_122 <- 0.6
theta_122relgsh <- c(all_phi_122, all_A_122, vec(W_122), lambdas_122, alpha1_122)
mod122relgsh <- STVAR(data=gdpdef, p=1, M=2, d=2, params=theta_122relgsh, weight_function="relative_dens", identification="heteroskedasticity")

## ind_Student

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
params322thresit <- c(0.53077, 0.003653, 1.299719, -0.004862, 0.276631, 0.045029, -0.029358, 0.427229, 0.296498, 0.027308,
                     -0.098684, 0.160642, -0.121114, 0.021843, -0.06632, 0.283899, 0.186097, -0.017753, -0.361943, 0.782694,
                     0.155553, 0.023192, 0.153155, 0.051293, -0.068252, 0.02092, -0.310782, 0.115128, 0.626432, 0.037041,
                     0.081263, -0.186364, -0.878845, -0.141764, 0.479942, -0.314881, 1, 4.7304, 9.933647)

mod322thresit <- STVAR(gdpdef, p=3, M=2, params=params322thresit, weight_function="threshold", weightfun_pars=c(2, 1),
                      cond_dist="ind_Student", identification="non-Gaussianity")


set.seed(1); p112 <- predict(mod112relg, nsteps=1, nsim=1, pred_type="mean")
set.seed(2); p122 <- predict(mod122relg, nsteps=3, nsim=3, pred_type="median")
set.seed(3); p222 <- predict(mod222relg, nsteps=2, nsim=10, pred_type="mean", pi=c(0.5))
set.seed(4); p123 <- predict(mod123relg, nsteps=4, nsim=7, pred_type="median", pi=c(0.3, 0.5, 0.7))
set.seed(5); p322logistic <- predict(mod322logistic, nsteps=3, nsim=5, pred_type="mean", pi=0.9)
set.seed(5); p322mlogit <- predict(mod322mlogit, nsteps=3, nsim=5, pred_type="median", pi=c(0.7, 0.8))
set.seed(5); p322exp <- predict(mod322exp, nsteps=3, nsim=5, pred_type="mean", pi=0.95)
set.seed(5); p322thres <- predict(mod322thres, nsteps=3, nsim=5, pred_type="median", pi=c(0.6, 0.81))

set.seed(5); p322logistict <- predict(mod322logistict, nsteps=3, nsim=5, pred_type="mean", pi=0.9)
set.seed(5); p322mlogitt <- predict(mod322mlogitt, nsteps=3, nsim=5, pred_type="median", pi=c(0.7, 0.8))
set.seed(5); p322expt <- predict(mod322expt, nsteps=3, nsim=5, pred_type="mean", pi=0.95)
set.seed(5); p322threst <- predict(mod322threst, nsteps=3, nsim=5, pred_type="median", pi=c(0.6, 0.81))

set.seed(6); p322thresr <- predict(mod322thres, nsteps=3, nsim=5, pred_type="mean", pi=c(0.7, 0.99))
set.seed(6); p122relgsh <- predict(mod122relgsh, nsteps=4, nsim=3, pred_type="median", pi=c(0.7, 0.98))

set.seed(7); p122exocit <- predict(mod122exocit, nsteps=3, nsim=5, pred_type="mean", pi=0.9,
                                   exo_weights=cbind(c(0.1, 0.6, 0.9), c(0.9, 0.4, 0.1)))
set.seed(8); p322thresit <- predict(mod322thresit, nsteps=2, nsim=3, pred_type="median", pi=0.49)

test_that("simulate.stvar works correctly", {
  # Relative_dens Gaussian STVAR
  expect_equal(unname(p112$pred[1,]), c(0.2393741, 0.4936729), tol=1e-4)
  expect_equal(unname(p112$pred_ints[1, 3,]), c(0.2393741, 0.4936729), tol=1e-4)
  expect_equal(unname(p112$trans_pred[1,]), c(1), tol=1e-4)
  expect_equal(unname(p112$trans_pred_ints[1, 1, 1]), c(1), tol=1e-4)
  expect_equal(unname(p122$pred[3,]), c(1.4321082, 0.4083084), tol=1e-4)
  expect_equal(unname(p122$pred_ints[3, 2,]), c(1.0924034, 0.2907614), tol=1e-4)
  expect_equal(unname(p222$trans_pred[2,]), c(0.95938141, 0.04061859), tol=1e-4)
  expect_equal(unname(p222$trans_pred_ints[2, 1,]), c(0.95597200, 0.03039205), tol=1e-4)
  expect_equal(unname(p123$pred[4,]), c(2.154655, 1.960712, 1.612034), tol=1e-4)
  expect_equal(unname(p123$pred_ints[4, 5,]), c(2.422710, 2.372279, 7.711421), tol=1e-4)
  expect_equal(unname(p123$trans_pred[4,]), c(0.0002149667, 0.9997850333), tol=1e-4)
  expect_equal(unname(p123$trans_pred_ints[3, 4,]), c(0.1428647, 1.0000000), tol=1e-4)

  expect_equal(unname(p322logistic$pred[3,]), c(1.042285, 0.327765), tol=1e-4)
  expect_equal(unname(p322logistic$pred_ints[3, 2,]), c(1.7293984, 0.6218133), tol=1e-4)
  expect_equal(unname(p322mlogit$pred[3,]), c(1.3200944, 0.3322151), tol=1e-4)
  expect_equal(unname(p322mlogit$pred_ints[3, 4,]), c(1.5454664, 0.5335837), tol=1e-4)

  expect_equal(unname(p322exp$pred[3,]), c(0.9682415, 0.3139115), tol=1e-4)
  expect_equal(unname(p322exp$pred_ints[3, 2,]), c(1.6870691, 0.5851358), tol=1e-4)
  expect_equal(unname(p322thres$pred[3,]), c(1.364332, 0.344308), tol=1e-4)
  expect_equal(unname(p322thres$pred_ints[3, 4,]), c(1.4781161, 0.5018161), tol=1e-4)

  # Student
  expect_equal(unname(p322logistict$pred[3,]), c(0.8981389, 0.5252586), tol=1e-4)
  expect_equal(unname(p322logistict$pred_ints[3, 2,]), c(1.6169003, 0.6830229), tol=1e-4)
  expect_equal(unname(p322mlogitt$pred[3,]), c(0.9551891, 0.5039055), tol=1e-4)
  expect_equal(unname(p322mlogitt$pred_ints[3, 4,]), c(1.409086, 0.624657), tol=1e-4)

  expect_equal(unname(p322expt$pred[3,]), c(1.0935312, 0.5088775), tol=1e-4)
  expect_equal(unname(p322expt$pred_ints[3, 2,]), c(1.7067186, 0.6818502), tol=1e-4)
  expect_equal(unname(p322threst$pred[3,]), c(0.7924443, 0.4580927), tol=1e-4)
  expect_equal(unname(p322threst$pred_ints[3, 4,]), c(0.8058259, 0.4715564), tol=1e-4)

  # Structural
  expect_equal(unname(p322thresr$pred[3,]), c(1.6001306, 0.3472043), tol=1e-4)
  expect_equal(unname(p322thresr$pred_ints[3, 3,]), c(1.965468, 0.650851), tol=1e-4)
  expect_equal(unname(p122relgsh$pred[4,]), c(0.8342902, 0.8843367), tol=1e-4)
  expect_equal(unname(p122relgsh$pred_ints[4, 1,]), c(0.3007967, 0.3969968), tol=1e-4)

  # ind_Student
  expect_equal(unname(p122exocit$pred[3,]), c(0.3939130, 0.4663528), tol=1e-4)
  expect_equal(c(unname(p122exocit$pred_ints[3, ,])), c(-0.3053796, 1.5096820, 0.3453159, 0.6648726), tol=1e-4)
  expect_equal(c(unname(p322thresit$pred[2,])), c(0.0383185, 0.4116402), tol=1e-4)
  expect_equal(c(unname(p122exocit$pred_ints[2, ,])), c(-0.6105301, 1.8325592, 0.3266453, 0.7264805), tol=1e-4)
})

