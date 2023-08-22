context("uncondMoments")
library(sstvars)


alt_pcovmat <- function(p, d, all_A, all_Omegas) {
  # Calculate the (dp x dp) covariance matrix Sigma_{m,p} (Lutkepohl 2005, eq. (2.1.39))
  all_boldA <- form_boldA(p=p, M=1, d=d, all_A=all_A)
  I_dp2 <- diag(nrow=(d*p)^2)
  ZER_lower <- matrix(0, nrow=d*(p - 1), ncol=d*p)
  ZER_right <- matrix(0, nrow=d, ncol=d*(p - 1))
  kronmat <- I_dp2 - kronecker(all_boldA[, , 1], all_boldA[, , 1])
  sigma_epsm <- rbind(cbind(all_Omegas[, , 1], ZER_right), ZER_lower)
  Sigma_m <- solve(kronmat, vec(sigma_epsm))
  matrix(Sigma_m, nrow=d*p, ncol=d*p)
}

# pMd
params112 <- c(-0.114, 4.704, 0.23, -0.522, 0.014, 0.217, 0.04, -0.047, 1.649)
all_A112 <- pick_allA(p=1, M=1, d=2, params=params112)
all_Omegas112 <- pick_Omegas(p=1, M=1, d=2, params=params112)
all_boldA112 <- form_boldA(p=1, M=1, d=2, all_A=all_A112)

params115 <- c(0.936, 3.738, 4.135, 9.663, 3.955, -0.579, -0.009, -0.187, 0.019, -0.241, -0.02, -0.128, -0.036,
               -0.167, -0.218, 0.042, -0.259, -0.029, -0.03, -0.133, 0.377, -0.321, 0.017, -0.205, -0.29, -0.01,
               0.483, 0.286, -0.035, -0.216, 2.373, -0.613, 1.99, -2.255, 0.08, 1.14, 0.319, 0.342, -0.723, 2.832,
               -1.437, -0.769, 4.711, 2.736, 6.54)
all_A115 <- pick_allA(p=1, M=1, d=5, params=params115)
all_Omegas115 <- pick_Omegas(p=1, M=1, d=5, params=params115)

params212 <- c(1.551, 4.498, 0.146, -0.079, -0.045, -0.487, 0.271, -0.087, -0.697, -0.221, 0.25, 0.512, 2.507)
all_A212 <- pick_allA(p=2, M=1, d=2, params=params212)
all_Omegas212 <- pick_Omegas(p=2, M=1, d=2, params=params212)

params213 <- c(1.076, 1.777, -1.985, -0.141, 0.234, 0.088, -0.311, 0.112, 0.234, -0.421, -0.217, 0.017, -0.067,
               -0.104, 0.067, 0.075, 0.454, -0.249, 0.034, -0.247, -0.231, 0.499, -0.218, -1.044, 0.422, 0.986, 3.175)
all_A213 <- pick_allA(p=2, M=1, d=3, params=params213)
all_Omegas213 <- pick_Omegas(p=2, M=1, d=3, params=params213)
all_boldA213 <- form_boldA(p=2, M=1, d=3, all_A=all_A213)

params312 <- c(0.712, 4.819, -0.147, 0.365, -0.159, -0.055, 0.52, -0.257, -0.057, -0.313, 0.373, 0.182, -0.109,
               -0.111, 1.326, -0.895, 0.66)
all_A312 <- pick_allA(p=3, M=1, d=2, params=params312)
all_Omegas312 <- pick_Omegas(p=3, M=1, d=2, params=params312)


params714 <- c(2.211, 3.498, 5.362, 4.747, -0.114, -0.073, -0.138, 0.186, -0.062, -0.015, -0.102, -0.276, 0.203,
               0.084, 0.087, 0.049, 0.257, 0.262, 0.037, 0.087, 0.1, -0.243, -0.109, -0.275, 0.199, -0.125, -0.081,
               -0.129, 0.323, 0.426, -0.021, 0.323, 0.062, 0.057, 0.076, -0.045, -0.303, 0.093, 0.072, -0.209, 0.219,
               0.068, -0.131, -0.381, -0.03, -0.101, -0.083, 0.26, -0.116, 0.011, -0.043, 0.099, 0.114, 0.152, -0.209,
               0.003, 0.258, 0.047, 0.457, -0.034, -0.291, 0.325, 0.033, -0.004, -0.017, 0.244, 0.027, 0.054, 0.095,
               0.079, 0.258, -0.101, 0.09, 0.014, -0.123, -0.079, 0.178, 0.27, 0.173, 0.161, -0.075, -0.118, -0.245,
               0.08, -0.172, 0.261, -0.074, -0.046, 0.137, -0.081, -0.225, -0.055, 0.136, -0.366, 0.04, 0.162, -0.137,
               0.292, 0.087, 0.074, -0.026, -0.137, -0.011, 0.068, 0.083, 0.11, 0.139, 0.091, -0.297, -0.116, 0.166,
               -0.064, 0.132, -0.062, 0.251, -0.042, 0.681, -0.328, 0.944, 0.585, 2.154, 0.684, -2.075, 4.601, -1.727, 3.392)
all_A714 <- pick_allA(p=7, M=1, d=4, params=params714)
all_Omegas714 <- pick_Omegas(p=7, M=1, d=4, params=params714)
all_boldA714 <- form_boldA(p=7, M=1, d=4, all_A=all_A714)

params516 <- c(0.225461, 0.2072, -0.143174, -0.089739, -0.018469, -0.050561, -0.170449, 0.056068, 0.008218, 0.104494,
               -0.054231, -0.065922, -0.101427, 0.303574, 0.09502, 0.173378, 0.154364, -0.050695, 0.125909, 0.423252,
               -0.090888, 0.008861, 0.033686, 0.229432, -0.204607, 0.12374, -0.246613, -0.272826, -0.079953, -0.155406,
               0.147952, -0.045087, 0.135298, 0.024628, 0.112773, -0.087763, -0.040376, -0.072829, 0.006114, -0.006132,
               -0.081802, -0.099939, -0.026151, -0.080008, -0.066702, 0.160827, 0.10961, 0.032163, 0.278354, -0.151524,
               -0.159157, 0.069496, 0.2282, 0.141799, -0.042404, -0.035857, 0.001668, 0.051633, 0.177808, -0.108081,
               -0.003776, 0.086418, -0.177295, -0.007459, 0.001096, -0.08415, 0.077412, -0.153839, -0.222861, -0.221952,
               0.029957, 0.081076, 0.187479, 0.236314, 0.069664, -0.034321, 0.138018, -0.090784, 0.078337, 0.11314,
               0.072923, -0.069141, 0.155786, 0.095607, -0.013958, 0.25297, 0.143652, 0.085351, -0.191761, 0.110038,
               0.122346, -0.068535, 0.090091, 0.069526, -0.131725, -0.173419, 0.080161, 0.013833, -0.163593, -0.001283,
               -0.011312, -0.020348, -0.004064, 0.089246, 0.115923, -0.001136, 0.108662, 0.044961, 0.150245, -0.355093,
               -0.065505, -0.08575, -0.084668, 0.024026, -0.114004, -0.00582, -0.014616, -0.11759, 0.217696, 0.088498,
               0.078428, -0.24494, 0.108445, -0.085994, 0.010226, -0.012086, 0.062526, 0.032096, 0.15084, 0.187635,
               0.125639, 0.05251, 0.111083, 0.091175, -0.042892, 0.016866, -0.035713, -0.028081, -0.018831, -0.174531,
               -0.084074, 0.047865, 0.122736, 0.11565, -0.131177, 0.187818, 0.123726, -0.065527, 0.019519, 0.037157,
               -0.368623, 0.195931, -0.015609, -0.162628, -0.024976, -0.138307, -0.112833, -0.01342, -0.152109, -0.127772,
               -0.246636, 0.073784, -0.110782, -0.007831, -0.003678, 0.048329, -0.294818, 0.267569, -0.020779, 0.251712,
               -0.087998, -0.010173, -0.012708, -0.014015, 0.206771, 0.188524, -0.057268, -0.050239, 0.098685, -0.18502,
               0.33091, -0.154566, 0.044506, -0.18042, -0.087658, 0.090132, 0.934282, -0.050769, 0.148077, 0.137113,
               -0.141831, 0.107824, 0.549019, -0.137383, 0.08692, -0.036844, 0.092005, 0.782197, 0.003268, 0.06118,
               0.088202, 0.701364, -0.142686, 0.224613, 0.735951, -0.054841, 0.979966)
all_A516 <- pick_allA(p=5, M=1, d=6, params=params516)
all_Omegas516 <- pick_Omegas(p=5, M=1, d=6, params=params516)

params412 <- c(3.201, 1.754, -2.525, 1.181, -1.138, 0.693, -3.132, 5.712, -1.882, 2.174, -1.832, 6.152, -1.5, 3.525,
               -0.297, 1.859, -0.786, 3.065, 2.544, 2.446, 2.805)
all_A412 <- pick_allA(p=4, M=1, d=2, params=params412)
all_Omegas412 <- pick_Omegas(p=4, M=1, d=2, params=params412)

params413 <- c(-0.756, 0.121, 2.396, -1.019, -3.577, -2.476, 1.372, 3.627, 1.099, -1.98, -3.385, 0.109, 0.723, 2.879,
               1.895, -2.923, -4.605, -1.746, 2.225, 3.212, 1.138, -2.113, -4.057, -1.349, 1.13, 1.123, -0.01, -0.357,
               -0.428, -0.677, 1.995, 2.943, 1.059, -0.368, 0.026, 0.268, -0.308, -0.642, -0.022, 1.339, -1.746, 0.607,
               2.511, -1.196, 1.17)
all_A413 <- pick_allA(p=4, M=1, d=3, params=params413)
all_Omegas413 <- pick_Omegas(p=4, M=1, d=3, params=params413)

params512 <- c(-1.413, 1.365, 1.333, 0.491, 0.073, -0.233, -1.222, 0.302, 0.046, 1.803, 1.099, -0.423, -1.833, 0.412,
               -0.238, -0.545, 2.679, -1.216, 0.006, -0.221, -0.975, -0.543, 1.86, 0.459, 1.203)
all_A512 <- pick_allA(p=5, M=1, d=2, params=params512)
all_Omegas512 <- pick_Omegas(p=5, M=1, d=2, params=params512)

params314 <- c(-0.888, 3.587, 1.846, 6.059, 0.54, 0.847, -0.409, -0.053, -1.328, 1.203, -2.687, 0.804, 0.005, -0.067,
               -1.156, 0.765, 0.662, 0.269, 0.584, 0.17, 0.377, -0.612, 1.724, -1.784, 1.445, -0.185, 2.275, -0.267,
               -1.149, 0.441, -1.468, 0.835, 0.361, -0.781, 1.687, 0.039, -0.113, 0.042, -0.81, -0.882, 0.102, -0.316,
               0.774, 0.458, 0.653, 0.406, -0.119, 0.939, -0.481, 0.169, -0.835, 0.549, 0.817, -0.224, -0.043, -0.638,
               3.028, 0.801, 0.083, 0.762, 1.942, 8.485)
all_A314 <- pick_allA(p=3, M=1, d=4, params=params314)
all_Omegas314 <- pick_Omegas(p=3, M=1, d=4, params=params314)



test_that("VAR_pcovmat works correctly", {
  expect_equal(VAR_pcovmat(p=1, d=2, all_Am=array(all_A112[, , , 1], dim=c(2, 2, 1)), Omega_m=all_Omegas112[, , 1]),
               alt_pcovmat(p=1, d=2, all_A=all_A112, all_Omegas=all_Omegas112), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=1, d=5, all_Am=array(all_A115[, , , 1], dim=c(5, 5, 1)), Omega_m=all_Omegas115[, , 1]),
               alt_pcovmat(p=1, d=5, all_A=all_A115, all_Omegas=all_Omegas115), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=2, d=2, all_Am=all_A212[, , , 1], Omega_m=all_Omegas212[, , 1]),
               alt_pcovmat(p=2, d=2, all_A=all_A212, all_Omegas=all_Omegas212), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=2, d=3, all_Am=all_A213[, , , 1], Omega_m=all_Omegas213[, , 1]),
               alt_pcovmat(p=2, d=3, all_A=all_A213, all_Omegas=all_Omegas213), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=3, d=2, all_Am=all_A312[, , , 1], Omega_m=all_Omegas312[, , 1]),
               alt_pcovmat(p=3, d=2, all_A=all_A312, all_Omegas=all_Omegas312), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=7, d=4, all_Am=all_A714[, , , 1], Omega_m=all_Omegas714[, , 1]),
               alt_pcovmat(p=7, d=4, all_A=all_A714, all_Omegas=all_Omegas714), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=5, d=6, all_Am=all_A516[, , , 1], Omega_m=all_Omegas516[, , 1]),
               alt_pcovmat(p=5, d=6, all_A=all_A516, all_Omegas=all_Omegas516), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=4, d=2, all_Am=all_A412[, , , 1], Omega_m=all_Omegas412[, , 1]),
               alt_pcovmat(p=4, d=2, all_A=all_A412, all_Omegas=all_Omegas412), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=4, d=3, all_Am=all_A413[, , , 1], Omega_m=all_Omegas413[, , 1]),
               alt_pcovmat(p=4, d=3, all_A=all_A413, all_Omegas=all_Omegas413), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=5, d=2, all_Am=all_A512[, , , 1], Omega_m=all_Omegas512[, , 1]),
               alt_pcovmat(p=5, d=2, all_A=all_A512, all_Omegas=all_Omegas512), tolerance=1e-5)

  expect_equal(VAR_pcovmat(p=3, d=4, all_Am=all_A314[, , , 1], Omega_m=all_Omegas314[, , 1]),
               alt_pcovmat(p=3, d=4, all_A=all_A314, all_Omegas=all_Omegas314), tolerance=1e-5)
})


test_that("get_Sigmas works correctly", {
  expect_equal(get_Sigmas(p=1, M=1, d=2, all_A=all_A112, all_Omegas=all_Omegas112, all_boldA=all_boldA112),
               array(alt_pcovmat(p=1, d=2, all_A=all_A112, all_Omegas=all_Omegas112), dim=c(2, 2, 1)), tolerance=1e-5)

  expect_equal(get_Sigmas(p=2, M=1, d=3, all_A=all_A213, all_Omegas=all_Omegas213, all_boldA=all_boldA213),
               array(alt_pcovmat(p=2, d=3, all_A=all_A213, all_Omegas=all_Omegas213), dim=c(6, 6, 1)), tolerance=1e-5)

  expect_equal(get_Sigmas(p=7, M=1, d=4, all_A=all_A714, all_Omegas=all_Omegas714, all_boldA=all_boldA714),
               array(alt_pcovmat(p=7, d=4, all_A=all_A714, all_Omegas=all_Omegas714), dim=c(28, 28, 1)), tolerance=1e-5)
})


# p=1, M=2, d=2
theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863, -0.3124, 0.505251, 0.298483, 0.030096,
                   -0.176925, 0.838898, 0.310863, 0.007512, 0.018244, 0.949533, -0.016941, 0.121403, 0.573269)

# p=2, M=2, d=2
theta_222relg <- c(0.356914, 0.107436, 0.356386, 0.08633, 0.13996, 0.035172, -0.164575, 0.386816, 0.451675, 0.013086,
                   0.227882, 0.336084, 0.239257, 0.024173, -0.021209, 0.707502, 0.063322, 0.027287, 0.009182, 0.197066,
                   0.205831, 0.005157, 0.025877, 1.092094, -0.009327, 0.116449, 0.592446)

# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2)
gamma1_222_12_2 <- c(0.1, 0.2, 0.11, 0.22, 0.33)
theta_222log_12_2 <- c(theta_222relg[-length(theta_222relg)], gamma1_222_12_2)

## weight_function == "threshold"

# p=2, M=2, d=2, weight_function="threshold", weighfun_pars=c(2, 1)
theta_222thres_2_1 <- c(theta_222relg[-length(theta_222relg)], 1)

## Constrained models
rbind_diags <- function(p, M, d) {
  I <- diag(p*d^2)
  Reduce(rbind, replicate(M, I, simplify=FALSE))
}

# p=2, M=2, d=2, C_222
C_222 <- rbind_diags(p=2, M=2, d=2)
theta_222relgc <- c(0.36, 0.12, 0.48, 0.07, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, 0.21, 0.01, 0.03, 1.1, 0.01, 0.11, 0.37)

# p=2, M=2, d=2, mean_constraints=list(1:2), AR_constraints=C_222
theta_222relgcm <- c(0.7209658, 0.810858, 0.22, 0.06, -0.15, 0.39, 0.41, -0.01, 0.08, 0.3, 0.21, 0.01,
                     0.03, 1.1, 0.01, 0.11, 0.37)

# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2),
# AR_constraints=C_222, parametrization="mean"
theta_222logcm_12_2 <- c(theta_222relgcm[-length(theta_222relgcm)], gamma1_222_12_2)

# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2),
# AR_constraints=C_222, parametrization="mean",
# weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5), r=c(0, 0.11, 0.12, 0.13, 0))
xi_222logcmw_12_2 <- c(0.22, 0.33)
theta_222logcmw_12_2 <- c(theta_222relgcm[-length(theta_222relgcm)], xi_222logcmw_12_2)

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean"
xi_222logisticcmw_2_1 <- c(0.33)
theta_222logisticcmw_2_1 <- c(theta_222relgcm[-length(theta_222relgcm)], xi_222logisticcmw_2_1)

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean"
xi_222expcmw_2_1 <- c(0.33)
theta_222expcmw_2_1 <- c(theta_222relgcm[-length(theta_222relgcm)], xi_222logisticcmw_2_1)


## Student

# p=2, M=2, d=2, weight_function="threshold", weighfun_pars=c(2, 1), cond_dist="Student"
theta_222threst_2_1 <- c(theta_222thres_2_1, 13)

# p=2, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2), cond_dist="Student", mean_constraints=list(1:2),
# AR_constraints=C_222, parametrization="mean"
theta_222logcmt_12_2 <- c(theta_222logcm_12_2, 2.13)

# p=2, M=2, d=2, weight_function="exponential", weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean"
theta_222expcmwt_2_1 <- c(theta_222expcmw_2_1, 4)

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2), AR_constraints=C_222,
# weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)), parametrization="mean"
theta_222logisticcmwt_2_1 <- theta_222expcmwt_2_1


### Structural models
# (recursively identified models use the same parametrization as reduced form models)

# p=1, M=2, d=2, weight_function="relative_dens", identification="heteroskedasticity"
all_phi_122 <- c(0.734054, 0.225598, 0.705744, 0.187897)
all_A_122 <- c(0.259626, -0.000863, -0.3124, 0.505251, 0.298483, 0.030096, -0.176925, 0.838898)
#all_Omega_122 <- c(0.310863, 0.007512, 0.018244, 0.949533, -0.016941, 0.121403)
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

# p=1, M=2, d=2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=1), identification="heteroskedasticity"
phi10_122 <- c(0.55, 0.11)
A11_122 <- matrix(c(0.34, 0.05, -0.01, 0.72), nrow=2, byrow=FALSE)
Omega1_122 <- matrix(c(0.58, 0.01, 0.01, 0.06), nrow=2, byrow=FALSE)
phi20_122 <- c(0.17, 0.25)
A21_122 <- A11_122
Omega2_122 <- matrix(c(0.50, -0.01, -0.01, 0.20), nrow=2, byrow=FALSE)
alpha1_122 <- 0.60
gamma1_122_12_1 <- c(0.1, 0.2, 0.3)
theta_122logsh_12_1 <- c(phi10_122, phi20_122, vec(A11_122), vec(A21_122), vec(W_122), lambdas_122, gamma1_122_12_1)

# p=1, M=2, d=3, weight_function="exponential", weightfun_pars=c(1, 1), identification="heteroskedasticity"
phi10_123 <- c(1, 2, 3)
A11_123 <- matrix(c(0.1, 0.02, 0.12, 0.3, 0.21, 0.11, 0.05, 0.03, 0.09), nrow=3)
A21_123 <- matrix(c(0.13, 0.03, 0.21, 0.03, 0.14, 0.15, 0.06, 0.07, 0.08), nrow=3)
phi20_123 <- c(0.1, 0.2, 0.3)
alpha1_123 <- 0.6
W_123 <- matrix(c(-0.47, -0.40, 1.25, 0.58, -1.01, 0.18, -0.66, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE)
lambdas_123 <- c(1.56, 1.44, 0.59)
c_and_gamma_123_1_1 <- c(0.1, 0.5)
theta_123expsh_1_1 <- c(phi10_123, phi20_123, vec(A11_123), vec(A11_123), vec(W_123), lambdas_123, c_and_gamma_123_1_1)


# p=2, M=3, d=2, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist="Student",
# identification="heteroskedasticity"
df_232_1_1 <- 10
phi10_232 <- phi10_122; phi20_232 <- phi20_122; phi30_232 <- c(12, 13)
A11_232 <- A11_122; A21_232 <- A21_122; A31_232 <- matrix(c(0.1, 0.0, 0.0, 0.1), nrow=2)
A12_232 <- -A11_122; A22_232 <- -A21_122; A32_232 <- -A31_232
W_232 <- W_122; lambdas2_232 <- lambdas_122; lambdas3_232 <- c(2.1, 0.62)
r1_232_1_1 <- -0.01; r2_232_1_1 <- 1.01
theta_232threstsh_1_1 <- c(phi10_232, phi20_232, phi30_232, vec(A11_232), vec(A12_232), vec(A21_232), vec(A22_232),
                           vec(A31_232), vec(A32_232), vec(W_232), lambdas2_232, lambdas3_232, r1_232_1_1, r2_232_1_1,
                           df_232_1_1)

## Structural models imposing constraints

# p=1, M=2, d=2, weight_function="relative_dens", identification="heteroskedasticity", AR_constraints=C_122
C_122 <- rbind_diags(p=1, M=2, d=2)
W_122 <- matrix(c(-0.03, 0.24, -0.76, -0.02), nrow=2, ncol=2, byrow=FALSE)
lambdas_122 <- c(3.36, 0.86)
theta_122relgshc <- c(phi10_122, phi20_122, vec(A11_122), vec(W_122), lambdas_122, alpha1_122)
theta_122relgshc_expanded <- c(phi10_122, phi20_122, vec(A11_122), vec(A11_122), vec(W_122), lambdas_122, alpha1_122)

# p=1, M=3, d=2, weight_function="relative_dens", identification="heteroskedasticity",
# B_constraints=matrix(c(0.1, NA, 0.3, 0), nrow=2)
W_132b <- matrix(c(0.11, 0.22, 0.33, 0), nrow=2);
phi10_132 <- phi10_232; phi20_132 <- phi20_232; phi30_132 <- phi30_232;
A11_132 <- A11_232; A21_132 <- A21_232; A31_132 <- A31_232;
lambdas2_132 <- lambdas2_232; lambdas3_132 <- lambdas3_232;
alpha1_132 <- 0.5; alpha2_132 <- 0.3
theta_132relgshb <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                      Wvec(W_132b), lambdas2_132, lambdas3_132, alpha1_132, alpha2_132)
theta_132relgshb_expanded <- c(phi10_132, phi20_132, phi30_132, vec(A11_132), vec(A21_132), vec(A31_132),
                               vec(W_132b), lambdas2_132, lambdas3_132, alpha1_132, alpha2_132)

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student",  parametrization="mean"
# identification="heteroskedasticity", mean_constraints=list(1:2), B_constraints=matrix(c(0.1, 0, 0, 0.3), nrow=2)
W_222b <- matrix(c(0.12, 0, 0, 0.31), nrow=2)
phi10_222 <- phi10_122
A11_222 <- A11_122; A12_222 <- -0.2*A11_122
A21_222 <- A21_122; A22_222 <- -0.2*A21_122
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
C_123 <- rbind_diags(p=1, M=2, d=3)
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



test_that("get_regime_means works correctly", {
  expect_equal(c(get_regime_means(p=2, M=2, d=2, params=theta_222thres_2_1, weight_function="threshold", weightfun_pars=c(2, 1))),
               c(0.9600324, 0.5549088, 0.4908410, 1.1693004), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=2, M=2, d=2, params=theta_222expcmw_2_1, weight_function="exponential", parametrization="mean",
                                  weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222,
                                  weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)))),
               c(0.7209658, 0.8108580, 0.7209658, 0.8108580), tolerance=1e-3)

  expect_equal(c(get_regime_means(p=2, M=2, d=2, params=theta_222logcm_12_2, weight_function="mlogit", parametrization="mean",
                                  weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222)),
               c(0.7209658, 0.8108580, 0.7209658, 0.8108580), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=2, M=2, d=2, params=theta_222logcmw_12_2, weight_function="mlogit", parametrization="mean",
                                  weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222,
                                  weight_constraints=list(R=matrix(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 1), nrow=5),
                                                          r=c(0, 0.11, 0.12, 0.13, 0)))),
               c(0.7209658, 0.8108580, 0.7209658, 0.8108580), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=2, M=2, d=2, params=theta_222logcmw_12_2, weight_function="logistic", parametrization="mean",
                                  weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222,
                                  weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)))),
               c(0.7209658, 0.8108580, 0.7209658, 0.8108580), tolerance=1e-3)

  expect_equal(c(get_regime_means(p=1, M=1, d=2, params=params112, weight_function="relative_dens")),
               c(-0.03835678, 6.03323402), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=1, M=1, d=5, params=params115, weight_function="relative_dens")),
               c(2.5744586, 0.3398819, 3.7756708, 7.9074659, 0.3825098), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=2, M=1, d=3, params=params213, weight_function="relative_dens")),
               c(0.2622567, 5.9641285, -1.6752983), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=1, M=2, d=2, params=theta_122relg, weight_function="relative_dens")),
               c(0.7996501, 0.4545899, 0.6798440, 1.2933271), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=2, M=2, d=2, params=theta_222relg, weight_function="relative_dens")),
               c(0.9600324, 0.5549088, 0.4908410, 1.1693004), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=2, M=2, d=2, params=theta_222relgc, weight_function="relative_dens", AR_constraints=C_222)),
               c(0.8730964, 0.5279188, 1.2174281, 0.4221658), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=2, M=2, d=2, params=theta_222relgcm, weight_function="relative_dens", parametrization="mean",
                                  AR_constraints=C_222, mean_constraints=list(1:2))),
               c(0.7209658, 0.8108580, 0.7209658, 0.8108580), tolerance=1e-3)

  expect_equal(c(get_regime_means(p=2, M=2, d=2, params=theta_222log_12_2, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2))),
               c(0.9600324, 0.5549088, 0.4908410, 1.1693004), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=2, M=2, d=2, params=theta_222relgcm, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=2),
                                  parametrization="mean", AR_constraints=C_222, mean_constraints=list(1:2))),
               c(0.7209658, 0.8108580, 0.7209658, 0.8108580), tolerance=1e-3)

  # Student
  expect_equal(c(get_regime_means(p=2, M=2, d=2, params=theta_222threst_2_1, weight_function="threshold", weightfun_pars=c(2, 1),
                                  cond_dist="Student")), c(0.9600324, 0.5549088, 0.4908410, 1.1693004), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=2, M=2, d=2, params=theta_222logcmt_12_2, weight_function="mlogit", cond_dist="Student",
                                  weightfun_pars=list(vars=1:2, lags=2), mean_constraints=list(1:2), AR_constraints=C_222,
                                  parametrization="mean")), c(0.7209658, 0.8108580, 0.7209658, 0.8108580), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=2, M=2, d=2, params=theta_222expcmwt_2_1, weight_function="exponential", cond_dist="Student",
                                  weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222, parametrization="mean",
                                  weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)))),
               c(0.7209658, 0.8108580, 0.7209658, 0.8108580), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=2, M=2, d=2, params=theta_222logisticcmwt_2_1, weight_function="logistic", cond_dist="Student",
                                  weightfun_pars=c(2, 1), mean_constraints=list(1:2), AR_constraints=C_222, parametrization="mean",
                                  weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.01, 0)))),
               c(0.7209658, 0.8108580, 0.7209658, 0.8108580), tolerance=1e-3)

  # Structural
  expect_equal(c(get_regime_means(p=2, M=2, d=2, params=theta_222threst_2_1, weight_function="threshold", weightfun_pars=c(2, 1),
                                  cond_dist="Student", identification="recursive")), c(0.9600324, 0.5549088, 0.4908410, 1.1693004), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=1, M=1, d=5, params=params115, weight_function="relative_dens", identification="recursive")),
               c(2.5744586, 0.3398819, 3.7756708, 7.9074659, 0.3825098), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=1, M=2, d=2, params=theta_122relgsh, weight_function="relative_dens", identification="heteroskedasticity")),
               c(0.7996501, 0.4545899, 0.6798440, 1.2933271), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=2, M=2, d=2, params=theta_222logistictsh_2_1, weight_function="logistic", weightfun_pars=c(2, 1),
                                  cond_dist="Student", identification="heteroskedasticity")),
               c(1.9771355, 1.1584648, 0.9689697, 2.4914094), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=1, M=2, d=2, params=theta_122logsh_12_1, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=1),
                                  identification="heteroskedasticity")),
               c(0.8251484, 0.5402051, 0.2433891, 0.9363195), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=1, M=2, d=3, params=theta_123expsh_1_1,weight_function="exponential", weightfun_pars=c(1, 1),
                                  identification="heteroskedasticity")),
               c(2.2415237, 2.7373745, 3.9231803, 0.2241524, 0.2737375, 0.3923180), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=1, M=2, d=3, params=theta_123expsh_1_1,weight_function="exponential", weightfun_pars=c(1, 1),
                                  identification="heteroskedasticity", parametrization="mean")),
               c(1.0, 2.0, 3.0, 0.1, 0.2, 0.3), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=2, M=3, d=2, params=theta_232threstsh_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                  cond_dist="Student", identification="heteroskedasticity")),
               c(0.55, 0.11, 0.17, 0.25, 12.00, 13.00), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=2, M=3, d=2, params=theta_232threstsh_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                  cond_dist="Student", identification="heteroskedasticity", parametrization="mean")),
               c(0.55, 0.11, 0.17, 0.25, 12.00, 13.00), tolerance=1e-3)

  expect_equal(c(get_regime_means(p=1, M=2, d=2, params=theta_122relgshc, weight_function="relative_dens",
                                  identification="heteroskedasticity", AR_constraints=C_122)),
               c(0.8251484, 0.5402051, 0.2433891, 0.9363195), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=1, M=3, d=2, params=theta_132relgshb, weight_function="relative_dens", identification="heteroskedasticity",
                                  B_constraints=matrix(c(0.1, NA, 0.3, 0), nrow=2))),
               c(0.8251484, 0.5402051, 0.2433891, 0.9363195, 13.3333333, 14.4444444), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=2, M=2, d=2, params=theta_222logistictshmb_2_1, weight_function="logistic", weightfun_pars=c(2, 1),
                                  cond_dist="Student",  parametrization="mean", identification="heteroskedasticity",
                                  mean_constraints=list(1:2), B_constraints=matrix(c(0.1, 0, 0, 0.3), nrow=2))),
               c(0.55, 0.11, 0.55, 0.11), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=1, M=2, d=2, params=theta_122logshwb_12_1, weight_function="mlogit", weightfun_pars=list(vars=1:2, lags=1),
                                  identification="heteroskedasticity", weight_constraints=list(R=0, r=c(0.1, 0.2, 0.3)),
                                  B_constraints=matrix(c(0.1, 0.2, 0.3, 0), nrow=2))),
               c(0.8251484, 0.5402051, 0.2433891, 0.9363195), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=1, M=2, d=3, params=theta_123expshcwb_1_1, weight_function="exponential", weightfun_pars=c(1, 1),
                                  identification="heteroskedasticity", AR_constraints=C_123,
                                  weight_constraints=list(R=matrix(c(1, 0.5), nrow=2), r=c(0, 0)),
                                  B_constraints=matrix(c(-0.47, -0.40, 0, 0.58, -1.01, -0.66, 0, -0.91, -1.19), nrow=3, ncol=3, byrow=FALSE))),
               c(2.2415237, 2.7373745, 3.9231803, 0.2241524, 0.2737375, 0.3923180), tolerance=1e-3)
  expect_equal(c(get_regime_means(p=2, M=3, d=2, params=theta_232threstshb_1_1, weight_function="threshold", weightfun_pars=c(1, 1),
                                  cond_dist="Student", identification="heteroskedasticity", B_constraints=matrix(c(0.1, 0.2, -0.3, 0), nrow=2))),
               c(0.55, 0.11, 0.17, 0.25, 12.00, 13.00), tolerance=1e-3)
})
