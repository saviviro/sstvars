context("standardErrors")
library(sstvars)

# p=1, M=1, d=2
theta_112relg <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103, 0.601786, -0.002945, 0.067224)

# p=1, M=2, d=2
theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863, -0.3124, 0.505251, 0.298483, 0.030096,
                   -0.176925, 0.838898, 0.310863, 0.007512, 0.018244, 0.949533, -0.016941, 0.121403, 0.573269)

pick_Omegas(p=1, M=2, d=2, params=theta_122relg)
diag_Omegas(Omega1=pick_Omegas(p=1, M=2, d=2, params=theta_122relg)[, , 1], Omega2=pick_Omegas(p=1, M=2, d=2, params=theta_122relg)[, , 2])

# p=2, M=2, d=2
theta_222relg <- c(0.356914, 0.107436, 0.356386, 0.08633, 0.13996, 0.035172, -0.164575, 0.386816, 0.451675, 0.013086,
                   0.227882, 0.336084, 0.239257, 0.024173, -0.021209, 0.707502, 0.063322, 0.027287, 0.009182, 0.197066,
                   0.205831, 0.005157, 0.025877, 1.092094, -0.009327, 0.116449, 0.592446)

## Structural
# p=1, M=2, d=2, weight_function="relative_dens", identification="heteroskedasticity"
theta_122relgh <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863, -0.3124, 0.505251, 0.298483, 0.030096,
                   -0.176925, 0.838898, -0.07700669, 0.13124778, -0.55220736, -0.03190641, 6.87154033, 2.98027648, 0.573269)

#mod122relgh <- STVAR(data=gdpdef, p=1, M=2, params=theta_122relgh, weight_function="relative_dens",
#                     identification="heteroskedasticity", calc_std_errors=TRUE)

test_that("standard_errors works correctly", {
  # Relative_dens Gaussian STVAR
  expect_equal(standard_errors(data=gdpdef, p=1, M=1, params=theta_112relg, weight_function="relative_dens"),
               c(0.10110444, 0.03379184, 0.06105941, 0.02040767, 0.08638925, 0.02887363, 0.05459461, 0.01290402, 0.00609869),
               tolerance=1e-2)
  expect_equal(standard_errors(data=gdpdef, p=1, M=2, params=theta_122relg, weight_function="relative_dens"),
               c(0.230509570, 0.064782472, 0.223722542, 0.074611696, 0.121600680, 0.035379322, 0.425862474, 0.121520745,
                 0.093879073, 0.033526832, 0.149737655, 0.050679108, 0.047216323, 0.007952795, 0.003093586, 0.136559303,
                 0.031946911, 0.018560223, 0.109767610),
               tolerance=1e-2)

#  expect_equal(standard_errors(data=gdpdef, p=1, M=2, params=theta_122relgh, weight_function="relative_dens",
#                               identification="heteroskedasticity"),
#               c(0.23063165, 0.06480410, 0.22373391, 0.07461276, 0.12161079, 0.03537993, 0.42608782, 0.12158182, 0.09388320,
#                 0.03352708, 0.14974192, 0.05068029, 0.08269308, 0.01209446, 0.04592398, 0.03088078, 1.58919345, 0.64870247, 0.10980437),
#               tolerance=1e-2)

})
