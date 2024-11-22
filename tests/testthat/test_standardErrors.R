context("standardErrors")
library(sstvars)

# p=1, M=1, d=2
theta_112relg <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103, 0.601786, -0.002945, 0.067224)

# p=1, M=2, d=2
theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863, -0.3124, 0.505251, 0.298483, 0.030096,
                   -0.176925, 0.838898, 0.310863, 0.007512, 0.018244, 0.949533, -0.016941, 0.121403, 0.573269)

test_that("standard_errors works correctly", {
  # Relative_dens Gaussian STVAR
  expect_equal(standard_errors(data=gdpdef, p=1, M=1, params=theta_112relg, weight_function="relative_dens"),
               c(0.10110451, 0.03379183, 0.06105937, 0.02040767, 0.08638933, 0.02887362, 0.05459424, 0.01290306, 0.00609464),
               tolerance=1e-2)
  expect_equal(standard_errors(data=gdpdef, p=1, M=2, params=theta_122relg, weight_function="relative_dens"),
               c(0.230462862, 0.064763161, 0.223722175, 0.074604006, 0.121597928, 0.035377659, 0.425770535, 0.121462854,
                 0.093878838, 0.033527038, 0.149737260, 0.050667792, 0.047213702, 0.007951471, 0.003068020, 0.136530338,
                 0.031946320, 0.018553981, 0.109743334),
               tolerance=1e-2)

})
