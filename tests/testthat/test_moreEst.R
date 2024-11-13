context("moreEst")
library(sstvars)

# The estimation process computationally demanding. Therefore, we only employ minimal
# tests checking that the functions work correctly without producing any errors.

# p=1, M=2, d=2 relative_dens Gaussian STVAR with the means and AR matrices constrained to be identical in both regimes
mod12cm <- STVAR(gdpdef, p=1, M=2, params=c(0.831093, 0.474375, 0.317932, 0.016297, -0.125438, 0.886249, 0.972206,
                                            -0.034546, 0.143727, 0.324064, 0.013734, 0.010002, 0.417756),
                 AR_constraints=rbind(diag(1*2^2), diag(1*2^2)), mean_constraints=list(1:2), parametrization="mean",
                 calc_std_errors=FALSE)
fit12cm <- iterate_more(mod12cm, maxit=30, calc_std_errors=FALSE) # Not found later if defined inside test_that only

test_that("iterate_more works correctly", {
  expect_equal(fit12cm$params, c(0.833884628, 0.475115791, 0.309214073, 0.019662000, -0.120540749, 0.886574549, 0.947106362,
                                 -0.034873197, 0.151487900, 0.319143635, 0.015522353, 0.009552297, 0.417772934), tolerance=1e-4)
})


test_that("fitSSTVAR works correctly", {
  fit12cm_rec <- fitSSTVAR(fit12cm, identification="recursive", calc_std_errors=FALSE)
  expect_equal(fit12cm_rec$params, c(0.833884628, 0.475115791, 0.309214073, 0.019662000, -0.120540749, 0.886574549, 0.947106362,
                                 -0.034873197, 0.151487900, 0.319143635, 0.015522353, 0.009552297, 0.417772934), tolerance=1e-4)

  fit12cm_hsked <- fitSSTVAR(fit12cm, identification="heteroskedasticity", calc_std_errors=FALSE)
  expect_equal(fit12cm_hsked$params, c(0.83388463, 0.47511579, 0.30921407, 0.01966200, -0.12054075, 0.88657455, 0.94258355,
                                       0.06173125, -0.24216237, 0.38428785, 0.35554455, 0.05550897, 0.41777293), tolerance=1e-4)

  fit12cm_hsked2 <- fitSSTVAR(fit12cm_hsked, identification="heteroskedasticity", maxit=100, maxit_robust=100, print_res=FALSE,
                              B_constraints=matrix(c(1, 0, NA, 1), nrow=2), robust_method="Nelder-Mead", calc_std_errors=FALSE)
  expect_equal(fit12cm_hsked2$params, c(0.837936908, 0.477839220, 0.303390364, 0.024522471, -0.122922610, 0.889240290, 0.982442765,
                                        -0.001023027, 0.386047861, 0.313316913, 0.069656755, 0.426723021), tolerance=1e-4)

})

test_that("estim_LS works correctly", {
  # M=1 models
  expect_equal(estim_LS(gdpdef, p=1, M=1, weight_function="threshold", weightfun_pars=c(1, 1), use_parallel=FALSE),
               c(0.64952757, 0.06650693, 0.28852542, 0.02176643, -0.14402552, 0.89710291), tolerance=1e-4)
  expect_equal(estim_LS(gdpdef, p=3, M=1, weight_function="threshold", weightfun_pars=c(1, 1), use_parallel=FALSE),
               c(0.56446992, -0.00112219, 0.24472216, 0.01715902, -0.13104438, 0.61071319, 0.19089925, 0.01448693, 0.15980970,
                 0.13639022, -0.03491110, 0.03055308, -0.17216331, 0.19668905), tolerance=1e-4)
  expect_equal(estim_LS(usamone, p=2, M=1, weight_function="threshold", weightfun_pars=c(1, 1), use_parallel=FALSE),
               c(0.210079734, 0.085210004, -0.029968781, 0.871229633, -0.021505139, 0.159748400, -0.224477388, 0.528035630, -0.042948559,
                 0.007392181, 0.075415536, 1.145718677, -0.083871763, 0.016815872, -0.040748637, 0.206061999, 0.357275942, 0.445073600,
                 -0.051452440, -0.072672532, -0.207935880), tolerance=1e-4)

  expect_equal(estim_LS(gdpdef, p=1, M=1, weight_function="threshold", weightfun_pars=c(1, 1), AR_constraints=diag(1*1*2^2),
                        use_parallel=FALSE),
               c(0.64952757, 0.06650693, 0.28852542, 0.02176643, -0.14402552, 0.89710291), tolerance=1e-4)
  expect_equal(estim_LS(usamone, p=2, M=1, weight_function="threshold", weightfun_pars=c(1, 1), AR_constraints=diag(1*2*3^2),
                        use_parallel=FALSE),
               c(0.210079734, 0.085210004, -0.029968781, 0.871229633, -0.021505139, 0.159748400, -0.224477388, 0.528035630, -0.042948559,
                 0.007392181, 0.075415536, 1.145718677, -0.083871763, 0.016815872, -0.040748637, 0.206061999, 0.357275942, 0.445073600,
                 -0.051452440, -0.072672532, -0.207935880), tolerance=1e-4)

  # M=2 models
  expect_equal(estim_LS(gdpdef, p=1, M=2, weight_function="threshold", weightfun_pars=c(1, 1), use_parallel=FALSE, prefer_stab=TRUE,
                        stab_tol=0.02),
               c(0.66873074, 0.04280673, 1.99561678, 0.05260033, 0.25296743, 0.03227000, -0.16681543, 0.92224400, -0.32268851, 0.07161651,
                 -0.08380163, 0.78673395, 1.67080000), tolerance=1e-4)
  expect_equal(estim_LS(gdpdef, p=1, M=2, weight_function="threshold", weightfun_pars=c(1, 1), AR_constraints=diag(2*1*2^2),
                        use_parallel=FALSE, prefer_stab=FALSE),
               c(0.65960789, 0.03526952, 2.66153644, -0.06036853, 0.25877903, 0.03617575, -0.15627433, 0.92984380, -0.56807458, 0.13839834,
                 -0.10073727, 0.72599353, 1.82058000), tolerance=1e-4)
  expect_equal(estim_LS(gdpdef[1:20,], p=1, M=2, weight_function="threshold", weightfun_pars=c(1, 1),
                        AR_constraints=rbind(diag(1*2^2), diag(1*2^2)), use_parallel=FALSE, prefer_stab=FALSE),
               c(1.1522529, 0.4571413, 1.8530941, 0.1888232, -0.1940042, 0.1311513, -1.3653882, -0.4211367, 0.6749000), tolerance=1e-4)

  expect_equal(estim_LS(usamone, p=3, M=2, weight_function="threshold", weightfun_pars=c(2, 3), use_parallel=FALSE,
                        prefer_stab=FALSE),
               c(0.652576, 0.175326, 0.002851, 0.083691, 0.032548, -0.078254, 0.567213, -0.042183, 0.000817, -1.243842, 0.226208, -0.076736,
                 0.462857, 0.0652, 1.553825, -0.10822, 0.019559, -0.038473, 0.368459, 0.146791, 0.223963, 0.115597, 0.005917, -0.592009,
                 0.028029, -0.046028, 0.048644, -0.985386, 0.392626, -0.161709, -0.531623, -0.06505, 0.041138, 1.076105, 0.00633, 0.262533,
                 0.312795, 0.527391, 0.208268, -0.029051, 0.060797, 1.161442, 0.060284, 0.004856, -0.088378, 0.080577, 0.341107, 0.661855,
                 -0.293046, -0.058348, -0.476601, -0.257079, 0.0168, 0.01489, -0.337263, 0.088686, -0.434576, 0.296081, -0.002405, 0.26017,
                 0.438687), tolerance=1e-4)
  expect_equal(estim_LS(usamone[1:30,], p=2, M=2, weight_function="threshold", weightfun_pars=c(1, 2),
                        AR_constraints=rbind(diag(2*3^2), diag(2*3^2)), use_parallel=FALSE, prefer_stab=TRUE, stab_tol=0.2),
               c(1.74880035, 0.43517156, 0.26252207, 1.67691778, 0.50802265, 0.56614854, 1.00087591, -0.09693911, 0.21331656, -0.09051708,
                 -0.02291946, 0.29522832, -0.50006675, 0.15132866, 1.05887188, -0.15210607, 0.06647288, -0.17334719, -0.17620033,
                 0.50939130, 0.31700377, -0.11856343, -0.23643659, -0.34043967, -0.63293553), tolerance=1e-4)

  expect_equal(estim_LS(gdpdef, p=2, M=2, weight_function="threshold", weightfun_pars=c(2, 1), weight_constraints=list(R=0, r=1),
                        use_parallel=FALSE, prefer_stab=FALSE),
               c(0.49015686, 0.05411437, 1.52344972, 0.09500570, 0.23879165, 0.05260864, -0.11313588, 0.50742467, 0.27536129, 0.02111511,
                 -0.13015735, 0.31240533, 0.15806872, -0.01213539, -0.66667993, 0.75003530, 0.02976018, 0.01282125, 0.07731179, 0.15617378),
               tolerance=1e-4)
  expect_equal(estim_LS(usamone, p=1, M=2, weight_function="threshold", weightfun_pars=c(3, 1), weight_constraints=list(R=0, r=4),
                        use_parallel=FALSE, prefer_stab=TRUE, stab_tol=0.3),
               c(0.216331, 0.267828, -0.048065, 0.427145, 0.055718, 0.184385, 0.719071, -0.008879, 0.097187, -0.077482, 0.480341,
                 0.127052, -0.077407, -0.000328, 1.007551, 0.871984, 0.020279, 0.208262, 0.036304, 0.902758, 0.530248, -0.072479,
                 0.00656, 0.896643),
               tolerance=1e-4)


  # M=3 models
  expect_equal(estim_LS(gdpdef[1:25,], p=1, M=3, weight_function="threshold", weightfun_pars=c(2, 1), use_parallel=FALSE, prefer_stab=FALSE),
               c(4.638511664, 0.850787779, 3.77347740, 0.114723636, 2.120097336, 0.577186587, -1.04435040, 0.056701620, -12.835751344,
                 -3.218735914, 0.415754993, 0.009862679, -10.805408192, 0.655210266, -1.330640729, -0.106618509, 1.256396612, -0.270837751,
                 0.215530000, 0.385060000), tolerance=1e-4)
  #expect_equal(estim_LS(gdpdef[1:25,], p=1, M=3, weight_function="threshold", weightfun_pars=c(2, 1), AR_constraints=diag(3*1*2^2),
  #                      use_parallel=FALSE),
  #             c(4.638511664, 0.850787779, 3.77347740, 0.114723636, 2.120097336, 0.577186587, -1.04435040, 0.056701620, -12.835751344,
  #               -3.218735914, 0.415754993, 0.009862679, -10.805408192, 0.655210266, -1.330640729, -0.106618509, 1.256396612, -0.270837751,
  #               0.215530000, 0.385060000), tolerance=1e-4)
  expect_equal(estim_LS(gdpdef[1:18,], p=1, M=3, weight_function="threshold", weightfun_pars=c(2, 1),
                        AR_constraints=rbind(diag(1*2^2), diag(1*2^2), diag(1*2^2)), use_parallel=FALSE, prefer_stab=TRUE, stab_tol=0.05),
               c(0.2296979897, 0.5741459197, 1.0920239986, 0.5927762583, -0.0005306983, 0.8614160723, 0.2342120451, -0.0227530898,
                 0.9539800003, -1.3059003282, 0.2060200000, 0.2946802878), tolerance=1e-4)

  expect_equal(estim_LS(usamone, p=1, M=3, weight_function="threshold", weightfun_pars=c(2, 1), weight_constraints=list(R=0, r=c(0.5, 1)),
                        use_parallel=FALSE, prefer_stab=FALSE),
               c(0.440583, 0.327383, 0.20408, 0.005529, 0.050251, -0.062963, 1.047048, -0.004801, 0.393968, 0.728018, -0.017939, 0.131222,
                 -1.021616, 0.214329, -0.519795, -0.014421, 0.008286, 0.952177, 0.933074, -0.005109, 0.136908, -0.028199, 0.74166, 0.300739,
                 -0.012535, 0.025099, 0.959378, 0.823539, 0.055133, 0.232099, -0.260387, 0.877925, 0.085019, -0.084918, 0.018238, 0.959502),
               tolerance=1e-4)
  expect_equal(estim_LS(usamone, p=1, M=3, weight_function="threshold", weightfun_pars=c(2, 1), weight_constraints=list(R=0, r=c(0.5, 1)),
                        AR_constraints=diag(3*1*3^2), use_parallel=FALSE, prefer_stab=FALSE),
               c(0.440583, 0.327383, 0.20408, 0.005529, 0.050251, -0.062963, 1.047048, -0.004801, 0.393968, 0.728018, -0.017939, 0.131222,
                 -1.021616, 0.214329, -0.519795, -0.014421, 0.008286, 0.952177, 0.933074, -0.005109, 0.136908, -0.028199, 0.74166, 0.300739,
                 -0.012535, 0.025099, 0.959378, 0.823539, 0.055133, 0.232099, -0.260387, 0.877925, 0.085019, -0.084918, 0.018238, 0.959502),
               tolerance=1e-4)
  expect_equal(estim_LS(usamone, p=2, M=3, weight_function="threshold", weightfun_pars=c(2, 1), weight_constraints=list(R=0, r=c(0.5, 1)),
                        AR_constraints=rbind(diag(2*3^2), diag(2*3^2), diag(2*3^2)), use_parallel=FALSE, prefer_stab=TRUE, stab_tol=0.02),
               c(0.313913403, 0.107614242, 0.001148562, 0.455225826, 0.078647378, 0.174906939, 1.024931371, 0.142042391, 0.477208777,
                 0.877915239, -0.020575558, 0.162885711, -0.646278559, 0.511470637, -0.333893312, 0.001718157, 0.074646405, 1.143012311,
                 -0.082142502, 0.017324091, -0.040528981, 0.187354299, 0.346558229, 0.454233334, -0.044714566, -0.071609548, -0.205052816),
               tolerance=1e-4)

  # M=4 models
  expect_equal(estim_LS(gdpdef, p=2, M=4, weight_function="threshold", weightfun_pars=c(1, 1), weight_constraints=list(R=0, r=c(0, 0.5, 1)),
                        use_parallel=FALSE, prefer_stab=FALSE),
               c(0.60569, 0.06846, 0.77101, -0.16723, 0.28109, 0.15959, 0.64877, 0.02163, 0.44802, 0.0563, -0.06267, 1.10369, 0.05106,
                 0.01208, 0.05355, -0.18523, -0.30673, 0.26092, -0.84912, 0.71243, 0.21501, 0.02838, 0.66451, 0.39614, 0.39249, -0.0088,
                 0.04541, 0.49218, 0.29722, -0.01916, -0.15529, 0.31985, 0.2107, 0.01915, 0.06456, 0.55056, 0.14858, 0.05162, -0.23531,
                 0.33567), tolerance=1e-4)
})
