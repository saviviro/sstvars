context("LSest")
library(sstvars)


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
  expect_equal(estim_LS(gdpdef, p=1, M=2, weight_function="threshold", weightfun_pars=c(1, 1), use_parallel=FALSE),
               c(0.65960789, 0.03526952, 2.66153644, -0.06036853, 0.25877903, 0.03617575, -0.15627433, 0.92984380, -0.56807458,
                 0.13839834, -0.10073727, 0.72599353, 1.82058000), tolerance=1e-4)
  expect_equal(estim_LS(gdpdef, p=1, M=2, weight_function="threshold", weightfun_pars=c(1, 1), AR_constraints=diag(2*1*2^2),
                        use_parallel=FALSE),
               c(0.65960789, 0.03526952, 2.66153644, -0.06036853, 0.25877903, 0.03617575, -0.15627433, 0.92984380, -0.56807458, 0.13839834,
                 -0.10073727, 0.72599353, 1.82058000), tolerance=1e-4)
  expect_equal(estim_LS(gdpdef[1:20,], p=1, M=2, weight_function="threshold", weightfun_pars=c(1, 1),
                        AR_constraints=rbind(diag(1*2^2), diag(1*2^2)), use_parallel=FALSE, min_obs_coef=2.4),
               c(1.1522529, 0.4571413, 1.8530941, 0.1888232, -0.1940042, 0.1311513, -1.3653882, -0.4211367, 0.6749000), tolerance=1e-4)

  expect_equal(estim_LS(usamone, p=3, M=2, weight_function="threshold", weightfun_pars=c(2, 3), use_parallel=FALSE, min_obs_coef=2),
               c(0.652576, 0.175326, 0.002851, 0.083691, 0.032548, -0.078254, 0.567213, -0.042183, 0.000817, -1.243842, 0.226208, -0.076736,
                 0.462857, 0.0652, 1.553825, -0.10822, 0.019559, -0.038473, 0.368459, 0.146791, 0.223963, 0.115597, 0.005917, -0.592009,
                 0.028029, -0.046028, 0.048644, -0.985386, 0.392626, -0.161709, -0.531623, -0.06505, 0.041138, 1.076105, 0.00633, 0.262533,
                 0.312795, 0.527391, 0.208268, -0.029051, 0.060797, 1.161442, 0.060284, 0.004856, -0.088378, 0.080577, 0.341107, 0.661855,
                 -0.293046, -0.058348, -0.476601, -0.257079, 0.0168, 0.01489, -0.337263, 0.088686, -0.434576, 0.296081, -0.002405, 0.26017,
                 0.438687), tolerance=1e-4)

  expect_equal(estim_LS(gdpdef, p=2, M=2, weight_function="threshold", weightfun_pars=c(2, 1), weight_constraints=list(R=0, r=1),
                        use_parallel=FALSE),
               c(0.49015686, 0.05411437, 1.52344972, 0.09500570, 0.23879165, 0.05260864, -0.11313588, 0.50742467, 0.27536129, 0.02111511,
                 -0.13015735, 0.31240533, 0.15806872, -0.01213539, -0.66667993, 0.75003530, 0.02976018, 0.01282125, 0.07731179, 0.15617378),
               tolerance=1e-4)
  expect_equal(estim_LS(usamone, p=1, M=2, weight_function="threshold", weightfun_pars=c(3, 1), weight_constraints=list(R=0, r=4),
                        use_parallel=FALSE),
               c(0.216331, 0.267828, -0.048065, 0.427145, 0.055718, 0.184385, 0.719071, -0.008879, 0.097187, -0.077482, 0.480341,
                 0.127052, -0.077407, -0.000328, 1.007551, 0.871984, 0.020279, 0.208262, 0.036304, 0.902758, 0.530248, -0.072479,
                 0.00656, 0.896643),
               tolerance=1e-4)

  # M=3 models
  expect_equal(estim_LS(gdpdef[1:25,], p=1, M=3, weight_function="threshold", weightfun_pars=c(2, 1), use_parallel=FALSE, min_obs_coef=1.5),
               c(2.36714379530292, 0.348740954662662, 3.94721381125168, -0.00284194993636204, 2.12009733634266, 0.57718658707643,
                 -0.70395526927183, 0.131940113253948, -2.15885122994856, -0.858790670221793, 0.412534469646031, 0.0120419726396386,
                 -11.3221301261482, 1.00487052179712, -1.33064072941913, -0.10661850904367, 1.25639661192596, -0.270837750765991,
                 0.230864371859296, 0.38506), tolerance=1e-4)

  expect_equal(estim_LS(gdpdef[1:25,], p=1, M=3, weight_function="threshold", weightfun_pars=c(2, 1), use_parallel=FALSE, min_obs_coef=1.5,
                        sparse_grid=TRUE),
               c(2.36714379530292, 0.348740954662662, 3.94721381125168, -0.00284194993636204, 2.12009733634266, 0.57718658707643,
                 -0.70395526927183, 0.131940113253948, -2.15885122994856, -0.858790670221793, 0.412534469646031, 0.0120419726396386,
                 -11.3221301261482, 1.00487052179712, -1.33064072941913, -0.10661850904367, 1.25639661192596, -0.270837750765991,
                 0.230864371859296, 0.38506), tolerance=1e-4) # Same as without sparse_grid though

  expect_equal(estim_LS(usamone, p=1, M=3, weight_function="threshold", weightfun_pars=c(2, 1), weight_constraints=list(R=0, r=c(0.5, 1)),
                        use_parallel=FALSE),
               c(0.440582833518931, 0.32738262842581, 0.204080051609945, 0.00552929395908314, 0.0502510939283292, -0.0629632827966964,
                 1.04704831560802, -0.00480099297375198, 0.39396837101859, 0.728017891338885, -0.0179392259744515, 0.131222048924701,
                 -1.02161564038501, 0.214329105514438, -0.519795427706489, -0.0144214859747971, 0.00828622075166788, 0.952177193264171,
                 0.933073770758613, -0.00510853816877957, 0.136907555493998, -0.0281992692985051, 0.741659582427115, 0.300738694026078,
                 -0.0125346313304714, 0.0250986869114289, 0.959378471795302, 0.823538570267156, 0.0551331646339559, 0.23209884715787,
                 -0.260387379191194, 0.877925094739561, 0.0850192474405208, -0.0849180323460721, 0.0182382202227526, 0.959501632761606),
               tolerance=1e-4)

  # M=4 models
  expect_equal(estim_LS(gdpdef, p=2, M=4, weight_function="threshold", weightfun_pars=c(1, 1), weight_constraints=list(R=0, r=c(0, 0.5, 1)),
                        use_parallel=FALSE),
               c(0.60569, 0.06846, 0.77101, -0.16723, 0.28109, 0.15959, 0.64877, 0.02163, 0.44802, 0.0563, -0.06267, 1.10369, 0.05106,
                 0.01208, 0.05355, -0.18523, -0.30673, 0.26092, -0.84912, 0.71243, 0.21501, 0.02838, 0.66451, 0.39614, 0.39249, -0.0088,
                 0.04541, 0.49218, 0.29722, -0.01916, -0.15529, 0.31985, 0.2107, 0.01915, 0.06456, 0.55056, 0.14858, 0.05162, -0.23531,
                 0.33567), tolerance=1e-4)
})


test_that("estim_NLS works correctly", {
  # M=1 models
  expect_equal(estim_NLS(gdpdef, p=1, M=1, weight_function="threshold", weightfun_pars=c(1, 1), use_parallel=FALSE),
               c(0.64952757, 0.06650693, 0.28852542, 0.02176643, -0.14402552, 0.89710291), tolerance=1e-4)
  expect_equal(estim_NLS(gdpdef, p=3, M=1, weight_function="threshold", weightfun_pars=c(1, 1), use_parallel=FALSE),
               c(0.56446992, -0.00112219, 0.24472216, 0.01715902, -0.13104438, 0.61071319, 0.19089925, 0.01448693, 0.15980970,
                 0.13639022, -0.03491110, 0.03055308, -0.17216331, 0.19668905), tolerance=1e-4)
  expect_equal(estim_NLS(usamone, p=2, M=1, weight_function="threshold", weightfun_pars=c(1, 1), use_parallel=FALSE),
               c(0.210079734, 0.085210004, -0.029968781, 0.871229633, -0.021505139, 0.159748400, -0.224477388, 0.528035630, -0.042948559,
                 0.007392181, 0.075415536, 1.145718677, -0.083871763, 0.016815872, -0.040748637, 0.206061999, 0.357275942, 0.445073600,
                 -0.051452440, -0.072672532, -0.207935880), tolerance=1e-4)

  expect_equal(estim_NLS(gdpdef, p=1, M=1, weight_function="threshold", weightfun_pars=c(1, 1), AR_constraints=diag(1*1*2^2),
                         use_parallel=FALSE),
               c(0.64952757, 0.06650693, 0.28852542, 0.02176643, -0.14402552, 0.89710291), tolerance=1e-4)
  expect_equal(estim_NLS(usamone, p=2, M=1, weight_function="threshold", weightfun_pars=c(1, 1), AR_constraints=diag(1*2*3^2),
                         use_parallel=FALSE),
               c(0.210079734, 0.085210004, -0.029968781, 0.871229633, -0.021505139, 0.159748400, -0.224477388, 0.528035630, -0.042948559,
                 0.007392181, 0.075415536, 1.145718677, -0.083871763, 0.016815872, -0.040748637, 0.206061999, 0.357275942, 0.445073600,
                 -0.051452440, -0.072672532, -0.207935880), tolerance=1e-4)

  # M=2 models
  expect_equal(estim_NLS(gdpdef[1:20,], p=1, M=2, weight_function="threshold", weightfun_pars=c(1, 1), use_parallel=FALSE, min_obs_coef=1.1),
               c(0.34473043, 0.43327746, 2.35107664, 0.13900657, -0.17596226, 0.09129371, 0.70432518, -0.24048340, -0.53479980, 0.19431477, -1.08831092,
                 -0.60493940, 0.48942000), tolerance=1e-4)
  expect_equal(estim_NLS(gdpdef[1:20,], p=1, M=2, weight_function="threshold", weightfun_pars=c(1, 1), AR_constraints=diag(2*1*2^2),
                         use_parallel=FALSE, min_obs_coef=1.5),
               c(1.20382493, 0.28795213, 3.21527602, -0.06658233, 0.24818082, -0.02408473, -1.40211728, 0.05933457, -0.97643737, 0.31875529, -1.03523136,
                 -0.75750429, 1.09157000), tolerance=1e-4)
  expect_equal(estim_NLS(gdpdef[1:20,], p=1, M=2, weight_function="threshold", weightfun_pars=c(1, 1),
                        AR_constraints=rbind(diag(1*2^2), diag(1*2^2)), use_parallel=FALSE, min_obs_coef=1.5),
               c(0.8010626, 0.5402487, 2.0314077, 0.1877976, -0.3801979, 0.1567008, -0.8243396, -0.5596906, 0.4894200), tolerance=1e-4)

  #expect_equal(estim_NLS(gdpdef[1:30,], p=1, M=2, weight_function="logistic", weightfun_pars=c(1, 1), use_parallel=FALSE),
  #             c(-0.17291360, 0.37857095, 1.83926953, -0.09261452, 0.02913881, 0.09136100, 3.14098530, -0.02978985, -0.42089036, 0.22668045,
  #               0.57120842, 0.20068665, 0.67490000, 38.46464646), tolerance=1e-4) # Takes too long to compute but can be uncommented for testing

  expect_equal(estim_NLS(gdpdef[1:30,], p=1, M=2, weight_function="logistic", weightfun_pars=c(1, 1), use_parallel=FALSE,
                         weight_constraints=list(R=0, r=c(0.5, 2)), min_obs_coef=1.4),
               c(-1.15590557, 0.82829333, 3.01236194, -0.69803392, -0.63134283, 0.47766326, 3.30057926, 0.02734005, -0.91632485, 0.48707751,
                 0.56006051, 0.17498810), tolerance=1e-4)
  expect_equal(estim_NLS(gdpdef[1:30,], p=1, M=2, weight_function="logistic", weightfun_pars=c(1, 1), AR_constraints=diag(2*1*2^2),
                         use_parallel=FALSE, weight_constraints=list(R=0, r=c(0.5, 2)), min_obs_coef=1.5),
               c(-1.15590557, 0.82829333, 3.01236194, -0.69803392, -0.63134283, 0.47766326, 3.30057926, 0.02734005, -0.91632485, 0.48707751,
                 0.56006051, 0.17498810), tolerance=1e-4)
  expect_equal(estim_NLS(gdpdef[1:20,], p=1, M=2, weight_function="logistic", weightfun_pars=c(1, 1),
                         AR_constraints=rbind(diag(1*2^2), diag(1*2^2)), use_parallel=FALSE, weight_constraints=list(R=0, r=c(0.5, 4)),
                         min_obs_coef=1.1),
               c(0.5727287, 0.7042587, 2.5227436, -0.1069938, -0.5779311, 0.3004377, -1.2132161, -0.4901011), tolerance=1e-4)

  # expect_equal(estim_NLS(gdpdef[1:20,], p=1, M=2, weight_function="logistic", weightfun_pars=c(1, 1),
  #                        AR_constraints=rbind(diag(1*2^2), diag(1*2^2)), use_parallel=FALSE,
  #                        min_obs_coef=1.1, sparse_grid=TRUE),
  #              c(1.42600611, 0.36827734, -0.07213367, 0.68287135, 0.50116985, -0.06466479, -2.07634991, -0.24236217,
  #                1.91955000, 12.90769231), tolerance=1e-4) # Tests sparse_grid but takes too long to be used in automatic testing; uncomment to test in development


  expect_equal(estim_NLS(gdpdef[1:30,], p=1, M=2, weight_function="exponential", weightfun_pars=c(1, 1), use_parallel=FALSE,
                         weight_constraints=list(R=0, r=c(0.5, 2)), min_obs_coef=1.3),
               c(0.17288528, 0.53099751, 0.81191988, 0.18943739, 0.97580964, -0.42119335, 1.72768163, 0.06092438, 0.11608078, 0.07478986,
                 0.33425175, 0.31072873), tolerance=1e-4)
  expect_equal(estim_NLS(gdpdef[1:30,], p=1, M=2, weight_function="exponential", weightfun_pars=c(1, 1), AR_constraints=diag(2*1*2^2),
                         use_parallel=FALSE, weight_constraints=list(R=0, r=c(0.5, 2)), min_obs_coef=1.6),
               c(0.17288528, 0.53099751, 0.81191988, 0.18943739, 0.97580964, -0.42119335, 1.72768163, 0.06092438, 0.11608078, 0.07478986,
                 0.33425175, 0.31072873), tolerance=1e-4)

  expect_equal(estim_NLS(gdpdef[1:30,], p=1, M=2, weight_function="exogenous", use_parallel=FALSE,
                         weightfun_pars=cbind(c(rep(0.2, times=15), rep(0.6, times=14)), c(rep(0.8, times=15), rep(0.4, times=14))),
                         min_obs_coef=2),
               c(1.85072178, 0.15252943, 1.17258260, 0.29669609, -0.79740161, 0.28395474, 3.81914999, 0.10263741, 0.37082979, -0.07181246,
                 -3.17175965, 0.06883980), tolerance=1e-4)
  expect_equal(estim_NLS(gdpdef[1:30,], p=1, M=2, weight_function="exogenous", use_parallel=FALSE, AR_constraints=diag(2*1*2^2),
                         weightfun_pars=cbind(c(rep(0.2, times=15), rep(0.6, times=14)), c(rep(0.8, times=15), rep(0.4, times=14)))),
               c(1.85072178, 0.15252943, 1.17258260, 0.29669609, -0.79740161, 0.28395474, 3.81914999, 0.10263741, 0.37082979, -0.07181246,
                 -3.17175965, 0.06883980), tolerance=1e-4)

  # M=3 models
  # expect_equal(suppressMessages(estim_NLS(gdpdef[1:50,], p=1, M=3, weight_function="mlogit", weightfun_pars=list(1, 1), use_parallel=FALSE)),
  #              c(-0.92276098, -0.47639232, 1.00711707, 0.11671792, 0.96729336, 0.29141893, -0.87444879, -0.22949753, 1.94699160,
  #                1.32570209, 0.22321649, -0.01105715, -0.64472324, 0.91148656, 0.85783938, -0.29305193, -0.37929860, 0.72152604,
  #                -11.11111111, -20.00000000, -15.55555556, 20.00000000), tolerance=1e-4) # Takes too long but can be used for testing

  expect_equal(estim_NLS(gdpdef[1:50,], p=1, M=3, weight_function="mlogit", weightfun_pars=list(1, 1), use_parallel=FALSE,
                                          weight_constraints=list(R=0, r=c(-0.1, 0.1, 0.2, -0.2)), min_obs_coef=1.01),
               c(3474364.6482, -155685.8412, -463167.5732, 22000.0581, -2578017.8168, 113999.8299, -126749.2284, 5620.2018, 195.0502, -1.6739, -28893.0525,
                 1397.3703, 92.7365, 0.6262, -277532.6599, 12668.6575, -290.4225, 3.5395), tolerance=1)

  # expect_equal(estim_NLS(gdpdef[1:50,], p=1, M=3, weight_function="mlogit", weightfun_pars=list(1, 1), use_parallel=FALSE,
  #                        sparse_grid=TRUE, min_obs_coef=1.01),
  #              c(0.357135700, 0.173148447, 2.690586016, -0.081892886, 2.034396871, -0.007423291, -0.527916687, -0.018411182, 0.274082625, 0.839871246,
  #                -0.678212109, 0.130056592, -0.135290813, 0.727828900, -0.772068788, 0.063968248, -0.691329060, 0.931194971, 5.000000000 ,-20.000000000,
  #                -15.000000000, 10.000000000), tolerance=0.01) # Tests sparse_grid but takes too long to be used in automatic testing; uncomment to test in development
})

