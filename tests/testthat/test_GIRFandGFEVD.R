context("GIRFandGFEVD")
library(sstvars)

# NOTE: predict uses simulate, so some of the setups are only tested in predict, which breaks if simulate breaks

# p=1, M=1, d=2, identification="recursive"
theta_112rec <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103, 0.601786, -0.002945, 0.067224)
mod112rec <- STVAR(data=gdpdef, p=1, M=1, params=theta_112rec, identification="recursive")

# p=3, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="Student", identification="recursive"
params322logt <- c(0.5959, 0.0447, 2.6279, 0.2897, 0.2837, 0.0504, -0.2188, 0.4008, 0.3128, 0.0271, -0.1194,
                   0.1559, -0.0972, 0.0082, -0.1118, 0.2391, 0.164, -0.0363, -1.073, 0.6759, 3e-04, 0.0069,
                   0.4271, 0.0533, -0.0498, 0.0355, -0.4686, 0.0812, 0.3368, 0.0035, 0.0325, 1.2289, -0.047,
                   0.1666, 1.2067, 7.2392, 11.6091)
mod322logt <- STVAR(gdpdef, p=3, M=2, params=params322logt, weight_function="logistic", weightfun_pars=c(2, 1),
                    cond_dist="Student", identification="recursive")
mod322logt2 <- STVAR(gdpdef[1:60,], p=3, M=2, params=params322logt, weight_function="logistic", weightfun_pars=c(2, 1),
                    cond_dist="Student", identification="recursive")

# data=usamone, p=1, M=2, d=3, weight_function="relative_dens", identification="heteroskedasticity"
theta_123relgh <- c(0.301922, 0.14077, 0.054866, 0.697999, 0.115812, -0.096371, 0.800912, -0.015313, 0.097591,
                    -0.680523, 0.481524, -0.017179, 0.061408, 0.070541, 1.010242, 0.023394, 0.011553, 0.063565,
                    -0.120632, 0.691059, 0.14614, 0.128499, 0.000553, 0.989699, 1.111114, -0.144676, -0.399165,
                    0.845657, 0.091607, 0.875818, -0.540155, -0.411713, 0.080451, 0.603213, 0.174956, 0.087595,
                    0.526162)
mod123relgh <- STVAR(data=usamone, p=1, M=2, params=theta_123relgh, weight_function="relative_dens",
                     identification="heteroskedasticity")
mod123relgh2 <- STVAR(data=usamone[1:50,], p=1, M=2, params=theta_123relgh, weight_function="relative_dens",
                      identification="heteroskedasticity")

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

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="ind_Student", identification="non-Gaussianity"
params222logit <- c(0.428924, 0.153884, 1.455217, 0.264533, 0.314288, 0.124903, -2.633578, -0.313929, 0.293621, 0.036355, -0.273429,
                    0.053141, 0.210305, -0.023328, -0.516816, 0.501481, 0.190278, 0.035209, -0.054981, 0.340776, 1.384283, -0.159169,
                    0.449031, 0.522074, -1.213738, -0.086115, 0.353828, -0.442951, 0.313206, 2.807629, 4.924918, 7.722209)
mod222logit <- STVAR(data=gdpdef, p=2, M=2, params=params222logit, weight_function="logistic", weightfun_pars=c(2, 1),
                     cond_dist="ind_Student", identification="non-Gaussianity")

# p=1, M=2, d=2, weight_function="exogenous", cond_dist="ind_skewed_t", weightfun_pars=twmat, AR_constraints=C_122,
params122exocikt <- c(params122exocit, -0.2, 0.3)
mod122exocikt <- STVAR(data=gdpdef, p=1, M=2, d=2, params=params122exocikt, weight_function="exogenous", cond_dist="ind_skewed_t",
                       weightfun_pars=twmat, AR_constraints=C_122)

# p=2, M=2, d=2, weight_function="logistic", weightfun_pars=c(2, 1), cond_dist="ind_skewed_t", identification="non-Gaussianity"
params222logikt <- c(params222logit, 0.12, -0.17)
mod222logikt <- STVAR(data=gdpdef, p=2, M=2, params=params222logikt, weight_function="logistic", weightfun_pars=c(2, 1),
                      cond_dist="ind_skewed_t", identification="non-Gaussianity")

test_that("GIRF works correctly", {
  girf1 <- GIRF(mod112rec, which_shocks=1:2, shock_size=1, N=3, R1=2, R2=4, init_regime=1, which_cumulative=2,
                use_parallel=FALSE, seeds=1:4, ci=0.9)
  girf2 <- GIRF(mod322logt, which_shocks=1, shock_size=1, N=1, R1=2, R2=3, init_regime=2, which_cumulative=numeric(0),
                ci=c(0.99, 0.5), scale_type="peak", scale=c(1, 1, 0.5), use_parallel=FALSE, seeds=1:3)
  girf3 <- GIRF(mod123relgh, which_shocks=2:3, shock_size=-2, N=2, R1=3, R2=4, init_regime=1, which_cumulative=1:2,
                ci=c(0.95, 0.2), scale_type="instant", scale=cbind(c(2, 1, 0.5), c(3, 3, -0.25)), use_parallel=FALSE,
                seeds=1:4, init_values=array(mod123relgh$data, dim=c(1, 3, 4)))
  girf4 <- GIRF(mod122exocit, which_shocks=2, shock_size=1.3, N=3, R1=2, R2=2, init_regime=1, which_cumulative=1:2,
                ci=c(0.95), scale_type="instant", scale=c(2, 1, 0.3), use_parallel=FALSE, exo_weights=twmat[3:6,],
                seeds=1:2)
  girf5 <- GIRF(mod222logit, which_shocks=1:2, shock_size=0.4, N=2, R1=2, R2=3, init_regime=1, which_cumulative=1,
                ci=c(0.99, 0.2), scale_type="peak", scale=cbind(c(2, 1, 0.5), c(1, 1, -0.25)), use_parallel=FALSE, seeds=1:3)
  girf6 <- GIRF(mod122exocikt, which_shocks=2, shock_size=1.3, N=3, R1=2, R2=2, init_regime=1, which_cumulative=1:2,
                ci=c(0.95), scale_type="instant", scale=c(2, 1, 0.3), use_parallel=FALSE, exo_weights=twmat[3:6,],
                seeds=1:2)
  girf7 <- GIRF(mod222logikt, which_shocks=1:2, shock_size=0.4, N=2, R1=2, R2=3, init_regime=1, which_cumulative=1,
                ci=c(0.99, 0.2), scale_type="peak", scale=cbind(c(2, 1, 0.5), c(1, 1, -0.25)), use_parallel=FALSE, seeds=1:3)
  girf8 <- GIRF(mod222logit, which_shocks=1, N=2, R1=2, init_regime=1, which_cumulative=1,
                ci=c(0.9), use_data_shocks=TRUE, data_girf_pars=c(2, 0.6, -1, 1, 0.8), scale_type="peak",
                scale=c(1, 1, 0.5), use_parallel=FALSE, seeds=1)
  girf9 <- GIRF(mod123relgh, which_shocks=2, N=2, R1=2, init_regime=1, which_cumulative=2,
                ci=c(0.9), use_data_shocks=TRUE, data_girf_pars=c(0, 0.6, 1, 2, 0.8), scale_type="instant",
                scale=c(2, 1, 0.5), use_parallel=FALSE, seeds=1)

  expect_equal(unname(girf1$girf_res$shock1$point_est[4, ]), c(0.00851552, 0.02350278, 0.00000000), tol=1e-4)
  expect_equal(c(unname(girf1$girf_res$shock1$conf_ints[4, , ])), c(0.0006007925, 0.0181154039, 0.0016581835, 0.0499984036,
                                                                    0.0000000000, 0.0000000000), tol=1e-4)
  expect_equal(unname(girf1$girf_res$shock2$point_est[4, ]), c(-0.04182505, 0.86688967, 0.00000000), tol=1e-4)

  expect_equal(unname(girf2$girf_res$shock1$point_est[2, ]), c(0.099155053, -0.028546980, 0.006293229, -0.006293229), tol=1e-4)
  expect_equal(c(unname(girf2$girf_res$shock1$conf_ints[2, , ]))[1:5],
               c(0.09623235, 0.09760175, 0.10063038, 0.10222903, -0.03096390), tol=1e-4)

  expect_equal(unname(girf3$girf_res$shock2$point_est[3, ]), c(1.20918609, 0.18394122, 0.58548529, -0.06394171, 0.06394171), tol=1e-4)
  expect_equal(unname(girf3$girf_res$shock3$point_est[3, ]), c(1.36969884, 2.20495769, -0.09185125, -0.30597982, 0.30597982), tol=1e-4)
  expect_equal(c(unname(girf3$girf_res$shock3$conf_ints[3, ,])),
               c(-0.2820474, 1.2519878, 2.0576334, 2.5367551, 1.5482976, 2.1344350, 2.2240372, 2.9053446, -0.7768956, -0.2943443,
                 0.1092289, 0.5943941, -0.7466761, -0.5207280, -0.3415312, 0.3474711, -0.3474711, 0.3415312, 0.5207280, 0.7466761),
               tol=1e-4)

  expect_equal(unname(girf4$girf_res$shock2$point_est[4, ]), c(0.3782161, -0.7688280, 0.0000000, 0.0000000), tol=1e-4)
  expect_equal(c(unname(girf4$girf_res$shock2$conf_ints[4, ,])), c(0.3782161, 0.3782161, -0.7688280, -0.7688280, 0.0000000,
                                                                   0.0000000, 0.0000000, 0.0000000), tol=1e-4)
  expect_equal(unname(girf5$girf_res$shock1$point_est[3, ]), c(-0.250000000, -0.006672366, 0.002277241, -0.002277241), tol=1e-4)
  expect_equal(c(unname(girf5$girf_res$shock2$conf_ints[3, ,]))[1:4],
               c(0.03261828, 0.40557945, 0.50000000, 0.50000000), tol=1e-4)

  expect_equal(unname(girf6$girf_res$shock2$point_est[4, ]), c(0.3782161, -0.7688280, 0.0000000, 0.0000000), tol=1e-4)
  expect_equal(c(unname(girf6$girf_res$shock2$conf_ints[4, ,])), c(0.3782161, 0.3782161, -0.7688280, -0.7688280, 0.0000000, 0.0000000,
                                                                   0.0000000, 0.0000000), tol=1e-4)
  expect_equal(unname(girf7$girf_res$shock1$point_est[3, ]), c(-0.250000000, -0.008856998, 0.003641080, -0.003641080), tol=1e-4)
  expect_equal(c(unname(girf7$girf_res$shock2$conf_ints[3, ,]))[1:5],
               c(-0.227935215, 0.117582058, 0.264043828, 0.497050548, 0.001396903), tol=1e-4)

  expect_equal(unname(girf8$girf_res$shock1$point_est[3, ]), c(0.48598914, 0.05717324, -0.02531947, 0.02531947), tol=1e-4)
  expect_equal(c(unname(girf8$girf_res$shock1$conf_ints[3, ,]))[1:5],
               c(0.378238636, 0.500000000, 0.005303634, 0.149496147, -0.096839469), tol=1e-4)

  expect_equal(unname(girf9$girf_res$shock2$point_est[3, ]), c(0.42682830, -0.01464495, 0.69076886, 0.01151584, -0.01151584), tol=1e-4)
  expect_equal(c(unname(girf9$girf_res$shock2$conf_ints[3, ,]))[1:5],
               c(0.01336784, 0.94680048, -0.44187445, 0.22429259, 0.37485890), tol=1e-4)
})


test_that("GFEVD works correctly", {
  gfevd1 <- GFEVD(mod112rec, shock_size=1, N=3, initval_type="data", R1=2, R2=4, which_cumulative=2,
                  use_parallel=FALSE, seeds=1:nrow(mod112rec$data))
  gfevd2 <- GFEVD(mod322logt, shock_size=-2.2, N=2, initval_type="random", R1=5, R2=3, init_regime=2,
                  which_cumulative=numeric(0), use_parallel=FALSE, seeds=1:3)
  gfevd3 <- GFEVD(mod123relgh, shock_size=2, N=1, initval_type="fixed", R1=1, R2=4, init_values=array(mod123relgh$data, dim=c(1, 3, 4)),
                  which_cumulative=c(1, 3), use_parallel=FALSE, seeds=1:4)
  gfevd4 <- GFEVD(mod322logt2, use_data_shocks=TRUE, R1=1, N=2, use_parallel=FALSE, seeds=1:(nrow(mod322logt2$data)-3))
  gfevd5 <- GFEVD(mod123relgh2, use_data_shocks=TRUE, R1=2, N=3, which_cumulative=1:2,
                  use_parallel=FALSE, seeds=1:(nrow(mod123relgh2$data)-1))
  gfevd6 <- GFEVD(mod122exocit, use_data_shocks=TRUE, R1=2, N=3, which_cumulative=2,
                  use_parallel=FALSE, seeds=1:(nrow(mod122exocit$data)-1), exo_weights=twmat[3:6,])
  gfevd7 <- GFEVD(mod222logit, use_data_shocks=FALSE, R1=3, R2=2, N=2, which_cumulative=1, init_regime=1,
                  initval_type="random", use_parallel=FALSE, seeds=1:2)
  gfevd8 <- GFEVD(mod222logikt, use_data_shocks=FALSE, R1=1, R2=2, N=2, which_cumulative=1, init_regime=1,
                  initval_type="random", use_parallel=FALSE, seeds=1:2)
  gfevd9 <- GFEVD(mod123relgh2, use_data_shocks=TRUE, data_gfevd_pars=c(1, 0.60), R1=2, N=3, which_cumulative=2,
                  use_parallel=FALSE, seeds=1:nrow(mod123relgh2$data))
  gfevd10 <- GFEVD(mod322logt, initval_type="data", data_gfevd_pars=c(2, 0.70), R1=2, N=2, which_cumulative=2,
                   shock_size=2, use_parallel=FALSE, seeds=1:nrow(mod322logt$data))

  expect_equal(c(unname(gfevd1$gfevd_res[4, , 1:2])), c(0.94476624, 0.05523376, 0.03315255, 0.96684745), tol=1e-4)
  expect_equal(c(unname(gfevd2$gfevd_res[3, ,]))[1:4], c(0.993966532, 0.006033468, 0.032458296, 0.967541704), tol=1e-4)
  expect_equal(c(unname(gfevd3$gfevd_res[2, ,])), c(0.695699319, 0.259120277, 0.045180404, 0.380183965, 0.100550201, 0.519265834, 0.249111917,
                                                    0.747783504, 0.003104579, 0.003353534, 0.895162332, 0.101484133, 0.003353534, 0.895162332,
                                                    0.101484133), tol=1e-4)
  expect_equal(c(unname(gfevd4$gfevd_res[2, 1:2,])), c(0.88593418, 0.11406582, 0.11669351, 0.88330649, 0.04721711, 0.95278289, 0.04721711,
                                                       0.95278289), tol=1e-4)
  expect_equal(c(unname(gfevd4$ind_gfevd_res[1, 2, 1, 1:2])), c(0.3099849, 0.5128602), tol=1e-4)
  expect_equal(c(unname(gfevd5$gfevd_res[1, 1:3, 1:3])), c(0.56943516, 0.27353994, 0.15702490, 0.26239241, 0.11191679, 0.62569081, 0.39344766,
                                                           0.57486996, 0.03168238), tol=1e-4)
  expect_equal(c(unname(gfevd5$ind_gfevd_res[2, 3, 2, 1:2])), c(0.9447989, 0.8897515), tol=1e-4)

  expect_equal(c(unname(gfevd6$gfevd_res[3, 1:2, 1:2])), c(0.6779563, 0.3220437, 0.3568694, 0.6431306), tol=1e-4)
  expect_equal(c(unname(gfevd6$ind_gfevd_res[4, 1:2, 1:2, 20])), c(0.32863707, 0.02956954, 0.67136293, 0.97043046), tol=1e-4)
  expect_equal(c(unname(gfevd7$gfevd_res[3, 1:2, 1:3])), c(0.8828139, 0.1171861, 0.1999002, 0.8000998, 0.1564331, 0.8435669), tol=1e-4)
  expect_equal(c(unname(gfevd8$gfevd_res[3, 1:2, 1:3])), c(0.6839871, 0.3160129, 0.2216477, 0.7783523, 0.1850926, 0.8149074), tol=1e-4)
  expect_equal(c(unname(gfevd9$gfevd_res[4, 1:3, 1:3])), c(0.63494934, 0.24851297, 0.11653769, 0.35639105, 0.14887077, 0.49473818, 0.28038887,
                                                           0.68776030, 0.03185084), tol=1e-4)
  expect_equal(c(unname(gfevd10$gfevd_res[3, 1:2, 1:4])), c(0.87008734, 0.12991266, 0.06337238, 0.93662762, 0.29883251, 0.70116749, 0.29883251,
                                                            0.70116749), tol=1e-4)
})
