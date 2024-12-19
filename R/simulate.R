#' @title Simulate method for class 'stvar' objects
#'
#' @description \code{simulate.stvar} is a simulate method for class 'stvar' objects.
#'
#' @param object an object of class \code{'stvar'}.
#' @param nsim number of observations to be simulated.
#' @param seed set seed for the random number generator?
#' @param ... currently not in use.
#' @param init_values a size \eqn{(p\times d)} matrix specifying the initial values, where d is the number
#'   of time series in the system. The \strong{last} row will be used as initial values for the first lag,
#'   the second last row for second lag etc. If not specified, initial values will be drawn from
#'   the regime specified in \code{init_regimes} (for Gaussian models only).
#' @param init_regime an integer in \eqn{1,...,M} specifying the regime from which
#'   the initial values should be generated from. The initial values will be generated
#'   from the stationary distribution of the specific regime. Due to the lack of
#'   knowledge of the stationary distribution, models with other than Gaussian conditional distribution
#'   uses a simulation procedure with a burn-in period. See the details section.
#' @param ntimes how many sets of simulations should be performed?
#' @param burn_in Burn-in period for simulating initial values from a regime when \code{cond_dist!="Gaussian"}.
#'  See the details section.
#' @param exo_weights if \code{weight_function="exogenous"}, provide a size \eqn{(nsim \times M)} matrix of exogenous
#'  transition weights for the regimes: \code{[t, m]} for a time \eqn{t} and regime \eqn{m} weight. Ignored
#'  if \code{weight_function!="exogenous"}.
#' @param drop if \code{TRUE} (default) then the components of the returned list are coerced to lower dimension if
#'  \code{ntimes==1}, i.e., \code{$sample} and \code{$transition_weights} will be matrices, and \code{$component}
#'   will be vector.
#' @param girf_pars This argument is used internally in the estimation of generalized impulse response functions
#'   (see \code{?GIRF}). You should ignore it (specifying something else than null to it will change how the function behaves).
#' @details The stationary distribution of each regime is not known when \code{cond_dist!="Gaussian"}. Therefore, when using
#'   \code{init_regime} to simulate the initial values from a given regime, we employ the following simulation procedure to
#'   obtain the initial values. First, we set the initial values to the unconditional mean of the specified regime. Then,
#'   we simulate a large number observations from that regime as specified in the argument \code{burn_in}. Then, we simulate
#'   \eqn{p + 100} observations more after the burn in period, and for the \eqn{100} observations calculate the transition
#'   weights for them and take the consecutive \eqn{p} observations that yield the highest transition weight for the given regime.
#'   For models with exogenous transition weights, takes just the last \eqn{p} observations after the burn-in period.
#'
#'   The argument \code{ntimes} is intended for forecasting, which is used by the predict method (see \code{?predict.stvar}).
#' @return Returns a list containing the simulation results. If \code{drop==TRUE} and \code{ntimes==1} (default),
#'   contains the following entries:
#'   \item{sample}{a size (\code{nsim}\eqn{\times d}) matrix containing the simulated time series.}
#'   \item{transition weights:}{a size (\code{nsim}\eqn{\times M}) matrix containing the transition weights corresponding
#'         to the simulated sample.}
#'   Otherwise, returns a list with the following entries:
#'   \item{\code{$sample}}{a size (\code{nsim}\eqn{\times d\times}\code{ntimes}) array containing the samples: the dimension
#'      \code{[t, , ]} is the time index, the dimension \code{[, d, ]} indicates the marginal time series, and the dimension
#'      \code{[, , i]} indicates the i:th set of simulations.}
#'   \item{\code{$transition_weights}}{a size (\code{nsim}\eqn{\times M \times}\code{ntimes}) array containing the transition weights
#'      corresponding to the sample: the dimension \code{[t, , ]} is the time index, the dimension \code{[, m, ]} indicates the
#'      regime, and the dimension \code{[, , i]} indicates the i:th set of simulations.}
#' @seealso \code{\link{predict.stvar}},\code{\link{GIRF}}, \code{\link{GFEVD}},  \code{\link{fitSTVAR}},
#'   \code{\link{fitSSTVAR}} \code{\link{STVAR}}
#' @inherit loglikelihood references
#' @examples
#'  # Gaussian STVAR(p=2, M=2) model with weighted relative stationary densities
#'  # of the regimes as the transition weight function:
#'  theta_222relg <- c(0.356914, 0.107436, 0.356386, 0.08633, 0.13996, 0.035172,
#'    -0.164575, 0.386816, 0.451675, 0.013086, 0.227882, 0.336084, 0.239257, 0.024173,
#'    -0.021209, 0.707502, 0.063322, 0.027287, 0.009182, 0.197066, 0.205831, 0.005157,
#'    0.025877, 1.092094, -0.009327, 0.116449, 0.592446)
#'  mod222relg <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222relg,
#'    weight_function="relative_dens")
#'
#'  # Simulate T=200 observations using given initial values:
#'  init_vals <- matrix(c(0.5, 1.0, 0.5, 1), nrow=2)
#'  sim1 <- simulate(mod222relg, nsim=200, seed=1, init_values=init_vals)
#'  plot.ts(sim1$sample) # Sample
#'  plot.ts(sim1$transition_weights) # Transition weights
#'
#'  # Simulate T=100 observations, with initial values drawn from the stationary
#'  # distribution of the 1st regime:
#'  sim2 <- simulate(mod222relg, nsim=200, seed=1, init_regime=1)
#'  plot.ts(sim2$sample) # Sample
#'  plot.ts(sim2$transition_weights) # Transition weights
#'
#' # Logistic Student's t STVAR with p=1, M=2, and the first lag of the second variable
#' # as the switching variable.
#' params12 <- c(0.62906848, 0.14245295, 2.41245785, 0.66719269, 0.3534745, 0.06041779, -0.34909745,
#'   0.61783824, 0.125769, -0.04094521, -0.99122586, 0.63805416, 0.371575, 0.00314754, 0.03440824,
#'   1.29072533, -0.06067807, 0.18737385, 1.21813844, 5.00884263, 7.70111672)
#' fit12 <- STVAR(data=gdpdef, p=1, M=2, params=params12, weight_function="logistic",
#'   weightfun_pars=c(2, 1), cond_dist="Student")
#'
#' # Simulate T=100 observations with initial values drawn from the second regime.
#' # Since the stationary distribution of the Student's regime is not known, we
#' # use a simulation procedure that starts from the unconditional mean of the regime,
#' # then simulates a number of observations from the regime for a "burn-in" period,
#' # and finally takes the last p observations generated from the regime as the initial
#' # values for the simulation from the STVAR model:
#' sim3 <- simulate(fit12, nsim=100, init_regime=1, burn_in=1000)
#' plot.ts(sim3$sample) # Sample
#' plot.ts(sim3$transition_weights) # Transition weights
#' @export

simulate.stvar <- function(object, nsim=1, seed=NULL, ..., init_values=NULL, init_regime, ntimes=1, burn_in=1000,
                           exo_weights=NULL, drop=TRUE, girf_pars=NULL) {
  # girf_pars$shock_numb - which shock?
  # girf_pars$shock_size - size of the structural shock?
  # Returns a size (N+1 x d+M) vector containing the estimated GIRFs for the variables and
  # the transition weights (column d+m for the m:th regime). The first row for response at impact.

  # Checks etc
  if(!is.null(seed)) set.seed(seed)
  check_stvar(object, object_name="object")
  epsilon <- round(log(.Machine$double.xmin) + 10)
  stvar <- object
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  weight_function <- stvar$model$weight_function
  weightfun_pars <- check_weightfun_pars(data=stvar$data, p=p, M=M, d=d, weight_function=weight_function,
                                         weightfun_pars=stvar$model$weightfun_pars)
  cond_dist <- stvar$model$cond_dist
  identification <- stvar$model$identification
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  weight_constraints <- stvar$model$weight_constraints
  B_constraints <- stvar$model$B_constraints
  allow_unstab <- stvar$allow_unstab # Always FALSE with relative_dens weight function

  # Check the exogenous weights given for simulation
  if(weight_function == "exogenous") {
    check_exoweights(M=M, exo_weights=exo_weights, how_many_rows=nsim, name_of_row_number="nsim")
  }

  if(is.null(init_values) & missing(init_regime)) {
    stop("Either init_values or init_regime needs to be specified")
  }
  if(!missing(init_regime)) {
    stopifnot(init_regime %in% 1:M)
  }
  if(!all_pos_ints(c(nsim, ntimes, burn_in))) stop("Arguments nsim, ntimes, and burn_in must be positive integers")
  if(!is.null(init_values)) {
    if(!is.matrix(init_values)) stop("init_values must be a numeric matrix")
    if(anyNA(init_values)) stop("init_values contains NA values")
    if(ncol(init_values) != d | nrow(init_values) < p) stop("init_values must contain d columns and at least p rows")
  }

  # Collect parameter values
  params <- reform_constrained_pars(p=p, M=M, d=d, params=stvar$params,
                                    weight_function=weight_function, weightfun_pars=weightfun_pars,
                                    cond_dist=cond_dist, identification=identification,
                                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                    weight_constraints=weight_constraints, B_constraints=B_constraints)
  if(stvar$model$parametrization == "mean") {
    params <- change_parametrization(p=p, M=M, d=d, params=params,
                                     weight_function=weight_function, weightfun_pars=weightfun_pars,
                                     identification=identification, cond_dist=cond_dist,
                                     AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL,
                                     B_constraints=NULL, change_to="intercept")
  }

  all_mu <- get_regime_means(p=p, M=M, d=d, params=params,
                             weight_function=weight_function, weightfun_pars=weightfun_pars,
                             cond_dist=cond_dist, parametrization="intercept",
                             identification=identification,
                             AR_constraints=NULL, mean_constraints=NULL,
                             weight_constraints=NULL, B_constraints=NULL) # Not necessarily valid if allow_unstab
  all_phi0 <- pick_phi0(M=M, d=d, params=params)
  all_A <- pick_allA(p=p, M=M, d=d, params=params)
  all_A2 <- array(all_A, dim=c(d, d*p, M)) # cbind coefficient matrices of each component: m:th component is obtained at [, , m]
  all_Omegas <- pick_Omegas(p=p, M=M, d=d, params=params, cond_dist=cond_dist, identification=identification)
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  weightpars <- pick_weightpars(p=p, M=M, d=d, params=params,
                                weight_function=weight_function, weightfun_pars=weightfun_pars,
                                cond_dist=cond_dist)
  distpars <- pick_distpars(d=d, params=params, cond_dist=cond_dist)

  # Structural pars (recursive just takes Cholesky and model identified by non-Gaussianity have B_m in all_Omegas)
  if(identification == "heteroskedasticity") {
    W <- pick_W(p=p, M=M, d=d, params=params, identification=identification)
    lambdas <- matrix(pick_lambdas(p=p, M=M, d=d, params=params, identification=identification), nrow=d, ncol=M-1)
  }

  # Calculate statistics that remain constant through the iterations
  if((cond_dist == "Gaussian" && !allow_unstab) || weight_function == "relative_dens") {
    # Initial regime Gaussian stat dist simu + relative_dens weight function uses this;
    # relative dens only has Gaussian cond dist, but it is used in simulate_from_regime with Student cond dist
     Sigmas <- get_Sigmas(p=p, M=M, d=d, all_A=all_A, all_boldA=all_boldA, all_Omegas=all_Omegas)
     inv_Sigmas <- array(NA, dim=c(d*p, d*p, M)) # Store inverses of the (dpxdp) covariance matrices
     det_Sigmas <- numeric(M) # Store determinants of the (dpxdp) covariance matrices
     chol_Sigmas <- array(dim=c(d*p, d*p, M)) # Cholesky decompositions of the  (dpxdp) covariance matrices
     for(m in 1:M) {
       chol_Sigmas[, , m] <- chol(Sigmas[, , m]) # Upper triangle
       inv_Sigmas[, , m] <- chol2inv(chol_Sigmas[, , m]) # Faster inverse
       det_Sigmas[m] <- prod(diag(chol_Sigmas[, , m]))^2 # Faster determinant
     }
  }
  if(weight_function == "mlogit") {
    all_gamma_m <- matrix(weightpars, ncol=M-1)
    vars <- weightfun_pars[[1]]
    lags <- weightfun_pars[[2]]
    lowers <- (1:lags - 1)*d # We want add vars to each of these
    # Indices of switching variables in cbind(1, Y): we add +1 to the indices since the column of ones on the left,
    # and then the index is added to always account for the constant term.
    inds_of_switching_vars <- c(1, as.vector(matrix(lowers, nrow=length(vars), ncol=length(lowers), byrow=TRUE) + vars + 1))
  }

  # GIRF stuff, particularly for reduced form models, which assume Cholesky identification
  if(!is.null(girf_pars)) {
    R1 <- girf_pars$R1
    all_Omegas_as_matrix <- t(matrix(all_Omegas, nrow=d^2, ncol=M)) # Used for reduced form model GIRF [,m]
  }

  # Set/generate initial values
  if(is.null(init_values)) {
    if(cond_dist == "Gaussian" && !allow_unstab) {
      # Generate the initial values from the stationary distribution of init_regime; Gaussian dist
      mu <- rep(all_mu[, init_regime], p)
      L <- t(chol_Sigmas[, , init_regime]) # Lower triangle
      init_values <- matrix(mu + L%*%rnorm(d*p), nrow=p, ncol=d, byrow=TRUE) # i:th row for the i:th length d random vector
    } else {
      # Generate the initial values using the simulation procedure with a burn-in period.
      init_values <- simulate_from_regime(stvar=stvar, regime=init_regime, nsim=burn_in, use_transweights=TRUE)
    }
  }
  # Take the last p rows of initial values as the initial values
  init_values <- init_values[(nrow(init_values) - p + 1):nrow(init_values), , drop=FALSE]

  # Container for the simulated values and initial values. First row row initial values vector, and t:th row for (y_{i-1},...,y_{i-p})
  Y <- matrix(nrow=nsim + 1, ncol=d*p)
  Y[1,] <- reform_data(init_values, p=p)
  if(!is.null(girf_pars)) Y2 <- Y # Storage for the second sample path in the GIRF algorithm

  # Initialize data structures
  sample <- array(dim=c(nsim, d, ntimes))
  component <- matrix(nrow=nsim, ncol=ntimes)
  transition_weights <- array(dim=c(nsim, M, ntimes))
  if(!is.null(girf_pars)) {
    sample2 <- array(dim=c(nsim, d, ntimes))
    transition_weights2 <- array(dim=c(nsim, M, ntimes))
  }
  if(weight_function == "exogenous") { # Fill in the exogenous weights
    transition_weights <- array(exo_weights, dim=c(nsim, M, ntimes))
  }

  # Some functions to be used
  if(weight_function == "relative_dens") {
    # Get log multivariate normal densities for calculating the transition weights
    get_logmvdvalues <- function(Y, i1) {
      vapply(1:M,
             function(m) -0.5*d*p*log(2*base::pi) - 0.5*log(det_Sigmas[m]) - 0.5*(crossprod(Y[i1,] - rep(all_mu[, m], p),
                                                                                  inv_Sigmas[, , m])%*%(Y[i1,] - rep(all_mu[, m], p))),
             numeric(1))
    } # Returns M x 1 vector; transformed into a matrix in get_alpha_mt
  } else if(weight_function %in% c("logistic", "exponential", "threshold")) {
     # Nothing to do here, can use get_alpha_mt directly
  } else if(weight_function == "mlogit") {
     # Calculate the sub model regressions for calculating the transition weights
     get_regression_values <- function(Y, i1) {
        regressions_mt <- numeric(M) # Vector of zeros
        for(i2 in 1:ncol(all_gamma_m)) {
          regressions_mt[i2] <- crossprod(all_gamma_m[,i2], cbind(1, Y)[i1, inds_of_switching_vars])
        }
        regressions_mt
     }
  } else if(weight_function == "exogenous") {
    # Nothing to do here, uses exo_weights directly
  }

  # Get bounding constants for acceptance-rejection sampling for skewed t-distribution
  if(cond_dist == "ind_skewed_t") {
    all_nu <- distpars[1:d] # df pars
    all_lambda <- distpars[(d + 1):(2*d)] # skewness pars, note that het.sked ident not possible here
    all_bc_M <- vapply(1:d, function(i1) bounding_const_M(nu=all_nu[i1], lambda=all_lambda[i1]), numeric(1)) # bounding constants
  }

  # Run through the time periods and repetitions
  for(j1 in seq_len(ntimes)) {
    if(cond_dist == "ind_skewed_t") { # Genererate sequences of structural shocks here for computational efficiency
      all_e_t <- matrix(nrow=nsim, ncol=d) # [nsim, d]
      for(i1 in 1:d) {
        all_e_t[,i1] <- generate_skewed_t(n=nsim, nu=all_nu[i1], lambda=all_lambda[i1],
                                          bc_M=all_bc_M[i1]) # nsim draws from the i1:th shock
      }
    }
    for(i1 in seq_len(nsim)) {
      # Calculate transition weights
      if(M == 1) {
        alpha_mt <- matrix(1)
      } else {
        if(weight_function == "relative_dens") {
          log_mvdvalues <- get_logmvdvalues(Y=Y, i1=i1)
          alpha_mt <- get_alpha_mt(M=M, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                   weightpars=weightpars, log_mvdvalues=log_mvdvalues, epsilon=epsilon)
        } else if(weight_function %in% c("logistic", "exponential", "threshold")) {
          alpha_mt <- get_alpha_mt(M=M, d=d, Y2=Y[i1, , drop=FALSE], weight_function=weight_function,
                                   weightfun_pars=weightfun_pars, weightpars=weightpars, epsilon=epsilon)
        } else if(weight_function == "mlogit") {
          regression_values <- get_regression_values(Y=Y, i1=i1) # Uses regressions as logmvd values in relative_dens fun
          alpha_mt <- get_alpha_mt(M=M, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                   weightpars=rep(1, times=M), log_mvdvalues=regression_values, epsilon=epsilon)
        } else if(weight_function == "exogenous") {
          alpha_mt <- exo_weights[i1, , drop=FALSE]
        }
      }
      transition_weights[i1, , j1] <- alpha_mt

      # Calculate conditional mean
      mu_yt <- get_mu_yt_Cpp(obs=matrix(Y[i1,], nrow=1), all_phi0=all_phi0, all_A=all_A2, alpha_mt=alpha_mt)

      # Calculate conditional covariance matrix
      if(cond_dist == "ind_Student" || cond_dist == "ind_skewed_t" || identification == "non-Gaussianity") {
        # Paramtrization with impact matrices. Omega_yt is not used anywhere so no need to calculate it.
      } else { # Parametrization with covariance matrices
        Omega_yt <- matrix(rowSums(vapply(1:M, function(m) alpha_mt[, m]*as.vector(all_Omegas[, , m]), numeric(d*d))),
                           nrow=d, ncol=d)
      }

      # Calculate B_t
      if(cond_dist == "ind_Student" || cond_dist == "ind_skewed_t" || identification == "non-Gaussianity") {
        B_t <- matrix(rowSums(vapply(1:M, function(m) alpha_mt[, m]*as.vector(all_Omegas[, , m]), numeric(d*d))),
                      nrow=d, ncol=d) # weighted sum of the impact matrices.
      } else if(identification == "reduced_form") {
        B_t <- matrix(get_symmetric_sqrt(Omega_yt), nrow=d, ncol=d)
      } else if(identification == "recursive") {
        B_t <- t(chol(Omega_yt))
      } else if(identification == "heteroskedasticity") {
        if(M == 1) {
          B_t <- W
        } else {
          tmp <- array(dim=c(d, d, M)) # Store alpha_mt[m]*Lambda_m
          tmp[, , 1] <- alpha_mt[1]*diag(d) # m=1, Lambda = I_d
          for(m in 2:M) {
            tmp[, , m] <- alpha_mt[m]*diag(lambdas[, m - 1])
          }
          B_t <- W%*%sqrt(apply(tmp, MARGIN=1:2, FUN=sum)) # Calculate B_t
        }
      }

      # Draw the structural error
      e_t <- rnorm(d) # This is used to create Student's t variables as well (but will be updated)

      if(cond_dist == "Student") {
        Z <- sqrt((distpars - 2)/distpars)*e_t # Sample from N(0, (df-2)/df*I_d) # df == distpars
        e_t <- Z*sqrt(distpars/rchisq(n=1, df=distpars)) # Sample from t_d(0, I_d, df)
      } else if(cond_dist == "ind_Student") {
        e_t <- sqrt((distpars - 2)/distpars)*rt(n=d, df=distpars) # Sample from independent Student's t distributions
      } else if(cond_dist == "ind_skewed_t") {
        e_t <- all_e_t[i1,] # Sample from independent skewed t distributions (drawn earlier)
      }

      # Calculate the current observation
      sample[i1, , j1] <- t(mu_yt) + B_t%*%e_t

      # Update storage matrix Y (overwrites when ntimes > 1)
      if(p == 1) {
        Y[i1 + 1,] <- sample[i1, , j1]
      } else {
        Y[i1 + 1,] <- c(sample[i1, , j1], Y[i1, 1:(d*p - d)])
      }

      ## For the second sample in GIRF (with a specific structural shock occurring)
      if(!is.null(girf_pars)) {
        # Calculate transition weights
        if(M == 1) {
          alpha_mt2 <- matrix(1)
        } else {
          if(weight_function == "relative_dens") {
            log_mvdvalues2 <- get_logmvdvalues(Y=Y2, i1=i1)
            alpha_mt2 <- get_alpha_mt(M=M, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                      weightpars=weightpars, log_mvdvalues=log_mvdvalues2, epsilon=epsilon)
          } else if(weight_function %in% c("logistic", "exponential", "threshold")) {
            alpha_mt2 <- get_alpha_mt(M=M, d=d, Y2=Y2[i1, , drop=FALSE], weight_function=weight_function,
                                     weightfun_pars=weightfun_pars, weightpars=weightpars, epsilon=epsilon)
          } else if(weight_function == "mlogit") {
            regression_values2 <- get_regression_values(Y=Y2, i1=i1) # Uses regressions as logmvd values in relative_dens fun
            alpha_mt2 <- get_alpha_mt(M=M, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                      weightpars=rep(1, times=M), log_mvdvalues=regression_values2, epsilon=epsilon)
          } else if(weight_function == "exogenous") {
            alpha_mt2 <- exo_weights[i1, , drop=FALSE]
          }
        }
        transition_weights2[i1, , j1] <- alpha_mt2

        # Calculate conditional mean and conditional covariance matrix
        mu_yt2 <- get_mu_yt_Cpp(obs=matrix(Y2[i1,], nrow=1), all_phi0=all_phi0, all_A=all_A2, alpha_mt=alpha_mt2)

        if(cond_dist == "ind_Student" || cond_dist == "ind_skewed_t" || identification == "non-Gaussianity") {
          # Omega_yt2 is not used anywhere so no need to calculate it
        } else { # Parametrization with covariance matrices
          Omega_yt2 <- matrix(rowSums(vapply(1:M, function(m) alpha_mt2[, m]*as.vector(all_Omegas[, , m]),
                                             numeric(d*d))), nrow=d, ncol=d)
        }

        # Calculate B_t
        if(cond_dist == "ind_Student" || cond_dist == "ind_skewed_t" || identification == "non-Gaussianity") {
          B_t2 <- matrix(rowSums(vapply(1:M, function(m) alpha_mt2[, m]*as.vector(all_Omegas[, , m]),
                                        numeric(d*d))), nrow=d, ncol=d) # weighted sum of the impact matrices.
        } else if(identification == "reduced_form") {
          B_t2 <- matrix(get_symmetric_sqrt(Omega_yt2), nrow=d, ncol=d)
        } else if(identification == "recursive") {
          B_t2 <- t(chol(Omega_yt2))
        } else if(identification == "heteroskedasticity") {
          if(M == 1) {
            B_t2 <- W
          } else {
            tmp2 <- array(dim=c(d, d, M)) # Store alpha_mt[m]*Lambda_m
            tmp2[, , 1] <- alpha_mt2[1]*diag(d) # m=1, Lambda = I_d
            for(m in 2:M) {
              tmp2[, , m] <- alpha_mt2[m]*diag(lambdas[, m - 1])
            }
            B_t2 <- W%*%sqrt(apply(tmp2, MARGIN=1:2, FUN=sum)) # Calculate B_t
          }
        }

        # At impact: impose a specific shock to the structural error e_t (which is drawn already for 1st path)
        if(i1 == 1) {
          e_t[girf_pars$shock_numb] <- girf_pars$shock_size
        }

        # Calculate the current observation
        sample2[i1, , j1] <- t(mu_yt2) + B_t2%*%e_t

        # Update storage matrix Y (overwrites when ntimes > 1)
        if(p == 1) {
          Y2[i1 + 1,] <- sample2[i1, , j1]
        } else {
          Y2[i1 + 1,] <- c(sample2[i1, , j1], Y2[i1, 1:(d*p - d)])
        }
      }
    }
  }

  # Calculate a single GIRF for the given structural shock: (N + 1 x d) matrix
  if(!is.null(girf_pars)) {
    one_girf <- apply(X=sample2 - sample, MARGIN=1:2, FUN=mean)
    if(!is.null(stvar$data)) {
      colnames(one_girf) <- colnames(stvar$data)
    } else {
      colnames(one_girf) <- paste("shock", 1:d)
    }
    tw_girf <- apply(X=transition_weights2 - transition_weights, MARGIN=1:2, FUN=mean)
    colnames(tw_girf) <- paste("tw reg.", 1:M)
    return(cbind(one_girf, tw_girf))
  }


  if(ntimes == 1 && drop) {
    sample <- matrix(sample, nrow=nsim, ncol=d)
    transition_weights <- matrix(transition_weights, nrow=nsim, ncol=M)
  }

  list(sample=sample,
       transition_weights=transition_weights)
}



#' @title Simulate observations from a regime of a STVAR model
#'
#' @description \code{simulate_from_regime} allows to simulate observations from a single
#'   regime of a STVAR model
#'
#' @inheritParams simulate.stvar
#' @param stvar an object of class \code{'stvar'}.
#' @param regime an integer in \eqn{1,...,M} determining the regime from which to simulate observations from
#' @param init_values a size \eqn{(p\times d)} matrix specifying the initial values, where d is the number
#'   of time series in the system. The \strong{last} row will be used as initial values for the first lag,
#'   the second last row for second lag etc. If not specified, initial values are set to the unconditional
#'   mean of the regime.
#' @param use_transweights if \code{TRUE} will calculate the transition weights of the provided model, simulate
#'   \eqn{p + 100} observations more, calculate the transition weights for the last \eqn{100} observations, and
#'   return the the consecutive \eqn{p} observations have the highest transition weight for the specified regime.
#' @details Does not take random number generator seed as an argument to avoid unwanted behavior,
#'    because \code{simulate_from_regime} is mostly called from \code{simulate.stvar}
#'    that takes a seed as its argument, and \code{simulate_from_regime} calls \code{simulate.stvar} to simulate the observations.
#'    Specifically, \code{simulate_from_regime} generates a STVAR model from the given regime, sets up the initial values to the
#'    (if not specified), and then calls \code{simulate.stvar} accordingly.
#' @return
#'   \describe{
#'     \item{If \code{use_transweights=FALSE}:}{Returns a \eqn{(nsim \times d)} matrix such that the \eqn{t}th row
#'       contains the \eqn{t}th simulated observation.}
#'     \item{If \code{use_transweights=TRUE}:}{Returns a \eqn{(p \times d)} such that the \eqn{t}th row constrains
#'       the \eqn{t}th observations.}
#'   }
#' @seealso \code{\link{simulate.stvar}}
#' @inherit simulate.stvar references
#' @keywords internal

simulate_from_regime <- function(stvar, regime=1, nsim=1, init_values=NULL, use_transweights=TRUE) {
  check_stvar(stvar)

  # Model specifications
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  weight_function <- stvar$model$weight_function
  weightfun_pars <- check_weightfun_pars(data=stvar$data, p=p, M=M, d=d, weight_function=weight_function,
                                         weightfun_pars=stvar$model$weightfun_pars)
  cond_dist <- stvar$model$cond_dist
  identification <- stvar$model$identification
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  weight_constraints <- stvar$model$weight_constraints
  B_constraints <- stvar$model$B_constraints
  allow_unstab <- stvar$allow_unstab

  # Collect parameter values
  params <- reform_constrained_pars(p=p, M=M, d=d, params=stvar$params,
                                    weight_function=weight_function, weightfun_pars=weightfun_pars,
                                    cond_dist=cond_dist, identification=identification,
                                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                    weight_constraints=weight_constraints, B_constraints=B_constraints)
  if(stvar$model$parametrization == "mean") {
    params <- change_parametrization(p=p, M=M, d=d, params=params,
                                     weight_function=weight_function, weightfun_pars=weightfun_pars,
                                     identification=identification, cond_dist=cond_dist,
                                     AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL,
                                     B_constraints=NULL, change_to="intercept")
  }
  all_phi0 <- pick_phi0(M=M, d=d, params=params)
  all_A <- pick_allA(p=p, M=M, d=d, params=params)
  all_Omegas <- pick_Omegas(p=p, M=M, d=d, params=params, cond_dist=cond_dist, identification=identification)
  distpars <- pick_distpars(d=d, params=params, cond_dist=cond_dist)

  # Create VAR params corresponding the given regime
  # Note that structural models not implemented here: no need here because just simulates initial values to the simulator function
  if(cond_dist == "ind_Student" || cond_dist == "ind_skewed_t" || identification == "non-Gaussianity") {
    new_params <- c(all_phi0[, regime], all_A[, , , regime],
                    vec(all_Omegas[, , regime]), distpars)
  } else {
    new_params <- c(all_phi0[, regime], all_A[, , , regime],
                    vech(all_Omegas[, , regime]), distpars)
  }

  # Note that weight_function does not matter because there is just one regime; also, all constraints are removed prior
  # to building the model.
  new_stvar <- STVAR(p=p, M=1, d=d, params=new_params, weight_function="threshold", weightfun_pars=c(1, 1), cond_dist=cond_dist,
                     parametrization="intercept", identification="reduced_form", AR_constraints=NULL, mean_constraints=NULL,
                     weight_constraints=NULL, B_constraints=NULL, calc_std_errors=FALSE, allow_unstab=allow_unstab)

  if(is.null(init_values)) {
    all_mu <- get_regime_means(p=p, M=M, d=d, params=params,
                               weight_function=weight_function, weightfun_pars=weightfun_pars,
                               cond_dist=cond_dist, parametrization="intercept",
                               identification=identification,
                               AR_constraints=NULL, mean_constraints=NULL,
                               weight_constraints=NULL, B_constraints=NULL)
    init_values <- matrix(rep(all_mu[, regime], times=p), nrow=p, ncol=d, byrow=TRUE)
  } # else: simulate.stvar takes care of hand-specified initial values

  # Simulate and return the sample
  tmp_sim <- simulate.stvar(new_stvar, nsim=ifelse(use_transweights, nsim+100+p, nsim), ntimes=1, init_values=init_values, drop=TRUE)
  ret <- tmp_sim$sample # Sample
  if(use_transweights && stvar$model$weight_function != "exogenous") {
    # Calculate the transition weights, nrow = nrow(ret) - p, as the first p values were used is initial values
    tw <- loglikelihood(data=ret, p=stvar$model$p, M=stvar$model$M, params=stvar$params,
                        weight_function=stvar$model$weight_function, weightfun_pars=stvar$model$weightfun_pars,
                        cond_dist=stvar$model$cond_dist, parametrization=stvar$model$parametrization,
                        identification=stvar$model$identification, AR_constraints=stvar$model$AR_constraints,
                        mean_constraints=stvar$model$mean_constraints, weight_constraints=stvar$model$weight_constraints,
                        B_constraints=stvar$model$B_constraints, to_return="tw", allow_unstab=allow_unstab)

    tw_both <- tw[(nsim + 1):(nrow(tw)), ]
    tw <- tw[(nsim + 1):(nrow(tw)), regime] # Take the transition weights of the last 100 observations (nsim+100+p obs simulated)
    twmax_ind <- which(abs(tw - max(tw)) < 0.001)[1] # Ind with highest tw, but not exactly to avoid overly skewed results
    samp <- ret[(nsim+1):nrow(ret), , drop=FALSE] # Take the last nsim obs
    ret <- samp[twmax_ind:(twmax_ind + p - 1), , drop=FALSE] # Return the previous p obs from the one with highest tw
  }
  ret
}
