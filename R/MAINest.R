#' @title Two-phase or three-phase (penalized) maximum likelihood estimation of a reduced form smooth transition VAR model
#'
#' @description \code{fitSTVAR} estimates a reduced form smooth transition VAR model in two phases
#'   or three phases. In additional ML estimation, also penalized ML estimation is available.
#'
#' @inheritParams GAfit
#' @param estim_method either \code{"two-phase"} or \code{"three-phase"} (the latter is the default
#'   option for threshold models and the former is currently the only option for other models). See details.
#' @param penalized should penalized log-likelihood function be used that penalizes the log-likelihood function when
#'   the parameter values are close the boundary of the stability region or outside it? If \code{TRUE}, estimates
#'   that do not satisfy the stability condition are allowed (except when \code{weight_function="relative_dens"}).
#'   The default is \code{TRUE} for three-phase estimation and \code{FALSE} for two-phase estimation.
#' @param penalty_params a numeric vector with two positive elements specifying the penalization parameters:
#'   the first element determined how far from the boundary of the stability region the penalization starts
#'   (a number between zero and one, smaller number starts penalization closer to the boundary) and the second element
#'   is a tuning parameter for the penalization (a positive real number, a higher value penalizes non-stability more).
#' @param min_obs_coef In the LS/NLS step of the three phase estimation, the smallest accepted number of observations
#'   (times variables) from each regime relative to the number of parameters in the regime.
#' @param nrounds the number of estimation rounds that should be performed. The default is \code{(M*ncol(data))^3}
#'   when \code{estim_method="two-phase"} and \code{(M*ncol(data))^2} when \code{estim_method="three-phase"}.
#' @param ncores the number CPU cores to be used in parallel computing.
#' @param maxit the maximum number of iterations in the variable metric algorithm.
#' @param seeds a length \code{nrounds} vector containing the random number generator seed
#'  for each call to the genetic algorithm, or \code{NULL} for not initializing the seed.
#' @param print_res should summaries of estimation results be printed?
#' @param use_parallel employ parallel computing? If \code{use_parallel=FALSE && print_res=FALSE},
#'  nothing is printed during the estimation process.
#' @param calc_std_errors Calculate approximate standard errors (based on standard asymptotics)?
#' @param ... additional settings passed to the function \code{GAfit} employing the genetic algorithm.
#' @details
#'  If you wish to estimate a structural model, estimate first the reduced form model and then use the
#'  use the function \code{fitSSTVAR} to create (and estimate if necessary) the structural model
#'  based on the estimated reduced form model.
#'
#'  \strong{three-phase estimation.} With \code{estim_method="three-phase"} (not available for models with \code{relative_dens}
#'  weight function), an extra phase is added to the beginning of the two-phase estimation procedure:
#'  the autoregressive and weight function parameters are first estimated by the method of (penalized) least squares. Then,
#'  the rest of the parameters are estimated by (penalized) ML with the genetic algorithm conditionally on the LS estimates.
#'  Finally, all the parameters are estimated by (penalized) ML by initializing a gradient based variable metric algorithm
#'  from initial estimates obtained from the first two phases. This allows to use substantially decrease the required
#'  number of estimation rounds, and thereby typically speeds up the estimation substantially. On the other hand, the three-phase
#'  procedure tends to produce estimates close to the initial (penalized) LS estimates, while the two-phase procedure explores
#'  the parameter space more thoroughly (when a large enough number of estimation rounds is ran).
#'
#'  \strong{Penalized estimation.} The penalized estimation (\code{penalized=TRUE}) maximizes the penalized log-likelihood function
#'  in which a penalty term is added. The penalty term becomes nonzero when the parameter values are close to the boundary of the
#'  stability region or outside it, it increases in the modulus of the eigenvalues of the companion form AR matrices of the regimes.
#'  With \code{allow_unstab=TRUE}, allowing for unstable estimates, it allows the estimation algorithms to explore the parameter space
#'  outside the stability region, but with high enough penalization, the optimization algorithm eventually converges back to the
#'  stability region. By default, penalized estimation (with unstable estimates allow) is used for \code{estim_method="three-phase"}
#'  and not used for \code{estim_method="two-phase"}.
#'
#'  \strong{The rest concerns both two-phase and three-phase procedures.}\\
#'  Because of complexity and high multimodality of the log-likelihood function, it is \strong{not certain}
#'  that the estimation algorithm will end up in the global maximum point. When \code{estim_method="two-phase"},
#'  it is expected that many of the estimation rounds will end up in some local maximum or a saddle point instead.
#'  Therefore, a (sometimes very large) number of estimation rounds is required for reliable results
#'  (when \code{estim_method="three-phase"} substantially smaller number should be sufficient). Due to
#'  identification problems and high complexity of the surface of the log-likelihood function, the estimation may
#'  fail especially in the cases where the number of regimes is chosen too large.
#'
#'  The estimation process is computationally heavy and it might take considerably long time for large models to
#'  estimate, particularly if \code{estim_method="two-phase"}. Note that reliable estimation of model with
#'  \code{cond_dist == "ind_Student"} or \code{"ind_skewed_t"} is more difficult than with Gaussian or Student's t
#'  models due to the increased complexity.
#'
#'  If the iteration limit \code{maxit} in the variable metric algorithm is reached, one can continue the estimation by
#'  iterating more with the function \code{iterate_more}. Alternatively, one may use the found estimates as starting values
#'  for the genetic algorithm and employ another round of estimation (see \code{??GAfit} how to set up an initial population
#'  with the dot parameters).
#'
#'  \strong{If the estimation algorithm performs poorly}, it usually helps to scale the individual series so that they
#'  vary roughly in the same scale. This makes it is easier to draw reasonable AR coefficients and (with some weight functions)
#'  weight parameter values in the genetic algorithm. Even if the estimation algorithm somewhat works, it should be preferred
#'  to scale the data so that most of the AR coefficients will not be very large, as the genetic algorithm works better with
#'  relatively small AR coefficients. If needed, another package can be used to fit linear VARs to the series to see which scaling
#'  of the series results in relatively small AR coefficients. You should avoid very small (or very high) variance in the data as
#'  well, so that the eigenvalues of the covariance matrices are in a reasonable range.
#'
#'  \strong{weight_constraints:} If you are using weight constraints other than restricting some of the weight parameters to known
#'  constants, make sure the constraints are sensible. Otherwise, the estimation may fail due to the estimation algorithm not being
#'  able to generate reasonable random guesses for the values of the constrained weight parameters.
#'
#'  \strong{Filtering inappropriate estimates:} \code{fitSTVAR} automatically filters through estimates
#'  that it deems "inappropriate". That is, estimates that are not likely solutions of interest.
#'  Specifically, solutions that incorporate a near-singular error term covariance matrix (any eigenvalue less than \eqn{0.002}),
#'  any modulus of the eigenvalues of the companion form AR matrices larger than $0.9985$ (indicating the necessary condition for
#'  stationarity is close to break), or transition weights such that they are close to zero for almost all \eqn{t} for at least
#'  one regime. You can also always find the solutions of interest yourself by using the function \code{alt_stvar} as well since
#'  results from all estimation rounds are saved).
#'
#' @inherit STVAR return
#' @section S3 methods:
#'   The following S3 methods are supported for class \code{'stvar'}: \code{logLik}, \code{residuals}, \code{print}, \code{summary},
#'    \code{predict}, \code{simulate}, and \code{plot}.
#' @seealso \code{\link{fitSSTVAR}}, \code{\link{STVAR}}, \code{\link{GAfit}}, \code{\link{iterate_more}}, \code{\link{filter_estimates}}
#' @references
#'  \itemize{
#'    \item Anderson H., Vahid F. 1998. Testing multiple equation systems for common nonlinear components.
#'      \emph{Journal of Econometrics}, \strong{84}:1, 1-36.
#'    \item Hubrich K., Teräsvirta. T. 2013. Thresholds and Smooth Transitions in Vector Autoregressive Models.
#'      \emph{CREATES Research Paper 2013-18, Aarhus University.}
#'    \item Koivisto T., Luoto J., Virolainen S. 2025. Unpublished working paper.
#'    \item Lanne M., Virolainen S. 2024. A Gaussian smooth transition vector autoregressive model:
#'       An application to the macroeconomic effects of severe weather shocks. Unpublished working
#'       paper, available as arXiv:2403.14216.
#'    \item Kheifets I.L., Saikkonen P.J. 2020. Stationarity and ergodicity of Vector STAR models.
#'      \emph{Econometric Reviews}, \strong{39}:4, 407-414.
#'    \item Tsay R. 1998. Testing and Modeling Multivariate Threshold Models.
#'      \emph{Journal of the American Statistical Association}, \strong{93}:443, 1188-1202.
#'    \item Virolainen S. 2024. Identification by non-Gaussianity in structural threshold and
#'       smooth transition vector autoregressive models. Unpublished working
#'       paper, available as arXiv:2404.19707.
#'  }
#' @examples
#' \donttest{
#' ## These are long running examples. Running all the below examples will take
#' ## approximately three minutes.
#' # When estimating the models in empirical applications, typically a large number
#' # of estimation rounds (set by the argument 'nrounds') should be used. These examples
#' # use only a small number of rounds to make the running time of the examples reasonable.
#'
#' # The below examples make use of the two-variate dataset 'gdpdef' containing
#' # the the quarterly U.S. GDP and GDP deflator from 1947Q1 to 2019Q4.
#'
#' # Estimate Gaussian STVAR model of autoregressive order p=3 and two regimes (M=2),
#' # with the weighted relative stationary densities of the regimes as the transition
#' # weight function. The estimation is performed with 2 rounds and 2 CPU cores, with
#' # the random number generator seeds set for reproducibility (two-phase estimation):
#' fit32 <- fitSTVAR(gdpdef, p=3, M=2, weight_function="relative_dens", cond_dist="Gaussian",
#'  nrounds=2, ncores=2, seeds=1:2)
#'
#' # Examine the results:
#' fit32 # Printout of the estimates
#' summary(fit32) # A more detailed summary printout
#' plot(fit32) # Plot the fitted transition weights
#' get_foc(fit32) # Gradient of the log-likelihood function about the estimate
#' get_soc(fit32) # Eigenvalues of the Hessian of the log-lik. fn. about the estimate
#' profile_logliks(fit32) # Profile log-likelihood functions about the estimate
#'
#' # Estimate a two-regime Student's t STVAR p=5 model with logistic transition weights
#' # and the first lag of the second variable as the switching variable, only two
#' # estimation rounds using two CPU cores (three-phase estimation):
#' fitlogistict32 <- fitSTVAR(gdpdef, p=3, M=2, weight_function="logistic", weightfun_pars=c(2, 1),
#'  cond_dist="Student", nrounds=2, ncores=2, seeds=1:2)
#' summary(fitlogistict32) # Summary printout of the estimates
#'
#' # Estimate a two-regime threshold VAR p=3 model with independent skewed t shocks
#' # using the three-phase estimation procedure.
#' # The first lag of the the second variable is specified as the switching variable,
#' # and the threshold parameter constrained to the fixed value 1 (three-phase estimation):
#' fitthres32wit <- fitSTVAR(gdpdef, p=3, M=2, weight_function="threshold", weightfun_pars=c(2, 1),
#'   cond_dist="ind_skewed_t", weight_constraints=list(R=0, r=1), nrounds=2, ncores=2, seeds=1:2)
#' plot(fitthres32wit) # Plot the fitted transition weights
#'
#' # Estimate a two-regime STVAR p=1 model with exogenous transition weights defined as the indicator
#' # of NBER based U.S. recessions (source: St. Louis Fed database). Moreover, restrict the AR matrices
#' # to be identical across the regimes (i.e., allowing for time-variation in the intercepts and the
#' # covariance matrix only), (three-phase estimation):
#'
#' # Step 1: Define transition weights of Regime 1 as the indicator of NBER based U.S. recessions
#' # (the start date of weights is start of data + p, since the first p values are used as the initial
#' # values):
#' tw1 <- c(0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#'  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1,
#'  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1,
#'  1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#'  1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#'  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#'  0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
#'  0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
#'
#' # Step 2: Define the transition weights of Regime 2 as one minus the weights of Regime 1, and
#' # combine the weights to matrix of transition weights:
#' twmat <- cbind(tw1, 1 - tw1)
#'
#' # Step 3: Create the appropriate constraint matrix:
#' C_122 <- rbind(diag(1*2^2), diag(1*2^2))
#'
#' # Step 4: Estimate the model by specifying the weights in the argument 'weightfun_pars'
#' # and the constraint matrix in the argument 'AR_constraints':
#' fitexo12cit <- fitSTVAR(gdpdef, p=1, M=2, weight_function="exogenous", weightfun_pars=twmat,
#'  cond_dist="ind_Student", AR_constraints=C_122, nrounds=2, ncores=2, seeds=1:2)
#' plot(fitexo12cit) # Plot the transition weights
#' summary(fitexo12cit) # Summary printout of the estimates
#'
#' # Estimate a two-regime Gaussian STVAR p=1 model with the weighted relative stationary densities
#' # of the regimes as the transition weight function; constrain AR matrices to be identical
#' # across the regimes and also constrain the off-diagonal elements of the AR matrices to be zero.
#' # Moreover, constrain the unconditional means of both regimes to be equal (two-phase estimation):
#' mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE)
#' C_222 <- rbind(mat0, mat0) # The constraint matrix
#' fit22cm <- fitSTVAR(gdpdef, p=2, M=2, weight_function="relative_dens", cond_dist="Gaussian",
#'  parametrization="mean", AR_constraints=C_222, mean_constraints=list(1:2), nrounds=2, seeds=1:2)
#' fit22cm # Print the estimates
#'
#' # Estimate a two-regime Student's t STVAR p=3 model with logistic transition weights and the
#' # first lag of the second variable as the switching variable. Constraint the location parameter
#' # to the fixed value 1 and leave the scale parameter unconstrained (three-phase estimation):
#' fitlogistic32w <- fitSTVAR(gdpdef, p=3, M=2, weight_function="logistic", weightfun_pars=c(2, 1),
#'  cond_dist="Student", weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(1, 0)), nrounds=2,
#'  seeds=1:2)
#' plot(fitlogistic32w) # Plot the fitted transition weights
#' }
#' @export

fitSTVAR <- function(data, p, M, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                     weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t"),
                     parametrization=c("intercept", "mean"), AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL,
                     estim_method, penalized, penalty_params=c(0.05, 0.2), allow_unstab, min_obs_coef=3, nrounds, ncores=2, maxit=2000, seeds=NULL,
                     print_res=TRUE, use_parallel=TRUE, calc_std_errors=TRUE, ...) {
  # Initial checks etc
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  if(weight_function == "relative_dens" && M == 1) { # Use threshold weights for a linear model, so enable all appropriate features
    weight_function <- "threshold"
    weightfun_pars <- c(1, 1)
  }
  if(missing(estim_method)) { # Determine the estimation method
    if(weight_function != "relative_dens" && is.null(mean_constraints)) {
      if(is.null(weight_constraints)) {
        estim_method <- "three-phase"
      } else {
        if(!all(weight_constraints$R == 0)) {
          estim_method <- "two-phase"
        } else {
          estim_method <- "three-phase"
        }
      }
    } else {
      estim_method <- "two-phase"
    }
  } else {
    stopifnot(estim_method %in% c("two-phase", "three-phase"))
    if(estim_method == "three-phase" && weight_function == "relative_dens") {
      stop("The three-phase estimation is not available for models with relative_dens weight function.")
    } else if(estim_method == "three-phase" && !is.null(mean_constraints)) {
      stop("The three-phase estimation does not support mean_constraints.")
    } else if(estim_method == "three-phase" && !is.null(weight_constraints) && !all(weight_constraints$R == 0)) {
      stop("The three-phase estimation only supports weight_constraints with R=0.")
    }
  }
  if(missing(penalized)) { # Penalize estimates close to the boundary of the stability region or outside it?
    penalized <- estim_method == "three-phase" # Penalize in three-phase estimation
  }
  stopifnot(is.numeric(penalty_params) && length(penalty_params) == 2 && all(penalty_params >= 0) && penalty_params[1] < 1)
  if(missing(allow_unstab)) {
    allow_unstab <- ifelse(weight_function == "relative_dens", FALSE, penalized) # Allow unstable estimates?

  }
  if(parametrization == "mean") {
    message("Unstable estimates are cannot be allowed with parametrization='mean' (how do you derive the mean?)")
    allow_unstab <- FALSE
  }

  if(missing(nrounds)) {
    if(estim_method == "two-phase") {
      nrounds <- (M*ncol(data))^3
    } else { # three-phase estimation, GA estimation of distribution params only
      nrounds <- (M*ncol(data))^2
    }
  } else {
    stopifnot(length(nrounds) == 1 && all_pos_ints(nrounds))
  }

  if(length(M) != 1 && !all_pos_ints(M)) stop("Argument M must be a positive integer")
  if(M == 1 && weight_function %in% c("logistic", "mlogit", "exponential")) {
    # Set to threshold if only regime (we assume two regimes for logistic and exponential weights)
    weight_function <- "threshold"
    weightfun_pars <- c(1, 1) # There have no affect with M == 1
  }
  check_pMd(p=p, M=M, weight_function=weight_function, identification="reduced_form")
  no_prints <- !use_parallel && !print_res
  if(!all_pos_ints(c(nrounds, ncores, maxit))) stop("Arguments nrounds, ncores, and maxit must be positive integers")
  stopifnot(length(nrounds) == 1)
  if(!is.null(seeds) && length(seeds) != nrounds) stop("The argument 'seeds' should be NULL or a vector of length 'nrounds'")
  data <- check_data(data=data, p=p)
  d <- ncol(data)
  n_obs <- nrow(data)
  T_obs <- n_obs - p
  weightfun_pars <- check_weightfun_pars(data=data, p=p, d=d, M=M, weight_function=weight_function,
                                         weightfun_pars=weightfun_pars, cond_dist=cond_dist)
  if(!is.null(mean_constraints) && parametrization == "intercept") {
    if(!no_prints) message("mean_constraints can be applied for mean-parametrized models only. Switching to parametrization = 'mean'.")
    parametrization <- "mean"
  }
  check_constraints(data=data, p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    parametrization=parametrization, identification="reduced_form",
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=NULL)
  npars <- n_params(p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                    B_constraints=NULL, identification="reduced_form")
  if(npars >= d*T_obs) stop(paste("There are at least as many parameters in the model as there are observations in the data,",
                                  "so you should use smaller p or M."))
  pars_per_reg <- ifelse(is.null(mean_constraints), d, length(mean_constraints)*d/M) + # mean/int params
    ifelse(is.null(AR_constraints), p*d^2, ncol(AR_constraints)/M) + # AR params
    ifelse(cond_dist %in% c("ind_Student", "ind_skewed_t"), d^2, d*(d + 1)/2) # Covmat params

  if(pars_per_reg > d*T_obs/M/2) {
    if(estim_method == "two-phase") { # For three phase, a message is given in LS_estim or NLS_estim
      message("The model has a large number of parameters per regime, which may lead to overfitting.")
    }
  }

  dot_params <- list(...)
  minval <- ifelse(is.null(dot_params$minval), get_minval(data), dot_params$minval)
  red_criteria <- ifelse(rep(is.null(dot_params$red_criteria), 2), c(0.05, 0.01), dot_params$red_criteria)

  # Parallel computing checks
  if(ncores > parallel::detectCores()) {
    ncores <- parallel::detectCores()
    message("ncores was set to be larger than the number of cores detected.")
  }

  if(use_parallel) {
    message(paste("Using", ncores, "cores for", nrounds, "estimations rounds..."))
  }

  if(estim_method == "three-phase") { # Three-phase estimation: initial AR and weight estimates by LS
    ###################################################################
    ### Phase 1: Estimate AR and weight parameters by least squares ##
    ###################################################################

    if(!no_prints && !use_parallel) {
      message(paste0("PHASE 1: Estimating the AR and weight parameters by ",
                     ifelse(weight_function == "threshold", "least squares", "nonlinear least squares"),
                     "...")) # parallel prints inside estim_LS/NLS
    }

    if(weight_function == "threshold") { # Use least squares
      LS_results <- estim_LS(data=data, p=p, M=M, weight_function=weight_function, weightfun_pars=weightfun_pars,
                             cond_dist=cond_dist, parametrization="intercept", AR_constraints=AR_constraints,
                             mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                             penalized=penalized, penalty_params=penalty_params,
                             ncores=ncores, use_parallel=use_parallel) # Always intercept parametrization used here
    } else { # Use nonlinear least squares
      LS_results <- estim_NLS(data=data, p=p, M=M, weight_function=weight_function, weightfun_pars=weightfun_pars,
                              cond_dist=cond_dist, parametrization="intercept", AR_constraints=AR_constraints,
                              mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                              penalized=penalized, penalty_params=penalty_params,
                              ncores=ncores, use_parallel=use_parallel) # Always intercept parametrization used here
    }

    # Check whether the least squares estimates satisfy the stability condition
    if(!is.null(AR_constraints)) { # Expand the AR constraints
      pars_to_check <- c(LS_results[1:(M*d)], AR_constraints%*%LS_results[(M*d + 1):(M*d + ncol(AR_constraints))])
    } else {
      pars_to_check <- LS_results
    }
    all_phi0 <- pick_phi0(M=M, d=d, params=pars_to_check)
    all_A <- pick_allA(p=p, M=M, d=d, params=pars_to_check)
    all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
    stab_ok <- logical(M)
    stab_tol_to_use <- penalty_params[1]
    any_eigen_large <- FALSE
    for(m in 1:M) { # Check stability conditon for each regime
      boldA_mods <- abs(eigen(all_boldA[, , m], symmetric=FALSE, only.values=TRUE)$values)
      stab_ok[m] <- all(boldA_mods < 1 - stab_tol_to_use)
      if(any(boldA_mods > 1.2)) {
        any_eigen_large <- TRUE
      }
    }

    if(!allow_unstab) {
      if(!all(stab_ok)) {
        message("The least squares estimates do not satisfy the stability condition (with a large enough margin)!")
        if(any_eigen_large) {
          message(paste("\nThe", ifelse(weight_function == "threshold", "LS", "NSL"),
                        "estimates a substantially far from satisfying the usual stability condition",
                        "(which will be imposed in the ML estimation)!\n"))
        }
        message(paste("Adjusting the", ifelse(weight_function == "threshold", "LS", "NSL"),
                      "estimates to satisfy the stability condition..."))
        Id <- diag(nrow=d)
        means_prior_to_adjustment <- vapply(1:M, function(m) solve(Id - rowSums(all_A[, , , m, drop=FALSE], dims=2), all_phi0[,m]),
                                            numeric(d)) # [,m]

        for(m in 1:M) {
          if(!stab_ok[m]) {
            # Eigenvalue adjustment approach to force stability to the estimates
            eigen_decomp <- eigen(all_boldA[, , m], symmetric=FALSE)
            eigenvals <- eigen_decomp$values # Eigenvalues, possibly imaginary
            eigenvecs <- eigen_decomp$vectors # Eigenvectors, possibly imaginary
            orig_eigenvals <- eigenvals # Original eigenvalues

            i1 <- 1
            while(TRUE) { # Iteratively adjust eigenvalues until the produced AR matrices are stable
              which_to_adjust <- which(abs(eigenvals) >= 1 - stab_tol_to_use) # Which eigenvalues do not satisfy stability condition

              # Adjust the eigenvalues iteratively:
              eigenvals[which_to_adjust] <- (0.98^i1)*orig_eigenvals[which_to_adjust] # Use original eigenvalues in the scaling
              # Vastly decreasing scaling factor is required to ensure that the method converges to stable estimates

              # Remove the any negligible imaginary part due to numerical error, and obtain the adjusted companion form AR matrix
              # (since conjugate pairs are both equally adjusted, the resulting matrix is real up to numerical error):
              new_boldA <- Re(eigenvecs%*%diag(eigenvals)%*%solve(eigenvecs))

              # Fill in the new AR matrices:
              all_A[, , , m] <- as.vector(new_boldA[1:d, ])

              # Reconstruct and refill the companion form AR matrix:
              all_boldA[, , m] <- form_boldA(p=p, M=1, d=d, all_A=all_A[, , , m, drop=FALSE])

              # Check the stability of the adjusted AR matrix:
              eigenvals <- eigen(all_boldA[, , m], symmetric=FALSE, only.values=TRUE)$values # Updated eigenvalues

              if(all(abs(eigenvals) < 1 - stab_tol_to_use)) { # Check whether the adjusted AR matrices are stable with large enough num tol
                break # Stable estimates are found; break the loop.
              }
              i1 <- i1 + 1
            }

            # Adjust the intercepts so that the means are preserved (but the AR matrices are adjusted):
            all_phi0[, m] <- (Id - apply(all_A[, , , m], MARGIN=1:2, sum))%*%means_prior_to_adjustment[, m]
          }

          # Fill in the adjusted AR matrices:
          if(is.null(AR_constraints)) {
            # We can use the adjusted AR matrices directly:
            LS_results[(M*d + 1):(M*d + M*p*d^2)] <- as.vector(all_A)
          } else { # AR_constraints employed
            # We use Moore-Penrose Pseudoinverse of the AR constraints matrix (which has full column rank)
            # to obtain a stable estimate of psi from the adjusted AR matrices:
            LS_results[(M*d + 1):(M*d + ncol(AR_constraints))] <- solve(crossprod(AR_constraints), AR_constraints)%*%as.vector(all_A)
          }

          # Fill in the adjusted intercepts
          LS_results[1:(M*d)] <- as.vector(all_phi0)
        }
      }
    }

    ## Obtain the LS results with the correct parametrization (to be included in the returned object):
    if(parametrization == "intercept") {
      LS_res_to_ret <- LS_results # Return LS estimates in the intercept paramatrization
    } else { # Parametrization = "mean" and thus, allow_unstab = FALSE
      # Change the LS_res_to_ret to mean parametrization:
      if(is.null(AR_constraints)) {
        mean_and_ar_pars <- LS_results[1:(M*d + M*p*d^2)]
      } else {
        # Expand the AR constraints:
        mean_and_ar_pars <- c(LS_results[1:(M*d)], AR_constraints%*%LS_results[(M*d + 1):(M*d + ncol(AR_constraints))])
      }
      # Change to mean parametrization
      Id <- diag(nrow=d)
      all_phi0 <- pick_phi0(M=M, d=d, params=mean_and_ar_pars)
      all_A <- pick_allA(p=p, M=M, d=d, params=mean_and_ar_pars)
      LS_res_to_ret <- LS_results
      LS_res_to_ret[1:(M*d)] <- vapply(1:M, function(m) solve(Id - rowSums(all_A[, , , m, drop=FALSE], dims=2), all_phi0[,m]), numeric(d))
    }
  }

  #####################################################################
  ## PHASE 1 or 2: estimate all parameters with ta genetic algorithm ##
  #####################################################################

  ### Optimization with the genetic algorithm ###
  GA_parametrization <- ifelse(estim_method == "three-phase", "intercept", parametrization) # parametrization to be used GAfit
  which_phase <- ifelse(estim_method == "three-phase", "PHASE 2", "PHASE 1")
  which_pars_est <- ifelse(estim_method == "three-phase", "the error distribution", "all the")

  if(estim_method == "three-phase") {
    fixed_params <- LS_results
  } else {
    fixed_params <- NULL
  }

  if(!no_prints) message(paste0(which_phase, ": Estimating ", which_pars_est, " parameters with a genetic algorithm..."))

  if(use_parallel) {
    cl <- parallel::makeCluster(ncores)
    on.exit(try(parallel::stopCluster(cl), silent=TRUE)) # Close the cluster on exit, if not already closed.
    parallel::clusterExport(cl, ls(environment(fitSTVAR)), envir=environment(fitSTVAR)) # assign all variables from package:sstvars
    parallel::clusterEvalQ(cl, c(library(pbapply), library(Rcpp), library(RcppArmadillo), library(sstvars)))

    GAresults <- pbapply::pblapply(1:nrounds, function(i1) GAfit(data=data, p=p, M=M,
                                                                 weight_function=weight_function,
                                                                 weightfun_pars=weightfun_pars,
                                                                 cond_dist=cond_dist,
                                                                 parametrization=GA_parametrization,
                                                                 AR_constraints=AR_constraints,
                                                                 mean_constraints=mean_constraints,
                                                                 weight_constraints=weight_constraints,
                                                                 fixed_params=fixed_params,
                                                                 penalized=penalized,
                                                                 penalty_params=penalty_params,
                                                                 allow_unstab=allow_unstab,
                                                                 seed=seeds[i1], ...), cl=cl)

  } else {
    tmpfunGA <- function(i1, ...) {
      if(!no_prints) message(i1, "/", nrounds, "\r")
      GAfit(data=data, p=p, M=M,
            weight_function=weight_function,
            weightfun_pars=weightfun_pars,
            cond_dist=cond_dist,
            parametrization=GA_parametrization,
            AR_constraints=AR_constraints,
            mean_constraints=mean_constraints,
            weight_constraints=weight_constraints,
            fixed_params=fixed_params,
            penalized=penalized,
            penalty_params=penalty_params,
            allow_unstab=allow_unstab,
            seed=seeds[i1], ...)
    }
    GAresults <- lapply(1:nrounds, function(i1) tmpfunGA(i1, ...))
  }

  ## Calculate the log-likelihoods
  loks <- vapply(1:nrounds, function(i1) loglikelihood(data=data, p=p, M=M, params=GAresults[[i1]],
                                                       weight_function=weight_function, weightfun_pars=weightfun_pars,
                                                       cond_dist=cond_dist, parametrization=GA_parametrization,
                                                       identification="reduced_form", AR_constraints=AR_constraints,
                                                       mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                                       B_constraints=NULL, to_return="loglik", check_params=TRUE,
                                                       penalized=penalized, penalty_params=penalty_params, allow_unstab=allow_unstab,
                                                       minval=minval, alt_par=TRUE), numeric(1)) # GA results always in alt_par

  if(print_res) {
    print_loks <- function(loks) {
      printfun <- function(txt, FUN) message(paste(txt, round(FUN(loks), 3)))
      printfun("The lowest loglik: ", min)
      printfun("The largest loglik:", max)
    }
    message("Results from the genetic algorithm:")
    print_loks(loks)
  }

  ## Change to intercept parametrization if parametrization == "intercept" and GA returns mean parametrization:
  if(parametrization == "intercept" && GA_parametrization == "mean") {
    GAresults <- lapply(GAresults, function(pars) change_parametrization(p=p, M=M, d=d, params=pars,
                                                                         weight_function=weight_function,
                                                                         weightfun_pars=weightfun_pars,
                                                                         cond_dist=cond_dist,
                                                                         identification="reduced_form",
                                                                         AR_constraints=AR_constraints,
                                                                         mean_constraints=mean_constraints,
                                                                         weight_constraints=weight_constraints,
                                                                         B_constraints=NULL,
                                                                         change_to="intercept"))
  }
  # Below this parametrization = parametrization can be used.

  #####################################################################
  ### PHASE 2 or 3: Optimization with the variable metric algorithm ###
  #####################################################################

  ## A function to calculate the log-likelihood function
  loglik_fn <- function(params) {
    tryCatch(loglikelihood(data=data, p=p, M=M, params=params,
                           weight_function=weight_function, weightfun_pars=weightfun_pars,
                           cond_dist=cond_dist, parametrization=parametrization,
                           identification="reduced_form", AR_constraints=AR_constraints,
                           mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                           B_constraints=NULL, to_return="loglik", check_params=TRUE, penalized=penalized,
                           penalty_params=penalty_params,
                           allow_unstab=allow_unstab, minval=minval, alt_par=TRUE), # alt_par used in the GA estimation
             error=function(e) minval)
  }

  ## A function to calculate the gradient of the log-likelihood function# using central difference approximation:
  h <- 1e-3
  I <- diag(rep(1, times=npars))
  loglik_grad <- function(params) {
    vapply(1:npars, function(i1) (loglik_fn(params + I[i1,]*h) - loglik_fn(params - I[i1,]*h))/(2*h),
           numeric(1))
  }

  which_phase <- ifelse(estim_method == "three-phase", "PHASE 3", "PHASE 2")
  if(!no_prints) message(paste0(which_phase, ": Estimating all the parameters with a variable metric algorithm..."))
  if(use_parallel) {
    NEWTONresults <- pbapply::pblapply(1:nrounds, function(i1) optim(par=GAresults[[i1]], fn=loglik_fn, gr=loglik_grad,
                                                                     method="BFGS", control=list(fnscale=-1, maxit=maxit)), cl=cl)
    parallel::stopCluster(cl=cl)
  } else {
    tmpfunNE <- function(i1) {
      if(!no_prints) message(i1, "/", nrounds, "\r")
      optim(par=GAresults[[i1]], fn=loglik_fn,  gr=loglik_grad,
            method="BFGS", control=list(fnscale=-1, maxit=maxit))
    }
    NEWTONresults <- lapply(1:nrounds, function(i1) tmpfunNE(i1))
  }

  loks <- vapply(1:nrounds, function(i1) NEWTONresults[[i1]]$value, numeric(1)) # Log-likelihoods
  converged <- vapply(1:nrounds, function(i1) NEWTONresults[[i1]]$convergence == 0, logical(1)) # Which coverged

  if(print_res) {
    message("Results from the variable metric algorithm:")
    print_loks(loks)
  }

  ### Obtain estimates, change back to original parametrization, and filter the inapproriate estimates
  all_estimates <- lapply(NEWTONresults, function(x) x$par)
  if(cond_dist == "ind_Student" || cond_dist == "ind_skewed_t") {
    all_estimates <- lapply(all_estimates, function(pars) change_parametrization(p=p, M=M, d=d, params=pars,
                                                                                 weight_function=weight_function,
                                                                                 weightfun_pars=weightfun_pars,
                                                                                 cond_dist=cond_dist,
                                                                                 identification="reduced_form",
                                                                                 AR_constraints=AR_constraints,
                                                                                 mean_constraints=mean_constraints,
                                                                                 weight_constraints=weight_constraints,
                                                                                 B_constraints=NULL,
                                                                                 change_to="orig"))
  }

  ### Filter inappropriate estimates ###
  if(!no_prints) message("Filtering inappropriate estimates...")
  ord_by_loks <- order(loks, decreasing=TRUE) # Ordering from the largest loglik to the smallest

  # Go through estimates, take the estimate that yield the higher likelihood
  # among estimates that are do not include wasted regimes or near-singular
  # error term covariance matrices.
  for(i1 in 1:length(all_estimates)) {
    which_round <- ord_by_loks[i1] # Est round with i1:th largest loglik
    pars <- all_estimates[[which_round]]
    pars_std <- reform_constrained_pars(p=p, M=M, d=d, params=pars,
                                        weight_function=weight_function, weightfun_pars=weightfun_pars,
                                        cond_dist=cond_dist, identification="reduced_form",
                                        AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                        weight_constraints=weight_constraints,
                                        B_constraints=NULL) # Pars in standard form for pick pars fns
    # Check Omegas
    Omega_eigens <- get_omega_eigens_par(p=p, M=M, d=d, params=pars_std,
                                         weight_function=weight_function, weightfun_pars=weightfun_pars,
                                         cond_dist=cond_dist, identification="reduced_form",
                                         AR_constraints=NULL, mean_constraints=NULL,
                                         weight_constraints=NULL, B_constraints=NULL)
    Omegas_ok <- !any(Omega_eigens < 0.002)

    # Checks AR matrices
    boldA_eigens <- get_boldA_eigens_par(p=p, M=M, d=d, params=pars_std,
                                         weight_function=weight_function, weightfun_pars=weightfun_pars,
                                         cond_dist=cond_dist, identification="reduced_form",
                                         AR_constraints=NULL, mean_constraints=NULL,
                                         weight_constraints=NULL, B_constraints=NULL)
    stat_ok <- !any(boldA_eigens > 0.9985)

    # Check weight parameters
    if(weight_function == "relative_dens") {
      alphas <- pick_weightpars(p=p, M=M, d=d, params=pars_std, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                cond_dist=cond_dist)
      weightpars_ok <- !any(alphas < 0.01)
    } else {
      weightpars_ok <- TRUE
    }

    # Check transition weights
    tweights <- loglikelihood(data=data, p=p, M=M, params=pars_std,
                              weight_function=weight_function, weightfun_pars=weightfun_pars,
                              cond_dist=cond_dist, parametrization=parametrization,
                              identification="reduced_form", AR_constraints=NULL,
                              mean_constraints=NULL, B_constraints=NULL, weight_constraints=NULL,
                              to_return="tw", check_params=TRUE, penalized=penalized, penalty_params=penalty_params,
                              allow_unstab=allow_unstab, minval=matrix(0, nrow=n_obs-p, ncol=M))
    tweights_ok <- !any(vapply(1:M, function(m) sum(tweights[,m] > red_criteria[1]) < red_criteria[2]*n_obs, logical(1)))
    if(Omegas_ok && stat_ok && tweights_ok && weightpars_ok) {
      which_best_fit <- which_round # The estimation round of the appropriate estimate with the largest loglik
      if(!no_prints) message(paste0("Filtered through ", i1-1, " inappropriate ", ifelse(allow_unstab, "(or unstable)", ""),
                                    " estimates with a larger", ifelse(penalized, " (penalized) ", " "), "log-likelihood"))
      break
    }
    if(i1 == length(all_estimates)) {
      message(paste0("No 'appropriate'", ifelse(allow_unstab, "(or unstable)", ""),
                     "estimates found! Check that all the variables are scaled to vary in similar magninutes, ",
                     "also not very small or large magnitudes."))
      message("Consider running more estimation rounds or study the obtained estimates one-by-one with the function alt_stvar.")
      if(allow_unstab) {
        message("Unstable estimates were allowed in the estimation. Consider also using penalized estimation with a higher value
                 for the tuning parameter in 'penalty_params'.")
      }
      if(M > 2) {
        message("Consider also using smaller M. Too large M leads to identification problems.")
      }

      which_best_fit <- which(loks == max(loks))[1]
    }
  }
  params <- all_estimates[[which_best_fit]]


  ## Obtain standard errors, calculate IC ###
  # Sort regimes if no constraints are employed (affects params only with specific weight_functions)
  if(is.null(AR_constraints) && is.null(mean_constraints) && is.null(weight_constraints)) {
    params <- sort_regimes(p=p, M=M, d=d, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
                           cond_dist=cond_dist, identification="reduced_form")
    all_estimates <- lapply(all_estimates, function(pars) sort_regimes(p=p, M=M, d=d, params=pars,
                                                                       weight_function=weight_function,
                                                                       weightfun_pars=weightfun_pars,
                                                                       cond_dist=cond_dist,
                                                                       identification="reduced_form"))
  }

  # Sort and sign change the columns of the impact matrices if cond_dist == "ind_Student" or "ind_skewed_t
  if(cond_dist == "ind_Student" || cond_dist == "ind_skewed_t") {
    params <- sort_impactmats(p=p, M=M, d=d, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
                              cond_dist=cond_dist, AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                              weight_constraints=weight_constraints)
    all_estimates <- lapply(all_estimates, function(pars) {
      sort_impactmats(p=p, M=M, d=d, params=pars,
                      weight_function=weight_function,
                      weightfun_pars=weightfun_pars,
                      cond_dist=cond_dist,
                      AR_constraints=AR_constraints,
                      mean_constraints=mean_constraints,
                      weight_constraints=weight_constraints)
    })
  }

  if(!converged[which_best_fit]) {
    if(!no_prints) message(paste("Iteration limit was reached when estimating the best fitting individual!",
                                 "Consider further estimation with the function 'iterate_more'"))
  }

  transition_weights <- loglikelihood(data=data, p=p, M=M, params=params,
                                      weight_function=weight_function, weightfun_pars=weightfun_pars,
                                      cond_dist=cond_dist, parametrization=parametrization,
                                      identification="reduced_form", AR_constraints=AR_constraints,
                                      mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                      B_constraints=NULL, to_return="tw", check_params=TRUE, penalized=penalized,
                                      penalty_params=penalty_params, allow_unstab=allow_unstab, minval=minval)
  if(any(vapply(1:M, function(i1) sum(transition_weights[,i1] > red_criteria[1]) < red_criteria[2]*n_obs, logical(1)))) {
    message("At least one of the regimes in the estimated model seems to be wasted in the best fitting individual!")
  }

  ### Wrap up ###
  ret <- STVAR(data=data, p=p, M=M, d=d, params=params,
               weight_function=weight_function,
               weightfun_pars=weightfun_pars,
               cond_dist=cond_dist,
               parametrization=parametrization,
               identification="reduced_form",
               AR_constraints=AR_constraints,
               mean_constraints=mean_constraints,
               weight_constraints=weight_constraints,
               B_constraints=NULL,
               penalized=penalized,
               penalty_params=penalty_params,
               allow_unstab=allow_unstab,
               calc_std_errors=calc_std_errors)
  ret$all_estimates <- all_estimates
  ret$all_logliks <- loks
  ret$which_converged <- converged
  ret$which_round <- which_best_fit # Which estimation round induced the largest log-likelihood?
  if(estim_method == "three-phase") {
    ret$LS_estimates <- LS_res_to_ret # Least squares estimates (mean or intercept parametrization)
  }
  warn_eigens(ret, allow_unstab=allow_unstab)
  if(!no_prints) message("Finished!\n")
  ret
}


