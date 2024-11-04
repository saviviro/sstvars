#' @title Two-phase or multiple-phase maximum likelihood estimation of a reduced form smooth transition VAR model
#'
#' @description \code{fitSTVAR} estimates a reduced form smooth transition VAR model in two phases
#'   or multiple phases (see details-section for the multiple phase estimation). In the two-phase
#'   procedure:
#'   in the first phase, it uses a genetic algorithm (GA) to find starting values for a gradient based
#'   variable metric algorithm (VA), which it then uses to finalize the estimation in the second phase.
#'   Parallel computing is utilized to perform multiple rounds of estimations in parallel.
#'
#' @inheritParams GAfit
#' @param estim_method either \code{"two-phase"} or \code{"multiple-phase"} (the latter is the default
#'   option for threshold models and the former is currently the only option for other models). See details.
#' @param nrounds the number of estimation rounds that should be performed. The default is \code{(M*ncol(data))^3}
#'   when \code{estim_method="two-phase"} and \code{(M*ncol(data))^2} when \code{estim_method="multiple-phase"}.
#' @param ncores the number CPU cores to be used in parallel computing.
#' @param maxit the maximum number of iterations in the variable metric algorithm.
#' @param seeds a length \code{nrounds} vector containing the random number generator seed
#'  for each call to the genetic algorithm, or \code{NULL} for not initializing the seed.
#' @param print_res should summaries of estimation results be printed?
#' @param use_parallel employ parallel computing? If \code{use_parallel=FALSE && print_res=FALSE},
#'  nothing is printed during the estimation process.
#' @param filter_estimates should the likely inappropriate estimates be filtered? See details.
#' @param ... additional settings passed to the function \code{GAfit} employing the genetic algorithm.
#' @details
#'  If you wish to estimate a structural model, estimate first the reduced form model and then use the
#'  use the function \code{fitSSTVAR} to create (and estimate if necessary) the structural model
#'  based on the estimated reduced form model.
#'
#'  \strong{Multiple-phase estimation:}\\
#'  With \code{estim_method="multiple-phase"} (currently only available for threshold VAR models),
#'  the following multiple-phase procedure proposed by Koivisto, Luoto, and Virolainen (2025) is employed:
#'  \enumerate{
#'    \item The AR and weight function parameters are estimated by the method of least squares.
#'    \item The rest of the parameters are ML estimated conditional on the LS estimates of the AR parameters with
#'      the two-phase procedure (GA + VA).
#'    \item The AR parameters are ML estimated with VA conditional on the estimates of the rest of the parameters,
#'      initializing the algorithm from the LS estimates.
#'    \item All parameters are ML estimated by initializing VA from the estimates obtained in the previous step.
#'  }
#'  Note that \code{mean_constraints} are not supported in the multiple-phase procedure and only such \code{weight_constraints}
#'  are supported that fix the values of the weight function parameters to known constants. Also, despite running multiple
#'  estimation rounds in Phase~2, the multiple-phase estimation procedure produces only one final estimate, since the
#'  best appropriate estimate is automatically selected after Phase~2 (see "Filtering inappropriate estimates" below).
#'
#'  If structural model identified by non-Gaussianity is estimated, note that the identification can be weak with
#'  respect to the ordering and signs of some of the columns of \eqn{B_2,...,B_M} (see Virolainen, 2024). You can
#'  estimate models with different orderings of the columns of \eqn{B_2,...,B_M} with the function FILL IN.
#'
#'  \strong{Related also to the two-phase estimation:}\\
#'  Because of complexity and high multimodality of the log-likelihood function, it is \strong{not certain}
#'  that the estimation algorithm will end up in the global maximum point. When \code{estim_method="two-phase"},
#'  it is expected that many of the estimation rounds will end up in some local maximum or a saddle point instead.
#'  Therefore, a (sometimes very large) number of estimation rounds is required for reliable results
#'  (when \code{estim_method="multiple-phase"} substantially smaller number should be sufficient). Due to
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
#'  \strong{Filtering inappropriate estimates:} If \code{filter_estimates == TRUE}, the code will automatically filter
#'  through estimates that it deems "inappropriate". That is, estimates that are not likely solutions of interest.
#'  Specifically, solutions that incorporate a near-singular error term covariance matrix (any eigenvalue less than \eqn{0.002}),
#'  any modulus of the eigenvalues of the companion form AR matrices larger than $0.9985$ (indicating the necessary condition for
#'  stationarity is close to break), or transition weights such that they are close to zero for almost all \eqn{t} for at least
#'  one regime. You can also set \code{filter_estimates=FALSE} and find the solutions of interest yourself by using the
#'  function \code{alt_stvar} (which can used with \code{filter_estimates=TRUE} as well since results from all estimation rounds
#'  are saved). If \code{estim_method = "multiple-phase"} filtering is always automatically applied after Phase~2 to facilitate
#'  reasonable initial estimates in the subsequent phases.
#'
#' @inherit STVAR return
#' @section S3 methods:
#'   The following S3 methods are supported for class \code{'stvar'}: \code{logLik}, \code{residuals}, \code{print}, \code{summary},
#'    \code{predict}, \code{simulate}, and \code{plot}.
#' @seealso \code{\link{fitSSTVAR}}, \code{\link{STVAR}}, \code{\link{GAfit}}, \code{\link{iterate_more}}
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
#' # the random number generator seeds set for reproducibility.
#' fit32 <- fitSTVAR(gdpdef, p=3, M=2, weight_function="relative_dens",
#'  cond_dist="Gaussian", nrounds=2, ncores=2, seeds=1:2)
#'
#' # Examine the results:
#' fit32 # Printout of the estimates
#' summary(fit32) # A more detailed summary printout
#' plot(fit32) # Plot the fitted transition weights
#' get_foc(fit32) # Gradient of the log-likelihood function about the estimate
#' get_soc(fit32) # Eigenvalues of the Hessian of the log-lik. fn. about the estimate
#' profile_logliks(fit32) # Profile log-likelihood functions about the estimate
#'
#' # Estimate a two-regime Student's t STVAR p=3 model with logistic transition weights
#' # and the first lag of the second variable as the switching variable, only two
#' # estimation rounds using two CPU cores:
#' fitlogistict32 <- fitSTVAR(gdpdef, p=3, M=2, weight_function="logistic", weightfun_pars=c(2, 1),
#'  cond_dist="Student", nrounds=2, ncores=2, seeds=1:2)
#' summary(fitlogistict32) # Summary printout of the estimates
#'
#' # Estimate a two-regime threshold VAR p=3 model with independent skewed t shocks.
#' # The first lag of the the second variable is specified as the switching variable,
#' # and the threshold parameter constrained to the fixed value 1.
#' fitthres32wit <- fitSTVAR(gdpdef, p=3, M=2, weight_function="threshold", weightfun_pars=c(2, 1),
#'   cond_dist="ind_skewed_t", weight_constraints=list(R=0, r=1), nrounds=2, ncores=2, seeds=1:2)
#' plot(fitthres32wit) # Plot the fitted transition weights
#'
#' # Estimate a two-regime STVAR p=1 model with exogenous transition weights defined as the indicator
#' # of NBER based U.S. recessions (source: St. Louis Fed database). Moreover, restrict the AR matrices
#' # to be identical across the regimes (i.e., allowing for time-variation in the intercepts and the
#' # covariance matrix only):
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
#'   cond_dist="ind_Student", AR_constraints=C_122, nrounds=2, ncores=2, seeds=1:2)
#' plot(fitexo12cit) # Plot the transition weights
#' summary(fitexo12cit) # Summary printout of the estimates
#'
#' # Estimate a two-regime Gaussian STVAR p=1 model with the weighted relative stationary densities
#' # of the regimes as the transition weight function, and the means of the regimes
#' # and AR matrices constrained to be identical across the regimes (i.e., allowing for time-varying
#' # conditional covariance matrix only):
#' fit12cm <- fitSTVAR(gdpdef, p=1, M=2, weight_function="relative_dens", cond_dist="Gaussian",
#'  AR_constraints=C_122, mean_constraints=list(1:2), parametrization="mean", nrounds=2, seeds=1:2)
#' fit12cm # Print the estimates
#'
#' # Estimate a two-regime Gaussian STVAR p=1 model with the weighted relative stationary densities
#' # of the regimes as the transition weight function; constrain AR matrices to be identical
#' # across the regimes and also constrain the off-diagonal elements of the AR matrices to be zero.
#' mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE)
#' C_222 <- rbind(mat0, mat0) # The constraint matrix
#' fit22c <- fitSTVAR(gdpdef, p=2, M=2, weight_function="relative_dens", cond_dist="Gaussian",
#'  AR_constraints=C_222, nrounds=2, seeds=1:2)
#' fit22c # Print the estimates
#'
#' # Estimate a two-regime Student's t STVAR p=3 model with logistic transition weights
#' # and the first lag of the second variable as the switching variable. Constraint the location
#' # parameter to the fixed value 1 and leave the scale parameter unconstrained.
#' fitlogistic32w <- fitSTVAR(gdpdef, p=3, M=2, weight_function="logistic", weightfun_pars=c(2, 1),
#'  weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(1, 0)), nrounds=2, seeds=1:2)
#' plot(fitlogistic32w) # Plot the fitted transition weights
#'
#' # Estimate a two-regime Gaussian STVAR p=3 model with multinomial logit transition weights
#' # using the second variable is the switching variable with two lags. Constrain the AR matrices
#' # identical across the regimes (allowing for time-variation in the intercepts and covariance
#' # matrix).
#' C_322 <- rbind(diag(3*2^2), diag(3*2^2)) # The constraint matrix
#' fitmlogit32c <- fitSTVAR(gdpdef, p=3, M=2, weight_function="mlogit", cond_dist="Gaussian",
#'  weightfun_pars=list(vars=2, lags=2), AR_constraints=C_322, nrounds=1, seeds=3, ncores=1)
#' plot(fitmlogit32c) # Plot the fitted transition weights
#'
#' # Estimate a two-regime Gaussian STVAR p=3 model with exponential transition weights and the first
#' # lag of the second variable as switching variable, and AR parameter constrained identical across
#' # the regimes, means constrained identical across the regimes, and the location parameter
#' # constrained to 0.5 (but scale parameter unconstrained).
#' fitexp32cmw <- fitSTVAR(gdpdef, p=3, M=2, weight_function="exponential", weightfun_pars=c(2, 1),
#'  cond_dist="Student", AR_constraints=C_322, mean_constraints=list(1:2),
#'  weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(0.5, 0)), nrounds=1, seeds=1, ncores=1)
#' summary(fitexp32cmw) # Summary printout of the estimates
#' }
#' @export

fitSTVAR <- function(data, p, M, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                     weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t"),
                     parametrization=c("intercept", "mean"), AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL,
                     estim_method, nrounds=(M + 1)^5, ncores=2, maxit=1000, seeds=NULL, print_res=TRUE, use_parallel=TRUE,
                     filter_estimates=TRUE, ...) {
  # Initial checks etc
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  if(missing(estim_method)) { # Determine the estimation method
    if(weight_function == "threshold" && is.null(mean_constraints)) {
      if(is.null(weight_constraints)) {
        estim_method <- "multiple-phase"
      } else {
        if(weight_constraints$R != 0) {
          estim_method <- "two-phase"
        } else {
          estim_method <- "multiple-phase"
        }
      }
    } else {
      estim_method <- "two-phase"
    }
  } else {
    stopifnot(estim_method %in% c("two-phase", "multiple-phase"))
    if(estim_method == "multiple-phase" && weight_function != "threshold") {
      stop("The multiple-phase estimation is currently only available for threshold models.")
    } else if(estim_method == "multiple-phase" && !is.null(mean_constraints)) {
      stop("The multiple-phase estimation does not support mean_constraints.")
    } else if(estim_method == "multiple-phase" && !is.null(weight_constraints) && weight_constraints$R != 0) {
      stop("The multiple-phase estimation only supports weight_constraints with R=0.")
    }
  }
  if(missing(nrounds)) {
    if(estim_method == "two-phase") {
      nrounds <- (M*ncol(data))^3
    } else { # Multiple-phase estimation, GA estimation of distribution params only
      nrounds <- (M*ncol(data))^2
    }
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
  if(npars >= d*nrow(data)) stop(paste("There are at least as many parameters in the model as there are observations in the data,",
                                       "so you should use smaller p or M."))
  dot_params <- list(...)
  minval <- ifelse(is.null(dot_params$minval), get_minval(data), dot_params$minval)
  red_criteria <- ifelse(rep(is.null(dot_params$red_criteria), 2), c(0.05, 0.01), dot_params$red_criteria)

  # Parallel computing checks
  if(ncores > parallel::detectCores()) {
    ncores <- parallel::detectCores()
    message("ncores was set to be larger than the number of cores detected.")
  }
  if(ncores > nrounds) {
    ncores <- nrounds
    message("ncores was set to be larger than the number of estimation rounds.")
  }

  ## A function to calculate the gradient of the log-likelihood function for the optimization algorithms
  ## using central difference approximation:
  h <- 6e-6 # The difference used in the gradient
  loglik_grad <- function(params, FUN, number_of_pars) {
    # FUN = function to calculate the log-likelihood, with the parameter vector as the only argument
    # number_of_pars = number of parameters in the parameter vector that is the only argument of FUN
    I <- diag(rep(1, number_of_pars))
    vapply(1:npars, function(i1) (FUN(params + I[i1,]*h) - FUN(params - I[i1,]*h))/(2*h), numeric(1))
  }

  ## Function to filter inappropriate estimates
  filter_estimates_fun <- function(all_estimates, loks) {
    # all_estimates = list of all estimates, loks = vector of the corresponding log-likelihoods
    # Returns the estimation round number with the best "appropriate" estimate
    if(!no_prints) message("Filtering inappropriate estimates...")
    ord_by_loks <- order(loks, decreasing=TRUE) # Ordering from largest loglik to smaller

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
                                to_return="tw", check_params=TRUE, minval=matrix(0, nrow=n_obs-p, ncol=M))
      tweights_ok <- !any(vapply(1:M, function(m) sum(tweights[,m] > red_criteria[1]) < red_criteria[2]*n_obs, logical(1)))
      if(Omegas_ok && stat_ok && tweights_ok && weightpars_ok) {
        which_best_fit <- which_round # The estimation round of the appropriate estimate with the largest loglik
        if(!no_prints && i1 != 1) message(paste("Filtered out", i1-1, "inappropriate estimates with a larger log-likelihood"))
        break
      }
      if(i1 == length(all_estimates)) {
        message(paste("No 'appropriate' estimates were found! Check that all the variables are scaled to vary in similar magninutes,",
                      "also not very small or large magnitudes."))
        if(estim_method == "two-phase") {
          message("Consider running more estimation rounds or study the obtained estimates one-by-one with the function alt_stvar.")
        }
        if(M > 2) {
          message("Consider also using smaller M. Too large M leads to identification problems.")
        }

        which_best_fit <- which(loks == max(loks))[1]
      }
    }
    which_best_fit # return the estimation round with the best "appropriate" estimate
  }

  ####################################################################
  ####################### TWO-PHASE ESTIMATION #######################
  ####################################################################

  if(estim_method == "two-phase") {
    if(use_parallel) {
      message(paste("Using", ncores, "cores for", nrounds, "estimations rounds..."))

      ### Optimization with the genetic algorithm ###
      cl <- parallel::makeCluster(ncores)
      on.exit(try(parallel::stopCluster(cl), silent=TRUE)) # Close the cluster on exit, if not already closed.
      parallel::clusterExport(cl, ls(environment(fitSTVAR)), envir=environment(fitSTVAR)) # assign all variables from package:sstvars
      parallel::clusterEvalQ(cl, c(library(pbapply), library(Rcpp), library(RcppArmadillo), library(sstvars)))

      message("Optimizing with a genetic algorithm...")
      GAresults <- pbapply::pblapply(1:nrounds, function(i1) GAfit(data=data, p=p, M=M,
                                                                   weight_function=weight_function,
                                                                   weightfun_pars=weightfun_pars,
                                                                   cond_dist=cond_dist,
                                                                   parametrization=parametrization,
                                                                   AR_constraints=AR_constraints,
                                                                   mean_constraints=mean_constraints,
                                                                   weight_constraints=weight_constraints,
                                                                   seed=seeds[i1], ...), cl=cl)

    } else {
      if(!no_prints) message("Optimizing with a genetic algorithm...")

      tmpfunGA <- function(i1, ...) {
        if(!no_prints) message(i1, "/", nrounds, "\r")
        GAfit(data=data, p=p, M=M,
              weight_function=weight_function,
              weightfun_pars=weightfun_pars,
              cond_dist=cond_dist,
              parametrization=parametrization,
              AR_constraints=AR_constraints,
              mean_constraints=mean_constraints,
              weight_constraints=weight_constraints,
              seed=seeds[i1], ...)
      }
      GAresults <- lapply(1:nrounds, function(i1) tmpfunGA(i1, ...))
    }

    loks <- vapply(1:nrounds, function(i1) loglik_fn(params=GAresults[[i1]]), numeric(1))

    if(print_res) {
      print_loks <- function(loks) {
        printfun <- function(txt, FUN) message(paste(txt, round(FUN(loks), 3)))
        printfun("The lowest loglik: ", min)
        printfun("The largest loglik:", max)
      }
      message("Results from the genetic algorithm:")
      print_loks(loks)
    }

    ### Optimization with the variable metric algorithm ###

    # A function to calculate the log-likelihood function
    loglik_fn <- function(params) {
      tryCatch(loglikelihood(data=data, p=p, M=M, params=params,
                             weight_function=weight_function, weightfun_pars=weightfun_pars,
                             cond_dist=cond_dist, parametrization=parametrization,
                             identification="reduced_form", AR_constraints=AR_constraints,
                             mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                             B_constraints=NULL, to_return="loglik", check_params=TRUE, minval=minval,
                             alt_par=TRUE), # alt_par used in the GA estimation
               error=function(e) minval)
    }

    if(!no_prints) message("Optimizing with a variable metric algorithm...")
    if(use_parallel) {
      NEWTONresults <- pbapply::pblapply(1:nrounds, function(i1) optim(par=GAresults[[i1]], fn=loglik_fn, gr=loglik_grad,
                                                                       FUN=loglik_fn, number_of_pars=npars, # Passed to loglik_grad
                                                                       method="BFGS", control=list(fnscale=-1, maxit=maxit)), cl=cl)
      parallel::stopCluster(cl=cl)
    } else {
      tmpfunNE <- function(i1) {
        if(!no_prints) message(i1, "/", nrounds, "\r")
        optim(par=GAresults[[i1]], fn=loglik_fn, gr=loglik_grad,
              FUN=loglik_fn, number_of_pars=npars, # Passed to loglik_grad
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

    if(filter_estimates) {
      which_best_fit <- filter_estimates_fun(all_estimates=all_estimates, loks=loks)
    } else {
      which_best_fit <- which(loks == max(loks))[1]
    }
    params <- all_estimates[[which_best_fit]] # The params to return
  } else {
    #########################################################################
    ####################### MULTIPLE-PHASE ESTIMATION #######################
    #########################################################################
    if(use_parallel) {
      message(paste("Using", ncores, "cores for estimation..."))
    }

    ### Phase 1: Estimate AR and weight parameters by least squares (always estimates with intercept parametrization)
    if(!no_prints && !use_parallel) {
      message("PHASE 1: Estimating AR and weight parameters by least squares...") # parallel prints inside estim_LS
    }
    LS_results <- estim_LS(data=data, p=p, M=M, weight_function=weight_function, weightfun_pars=weightfun_pars,
                           cond_dist=cond_dist, parametrization=parametrization, AR_constraints=AR_constraints,
                           mean_constraints=mean_constraints, weight_constraints=weight_constraints, ncores=ncores,
                           use_parallel=use_parallel)

    ### Phase 2: Estimate the remaining parameters conditional on the LS estimates of the AR and weight parameters

    ## Part 1: Estimation by genetic algorithm ##
    if(!no_prints) message(paste0("PHASE 2a: Estimating error distribution parameters with a genetic algorithm (",
                                  nrounds, " rounds)..."))
    if(use_parallel) {
      cl <- parallel::makeCluster(ncores)
      on.exit(try(parallel::stopCluster(cl), silent=TRUE)) # Close the cluster on exit, if not already closed.
      parallel::clusterExport(cl, ls(environment(fitSTVAR)), envir=environment(fitSTVAR)) # assign all variables from package:sstvars
      parallel::clusterEvalQ(cl, c(library(pbapply), library(Rcpp), library(RcppArmadillo), library(sstvars)))

      GAresults <- pbapply::pblapply(1:nrounds, function(i1) GAfit(data=data, p=p, M=M,
                                                                   weight_function=weight_function,
                                                                   weightfun_pars=weightfun_pars,
                                                                   cond_dist=cond_dist,
                                                                   parametrization="intercept", # LS_results always intercept paramtrz
                                                                   AR_constraints=AR_constraints,
                                                                   mean_constraints=mean_constraints,
                                                                   weight_constraints=weight_constraints,
                                                                   fixed_params=LS_results,
                                                                   seed=seeds[i1], ...), cl=cl)
    } else {
      tmpfunGA <- function(i1, ...) {
        if(!no_prints) message(i1, "/", nrounds, "\r")
        GAfit(data=data, p=p, M=M,
              weight_function=weight_function,
              weightfun_pars=weightfun_pars,
              cond_dist=cond_dist,
              parametrization="intercept", # LS_results always intercept parametrized params
              AR_constraints=AR_constraints,
              mean_constraints=mean_constraints,
              weight_constraints=weight_constraints,
              fixed_params=LS_results,
              seed=seeds[i1], ...)
      }
      GAresults <- lapply(1:nrounds, function(i1) tmpfunGA(i1, ...))
    }

    ## Change to mean parametrization if parametrization == "mean"
    if(parametrization == "mean") {
      GAresults <- lapply(GAresults, function(pars) change_parametrization(p=p, M=M, d=d, params=pars,
                                                                           weight_function=weight_function,
                                                                           weightfun_pars=weightfun_pars,
                                                                           cond_dist=cond_dist,
                                                                           identification="reduced_form",
                                                                           AR_constraints=AR_constraints,
                                                                           mean_constraints=mean_constraints,
                                                                           weight_constraints=weight_constraints,
                                                                           B_constraints=NULL,
                                                                           change_to="mean"))
    }

    # Print results from GA estimation
    if(print_res) {
      loks <- vapply(1:nrounds, function(i1) loglikelihood(data=data, p=p, M=M, params=GAresults[[i1]],
                                                           weight_function=weight_function, weightfun_pars=weightfun_pars,
                                                           cond_dist=cond_dist, parametrization=parametrization,
                                                           identification="reduced_form", AR_constraints=AR_constraints,
                                                           mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                                           B_constraints=NULL, to_return="loglik", check_params=TRUE, minval=minval,
                                                           alt_par=TRUE), numeric(1)) # GA results always in alt_par
      print_loks <- function(loks) {
        printfun <- function(txt, FUN) message(paste(txt, round(FUN(loks), 3)))
        printfun("The lowest loglik:  ", min)
        printfun("The mean loglik:    ", mean)
        printfun("The largest loglik: ", max)
      }
      message("Results from the genetic algorithm:")
      print_loks(loks)
    }

    ## Part 2: Estimation by variable metric algorithm ##
    n_covmatpars <- ifelse(cond_dist %in% c("ind_Student", "ind_skewed_t"), M*d^2, M*d*(d+1)/2)
    n_arpars <- ifelse(is.null(AR_constraints), M*d + M*p*d^2, M*d + ncol(AR_constraints))
    weight_pars_are_fixed <- ifelse(!is.null(weight_constraints) && weight_constraints[[1]] == 0, TRUE, FALSE)
    if(weight_function == "exogenous" || M == 1 || weight_pars_are_fixed) {
      n_weightpars <- 0
    } else if(weight_function %in% c("threshold", "relative_dens")) {
      n_weightpars <- M - 1
    } else if(weight_function %in% c("logistic", "exponential")) {
      n_weightpars <- 2
    } else if(weight_function == "mlogit") {
      n_weightpars <- (M - 1)*(1 + length(weightfun_pars[[1]])*weightfun_pars[[2]])
    }
    if(cond_dist == "Gaussian") {
      n_distpars <- 0
    } else if(cond_dist == "Student") {
      n_distpars <- 1
    } else if(cond_dist == "ind_Student") {
      n_distpars <- d
    } else if(cond_dist == "ind_skewed_t") {
      n_distpars <- 2*d
    }

    # A function to calculate the log-likelihood function
    loglik_fn2 <- function(covmat_and_dist_pars) {
      if(n_weightpars == 0) {
        params <- c(LS_results[1:n_arpars], covmat_and_dist_pars)
      } else {
        params <- c(LS_results[1:n_arpars], covmat_and_dist_pars[1:n_covmatpars],
                    LS_results[(n_arpars + 1):(n_arpars + n_weightpars)],
                    covmat_and_dist_pars[(n_covmatpars + 1):(n_covmatpars + n_distpars)])
      }
      tryCatch(loglikelihood(data=data, p=p, M=M, params=params,
                             weight_function=weight_function, weightfun_pars=weightfun_pars,
                             cond_dist=cond_dist, parametrization=parametrization,
                             identification="reduced_form", AR_constraints=AR_constraints,
                             mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                             B_constraints=NULL, to_return="loglik", check_params=TRUE, minval=minval,
                             alt_par=TRUE), # alt_par params procuded by the GA estimation
               error=function(e) minval)
    }

    # Function that picks covmat and dist params from the full parameter vector
    get_covmat_and_dist_pars <- function(params) {
      if(n_weightpars == 0) {
        return(params[(n_arpars + 1):(n_arpars + n_covmatpars + n_distpars)])
      } else {
        return(c(params[(n_arpars + 1):(n_arpars + n_covmatpars)],
                 params[(n_arpars + n_covmatpars + n_weight_pars + 1):(n_arpars + n_covmatpars + n_weightpars + n_distpars)]))
      }
    }

    if(!no_prints) message(paste0("PHASE 2b: Estimating error distribution parameters with a variable algorithm (",
                                  nrounds, " rounds)..."))
    if(use_parallel) {
      NEWTONresults <- pbapply::pblapply(1:nrounds,
                                         function(i1) optim(par=get_covmat_and_dist_pars(GAresults[[i1]]),
                                                            fn=loglik_fn2, gr=loglik_grad,
                                                            FUN=loglik_fn2, number_of_pars=n_covmatpars+n_distpars, # passed to loglik_grad
                                                            method="BFGS", control=list(fnscale=-1, maxit=maxit)), cl=cl)
      parallel::stopCluster(cl=cl)
    } else {
      tmpfunNE <- function(i1) {
        if(!no_prints) message(i1, "/", nrounds, "\r")
        optim(par=get_covmat_and_dist_pars(GAresults[[i1]]),
              fn=loglik_fn2, gr=loglik_grad,
              FUN=loglik_fn2, number_of_pars=n_covmatpars+n_distpars, # passed to loglik_grad
              method="BFGS", control=list(fnscale=-1, maxit=maxit))
      }
      NEWTONresults <- lapply(1:nrounds, function(i1) tmpfunNE(i1))
    }

    loks <- vapply(1:nrounds, function(i1) NEWTONresults[[i1]]$value, numeric(1)) # Log-likelihoods

    if(print_res) {
      message("Results from the variable metric algorithm:")
      print_loks()
    }

    ## Obtain estimates, change back to original parametrization, and filter the inapproriate estimates
    all_estimates <- lapply(NEWTONresults, function(x) x$par)

    ## TÄÄLLÄ PITÄÄ RAKENTAA TAKAISIN PARAMETRIVEKTORI KOKONAISEKSI KAIKILLE YKSILÖILLE

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

    ## Filter inappropriate estimates
    which_best_fit <- filter_estimates_fun(all_estimates=all_estimates, loks=loks)
    phase2_estimate <- all_estimates[[which_best_fit]] # The params to initialize the next phase with

    ### Phase 3: estimate the AR and weight parameters by ML with VA with the distribution parameters fixed
    fixed_covmatpars <- phase2_estimate[(n_arpars + 1):(n_arpars + n_covmatpars)]
    fixed_distpars <- phase2_estimate[(n_arpars + n_covmatpars + n_weightpars + 1):
                                        (n_arpars + n_covmatpars + n_weightpars + n_distpars)]
    if(n_weightpars == 0) {
      phase2_ar_and_weight_estim <- phase2_estimate[1:n_arpars]
    } else {
      phase2_ar_and_weight_estim <- c(phase2_estimate[1:n_arpars],
                                      phase2_estimate[(n_arpars + n_covmatpars + 1):(n_arpars + n_covmatpars + n_weightpars)])
    }

    ## Log-lik function and gradient with fixed covmat and distpars:
    loglik_fn2 <- function(ar_and_weight_pars) {
      # ar and weight_pars assumed to be in the form (phi_{1,0},...,phi_{M,9},varphi_{1},...,\varphi_{M},alpha)

      ### TÄÄLLÄ HUOMIOI JÄLLEEN WEIGHTPARSSEJA OLLENKAAN!
      params <- c(ar_and_weight_pars[1:n_arpars], fixed_covmatpars, ar_and_weight_pars[(n_arpars + 1):(n_arpars + n_weightpars)],
                  fixed_distpars) # Create the full parameter vector with the fixed covmat and dist pars
      tryCatch(loglikelihood(data=data, p=p, M=M, params=params,
                             weight_function=weight_function, weightfun_pars=weightfun_pars,
                             cond_dist=cond_dist, parametrization=parametrization,
                             identification="reduced_form", AR_constraints=AR_constraints,
                             mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                             B_constraints=NULL, to_return="loglik", check_params=TRUE, minval=minval,
                             alt_par=FALSE), # We have swapped to orig parametrization here.
               error=function(e) minval)
    }
    n_ar_and_weightpars <- n_arpars + n_weightpars
    # I2 <- diag(rep(1, n_ar_and_weightpars))
    # loglik_grad2 <- function(ar_and_weight_pars) {
    #   # ar and weight_pars assumed to be inthe form (phi_{1,0},...,phi_{M,9},varphi_{1},...,\varphi_{M},alpha)
    #   #params <- c(ar_and_weight_pars[1:n_arpars], fixed_covmatpars, ar_and_weight_pars[(n_arpars + 1):(n_arpars + n_weightpars)],
    #   #            fixed_distpars) # Create the full parameter vector with the fixed covmat and dist pars
    #   vapply(1:n_ar_and_weightpars, function(i1) {
    #     # Determine the index of the i1:th param in ar_and_weight_pars in the full parameter vector
    #     #if(i1 <= n_arpars) { # AR or mean param
    #     #  i2 <- i1
    #     #} else { # Weight param
    #     #  i2 <- i1 + n_covmatpars
    #     #}
    #     #(loglik_fn(params + I[i2,]*h) - loglik_fn(params - I[i2,]*h))/(2*h)
    #     (loglik_fn2(ar_and_weight_pars + I2[i1,]*h) - loglik_fn(ar_and_weight_pars - I2[i1,]*h))/(2*h)
    #     }, numeric(1))
    # }

    if(!no_prints) message("PHASE 3: Estimating AR and weight parameters with a variable metric algorithm...")
    phase3_res <- optim(par=phase2_ar_and_weight_estim, fn=loglik_fn2, gr=loglik_grad2,
                        method="BFGS", control=list(fnscale=-1, maxit=maxit))
    phase3_ar_estimate <- phase3_res$par
    phase3_loglik <- phase3_res$value
    if(print_res) message(paste("The log-likelihood from PHASE 3:", round(phase3_res$value, 3)))
    if(n_weightpars == 0) { # No weight pars
      phase3_estimate <- c(phase3_ar_estimate[1:n_arpars], fixed_covmatpars, fixed_distpars) # Full Phase 3 estimate
    } else {
      phase3_estimate <- c(phase3_ar_estimate[1:n_arpars], fixed_covmatpars,
                           phase3_ar_estimate[(n_arpars + 1):(n_arpars + n_weightpars)],
                           fixed_distpars) # Full Phase 3 estimate
    }
    print(loglikelihood(data=data, p=p, M=M, params=phase3_estimate,
                        weight_function=weight_function, weightfun_pars=weightfun_pars,
                        cond_dist=cond_dist, parametrization=parametrization,
                        identification="reduced_form", AR_constraints=AR_constraints,
                        mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                        B_constraints=NULL, to_return="loglik", check_params=TRUE, minval=minval))

    ### Phase 4: estimate all parameters by ML initializing VA from Phase 3 estimates:
    if(!no_prints) message("PHASE 4: Estimating all parameters with a variable metric algorithm...")
    phase4_res <- optim(par=phase3_estimate, fn=loglik_fn, gr=loglik_grad,
                        method="BFGS", control=list(fnscale=-1, maxit=maxit))
    if(print_res) message(paste("The log-likelihood from PHASE 4:", round(phase4_res$value, 3)))
    print(loglikelihood(data=data, p=p, M=M, params=phase4_res$par,
                        weight_function=weight_function, weightfun_pars=weightfun_pars,
                        cond_dist=cond_dist, parametrization=parametrization,
                        identification="reduced_form", AR_constraints=AR_constraints,
                        mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                        B_constraints=NULL, to_return="loglik", check_params=TRUE, minval=minval))

    ## Obtain the results:
    converged <- phase4_res$convergence == 0 # Did the final result converge?
    which_best_fit <- 1 # Only one estimation round for the final estimate in multiple-phase estimation
    all_estimates <- list(phase4_res$par) # Only one final estimate in multiple-phase estimation
    loks <- phase4_res$value # The final log-likelihood, only one final estimate
    params <- phase4_res$par # The final estimate
  }

  ### Obtain standard errors, calculate IC ###
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
  # Sort and sign change the columns of the impact matrices if cond_dist == "ind_Student"
  if(cond_dist == "ind_Student" || cond_dist == "ind_skewed_t") {
    params <- sort_impactmats(p=p, M=M, d=d, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
                              cond_dist=cond_dist, AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                              weight_constraints=weight_constraints)
    all_estimates <- lapply(all_estimates, function(pars) sort_impactmats(p=p, M=M, d=d, params=pars,
                                                                          weight_function=weight_function,
                                                                          weightfun_pars=weightfun_pars,
                                                                          cond_dist=cond_dist,
                                                                          AR_constraints=AR_constraints,
                                                                          mean_constraints=mean_constraints,
                                                                          weight_constraints=weight_constraints))
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
                                      B_constraints=NULL, to_return="tw", check_params=TRUE, minval=minval)
  if(any(vapply(1:M, function(i1) sum(transition_weights[,i1] > red_criteria[1]) < red_criteria[2]*n_obs, logical(1)))) {
    message("At least one of the regimes in the estimated model seems to be wasted in the best fitting individual!")
  }

  ### Wrap up ###
  if(!no_prints) message("Calculating approximate standard errors...")
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
               calc_std_errors=TRUE)
  ret$all_estimates <- all_estimates
  ret$all_logliks <- loks
  ret$which_converged <- converged
  ret$which_round <- which_best_fit # Which estimation round induced the largest log-likelihood?
  warn_eigens(ret)
  if(!no_prints) message("Finished!\n")
  ret
}


