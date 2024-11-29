#' @title Create a class 'stvar' object defining a reduced form or structural smooth transition VAR model
#'
#' @description \code{STVAR} creates a class \code{'stvar'} object that defines
#'  a reduced form or structural smooth transition VAR model
#'
#' @inheritParams loglikelihood
#' @param data a matrix or class \code{'ts'} object with \code{d>1} columns. Each column is taken to represent
#'  a single times series. \code{NA} values are not supported. Ignore if defining a model without data is desired.
#' @param d number of times series in the system, i.e. \code{ncol(data)}. This can be
#'   used to define STVAR models without data and can be ignored if \code{data} is provided.
#' @param calc_std_errors should approximate standard errors be calculated?
#' @details If data is provided, then also residuals are computed and included in the returned object.
#' @return Returns an S3 object of class \code{'stvar'} defining a smooth transition VAR model. The returned list
#'  contains the following components (some of which may be \code{NULL} depending on the use case):
#'    \item{data}{The input time series data.}
#'    \item{model}{A list describing the model structure.}
#'    \item{params}{The parameters of the model.}
#'    \item{std_errors}{Approximate standard errors of the parameters, if calculated.}
#'    \item{transition_weights}{The transition weights of the model.}
#'    \item{regime_cmeans}{Conditional means of the regimes, if data is provided.}
#'    \item{total_cmeans}{Total conditional means of the model, if data is provided.}
#'    \item{total_ccovs}{Total conditional covariances of the model, if data is provided.}
#'    \item{uncond_moments}{A list of unconditional moments including regime autocovariances, variances, and means.}
#'    \item{residuals_raw}{Raw residuals, if data is provided.}
#'    \item{residuals_std}{Standardized residuals, if data is provided.}
#'    \item{structural_shocks}{Recovered structural shocks, if applicable.}
#'    \item{loglik}{Log-likelihood of the model, if data is provided.}
#'    \item{IC}{The values of the information criteria (AIC, HQIC, BIC) for the model, if data is provided.}
#'    \item{all_estimates}{The parameter estimates from all estimation rounds, if applicable.}
#'    \item{all_logliks}{The log-likelihood of the estimates from all estimation rounds, if applicable.}
#'    \item{which_converged}{Indicators of which estimation rounds converged, if applicable.}
#'    \item{which_round}{Indicators of which round of optimization each estimate belongs to, if applicable.}
#'    \item{LS_estimates}{The least squares estimates of the parameters in the form
#'      \eqn{(\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M,\alpha} (intercepts replaced by unconditional means
#'      if mean parametrization is used), if applicable.}
#' @section About S3 methods:
#'   If data is not provided, only the \code{print} and \code{simulate} methods are available.
#'   If data is provided, then in addition to the ones listed above, \code{predict} method is also available.
#'   See \code{?simulate.stvar} and \code{?predict.stvar} for details about the usage.
#' @seealso \code{\link{fitSTVAR}}, \code{\link{swap_parametrization}}, \code{\link{alt_stvar}}
#' @references
#'  \itemize{
#'    \item Anderson H., Vahid F. 1998. Testing multiple equation systems for common nonlinear components.
#'      \emph{Journal of Econometrics}, \strong{84}:1, 1-36.
#'    \item Hubrich K., Teräsvirta. T. 2013. Thresholds and Smooth Transitions in Vector Autoregressive Models.
#'      \emph{CREATES Research Paper 2013-18, Aarhus University.}
#'    \item Lanne M., Virolainen S. 2024. A Gaussian smooth transition vector autoregressive model:
#'       An application to the macroeconomic effects of severe weather shocks. Unpublished working
#'       paper, available as arXiv:2403.14216.
#'    \item Kheifets I.L., Saikkonen P.J. 2020. Stationarity and ergodicity of Vector STAR models.
#'      \emph{Econometric Reviews}, \strong{39}:4, 407-414.
#'    \item Lütkepohl H., Netšunajev A. 2017. Structural vector autoregressions with smooth transition in variances.
#'      \emph{Journal of Economic Dynamics & Control}, \strong{84}, 43-57.
#'    \item Tsay R. 1998. Testing and Modeling Multivariate Threshold Models.
#'      \emph{Journal of the American Statistical Association}, \strong{93}:443, 1188-1202.
#'    \item Virolainen S. 2024. Identification by non-Gaussianity in structural threshold and
#'       smooth transition vector autoregressive models. Unpublished working
#'       paper, available as arXiv:2404.19707.
#'  }
#' @examples
#' # Below examples use the example data "gdpdef", which is a two-variate quarterly data
#' # of U.S. GDP and GDP implicit price deflator covering the period from 1959Q1 to 2019Q4.
#'
#' # Gaussian STVAR p=1, M=2, model with the weighted relative stationary densities
#' # of the regimes as the transition weight function:
#' theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863,
#'   -0.3124, 0.505251, 0.298483, 0.030096, -0.176925, 0.838898, 0.310863, 0.007512,
#'   0.018244, 0.949533, -0.016941, 0.121403, 0.573269)
#' mod122 <- STVAR(data=gdpdef, p=1, M=2, params=theta_122relg)
#' print(mod122) # Printout of the model
#' summary(mod122) # Summary printout
#' plot(mod122) # Plot the transition weights
#' plot(mod122, plot_type="cond_mean") # Plot one-step conditional means
#'
#' # Logistic Student's t STVAR with p=1, M=2, and the first lag of the second variable
#' # as the switching variable:
#' params12 <- c(0.62906848, 0.14245295, 2.41245785, 0.66719269, 0.3534745, 0.06041779, -0.34909745,
#'   0.61783824, 0.125769, -0.04094521, -0.99122586, 0.63805416, 0.371575, 0.00314754, 0.03440824,
#'   1.29072533, -0.06067807, 0.18737385, 1.21813844, 5.00884263, 7.70111672)
#' fit12 <- STVAR(data=gdpdef, p=1, M=2, params=params12, weight_function="logistic",
#'   weightfun_pars=c(2, 1), cond_dist="Student")
#' summary(fit12) # Summary printout
#' plot(fit12) # Plot the transition weights
#'
#' # Threshold STVAR with p=1, M=2, the first lag of the second variable as switching variable:
#' params12thres <- c(0.5231, 0.1015, 1.9471, 0.3253, 0.3476, 0.0649, -0.035, 0.7513, 0.1651,
#'  -0.029, -0.7947, 0.7925, 0.4233, 5e-04, 0.0439, 1.2332, -0.0402, 0.1481, 1.2036)
#' mod12thres <- STVAR(data=gdpdef, p=1, M=2, params=params12thres, weight_function="threshold",
#'   weightfun_pars=c(2, 1))
#' mod12thres # Printout of the model
#'
#' # Student's t logistic STVAR with p=2, M=2 with the second lag of the second variable
#' # as the switching variable and structural shocks identified by heteroskedasticity;
#' # the model created without data:
#' params22log <- c(0.357, 0.107, 0.356, 0.086, 0.14, 0.035, -0.165, 0.387, 0.452,
#'  0.013, 0.228, 0.336, 0.239, 0.024, -0.021, 0.708, 0.063, 0.027, 0.009, 0.197,
#'   -0.03, 0.24, -0.76, -0.02, 3.36, 0.86, 0.1, 0.2, 7)
#' mod222logtsh <- STVAR(p=2, M=2, d=2, params=params22log, weight_function="logistic",
#'  weightfun_pars=c(2, 2), cond_dist="Student", identification="heteroskedasticity")
#' print(mod222logtsh) # Printout of the model
#'
#' # STVAR p=2, M=2, model with exogenous transition weights and mutually independent
#' # Student's t shocks:
#' set.seed(1); tw1 <- runif(nrow(gdpdef)-2) # Transition weights of Regime 1
#' params22exoit <- c(0.357, 0.107, 0.356, 0.086, 0.14, 0.035, -0.165, 0.387, 0.452,
#'  0.013, 0.228, 0.336, 0.239, 0.024, -0.021, 0.708, 0.063, 0.027, 0.009, 0.197,
#'  -0.1, 0.2, -0.15, 0.13, 0.21, 0.15, 0.11, -0.09, 3, 4)
#' mod222exoit <- STVAR(p=2, M=2, d=2, params=params22exoit, weight_function="exogenous",
#'  weightfun_pars=cbind(tw1, 1-tw1), cond_dist="ind_Student")
#' print(mod222exoit) # Printout of the model
#'
#' # Linear Gaussian VAR(p=1) model:
#' theta_112 <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103,
#'   0.601786, -0.002945, 0.067224)
#' mod112 <- STVAR(data=gdpdef, p=1, M=1, params=theta_112)
#' summary(mod112) # Summary printout
#' @export

STVAR <- function(data, p, M, d, params,
                  weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                  weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t"),
                  parametrization=c("intercept", "mean"),
                  identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity"),
                  AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL,
                  penalized=FALSE, penalty_params=c(0.05, 1), allow_unstab=FALSE, calc_std_errors=FALSE) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  identification <- match.arg(identification)
  if(cond_dist %in% c("ind_Student", "ind_skewed_t") && !(identification %in% c("reduced_form", "non-Gaussianity"))) {
    stop(paste("If cond_dist='ind_Student' or 'ind_skewed_t', identification must be 'reduced_form' or 'non-Gaussianity'",
               "(short-run restrictions can be imposed by specifying the argument 'B_constraints')"))
  } else if(!(cond_dist %in% c("ind_Student", "ind_skewed_t")) && identification == "non-Gaussianity") {
    stop("Identification by 'non-Gaussianity' is not available for models with cond_dist='Gaussian' or 'Student'.")
  }
  stopifnot(is.logical(penalized))
  stopifnot(is.numeric(penalty_params) && length(penalty_params) == 2 && all(penalty_params >= 0) && penalty_params[1] < 1)
  if(weight_function == "relative_dens" && allow_unstab) {
    message("allow_unstab cannot be used with the relative_dens weight function")
    allow_unstab <- FALSE
  }
  if(missing(data) & missing(d)) stop("data or d must be provided")
  if(missing(data) || is.null(data)) {
    data <- NULL
  } else {
    data <- check_data(data=data, p=p)
    if(missing(d)) {
      d <- ncol(data)
    } else if(ncol(data) != d) {
      warning("ncol(data) does not equal d. Using d = ncol(data)")
      d <- ncol(data)
    }
  }
  if(length(M) != 1 && !all_pos_ints(M)) stop("Argument M must be a positive integer")
  if(M == 1 && weight_function %in% c("logistic", "mlogit", "exponential")) {
    # Set to threshold if only regime (we assume two regimes for logistic and exponential weights)
    weight_function <- "threshold"
    weightfun_pars <- c(1, 1) # There have no affect with M == 1
  }
  check_pMd(p=p, M=M, d=d, weight_function=weight_function, identification=identification)
  weightfun_pars <- check_weightfun_pars(data=data, p=p, d=d, M=M, weight_function=weight_function,
                                         weightfun_pars=weightfun_pars, cond_dist=cond_dist)
  check_constraints(data=data, p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    parametrization=parametrization, identification=identification, AR_constraints=AR_constraints,
                    mean_constraints=mean_constraints, weight_constraints=weight_constraints, B_constraints=B_constraints)
  check_params(data=data, p=p, M=M, d=d, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
               cond_dist=cond_dist, parametrization=parametrization, identification=identification, AR_constraints=AR_constraints,
               mean_constraints=mean_constraints, weight_constraints=weight_constraints, B_constraints=B_constraints,
               allow_unstab=allow_unstab)
  npars <- n_params(p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    cond_dist=cond_dist, identification=identification, AR_constraints=AR_constraints,
                    mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                    B_constraints=B_constraints)

  # Log-likelihood, transition weights, residuals, IC
  if(is.null(data)) {
    lok_and_tw <- list(loglik=NA, mw=NA)
    IC <- data.frame(AIC=NA, HQIC=NA, BIC=NA)
    residuals_raw <- residuals_std <- NA
    structural_shocks <- NA
  } else {
    if(npars >= d*nrow(data)) warning("There are at least as many parameters in the model as there are observations in the data")
    lok_and_tw <- loglikelihood(data=data, p=p, M=M, params=params,
                                weight_function=weight_function, weightfun_pars=weightfun_pars,
                                cond_dist=cond_dist, parametrization=parametrization,
                                identification=identification, AR_constraints=AR_constraints,
                                mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                B_constraints=B_constraints, to_return="loglik_and_tw",
                                penalized=penalized, penalty_params=penalty_params,
                                allow_unstab=allow_unstab, check_params=TRUE, minval=NA)
    residuals_raw <- get_residuals(data=data, p=p, M=M, params=params,
                                   weight_function=weight_function, weightfun_pars=weightfun_pars,
                                   cond_dist=cond_dist, parametrization=parametrization,
                                   identification=identification, AR_constraints=AR_constraints,
                                   mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                   B_constraints=B_constraints, standardize=FALSE,
                                   penalized=penalized, penalty_params=penalty_params,
                                   allow_unstab=allow_unstab)
    residuals_std <- get_residuals(data=data, p=p, M=M, params=params,
                                   weight_function=weight_function, weightfun_pars=weightfun_pars,
                                   cond_dist=cond_dist, parametrization=parametrization,
                                   identification=identification, AR_constraints=AR_constraints,
                                   mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                   B_constraints=B_constraints, standardize=TRUE,
                                   penalized=penalized, penalty_params=penalty_params, allow_unstab=allow_unstab)
    IC <- get_IC(loglik=lok_and_tw$loglik, npars=npars, T_obs=nrow(data) - p)
    if(identification == "reduced_form" && cond_dist != "ind_Student" && cond_dist != "ind_skewed_t") {
      # Struct shocks are available for ind_Student and skewed_t mods
      structural_shocks <- NA
    } else {
      structural_shocks <- get_residuals(data=data, p=p, M=M, params=params,
                                         weight_function=weight_function, weightfun_pars=weightfun_pars,
                                         cond_dist=cond_dist, parametrization=parametrization,
                                         identification=identification, AR_constraints=AR_constraints,
                                         mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                         B_constraints=B_constraints, structural_shocks=TRUE,
                                         penalized=penalized, penalty_params=penalty_params, allow_unstab=allow_unstab)
    }
  }

  # Standard errors
  if(calc_std_errors) {
    if(is.null(data)) {
      warning("Approximate standard errors can't be calculated without data")
      std_errors <- rep(NA, npars)
    } else {
      std_errors <- tryCatch(standard_errors(data=data, p=p, M=M, params=params, weight_function=weight_function,
                                             weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                                             parametrization=parametrization, identification=identification,
                                             AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                             weight_constraints=weight_constraints, B_constraints=B_constraints,
                                             penalized=penalized, penalty_params=penalty_params, allow_unstab=allow_unstab),
                             error=function(e) {
                               warning("Approximate standard errors can't be calculated:")
                               warning(e)
                               std_errors=rep(NA, npars)
                             })
    }
  } else {
    std_errors <- rep(NA, npars)
  }

  # Conditional moments
  if(is.null(data)) {
    regime_cmeans <- total_cmeans <- total_ccovs <- NA
  } else {
    get_cm <- function(to_return) loglikelihood(data=data, p=p, M=M, params=params,
                                                weight_function=weight_function, weightfun_pars=weightfun_pars,
                                                cond_dist=cond_dist, parametrization=parametrization,
                                                identification=identification, AR_constraints=AR_constraints,
                                                mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                                B_constraints=B_constraints, to_return=to_return,
                                                check_params=TRUE, penalized=penalized, penalty_params=penalty_params,
                                                allow_unstab=allow_unstab, minval=NA)
    regime_cmeans <- get_cm("regime_cmeans")
    total_cmeans <- get_cm("total_cmeans")
    total_ccovs <- get_cm("total_ccovs")
  }

  # Some unconditional moments
  if(allow_unstab) {
    # Check whether the stability condition is satisfied
    params_std <- reform_constrained_pars(p=p, M=M, d=d, params=params, weight_function=weight_function,
                                          cond_dist=cond_dist, identification=identification,
                                          AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                          weight_constraints=weight_constraints, B_constraints=B_constraints,
                                          weightfun_pars=weightfun_pars)
    is_stable <- stab_conds_satisfied(p=p, M=M, d=d, params=params_std)
  } else { # Always stable if unstable is not allowed
    is_stable <- TRUE
  }
  if(is_stable) {
    regime_autocovs <- get_regime_autocovs(p=p, M=M, d=d, params=params, weight_function=weight_function,
                                           weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                                           identification=identification, AR_constraints=AR_constraints,
                                           mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                           B_constraints=B_constraints)
    regime_vars <- vapply(1:M, function(m) diag(regime_autocovs[, , 1, m]), numeric(d))
    regime_means <- get_regime_means(p=p, M=M, d=d, params=params, weight_function=weight_function,
                                     weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                                     parametrization=parametrization, identification=identification,
                                     AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                     weight_constraints=weight_constraints, B_constraints=B_constraints)
  } else { # Cannot calculate unconditional moments if the
    regime_autocovs <- array(NA, dim=c(d, d, p + 1, M))
    regime_vars <- matrix(NA, nrow=d, ncol=M)
    regime_means <- matrix(NA, nrow=d, ncol=M)
  }


  # Return
  structure(list(data=data,
                 model=list(p=p,
                            M=M,
                            d=d,
                            weight_function=weight_function,
                            weightfun_pars=weightfun_pars,
                            cond_dist=cond_dist,
                            parametrization=parametrization,
                            identification=identification,
                            AR_constraints=AR_constraints,
                            mean_constraints=mean_constraints,
                            weight_constraints=weight_constraints,
                            B_constraints=B_constraints),
                 params=params,
                 std_errors=std_errors,
                 transition_weights=lok_and_tw$tw,
                 regime_cmeans=regime_cmeans,
                 total_cmeans=total_cmeans,
                 total_ccovs=total_ccovs,
                 uncond_moments=list(regime_autocovs=regime_autocovs,
                                     regime_vars=regime_vars,
                                     regime_means=regime_means),
                 residuals_raw=residuals_raw,
                 residuals_std=residuals_std,
                 structural_shocks=structural_shocks,
                 loglik=structure(lok_and_tw$loglik,
                                  class="logLik",
                                  df=npars),
                 IC=IC,
                 penalized=penalized,
                 penalty_params=penalty_params,
                 allow_unstab=allow_unstab,
                 all_estimates=NULL,
                 all_logliks=NULL,
                 which_converged=NULL,
                 which_round=NULL),
            class="stvar")
}



#' @title Construct a STVAR model based on results from an arbitrary estimation round of \code{fitSTVAR}
#'
#' @description \code{alt_stvar} constructs a STVAR model based on results from an arbitrary estimation
#'   round of \code{fitSTVAR}
#'
#' @inheritParams get_boldA_eigens
#' @inheritParams STVAR
#' @param which_largest based on estimation round with which largest log-likelihood should the model be constructed?
#'   An integer value in 1,...,\code{nrounds}. For example, \code{which_largest=2} would take the second largest log-likelihood
#'   and construct the model based on the corresponding estimates.
#' @param which_round based on which estimation round should the model be constructed? An integer value in 1,...,\code{nrounds}.
#'   If specified, then \code{which_largest} is ignored.
#' @details It's sometimes useful to examine other estimates than the one with the highest log-likelihood. This function
#'   is wrapper around \code{STVAR} that picks the correct estimates from an object returned by \code{fitSTVAR}.
#' @seealso \code{\link{STVAR}}
#' @inherit STVAR references return
#' @examples
#' \donttest{
#' ## These are long-running examples that take approximately 10 seconds to run.
#'
#' # Estimate a Gaussian STVAR p=1, M=2 model with threshold weight function and
#' # the first lag of the second variable as the switching variables. Run only two
#' # estimation rounds and use the two-phase estimation method:
#' fit12 <- fitSTVAR(gdpdef, p=1, M=2, weight_function="threshold", weightfun_pars=c(2, 1),
#'  nrounds=2, seeds=c(1, 4), estim_method="two-phase")
#' fit12$loglik # Log-likelihood of the estimated model
#'
#' # Print the log-likelihood obtained from each estimation round:
#' fit12$all_logliks
#'
#' # Construct the model based on the second largest log-likelihood found in the
#' # estimation procedure:
#' fit12_alt <- alt_stvar(fit12, which_largest=2, calc_std_errors=FALSE)
#' fit12_alt$loglik # Log-likelihood of the alternative solution
#'
#' # Construct a model based on a specific estimation round, the first round:
#' fit12_alt2 <- alt_stvar(fit12, which_round=1, calc_std_errors=FALSE)
#' fit12_alt2$loglik # Log-likelihood of the alternative solution
#' }
#' @export

alt_stvar <- function(stvar, which_largest=1, which_round, calc_std_errors=FALSE) {
  check_stvar(stvar)
  stopifnot(!is.null(stvar$all_estimates))
  if(missing(which_round)) {
    stopifnot(which_largest >= 1 && which_largest <= length(stvar$all_estimates))
    which_round <- order(stvar$all_logliks, decreasing=TRUE)[which_largest]
  } else {
    stopifnot(which_round >= 1 && which_round <= length(stvar$all_estimates))
  }
  ret <- STVAR(data=stvar$data, p=stvar$model$p, M=stvar$model$M, d=stvar$model$d,
               params=stvar$all_estimates[[which_round]],
               weight_function=stvar$model$weight_function,
               weightfun_pars=stvar$model$weightfun_pars,
               cond_dist=stvar$model$cond_dist,
               parametrization=stvar$model$parametrization,
               identification=stvar$model$identification,
               AR_constraints=stvar$model$AR_constraints,
               mean_constraints=stvar$model$mean_constraints,
               weight_constraints=stvar$model$weight_constraints,
               B_constraints=stvar$model$B_constraints,
               penalized=stvar$penalized,
               penalty_params=stvar$penalty_params,
               allow_unstab=stvar$allow_unstab,
               calc_std_errors=calc_std_errors)

  # Pass the estimation results to the new object
  ret$all_estimates <- stvar$all_estimates
  ret$all_logliks <- stvar$all_logliks
  ret$which_converged <- stvar$which_converged
  if(!is.null(stvar$which_round)) {
    ret$which_round <- which_round
  }
  warn_eigens(ret, allow_unstab=stvar$allow_unstab)
  ret
}



#' @title Swap the parametrization of a STVAR model
#'
#' @description \code{swap_parametrization} swaps the parametrization of a STVAR model
#'  to \code{"mean"} if the current parametrization is \code{"intercept"}, and vice versa.
#'
#' @inheritParams diagnostic_plot
#' @inheritParams STVAR
#' @details \code{swap_parametrization} is a convenient tool if you have estimated the model in
#'  "intercept" parametrization but wish to work with "mean" parametrization in the future, or vice versa.
#' @inherit STVAR references return
#' @examples
#' ## Create a Gaussian STVAR p=1, M=2 model with the weighted relative stationary densities
#' # of the regimes as the transition weight function; use the intercept parametrization:
#' theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863,
#'   -0.3124, 0.505251, 0.298483, 0.030096, -0.176925, 0.838898, 0.310863, 0.007512,
#'   0.018244, 0.949533, -0.016941, 0.121403, 0.573269)
#' mod122 <- STVAR(p=1, M=2, d=2, params=theta_122relg, parametrization="intercept")
#' mod122$params[1:4] # The intercept parameters
#'
#' # Swap from the intercept parametrization to mean parametrization:
#' mod122mu <- swap_parametrization(mod122)
#' mod122mu$params[1:4] # The mean parameters
#'
#' # Swap back to the intercept parametrization:
#' mod122int <- swap_parametrization(mod122mu)
#' mod122int$params[1:4] # The intercept parameters
#'
#' ## Create a linear VAR(p=1) model with the intercept parametrization, include
#' # the two-variate data gdpdef to the model and calculate approximate standard errors:
#' theta_112 <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103,
#'   0.601786, -0.002945, 0.067224)
#' mod112 <- STVAR(data=gdpdef, p=1, M=1, params=theta_112, parametrization="intercept",
#'   calc_std_errors=TRUE)
#' print(mod112, standard_error_print=TRUE) # Standard errors are printed for the intercepts
#'
#' # To obtain standard errors for the unconditional means instead of the intercepts,
#' # swap to mean parametrization:
#' mod112mu <- swap_parametrization(mod112, calc_std_errors=TRUE)
#' print(mod112mu, standard_error_print=TRUE) # Standard errors are printed for the means
#' @export

swap_parametrization <- function(stvar, calc_std_errors=FALSE) {
  check_stvar(stvar)
  if(!is.null(stvar$model$mean_constraints)) {
    stop("Cannot change parametrization to intercept if the mean parameters are constrained")
  }
  change_to <- ifelse(stvar$model$parametrization == "intercept", "mean", "intercept")
  if(stvar$allow_unstab) {
    # Check whether the stability condition is satisfied
    params_std <- reform_constrained_pars(p=stvar$model$p, M=stvar$model$M, d=stvar$model$d, params=stvar$params,
                                          weight_function=stvar$model$weight_function, weightfun_pars=stvar$model$weightfun_pars,
                                          cond_dist=stvar$model$cond_dist, identification=stvar$model$identification,
                                          AR_constraints=stvar$model$AR_constraints, mean_constraints=stvar$model$mean_constraints,
                                          weight_constraints=stvar$model$weight_constraints, B_constraints=stvar$model$B_constraints)
    is_stable <- stab_conds_satisfied(p=stvar$model$p, M=stvar$model$M, d=stvar$model$d, params=params_std)
  } else { # Always stable if unstable is not allowed
    is_stable <- TRUE
  }
  if(!is_stable) {
    stop("Cannot swap parametrization if the model does not satisfy the usual stability condition (in all regimes)")
  }
  new_params <- change_parametrization(p=stvar$model$p, M=stvar$model$M, d=stvar$model$d, params=stvar$params,
                                       weight_function=stvar$model$weight_function, weightfun_pars=stvar$model$weightfun_pars,
                                       cond_dist=stvar$model$cond_dist, identification=stvar$model$identification,
                                       AR_constraints=stvar$model$AR_constraints, mean_constraints=stvar$model$mean_constraints,
                                       weight_constraints=stvar$model$weight_constraints, B_constraints=stvar$model$B_constraints,
                                       change_to=change_to)
  STVAR(data=stvar$data, p=stvar$model$p, M=stvar$model$M, d=stvar$model$d, params=new_params, parametrization=change_to,
        weight_function=stvar$model$weight_function, weightfun_pars=stvar$model$weightfun_pars,
        cond_dist=stvar$model$cond_dist, identification=stvar$model$identification,
        AR_constraints=stvar$model$AR_constraints, mean_constraints=stvar$model$mean_constraints,
        weight_constraints=stvar$model$weight_constraints, B_constraints=stvar$model$B_constraints,
        penalized=stvar$penalized, penalty_params=stvar$penalty_params, allow_unstab=stvar$allow_unstab,
        calc_std_errors=calc_std_errors)
}



#' @title Switch from two-regime reduced form STVAR model to a structural model identified by heteroskedasticity
#'
#' @description \code{get_hetsked_sstvar} constructs structural STVAR model identified by heteroskedasticity
#'   based on a reduced form STVAR model.
#'
#' @inheritParams fitSSTVAR
#' @details The switch is made by simultaneously diagonalizing the two error term covariance matrices
#'   with a well known matrix decomposition (Muirhead, 1982, Theorem A9.9) and then normalizing the
#'   diagonal of the matrix W positive (which implies positive diagonal of the impact matrix). Models with
#'   more that two regimes are not supported because the matrix decomposition does not generally
#'   exists for more than two covariance matrices.
#' @return Returns an object of class \code{'stvar'} defining a structural STVAR model identified by heteroskedasticity,
#'   with the main diagonal of the impact matrix normalized to be positive.
#' @seealso \code{\link{fitSSTVAR}}, \code{\link{STVAR}}, \code{\link{fitSTVAR}}
#'  \itemize{
#'    \item Muirhead R.J. 1982. Aspects of Multivariate Statistical Theory, \emph{Wiley}.
#'  }

get_hetsked_sstvar <- function(stvar, calc_std_errors=FALSE) {
  check_stvar(stvar)
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  data <- stvar$data
  pars_orig <- stvar$params
  weight_function <- stvar$model$weight_function
  cond_dist <- stvar$model$cond_dist
  parametrization <- stvar$model$parametrization
  identification <- stvar$model$identification
  if(identification != "reduced_form") {
    stop("Only reduced form models are supported")
  } else if(M != 2) {
    stop("Only two-regime models are supported")
  } else if(cond_dist == "ind_Student") {
    stop("Models with independent Student's t errors are not supported (they are readily statistically identified)")
  } else if(cond_dist == "ind_skewed_t") {
    stop("Models with independent skewed t errors are not supported (they are readily statistically identified)")
  }
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  weight_constraints <- stvar$model$weight_constraints
  B_constraints <- stvar$model$B_constraints
  weightfun_pars <- check_weightfun_pars(data=data, p=p, M=M, d=d, weight_function=weight_function,
                                         weightfun_pars=stvar$model$weightfun_pars, cond_dist=cond_dist)
  params <- reform_constrained_pars(p=p, M=M, d=d, params=pars_orig,
                                    weight_function=weight_function, weightfun_pars=weightfun_pars,
                                    cond_dist=cond_dist, identification=identification,
                                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                    weight_constraints=weight_constraints, B_constraints=B_constraints)

  # Obtain and simultaneously diagonalize the covariance matrices
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, cond_dist=cond_dist, identification=identification)
  tmp <- diag_Omegas(Omega1=all_Omega[, , 1], Omega2=all_Omega[, , 2])
  W <- matrix(tmp[1:(d^2)], nrow=d, ncol=d, byrow=FALSE)
  lambdas <- tmp[(d^2 + 1):(d^2 + d)]

  # Normalize the main diagonal of W to be positive
  for(i1 in 1:d) {
    if(W[i1, i1] < 0) W[,i1] <- -W[,i1]
  }

  ## Create the structural model parameter vector

  # First determine the numbers of each type of parameters in the original parameter vector,
  # so the new parameter vector can be constructed.
  if(is.null(mean_constraints)) {
    n_mean_pars <- M*d
  } else { # Means constrained
    n_mean_pars <- d*length(mean_constraints)
  }
  if(is.null(AR_constraints)) {
    n_ar_pars <- M*p*d^2
  } else { # AR matrices constrained
    n_ar_pars <- ncol(AR_constraints)
  }
  n_covmat_pars <- M*d*(d + 1)/2 # No B_constraints available here, as reduced form identification assumed
  if(is.null(weight_constraints)) {
    if(weight_function == "relative_dens" || weight_function == "threshold") {
      n_weight_pars <- M - 1
    } else if(weight_function == "logistic" || weight_function == "exponential") {
      n_weight_pars <- 2
    } else if(weight_function == "mlogit") {
      n_weight_pars <- (M - 1)*(1 + length(weightfun_pars[[1]])*weightfun_pars[[2]])
    }
  } else { # Constraints on the weight parameters
    if(all(weight_constraints[[1]] == 0)) {
      n_weight_pars <- 0 # alpha = r, not in the parameter vector
    } else {
      n_weight_pars <- ncol(weight_constraints[[1]]) # The dimension of xi
    }
  }

  mean_pars <- pars_orig[1:n_mean_pars]
  ar_pars <- pars_orig[(n_mean_pars + 1):(n_mean_pars + n_ar_pars)]
  weight_pars <- pars_orig[(n_mean_pars + n_ar_pars + n_covmat_pars
                            + 1):(n_mean_pars + n_ar_pars + n_covmat_pars + n_weight_pars)]
  if(cond_dist == "Gaussian") {
    dist_pars <- numeric(0)
  } else { # cond_dist == "Student"
    dist_pars <- pars_orig[length(pars_orig)] # df is the last param
  }
  new_params <- c(mean_pars, ar_pars, vec(W), lambdas, weight_pars, dist_pars)

  # Return the structural model identified by heteroskedasticity
  STVAR(data=data, p=p, M=M, d=d, params=new_params, cond_dist=cond_dist,
        weight_function=weight_function, weightfun_pars=weightfun_pars,
        parametrization=parametrization, identification="heteroskedasticity",
        AR_constraints=AR_constraints, mean_constraints=mean_constraints,
        weight_constraints=weight_constraints, B_constraints=NULL,
        penalized=stvar$penalized, penalty_params=stvar$penalty_params,
        allow_unstab=stvar$allow_unstab, calc_std_errors=calc_std_errors)
}


#' @title Reorder columns of impact matrix B (and lambda parameters if any) of a structural STVAR model
#'   that is identified by heteroskedasticity or non-Gaussianity.
#'
#' @description \code{reorder_B_columns} reorder columns of impact matrix B (and lambda parameters if any) of
#'   a structural STVAR model that is identified by heteroskedasticity or non-Gaussianity.
#'
#' @inheritParams STVAR
#' @param stvar a class 'stvar' object defining a structural STVAR model that is identified by heteroskedasticity
#'   or non-Gaussianity, typically created with \code{fitSSTVAR}.
#' @param perm an integer vector of length \eqn{d} specifying the new order of the columns of the impact matrix.
#'   For model identified by...
#'   \describe{
#'     \item{heteroskedasticity}{also lambda parameters of each regime will be reordered accordingly.}
#'     \item{non-Gaussianity}{the columns of the impact matrices of all the regimes and the component specific distribution
#'       parameters (degrees of freedom parameters) are reordered accordingly.}
#'   }
#' @details The order of the columns of the impact matrix can be changed without changing the implied reduced
#'   form model (as long as, for models identified by heteroskedasticity, the order of lambda parameters is also changed accordingly;
#'   and for model identified by non-Gaussianity, ordering of the columns of all the impact matrices and the component specific
#'   distribution  parameters is also changed accordingly). Note that constraints imposed on the impact matrix via \code{B_constraints}
#'   will also be modified accordingly.
#'
#'   Also all signs in any column of impact matrix can be swapped (without changing the implied reduced form model)
#'   with the function \code{swap_B_signs}. This obviously also swaps the sign constraints (if any) in the corresponding columns of
#'   the impact matrix.
#' @inherit STVAR return
#' @seealso \code{\link{GIRF}}, \code{\link{fitSSTVAR}}, \code{\link{swap_B_signs}}
#' @references
#'  \itemize{
#'    \item Lütkepohl H., Netšunajev A. 2018. Structural vector autoregressions with smooth transition in variances.
#'      \emph{Journal of Economic Dynamics & Control}, \strong{84}, 43-57.
#'  }
#' @examples
#' # Create a structural two-variate Student's t STVAR p=2, M=2 model with logistic transition
#' # weights and the first lag of the second variable  as the switching variable, and shocks
#' # identified by heteroskedasticity:
#' theta_222logt <- c(0.356914, 0.107436, 0.356386, 0.086330, 0.139960, 0.035172, -0.164575,
#'   0.386816, 0.451675, 0.013086, 0.227882, 0.336084, 0.239257, 0.024173, -0.021209, 0.707502,
#'   0.063322, 0.027287, 0.009182, 0.197066, -0.03, 0.24, -0.76, -0.02, 3.36, 0.86, 0.1, 0.2, 7)
#' mod222logt <- STVAR(p=2, M=2, d=2, params=theta_222logt, weight_function="logistic",
#'   weightfun_pars=c(2, 1), cond_dist="Student", identification="heteroskedasticity")
#'
#' # Print the parameter values, W and lambdas are printed in the bottom:
#' mod222logt
#'
#' # Reverse the ordering of the columns of W (or equally the impact matrix):
#' mod222logt_rev <- reorder_B_columns(mod222logt, perm=c(2, 1))
#' mod222logt_rev # The columns of the impact matrix are in a reversed order
#'
#' # Swap the ordering of the columns of the impact matrix back to the original:
#' mod222logt_rev2 <- reorder_B_columns(mod222logt_rev, perm=c(2, 1))
#' mod222logt_rev2 # The columns of the impact matrix are back in the original ordering
#'
#' # Below code does not do anything, as perm=1:2, so the ordering does not change:
#' mod222logt3 <- reorder_B_columns(mod222logt, perm=c(1, 2))
#' mod222logt3 # The ordering of the columns did not change from the original
#' @export

reorder_B_columns <- function(stvar, perm, calc_std_errors=FALSE) {
  check_stvar(stvar)
  cond_dist <- stvar$model$cond_dist
  identification <- stvar$model$identification
  if(cond_dist %in% c("ind_Student", "ind_skewed_t")) identification <- "non-Gaussianity" # ind_Stud and skewed_t mods are readily identified
  if(identification != "heteroskedasticity" && identification != "non-Gaussianity") {
    stop("Only model identified by heteroskedasticity or non-Gaussianity are supported!")
  }
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  weight_function <- stvar$model$weight_function
  weightfun_pars <- check_weightfun_pars(data=stvar$data, p=p, M=M, d=d, weight_function=weight_function,
                                         weightfun_pars=stvar$model$weightfun_pars)

  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  weight_constraints <- stvar$model$weight_constraints
  B_constraints <- stvar$model$B_constraints
  params <- reform_constrained_pars(p=p, M=M, d=d, params=stvar$params,
                                    weight_function=weight_function, weightfun_pars=weightfun_pars,
                                    cond_dist=cond_dist, identification=identification,
                                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                    weight_constraints=weight_constraints, B_constraints=B_constraints)

  ## The impact matrix and distribution parameters
  if(identification == "heteroskedasticity") {
    W <- pick_W(p=p, M=M, d=d, params=params, identification=identification)
    lambdas <- pick_lambdas(p=p, M=M, d=d, params=params, identification=identification)
  } else { # identification == "non-Gaussianity"
    all_B <- pick_Omegas(p=p, M=M, d=d, params=params, cond_dist=cond_dist, identification=identification)
    distpars <- pick_distpars(d=d, params=params, cond_dist=cond_dist)
  }
  n_distpars <- length(pick_distpars(d=d, params=params, cond_dist=cond_dist))

  ## Calculate the number of weight parameters
  if(weight_function == "exogenous") {
    n_weight_pars <- 0
  } else {
    if(is.null(weight_constraints)) {
      if(weight_function == "relative_dens" || weight_function == "threshold") {
        n_weight_pars <- M - 1
      } else if(weight_function == "logistic" || weight_function == "exponential") {
        n_weight_pars <- 2
      } else if(weight_function == "mlogit") {
        n_weight_pars <- (M - 1)*(1 + length(weightfun_pars[[1]])*weightfun_pars[[2]])
      }
    } else { # Constraints on the weight parameters
      if(all(weight_constraints[[1]] == 0)) {
        n_weight_pars <- 0 # alpha = r, not in the parameter vector
      } else {
        n_weight_pars <- ncol(weight_constraints[[1]]) # The dimension of xi
      }
    }
  }

  ## Create the new parameter vector
  if(identification == "heteroskedasticity") {
    W <- W[, perm]
    if(is.null(B_constraints)) {
      W <- vec(W) # If any zeros, they are specified by hand and included in the param vector
    } else {
      W <- Wvec(W) # Zeros removed
    }
    if(M > 1) {
      lambdas <- matrix(lambdas, nrow=d, ncol=M - 1, byrow=FALSE)
      lambdas <- vec(lambdas[perm,])
    }
    new_params <- stvar$params
    new_params[(length(new_params) - (n_weight_pars + length(W) + length(lambdas) + n_distpars)
                + 1):(length(new_params) - (n_weight_pars + n_distpars))] <- c(W, lambdas)
  } else { # identification == "non-Gaussianity"
    new_all_B <- numeric(0)
    for(m in 1:M) { # Reorder the columns of the impact matrices
      if(is.null(B_constraints)) {
        new_all_B <- c(new_all_B, vec(all_B[, perm, m])) # If any zeros, they are specified by hand and included in the param vector
      } else {
        new_all_B <- c(new_all_B, Wvec(all_B[, perm, m])) # Zeros removed (assumes no exact zeros in the estimates)
      }
    }
    new_params <- stvar$params
    new_params[(length(new_params) - (n_weight_pars + length(new_all_B) + n_distpars)
                + 1):(length(new_params) - (n_weight_pars + n_distpars))] <- new_all_B # New impact matrix params

    if(cond_dist == "ind_Student") {
      new_distpars <- distpars[perm] # Reorder the degrees of freedom parameters
    } else if(cond_dist == "ind_skewed_t") { # cond_dist == "ind_skewed_t"
      new_distpars <- c(distpars[1:d][perm], distpars[(d + 1):length(distpars)][perm]) # reorder df and skewness params
    } # Other options do not end up here as, we are inside identificatiom by non-Gaussianity

    new_params[(length(new_params) - n_distpars + 1):length(new_params)] <- new_distpars # Insert new distribution parameters
  }

  ## Reorder the columns of the B_constraints accordingly
  if(!is.null(B_constraints)) {
    new_B_constraints <- B_constraints[, perm]
  } else {
    new_B_constraints <- NULL
  }

  ## Construct the STVAR model based on the obtained structural parameters
  STVAR(data=stvar$data, p=p, M=M, d=d, params=new_params, weight_function=weight_function,
        weightfun_pars=weightfun_pars, cond_dist=cond_dist, parametrization=stvar$model$parametrization,
        identification=identification, AR_constraints=AR_constraints, mean_constraints=mean_constraints,
        weight_constraints=weight_constraints, B_constraints=new_B_constraints,
        penalized=stvar$penalized, penalty_params=stvar$penalty_params, allow_unstab=stvar$allow_unstab,
        calc_std_errors=calc_std_errors)
}



#' @title Swap all signs in pointed columns of the impact matrix of a structural STVAR model
#'   that is identified by heteroskedasticity or non-Gaussianity
#'
#' @description \code{swap_B_signs} swaps all signs in pointed columns of the impact matrix of
#'   a structural STVAR model that is identified by heteroskedasticity or non-Gaussianity.
#'
#' @inheritParams reorder_B_columns
#' @param which_to_swap a numeric vector of length at most \eqn{d} and elemnts in \eqn{1,..,d}
#'   specifying the columns of the impact matrix whose sign should be swapped.
#' @details All signs in any column of the impact matrix can be swapped without changing the implied reduced form model.
#'   For model identified by non-Gaussianity, the signs of the columns of the impact matrices of all the regimes are
#'   swapped accordingly. Note that the sign constraints imposed on the impact matrix via \code{B_constraints} are also
#'   swapped in the corresponding columns accordingly.
#'
#'   Also the order of the columns of the impact matrix can be changed (without changing the implied reduced
#'   form model) as long as the ordering of other related parameters is also changed accordingly. This can be
#'   done with the function \code{reorder_B_columns}.
#' @inherit STVAR return
#' @seealso \code{\link{GIRF}}, \code{\link{fitSSTVAR}}, \code{\link{reorder_B_columns}}
#' @inherit reorder_B_columns references
#' @examples
#' # Create a structural two-variate Student's t STVAR p=2, M=2, model with logistic transition
#' # weights and the first lag of the second variable as the switching variable, and shocks
#' # identified by heteroskedasticity:
#' theta_222logt <- c(0.356914, 0.107436, 0.356386, 0.086330, 0.139960, 0.035172, -0.164575,
#'   0.386816, 0.451675, 0.013086, 0.227882, 0.336084, 0.239257, 0.024173, -0.021209, 0.707502,
#'   0.063322, 0.027287, 0.009182, 0.197066, -0.03, 0.24, -0.76, -0.02, 3.36, 0.86, 0.1, 0.2, 7)
#' mod222logt <- STVAR(p=2, M=2, d=2, params=theta_222logt, weight_function="logistic",
#'   weightfun_pars=c(2, 1), cond_dist="Student", identification="heteroskedasticity")
#'
#' # Print the parameter values, W and lambdas are printed in the bottom:
#' mod222logt
#'
#' # Swap the signs of the first column of W (or equally the impact matrix):
#' mod222logt2 <- swap_B_signs(mod222logt, which_to_swap=1)
#' mod222logt2 # The signs of the first column of the impact matrix are swapped
#'
#' # Swap the signs of the second column of the impact matrix:
#' mod222logt3 <- swap_B_signs(mod222logt, which_to_swap=2)
#' mod222logt3 # The signs of the second column of the impact matrix are swapped
#'
#' # Swap the signs of both columns of the impact matrix:
#' mod222logt4 <- swap_B_signs(mod222logt, which_to_swap=1:2)
#' mod222logt4 # The signs of both columns of the impact matrix are swapped
#' @export

swap_B_signs <- function(stvar, which_to_swap, calc_std_errors=FALSE) {
  check_stvar(stvar)
  cond_dist <- stvar$model$cond_dist
  identification <- stvar$model$identification
  if(cond_dist %in% c("ind_Student", "ind_skewed_t")) identification <- "non-Gaussianity" # ind_Stud and skewed_t mods are readily identified
  if(identification != "heteroskedasticity" && identification != "non-Gaussianity") {
    stop("Only model identified by heteroskedasticity or non-Gaussianity are supported!")
  }
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  weight_function <- stvar$model$weight_function
  weightfun_pars <- check_weightfun_pars(stvar$data, p=p, M=M, d=d, weight_function=weight_function,
                                         weightfun_pars=stvar$model$weightfun_pars)
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  weight_constraints <- stvar$model$weight_constraints
  B_constraints <- stvar$model$B_constraints
  params <- reform_constrained_pars(p=p, M=M, d=d, params=stvar$params,
                                    weight_function=weight_function, weightfun_pars=weightfun_pars,
                                    cond_dist=cond_dist, identification=identification,
                                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                    weight_constraints=weight_constraints, B_constraints=B_constraints)

  ## The impact matrix and distribution parameters
  if(identification == "heteroskedasticity") {
    W <- pick_W(p=p, M=M, d=d, params=params, identification=identification)
    lambdas <- pick_lambdas(p=p, M=M, d=d, params=params, identification=identification)
  } else { # identification == "non-Gaussianity"
    all_B <- pick_Omegas(p=p, M=M, d=d, params=params, cond_dist=cond_dist, identification=identification)
    distpars <- pick_distpars(d=d, params=params, cond_dist=cond_dist)
  }
  n_distpars <- length(pick_distpars(d=d, params=params, cond_dist=cond_dist))

  ## Calculate the number of weight parameters
  if(weight_function == "exogenous") {
    n_weight_pars <- 0
  } else {
    if(is.null(weight_constraints)) {
      if(weight_function == "relative_dens" || weight_function == "threshold") {
        n_weight_pars <- M - 1
      } else if(weight_function == "logistic" || weight_function == "exponential") {
        n_weight_pars <- 2
      } else if(weight_function == "mlogit") {
        n_weight_pars <- (M - 1)*(1 + length(weightfun_pars[[1]])*weightfun_pars[[2]])
      }
    } else { # Constraints on the weight parameters
      if(all(weight_constraints[[1]] == 0)) {
        n_weight_pars <- 0 # alpha = r, not in the parameter vector
      } else {
        n_weight_pars <- ncol(weight_constraints[[1]]) # The dimension of xi
      }
    }
  }

  ## Create the new parameter vector
  if(identification == "heteroskedasticity") {
    W[, which_to_swap] <- -W[, which_to_swap]
    if(is.null(B_constraints)) {
      W <- vec(W) # If any zeros, they are specified by hand and included in the param vector
    } else {
      W <- Wvec(W) # Zeros removed
    }
    r <- d*(M - 1) # The number of lambda parameters
    new_params <- stvar$params
    new_params[(length(new_params) - (n_weight_pars + length(W) + r
                                      + n_distpars) + 1):(length(new_params) - (n_weight_pars + r + n_distpars))] <- W
  } else { # identification by non-Gaussianity
    new_all_B <- numeric(0)
    for(m in 1:M) { # Swap the signs of the columns of the impact matrices
      all_B[, which_to_swap, m] <- -all_B[, which_to_swap, m]
      if(is.null(B_constraints)) {
        new_all_B <- c(new_all_B, vec(all_B[, , m])) # If any zeros, they are specified by hand and included in the param vector
      } else {
        new_all_B <- c(new_all_B, Wvec(all_B[, , m])) # Zeros removed (assumes no exact zeros in the estimates)
      }
    }
    new_params <- stvar$params
    new_params[(length(new_params) - (n_weight_pars + length(new_all_B) + n_distpars)
                + 1):(length(new_params) - (n_weight_pars + n_distpars))] <- new_all_B # New impact matrix params
  }


  ## Swap the sign constraints accordingly
  if(!is.null(B_constraints)) {
    new_B_constraints <- B_constraints
    new_B_constraints[, which_to_swap] <- -new_B_constraints[, which_to_swap]
  } else {
    new_B_constraints <- NULL
  }

  # Construct the SSTVAR model based on the obtained structural parameters
  STVAR(data=stvar$data, p=p, M=M, d=d, params=new_params, weight_function=weight_function,
        weightfun_pars=weightfun_pars, cond_dist=cond_dist, parametrization=stvar$model$parametrization,
        identification=identification, AR_constraints=AR_constraints, mean_constraints=mean_constraints,
        weight_constraints=weight_constraints, B_constraints=new_B_constraints,
        penalized=stvar$penalized, penalty_params=stvar$penalty_params, allow_unstab=stvar$allow_unstab,
        calc_std_errors=calc_std_errors)
}


#' @title Filter inappropriate the estimates produced by fitSTVAR
#'
#' @description \code{filter_estimates} filters out inappropriate estimates produced by \code{fitSTVAR}:
#'   can be used to obtain the (possibly) appropriate estimate with the largest found log-likelihood
#'   (among possibly appropriate estimates) as well as (possibly) appropriate estimates based on smaller
#'   log-likelihoods.
#'
#' @inheritParams reorder_B_columns
#' @param which_largest an integer at least one specifying the (possibly) appropriate estimate corresponding
#'  to which largest log-likelihood should be returned. E.g., if \code{which_largest=2}, the function will
#'  return among the estimates that it does not deem inappropriate the one that has the second largest log-likelihood.
#' @param filter_stab Should estimates close to breaking the usual stability condition be filtered out?
#' @details The function goes through the estimates produced by \code{fitSTVAR} and checks which estimates are
#'  deemed inappropriate. That is, estimates that are not likely solutions of interest. Specifically, solutions
#'  that incorporate a near-singular error term covariance matrix (any eigenvalue less than \eqn{0.002}),
#'  any modulus of the eigenvalues of the companion form AR matrices larger than $0.9985$ (indicating the
#'  necessary condition for stationarity is close to break), or transition weights such that they are close to zero
#'  for almost all \eqn{t} for at least one regime. Then, among the solutions are not deemed inappropriate, it
#'  returns a STVAR models based on the estimate that has the \code{which_largest} largest log-likelihood.
#'
#'  The function \code{filter_estimates} is kind of a version of \code{alt_stvar} that only considers estimates
#'  that are not deemed inappropriate
#' @inherit STVAR return
#' @seealso \code{\link{fitSTVAR}}, \code{\link{alt_stvar}}
#' @examples
#' \donttest{
#'  # Fit a two-regime STVAR model with logistic transition weights and Student's t errors,
#'  # and use two-phase estimation method:
#'  fit12 <- fitSTVAR(gdpdef, p=1, M=2, weight_function="logistic", weightfun_pars=c(2, 1),
#'   cond_dist="Student", nrounds=2, ncores=2, seeds=1:2, estim_method="two-phase")
#'  fit12
#'
#'  # Filter through inappropriate estimates and obtain the second best appropriate solution:
#'  fit12_2 <- filter_estimates(fit12, which_largest=2)
#'  fit12_2 # The same model since the two estimation rounds yielded the same estimate
#' }
#' @export

filter_estimates <- function(stvar, which_largest=1, filter_stab=TRUE, calc_std_errors=FALSE) {
  check_stvar(stvar)
  stopifnot(all_pos_ints(which_largest))
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  data <- stvar$data
  pars_orig <- stvar$params
  weight_function <- stvar$model$weight_function
  cond_dist <- stvar$model$cond_dist
  parametrization <- stvar$model$parametrization
  identification <- stvar$model$identification
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  weight_constraints <- stvar$model$weight_constraints
  B_constraints <- stvar$model$B_constraints
  weightfun_pars <- stvar$model$weightfun_pars
  all_estimates <- stvar$all_estimates
  all_logliks <- stvar$all_logliks
  penalized <- stvar$penalized
  penalty_params <- stvar$penalty_params
  allow_unstab <- stvar$allow_unstab
  red_criteria <- c(0.05, 0.01)
  n_obs <- nrow(data)
  if(identification != "reduced_form") stop("Only reduced form models are supported!")
  if(is.null(all_estimates)) stop("No multiple estimates found in the model object; was it created by fitSTVAR?")
  if(is.null(all_logliks)) stop("No multiple log-likelihoods found in the model object; was it created by fitSTVAR?")
  if(length(all_estimates) != length(all_logliks)) stop("The number of estimates and log-likelihoods do not match")

  # Ordering from the largest loglik to the smallest:
  ord_by_loks <- order(all_logliks, decreasing=TRUE)

  # Go through estimates, take the estimate that yield the higher likelihood among estimates that are do not
  # include wasted regimes or near-singular error term covariance matrices.
  where_at_which_largest <- 0
  for(i1 in 1:length(all_estimates)) {
    which_round <- ord_by_loks[i1] # Est round with i1:th largest loglik
    pars <- all_estimates[[which_round]]
    pars_std <- reform_constrained_pars(p=p, M=M, d=d, params=pars,
                                        weight_function=weight_function, weightfun_pars=weightfun_pars,
                                        cond_dist=cond_dist, identification=identification,
                                        AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                        weight_constraints=weight_constraints,
                                        B_constraints=NULL) # Pars in standard form for pick pars fns
    # Check Omegas
    Omega_eigens <- get_omega_eigens_par(p=p, M=M, d=d, params=pars_std,
                                         weight_function=weight_function, weightfun_pars=weightfun_pars,
                                         cond_dist=cond_dist, identification=identification,
                                         AR_constraints=NULL, mean_constraints=NULL,
                                         weight_constraints=NULL, B_constraints=NULL)
    Omegas_ok <- !any(Omega_eigens < 0.002)

    # Checks AR matrices
    if(filter_stab) {
      boldA_eigens <- get_boldA_eigens_par(p=p, M=M, d=d, params=pars_std,
                                           weight_function=weight_function, weightfun_pars=weightfun_pars,
                                           cond_dist=cond_dist, identification=identification,
                                           AR_constraints=NULL, mean_constraints=NULL,
                                           weight_constraints=NULL, B_constraints=NULL)
      stat_ok <- !any(boldA_eigens > 0.9985)
    } else {
      stat_ok <- TRUE
    }

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
                              identification=identification, AR_constraints=NULL,
                              mean_constraints=NULL, B_constraints=NULL, weight_constraints=NULL,
                              penalized=penalized, penalty_params=penalty_params, allow_unstab=allow_unstab,
                              to_return="tw", check_params=TRUE, minval=matrix(0, nrow=n_obs-p, ncol=M))
    tweights_ok <- !any(vapply(1:M, function(m) sum(tweights[,m] > red_criteria[1]) < red_criteria[2]*n_obs, logical(1)))
    if(Omegas_ok && stat_ok && tweights_ok && weightpars_ok) {
      where_at_which_largest <- where_at_which_largest + 1
      if(where_at_which_largest == which_largest) {
        which_round <- which_round # The estimation round of the appropriate estimate with the which_largest largest loglik
        message(paste("Filtered through", i1-1, "'estimates' with a larger log-likelihood"))
        break
      }
    }
    if(i1 == length(all_estimates)) {
      message(paste("No (more) 'appropriate' estimates found! Returing the supplied model."))
      return(stvar)
    }
  }

  # Build a STVAR model from the obtained estimates
  ret <- STVAR(data=data, p=p, M=M, d=d,
               params=all_estimates[[which_round]],
               weight_function=weight_function,
               weightfun_pars=weightfun_pars,
               cond_dist=cond_dist,
               parametrization=parametrization,
               identification=identification,
               AR_constraints=AR_constraints,
               mean_constraints=mean_constraints,
               weight_constraints=weight_constraints,
               B_constraints=B_constraints,
               penalized=penalized,
               penalty_params=penalty_params,
               allow_unstab=allow_unstab,
               calc_std_errors=calc_std_errors)

  # Pass the estimation results to the new object
  ret$all_estimates <- all_estimates
  ret$all_logliks <- all_logliks
  ret$which_converged <- stvar$which_converged
  if(!is.null(stvar$which_round)) {
    ret$which_round <- which_round
  }
  warn_eigens(ret, allow_unstab=allow_unstab)
  ret
}



#' @title Update STVAR model estimated with a version of the package <1.1.0 to
#'   be compatible with the versions >=1.1.0.
#'
#' @description \code{update_stvar_to_sstvar110} updates a STVAR model estimated with a version of
#'   the package <1.1.0 to be compatible with the versions >=1.1.0 by adding the elements
#'   \code{$penalized}, \code{$penalty_params}, and \code{$allow_unstab} to the model object.
#'
#' @inheritParams diagnostic_plot
#' @details The function is useful when a STVAR model estimated with a version of the package <1.1.0.
#'   Does not do anything if the elements \code{$penalized}, \code{$penalty_params}, and \code{$allow_unstab}
#'   are already containing in the model object.
#' @return Returns an object of class \code{'stvar'} with the elements \code{$penalized},
#'  \code{$penalty_params}, and \code{$allow_unstab} added to it if they were missing.
#' @examples
#' # Linear Gaussian VAR(p=1) model:
#' theta_112 <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103,
#'   0.601786, -0.002945, 0.067224)
#' mod112 <- STVAR(data=gdpdef, p=1, M=1, params=theta_112)
#'
#' # Update to include the new elements (does not do anything they are already
#' # included):
#' mod112 <- stvar_to_sstvars110(mod112)
#' @export

stvar_to_sstvars110 <- function(stvar) {
  check_stvar(stvar)
  if(!is.null(stvar$penalized) && !is.null(stvar$penalty_params) && !is.null(stvar$allow_unstab)) {
    return(stvar)
  }
  if(!is.null(stvar$penalized)) {
    stvar$penalized <- FALSE
  }
  if(!is.null(stvar$penalty_params)) {
    stvar$penalty_params <- c(0.05, 0)
  }
  if(!is.null(stvar$allow_unstab)) {
    stvar$allow_unstab <- FALSE
  }
  stvar
}
