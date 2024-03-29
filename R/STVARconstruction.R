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
#' @return Returns an object of class \code{'stvar'} defining the specified reduced form or structural
#'  smooth transition VAR model. Can be used to work with other functions provided in \code{sstvars}.
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
#'    \item Kheifets I.L., Saikkonen P.J. 2020. Stationarity and ergodicity of Vector STAR models.
#'      \emph{Econometric Reviews}, \strong{39}:4, 407-414.
#'    \item Lütkepohl H., Netšunajev A. 2017. Structural vector autoregressions with smooth transition in variances.
#'      \emph{Journal of Economic Dynamics & Control}, \strong{84}, 43-57.
#'    \item Tsay R. 1998. Testing and Modeling Multivariate Threshold Models.
#'      \emph{Journal of the American Statistical Association}, \strong{93}:443, 1188-1202.
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
#' mod222logtsh <- STVAR(p=2, M=2, d=2, params=params22log,
#' weight_function="logistic", weightfun_pars=c(2, 2), cond_dist="Student",
#' identification="heteroskedasticity")
#' print(mod222logtsh) # Printout of the model
#'
#' # Linear Gaussian VAR(p=1) model:
#' theta_112 <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103,
#'   0.601786, -0.002945, 0.067224)
#' mod112 <- STVAR(data=gdpdef, p=1, M=1, params=theta_112)
#' summary(mod112) # Summary printout
#' @export

STVAR <- function(data, p, M, d, params, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold"),
                  weightfun_pars=NULL, cond_dist=c("Gaussian", "Student"), parametrization=c("intercept", "mean"),
                  identification=c("reduced_form", "recursive", "heteroskedasticity"),
                  AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL,
                  calc_std_errors=FALSE) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  identification <- match.arg(identification)
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
  check_pMd(p=p, M=M, d=d, weight_function=weight_function, identification=identification)
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                         cond_dist=cond_dist)
  check_constraints(p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    parametrization=parametrization, identification=identification,
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=B_constraints)
  check_params(p=p, M=M, d=d, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
               cond_dist=cond_dist, parametrization=parametrization, identification=identification,
               AR_constraints=AR_constraints, mean_constraints=mean_constraints,
               weight_constraints=weight_constraints, B_constraints=B_constraints)
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
                                check_params=TRUE, minval=NA)
    residuals_raw <- get_residuals(data=data, p=p, M=M, params=params,
                                   weight_function=weight_function, weightfun_pars=weightfun_pars,
                                   cond_dist=cond_dist, parametrization=parametrization,
                                   identification=identification, AR_constraints=AR_constraints,
                                   mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                   B_constraints=B_constraints, standardize=FALSE)
    residuals_std <- get_residuals(data=data, p=p, M=M, params=params,
                                   weight_function=weight_function, weightfun_pars=weightfun_pars,
                                   cond_dist=cond_dist, parametrization=parametrization,
                                   identification=identification, AR_constraints=AR_constraints,
                                   mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                   B_constraints=B_constraints, standardize=TRUE)
    IC <- get_IC(loglik=lok_and_tw$loglik, npars=npars, T_obs=nrow(data) - p)
    if(identification == "reduced_form") {
      structural_shocks <- NA
    } else {
      structural_shocks <- get_residuals(data=data, p=p, M=M, params=params,
                                         weight_function=weight_function, weightfun_pars=weightfun_pars,
                                         cond_dist=cond_dist, parametrization=parametrization,
                                         identification=identification, AR_constraints=AR_constraints,
                                         mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                         B_constraints=B_constraints, structural_shocks=TRUE)
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
                                             weight_constraints=weight_constraints, B_constraints=B_constraints),
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
                                                check_params=TRUE, minval=NA)
    regime_cmeans <- get_cm("regime_cmeans")
    total_cmeans <- get_cm("total_cmeans")
    total_ccovs <- get_cm("total_ccovs")
  }

  # Some unconditional moments
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
#' @inherit STVAR references return
#' @examples
#' \donttest{
#' ## These are long-running examples that take approximately 10 seconds to run.
#'
#' # Estimate a Gaussian STVAR p=1, M=2 model with exponential weight function and
#' # the first lag of the second variable as the switching variables. Run only two
#' # estimation rounds:
#' fit12 <- fitSTVAR(gdpdef, p=1, M=2, weight_function="exponential", weightfun_pars=c(2, 1),
#'  nrounds=2, seeds=c(1, 7))
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
#' # Construct a model based on a specific estimation round, the second round:
#' fit12_alt2 <- alt_stvar(fit12, which_round=2, calc_std_errors=FALSE)
#' fit12_alt2$loglik # Log-likelihood of the alternative solution
#' }
#' @export

alt_stvar <- function(stvar, which_largest=1, which_round, calc_std_errors=TRUE) {
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
               calc_std_errors=calc_std_errors)

  # Pass the estimation results to the new object
  ret$all_estimates <- stvar$all_estimates
  ret$all_logliks <- stvar$all_logliks
  ret$which_converged <- stvar$which_converged
  if(!is.null(stvar$which_round)) {
    ret$which_round <- which_round
  }
  warn_eigens(ret)
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
#' print_std_errors(mod112) # Standard errors are printed for the intercepts
#'
#' # To obtain standard errors for the unconditional means instead of the intercepts,
#' # swap to mean parametrization:
#' mod112mu <- swap_parametrization(mod112, calc_std_errors=TRUE)
#' print_std_errors(mod112mu) # Standard errors are printed for the means
#' @export

swap_parametrization <- function(stvar, calc_std_errors=FALSE) {
  check_stvar(stvar)
  if(!is.null(stvar$model$mean_constraints)) {
    stop("Cannot change parametrization to intercept if the mean parameters are constrained")
  }
  change_to <- ifelse(stvar$model$parametrization == "intercept", "mean", "intercept")
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
#' @seealso \code{\link{fitSSTVAR}}, \code{\link{STVAR}}, code{\link{fitSTVAR}}
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
  stopifnot(identification == "reduced_form")
  stopifnot(M == 2)
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  weight_constraints <- stvar$model$weight_constraints
  B_constraints <- stvar$model$B_constraints
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function,
                                         weightfun_pars=stvar$model$weightfun_pars, cond_dist=cond_dist)
  params <- reform_constrained_pars(p=p, M=M, d=d, params=pars_orig,
                                    weight_function=weight_function, weightfun_pars=weightfun_pars,
                                    cond_dist=cond_dist, identification=identification,
                                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                    weight_constraints=weight_constraints, B_constraints=B_constraints)

  # Obtain and simultaneously diagonalize the covariance matrices
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, identification=identification)
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
        calc_std_errors=calc_std_errors)
}


#' @title Reorder columns of the W-matrix and lambda parameters of a structural STVAR model
#'   that is identified by heteroskedasticity.
#'
#' @description \code{reorder_W_columns} reorder columns of the W-matrix and lambda parameters
#'   of a structural STVAR model that is identified by heteroskedasticity.
#'
#' @inheritParams STVAR
#' @param stvar a class 'stvar' object defining a structural STVAR model that is identified by heteroskedasticity,
#'   typically created with \code{fitSSTVAR}.
#' @param perm an integer vector of length \eqn{d} specifying the new order of the columns of \eqn{W}.
#'   Also lambda parameters of each regime will be reordered accordingly.
#' @details The order of the columns of \eqn{W} can be changed without changing the implied reduced
#'   form model as long as the order of lambda parameters is also changed accordingly. Note that the
#'   constraints imposed on \eqn{W} (or the  impact matrix) will also be modified accordingly.
#'
#'   This function does not support models with constraints imposed on the lambda parameters!
#'
#'   Also all signs in any column of \eqn{W} can be swapped (without changing the implied reduced form model)
#'   with the function \code{swap_W_signs} but this obviously also swaps the sign constraints in the
#'   corresponding columns of \eqn{W}.
#' @return Returns an object of class \code{'stvar'} defining a structural STVAR model with the columns of \eqn{W}
#'   reordered.
#' @seealso \code{\link{GIRF}}, \code{\link{fitSSTVAR}}, \code{\link{swap_W_signs}}
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
#' mod222logt_rev <- reorder_W_columns(mod222logt, perm=c(2, 1))
#' mod222logt_rev # The columns of W are in a reversed order
#'
#' # Swap the ordering of the columns of W back to the original:
#' mod222logt_rev2 <- reorder_W_columns(mod222logt_rev, perm=c(2, 1))
#' mod222logt_rev2 # The columns of W are back in the original ordering
#'
#' # Below code does not do anything, as perm=1:2, so the ordering does not change:
#' mod222logt3 <- reorder_W_columns(mod222logt, perm=c(1, 2))
#' mod222logt3 # The ordering of the columns did not change from the original
#' @export

reorder_W_columns <- function(stvar, perm, calc_std_errors=FALSE) {
  check_stvar(stvar)
  if(stvar$model$identification != "heteroskedasticity") stop("Only model identified by heteroskedasticity are supported!")
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  weight_function <- stvar$model$weight_function
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function, weightfun_pars=stvar$model$weightfun_pars)
  cond_dist <- stvar$model$cond_dist
  identification <- stvar$model$identification
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  weight_constraints <- stvar$model$weight_constraints
  B_constraints <- stvar$model$B_constraints
  params <- reform_constrained_pars(p=p, M=M, d=d, params=stvar$params,
                                    weight_function=weight_function, weightfun_pars=weightfun_pars,
                                    cond_dist=cond_dist, identification=identification,
                                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                    weight_constraints=weight_constraints, B_constraints=B_constraints)
  W <- pick_W(p=p, M=M, d=d, params=params, identification=identification)
  lambdas <- pick_lambdas(p=p, M=M, d=d, params=params, identification=identification)
  n_distpars <- length(pick_distpars(params=params, cond_dist=cond_dist))

  # Calculate the number of weight parameters
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

  # Create the new parameter vector
  W <- W[, perm]
  W <- Wvec(W) # Zeros removed
  if(M > 1) {
    lambdas <- matrix(lambdas, nrow=d, ncol=M - 1, byrow=FALSE)
    lambdas <- vec(lambdas[perm,])
  }
  new_params <- stvar$params
  new_params[(length(new_params) - (n_weight_pars + length(W) + length(lambdas)
                                    + n_distpars) + 1):(length(new_params) - (n_weight_pars + n_distpars))] <- c(W, lambdas)
  if(!is.null(B_constraints)) {
    new_B_constraints <- B_constraints[, perm]
  } else {
    new_B_constraints <- NULL
  }

  # Construct the Sstvar model based on the obtained structural parameters
  STVAR(data=stvar$data, p=p, M=M, d=d, params=new_params, weight_function=weight_function,
        weightfun_pars=weightfun_pars, cond_dist=cond_dist, parametrization=stvar$model$parametrization,
        identification=identification, AR_constraints=AR_constraints, mean_constraints=mean_constraints,
        weight_constraints=weight_constraints, B_constraints=new_B_constraints, calc_std_errors=calc_std_errors)
}




#' @title Swap all signs in pointed columns a the \eqn{W} matrix of a structural STVAR model
#'   that is identified by heteroskedasticity.
#'
#' @description \code{swap_W_signs} swaps all signs in pointed columns a the \eqn{W} matrix
#'  a structural STVAR model that is identified by heteroskedasticity. Consequently, signs in
#'   the columns of the impact matrix are also swapped accordingly.
#'
#' @inheritParams reorder_W_columns
#' @param which_to_swap a numeric vector of length at most \eqn{d} and elemnts in \eqn{1,..,d}
#'   specifying the columns of \eqn{W} whose sign should be swapped.
#' @details All signs in any column of \eqn{W} can be swapped without changing the implied reduced form model.
#'   Consequently, also the signs in the columns of the impact matrix are swapped. Note that the sign constraints
#'   imposed on \eqn{W} (or the impact matrix) are also swapped in the corresponding columns accordingly.
#'
#'   Also the order of the columns of \eqn{W} can be changed (without changing the implied reduced
#'   form model) as long as the order of lambda parameters is also changed accordingly. This can be
#'   done with the function \code{reorder_W_columns}.
#' @return Returns an object of class \code{'stvar'} defining a structural STVAR model with the signs of the
#'   columns of \eqn{W} swapped.
#' @seealso \code{\link{GIRF}}, \code{\link{fitSSTVAR}}, \code{\link{reorder_W_columns}}
#' @inherit reorder_W_columns references
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
#' mod222logt2 <- swap_W_signs(mod222logt, which_to_swap=1)
#' mod222logt2 # The signs of the first column of W are swapped
#'
#' # Swap the signs of the second column of W:
#' mod222logt3 <- swap_W_signs(mod222logt, which_to_swap=2)
#' mod222logt3 # The signs of the second column of W are swapped
#'
#' # Swap the signs of both columns of W:
#' mod222logt4 <- swap_W_signs(mod222logt, which_to_swap=1:2)
#' mod222logt4 # The signs of both columns of W are swapped
#' @export

swap_W_signs <- function(stvar, which_to_swap, calc_std_errors=FALSE) {
  check_stvar(stvar)
  if(stvar$model$identification != "heteroskedasticity") stop("Only model identified by heteroskedasticity are supported!")
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  weight_function <- stvar$model$weight_function
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function, weightfun_pars=stvar$model$weightfun_pars)
  cond_dist <- stvar$model$cond_dist
  identification <- stvar$model$identification
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  weight_constraints <- stvar$model$weight_constraints
  B_constraints <- stvar$model$B_constraints
  params <- reform_constrained_pars(p=p, M=M, d=d, params=stvar$params,
                                    weight_function=weight_function, weightfun_pars=weightfun_pars,
                                    cond_dist=cond_dist, identification=identification,
                                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                    weight_constraints=weight_constraints, B_constraints=B_constraints)
  W <- pick_W(p=p, M=M, d=d, params=params, identification=identification)
  lambdas <- pick_lambdas(p=p, M=M, d=d, params=params, identification=identification)
  n_distpars <- length(pick_distpars(params=params, cond_dist=cond_dist))

  # Calculate the number of weight parameters
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

  # Create the new parameter vector
  W[, which_to_swap] <- -W[, which_to_swap]
  W <- Wvec(W) # Zeros removed
  r <- d*(M - 1) # The number of lambda parameters

  new_params <- stvar$params
  new_params[(length(new_params) - (n_weight_pars + length(W) + r
                                    + n_distpars) + 1):(length(new_params) - (n_weight_pars + r + n_distpars))] <- W
  if(!is.null(B_constraints)) {
    new_B_constraints <- B_constraints
    new_B_constraints[, which_to_swap] <- -new_B_constraints[, which_to_swap]
  } else {
    new_B_constraints <- NULL
  }

  # Construct the Sstvar model based on the obtained structural parameters
  STVAR(data=stvar$data, p=p, M=M, d=d, params=new_params, weight_function=weight_function,
        weightfun_pars=weightfun_pars, cond_dist=cond_dist, parametrization=stvar$model$parametrization,
        identification=identification, AR_constraints=AR_constraints, mean_constraints=mean_constraints,
        weight_constraints=weight_constraints, B_constraints=new_B_constraints, calc_std_errors=calc_std_errors)
}
