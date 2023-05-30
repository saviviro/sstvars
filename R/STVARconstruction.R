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
#' @seealso \code{\link{fitSTVAR}}, \code{\link{swap_parametrization}}, \code{\link{GIRF}}
#' @references
#'  \itemize{
#'    \item TO BE FILLED IN
#'  }
#' @examples
#' # p=1, M=1, d=2, linear VAR
#' theta_112relg <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103,
#'   0.601786, -0.002945, 0.067224)
#' mod112 <- STVAR(data=gdpdef, p=1, M=1, params=theta_112relg)
#'
#' # p=1, M=2, d=2, relative dens weight function
#' theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863,
#'   -0.3124, 0.505251, 0.298483, 0.030096, -0.176925, 0.838898, 0.310863, 0.007512,
#'   0.018244, 0.949533, -0.016941, 0.121403, 0.573269)
#' mod122 <- STVAR(data=gdpdef, p=1, M=2, params=theta_122relg)
#' @export

STVAR <- function(data, p, M, d, params, weight_function=c("relative_dens", "logit"), weightfun_pars=NULL,
                  cond_dist=c("Gaussian", "Student"), parametrization=c("intercept", "mean"),
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
  check_pMd(p=p, M=M, d=d)
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars)
  check_constraints(p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
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
                            B_constraints),
                 params=params,
                 std_errors=std_errors,
                 transition_weights=lok_and_tw$tw,
                 regime_cmeans=regime_cmeans,
                 total_cmeans=total_cmeans,
                 total_ccovs=total_ccovs,
                 residuals_raw=residuals_raw,
                 residuals_std=residuals_std,
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
#' # STVAR p=1, M=2 model
#' fit12 <- fitSTVAR(gdpdef, p=1, M=2, nrounds=2, seeds=1:2, ngen=20)
#' fit12
#' fit12_alt <- alt_stvar(fit12, which_largest=2, calc_std_errors=FALSE)
#' fit12_alt # Estimate from an alternative local maximum
#' }
#' @export

alt_stvar <- function(stvar, which_largest=1, which_round, calc_std_errors=TRUE) {
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
