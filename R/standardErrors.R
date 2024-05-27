#' @title Calculate standard errors for estimates of a smooth transition VAR model
#'
#' @description \code{standard_errors} calculates approximate standard errors for the smooth transition
#'   VAR model using square roots of the diagonal of inverse of observed information matrix
#'   and central-difference approximation for the differentiation.
#'
#' @inheritParams loglikelihood
#' @details This function assumes the standard asymptotic distribution of the estimator
#' @return A vector containing the approximate standard errors of the estimates.
#' @inherit in_paramspace references
#' @keywords internal

standard_errors <- function(data, p, M, params,
                            weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                            weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student"), parametrization=c("intercept", "mean"),
                            identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity"),
                            AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL, minval) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  identification <- match.arg(identification)
  d <- ncol(data)
  check_pMd(p=p, M=M, d=d, weight_function=weight_function, identification=identification)
  weightfun_pars <- check_weightfun_pars(data=data, p=p, M=M, d=d, weight_function=weight_function,
                                         weightfun_pars=weightfun_pars, cond_dist=cond_dist)
  check_constraints(data=data, p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    parametrization=parametrization, identification=identification, AR_constraints=AR_constraints,
                    mean_constraints=mean_constraints, weight_constraints=weight_constraints, B_constraints=B_constraints)
  if(missing(minval)) {
    minval <- get_minval(data)
  }

  # The log-likelihood function to differentiate
  loglik_fn <- function(params) {
    tryCatch(loglikelihood(data=data, p=p, M=M, params=params,
                           weight_function=weight_function, weightfun_pars=weightfun_pars,
                           cond_dist=cond_dist, parametrization=parametrization,
                           identification=identification, AR_constraints=AR_constraints,
                           mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                           B_constraints=B_constraints, check_params=TRUE, to_return="loglik",
                           minval=minval),
             error=function(e) NA)
  }

  # Calculate Hessian
  Hess <- calc_hessian(x=params, fn=loglik_fn, h=6e-6)

  # Inverse of the observed information matrix
  inv_obs_inf <- tryCatch(solve(-Hess), error=function(e) matrix(NA, nrow=length(params), ncol=length(params)))

  # Calculate the standard errors
  unlist(lapply(diag(inv_obs_inf), function(x) ifelse(is.na(x) | x < 0, NA, sqrt(x))))
}
