#' @title Calculate residuals of a smooth transition VAR
#'
#' @description \code{get_residuals} calculates residuals of a smooth transition VAR
#'
#' @inheritParams loglikelihood
#' @param standardize standardize the residuals to identity matrix covariance matrix?
#' @return Returns a \eqn{(T \times d)} matrix containing...
#'    \describe{
#'      \item{If \code{standardize == TRUE}:}{the standardized Pearson residuals.}
#'      \item{If \code{standardize == FALSE}:}{the nonstandardized residuals.}
#'    }
#'   Note that the starting time is the start time of data plus \eqn{p + 1},
#'   as the first \eqn{p} observations are used as initial values.
#' @inherit loglikelihood references
#' @keywords internal

get_residuals <- function(data, p, M, params, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold"),
                          weightfun_pars=NULL, cond_dist=c("Gaussian", "Student"), parametrization=c("intercept", "mean"),
                          identification=c("reduced_form", "impact_responses", "heteroskedasticity", "other"),
                          AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL,
                          standardize=TRUE) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  identification <- match.arg(identification)
  T_obs <- nrow(data) - p
  d <- ncol(data)
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                         cond_dist=cond_dist)

  mu_t <- loglikelihood(data=data, p=p, M=M, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
                        cond_dist=cond_dist, parametrization=parametrization, identification=identification,
                        AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                        weight_constraints=weight_constraints, B_constraints=B_constraints,
                        check_params=TRUE, to_return="total_cmeans")

  y_minus_mu <- data[(p + 1):nrow(data),] - mu_t # nonstandardized residuals [T_obs, d]
  if(!standardize) {
    return(y_minus_mu) # Nonstandardized residuals
  }

  Omega_t <- loglikelihood(data=data, p=p, M=M, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
                           cond_dist=cond_dist, parametrization=parametrization, identification=identification,
                           AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                           weight_constraints=weight_constraints, B_constraints=B_constraints,
                           check_params=TRUE, to_return="total_ccovs")

  all_residuals <- matrix(nrow=T_obs, ncol=d)

  # Calculate the Pearson residuals
  for(i1 in 1:T_obs) {
    all_residuals[i1,] <- solve(unvec(d=d, a=get_symmetric_sqrt(Omega_t[, , i1])), y_minus_mu[i1,]) # Standardized Pearson residual
  }
  all_residuals
}
