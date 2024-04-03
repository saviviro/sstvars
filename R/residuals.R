#' @title Calculate residuals of a smooth transition VAR
#'
#' @description \code{get_residuals} calculates residuals of a smooth transition VAR
#'
#' @inheritParams loglikelihood
#' @param standardize standardize the residuals to identity matrix covariance matrix?
#' @param structural_shocks If \code{TRUE}, returns structural shocks instead of residuals
#'  (not available if \code{identification == "reduced_form"}, argument \code{standardize}
#'  is if structural shocks are to be returned).
#' @return Returns a \eqn{(T \times d)} matrix containing...
#'    \describe{
#'      \item{If \code{standardize == TRUE}:}{the standardized Pearson residuals.}
#'      \item{If \code{standardize == FALSE}:}{the nonstandardized residuals.}
#'      \item{If \code{structural_shocks == TRUE}:}{the structural shocks.}
#'    }
#'   Note that the starting time is \eqn{p + 1} counted from the beginning of the starting time of the data,
#'   as the first \eqn{p} observations are used as initial values.
#' @inherit loglikelihood references
#' @keywords internal

get_residuals <- function(data, p, M, params, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                          weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student"), parametrization=c("intercept", "mean"),
                          identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity"),
                          AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL,
                          standardize=TRUE, structural_shocks=FALSE) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  identification <- match.arg(identification)
  if(structural_shocks && identification == "reduced_form" && cond_dist != "ind_Student") {
    stop("Structural shocks are not available if identification == \"reduced_form\" and cond_dist != \"ind_Student\".")
  }
  T_obs <- nrow(data) - p
  d <- ncol(data)
  check_pMd(p=p, M=M, d=d, weight_function=weight_function, identification=identification)
  weightfun_pars <- check_weightfun_pars(data=data, p=p, M=M, d=d, weight_function=weight_function,
                                         weightfun_pars=weightfun_pars, cond_dist=cond_dist)
  check_constraints(data=data, p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    parametrization=parametrization, identification=identification,
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=B_constraints)

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
  transition_weights <- loglikelihood(data=data, p=p, M=M, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                      cond_dist=cond_dist, parametrization=parametrization, identification=identification,
                                      AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                      weight_constraints=weight_constraints, B_constraints=B_constraints,
                                      check_params=TRUE, to_return="tw")

  if(structural_shocks && identification == "heteroskedasticity") { # Obtain W and lambdas
     params_std <- reform_constrained_pars(p=p, M=M, d=d, params=params, weight_function=weight_function,
                                           weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                                           identification=identification, AR_constraints=AR_constraints,
                                           mean_constraints=mean_constraints,weight_constraints=weight_constraints,
                                           B_constraints=B_constraints)
     W <- pick_W(p=p, M=M, d=d, params=params_std, identification=identification)
     lambdas <- pick_lambdas(p=p, M=M, d=d, params=params_std, identification=identification)
     if(M > 1) lambdas <- cbind(1, matrix(lambdas, nrow=d, ncol=M-1)) # First column is column of ones for the first regime
  }
  if(cond_dist == "ind_Student") {
    all_Omegas <- pick_Omegas(p=p, M=M, d=d, params=reform_constrained_pars(p=p, M=M, d=d, params=params,
                                                                            weight_function=weight_function,
                                                                            weightfun_pars=weightfun_pars,
                                                                            cond_dist=cond_dist,
                                                                            identification=identification,
                                                                            AR_constraints=AR_constraints,
                                                                            mean_constraints=mean_constraints,
                                                                            weight_constraints=weight_constraints,
                                                                            B_constraints=B_constraints),
                              cond_dist=cond_dist, identification=identification)
    alpha_mt <- loglikelihood(data=data, p=p, M=M, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
                              cond_dist=cond_dist, parametrization=parametrization, identification=identification,
                              AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                              weight_constraints=weight_constraints, B_constraints=B_constraints,
                              check_params=TRUE, to_return="tw")
    all_Bt <- get_Bt_Cpp(all_Omegas=all_Omegas, alpha_mt=alpha_mt)
  }

  # Go through each point of time and calculate the residuals/shocks
  for(i1 in 1:T_obs) {
    if(cond_dist == "ind_Student") {
      # Structural shocks and standardized errors are the same thing here due to statistical identification
      all_residuals[i1,] <- solve(all_Bt[, , i1], y_minus_mu[i1, ])
    } else {
      if(structural_shocks) { # Structural shock:
        # Calculate the impact matrix
        if(identification == "recursive") {
          B_t <- t(chol(Omega_t[, , i1])) # Lower triangular Cholesky decomposition of Omega_t
        } else if(identification == "heteroskedasticity") {
          if(M == 1) {
            B_t <- W
          } else {
            tmp <- array(dim=c(d, d, M))
            for(m in 1:M) {
              tmp[, , m] <- transition_weights[i1, m]*diag(lambdas[, m])
            }
            B_t <- W%*%sqrt(apply(tmp, MARGIN=1:2, FUN=sum))
          }
        }
        # Recover the structural shock
        all_residuals[i1,] <- solve(B_t, y_minus_mu[i1, ])
      } else { # Standardized Pearson residual:
        all_residuals[i1,] <- solve(unvec(d=d, a=get_symmetric_sqrt(Omega_t[, , i1])), y_minus_mu[i1,])
      }
    }
  }
  all_residuals
}


