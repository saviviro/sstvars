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
#'    \item TO BE FILLED IN
#'  }
#' @examples
#' # p=1, M=1, d=2, linear VAR
#' theta_112relg <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103,
#'   0.601786, -0.002945, 0.067224)
#' mod112 <- STVAR(data=gdpdef, p=1, M=1, params=theta_112relg)
#' summary(mod112)
#'
#' # p=1, M=2, d=2, relative dens weight function
#' theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863,
#'   -0.3124, 0.505251, 0.298483, 0.030096, -0.176925, 0.838898, 0.310863, 0.007512,
#'   0.018244, 0.949533, -0.016941, 0.121403, 0.573269)
#' mod122 <- STVAR(data=gdpdef, p=1, M=2, params=theta_122relg)
#' summary(mod122)
#' plot(mod122)
#' plot(mod122, plot_type="cond_mean")
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
#' ## p=1, M=1, d=2, linear VAR
#' theta_112relg <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103,
#'   0.601786, -0.002945, 0.067224)
#' mod112 <- STVAR(data=gdpdef, p=1, M=1, params=theta_112relg, calc_std_errors=TRUE)
#' print_std_errors(mod112) # Standard errors for the intercepts
#'
#' # Swap to mean parametrization:
#' mod112mu <- swap_parametrization(mod112, calc_std_errors=TRUE)
#' print_std_errors(mod112mu) # Standard errors for the stationary means
#'
#' ## p=1, M=2, d=2, relative dens weight function; not calculating standard errors:
#' theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863,
#'   -0.3124, 0.505251, 0.298483, 0.030096, -0.176925, 0.838898, 0.310863, 0.007512,
#'   0.018244, 0.949533, -0.016941, 0.121403, 0.573269)
#' mod122 <- STVAR(data=gdpdef, p=1, M=2, params=theta_122relg, calc_std_errors=FALSE)
#' mod122mu <- swap_parametrization(mod122, calc_std_errors=FALSE) # mean parametrization
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
#' @inherit STVAR references
#' @examples
#' # p=2, M=2, d=2, Student's t logistic STVAR model with the first lag of the second
#'  # variable  as the switching variable and shocks identified by heteroskedasticity.
#' theta_222logt <- c(0.356914, 0.107436, 0.356386, 0.086330, 0.139960, 0.035172, -0.164575,
#'   0.386816, 0.451675, 0.013086, 0.227882, 0.336084, 0.239257, 0.024173, -0.021209, 0.707502,
#'   0.063322, 0.027287, 0.009182, 0.197066, -0.03, 0.24, -0.76, -0.02, 3.36, 0.86, 0.1, 0.2, 7)
#' mod222logt <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222logt,
#'                                weight_function="logistic", weightfun_pars=c(2, 1),
#'                                cond_dist="Student", identification="heteroskedasticity")
#' mod222logt
#'
#' # Reverse the ordering of the columns of W:
#' mod222logt_rev <- reorder_W_columns(mod222logt, perm=c(2, 1))
#' mod222logt_rev
#'
#' # Swap the ordering back to the original:
#' mod222logt_rev2 <- reorder_W_columns(mod222logt_rev, perm=c(2, 1))
#' mod222logt_rev2
#'
#' # Does not do anything, as perm=1:2 so the ordering does not change:
#' reorder_W_columns(mod222logt, perm=c(1, 2))
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
  distpars <- pick_distpars(params=params, cond_dist=cond_dist)

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
                                    + length(distpars)) + 1):(length(new_params) - (n_weight_pars + length(distpars)))] <- c(W, lambdas)
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
