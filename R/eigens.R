#' @title Calculate absolute values of the eigenvalues of the "bold A" matrices containing the AR coefficients
#'
#' @description \code{get_boldA_eigens} calculates absolute values of the eigenvalues of
#'   the "bold A" matrices containing the AR coefficients for each regime.
#'
#' @param stvar object of class \code{"stvar"}
#' @return Returns a matrix with \eqn{d*p} rows and \eqn{M} columns - one column for each regime.
#'  The \eqn{m}th column contains the absolute values (or modulus) of the eigenvalues of the "bold A" matrix containing
#'  the AR coefficients correspinding to regime \eqn{m}.
#' @inherit stab_conds_satisfied references
#' @keywords internal

get_boldA_eigens <- function(stvar) {
  check_stvar(stvar)
  get_boldA_eigens_par(p=stvar$model$p, M=stvar$model$M, d=stvar$model$d, params=stvar$params,
                       weight_function=stvar$model$weight_function, weightfun_pars=stvar$model$weightfun_pars,
                       cond_dist=stvar$model$cond_dist, parametrization=stvar$model$parametrization,
                       identification=stvar$model$identification, AR_constraints=stvar$model$AR_constraints,
                       mean_constraints=stvar$model$mean_constraints, weight_constraints=stvar$model$weight_constraints,
                       B_constraints=stvar$model$B_constraints)
}


#' @title Calculate the eigenvalues of the "Omega" error term covariance matrices
#'
#' @description \code{get_omega_eigens} calculates the eigenvalues of the "Omega" error
#'  term covariance matrices for each regime
#'
#' @inheritParams get_boldA_eigens
#' @return Returns a matrix with \eqn{d} rows and \eqn{M} columns - one column for each regime.
#'  The \eqn{m}th column contains the eigenvalues of the "Omega" error term covariance matrix
#'  of the \eqn{m}th regime.
#' @inherit get_boldA_eigens references
#' @keywords internal

get_omega_eigens <- function(stvar) {
  check_stvar(stvar)
  get_omega_eigens_par(p=stvar$model$p, M=stvar$model$M, d=stvar$model$d, params=stvar$params,
                       weight_function=stvar$model$weight_function, weightfun_pars=stvar$model$weightfun_pars,
                       cond_dist=stvar$model$cond_dist, parametrization=stvar$model$parametrization,
                       identification=stvar$model$identification, AR_constraints=stvar$model$AR_constraints,
                       mean_constraints=stvar$model$mean_constraints, weight_constraints=stvar$model$weight_constraints,
                       B_constraints=stvar$model$B_constraints)
}


#' @title Calculate absolute values of the eigenvalues of the "bold A" matrices containing the AR coefficients
#'
#' @description \code{get_boldA_eigens_par} calculates absolute values of the eigenvalues of
#'   the "bold A" matrices containing the AR coefficients for each regime.
#'
#' @inheritParams loglikelihood
#' @return Returns a matrix with \eqn{d*p} rows and \eqn{M} columns - one column for each regime.
#'  The \eqn{m}th column contains the absolute values (or modulus) of the eigenvalues of the "bold A" matrix containing
#'  the AR coefficients corresponding to regime \eqn{m}.
#' @inherit stab_conds_satisfied references
#' @keywords internal

get_boldA_eigens_par <- function(p, M, d, params,
                                 weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                                 weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student"), parametrization=c("intercept", "mean"),
                                 identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity"),
                                 AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL) {
  parametrization <- match.arg(parametrization)
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  identification <- match.arg(identification)
  weightfun_pars <- check_weightfun_pars(p=p, d=d, M=M, weight_function=weight_function,
                                         weightfun_pars=weightfun_pars, cond_dist=cond_dist)
  check_constraints(p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    parametrization=parametrization, identification=identification,
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=B_constraints)
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params,
                                    weight_function=weight_function, weightfun_pars=weightfun_pars,
                                    cond_dist=cond_dist, identification=identification,
                                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                    weight_constraints=weight_constraints, B_constraints=B_constraints)
  all_A <- pick_allA(p=p, M=M, d=d, params=params)
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  matrix(vapply(1:M, function(m) abs(eigen(all_boldA[, , m], symmetric=FALSE, only.values=TRUE)$'values'), numeric(d*p)),
         nrow=d*p, ncol=M, byrow=FALSE)
}



#' @title Calculate the eigenvalues of the "Omega" error term covariance matrices
#'
#' @description \code{get_omega_eigens_par} calculates the eigenvalues of the "Omega" error
#'  term covariance matrices for each regime
#'
#' @inheritParams loglikelihood
#' @return Returns a matrix with \eqn{d} rows and \eqn{M} columns - one column for each regime.
#'  The \eqn{m}th column contains the eigenvalues of the "Omega" error term covariance matrix
#'  of the \eqn{m}th regime.
#' @inherit get_boldA_eigens references
#' @keywords internal

get_omega_eigens_par <- function(p, M, d, params,
                                 weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                                 weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student"), parametrization=c("intercept", "mean"),
                                 identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity"),
                                 AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL) {
  parametrization <- match.arg(parametrization)
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  identification <- match.arg(identification)
  weightfun_pars <- check_weightfun_pars(p=p, d=d, M=M, weight_function=weight_function,
                                         weightfun_pars=weightfun_pars, cond_dist=cond_dist)
  check_constraints(p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    parametrization=parametrization, identification=identification,
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=B_constraints)
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params,
                                    weight_function=weight_function, weightfun_pars=weightfun_pars,
                                    cond_dist=cond_dist, identification=identification,
                                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                    weight_constraints=weight_constraints, B_constraints=B_constraints)
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, cond_dist=cond_dist, identification=identification)
  if(cond_dist == "ind_Student") {
    all_Omega <- array(vapply(1:M, function(m) tcrossprod(all_Omega[,,m]), numeric(d^2)), dim=c(d, d, M))
  }
  matrix(vapply(1:M, function(m) eigen(all_Omega[, , m], symmetric=TRUE, only.values=TRUE)$'values', numeric(d)),
           nrow=d, ncol=M, byrow=FALSE)
}


#' @title Warn about near-unit-roots in some regimes
#'
#' @description \code{warn_eigens} warns if the model contains near-unit-roots in some regimes
#'
#' @inheritParams get_boldA_eigens
#' @param tol if eigenvalue is closer than \code{tol} to its bound, a warning is thrown
#' @details Warns if, for some regime, some moduli of "bold A" eigenvalues are larger than \code{1 - tol} or
#'  some eigenvalue of the error term covariance matrix is smaller than \code{tol}.
#' @return Doesn't return anything.
#' @keywords internal

warn_eigens <- function(stvar, tol=0.002) {
  boldA_eigens <- get_boldA_eigens(stvar)
  omega_eigens <- get_omega_eigens(stvar)
  M <- stvar$model$M
  near_nonstat <- vapply(1:M, function(i1) any(abs(boldA_eigens[,i1]) > 1 - tol), logical(1))
  near_singular <- vapply(1:M, function(i1) any(abs(omega_eigens[,i1]) < tol), logical(1))
  if(any(near_nonstat)) {
    my_string1 <- ifelse(sum(near_nonstat) == 1,
                         paste("Regime", which(near_nonstat),"has near-unit-roots! "),
                         paste("Regimes", paste(which(near_nonstat), collapse=" and ") ,"have near-unit-roots! "))
  } else {
    my_string1 <- NULL
  }
  if(any(near_singular)) {
    my_string2 <- ifelse(sum(near_singular) == 1,
                         paste("Regime", which(near_singular),"has near-singular error term covariance matrix! "),
                         paste("Regimes", paste(which(near_singular), collapse=" and "),
                               "have near-singular error term covariance matrices! "))
  } else {
    my_string2 <- NULL
  }
  if(any(near_nonstat) || any(near_singular)) {
    warning(paste0(my_string1, my_string2,
                   "Consider building a model from the next-largest local maximum with the function 'alt_stvar'",
                   "by adjusting its argument 'which_largest'."))
  }
}

