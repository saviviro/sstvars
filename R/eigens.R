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
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  params <- stvar$params
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  # REFORM CONSTRAINED PARS OR IN _PAR-FUNCTION MAYBE
  if(!is.null(AR_constraints) || !is.null(mean_constraints)) {
    stop("Constrained models are not yet implemented to get_boldA_eigens")
  }
  get_boldA_eigens_par(p=p, M=M, d=d, params=params, AR_constraints=AR_constraints, mean_constraints=mean_constraints)
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
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  params <- stvar$params
  # REFORM CONSTRAINED PARS
  if(!is.null(stvar$model$AR_constraints) || !is.null(stvar$model$mean_constraints)) {
    stop("Constrained models are not yet implemented to get_omega_eigens")
  }
  if(stvar$model$identification != "reduced_form") {
    stop("Structural models not yet implemented to get_omega_eigens")
  }
  get_omega_eigens_par(p=p, M=M, d=d, params=params, identification=identification,
                       AR_constraints=AR_constraints, mean_constraints=mean_constraints)
}


#' @title Calculate absolute values of the eigenvalues of the "bold A" matrices containing the AR coefficients
#'
#' @description \code{get_boldA_eigens_par} calculates absolute values of the eigenvalues of
#'   the "bold A" matrices containing the AR coefficients for each regime.
#'
#' @inheritParams loglikelihood
#' @return Returns a matrix with \eqn{d*p} rows and \eqn{M} columns - one column for each regime.
#'  The \eqn{m}th column contains the absolute values (or modulus) of the eigenvalues of the "bold A" matrix containing
#'  the AR coefficients correspinding to regime \eqn{m}.
#' @inherit stab_conds_satisfied references
#' @keywords internal

get_boldA_eigens_par <- function(p, M, d, params, AR_constraints=NULL, mean_constraints=NULL) {
  # REFORM CONSTRAINED PARS
  if(!is.null(AR_constraints) || !is.null(mean_constraints)) {
    stop("Constrained models are not yet implemented to get_boldA_eigens_par")
  }
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

get_omega_eigens_par <- function(p, M, d, params, identification=c("reduced_form", "recursive", "heteroskedasticity"),
                                 AR_constraints=NULL, mean_constraints=NULL) {
  identification <- match.arg(identification)
  # REFORM CONSTRAINED PARS
  if(!is.null(AR_constraints) || !is.null(mean_constraints)) {
    stop("Constrained models are not yet implemented to get_omega_eigens_par")
  }
  if(identification != "reduced_form") {
    stop("Structural models not yet implemented to get_omega_eigens_par")
  }
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params)
  matrix(vapply(1:M, function(m) eigen(all_Omega[, , m], symmetric=TRUE, only.values=TRUE)$'values', numeric(d)),
         nrow=d, ncol=M, byrow=FALSE)
}
