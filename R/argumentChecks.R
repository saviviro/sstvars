#' @title Check the stability condition for each of the regimes
#'
#' @description \code{stab_conds_satisfied} checks whether the stability condition is satisfied
#'   for each of the regimes.
#'
#' @inheritParams pick_allA
#' @param all_boldA 3D array containing the \eqn{((dp)x(dp))} "bold A" (companion form) matrices of each regime,
#'   obtained from \code{form_boldA}. Will be computed if not given.
#' @param tolerance Returns \code{FALSE} if modulus of any eigenvalue of "bold A" is larger or equal to \code{1-tolerance}.
#' @inherit pick_phi0 details
#' @return Returns \code{TRUE} if the stability condition is satisfied for all regimes and \code{FALSE} if not.
#'   According to the argument \code{tolerance}, \code{stab_conds_satisfied} may return \code{FALSE} when the parameter
#'   vector satisfies the stability conditions but is very close to the boundary (this is used to ensure numerical stability
#'   in the estimation of the model parameters).
#' @section Warning:
#'  No argument checks!
#' @inherit pick_phi0 details
#' @inherit form_boldA references
#' @keywords internal

stab_conds_satisfied <- function(p, M, d, params, all_boldA=NULL, tolerance=1e-3) {
  if(is.null(all_boldA)) {
    all_A <- pick_allA(p=p, M=M, d=d, params=params)
    all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  }
  for(m in 1:M) {
    if(any(abs(eigen(all_boldA[, , m], symmetric=FALSE, only.values=TRUE)$'values') >= 1 - tolerance)) {
      return(FALSE)
    }
  }
  TRUE
}



#' @title Determine whether the parameter vector is in the parameter space
#'
#' @description \code{in_paramspace} checks whether the parameter vector is in the parameter
#'   space.
#'
#' @inheritParams loglikelihood
#' @inheritParams stab_conds_satisfied
#' @param weightpars numerical vector containing the transition weight parameters, obtained from \code{pick_weightpars}.
#' @param all_Omega 3D array containing all covariance matrices \eqn{\Omega_{m}}, obtained from \code{pick_Omegas}.
#' @param stab_tol numerical tolerance for the stability condition of each regime: if the "bold A" matrix of any regime
#'   has eigenvalues larger that \code{1 - stab_tol} the parameter vector will be labeled as not in the parameter space.
#'   Note that if the tolerance is too small, numerical evaluation of the log-likelihood might fail.
#' @param posdef_tol numerical tolerance for positive definiteness of the regime-specific covariance matrices: if
#'   the error term covariance matrix of any regime has eigenvalues smaller than this, the model is classified
#'   as not satisfying positive definiteness assumption. Note that if the tolerance is too small, numerical
#'   evaluation of the log-likelihood might fail and cause error.
#' @param df_tol the parameter vector is considered to be outside the parameter space the degrees of
#'   freedom parameters is not larger than \code{2 + df_tol}.
#' @details The parameter vector in the argument \code{params} should be unconstrained and reduced form.
#' @return Returns \code{TRUE} if the given parameter values are in the parameter space and \code{FALSE} otherwise.
#'   This function does NOT consider the identifiability condition!
#' @references
#'  \itemize{
#'    \item TO BE FILLED IN
#'  }
#'  @keywords internal

in_paramspace <- function(p, M, d, weight_function, cond_dist, all_boldA, all_Omegas, weightpars, df,
                          stab_tol=1e-3, posdef_tol=1e-8, df_tol=1e-8) {

  if(cond_dist == "Student") { # Check degrees of freedom parameter
    if(df <= 2 + df_tol) {
      return(FALSE)
    }
  }
  if(weight_function == "relative_dens") {
    if(M >= 2 & sum(weightpars[-M]) >= 1) {
      return(FALSE)
    } else if(any(weightpars <= 0)) {
      return(FALSE)
    }
  }
  if(!stab_conds_satisfied(p=p, M=M, d=d, all_boldA=all_boldA, tolerance=stab_tol)) {
    return(FALSE)
  }
  for(m in 1:M) {
    if(any(eigen(all_Omegas[, , m], symmetric=TRUE, only.values=TRUE)$values < posdef_tol)) {
      return(FALSE)
    }
  }
  TRUE
}


#' @title Check whether the parameter vector is in the parameter space and throw error if not
#'
#' @description \code{check_parameters} checks whether the parameter vector is in the parameter
#'   space.
#'
#' @inheritParams loglikelihood
#' @inheritParams in_paramspace
#' @return Throws an informative error if there is something wrong with the parameter vector.
#' @inherit in_paramspace references
#' @examples
#'  \dontrun{
#'  # There examples will cause an informative error
#'  params112relg_notpd <- c(6.5e-01, 7.0e-01, 2.9e-01, 2.0e-02, -1.4e-01,
#'   9.0e-01, 6.0e-01, -1.0e-02, 1.0e-07)
#'
#'  # FILL IN!
#'  }
#' @export

check_parameters <- function(p, M, d, params, weight_function=c("relative_dens", "logit"), cond_dist=c("Gaussian", "Student"),
                             parametrization=c("intercept", "mean"),
                             identification=c("reduced_form", "recursive", "heteroskedasticity"),
                             AR_constraints=NULL, mean_constraints=NULL, B_constraints=NULL,
                             stab_tol=1e-3, posdef_tol=1e-8, df_tol=1e-8) {

  # FILL IN!
  # FILL IN!
}


#' @title Check whether all arguments are positive integers
#'
#' @description \code{all_pos_ints} checks whether all the elements in a vector
#'   are positive integers.
#'
#' @param x a vector containing the elements to be tested.
#' @return Returns \code{TRUE} or \code{FALSE} accordingly.
#' @keywords internal

all_pos_ints <- function(x) {
  all(vapply(x, function(x1) x1 %% 1 == 0 && length(x1) == 1 && x1 >= 1, logical(1)))
}

#' @title Check that p, M, and d are correctly set
#'
#' @description \code{check_pMd} checks the arguments p, M, and d.
#'
#' @inheritParams stab_conds_satisfied
#' @return Throws an error if something is wrong.
#' @keywords internal

check_pMd <- function(p, M, d) {
  if(!all_pos_ints(M) || length(M) != 1) {
    stop("The argument M must be a positive integer!")
  }
  if(!all_pos_ints(p) || length(p) != 1) {
    stop("The argument p must be a positive integer!")
  }
  if(!missing(d)) {
    if(d < 2 | d%%1 != 0) {
      stop("The argument d, the number of columns in the data matrix, has to be a positive integer larger than one!")
    }
  }
}
