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
  # in_paramspace is internal function that always takes in non-constrained reduced form parameter vector
  # Reform the parameter vectors before checking with in_paramspace

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
  } else {
    stop("Other weight functions are not yet implemented!")
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
#' @description \code{check_params} checks whether the parameter vector is in the parameter
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
#'  check_params(p=1, M=1, d=2, params=params112relg_notpd)
#'  }
#' @export

check_params <- function(p, M, d, params, weight_function=c("relative_dens", "logit"), cond_dist=c("Gaussian", "Student"),
                         parametrization=c("intercept", "mean"),
                         identification=c("reduced_form", "recursive", "heteroskedasticity"),
                         AR_constraints=NULL, mean_constraints=NULL, B_constraints=NULL,
                         stab_tol=1e-3, posdef_tol=1e-8, df_tol=1e-8) {
  check_pMd(p=p, M=M, d=d)
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  identification <- match.arg(identification)
  if(n_params(p=p, M=M, d=d, weight_function=weight_function, cond_dist=cond_dist,
              identification=identification, AR_constraints=AR_constraints,
              mean_constraints=mean_constraints, B_constraints=B_constraints) != length(params)) {
    stop("The parameter vector has wrong dimension!")
  }
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, weight_function=weight_function, cond_dist=cond_dist,
                                    identification=identification, AR_constraints=AR_constraints,
                                    mean_constraints=mean_constraints, B_constraints=B_constraints)

  # Pick params
  all_phi0 <- pick_phi0(M=M, d=d, params=params) # phi0 or mean parameters
  all_A <- pick_allA(p=p, M=M, d=d, params=params) # [d, d, p, M]
  all_Omegas <- pick_Omegas(p=p, M=M, d=d, params=params) # [d, d, M]
  weightpars <- pick_weightpars(p=p, M=M, d=d, params=params, weight_function=weight_function, cond_dist=cond_dist)
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  df <- numeric(0) # FILL IN WHEN STUDENT IS IMPLEMENTED
  if(!is.null(B_constraints)) {
    stop("B_constraints are not yet implemented to check_params!")
  }
  if(identification != "reduced_form") {
    stop("Only reduced form models are currently implemented to check_params!")
  }

  if(cond_dist == "Student") { # Check degrees of freedom parameter
    stop("cond_dist = Student is not yet implented to check_params!")
    if(df <= 2 + df_tol) {
      stop("The degrees of freedom parameter needs to be strictly larger than two (with large enough numerical tolerance)!")
    }
  }
  if(weight_function == "relative_dens") {
    if(M >= 2 & sum(weightpars[-M]) >= 1) {
      stop("The transition weight parameter alphas must sum to one!")
    } else if(any(weightpars <= 0)) {
      stop("The transition weight parameter alphas must be strictly larger than zero!")
    }
  } else {
    stop("Other weight functions are not yet implemented!")
  }
  if(!stab_conds_satisfied(p=p, M=M, d=d, all_boldA=all_boldA, tolerance=stab_tol)) {
    stop("At least one of the regimes does not satisfy the stability condition (with large enough numerical tolerance)!")
  }
  for(m in 1:M) {
    if(any(eigen(all_Omegas[, , m], symmetric=TRUE, only.values=TRUE)$values < posdef_tol)) {
      stop(paste0("The conditional covariance matrix of Regime ", m, " is not positive definite (with large enough numerical tolerance)!"))
    }
  }
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


#' @title Check the data is in the correct form
#'
#' @description \code{check_data} checks the data.
#'
#' @inheritParams loglikelihood
#' @return Checks the data and tries to correct it. Throws an error if something is wrong and
#'   returns the corrected data otherwise.
#' @keywords internal

check_data <- function(data, p) {
  if(is.data.frame(data)) {
    data <- as.matrix(data)
  }
  if(!is.matrix(data)) {
    stop("The data must be numeric matrix (possibly a class 'ts' object)!")
  } else {
    if(anyNA(data)) stop("The data contains NA values!")
    if(!is.numeric(data)) stop("The data must be numeric!")
    if(ncol(data) < 2) stop("The data matrix must contain at least two columns!")
    if(nrow(data) < p + 1) stop("The data must contain at least p+1 observations!")
  }
  data
}


#' @title Calculate the number of (freely estimaed) parameters in the model
#'
#' @description \code{n_params} calculates the number of (freely estimaed) parameters in the model.
#'
#' @inheritParams check_params
#' @return Returns the number of parameters in the parameter vector of the specified model.
#' @section Warning:
#'  No argument checks!
#' @inherit in_paramspace references
#' @keywords internal

n_params <- function(p, M, d, weight_function=c("relative_dens", "logit"), cond_dist=c("Gaussian", "Student"),
                     identification=c("reduced_form", "recursive", "heteroskedasticity"),
                     AR_constraints=NULL, mean_constraints=NULL, B_constraints=NULL) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  identification <- match.arg(identification)
  if(identification != "reduced_form") stop("Only reduced form models are currently supported")

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
  if(is.null(B_constraints)) {
    n_covmat_pars <- M*d*(d + 1)/2
  } else { # Constraints on the impact matrix
    stop("B_constraints not yet implemented to n_parms")
    n_covmat_pars <- NULL
  }
  if(weight_function == "relative_dens") {
    n_weight_pars <- M - 1
  } else { # weight_function == "logit"
    stop("only relative_dens weight fn is implemented to n_params")
    n_weight_pars <- NULL
  }
  if(cond_dist == "Gaussian") {
    n_dist_pars <- 0
  } else { # cond_dist == "Student"
    n_dist_pars <- 1 # degrees of freedom param
  }
  n_mean_pars + n_ar_pars + n_covmat_pars + n_weight_pars + n_dist_pars
}
