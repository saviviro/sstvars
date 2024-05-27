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
#' @param all_Omegas A 3D array containing the covariance matrix parameters obtain from \code{pick_Omegas}...
#'   \describe{
#'     \item{If \code{cond_dist \%in\% c("Gaussian", "Student")}:}{all covariance matrices \eqn{\Omega_{m}} in \code{[, , m]}.}
#'     \item{If \code{cond_dist=="ind_Student"}:}{all impact matrices \eqn{B_m} of the regimes in \code{[, , m]}.}
#'   }
#' @param distpars A numeric vector containing the distribution parameters...
#'   \describe{
#'     \item{If \code{cond_dist=="Gaussian"}:}{Not used, i.e., a numeric vector of length zero.}
#'     \item{If \code{cond_dist=="Student"}:}{The degrees of freedom parameter, i.e., a numeric vector of length one.}
#'   }
#' @param transition_weights (optional; only for models with \code{cond_dist="ind_Student"} or \code{identification="non-Gaussianity"})
#'   A \eqn{T \times M} matrix containing the transition weights. If \code{cond_dist="ind_Student"} checks that the impact matrix
#'   \eqn{\sum_{m=1}^M\alpha_{m,t}^{1/2}B_m} is invertible for all \eqn{t=1,...,T}.
#' @details The parameter vector in the argument \code{params} should be unconstrained and reduced form.
#' @return Returns \code{TRUE} if the given parameter values are in the parameter space and \code{FALSE} otherwise.
#'   This function does NOT consider identification conditions!
#' @references
#'  \itemize{
#'    \item Kheifets I.L., Saikkonen P.J. 2020. Stationarity and ergodicity of Vector STAR models.
#'      \emph{Econometric Reviews}, \strong{39}:4, 407-414.
#'    \item Lütkepohl H. 2005. New Introduction to Multiple Time Series Analysis, \emph{Springer}.
#'    \item Lanne M., Virolainen S. 2024. A Gaussian smooth transition vector autoregressive model:
#'       An application to the macroeconomic effects of severe weather shocks. Unpublished working
#'       paper, available as arXiv:2403.14216.
#'    \item Virolainen S. 2024. Identification by non-Gaussianity in structural threshold and
#'       smooth transition vector autoregressive models. Unpublished working
#'       paper, available as arXiv:2404.19707.
#'  }
#'  @keywords internal

in_paramspace <- function(p, M, d, params,
                          weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                          weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student"),
                          identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity"),
                          B_constraints=NULL, other_constraints=NULL, all_boldA, all_Omegas, weightpars, distpars,
                          transition_weights, stab_tol=1e-3, posdef_tol=1e-8, distpar_tol=1e-8, weightpar_tol=1e-8) {
  # in_paramspace is internal function that always takes in non-constrained reduced form parameter vector
  # Reform the parameter vectors before checking with in_paramspace
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  identification <- match.arg(identification)

  # Check distribution parameters
  if(cond_dist == "Student") {
    if(distpars <= 2 + distpar_tol) {
      return(FALSE)
    }
  } else if(cond_dist == "ind_Student") {
    if(any(distpars <= 2 + distpar_tol)) {
      return(FALSE)
    }
  }

  # Check weight_function parameters
  if(weight_function == "relative_dens") {
    if(M >= 2 & sum(weightpars[-M]) >= 1) {
      return(FALSE)
    } else if(any(weightpars <= 0)) {
      return(FALSE)
    }
  } else if(weight_function == "logistic" | weight_function == "exponential") {
    if(weightpars[2] <= 0 + weightpar_tol) {
      return(FALSE) # We assume strictly positive scale parameter gamma
    }
  } else if(weight_function == "mlogit") {
    # All real numbers are ok, so nothing to check
  } else if(weight_function == "threshold") {
    if(!all(order(weightpars, decreasing=FALSE) == seq_len(M - 1))) {
      return(FALSE) # We assume increasing ordering of the thresholds
    }
  } else if(weight_function == "exogenous") {
    # No weightpars to test
  }

  # Check stability conditions of a linear VAR
  if(!stab_conds_satisfied(p=p, M=M, d=d, all_boldA=all_boldA, tolerance=stab_tol)) {
    return(FALSE)
  }

  # Check positive definiteness of the covariance matrices / invertibility of the impact matrices and matrix
  if(cond_dist == "Gaussian" || cond_dist == "Student") {
    for(m in 1:M) {
      if(any(eigen(all_Omegas[, , m], symmetric=TRUE, only.values=TRUE)$values < posdef_tol)) {
        return(FALSE)
      }
    }
  } else { # cond_dist == "ind_Student"
    for(m in 1:M) { # Check the invertibility of the impact matrices of the regimes
      if(abs(det(all_Omegas[, , m])) < posdef_tol) {
        return(FALSE)
      }
    }
    if(!missing(transition_weights)) { # Check the invertibility of the time-varying impact matrix for all t
      for(i1 in 1:nrow(transition_weights)) { # Impact matrices in all_Omegas
        if(abs(det(apply(X=all_Omegas, MARGIN=c(1, 2), FUN=function(mat) sum(mat*sqrt(transition_weights[i1,]))))) < posdef_tol) {
          return(FALSE)
        }
      }
    }
  }
  # Check that sign constraints in B_constraints are satisfied, if imposed
  if(!is.null(B_constraints)) {
    if(identification == "non-Gaussianity" || cond_dist == "ind_Student") {
      for(m in 1:M) {
        if(any(all_Omegas[, , m][B_constraints > 0] <= 0, na.rm=TRUE) || any(all_Omegas[B_constraints < 0] >= 0, na.rm=TRUE)) {
          return(FALSE)
        }
      }
    } else if(identification == "heteroskedasticity") {
      W <- pick_W(p=p, M=M, d=d, params=params, identification=identification)
      if(any(W[B_constraints > 0] <= 0, na.rm=TRUE) || any(W[B_constraints < 0] >= 0, na.rm=TRUE)) {
        return(FALSE)
      }
    }
  }
  # Check other constraints if applied
  if(identification == "non-Gaussianity" || cond_dist == "ind_Student") {
    if(!is.null(other_constraints$B1_constraints)) {
      # The first non-zero entry in each column of B_1 is constrained strictly positive and they
      # are constrained to be in a decreasing ordering.
      first_non_zero <- vapply(1:d, function(i1) which(all_Omegas[, i1, 1] != 0)[1], numeric(1))
      first_non_zero_entries <- vapply(1:d, function(i1) all_Omegas[first_non_zero[i1], i1, 1], numeric(1))
      if(any(first_non_zero_entries <= 0)) {
        return(FALSE)
      } else if(!all(order(first_non_zero_entries, decreasing=TRUE) == seq_len(d))) {
        return(FALSE)
      }
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
#'  # There examples will cause an informative error
#'  params112_notpd <- c(6.5e-01, 7.0e-01, 2.9e-01, 2.0e-02, -1.4e-01,
#'   9.0e-01, 6.0e-01, -1.0e-02, 1.0e-07)
#'  try(check_params(p=1, M=1, d=2, params=params112_notpd))
#'
#'  params112_notstat <- c(6.5e-01, 7.0e-01, 10.9e-01, 2.0e-02, -1.4e-01,
#'   9.0e-01, 6.0e-01, -1.0e-02, 1.0e-07)
#'  try(check_params(p=1, M=1, d=2, params=params112_notstat))
#'
#'  params112_wronglength <- c(6.5e-01, 7.0e-01, 2.9e-01, 2.0e-02, -1.4e-01,
#'   9.0e-01, 6.0e-01, -1.0e-02)
#'  try(check_params(p=1, M=1, d=2, params=params112_wronglength))
#' @export

check_params <- function(data, p, M, d, params, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                         weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student"), parametrization=c("intercept", "mean"),
                         identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity"),
                         AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL, transition_weights,
                         stab_tol=1e-3, posdef_tol=1e-8, distpar_tol=1e-8, weightpar_tol=1e-8) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  identification <- match.arg(identification)
  check_pMd(p=p, M=M, d=d, weight_function=weight_function, identification=identification)
  weightfun_pars <- check_weightfun_pars(data=data, p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                         cond_dist=cond_dist)
  if(n_params(p=p, M=M, d=d, weight_function=weight_function, cond_dist=cond_dist,
              identification=identification, AR_constraints=AR_constraints,
              mean_constraints=mean_constraints, weight_constraints=weight_constraints,
              B_constraints=B_constraints, weightfun_pars=weightfun_pars) != length(params)) {
    stop("The parameter vector has wrong dimension!")
  }
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, weight_function=weight_function, cond_dist=cond_dist,
                                    identification=identification, AR_constraints=AR_constraints,
                                    mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                    B_constraints=B_constraints, weightfun_pars=weightfun_pars)

  # Pick params
  all_phi0 <- pick_phi0(M=M, d=d, params=params) # phi0 or mean parameters
  all_A <- pick_allA(p=p, M=M, d=d, params=params) # [d, d, p, M]
  all_Omegas <- pick_Omegas(p=p, M=M, d=d, params=params, cond_dist=cond_dist, identification=identification) # [d, d, M]
  weightpars <- pick_weightpars(p=p, M=M, d=d, params=params, weight_function=weight_function, cond_dist=cond_dist,
                                weightfun_pars=weightfun_pars)
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  distpars <- pick_distpars(d=d, params=params, cond_dist=cond_dist)

  if(identification == "heteroskedasticity") { # Check W lambdas
    W <- pick_W(p=p, M=M, d=d, params=params, identification=identification)
    lambdas <- pick_lambdas(p=p, M=M, d=d, params=params, identification=identification)
    if(any(W[B_constraints == 0] != 0, na.rm=TRUE)) {
      stop("The matrix W doesn't satisfy the zero constraints in B_constraints")
    } else if(any(W[B_constraints > 0] <= 0, na.rm=TRUE) || any(W[B_constraints < 0] >= 0, na.rm=TRUE)) {
      stop("The matrix W doesn't satisfy the strict sign constraints in B_constraints")
    }
    if(any(lambdas <= 0)) {
      stop("The eigenvalues 'lambdas' should be strictly positive")
    }
  } else if(identification == "non-Gaussianity" || cond_dist == "ind_Student") { # Check B matrices
    for(i1 in 1:M) {
      if(any(all_Omegas[, , i1][B_constraints == 0] != 0, na.rm=TRUE)) {
        stop(paste0("The impact matrix of Regime ", i1, " doesn't satisfy the zero constraints in B_constraints"))
      } else if(any(all_Omegas[, , i1][B_constraints > 0] <= 0, na.rm=TRUE) || any(all_Omegas[B_constraints < 0] >= 0, na.rm=TRUE)) {
        stop(paste0("The impact matrix of Regime ", i1, " doesn't satisfy the strict sign constraints in B_constraints"))
      } else if(abs(det(all_Omegas[, , i1])) < posdef_tol) {
        stop(paste0("The impact matrix of Regime ", i1, " is not invertible"))
      }
    }
    if(!missing(transition_weights)) { # Check the invertibility of the time-varying impact matrix for all t
      for(i1 in 1:nrow(transition_weights)) { # Impact matrices in all_Omegas
        if(abs(det(apply(X=all_Omegas, MARGIN=c(1, 2), FUN=function(mat) sum(mat*transition_weights[i1,])))) < posdef_tol) {
          stop("The impact matrix B_t is not invertible for all t")
        }
      }
    }
  }

  if(cond_dist == "Student") { # Check degrees of freedom parameter
    if(distpars <= 2 + distpar_tol) {
      stop("The degrees of freedom parameter needs to be strictly larger than two (with large enough numerical tolerance)!")
    }
  } else if(cond_dist == "ind_Student") {
    if(any(distpars <= 2 + distpar_tol)) {
      stop("The degrees of freedom parameters need to be strictly larger than two (with large enough numerical tolerance)!")
    }
  }
  if(weight_function == "relative_dens") {
    if(M >= 2 & sum(weightpars[-M]) >= 1) {
      stop("The transition weight parameter alphas must sum to one!")
    } else if(any(weightpars <= 0)) {
      stop("The transition weight parameter alphas must be strictly larger than zero!")
    }
  } else if(weight_function == "logistic" || weight_function == "exponential") {
    if(weightpars[2] <= 0 + weightpar_tol) {
      stop(paste0("The scale parameter of ", weight_function,
                  " transition weights needs to be strictly positive (with large enough numerical tolerance)"))
    }
  } else if(weight_function == "mlogit") {
    if(!is.numeric(weightpars)) stop("Transition weight parameters need to be numeric") # All real numbers are ok
  } else if(weight_function == "threshold") {
    if(!all(order(weightpars, decreasing=FALSE) == seq_len(M - 1))) {
      stop("The threshold parameters need be in an increasing ordering with threshold weight function")
    }
  } else if(weight_function == "exogenous") {
    # No weightpars to test
  }
  if(!stab_conds_satisfied(p=p, M=M, d=d, all_boldA=all_boldA, tolerance=stab_tol)) {
    stop("At least one of the regimes does not satisfy the stability condition (with large enough numerical tolerance)!")
  }

  if(cond_dist != "ind_Student" && identification != "non-Gaussianity") {
    for(m in 1:M) {
      if(any(eigen(all_Omegas[, , m], symmetric=TRUE, only.values=TRUE)$values < posdef_tol)) {
        stop(paste0("The conditional covariance matrix of Regime ", m,
                    " is not positive definite (with large enough numerical tolerance)!"))
      }
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

check_pMd <- function(p, M, d, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                      identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity")) {
  weight_function <- match.arg(weight_function)
  identification <- match.arg(identification)
  if(!all_pos_ints(M) || length(M) != 1) {
    stop("The argument M must be a positive integer!")
  }
  if(identification == "heteroskedasticity") {
    if(M == 1) {
      stop("M=1 models cannot be identified by heteroskedasticty")
    }
  }
  if(!all_pos_ints(p) || length(p) != 1) {
    stop("The argument p must be a positive integer!")
  }
  if(!missing(d)) {
    if(d < 2 | d%%1 != 0) {
      stop("The argument d, the number of columns in the data matrix, has to be a positive integer larger than one!")
    }
  }
  if(weight_function == "logistic" && M != 2) {
    stop("Only two regime (M=2) models are accommodated by the logistic weight function. Use weight_function = 'mlogit' for M>2.")
  } else if(weight_function == "exponential" && M != 2) {
    stop("Only two regime (M=2) models are accommodated by the exponential weight function.")
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
    if(nrow(data) < p + 2) stop("The data must contain at least p+2 observations!")
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

n_params <- function(p, M, d, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                     weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student"),
                     identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity"),
                     AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  identification <- match.arg(identification)

  # Mean pars
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

  # Covmat pars
  if(identification == "non-Gaussianity" || cond_dist == "ind_Student") { # impact matrices parametrized directly
    if(is.null(B_constraints)) {
      n_zeros <- 0
    } else {
      n_zeros <- M*sum(B_constraints == 0, na.rm=TRUE) # Same number zeros for all the regimes
    }
    n_covmat_pars <- M*d^2 - n_zeros
  } else if(identification == "heteroskedasticity") {
    if(is.null(B_constraints)) {
      n_zeros <- 0
    } else {
      n_zeros <- sum(B_constraints == 0, na.rm=TRUE)
    }
    n_covmat_pars <- d^2 + d*(M - 1) - n_zeros
  } else { # identification == "reduced_form" or "recursive"
    n_covmat_pars <- M*d*(d + 1)/2 # No B_constraints available here
  }

  # Weight pars
  if(weight_function == "exogenous" || M == 1) {
    n_weight_pars <- 0
  } else {
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
  }

  # dist pars
  if(cond_dist == "Gaussian") {
    n_dist_pars <- 0
  } else if(cond_dist == "Student") {
    n_dist_pars <- 1 # degrees of freedom param
  } else { # cond_dist == "ind_Student"
    n_dist_pars <- d # Degrees of freedom for each component
  }

  n_mean_pars + n_ar_pars + n_covmat_pars + n_weight_pars + n_dist_pars
}


#' @title Check the constraint matrix has the correct form
#'
#' @description \code{check_constraints} checks that the constraints are correctly set.
#'
#' @inheritParams loglikelihood
#' @inheritParams stab_conds_satisfied
#' @return Does return anything but checks the constraints and throws an error if something is wrong.
#' @keywords internal

check_constraints <- function(data, p, M, d, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                              weightfun_pars=NULL, parametrization=c("intercept", "mean"),
                              identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity"),
                              AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL) {
  weight_function <- match.arg(weight_function)
  parametrization <- match.arg(parametrization)
  identification <- match.arg(identification)
  weightfun_pars <- check_weightfun_pars(data=data, p=p, M=M, d=d, weight_function=weight_function,
                                         weightfun_pars=weightfun_pars, cond_dist="Gaussian")
  # cond_dist="Gaussian" used above since we dont want to check rel_dens distribution here

  # Check AR_constraints
  if(!is.null(AR_constraints)) {
    if(!is.matrix(AR_constraints) | !is.numeric(AR_constraints)) {
      stop("The argument AR_constraints should be a numeric matrix (or NULL if no constraints should be employed).")
    } else if(nrow(AR_constraints) != M*p*d^2) {
      stop("The AR constraint matrix should have M*p*d^2 rows.")
    } else if(ncol(AR_constraints) > nrow(AR_constraints)) {
      stop("The AR constraint matrix has more columns than rows! What are you doing??")
    } else if(qr(AR_constraints)$rank != ncol(AR_constraints)) {
      stop("The AR constraint matrix should have full column rank.")
    }
  }

  # Check mean_constraints
  if(!is.null(mean_constraints)) {
    if(parametrization != "mean") {
      stop("mean_constraints are available only for models with parametrization = 'mean'.")
    }
    if(!is.list(mean_constraints)) {
      stop("The argument mean_constraints should a list (or null if mean parameters are not constrained).")
    } else if(length(mean_constraints) == 0) {
      stop("The argument mean_constraints should not of length zero.")
    }
    for(i1 in 1:length(mean_constraints)) {
      if(!is.numeric(mean_constraints[[i1]]) || length(mean_constraints[[i1]]) == 0) {
        stop("The elements of mean_constraints should be numeric vectors with strictly positive length.")
      }
    }
    tmp <- sort(unlist(mean_constraints), decreasing=FALSE)
    if(length(tmp) != M || !all(tmp == 1:M)) {
      stop("The argument mean_constraints should contains all regimes in some group exactly once.")
    }
  }

  # Check weight_constraints
  if(!is.null(weight_constraints)) {
    if(M == 1) {
      stop("weight_constraints cannot be employed for models with M=1 (because there are no weight parameters).")
    }
    if(!is.list(weight_constraints) || length(weight_constraints) != 2) {
      stop("The argument weight_constraints should be a list of length two.")
    }
    if(weight_function == "relative_dens" || weight_function == "threshold") {
      n_nonconstr_weightpars <- M - 1
    } else if(weight_function == "logistic" || weight_function == "exponential") {
      n_nonconstr_weightpars <- 2
    } else if(weight_function == "mlogit") {
      n_nonconstr_weightpars <- (M - 1)*(1 + length(weightfun_pars[[1]])*weightfun_pars[[2]])
    }
    if(weight_function != "exogenous") {
      if(!all(weight_constraints[[1]] == 0)) { # R != 0
        if(!is.matrix(weight_constraints[[1]]) || !is.numeric(weight_constraints[[1]])) {
          stop("The first element of the argument weight_constraints should be a numeric matrix R (or 0 or NULL).")
        } else if(nrow(weight_constraints[[1]]) != n_nonconstr_weightpars) {
          stop("The first element of weight_constraints (matrix R) has incorrect number of rows.")
        } else if(ncol(weight_constraints[[1]]) > nrow(weight_constraints[[1]])) {
          stop("The first element of weight_constraints (matrix R) has more columns than rows! What are you doing??")
        } else if(qr(weight_constraints[[1]])$rank != ncol(weight_constraints[[1]])) {
          stop("The first element of weight_constraints (matrix R) should have full column rank (or it should equal to zero).")
        }
      }
      # Check r
      if(!is.numeric(weight_constraints[[2]]) || !is.vector(weight_constraints[[2]])) {
        stop("The second element of the argument weight_constraints should be a numeric vector r.")
      } else if(length(weight_constraints[[2]]) != n_nonconstr_weightpars) {
        stop("The second element of the argument weight_constraints (vector r) has wrong dimension.")
      }
    }

    # Warnings if no errors
    if(weight_function == "logistic" || weight_function == "exponential") {
      if(weight_constraints[[2]][2] < 0) {
        warning(paste0("When weight_function=", weight_function,
                       "the scale parameter needs to be strictly positive, and there is a negative",
                       "constraint in r for the scale parameter, implying that the estimation may fail."))
      }
    }
    if(weight_function == "relative_dens") {
      if(any(weight_constraints[[2]] < 0)) {
        warning(paste("When weight_function='relative dens', the weight parameters need to be strictly positive,",
                      "and there is a negative constraint in r for a weight parameter, implying that the estimation may fail."))
      }
    }
    if(weight_function == "threshold") {
      if(!all(order(weight_constraints[[2]]) == seq_len(M - 1))) {
        warning(paste("When weight_function='threshold', the thresholds need to be in an increasing ordering,",
                      "and the constraints imposed in r for the thresholds are not in an increasing ordering,",
                      "implying that the estimation may fail."))
      }
    }
    if(weight_function == "exogenous" && !is.null(weight_constraints)) {
      warning("weight_constraints cannot be used with the exogenous weight function.")
    }
  }

  # Check B_constraints
  if(!is.null(B_constraints)) {
    if(identification != "heteroskedasticity" && identification != "non-Gaussianity") {
      stop("B_constraints are currently available only for models with identication = 'heteroskedasticity' or 'non-Gaussianity'.")
    }
    if(!is.matrix(B_constraints) || any(dim(B_constraints) != d)) {
      stop("B_constraints should be a (d x d) matrix")
    }
    n_zeros1 <- vapply(1:d, function(i1) sum(B_constraints[i1,] == 0, na.rm=TRUE), numeric(1))
    n_zeros2 <- vapply(1:d, function(i1) sum(B_constraints[,i1] == 0, na.rm=TRUE), numeric(1))
    if(any(n_zeros1 == d) || any(n_zeros2 == d)) {
      stop("The impact matrix/matrices should be invertible, so you cannot constrain it/them to be singular via B_constraints.")
    }
  }
}


#' @title Check the argument \code{weightfun_pars}
#'
#' @description \code{check_weightfun_pars} checks that the argument \code{weightfun_pars}.
#'   is correctly set, if not, tries to correct them.
#'
#' @inheritParams loglikelihood
#' @return Does checks the argument \code{weightfun_pars} and throws an error if something is wrong; returns
#'   a corrected version of the argument if possible.
#' @keywords internal

check_weightfun_pars <- function(data, p, M, d, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                                 weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student")) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)

  if(weight_function == "relative_dens") {   # weightfun_pars are not used in weight_function == "relative_dens"
    weightfun_pars <- NULL
    if(cond_dist != "Gaussian") {
      stop("When weight_function == 'relative_dens', only Gaussian conditional distribution is supported (cond_dist='Gaussian').")
    }
  } else if(weight_function %in% c("logistic", "exponential", "threshold")) {
    if(!is.numeric(weightfun_pars) || !is.vector(weightfun_pars) || length(weightfun_pars) != 2) {
      stop(paste0("When weight_function == ", weight_function,
                  " the argument weightfun_pars should be be a length two numeric vector."))
    }
    if(!weightfun_pars[1] %in% 1:d) {
      stop(paste0("When weight_function == ", weight_function,
                  " the first element of argument weightfun_pars, i.e., the switching variable,",
                  " should be an integer in 1,...,ncol(data)."))
    } else if(!weightfun_pars[2] %in% 1:p) {
      stop(paste0("When weight_function == ", weight_function,
           " the second element of argument weightfun_pars, i.e., the lag of the switching variable,
           should be an integer in 1,...,p."))
    }
  } else if(weight_function == "mlogit") {
    if(!is.list(weightfun_pars) || length(weightfun_pars) != 2) {
      stop("When weight_function == 'mlogit', the argument weightfun_pars should be be a list of length two..")
    }
    if(!is.numeric(weightfun_pars[[1]]) || !all(weightfun_pars[[1]] %in% 1:d) ||
       length(unique(weightfun_pars[[1]])) != length(weightfun_pars[[1]]) ) {
      stop("When weight_function == 'mlogit', the first element of the argument weightfun_pars should be a numeric
           vector containing the switching variables, i.e., it should have unique elements in 1,...,d where d is
           the number of variables in the data.")
    }
    if(!is.numeric(weightfun_pars[[2]]) || !all(weightfun_pars[[2]] %in% 0:p) || length(weightfun_pars[[2]]) != 1) {
      stop("When weight_function == 'mlogit', the second element of the argument weightfun_pars should be an integer
            in 0,1,...,p determining the number of lags to be used in the weight function.")
    }
    weightfun_pars[[1]] <- sort(weightfun_pars[[1]], decreasing=FALSE)

    if(is.null(names(weightfun_pars))) {
      names(weightfun_pars) <- c("vars", "lags")
    } else if(!all(names(weightfun_pars) == c("vars", "lags"))) {
      names(weightfun_pars) <- c("vars", "lags")
    }
  } else if(weight_function == "exogenous") {
    if(missing(data) || is.null(data)) { # Does not check data
      if(!is.numeric(weightfun_pars) || !is.matrix(weightfun_pars) || ncol(weightfun_pars) != M) {
        stop("When weight_function == 'exogenous' and there is no data, the argument weightfun_pars should be a numeric matrix with M columns.")
      }
    } else {
      if(!is.numeric(weightfun_pars) || !is.matrix(weightfun_pars) || ncol(weightfun_pars) != M || nrow(weightfun_pars) != nrow(data) - p) {
        stop("When weight_function == 'exogenous', the argument weightfun_pars should be a numeric (nrow(data)-p x M) matrix.")
      }
    }
    if(any(weightfun_pars < 0)) {
       stop("When weight_function == 'exogenous', the argument weightfun_pars should not contain strictly negative elements.")
     } else if(!isTRUE(all.equal(rowSums(weightfun_pars), rep(1, times=nrow(weightfun_pars))))) {
       warning("The each row of the matrix weightfun_pars should sum to one, normalizing them to sum to one.")
       weightfun_pars <- weightfun_pars/rowSums(weightfun_pars)
     }
  }
  weightfun_pars
}


#' @title Checks whether the given object has class attribute 'stvar'
#'
#' @description \code{check_stvar} checks that the object has class attribute 'stvar'.
#'
#' @param object S3 object to be tested
#' @param object_name what is the name of the object that should of class 'stvar'?
#' @return Throws an error if the object doesn't have the class attribute 'stvar'.
#' @keywords internal

check_stvar <- function(object, object_name) {
  if(missing(object_name)) object_name <- "stvar"
  if(!any(class(object) == "stvar")) {
    stop(paste("The object", object_name,
               "has to be of class 'stvar',",
               "typically created with the function 'STVAR', 'fitSTVAR', or 'fitSSTVAR'"))
  }
}


#' @title Checks whether the given exogenous transition weights for simulation are correctly specified.
#'
#' @description \code{check_exoweights} checks whether the given exogenous transition weights for
#'  simulation are correctly specified.
#'
#' @param exo_weights Exogenous transition weights weights given for simulation in some context.
#' @param how_many_rows how many rows the exogenous weights should have?
#' @param name_of_row_number what is the name of the object whose value should determine the
#'   the number of rows in the exogenous weights?
#' @details Used by simulate.stvar, predict.stvar, GIRF, and GFEVD.
#' @return Throws an error if the exogenous weights are incorrectly specified.
#' @keywords internal

check_exoweights <- function(M, exo_weights, how_many_rows, name_of_row_number) {
  # Check the exogenous weights given for simulation
  if(is.null(exo_weights)) stop("Exogenous weights must be provided in the argument 'exo_weights' when weight_function is 'exogenous'.")
  if(!is.matrix(exo_weights)) stop("Exogenous weights 'exo_weights' must be a matrix.")
  if(nrow(exo_weights) != how_many_rows) stop(paste0("Exogenous weights 'exo_weights' must have '", name_of_row_number, "' rows."))
  if(ncol(exo_weights) != M) stop("Exogenous weights 'exo_weights' must have M columns.")
  if(any(exo_weights < 0) || any(exo_weights > 1)) stop("Exogenous weights 'exo_weights' must be in [0, 1].")
  # Check that exogenous weights sum to one at each row withing a numerical accuracy:
  if(!all(abs(rowSums(exo_weights) - 1) < 1e-10)) stop("Exogenous weights 'exo_weights' must sum to one at each row.")
}
