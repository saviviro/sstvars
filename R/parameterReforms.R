#' @title Reform data
#'
#' @description \code{reform_data} reforms the data into a form that is
#'   easier to use when calculating log-likelihood values etc.
#'
#' @inheritParams loglikelihood
#' @details Assumes the observed data is \eqn{y_{-p+1},...,y_0,y_1,...,y_{T}}.
#' @return Returns the data reformed into a \eqn{((n_{obs}-p+1)\times dp)} matrix. The i:th row
#'   of the matrix contains the vector \eqn{(y_{i-1},...,y_{i-p})} \eqn{(dp\times 1)}, where
#'   \eqn{y_{i}=(y_{1i},...,y_{di})} \eqn{(d \times 1)}.
#' @section Warning:
#'  No argument checks!
#' @keywords internal

reform_data <- function(data, p) {
  d <- ncol(data)
  n_obs <- nrow(data)
  T_obs <- n_obs - p
  matrix(vapply(1:p, function(i1) as.vector(data[(p - i1 + 1):(T_obs + p - i1 + 1),]),
                numeric((n_obs - p + 1)*d)), nrow=n_obs - p + 1, byrow=FALSE)
}


#' @title Form the \eqn{(dp\times dp)} "bold A" matrices related to the VAR processes
#'
#' @description \code{form_boldA} creates the "bold A" (companien form) coefficient matrices related to
#'   VAR processes.
#'
#' @inheritParams pick_allA
#' @param all_A 4D array containing all coefficient matrices \eqn{A_{m,i}}, obtained from \code{pick_allA}.
#' @details The "bold A" (companion form) matrix is given, for instance, in Lütkepohl (2005, p. 15).
#' @return Returns 3D array containing the \eqn{(dp \times dp)} "bold A" matrices related to each component VAR-process.
#'  The matrix \strong{\eqn{A_{m}}} can be obtained by choosing \code{[, , m]}.
#' @section Warning:
#'  No argument checks!
#' @references
#'  \itemize{
#'    \item Lütkepohl H. 2005. New Introduction to Multiple Time Series Analysis, \emph{Springer}.
#'  }
#' @keywords internal

form_boldA <- function(p, M, d, all_A) {
  I_all <- diag(nrow=d*(p - 1))
  ZER_all <- matrix(0, nrow=d*(p - 1), ncol=d)
  array(vapply(1:M, function(m) rbind(matrix(all_A[, , 1:p, m], nrow=d, byrow=FALSE),
                                      cbind(I_all, ZER_all)), numeric((d*p)^2)), dim=c(d*p, d*p, M))
}


#' @title Change parametrization of a parameter vector
#'
#' @description \code{change_parametrization} changes the parametrization of the given parameter
#'   vector to \code{change_to}.
#'
#' @inheritParams loglikelihood
#' @inheritParams form_boldA
#' @param change_to
#'   \describe{
#'     \item{If you want to switch between mean and intercept parametrizations:}{either "intercept" or "mean" specifying
#'       to which parametrization it should be switched to. If set to \code{"intercept"}, it's assumed that \code{params}
#'       is mean parametrized, and if set to \code{"mean"} it's assumed that \code{params} is intercept parametrized.}
#'     \item{If you want to switch between the paramterizations \eqn{B_{y,t}=\sum_{m=1}^M\alpha_{m,t}B_m} and
#'           \eqn{B_{y,t}=B_1 + \sum_{m=2}^M\alpha_{m,t}B_m^{*}}, \eqn{B_m^{*} = B_m - B_1}:}{either "orig"
#'            (with \eqn{B_m}) or "alt" (with \eqn{B_m^{*}}). It is assumed that the parameter vector is in the
#'            other parametrization than the one specified in \code{change_to}.}
#'   }
#' @details Parametrization cannot be changed for models with mean constraints! Note that changing between "orig" and "alt"
#'   changes the meaning of sign constraints in \code{B_constraints} (sign constraints imposed on "alt" is very different to
#'   those imposed on "orig"). Thus, this function should not be used to switch between "orig" and "alt" when sign constraints
#'   are imposed!
#' @return Returns parameter vector described in \code{params}, but with parametrization changed according
#'   to \code{change_to}.
#' @section Warning:
#'  No argument checks!
#' @keywords internal

change_parametrization <- function(p, M, d, params,
                                   weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                                   weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t"),
                                   identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity"),
                                   AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL,
                                   change_to=c("intercept", "mean", "orig", "alt")) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  identification <- match.arg(identification)
  change_to <- match.arg(change_to)
  re_params <- params
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, weight_function=weight_function, cond_dist=cond_dist,
                                    identification=identification, AR_constraints=AR_constraints,
                                    mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                    B_constraints=B_constraints, weightfun_pars=weightfun_pars)
  Id <- diag(nrow=d)
  all_A <- pick_allA(p=p, M=M, d=d, params=params)
  all_phi0_or_mu <- pick_phi0(M=M, d=d, params=params)
  # Nothing to change in distpars or weightpars, already in place

  if(change_to == "mean" || change_to == "intercept") {
    stopifnot(is.null(mean_constraints))
    # Calculate means/intercepts and insert them to re_params
    if(change_to == "mean") { # params has original parametrization with intercept
      re_params[1:(M*d)] <- vapply(1:M, function(m) solve(Id - rowSums(all_A[, , , m, drop=FALSE], dims=2), all_phi0_or_mu[,m]),
                                   numeric(d))
    } else if(change_to == "intercept") { # params has mean parameters instead of phi0
      re_params[1:(M*d)] <- vapply(1:M, function(m) (Id - rowSums(all_A[, , , m, drop=FALSE], dims=2))%*%all_phi0_or_mu[,m],
                                   numeric(d))
    }
  } else { # Change between nong_orig and nong_alt
    if(M == 1) { # B_1 = B_1*, so nothing to change
      return(re_params)
    } else if(!cond_dist %in% c("ind_Student", "ind_skewed_t")) { # Other cond dists do not apply to B matrices below
      return(re_params)
    }
    # All B_m or B_m^* matrices:
    all_Omegas <- pick_Omegas(p=p, M=M, d=d, params=params, cond_dist=cond_dist,
                              identification=identification)

    if(change_to == "orig") {
      # Matrices in all_Omegas are B_m* (B_1=B_1*). Obtain B_m from them:
      for(m in 2:M) {
        all_Omegas[, , m] <- all_Omegas[, , m] + all_Omegas[, , 1]
      }
    } else { # change_to == "alt"
      # Matrices in all_Omegas are B_m (B_1=B_1*). Obtain B_m* from them:
      for(m in 2:M) {
        all_Omegas[, , m] <- all_Omegas[, , m] - all_Omegas[, , 1]
      }
    }
    # Removed elements from the B matrices that nonparametrized zeros constrainted
    # to zero by B_constraints, and collect the remaining elements to a vector:
    if(is.null(B_constraints)) {
      new_B_pars <- c(all_Omegas)
    } else {
      # Removed elements from the B matrices that nonparametrized zeros constrainted
      # to zero by B_constraints, and collect the remaining elements to a vector.
      n_zeros <- sum(B_constraints == 0, na.rm=TRUE) # n zeros in each B matrix
      new_B_pars <- numeric(M*(d^2 - n_zeros))
      for(m in 1:M) {
        new_B_pars[((m - 1)*(d^2 - n_zeros) + 1):(m*(d^2 - n_zeros))] <- all_Omegas[, , m][B_constraints != 0 | is.na(B_constraints)]
      }
    }

    # Calculate the number of mean and AR parameters:
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

    # Insert the new B_pars to re_params:
    re_params[(n_mean_pars + n_ar_pars + 1):(n_mean_pars + n_ar_pars + length(new_B_pars))] <- new_B_pars
  }

  re_params
}


#' @title Sort regimes in parameter vector according to transition weights into a decreasing order
#'
#' @description \code{sort_regimes} sorts regimes in the parameter vector according to
#'   the transition weight parameters.
#'
#' @inheritParams loglikelihood
#' @details Constrained parameter vectors are not supported (except \code{B_constraints} for structural models identified
#'   by heteroskedasticity).
#' @return Returns sorted parameter vector of the form described for the argument \code{params},
#'   with the regimes sorted so that...
#'   \describe{
#'     \item{If \code{weight_function == "relative_dens"}:}{\eqn{\alpha_{1}>...>\alpha_{M}}.}
#'     \item{If \code{weight_function == "logistic"}:}{Nothing to sort, so returns the original parameter vector given in \code{param}.}
#'     \item{If \code{weight_function == "mlogit"}:}{Does not currently sort, so returns the original parameter vector given in
#'       \code{param}.}
#'     \item{If \code{weight_function == "exponential"}:}{Nothing to sort, so returns the original parameter vector given in
#'       \code{param}.}
#'     \item{If \code{weight_function == "threshold"}:}{The increasing ordering of the thresholds is imposed in the parameter space,
#'           so nothing to sort and thereby returns the original parameter vector given in \code{param}.}
#'     \item{If \code{weight_function == "exogenous"}:}{Does not sort but returns the original parameter vector.}
#'   }
#' @keywords internal

sort_regimes <- function(p, M, d, params, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                         weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t"),
                         identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity"), B_constraints=NULL) {
  weight_function <- match.arg(weight_function)
  if(M == 1 || weight_function %in% c("logistic", "mlogit", "exponential", "threshold", "exogenous")) {
    return(params) # Does not sort / nothing to sort
  }
  cond_dist <- match.arg(cond_dist)
  identification <- match.arg(identification)

  all_weightpars <- pick_weightpars(p=p, M=M, d=d, params=params, weight_function=weight_function, cond_dist=cond_dist)
  if(weight_function == "relative_dens") {
    new_order <- order(all_weightpars, decreasing=TRUE)
    if(all(new_order == 1:M)) {
      return(params)
    }
    new_weightpars <- all_weightpars[new_order][-M]
  } else {
    return(params) # Other weight functions do not have sorting implemented
  }
  if(cond_dist %in% c("Student", "ind_Student", "ind_skewed_t") || identification == "non-Gaussianity") {
    # Should never end up here, but exists to detect errors in a flawed update
    stop(paste("cond_dist or identification in sort regime implies a weight function that",
               "does not support sorting ('alt parametrization' assumed)!"))
  }

  all_phi0 <- pick_phi0(M=M, d=d, params=params)
  all_A <- matrix(pick_allA(p=p, M=M, d=d, params=params), ncol=M) #matrix(params[(d*M + 1):(d*M + d^2*p*M)], ncol=M)

  # Covmat pars
  if(identification == "heteroskedasticity") {
    if(is.null(B_constraints)) {
      n_zeros <- 0
    } else {
      n_zeros <- sum(B_constraints == 0, na.rm=TRUE)
    }
    W_pars <- params[(d*M + M*p*d^2 + 1):(d*M + M*p*d^2 + d^2 - n_zeros)] # Does not include non parametrized zeros
    if(is.null(B_constraints)) {
      old_W <- W_pars
    } else {
      old_W <- numeric(d^2) # We include the unparametrized zeros here
      old_W[B_constraints != 0 | is.na(B_constraints)] <- W_pars
    }
    lambdas <- params[(d*M + M*p*d^2 + d^2 - n_zeros + 1):(d*M + M*p*d^2 + d^2 - n_zeros + d*(M - 1))]
    new_covmatpars <- Wvec(redecompose_Omegas(M=M, d=d, W=old_W, lambdas=lambdas, perm=new_order)) # Wvec removes the unparametrized zeros
  } else { # Identification %in% c("reduced_form", "recursive")
    all_Omega <- matrix(params[(d*M*(1 + p*d) + 1):(d*M*(1 + p*d) + M*d*(d + 1)/2)], nrow=d*(d + 1)/2, ncol=M)
    new_covmatpars <- all_Omega[,new_order]
  }

  # Weight and dist pars
  all_weightpars <- pick_weightpars(p=p, M=M, d=d, params=params, weight_function=weight_function,
                                    weightfun_pars=weightfun_pars, cond_dist=cond_dist)
  all_distpars <- pick_distpars(d=d, params=params, cond_dist=cond_dist)

  c(all_phi0[,new_order], all_A[,new_order], new_covmatpars, new_weightpars, all_distpars)
}


#' @title Change the parameters of a specific regime of the given parameter vector
#'
#' @description \code{change_regime} changes the regime parameters
#'   \eqn{(\phi_{m,0},vec(A_{m,1}),...,\vec(A_{m,p}),vech(\Omega_m))} (replace \eqn{vech(\Omega_m)}
#'   by \eqn{vec(B_m)} for \code{cond_dist="ind_Student"})
#'   of the given regime to the new given parameters.
#'
#' @inheritParams pick_regime
#' @param regime_pars
#'   \describe{
#'     \item{If \code{cond_dist="Gaussian"} or \code{cond_dist="Student"}:}{rhe \eqn{(dp + pd^2 + d(d+1)/2)} vector
#'      \eqn{(\phi_{m,0},vec(A_{m,1}),...,\vec(A_{m,p}),vech(\Omega_m))}.}
#'     \item{If \code{cond_dist="ind_Student"} or \code{"ind_skewed_t"}:}{the \eqn{(dp + pd^2 + d^2 + 1)} vector
#'      \eqn{(\phi_{m,0},vec(A_{m,1}),...,\vec(A_{m,p}),vec(B_m))}.}
#'   }
#' @return Returns parameter vector with \code{m}:th regime changed to \code{regime_pars}.
#' @details Does not support constrained models or structural models. Weight parameters and distribution
#'   parameters are not changed.
#' @inherit in_paramspace references
#' @keywords internal

change_regime <- function(p, M, d, params, m, regime_pars, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t")) {
  cond_dist <- match.arg(cond_dist)
  new_pars <- params
  new_pars[((m - 1)*d + 1):(m*d)] <- regime_pars[1:d]
  new_pars[(M*d + (m - 1)*p*d^2 + 1):(M*d + m*p*d^2)] <- regime_pars[(d + 1):(d + p*d^2)]
  if(cond_dist %in% c("ind_Student", "ind_skewed_t")) {
    new_pars[(M*d + M*p*d^2 + (m - 1)*d^2 + 1):(M*d + M*p*d^2 + m*d^2)] <- regime_pars[(d + p*d^2 + 1):(d + p*d^2 + d^2)]
  } else {
    new_pars[(M*d + M*p*d^2 + (m - 1)*d*(d + 1)/2 + 1):(M*d + M*p*d^2 + m*d*(d + 1)/2)] <- regime_pars[(d + p*d^2 + 1):length(regime_pars)]
  }
  new_pars
}


#' @title Reform constrained parameter vector into the "standard" form
#'
#' @description \code{reform_constrained_pars} reforms constrained parameter vector
#'   into the form that corresponds to unconstrained parameter vectors.
#'
#' @inheritParams loglikelihood
#' @param change_na change NA parameter values of constrained models to -9.999?
#' @return Returns "regular model" parameter vector corresponding to the constraints.
#' @section Warning:
#'  No argument checks!
#' @inherit in_paramspace references
#' @keywords internal

reform_constrained_pars <- function(p, M, d, params,
                                    weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                                    weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t"),
                                    identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity"),
                                    AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL,
                                    other_constraints=NULL, change_na=FALSE) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  identification <- match.arg(identification)

  if(is.null(AR_constraints) && is.null(mean_constraints) && is.null(weight_constraints) &&
     is.null(B_constraints) && is.null(other_constraints)) {
    return(params)
  } else if(is.null(AR_constraints) && is.null(mean_constraints) && is.null(weight_constraints) &&
            is.null(B_constraints) && identification == "non-Gaussianity") {
    return(params) # Other constraints do not affect the parameter vector here
  }

  ## Obtain the mean parameters ##
  if(is.null(mean_constraints)) {
    less_pars <- 0 # Number of parameters less compared to models without same mean constraints
  } else {
    g <- length(mean_constraints) # Number groups with the same mean parameters
    less_pars <- d*(M - g) # Number of parameters less compared to models without same mean constraints
  }
  if(is.null(mean_constraints)) {
    all_phi0 <- matrix(params[1:(d*M)], nrow=d, ncol=M) # params[((m - 1)*d + 1):(m*d)]
  } else { # mean_constraints employed
    group_phi0 <- matrix(params[1:(d*g)], nrow=d, ncol=g) # Column for each group
    all_phi0 <- matrix(NA, nrow=d, ncol=M) # Storage for all phi0 (=mean parameters in this case)
    for(i1 in 1:g) {
      all_phi0[,mean_constraints[[i1]]] <- group_phi0[,i1]
    }
  }

  ## Obtain the AR matrix parameters ##
  if(is.null(AR_constraints)) { # For structural model with constrained structural parameters but no AR constraints
    q <- M*p*d^2
    psi_expanded <- params[(d*M + 1 - less_pars):(d*M + d^2*p*M - less_pars)] # AR coefficients (without constraints)
    psiNA <- FALSE
  } else {
    q <- ncol(AR_constraints)
    psi <- params[(M*d + 1 - less_pars):(M*d + q - less_pars)]
    if(change_na) {
      if(anyNA(psi)) {
        warning("Replaced some NA values with -9.999")
        psiNA <- TRUE
      } else {
        psiNA <- FALSE
      }
      psi[is.na(psi)] <- -9.999
    }
    psi_expanded <- AR_constraints%*%psi
  }

  ## Obtain the covariance matrix parameters ##
  if(identification == "non-Gaussianity" || cond_dist == "ind_Student" || cond_dist == "ind_skewed_t") {
    if(is.null(B_constraints)) {
      covmatpars <- params[(d*M - less_pars + q + 1):(d*M - less_pars + q + M*d^2)]
      less_covmatpars <- 0
    } else {
      n_zeros <- sum(B_constraints == 0, na.rm=TRUE)
      less_covmatpars <- M*n_zeros # The same zero constraints are imposed on each impact matrix
      B_constraints[is.na(B_constraints)] <- 1 # Insert some nonzero value to replace the NA values
      all_B_m <- array(0, dim=c(d, d, M))
      for(m in 1:M) { # Fill the non-zero parameters to the correct places
        all_B_m[, , m][B_constraints != 0] <- params[(d*M - less_pars + q + (m - 1)*(d^2 - n_zeros)
                                                      + 1):(d*M - less_pars + q + m*(d^2 - n_zeros))]
      }
      covmatpars <- as.vector(all_B_m)
    }
  } else if(identification == "heteroskedasticity") {
    if(is.null(B_constraints)) {
      if(is.null(other_constraints$fixed_lambdas)) {
        covmatpars <- params[(d*M - less_pars + q + 1):(d*M - less_pars + q + d^2 + d*(M - 1))]
        less_covmatpars <- 0
      } else { # lambdas defined by fixed_lambdas
        covmatpars <- c(params[(d*M - less_pars + q + 1):(d*M - less_pars + q + d^2)], other_constraints$fixed_lambdas)
        less_covmatpars <- d*(M - 1)
      }
    } else { # B_constraints imposed on W
      n_zeros <- sum(B_constraints == 0, na.rm=TRUE)
      less_covmatpars <- n_zeros
      new_W <- numeric(d^2)
      W_pars <- params[(d*M - less_pars + q + 1):(d*M - less_pars + q + d^2 - n_zeros)] # Does not include non parametrized zeros
      new_W[B_constraints != 0 | is.na(B_constraints)] <- W_pars
      if(M > 1) {
        if(is.null(other_constraints$fixed_lambdas)) {
          lambdas <- params[(d*M - less_pars + q + d^2 - n_zeros + 1):(d*M - less_pars + q + d^2 - n_zeros + d*(M - 1))]
        } else { # lambdas defined by fixed_lambdas
          lambdas <- other_constraints$fixed_lambdas
          less_covmatpars <- n_zeros + d*(M - 1)
        }
      } else {
        lambdas <- numeric(0)
      }
      covmatpars <- c(new_W, lambdas)
    }
  } else { # identification == "reduced_form" || identification == "recursive" and cond_dist is Gaussian or Student
    covmatpars <- params[(d*M - less_pars + q + 1):(d*M - less_pars + q + M*d*(d + 1)/2)]
    less_covmatpars <- 0 # How much covmats pars less are there in the constrained param vector relative to the expanded one
  }
  n_covmatpars <- length(covmatpars) - less_covmatpars

  ## Obtain the weight parameters ##
  if(M > 1 && weight_function != "exogenous") {
    if(weight_function == "relative_dens" || weight_function == "threshold") {
      n_nonconstr_weightpars <- M - 1
    } else if(weight_function == "logistic" || weight_function == "exponential") {
      n_nonconstr_weightpars <- 2
    } else if(weight_function == "mlogit") {
      n_nonconstr_weightpars <- (M - 1)*(1 + length(weightfun_pars[[1]])*weightfun_pars[[2]])
    }
    if(is.null(weight_constraints)) {
      weightpars <- params[(d*M - less_pars + q + n_covmatpars + 1):(d*M - less_pars + q + n_covmatpars + n_nonconstr_weightpars)]
    } else { # Obtain unconstrained weightpars
      if(all(weight_constraints[[1]] == 0)) {
        weightpars <- weight_constraints[[2]] # alpha = r if R=0
      } else {
        l <- ncol(weight_constraints[[1]])
        xi <- params[(d*M - less_pars + q + n_covmatpars + 1):(d*M - less_pars + q + n_covmatpars + l)]
        weightpars <- weight_constraints[[1]]%*%xi + weight_constraints[[2]] # alpha = R%*%xi + r
      }
    }
  } else { # No weightpars if M == 1 or weight_function == "exogenous"
    weightpars <- numeric(0)
  }

  ## Obtain the distribution parameters ##
  if(cond_dist == "Gaussian") {
    distpars <- numeric(0)
  } else if(cond_dist == "Student") {
    distpars <- params[length(params)]
  } else if(cond_dist == "ind_Student") {
    distpars <- params[(length(params) - d + 1):length(params)]
  } else if(cond_dist == "ind_skewed_t") {
    distpars <- params[(length(params) - 2*d + 1):length(params)]
  }

  c(all_phi0, psi_expanded, covmatpars, weightpars, distpars)
}



#' @title Sort and sign change the columns of the impact matrices of the regimes so that the first element in each column of \eqn{B_1}
#'   is positive and in a decreasing order.
#'
#' @description \code{sort_impactmats} sorts and sign changes the columns of the impact matrices of the regimes so that the first element
#'   in each column of \eqn{B_1} is positive and in a decreasing order (for skewed distributions, the signs are not normalized).
#'   The same reordering and sign changes performed to the columns of \eqn{B_1} are applied to the rest of the impact matrices to
#'   obtain an observationally equivalent model.
#'
#' @inheritParams loglikelihood
#' @details This function is internally used by \code{GAfit} and \code{fitSTVAR}, so structural models or \code{B_constraints}
#'   are not supported.
#' @return Returns sorted parameter vector of the form described for the argument \code{params},
#'   with the regimes sorted so that...
#'   \describe{
#'     \item{If \code{cond_dist == "ind_Student"}:}{The parameter vector with the columns of the impact
#'           matrices sorted and sign changed so that the first element in each column of \eqn{B_1} is
#'           positive and in a decreasing order. Sorts also the degrees of freedom and skewness parameters
#'           (if any) accordingly.}
#'     \item{If \code{cond_dist == "ind_skewed_t"}:}{The parameter vector with the columns of the impact
#'           matrices sorted so that the first element in each column of \eqn{B_1} are in a decreasing order.}
#'     \item{Otherwise:}{Nothing to sort, so returns the original parameter vector given in \code{param}.}
#'   }
#' @keywords internal

sort_impactmats <- function(p, M, d, params,
                            weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                            weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t"),
                            AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  if(!cond_dist %in% c("ind_Student", "ind_skewed_t")) {
    return(params) # No impact matrices whose columns to sort
  }

  # The number of weight pars
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
  # The number of degrees of freedom parameters is d for cond_dist == "ind_Student",
  # whereas the number of df and skewness parameters is 2d for cond_dist == "ind_skewed_t".
  n_distpars <- ifelse(cond_dist == "ind_Student", d, 2*d) # Models with other cond dists do not arrive here

  # Obtain the impact matrices
  all_B_m <- array(params[(length(params) - n_weight_pars - n_distpars - M*d^2 + 1):(length(params) - n_weight_pars - n_distpars)],
                   dim=c(d, d, M))

  # Determine which columns should go through sign change, and change the signs of those columns
  if(cond_dist == "ind_Student") {
    for(i1 in 1:d) {
      if(all_B_m[1, i1, 1] < 0) {
        all_B_m[, i1, 1] <- -all_B_m[, i1, 1]
        # Change the signs of the corresponding columns in the rest of the impact matrices
        if(M > 1) {
          for(m in 2:M) {
            all_B_m[, i1, m] <- -all_B_m[, i1, m]
          }
        }
      }
    }
  }

  # Sort the columns of the impact matrices so that the first element in each column of B_1 is in a decreasing order,
  # i.e., the first row is in a decreasing order (note: the first element is positive due to the above sign changes).
  new_ordering <- order(all_B_m[1, , 1], decreasing=TRUE)
  for(m in 1:M) {
    all_B_m[, , m] <- all_B_m[, new_ordering, m]
  }

  # Fill in the new impact matrices to the parameter vector
  params[(length(params) - n_weight_pars - n_distpars - M*d^2 + 1):(length(params) - n_weight_pars - n_distpars)] <- as.vector(all_B_m)

  # Sort the degrees of freedom parameters accordingly:
  if(cond_dist == "ind_Student") { # Only degrees of freedom parameters to sort
    distpars <- params[(length(params) - n_distpars + 1):length(params)]
    params[(length(params) - n_distpars + 1):length(params)] <- distpars[new_ordering]
  } else { # ind_skewed_t, sort both degrees of freedom and skewness parameters
    distpars <- matrix(params[(length(params) - n_distpars + 1):length(params)], ncol=2) # [d,]
    distpars <- distpars[new_ordering,]
    params[(length(params) - n_distpars + 1):length(params)] <- as.vector(distpars)
  }

  # Return the new parameter vector
  params
}
