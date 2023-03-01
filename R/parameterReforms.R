#' @title Reform data
#'
#' @description \code{reform_data} reforms the data into a form that is
#'   easier to use when calculating log-likelihood values etc.
#'
#' @inheritParams loglikelihood
#' @return Returns the data reformed into a \eqn{((n_{obs}-p+1)x(dp))} matrix. The i:th row
#'   of the matrix contains the vector \eqn{(y_{i-1}',...,y_{i-p}')} \eqn{((dp)x1)}, where
#'   \eqn{y_{i}=(y_{1i},...,y_{di})} \eqn{(dx1)}.
#' @section Warning:
#'  No argument checks!
#' @keywords internal

reform_data <- function(data, p) {
  d <- ncol(data)
  n_obs <- nrow(data)
  T_obs <- n_obs - p
  matrix(vapply(1:p, function(i1) as.vector(data[(p - i1 + 1):(T_obs + p - i1 + 1),]), numeric((n_obs - p + 1)*d)), nrow=n_obs - p + 1, byrow=FALSE)
}


#' @title Form the \eqn{((dp)x(dp))} "bold A" matrices related to the VAR processes
#'
#' @description \code{form_boldA} creates the "bold A" (companien form) coefficient matrices related to
#'   VAR processes.
#'
#' @inheritParams pick_allA
#' @param all_A 4D array containing all coefficient matrices \eqn{A_{m,i}}, obtained from \code{pick_allA}.
#' @details The "bold A" (companion form) matrix is given, for instance, in Lütkepohl (2005, p. 15).
#' @return Returns 3D array containing the \eqn{((dp)x(dp))} "bold A" matrices related to each component VAR-process.
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
  array(vapply(1:M, function(m) rbind(matrix(all_A[, , 1:p, m], nrow=d, byrow=FALSE), cbind(I_all, ZER_all)), numeric((d*p)^2)), dim=c(d*p, d*p, M))
}


#' @title Change parametrization of a parameter vector
#'
#' @description \code{change_parametrization} changes the parametrization of the given parameter
#'   vector to \code{change_to}.
#'
#' @inheritParams loglikelihood
#' @inheritParams form_boldA
#' @param change_to either "intercept" or "mean" specifying to which parametrization it should be switched to.
#'   If set to \code{"intercept"}, it's assumed that \code{params} is mean parametrized, and if set to \code{"mean"}
#'   it's assumed that \code{params} is intercept parametrized.
#' @return Returns parameter vector described in \code{params}, but with parametrization changed from intercept to mean
#'   (when \code{change_to == "mean"}) or from mean to intercept (when \code{change_to == "intercept"}).
#' @details Parametrization cannot be changed for models with mean constraints constraints!
#' @section Warning:
#'  No argument checks!
#' @keywords internal

change_parametrization <- function(p, M, d, params, AR_constraints=NULL, mean_constraints=NULL,
                                   change_to=c("intercept", "mean")) {
  stopifnot(is.null(mean_constraints))
  change_to <- match.arg(change_to)
  re_params <- params
  if(!is.null(AR_constraints)) {
    stop("AR_constraints not yet implemented to change_parametrization")
    # Create reform_constrained pars and call if here
  }
  Id <- diag(nrow=d)
  all_A <- pick_allA(p=p, M=M, d=d, params=params)
  all_phi0_or_mu <- pick_phi0(M=M, d=d, params=params)

  # Calculate means/intercepts and insert them to re_params
  if(change_to == "mean") { # params has original parametrization with intercept
    re_params[1:(M*d)] <- vapply(1:M, function(m) solve(Id - rowSums(all_A[, , , m, drop=FALSE], dims=2), all_phi0_or_mu[,m]),
                                 numeric(d))
  } else { # mean parameters instead of phi0
    re_params[1:(M*d)] <- vapply(1:M, function(m) (Id - rowSums(all_A[, , , m, drop=FALSE], dims=2))%*%all_phi0_or_mu[,m],
                                 numeric(d))
  }
  re_params
}


#' @title Sort regimes in parameter vector according to mixing weights into a decreasing order
#'
#' @description \code{sort_regimes} sorts regimes in the parameter vector according to
#'   the transition weight parameters.
#'
#' @inheritParams loglikelihood
#' @details Constrained parameter vectors are not supported. Currently only reduced form models
#'   are supported.
#' @return Returns sorted parameter vector of the form described for the argument \code{params},
#'   with the regimes sorted so that...
#'   \describe{
#'     \item{If \code{weight_function == "relative_dens"}:}{\eqn{\alpha_{1}>...>\alpha_{M}}.}
#'     \item{If \code{weight_function == "logit"}:}{NOT YET IMPLEMENTED}
#'   }
#' @keywords internal

sort_regimes <- function(p, M, d, params, weight_function=c("relative_dens", "logit"), cond_dist=c("Gaussian", "Student"),
                         identification=c("reduced_form", "recursive", "heteroskedasticity")) {
  if(M == 1) return(params) # Nothing to sort
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  identification <- match.arg(identification)
  if(identification != "reduced_form") stop("Structural models not yet implemented to sort_regimes!")
  if(cond_dist != "Gaussian") stop("Only Gaussian cond_dist is implemented to sort_regimes")

  all_weightpars <- pick_weightpars(p=p, M=M, d=d, params=params, weight_function=weight_function,
                                    cond_dist=cond_dist)
  if(weight_function == "relative_dens") {
    new_order <- order(all_weightpars, decreasing=TRUE)
    if(all(new_order == 1:M)) {
      return(params)
    }
    new_weightpars <- all_weightpars[new_order][-M]
  } else {
    stop("Only relative_dens weight function is implementesd to sort_regimes!")
  }


  all_phi0 <- pick_phi0(M=M, d=d, params=params)
  all_A <- matrix(pick_allA(p=p, M=M, d=d, params=params), ncol=M) #matrix(params[(d*M + 1):(d*M + d^2*p*M)], ncol=M)
  all_Omega <- matrix(params[(d*M*(1 + p*d) + 1):(d*M*(1 + p*d) + M*d*(d + 1)/2)], nrow=d*(d + 1)/2, ncol=M)
  all_weightpars <- pick_weightpars(p=p, M=M, d=d, params=params, weight_function=weight_function,
                                    cond_dist=cond_dist)

  c(all_phi0[,new_order], all_A[,new_order], all_Omega[,new_order], new_weightpars) # new_distpars
}
