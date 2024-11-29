#' @title Maximum likelihood estimation of a reduced form or structural STVAR model based on preliminary estimates
#'
#' @description \code{iterate_more} uses a variable metric algorithm to estimate a reduced form or structural STVAR model
#'  (object of class \code{'stvar'}) based on preliminary estimates.
#'
#' @inheritParams fitSTVAR
#' @inheritParams STVAR
#' @param stvar an object of class \code{'stvar'}, created by, e.g., \code{fitSTVAR} or \code{fitSSTVAR}.
#' @param h the step size used in the central difference approximation of the gradient of the log-likelihood function, so
#'   \code{h} should be a small positive real number.
#' @param print_trace should the trace of the optimization algorithm be printed?
#' @details The purpose of \code{iterate_more} is to provide a simple and convenient tool to finalize
#'   the estimation when the maximum number of iterations is reached when estimating a STVAR model
#'   with the main estimation function \code{fitSTVAR} or \code{fitSSTVAR}.
#' @inherit STVAR return
#' @seealso \code{\link{fitSTVAR}}, \code{\link{STVAR}}, \code{\link[stats]{optim}},
#'  \code{\link{swap_B_signs}}, \code{\link{reorder_B_columns}}
#' @inherit fitSTVAR references
#' @examples
#' \donttest{
#' ## These are long running examples that take approximately 20 seconds to run.
#'
#' # Estimate two-regime Gaussian STVAR p=1 model with the weighted relative stationary densities
#' # of the regimes as the transition weight function, but only 5 iterations of the variable matrix
#' # algorithm:
#' fit12 <- fitSTVAR(gdpdef, p=1, M=2, nrounds=1, seeds=1, ncores=1, maxit=5)
#'
#' # The iteration limit was reached, so the estimate is not local maximum.
#' # The gradient of the log-likelihood function:
#' get_foc(fit12) # Not close to zero!
#'
#' # So, we run more iterations of the variable metric algorithm:
#' fit12 <- iterate_more(fit12)
#'
#' # The gradient of the log-likelihood function after iterating more:
#' get_foc(fit12) # Close (enough) to zero!
#' }
#' @export

iterate_more <- function(stvar, maxit=100, h=1e-3, penalized, penalty_params, allow_unstab, calc_std_errors=TRUE, print_trace=TRUE) {
  check_stvar(stvar)
  stopifnot(maxit %% 1 == 0 & maxit >= 1)
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  data <- stvar$data
  weight_function <- stvar$model$weight_function
  cond_dist <- stvar$model$cond_dist
  parametrization <- stvar$model$parametrization
  identification <- stvar$model$identification
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  weight_constraints <- stvar$model$weight_constraints
  B_constraints <- stvar$model$B_constraints
  if(missing(penalized)) {
    penalized <- stvar$penalized
  } else {
    stopifnot(is.logical(penalized))
  }
  if(missing(penalty_params)) {
    penalty_params <- stvar$penalty_params
  } else {
    stopifnot(is.numeric(penalty_params) && length(penalty_params) == 2 && all(penalty_params >= 0) && penalty_params[1] < 1)
  }
  if(missing(allow_unstab)) {
    allow_unstab <- stvar$allow_unstab
  } else {
    stopifnot(is.logical(allow_unstab))
  }
  weightfun_pars <- check_weightfun_pars(data=data, p=p, M=M, d=d, weight_function=weight_function,
                                         weightfun_pars=stvar$model$weightfun_pars, cond_dist=cond_dist)
  minval <- get_minval(stvar$data)
  npars <- length(stvar$params)
  stopifnot(h > 0)

  # Function to optimize
  loglik_fn <- function(params) {
    tryCatch(loglikelihood(data=data, p=p, M=M, params=params,
                           weight_function=weight_function, weightfun_pars=weightfun_pars,
                           cond_dist=cond_dist, parametrization=parametrization,
                           identification=identification, AR_constraints=AR_constraints,
                           mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                           B_constraints=B_constraints, to_return="loglik", check_params=TRUE,
                           penalized=penalized, penalty_params=penalty_params,
                           allow_unstab=allow_unstab, minval=minval),
             error=function(e) minval)
  }

  # Gradient of the log-likelihood function using central difference approximation
  I <- diag(rep(1, times=npars))
  loglik_grad <- function(params) {
    vapply(1:npars, function(i1) (loglik_fn(params + I[i1,]*h) - loglik_fn(params - I[i1,]*h))/(2*h), numeric(1))
  }

  res <- optim(par=stvar$params, fn=loglik_fn, gr=loglik_grad, method=c("BFGS"),
               control=list(fnscale=-1, maxit=maxit, trace=ifelse(print_trace, 6, 0)))
  if(res$convergence == 1) message("The maximum number of iterations was reached! Consired iterating more.")

  ret <- STVAR(data=data, p=p, M=M, d=d, params=res$par,
               weight_function=weight_function,
               weightfun_pars=weightfun_pars,
               cond_dist=cond_dist,
               parametrization=parametrization,
               identification=identification,
               AR_constraints=AR_constraints,
               mean_constraints=mean_constraints,
               weight_constraints=weight_constraints,
               B_constraints=B_constraints,
               penalized=penalized,
               penalty_params=penalty_params,
               allow_unstab=allow_unstab,
               calc_std_errors=calc_std_errors)

  ret$all_estimates <- stvar$all_estimates
  ret$all_logliks <- stvar$all_logliks
  ret$which_converged <- stvar$which_converged
  if(!is.null(stvar$which_round)) {
    ret$which_round <- stvar$which_round
    ret$all_estimates[[stvar$which_round]] <- ret$params
    ret$all_logliks[stvar$which_round] <- ret$loglik
    ret$which_converged[stvar$which_round] <- res$convergence == 0
  }
  warn_eigens(ret, allow_unstab=allow_unstab)
  ret
}



#' @title Maximum likelihood estimation of a structural STVAR model based on preliminary estimates from
#'   a reduced form model.
#'
#' @description \code{fitSSTVAR} uses a robust method and a variable metric algorithm to estimate
#'   a structural STVAR model based on preliminary estimates from a reduced form model.
#'
#' @inheritParams fitSTVAR
#' @inheritParams STVAR
#' @param maxit_robust the maximum number of iterations on the first phase robust estimation, if employed.
#' @param stvar a an object of class \code{'stvar'}, created by, e.g., \code{fitSTVAR},
#'   specifying a reduced form or a structural model
#' @param robust_method Should some robust estimation method be used in the estimation before switching
#'   to the gradient based variable metric algorithm? See details.
#' @param identification Which identification should the structural model use?
#'  (see the vignette or the references for details)
#'   \describe{
#'     \item{\code{"recursive"}:}{The usual lower-triangular recursive identification of the shocks via their impact responses.}
#'     \item{\code{"heteroskedasticity"}:}{Identification by conditional heteroskedasticity, which imposes constant relative
#'       impact responses for each shock.}
#'   }
#' @param B_constraints Employ further constraints on the impact matrix?
#'   A \eqn{(d \times d)} matrix with its entries imposing constraints on the impact matrix \eqn{B_t}:
#'   \code{NA} indicating that the element is unconstrained, a positive value indicating strict positive sign constraint,
#'   a negative value indicating strict negative sign constraint, and zero indicating that the element is constrained to zero.
#'   Currently only available for models with \code{identification="heteroskedasticity"} due to the (in)availability of appropriate
#'   parametrizations that allow such constraints to be imposed.
#' @param B_pm_reg an integer between \eqn{1} and \eqn{M} specifying the regime the permutations and sign changes of \eqn{B_m}
#'   specified in the arguments \code{B_perm} and \code{B_signs} are applied to.
#' @param B_perm a numeric vector of length \eqn{d} specifying the permutation of the columns of the impact matrix \eqn{B_m}
#'   of a single regime specified in the argument \code{B_pm_reg} prior to re-estimating the model. Applicable only for models
#'   with \code{cond_dist = "ind_Student"} or \code{"ind_skewed_t"}.
#' @param B_signs a numeric vector specifying the columns of the impact matrix of a single regime specified in the argument
#'   \code{B_pm_reg} that should be multiplied by -1 \strong{prior} to reordering them according to \code{B_perm} (if specified).
#'    Applicable only for models with \code{cond_dist = "ind_Student"} or \code{"ind_skewed_t"}.
#' @details When the structural model does not impose overidentifying constraints, it is directly
#'   obtained from the reduced form model, and estimation is not required. When overidentifying constraints
#'   are imposed, the model is estimated subject to the constraints.
#'
#'   Using the robust estimation method before switching to the variable metric can be useful if the initial
#'   estimates are not very close to the ML estimate of the structural model, as the variable metric algorithm
#'   (usually) converges to a nearby local maximum or saddle point. However, if the initial estimates are far from
#'   the ML estimate, the resulting solution is likely local only due to the complexity of the model. Note that
#'   Nelder-Mead algorithm is much faster than SANN but can get stuck at a local solution.
#'   This is particularly the case when the imposed overidentifying restrictions are such that the unrestricted
#'   estimate is not close to satisfying them. Nevertheless, in most practical cases, the model is just identified
#'   and estimation is not required, and often reasonable overidentifying constraints are close to the unrestricted estimate.
#'
#'   Employs the estimation function \code{optim} from the package \code{stats} that implements the optimization
#'   algorithms. See \code{?optim} for the documentation on the optimization methods.
#'
#'   The arguments \code{B_pm_reg}, \code{B_perm}, and \code{B_signs} can be used to explore estimates based various orderings
#'   and sign changes of the columns of the impact matrices \eqn{B_m} of specific regimes. This can be useful in the presence
#'   of weak identification with respect to the ordering or signs of the columns \eqn{B_2,...,B_M} (see Virolainen 2024).
#' @inherit STVAR return
#' @seealso \code{\link{fitSTVAR}}, \code{\link{STVAR}}, \code{\link[stats]{optim}}
#' @references
#' \itemize{
#'    \item Kilian L., Lütkepohl H. 20017. Structural Vector Autoregressive Analysis. 1st edition.
#'    \emph{Cambridge University Press}, Cambridge.
#'    \item Lütkepohl H., Netšunajev A. 2017. Structural vector autoregressions with smooth transition in variances.
#'      \emph{Journal of Economic Dynamics & Control}, \strong{84}, 43-57.
#'    \item Virolainen S. 2024. Identification by non-Gaussianity in structural threshold and
#'       smooth transition vector autoregressive models. Unpublished working
#'       paper, available as arXiv:2404.19707.
#'  }
#' @examples
#' \donttest{
#' ## These are long running examples that take approximately 1 minute to run.
#'
#' ## Estimate first a reduced form Gaussian STVAR p=3, M=2 model with the weighted relative
#' # stationary densities of the regimes as the transition weight function, and the means and
#' # AR matrices constrained to be identical across the regimes:
#' fit32cm <- fitSTVAR(gdpdef, p=3, M=2, AR_constraints=rbind(diag(3*2^2), diag(3*2^2)),
#'   weight_function="relative_dens", mean_constraints=list(1:2), parametrization="mean",
#'   nrounds=1, seeds=1, ncores=1)
#'
#' # Then, we estimate/create various structural models based on the reduced form model.
#' # Create a structural model with the shocks identified recursively:
#' fit32cms_rec <- fitSSTVAR(fit32cm, identification="recursive")
#'
#' # Create a structural model with the shocks identified by conditional heteroskedasticity:
#' fit32cms_hetsked <- fitSSTVAR(fit32cm, identification="heteroskedasticity")
#' fit32cms_hetsked # Print the estimates
#'
#' # Estimate a structural model with the shocks identified by conditional heteroskedasticity
#' # and overidentifying constraints imposed on the impact matrix: positive diagonal element
#' # and zero upper right element:
#' fit32cms_hs2 <- fitSSTVAR(fit32cm, identification="heteroskedasticity",
#'  B_constraints=matrix(c(1, NA, 0, 1), nrow=2))
#'
#' # Estimate a structural model with the shocks identified by conditional heteroskedasticity
#' # and overidentifying constraints imposed on the impact matrix: positive diagonal element
#' # and zero off-diagonal elements:
#' fit32cms_hs3 <- fitSSTVAR(fit32cms_hs2, identification="heteroskedasticity",
#'  B_constraints=matrix(c(1, 0, 0, 1), nrow=2))
#'
#' # Estimate first a reduced form two-regime Threshold VAR p=1 model with
#' # with independent skewed t shocks, and the first lag of the second variable
#' # as the switching variable, and AR matrices constrained to be identical
#' # across the regimes:
#' fit12c <- fitSTVAR(gdpdef, p=1, M=2, cond_dist="ind_skewed_t",
#'  AR_constraints=rbind(diag(1*2^2), diag(1*2^2)), weight_function="threshold",
#'  weightfun_pars=c(2, 1), nrounds=1, seeds=1, ncores=1)
#'
#' # Due to the independent non-Gaussian shocks, the structural shocks are readily
#' # identified. The following returns the same model but marked as structural
#' # with the shocks identified by non-Gaussianity:
#' fit12c <- fitSSTVAR(fit12c)
#'
#' # Estimate a model based on a reversed ordering of the columns of the impact matrix B_2:
#' fit12c2 <- fitSSTVAR(fit12c, B_pm_reg=2, B_perm=c(2, 1))
#'
#' # Estimate a model based on reversed signs of the second column of B_2 and reversed
#' # ordering of the columns of B_2:
#' fit12c3 <- fitSSTVAR(fit12c, B_pm_reg=2, B_perm=c(2, 1), B_signs=2)
#' }
#' @export

fitSSTVAR <- function(stvar, identification=c("recursive", "heteroskedasticity", "non-Gaussianity"), B_constraints=NULL,
                      B_pm_reg=NULL, B_perm=NULL, B_signs=NULL, maxit=1000, maxit_robust=1000,
                      robust_method=c("Nelder-Mead", "SANN", "none"), print_res=TRUE, calc_std_errors=TRUE) {
  check_stvar(stvar)
  stopifnot(maxit %% 1 == 0 & maxit >= 1)
  stopifnot(maxit_robust %% 1 == 0 & maxit_robust >= 1)
  robust_method <- match.arg(robust_method)
  if(length(identification > 1) && match.arg(identification) == "recursive"
     && stvar$model$cond_dist %in% c("ind_Student", "ind_skewed_t")) {
    identification <- "non-Gaussianity"
  } else {
    identification <- match.arg(identification)
  }
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  data <- stvar$data
  params <- stvar$params
  weight_function <- stvar$model$weight_function
  cond_dist <- stvar$model$cond_dist
  parametrization <- stvar$model$parametrization
  old_identification <- stvar$model$identification
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  weight_constraints <- stvar$model$weight_constraints
  old_B_constraints <- stvar$model$B_constraints
  penalized <- stvar$penalized
  penalty_params <- stvar$penalty_params
  allow_unstab <- stvar$allow_unstab
  weightfun_pars <- check_weightfun_pars(data=data, p=p, M=M, d=d, weight_function=weight_function,
                                         weightfun_pars=stvar$model$weightfun_pars, cond_dist=cond_dist)

  # Check B_pm_reg, B_perm and B_signs
  if(!is.null(B_pm_reg)) {
    if(length(B_pm_reg) != 1 || !B_pm_reg %in% 1:M) {
      stop("B_pm_reg should be an integer between 1 and M")
    }
  }
  if(!is.null(B_perm)) {
    if(cond_dist != "ind_Student" && cond_dist != "ind_skewed_t") {
      stop("B_perm is only available for models with cond_dist='ind_Student' or 'ind_skewed_t'")
    } else if(length(B_perm) != d) {
      stop("The length of B_perm should be equal to the number of variables")
    } else if(!all(sort(B_perm, decreasing=FALSE) == 1:d)) {
      stop("B_perm should contain the integers from 1 to d, each exactly once")
    }
  }
  if(!is.null(B_signs)) {
    if(cond_dist != "ind_Student" && cond_dist != "ind_skewed_t") {
      stop("B_signs are only available for models with cond_dist='ind_Student' or 'ind_skewed_t'")
    } else if(length(B_signs) > d) {
      stop("The length of B_signs should be at most equal to the number of variables")
    }
    B_signs <- unique(B_signs) # Take only unique elements in B_signs
    if(!all(B_signs %in% 1:d)) {
      stop("B_signs should only contain integers between 1 and d")
    }
  }
  if(!is.null(B_constraints) || !is.null(old_B_constraints)) {
    if(!is.null(B_perm) || !is.null(B_signs)) {
      stop(paste("B_perm and B_signs cannot be used with B_constraints! You should first use B_perm and B_signs,",
                 "and then you can impose further constraints on the impact matrix by running fitSSTVAR again."))
    }
  }

  minval <- get_minval(stvar$data)
  old_npars <- length(stvar$params) # old npars, B_constraints not adjusted!

  # Check identification
  if(identification != "non-Gaussianity" && cond_dist %in% c("ind_Student", "ind_skewed_t")) {
     message(paste("Only identification by non-Gaussianity is supported with independent",
                   ifelse(cond_dist == "ind_Student", "Student's", "skewed"), "t-distributed errors.",
                   "Using identification by non-Gaussianity. Note that recursive structure and other overidentifying",
                   "constraints can by imposed by the argument B_constraints."))
    identification <- "non-Gaussianity"
  }

  if(old_identification != "reduced_form" && cond_dist != "ind_Student" && cond_dist != "ind_skewed_t") {
    if(old_identification == "recursive") {
      warning(paste("Recursively identified model supplied (not possible to impose overidentifying restrictions,",
                    "so the original model is returned"))

      return(stvar)
    } else if(old_identification != identification) {
      # Old model structural but identification changes
      message(paste("The type of the identification of a structural model cannot be changed.",
                    "Using the old identification as 'identification'."))
      identification <- old_identification
    }
  }

  if(identification == "recursive") {
    # Old identification reduced form, recursively identified model should be returned.
    if(!is.null(B_constraints)) warning("B_constraints cannot be imposed on recursively identified models")

    return(STVAR(data=data, p=p, M=M, d=d, params=params, cond_dist=cond_dist,
                 weight_function=weight_function, weightfun_pars=weightfun_pars,
                 parametrization=parametrization, identification=identification,
                 AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                 weight_constraints=weight_constraints, B_constraints=NULL,
                 penalized=penalized, penalty_params=penalty_params, allow_unstab=allow_unstab,
                 calc_std_errors=calc_std_errors))
  }

  # Check the new B_constraints
  check_constraints(data=data, p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    parametrization=parametrization, identification=identification,
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=B_constraints)

  ## The numbers of parameters (needed to create the new initial param vector)
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
  # Covmatpars
  if(cond_dist == "ind_Student" || cond_dist == "ind_skewed_t" || identification == "non-Gaussiniaty") {
    if(is.null(old_B_constraints)) {
      n_zeros_in_B <- 0
    } else {
      n_zeros_in_B <- sum(old_B_constraints == 0, na.rm=TRUE)
    }
    n_covmat_pars <- M*d^2 - n_zeros_in_B
  } else if(old_identification %in% c("reduced_form", "recursive")) {
    n_covmat_pars <- M*d*(d + 1)/2 # No B_constraints available here
    n_zeros_in_W <- 0
  } else { # identification == "heteroskedasticity
    if(is.null(old_B_constraints)) {
      n_zeros_in_W <- 0
    } else {
      n_zeros_in_W <- sum(old_B_constraints == 0, na.rm=TRUE)
    }
    n_covmat_pars <- d^2 + d*(M - 1) - n_zeros_in_W
  }
  if(weight_function == "exogenous") {
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

  if(cond_dist == "Gaussian") {
    n_dist_pars <- 0
  } else if(cond_dist == "Student") {
    n_dist_pars <- 1
  } else if(cond_dist == "ind_Student") {
    n_dist_pars <- d
  } else { # cond_dist == "ind_skewed_t"
    n_dist_pars <- 2*d
  }

  n_pars <- n_mean_pars + n_ar_pars + n_covmat_pars + n_weight_pars + n_dist_pars

  if(cond_dist == "ind_Student" || cond_dist == "ind_skewed_t" || identification == "non-Gaussianity") {
    if(!is.null(B_perm) || !is.null(B_signs)) {
      alt_par <- FALSE # Columns are permutated or signs changed, use normal parametrization
    } else if(!is.null(B_constraints) && length(which(B_constraints != 0)) != 0) {
      alt_par <- FALSE # Sign constraints employed, use normal parametrization
    } else { # no B_perm or B_signs, B_constraints not used or only used for zero constraints
      alt_par <- TRUE # Switch to parametrization with B_1*,...,B_M*
      params <- change_parametrization(p=p, M=M, d=d, params=params, cond_dist=cond_dist,
                                       weight_function=weight_function, weightfun_pars=weightfun_pars,
                                       identification=identification, AR_constraints=AR_constraints,
                                       mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                       B_constraints=B_constraints, change_to="alt")
    }
  } else {
    alt_par <- FALSE
  }

  # Set up initial values using the obtained parameters
  if(identification == "heteroskedasticity") {
    if(M == 1) {
      stop("Identification by heteroskedasticity requires at least two regimes")
    } else if(M == 2 && is.null(B_constraints) && old_identification == "reduced_form") {
      # Nothing to estimate, decompose the error term covariance matrices and create new params
      # and return the new model
      return(get_hetsked_sstvar(stvar, calc_std_errors=calc_std_errors))
    } else if(is.null(old_B_constraints) && is.null(B_constraints) && old_identification == "heteroskedasticity") {
      # The model does not change, return the original model
      return(stvar)
    }
    ## If arrived here, overidentifying constraints are imposed, either via M>2 or B_constraints.
    if(old_identification == "reduced_form") {
      if(M == 2) {
        params <- get_hetsked_sstvar(stvar, calc_std_errors=FALSE)$params # Change params to hetsked params
      } else {
        # M > 2, so needs to obtain initial estimates some other way. We decompose the first two Omegas
        # to obtain W and the first lambdas; then we obtain the lambda_m by decomposing the Omega_1 and Omega_m.
        # Maybe not the optimal solution, but easily scales to any number of regimes and gives the most weight
        # to the first regime, and the second most to the third regime.

        # First, pick reduced form error covariance matrices:
        params_std <- reform_constrained_pars(p=p, M=M, d=d, params=params, cond_dist=cond_dist,
                                              weight_function=weight_function, weightfun_pars=weightfun_pars,
                                              identification=old_identification,
                                              AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                              weight_constraints=weight_constraints, B_constraints=old_B_constraints)
        all_Omegas <- pick_Omegas(p=p, M=M, d=d, params=params_std, cond_dist=cond_dist, identification=old_identification)

        # W and lambda_2:
        W_and_lambdas <- numeric(d^2 + (M - 1)*d)
        W_and_lambdas[1:(d^2 + d)] <- diag_Omegas(Omega1=all_Omegas[, , 1], Omega2=all_Omegas[, , 2])

        # lambda_3,...,lambda_M:
        for(i1 in 3:M) {
          tmp <- diag_Omegas(Omega1=all_Omegas[, , 1], Omega2=all_Omegas[, , i1])
          W_and_lambdas[(d^2 + d*(i1 - 2) + 1):(d^2 + d*(i1 - 1))] <- tmp[(d^2 + 1):(d^2 + d)]
        }
      }
    }

    # Obtain the W pars
    if(old_identification == "reduced_form" && M > 2) {
      oldWpars <- W_and_lambdas[1:d^2]
    } else {
      oldWpars <- params[(n_mean_pars + n_ar_pars + 1):(n_mean_pars + n_ar_pars + d^2 - n_zeros_in_W)]
    }
    if(is.null(old_B_constraints)) {
      old_B_constraints <- matrix(NA, nrow=d, ncol=d) # No constraints imposed here but avoids null problems below
    }
    if(is.null(B_constraints)) {
      B_constraints <- matrix(NA, nrow=d, ncol=d) # No constraints imposed here but avoids null problems below
    }

    n_sign_changes <- 0
    oldWpars_inW <- matrix(0, nrow=d, ncol=d) # Fill in the parameters including non-parametrized zeros:
    oldWpars_inW[old_B_constraints != 0 | is.na(old_B_constraints)] <- oldWpars

    # Go through the matrices old_B_constraints, B_constraints, and change the oldWpars accordingly
    newWpars <- matrix(NA, nrow=d, ncol=d)
    for(i1 in 1:d) { # i1 marks the row
      for(i2 in 1:d) { # i2 marks the column
        so <- sign(old_B_constraints[i1, i2])
        sn <- sign(B_constraints[i1, i2])
        if(is.na(so)) { # If no constraint, just handle it as if there was sign constraint of the estimate's sign
          so <- sign(oldWpars_inW[i1, i2])
        }
        if(is.na(sn) && so != 0) { # No new constraint, old constraint is not zero
          newWpars[i1, i2] <- oldWpars_inW[i1, i2] # Insert old param
        } else if(is.na(sn) && so == 0) { # No new constraints, old constraint is zero
          newWpars[i1, i2] <- 1e-6  # Insert close to zero value: not exactly zero to avoid the parameter disappearing in Wvec
        } else if(so == sn) { # Same constraint in new and old W
          newWpars[i1, i2] <- oldWpars_inW[i1, i2] # Insert old param
        } else if(sn == 0) { # Zero constraint in new W
          newWpars[i1, i2] <- 0 # Insert zero (which will be removed as it is not parametrized)
        } else if(so > sn) { # sn must be negative, so could be zero or positive
          newWpars[i1, i2] <- -0.05 # Insert small negative number
          n_sign_changes <- n_sign_changes + 1
        } else { # It must be that so < sn, which implies sn > 0, while so could be zero or negative
          newWpars[i1, i2] <- 0.05 # Insert s mall positive number
          n_sign_changes <- n_sign_changes + 1
        }
      }
    }
    newWpars <- Wvec(newWpars) # New initial parameters for W


    # New initial params
    if(old_identification == "reduced_form" && M > 2) {
      new_params <- c(params[1:(n_mean_pars + n_ar_pars)], newWpars, W_and_lambdas[(d^2 + 1):(d^2 + (M - 1)*d)],
                      params[(n_mean_pars + n_ar_pars + n_covmat_pars + 1):length(params)])
    } else { # Het.sked model already originally or M=2
      new_params <- c(params[1:(n_mean_pars + n_ar_pars)], newWpars,
                      params[(n_mean_pars + n_ar_pars + length(oldWpars) + 1):length(params)])
    }
    if(n_sign_changes > 0) {
      message(paste0("There was ", n_sign_changes,
              " sign changes in the impact matrix when creating preliminary estimates for the new model. ",
              "The sign changes make the results more unrealiable. To obtain more reliable results, consider using ",
              "the function 'swap_B_signs' to create a model that has the sign constraints readily satisfied and ",
              "then applying this function.\n"))
    }
  } else { # identification == "non-Gaussianity"
    if(!is.null(B_perm) || !is.null(B_signs)) { # Permutate or swap signs of the columns of B_1,...,B_M
      # Construct initial estimates after reordering the columns or swapping the signs

      # First obtain the old estimates (each B_m has d^2 parameters here, as no B_constraints are allowed)
      stopifnot(n_covmat_pars == M*d^2)
      all_Bm <- array(params[(n_mean_pars + n_ar_pars + 1):(n_mean_pars + n_ar_pars + n_covmat_pars)], dim=c(d, d, M))

      # Take the B_m alter:
      Bm <- all_Bm[, , B_pm_reg]

      # Change the signs of the columns of B_m
      if(!is.null(B_signs)) {
        Bm[, B_signs] <- -Bm[, B_signs]
      }

      # Permutate the columns of B_m
      if(!is.null(B_perm)) {
        Bm <- Bm[, B_perm]
      }

      # Fill in the new B_m to the new_params, which is new initial parameter vector for estimation
      all_Bm[, , B_pm_reg] <- Bm
      new_params <- c(params[1:(n_mean_pars + n_ar_pars)], all_Bm, params[(n_mean_pars + n_ar_pars + n_covmat_pars + 1):length(params)])
    } else { # B_constraints imposed or changed
      if(is.null(old_B_constraints) && is.null(B_constraints)) {
        return(stvar) # Nothing to do here, return the original model
      } else {
        ## If arrived here, there are new B_constraints or the old ones are relaxed.
        if(is.null(old_B_constraints)) {
          old_B_constraints <- matrix(NA, nrow=d, ncol=d) # No constraints imposed here but avoids null problems below
        }
        if(is.null(B_constraints)) {
          B_constraints <- matrix(NA, nrow=d, ncol=d) # No constraints imposed here but avoids null problems below
        }
      }
      ## Construct initial estimates for the new model
      # First obtain the old estimates, including the non-parametrized zeros:
      params_std <- reform_constrained_pars(p=p, M=M, d=d, params=params, cond_dist=cond_dist,
                                            weight_function=weight_function, weightfun_pars=weightfun_pars,
                                            identification=old_identification,
                                            AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                            weight_constraints=weight_constraints, B_constraints=old_B_constraints)
      old_all_B <- pick_Omegas(p=p, M=M, d=d, params=params_std, cond_dist=cond_dist, identification=old_identification)
      new_all_B <- array(NA, dim=c(d, d, M))

      ## Construct the new B matrices based in the new set of constraints in B_constraints.
      n_sign_changes <- 0 # Track the number of sign changes

      # Go through the matrices old_B_constraints, B_constraints, and change the impact matrices accordingly
      for(m in 1:M) {
        for(i1 in 1:d) { # i1 marks the row
          for(i2 in 1:d) { # i2 marks the column
            so <- sign(old_B_constraints[i1, i2])
            sn <- sign(B_constraints[i1, i2])
            if(is.na(so)) { # If no constraint, just handle it as if there was sign constraint of the estimate's sign
              so <- sign(old_all_B[i1, i2, m])
            }
            if(is.na(sn) && so != 0) { # No new constraint, old constraint is not zero
              new_all_B[i1, i2, m] <- old_all_B[i1, i2, m] # Insert old param
            } else if(is.na(sn) && so == 0) { # No new constraints, old constraint is zero
              new_all_B[i1, i2, m] <- 1e-6  # Insert close to zero value: not exactly zero to avoid the parameter disappearing in Wvec
            } else if(so == sn) { # Same constraint in new and old W
              new_all_B[i1, i2, m] <- old_all_B[i1, i2, m] # Insert old param
            } else if(sn == 0) { # Zero constraint in new W
              new_all_B[i1, i2, m] <- 0 # Insert zero (which will be removed as it is not parametrized)
            } else if(so > sn) { # sn must be negative, so could be zero or positive
              new_all_B[i1, i2, m] <- -0.05 # Insert small negative number
              n_sign_changes <- n_sign_changes + 1
            } else { # It must be that so < sn, which implies sn > 0, while so could be zero or negative
              new_all_B[i1, i2, m] <- 0.05 # Insert s mall positive number
              n_sign_changes <- n_sign_changes + 1
            }
          }
        }
      }
      new_B_pars <- numeric(0)
      for(m in 1:M) {
        new_B_pars <- c(new_B_pars, Wvec(new_all_B[, , m])) # Remove non-parametrized zeros
      }

      new_params <- c(params[1:(n_mean_pars + n_ar_pars)], new_B_pars,
                      params[(n_mean_pars + n_ar_pars + M*d^2 - M*n_zeros_in_B + 1):length(params)])

      if(n_sign_changes > 0) {
        message(paste0("There was in total of ", n_sign_changes,
                       " sign changes in the impact matrices of the regimes when creating preliminary estimates for the new model. ",
                       "The sign changes make the results more unrealiable. To obtain more reliable results, consider using ",
                       "the function 'swap_B_signs' to create a model that has the sign constraints readily satisfied and ",
                       "then applying this function.\n"))
      }
    }
  }

  ### Check ups for the initial estimates of the new model
  new_loglik <- loglikelihood(data=data, p=p, M=M, params=new_params, cond_dist=cond_dist,
                              weight_function=weight_function, weightfun_pars=weightfun_pars,
                              parametrization=parametrization, identification=identification,
                              AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                              weight_constraints=weight_constraints, B_constraints=B_constraints,
                              penalized=penalized, penalty_params=penalty_params,
                              allow_unstab=allow_unstab, to_return="loglik", check_params=TRUE, minval=NULL,
                              alt_par=alt_par)

  if(is.null(new_loglik)) {
    message("Problem with the new parameter vector - try different B_constraints?\n See:")
    check_params(data=data, p=p, M=M, d=d, params=new_params, cond_dist=cond_dist,
                 weight_function=weight_function, weightfun_pars=weightfun_pars,
                 parametrization=parametrization, identification=identification,
                 AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                 weight_constraints=weight_constraints, B_constraints=B_constraints,
                 allow_unstab=allow_unstab)
  }

  if(print_res) {
    message("The log-likelihood of the supplied model:    ", round(c(stvar$loglik), 3),
            "\nConstrained log-likelihood prior estimation: ", round(new_loglik, 3))
  }


  ### Estimation
  if(is.null(data)) {
    if(alt_par) {
      new_params <- change_parametrization(p=p, M=M, d=d, params=new_params, cond_dist=cond_dist,
                                           weight_function=weight_function, weightfun_pars=weightfun_pars,
                                           identification=identification, AR_constraints=AR_constraints,
                                           mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                           B_constraints=B_constraints, change_to="orig") # Change back to orig
    }

    # No data, so returns the model with preliminary param values without estimation.
    message("There is no data for the model, so no estimation was done")
    return(STVAR(data=data, p=p, M=M, d=d, params=new_params, cond_dist=cond_dist,
                 weight_function=weight_function, weightfun_pars=weightfun_pars,
                 parametrization=parametrization, identification=identification,
                 AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                 weight_constraints=weight_constraints, B_constraints=B_constraints,
                 penalized=penalized, penalty_params=penalty_params, allow_unstab=allow_unstab,
                 calc_std_errors=FALSE))
  }


  # Function to optimize
  loglik_fn <- function(params) {
    tryCatch(loglikelihood(data=data, p=p, M=M, params=params,
                           weight_function=weight_function, weightfun_pars=weightfun_pars,
                           cond_dist=cond_dist, parametrization=parametrization,
                           identification=identification, AR_constraints=AR_constraints,
                           mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                           B_constraints=B_constraints, to_return="loglik", check_params=TRUE,
                           penalized=penalized, penalty_params=penalty_params, allow_unstab=allow_unstab,
                           minval=minval, alt_par=alt_par),
             error=function(e) minval)
  }

  ## Gradient of the log-likelihood function using central difference approximation
  npars <- length(new_params)
  h <- 1e-3
  I <- diag(rep(1, times=npars))
  loglik_grad <- function(params) {
    vapply(1:npars, function(i1) (loglik_fn(params + I[i1,]*h) - loglik_fn(params - I[i1,]*h))/(2*h), numeric(1))
  }

  ## Estimation by robust method
  if(robust_method != "none") {
    robust_results <- stats::optim(par=new_params, fn=loglik_fn, gr=loglik_grad, method=robust_method,
                                   control=list(fnscale=-1, maxit=maxit_robust))
    new_params <- robust_results$par
    if(print_res) message(paste0("The log-likelihood after robust estimation:  ", round(robust_results$value, 3)))
  }

  ## Estimation by the variable metric algorithm...
  final_results <- stats::optim(par=new_params, fn=loglik_fn, gr=loglik_grad, method="BFGS",
                                control=list(fnscale=-1, maxit=maxit))
  new_params <- final_results$par
  if(print_res) message(paste0("The log-likelihood after final estimation:   ", round(final_results$value, 3)))

  if(alt_par) {
    new_params <- change_parametrization(p=p, M=M, d=d, params=new_params, cond_dist=cond_dist,
                                         weight_function=weight_function, weightfun_pars=weightfun_pars,
                                         identification=identification, AR_constraints=AR_constraints,
                                         mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                         B_constraints=B_constraints, change_to="orig") # Change back to orig
  }

  # Return the estimated model
  STVAR(data=data, p=p, M=M, d=d, params=new_params, cond_dist=cond_dist,
        weight_function=weight_function, weightfun_pars=weightfun_pars,
        parametrization=parametrization, identification=identification,
        AR_constraints=AR_constraints, mean_constraints=mean_constraints,
        weight_constraints=weight_constraints, B_constraints=B_constraints,
        calc_std_errors=calc_std_errors)
}



#' @title Internal estimation function for estimating STVAR model when bootstrapping confidence
#'   bounds for IRFs in \code{linear_IRF}
#'
#' @description \code{fitbsSSTVAR} uses a robust method and a variable metric algorithm to estimate
#'   a structural STVAR model based on preliminary estimates.
#'
#' @inheritParams loglikelihood
#' @inheritParams fitSSTVAR
#' @param seed the seed for the random number generator (relevant when using SANN).
#' @details Used internally in the functions \code{linear_IRF} for estimating the model in each bootstrap replication.
#'
#'   Employs the estimation function \code{optim} from the package \code{stats} that implements the optimization
#'   algorithms.
#' @inherit STVAR return
#' @section warning: No argument checks!
#' @seealso \code{\link{linear_IRF}}, \code{\link[stats]{optim}}
#' @inherit fitSSTVAR references
#' @keywords internal

fitbsSSTVAR <- function(data, p, M, params,
                        weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                        weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t"),
                        parametrization=c("intercept", "mean"),
                        identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity"),
                        AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL,
                        other_constraints=NULL, robust_method=c("Nelder-Mead", "SANN", "none"),
                        penalized=FALSE, penalty_params=c(0.05, 0.2), allow_unstab=FALSE, minval=NULL,
                        maxit=1000, maxit_robust=1000, seed=NULL) {
  set.seed(seed)
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  identification <- match.arg(identification)
  robust_method <- match.arg(robust_method)
  minval <- get_minval(data)
  d <- ncol(data)

  if(cond_dist == "ind_Student" || cond_dist == "ind_skewed_t" || identification == "non-Gaussianity") {
    if(!is.null(B_constraints) && length(which(B_constraints != 0)) != 0) {
      alt_par <- FALSE # Sign constraints employed, use normal parametrization
    } else { # B_constraints not used or only used for zero constraints
      alt_par <- TRUE # Switch to parametrization with B_1*,...,B_M*
      params <- change_parametrization(p=p, M=M, d=d, params=params, cond_dist=cond_dist,
                                       weight_function=weight_function, weightfun_pars=weightfun_pars,
                                       identification=identification, AR_constraints=AR_constraints,
                                       mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                       B_constraints=B_constraints, change_to="alt")
    }
  } else {
    alt_par <- FALSE
  }

  # Function to optimize
  loglik_fn <- function(params) {
    tryCatch(loglikelihood(data=data, p=p, M=M, params=params,
                           weight_function=weight_function, weightfun_pars=weightfun_pars,
                           cond_dist=cond_dist, parametrization=parametrization,
                           identification=identification, AR_constraints=AR_constraints,
                           mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                           B_constraints=B_constraints, other_constraints=other_constraints, to_return="loglik",
                           check_params=TRUE, penalized=penalized, penalty_params=penalty_params,
                           allow_unstab=allow_unstab, minval=minval, alt_par=alt_par),
             error=function(e) minval)
  }

  ## Gradient of the log-likelihood function using central difference approximation
  npars <- length(params)
  I <- diag(rep(1, times=npars))
  h <- 1e-3
  loglik_grad <- function(params) {
    vapply(1:npars, function(i1) (loglik_fn(params + I[i1,]*h) - loglik_fn(params - I[i1,]*h))/(2*h), numeric(1))
  }

  ## Estimation by robust method
  if(robust_method != "none") {
    robust_results <- stats::optim(par=params, fn=loglik_fn, gr=loglik_grad, method=robust_method,
                                   control=list(fnscale=-1, maxit=maxit_robust))
   params <- robust_results$par
  }

  ## Estimation by the variable metric algorithm..
  final_results <- stats::optim(par=params, fn=loglik_fn, gr=loglik_grad, method="BFGS",
                                control=list(fnscale=-1, maxit=maxit))
  params <- final_results$par

  # Change back to original parametrization
  if(alt_par) {
    params <- change_parametrization(p=p, M=M, d=d, params=params, cond_dist=cond_dist,
                                     weight_function=weight_function, weightfun_pars=weightfun_pars,
                                     identification=identification, AR_constraints=AR_constraints,
                                     mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                     B_constraints=B_constraints, change_to="orig") # Change back to orig
  }

  ## Return the estimate
  params
}


#' @title Internal estimation function for estimating autoregressive and threshold parameters of
#'   TVAR models by the method of least squares.
#'
#' @description \code{estim_LS} estimates the autoregressive and threshold parameters of TVAR models
#'   by the method of least squares.
#'
#' @inheritParams loglikelihood
#' @param ncores the number CPU cores to be used in parallel computing.
#' @param use_parallel employ parallel computing? If \code{FALSE}, does not print anything.
#' @details Used internally in the multiple phase estimation procedure proposed by Koivisto,
#'  Luoto, and Virolainen (2025). Mean constraints are not supported. Only weight constraints that
#'  specify the threshold parameters as fixed values are supported. Only intercept parametrization is
#'  supported.
#' @return Returns the estimated parameters in a vector of the form
#'  \eqn{(\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M,\alpha}, where
#'  \itemize{
#'     \item{\eqn{\phi_{m,0} = } the \eqn{(d \times 1)} intercept vector of the \eqn{m}th regime.}
#'     \item{\eqn{\varphi_m = (vec(A_{m,1}),...,vec(A_{m,p}))} \eqn{(pd^2 \times 1)}.}
#'     \item{\eqn{\alpha = (r_1,...,r_{M-1})} the \eqn{(M-1\times 1)} vector of the threshold parameters.}
#'  }
#'  For models with...
#'   \describe{
#'     \item{AR_constraints:}{Replace \eqn{\varphi_1,...,\varphi_M} with \eqn{\psi} as described in the argument \code{AR_constraints}.}
#'     \item{weight_constraints:}{If linear constraints are imposed, replace \eqn{\alpha} with \eqn{\xi} as described in the
#'      argument \code{weigh_constraints}. If weight functions parameters are imposed to be fixed values, simply drop \eqn{\alpha}
#'      from the parameter vector.}
#'   }
#' @references
#'  \itemize{
#'    \item Hubrich K., Teräsvirta. T. 2013. Thresholds and Smooth Transitions in Vector Autoregressive Models.
#'      \emph{CREATES Research Paper 2013-18, Aarhus University.}
#'    \item Koivisto T., Luoto J., Virolainen S. 2025. Unpublished working paper.
#'    \item Tsay R. 1998. Testing and Modeling Multivariate Threshold Models.
#'      \emph{Journal of the American Statistical Association}, \strong{93}:443, 1188-1202.
#'    \item Virolainen S. 2024. Identification by non-Gaussianity in structural threshold and
#'       smooth transition vector autoregressive models. Unpublished working
#'       paper, available as arXiv:2404.19707.
#'  }
#' @keywords internal

estim_LS <- function(data, p, M, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                     weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t"),
                     parametrization=c("intercept", "mean"), AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL,
                     penalized=TRUE, penalty_params=c(0.05, 0.2), use_parallel=TRUE, ncores=2) {
  # Checks
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  if(!is.null(mean_constraints)) {
    stop("Mean constraints are not supported by the LS estimation")
  } else if(weight_function != "threshold") {
    stop("Only threshold weight function is supported by the LS estimation")
  } else if(parametrization != "intercept") {
    stop("Only the intercept parametrization is supported by the LS estimation")
  }
  stopifnot(is.numeric(penalty_params) && length(penalty_params) == 2 && all(penalty_params >= 0) && penalty_params[1] < 1)
  stab_tol <- penalty_params[1]
  tuning_par <- penalty_params[2]

  # Check the weight constraints
  if(!is.null(weight_constraints)) {
    if(!all(weight_constraints[[1]] == 0)) {
      stop(paste("Only such weight_constraints that specify the threshold parameters some known fixed values",
                 "are supported in the least squares estimation."))
    }
  }
  data <- check_data(data, p=p)
  d <- ncol(data)
  weightfun_pars <- check_weightfun_pars(data=data, p=p, d=d, M=M, weight_function=weight_function,
                                         weightfun_pars=weightfun_pars, cond_dist=cond_dist)
  check_constraints(data=data, p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    parametrization=parametrization, identification="reduced_form",
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=NULL)

  # Obtain relevant statistics
  n_obs <- nrow(data)
  T_obs <- n_obs - p
  pars_per_regime <- d + ifelse(is.null(AR_constraints), p*d^2, ncol(AR_constraints)/M) + d^2
  T_min <- 2/d*pars_per_regime # Minimum number of obs in each regime
  if(T_obs/M < T_min) { # Try smaller T_min
    T_min <- 1.5/d*pars_per_regime
    message("The number of observations relative to the number of parameters is small, consider decreasing p or M.")
    if(T_obs/M < T_min) {
      stop("The number of observations is too small for reasonable estimation. Decrease the order p or the number of regimes M.")
    }
  }

  ################################
  ## Create estim etc functions ##
  ################################

  ## Least squares estimation function given thresholds for models without AR constraints
  LS_without_AR_constraints <- function(thresholds) {
    # threshold = length M-1 vector of the thresholds r_1,...,r_{M-1}; if M=1 anything is ok
    # Other arguments are taken from the parent environment.

    # Storage for the estimates
    all_intercepts <- matrix(NA, nrow=d, ncol=M) # (d x M), [, m] for regime m
    all_AR_mats <- array(NA, dim=c(d, d*p, M)) # [, , m] for A_{m,1} : ... : A_{m,p}

    # Storage for the sums of squares of residuals
    all_rss <- numeric(M) # The sum of squares of residuals for each regime

    # In Y, i:th row denotes the vector \bold{y_{i-1}} = (y_{i-1},...,y_{i-p}) (dpx1),
    # assuming the observed data is y_{-p+1},...,y_0,y_1,...,y_{T}. The last row is for
    # the vector (y_{T},...,y_{T-p}).
    Y <- reform_data(data, p) # (T_obs + 1 x dp)
    Y2 <- Y[1:T_obs, , drop=FALSE] # Last row removed; not needed when only lagged observations used

    # The transition weights: (T x M) matrix, the t:th row is for the time point t and the m:th column is for the regime m.
    alpha_mt <- get_alpha_mt(data, Y2=Y2, p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                             weightpars=thresholds) # Transition weights (T x M), t:th row

    # Go through the regimes
    for(m in 1:M) {
      # Obtain the time periods during which regime m prevails
      m_periods <- which(abs(alpha_mt[, m] - 1) < 1e-6) # The time periods t when regime m prevails

      # Create the matrix Y for regime m; the first p values in data are the initial valus
      Y_m <- data[m_periods + p, , drop=FALSE] # (T_m x d)]

      # Create the matrix X for regime m
      X_m <- cbind(1, Y2[m_periods, , drop=FALSE]) # (T_m x d(p+1))

      # Calculate the lest squares estimates
      C_m <- tryCatch(solve(crossprod(X_m, X_m), crossprod(X_m, Y_m)), # (d x d(p+1)), [\phi_{m,0} : A_{m,1} : ... : A_{m,p}]
                      error=function(e) matrix(0, nrow=d*p + 1, ncol=d)) # zero estimates are legal but bad, dummy estimates
      tC_m <- t(C_m)

      # Store the estimates
      all_intercepts[, m] <- tC_m[, 1]
      all_AR_mats[, , m] <- tC_m[, -1]

      # Calculate the sum of squares of residuals
      U_m <- Y_m - X_m%*%C_m # (T_m x d), residuals
      all_rss[m] <- sum(diag((crossprod(U_m, U_m)))) # The sum of squares of residuals
    }

    # Obtain the estimates in the vector (\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M)
    estims <- c(all_intercepts, all_AR_mats)

    # Return the estimates and the sum of squares of residuals, the last element is the for rss
    c(estims, sum(all_rss))
  }

  ## Least squares estimation function given thresholds for models with AR constraints
  LS_with_AR_constraints <- function(thresholds) {
    # threshold = length M-1 vector of the thresholds r_1,...,r_{M-1}; if M=1 anything is ok
    # Other arguments are taken from the parent environment.

    # Create the constraint matrix \tilde{C} with dummy constraints for the intercepts as well
    q <- ncol(AR_constraints) # The number of intercept and AR parameters under constraints
    C_tilde <- rbind(cbind(diag(rep(1, times=M*d)), matrix(0, nrow=M*d, ncol=q)),
                     cbind(matrix(0, nrow=M*p*d^2, ncol=M*d), AR_constraints))

    # Storage for the estimates
    estims <- array(NA, dim=c(d, d*p + 1, M)) # [, , m] for [\phi_{m,0} : A_{m,1} : ... : A_{m,p}]

    # In Y, i:th row denotes the vector \bold{y_{i-1}} = (y_{i-1},...,y_{i-p}) (dpx1),
    # assuming the observed data is y_{-p+1},...,y_0,y_1,...,y_{T}. The last row is for
    # the vector (y_{T},...,y_{T-p}).
    Y <- reform_data(data, p) # (T_obs + 1 x dp)
    Y2 <- Y[1:T_obs, , drop=FALSE] # Last row removed; not needed when only lagged observations used

    # The transition weights: (T x M) matrix, the t:th row is for the time point t and the m:th column is for the regime m.
    alpha_mt <- get_alpha_mt(data, Y2=Y2, p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                             weightpars=thresholds) # Transition weights (T x M), t:th row

    # Now the observations are stored to a (Td x 1) vector defining the matrix Y.
    Y <- as.matrix(vec(t(data[(p+1):nrow(data),]))) # (Td x 1), we don't take the first p observations (init values)

    # Initialize the matrix of regressors X.
    X <- matrix(0, nrow=d*T_obs, ncol=M*d + M*p*d^2) # (Td x Md + Mpd^2)

    # We construct the (Td x Md + Mpd^2) matrix X consisting of (d x Md + Mpd^2) blocks X_t for each time period t.
    # Each block X_t can be partioned to the intercept part X_t_phi (d x Md) and the AR part X_t_A (d x Mpd^2).
    I_d <- diag(rep(1, times=d)) # The (d x d) identity matrix
    for(t in 1:T_obs) { # We need to loop through all time periods
      m <- which(abs(alpha_mt[t, ] - 1) < 1e-6) # The regime m that prevails at time t
      X_t_phi <- cbind(matrix(0, nrow=d, ncol=(m - 1)*d), I_d, matrix(0, nrow=d, ncol=(M - m)*d)) # (d x Md)
      X_t_A <- matrix(0, nrow=d, ncol=M*p*d^2) # (d x Mpd^2), initialize the AR part
      for(j in 1:p) { # Go through the lags and fill in the blocks
        # The Block_{A_{m,j}} is the columns: (dp*(m - 1)+ d*(j - 1) + 1):(dp*(m - 1) + d*j)
        X_t_A[, ((m - 1)*p*d^2 + (j - 1)*d^2 + 1):((m - 1)*p*d^2 + j*d^2)] <- kronecker(t(data[t + p - j,]), I_d) # (d x d^2) block
      }
      # Will in the X_t block to the X matrix
      X[((t - 1)*d + 1):((t - 1)*d + d), ] <- cbind(X_t_phi, X_t_A) # (d x Md + Mpd^2)
    }

    # Integrate the constraints into the design matrix
    X_bold <- X%*%C_tilde # (Td x q)

    estims <- tryCatch(solve(crossprod(X_bold, X_bold), crossprod(X_bold, Y)), # (M*d + q x 1)
                       error=function(e) matrix(0, nrow=M*d + q, ncol=1)) # zero estimates are legal but bad, dummy estimates

    # Calculate the sum of squares of residuals
    U <- Y - X_bold%*%estims # (Td x 1), residuals
    rss <- crossprod(U, U) # The sum of squares of residuals

    # Return the estimates and the sum of squares of residuals, the last element is the for rss
    c(estims, rss)
  }

  ## A function to check whether the stability condition is satisfied for the AR matrices, and
  ## if not, to what extend it is not satisfied.
  stab_exceeded <- function(estims) {
    # Estims should be a vector of the form (\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M)
    if(!is.null(AR_constraints)) { # Expand the AR constraints
      pars_to_check <- c(estims[1:(M*d)], AR_constraints%*%estims[(M*d + 1):(M*d + ncol(AR_constraints))])
    } else {
      pars_to_check <- estims[1:(M*d + M*p*d^2)]
    }
    all_phi0 <- pick_phi0(M=M, d=d, params=pars_to_check)
    all_A <- pick_allA(p=p, M=M, d=d, params=pars_to_check)
    all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
    all_stab_exceeds <- matrix(nrow=nrow(all_boldA[, , 1]), ncol=M) # Square of how much modulus of eigenvalues exceed 1 - stab_tol
    for(m in 1:M) { # Check stability condition for each regime
      abs_eigs <- abs(eigen(all_boldA[, , m], symmetric=FALSE, only.values=TRUE)$values)
      all_stab_exceeds[, m] <- pmax(0, abs_eigs - (1 - stab_tol))^2 # How much abs eigens exceed 1 - stab_tol, squared
    }
    sum(all_stab_exceeds) # Sum of the squared exceeded values of stab cond
  }

  ###########################################
  ## Create threshold vectors and estimate ##
  ###########################################

  ## Create the set of threshold vectors for the optimization; M=1 will use numeric(0) and run the LS only once
  if(is.null(weight_constraints)) {
    switch_var_series <- data[,weightfun_pars[1]] # The switching variable time series
    sv_sorted_full <- sort(switch_var_series, decreasing=FALSE) # The sorted switch variable series

    # Remove the values from the sorted switch variable series that would leave less than T_min observations
    # smaller or larger than the threshold.
    switch_var_sorted <- sv_sorted_full[T_min:(length(switch_var_series) - T_min)]
    min_switchvar <- min(switch_var_sorted)
    max_switchvar <- max(switch_var_sorted)

    # The maximum number of grid points is calculated so that the number M-1 dimensional of multisets
    # is at most 20000 for M <= 4.
    if(M >= 2) { # M=1 case is separately handled
      max_thresholds <- 1000 # The maximum number of threshold values to considered
      if(M == 2) {
        if(length(switch_var_sorted) < max_thresholds) {
          grid_points <- switch_var_sorted
        } else {
          grid_points <- seq(from=min_switchvar, to=max_switchvar, length.out=max_thresholds)
        }
      } else if(M == 3) {
        grid_points <- seq(from=min_switchvar, to=max_switchvar, length.out=min(200, max_thresholds))
      } else if(M == 4) {
        grid_points <- seq(from=min_switchvar, to=max_switchvar, length.out=50)
      } else {
        grid_points <- seq(from=min_switchvar, to=max_switchvar, length.out=30)
      }
      thresholds <- t(utils::combn(x=grid_points, m=M - 1, simplify=TRUE)) # M-1 dim multisets of lexically ordered grid points
    }
    if(M == 2) {
      thresvecs <- thresholds # Each row for each threshold vector (scalar in this case)
    } else if(M > 2) {
      obs_between_thresholds <-  matrix(NA, nrow=nrow(thresholds), ncol=M - 2) # The number of observations between the thresholds
      for(m in 1:(M - 2)) {
        n_at_most_upper <- findInterval(x=thresholds[, m + 1], vec=sv_sorted_full, left.open=TRUE, rightmost.closed=TRUE)
        n_at_most_lower <- findInterval(x=thresholds[, m], vec=sv_sorted_full, left.open=TRUE, rightmost.closed=TRUE)
        obs_between_thresholds[,m] <- n_at_most_upper - n_at_most_lower # Number of observations between lower and upper threshold
      }
      # The threshold vectors with enough observations in all regimes
      #print(which(rowSums(obs_between_thresholds >= T_min) == ncol(obs_between_thresholds)))
      thresvecs <- thresholds[which(rowSums(obs_between_thresholds >= T_min) == ncol(obs_between_thresholds)), , drop=FALSE]
    }
  } else { # thresholds fixed to known numbers
    thresvecs <- matrix(weight_constraints[[2]], nrow=1, ncol=M - 1)
  }
  # Each row in thresvecs of a threshold vector

  ## Estimate the model for all thresholds vectors in thresvecs
  estim_fun <- if(is.null(AR_constraints)) LS_without_AR_constraints else LS_with_AR_constraints
  estim_length <- if(is.null(AR_constraints)) M*d + M*p*d^2 + 1 else M*d + ncol(AR_constraints) + 1

  if(M == 1) {
    if(use_parallel) message(paste("PHASE 1: Estimating the AR and weight parameters by least squares..."))
    estims <- as.matrix(estim_fun(numeric(0)))
    all_stab_ex <- stab_exceeded(estims[,1])
  } else {
    if(use_parallel) {
      if(ncores > parallel::detectCores()) {
        ncores <- parallel::detectCores()
      }
      n_thresvecs <- ifelse(M == 1 || !is.null(weight_constraints), 1, nrow(thresvecs))
      message(paste0("PHASE 1: Estimating the AR and weight parameters by least squares for ", n_thresvecs,
                     " vectors of thresholds...")) # "PHASE 1" print i related to the multiple-phase estimation procedure
      cl <- parallel::makeCluster(ncores)
      on.exit(try(parallel::stopCluster(cl), silent=TRUE)) # Close the cluster on exit, if not already closed.
      parallel::clusterExport(cl, ls(environment(estim_LS)), envir=environment(estim_LS)) # assign all variables from package:sstvars
      parallel::clusterEvalQ(cl, c(library(pbapply), library(sstvars)))
      estims <- as.matrix(simplify2array(pbapply::pblapply(1:nrow(thresvecs), FUN=function(i1) estim_fun(thresvecs[i1,]), cl=cl)))

      if(penalized) {
        if(M > 2) {
          message(paste0("Checking the stability condition for all the LS estimates..."))
          all_stab_ex <- simplify2array(pbapply::pblapply(1:nrow(thresvecs), FUN=function(i1) stab_exceeded(estims[,i1]), cl=cl))
        } else { # Less prints, since the calculations are fast enough
          all_stab_ex <- simplify2array(pbapply::pblapply(1:nrow(thresvecs), FUN=function(i1) stab_exceeded(estims[,i1]), cl=cl))
        }
      }
      parallel::stopCluster(cl=cl)
    } else { # No parallel computing
      estims <- as.matrix(vapply(1:nrow(thresvecs), FUN=function(i1) estim_fun(thresvecs[i1,]),
                                 FUN.VALUE=numeric(estim_length)))

      if(penalized) {
        all_stab_ex <- vapply(1:nrow(thresvecs), FUN=function(i1) stab_exceeded(estims[,i1]), FUN.VALUE=numeric(1))
      }
    }
  }
  # Each column in estims corresponds to each vector of thresholds

  ## Obtain the LS estimates, possibly among stable estimates
  if(penalized) {
    # Determine the tuning parameter value that controls the extend of the penalization
    all_rss <- estims[nrow(estims),]
    min_rss <- min(all_rss) # The smallest residual sum of squares
    penalty_coef <- tuning_par*min_rss # The penalty coefficient

    # Obtain the index with the smallest penalized sum of squares
    penalized_stab_ex <- penalty_coef*all_stab_ex # Penalization for non-stable estimates
    all_pen_rss <- all_rss + penalized_stab_ex # Penalized sum of squares of residuals
    min_rss_index <- which.min(all_pen_rss)[1] # The index for which the penalized sum of squares is the smallest

  } else {
    # Find the index for which the sum of squares of residuals is the smallest (regardless of stability)
    min_rss_index <- which.min(estims[nrow(estims),])[1]
  }


  ## Obtain and return the estimates corresponding the smallest sum of squares of residuals
  int_and_ar_estims <- estims[1:(nrow(estims) - 1), min_rss_index]
  if(M == 1 || !is.null(weight_constraints)) {
    threshold_estims <- numeric(0) # No threshold estimates
  } else {
    threshold_estims <- thresvecs[min_rss_index,]
  }
  c(int_and_ar_estims, threshold_estims) # Return the estimates
}


#' @title Internal estimation function for estimating autoregressive and weight parameters of
#'   STVAR models by the method of nonlinear least squares.
#'
#' @description \code{estim_NLS} estimates the autoregressive and weight parameters of STVAR models
#'   by the method of least squares (\code{relative_dens} weight function is not supported).
#'
#' @inheritParams estim_LS
#' @details Used internally in the multiple phase estimation procedure proposed by Virolainen (2025).
#'  The weight function \code{relative_dens} is not supported.  Mean constraints are not supported.
#'  Only weight constraints that specify the weight parameters as fixed values are supported.
#'  Only intercept parametrization is supported.
#' @return Returns the estimated parameters in a vector of the form
#'  \eqn{(\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M,\alpha}, where
#'  \itemize{
#'    \item{\eqn{\phi_{m,0} = } the \eqn{(d \times 1)} intercept vector of the \eqn{m}th regime.}
#'    \item{\eqn{\varphi_m = (vec(A_{m,1}),...,vec(A_{m,p}))} \eqn{(pd^2 \times 1)}.}
#'    \item{\eqn{\alpha}} is the vector of the weight parameters: \describe{
#'      \item{\code{weight_function="relative_dens"}:}{\eqn{\alpha = (\alpha_1,...,\alpha_{M-1})}
#'        \eqn{(M - 1 \times 1)}, where \eqn{\alpha_m} \eqn{(1\times 1), m=1,...,M-1} are the transition weight parameters.}
#'      \item{\code{weight_function="logistic"}:}{\eqn{\alpha = (c,\gamma)}
#'        \eqn{(2 \times 1)}, where \eqn{c\in\mathbb{R}} is the location parameter and \eqn{\gamma >0} is the scale parameter.}
#'      \item{\code{weight_function="mlogit"}:}{\eqn{\alpha = (\gamma_1,...,\gamma_M)} \eqn{((M-1)k\times 1)},
#'        where \eqn{\gamma_m} \eqn{(k\times 1)}, \eqn{m=1,...,M-1} contains the multinomial logit-regression coefficients
#'        of the \eqn{m}th regime. Specifically, for switching variables with indices in \eqn{I\subset\lbrace 1,...,d\rbrace}, and with
#'        \eqn{\tilde{p}\in\lbrace 1,...,p\rbrace} lags included, \eqn{\gamma_m} contains the coefficients for the vector
#'        \eqn{z_{t-1} = (1,\tilde{z}_{\min\lbrace I\rbrace},...,\tilde{z}_{\max\lbrace I\rbrace})}, where
#'        \eqn{\tilde{z}_{i} =(y_{it-1},...,y_{it-\tilde{p}})}, \eqn{i\in I}. So \eqn{k=1+|I|\tilde{p}}
#'        where \eqn{|I|} denotes the number of elements in \eqn{I}.}
#'      \item{\code{weight_function="exponential"}:}{\eqn{\alpha = (c,\gamma)}
#'        \eqn{(2 \times 1)}, where \eqn{c\in\mathbb{R}} is the location parameter and \eqn{\gamma >0} is the scale parameter.}
#'      \item{\code{weight_function="threshold"}:}{\eqn{\alpha = (r_1,...,r_{M-1})}
#'         \eqn{(M-1 \times 1)}, where \eqn{r_1,...,r_{M-1}} are the thresholds.}
#'      \item{\code{weight_function="exogenous"}:}{Omit \eqn{\alpha} from the parameter vector.}
#'    }
#'  }
#'  For models with...
#'   \describe{
#'     \item{AR_constraints:}{Replace \eqn{\varphi_1,...,\varphi_M} with \eqn{\psi} as described in the argument \code{AR_constraints}.}
#'     \item{weight_constraints:}{If linear constraints are imposed, replace \eqn{\alpha} with \eqn{\xi} as described in the
#'      argument \code{weigh_constraints}. If weight functions parameters are imposed to be fixed values, simply drop \eqn{\alpha}
#'      from the parameter vector.}
#'   }
#' @references
#'  \itemize{
#'    \item Hubrich K., Teräsvirta. T. 2013. Thresholds and Smooth Transitions in Vector Autoregressive Models.
#'      \emph{CREATES Research Paper 2013-18, Aarhus University.}
#'    \item Virolainen S. 2024. Identification by non-Gaussianity in structural threshold and
#'       smooth transition vector autoregressive models. Unpublished working
#'       paper, available as arXiv:2404.19707.
#'    \item Virolainen S. 2025. Unpublished working paper.
#'  }
#' @keywords internal

estim_NLS <- function(data, p, M, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                      weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t"),
                      parametrization=c("intercept", "mean"), AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL,
                      penalized=TRUE, penalty_params=c(0.05, 0.2), use_parallel=TRUE, ncores=2) {
  # Checks
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  if(!is.null(mean_constraints)) {
    stop("Mean constraints are not supported by the NLS estimation")
  } else if(weight_function == "relative_dens") {
    stop("The relative_dens weight function is not supported by the NLS estimation")
  } else if(parametrization != "intercept") {
    stop("Only the intercept parametrization is supported by the LS estimation")
  }
  # Check the weight constraints
  if(!is.null(weight_constraints)) {
    if(!all(weight_constraints[[1]] == 0)) {
      stop(paste("Only such weight_constraints that specify the threshold parameters some known fixed values",
                 "are supported in the least squares estimation."))
    }
  }
  data <- check_data(data, p=p)
  d <- ncol(data)
  weightfun_pars <- check_weightfun_pars(data=data, p=p, d=d, M=M, weight_function=weight_function,
                                         weightfun_pars=weightfun_pars, cond_dist=cond_dist)
  check_constraints(data=data, p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    parametrization=parametrization, identification="reduced_form",
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=NULL)

  stopifnot(is.numeric(penalty_params) && length(penalty_params) == 2 && all(penalty_params >= 0) && penalty_params[1] < 1)
  stab_tol <- penalty_params[1]
  tuning_par <- penalty_params[2]

  # Obtain relevant statistics
  n_obs <- nrow(data)
  T_obs <- n_obs - p
  if(weight_function != "exogenous") {
    T_min <- 2/d*ifelse(is.null(AR_constraints), (p + 1)*d^2 + d,
                        ncol(AR_constraints)/M + d + d^2) # Minimum number of obs in each regime
    if(T_obs/M < T_min) { # Try smaller T_min
      T_min <- 1.5/d*ifelse(is.null(AR_constraints), (p + 1)*d^2 + d,
                            ncol(AR_constraints)/M + d + d^2)
      message("The number of observations in relative to the number of parameters is small, consider decreasing p or M.")
      if(T_obs/M < T_min) {
        stop("The number of observations is too small for reasonable estimation. Decrease the order p or the number of regimes M.")
      }
    }
  }

  ################################
  ## Create estim etc functions ##
  ################################

  # Logarithm of the smallest value that can be handled normally (used in get_alpha_mt)
  epsilon <- round(log(.Machine$double.xmin) + 10)

  ## Nonlinear least squares estimation function given weight parameters
  NLS_est <- function(weightpars, AR_constraints) {
    # weightpars = vector of the weight parameters, AR_constraints = as usual
    # Other arguments are taken from the parent env

    if(!is.null(AR_constraints)) {
      # Create the matrix C_tilde that sets dummy constraints for the intercepts:
      q <- ncol(AR_constraints) # The number of intercept and AR parameters under constraints
      C_tilde <- rbind(cbind(diag(rep(1, times=M*d)), matrix(0, nrow=M*d, ncol=q)),
                       cbind(matrix(0, nrow=M*p*d^2, ncol=M*d), AR_constraints))

      # Create the permutation matrix that expand \tilde{C}\tilde{\psi} to \beta for AR constraints:
      P <- matrix(0, nrow=M*d + M*p*d^2, ncol=M*d + M*p*d^2)
      I_d <- diag(nrow=d)
      I_pd2 <- diag(nrow=p*d^2)
      for(m in 1:M) { # Go through the regimes
        # Insert the identity matrix for the intercepts
        P[((m - 1)*d + (m - 1)*p*d^2 + 1):(m*d + (m - 1)*p*d^2), ((m - 1)*d + 1):(m*d)] <- I_d

        # Insert the identity matrix for the AR matrices
        P[(m*d + (m - 1)*p*d^2 + 1):(m*d + m*p*d^2), (M*d + (m - 1)*p*d^2 + 1):(M*d + m*p*d^2)] <- I_pd2
      }

      # Calculate the multiplication of P and C_tilde that expands psi_tilde to the usual parameter vector:
      PC <- P%*%C_tilde
    }

    # In Y, i:th row denotes the vector \bold{y_{i-1}} = (y_{i-1},...,y_{i-p}) (dpx1),
    # assuming the observed data is y_{-p+1},...,y_0,y_1,...,y_{T}. The last row is for
    # the vector (y_{T},...,y_{T-p}).
    Y <- reform_data(data, p) # (T_obs + 1 x dp)
    Y2 <- Y[1:T_obs, , drop=FALSE] # Last row removed; not needed when only lagged observations used

    # The transition weights: (T x M) matrix, the t:th row is for the time point t and the m:th column is for the regime m.
    alpha_mt <- get_alpha_mt(data, Y2=Y2, p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                             weightpars=weightpars, epsilon=epsilon) # Transition weights (T x M), t:th row

    # Storage for the data matrices \Psi_t in the form y_t = \Psi_t\beta + u_t,
    # \beta = (vec(\Phi_1),...,vec(\Phi_M)), \Phi_m = [\phi_{m,0} : A_{m,1} : ... : A_{m,p}].
    all_Psi <- array(NA, dim=c(d, M*d + M*p*d^2, T_obs)) # [, , t] for \Psi_t
    all_cPsi <- array(NA, dim=c(M*d + M*p*d^2, M*d + M*p*d^2, T_obs)) # [, , t] for crossprod(\Psi_t)
    all_tPsiy <- array(NA, dim=c(M*d + M*p*d^2, 1, T_obs)) # [, , t] for \Psi_t'y_t

    # Storages for models with AR constraints
    if(!is.null(AR_constraints)) {
      all_PsiPC <- array(NA, dim=c(d, ncol(C_tilde), T_obs)) # [, , t] for \Psi_tPC
      all_cPsiPC <- array(NA, dim=c(ncol(C_tilde), ncol(C_tilde), T_obs)) # [, , t] for crossprod(\Psi_tPC)
      all_tPsiPCy <- array(NA, dim=c(ncol(C_tilde), 1, T_obs)) # [,  ,t] for (\Psi_tPC)'y_t
    }

    # Go through the time periods
    for(t in 1:T_obs) {

      # Create the vector Z_t = (1, y_{t-1},...,y_{t-p}) (1 + dp x 1), and calculate Knonecker product between Z_t' and I_d
      Z_kron <- kronecker(matrix(c(1, Y[t,]), nrow=1), diag(d)) # (d x d(1 + dp)), the Kronecker product

      # Go through the regimes
      for(m in 1:M) {
        all_Psi[, ((m - 1)*(d + p*d^2) + 1):(m*(d + p*d^2)), t] <- alpha_mt[t, m]*Z_kron # Fill in Psi_t, also for AR_constraints
      }

      if(!is.null(AR_constraints)) {
        all_PsiPC[, , t] <- all_Psi[, , t]%*%PC # Multiply Psi_t with PC
        all_cPsiPC[, , t] <- crossprod(all_PsiPC[, , t]) # Fill in Psi_t'Psi_t
        all_tPsiPCy[, , t] <- crossprod(all_PsiPC[, , t], data[p + t,]) # Fill in \Psi_t'y_t
      } else {
        all_cPsi[, , t] <- crossprod(all_Psi[, , t]) # Fill in Psi_t'Psi_t
        all_tPsiy[, , t] <- crossprod(all_Psi[, , t], data[p + t,]) # Fill in \Psi_t'y_t
      }
    }

    # Sum over the time periods in cPsi and tPsiy
    if(!is.null(AR_constraints)) {
      sum_cPsi <- apply(all_cPsiPC, MARGIN=c(1, 2), FUN=sum)
      sum_tPsiy <- apply(all_tPsiPCy, MARGIN=c(1, 2), FUN=sum)
    } else {
      sum_cPsi <- apply(all_cPsi, MARGIN=c(1, 2), FUN=sum)
      sum_tPsiy <- apply(all_tPsiy, MARGIN=c(1, 2), FUN=sum)
    }
    ## Obtain the estimates in the vector \beta = (vec(\Phi_1),...,vec(\Phi_M)), where \Phi_m = [\phi_{m,0} : A_{m,1} : ... : A_{m,p}]
    # (with AR_constraints \beta = (\phi_{1,0},...,\phi_{M,0},\psi)
    #estims <- solve(sum_cPsi, sum_tPsiy) # (M*d + M*p*d^2 x 1)
    estims <- tryCatch(solve(sum_cPsi, sum_tPsiy), # (M*d + M*p*d^2 x 1), fails if the system is singular
                       error=function(e) matrix(0, nrow=M*d + M*p*d^2, ncol=1)) # zero estimates are legal but bad, dummy estimates

    ## Calculate the residual sums of squares
    if(is.null(AR_constraints)) {
      est_to_use <- estims
    } else {
      est_to_use <- PC%*%estims # Estimates in the standard form
    }
    all_rss <- numeric(T_obs) # Storage for all residual sums of squares
    for(t in 1:T_obs) { # Go through the time periods again
      u_t <- data[p + t,] - all_Psi[, , t]%*%est_to_use # The residual for time period t
      all_rss[t] <- crossprod(u_t, u_t) # The residual sum of squares for time period t
    }

    ## Obtain the estimates in the vector (\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M)
    # (for models with AR constraints (\phi_{1,0},...,\phi_{M,0},\psi), already in the correct form):
    if(is.null(AR_constraints)) {
      estims <- matrix(estims, nrow=d) # estims in the form [Phi_1,...,Phi_M]
      int_cols <- 0:(M - 1)*(d*p + 1) + 1 # which columns have the intercept parameters
      estims <- c(estims[, int_cols], estims[, -int_cols]) # Estimates in the form (\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M)
    }

    ## Return the estimates and the sum of squares of residuals, the last element is the for rss
    c(estims, sum(all_rss))
  }

  ## A function to check whether the stability condition is satisfied for the AR matrices, and
  ## if not, to what extend it is not satisfied.
  stab_exceeded <- function(estims) {
    # Estims should be a vector of the form (\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M)
    if(!is.null(AR_constraints)) { # Expand the AR constraints
      pars_to_check <- c(estims[1:(M*d)], AR_constraints%*%estims[(M*d + 1):(M*d + ncol(AR_constraints))])
    } else {
      pars_to_check <- estims[1:(M*d + M*p*d^2)]
    }
    all_phi0 <- pick_phi0(M=M, d=d, params=pars_to_check)
    all_A <- pick_allA(p=p, M=M, d=d, params=pars_to_check)
    all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
    all_stab_exceeds <- matrix(nrow=nrow(all_boldA[, , 1]), ncol=M) # Square of how much modulus of eigenvalues exceed 1 - stab_tol
    for(m in 1:M) { # Check stability condition for each regime
      abs_eigs <- abs(eigen(all_boldA[, , m], symmetric=FALSE, only.values=TRUE)$values)
      all_stab_exceeds[, m] <- pmax(0, abs_eigs - (1 - stab_tol))^2 # How much abs eigens exceed 1 - stab_tol, squared
    }
    sum(all_stab_exceeds) # Sum of the squared exceeded values of stab cond
  }

  ############################################
  ## Create weight par vectors and estimate ##
  ############################################

  ## Create the set of weight parameters for the optimization; M=1 will use numeric(0) and run the NLS only once
  if(is.null(weight_constraints) && weight_function != "exogenous") {
    if(weight_function != "mlogit") {
      switch_var_series <- data[,weightfun_pars[1]] # The switching variable time series
      sv_sorted_full <- sort(switch_var_series, decreasing=FALSE) # The sorted switch variable series
      switch_var_sorted <- sv_sorted_full[T_min:(length(switch_var_series) - T_min)] # Switch var vals >=T_min vals below and above
      min_switchvar <- min(switch_var_sorted)
      max_switchvar <- max(switch_var_sorted)
    }
    if(weight_function %in% c("logistic", "exponential")) {
      # Here always M=2, the first weight parameter is location parameter, and the second one is strictly positive scale parameter
      c_grid <- seq(from=min_switchvar, to=max_switchvar, length.out=100) # The grid for the location parameter
      gamma_grid <- seq(from=0.1, to=100, length.out=100) # The grid for the scale parameter
      weightparvecs <- unname(simplify2array(expand.grid(c_grid, gamma_grid))) # Each row for a vector of weight parameters
    } else if(weight_function == "mlogit") {
      n_weight_pars <- (M - 1)*(1 + length(weightfun_pars[[1]])*weightfun_pars[[2]])
      weightparvecs <- unname(simplify2array(expand.grid(replicate(n=n_weight_pars,
                                                                   expr=seq(from=-20, to=20,
                                                                            length.out=max(2, ceiling(10000^(1/n_weight_pars)))),
                                                                   simplify=FALSE)))) # Roughly 10000-100000 grid points
    } else if(weight_function == "threshold") {
      # Remove the values from the sorted switch variable series that would leave less than T_min observations
      # smaller or larger than the threshold.
      switch_var_sorted <- sv_sorted_full[T_min:(length(switch_var_series) - T_min)]
      min_switchvar <- min(switch_var_sorted)
      max_switchvar <- max(switch_var_sorted)

      # The maximum number of grid points is calculated so that the number M-1 dimensional of multisets
      # is at most 20000 for M <= 4.
      if(M >= 2) { # M=1 case is separately handled
        max_thresholds <- 1000 # The maximum number of threshold values to considered
        if(M == 2) {
          if(length(switch_var_sorted) < max_thresholds) {
            grid_points <- switch_var_sorted
          } else {
            grid_points <- seq(from=min_switchvar, to=max_switchvar, length.out=max_thresholds)
          }
        } else if(M == 3) {
          grid_points <- seq(from=min_switchvar, to=max_switchvar, length.out=min(200, max_thresholds))
        } else if(M == 4) {
          grid_points <- seq(from=min_switchvar, to=max_switchvar, length.out=50)
        } else {
          grid_points <- seq(from=min_switchvar, to=max_switchvar, length.out=30)
        }
        thresholds <- t(utils::combn(x=grid_points, m=M - 1, simplify=TRUE)) # M-1 dim multisets of lexically ordered grid points
      }
      if(M == 2) {
        weightparvecs <- thresholds # Each row for each threshold vector (scalar in this case)
      } else if(M > 2) {
        obs_between_thresholds <-  matrix(NA, nrow=nrow(thresholds), ncol=M - 2) # The number of observations between the thresholds
        for(m in 1:(M - 2)) {
          n_at_most_upper <- findInterval(x=thresholds[, m + 1], vec=sv_sorted_full, left.open=TRUE, rightmost.closed=TRUE)
          n_at_most_lower <- findInterval(x=thresholds[, m], vec=sv_sorted_full, left.open=TRUE, rightmost.closed=TRUE)
          obs_between_thresholds[,m] <- n_at_most_upper - n_at_most_lower # Number of observations between lower and upper threshold
        }
        # The threshold vectors with enough observations in all regimes
        weightparvecs <- thresholds[which(rowSums(obs_between_thresholds >= T_min) == ncol(obs_between_thresholds)), , drop=FALSE]
      }
    }

  } else if(weight_function == "exogenous") { # Exogenous weights
    weightparvecs <- weightfun_pars
  } else { # weight parameters fixed to known numbers
    weightparvecs <- matrix(weight_constraints[[2]], nrow=1)
  }
  # Each row in weightparvecs correspond to one vector of weight parameters

  ## Estimate the model for all weight pars in weight_pars
  estim_length <- if(is.null(AR_constraints)) M*d + M*p*d^2 + 1 else M*d + ncol(AR_constraints) + 1

  if(M == 1) {
    if(use_parallel) message(paste("PHASE 1: Estimating the AR and weight parameters by nonlinear least squares..."))
    estims <- as.matrix(NLS_est(numeric(0), AR_constraints=AR_constraints))
    all_stab_ex <- stab_exceeded(estims[,1])
  } else {
    if(use_parallel) {
      if(ncores > parallel::detectCores()) {
        ncores <- parallel::detectCores()
      }
      n_weightvecs <- ifelse(M == 1 || !is.null(weight_constraints), 1, nrow(weightparvecs))
      message(paste0("PHASE 1: Estimating the AR and weight parameters by least squares for ", n_weightvecs,
                     " vectors of thresholds...")) # "PHASE 1" print i related to the multiple-phase estimation procedure
      cl <- parallel::makeCluster(ncores)
      on.exit(try(parallel::stopCluster(cl), silent=TRUE)) # Close the cluster on exit, if not already closed.
      parallel::clusterExport(cl, ls(environment(estim_LS)), envir=environment(estim_LS)) # assign all variables from package:sstvars
      parallel::clusterEvalQ(cl, c(library(pbapply), library(sstvars)))
      estims <- as.matrix(simplify2array(pbapply::pblapply(1:nrow(weightparvecs),
                                                           FUN=function(i1) NLS_est(weightparvecs[i1,],
                                                                                    AR_constraints=AR_constraints), cl=cl)))

      if(penalized) {
        if(M > 2) {
          message(paste0("Checking the stability condition for all the LS estimates..."))
          all_stab_ex <- simplify2array(pbapply::pblapply(1:nrow(weightparvecs), FUN=function(i1) stab_exceeded(estims[,i1]), cl=cl))
        } else { # Less prints, since the calculations are fast enough
          all_stab_ex <- simplify2array(pbapply::pblapply(1:nrow(weightparvecs), FUN=function(i1) stab_exceeded(estims[,i1]), cl=cl))
        }
      }
      parallel::stopCluster(cl=cl)
    } else { # No parallel computing
      estims <- as.matrix(vapply(1:nrow(weightparvecs), FUN=function(i1) NLS_est(weightparvecs[i1,], AR_constraints=AR_constraints),
                                 FUN.VALUE=numeric(estim_length)))

      if(penalized) {
        all_stab_ex <- vapply(1:nrow(weightparvecs), FUN=function(i1) stab_exceeded(estims[,i1]), FUN.VALUE=numeric(1))
      }
    }
  }
  # Each column in estims corresponds to each vector of thresholds

  ## Obtain the LS estimates, possibly among stable estimates
  if(penalized) {
    # Determine the tuning parameter value that controls the extend of the penalization
    all_rss <- estims[nrow(estims),]
    min_rss <- min(all_rss) # The smallest residual sum of squares
    penalty_coef <- tuning_par*min_rss # The penalty coefficient

    # Obtain the index with the smallest penalized sum of squares
    penalized_stab_ex <- penalty_coef*all_stab_ex # Penalization for non-stable estimates
    all_pen_rss <- all_rss + penalized_stab_ex # Penalized sum of squares of residuals
    min_rss_index <- which.min(all_pen_rss)[1] # The index for which the penalized sum of squares is the smallest

  } else {
    # Find the index for which the sum of squares of residuals is the smallest (regardless of stability)
    min_rss_index <- which.min(estims[nrow(estims),])[1]
  }


  ## Obtain and return the estimates corresponding the smallest sum of squares of residuals
  int_and_ar_estims <- estims[1:(nrow(estims) - 1), min_rss_index]
  if(M == 1 || !is.null(weight_constraints) || weight_function == "exogenous") {
    weightpar_estims <- numeric(0) # No weightpar estimates
  } else {
    weightpar_estims <- weightparvecs[min_rss_index,]
  }
  c(int_and_ar_estims, weightpar_estims) # Return the estimates
}
