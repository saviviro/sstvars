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

iterate_more <- function(stvar, maxit=1000, h=1e-3, penalized, penalty_params, allow_unstab, calc_std_errors=TRUE, print_trace=TRUE) {
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
  ret$LS_estimates <- stvar$LS_estimates
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

