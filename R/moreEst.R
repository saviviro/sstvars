#' @title Maximum likelihood estimation of a reduced form or structural STVAR model based on preliminary estimates
#'
#' @description \code{iterate_more} uses a variable metric algorithm to estimate a reduced form or structural STVAR model
#'  (object of class \code{'stvar'}) based on preliminary estimates.
#'
#' @inheritParams fitSTVAR
#' @inheritParams STVAR
#' @param stvar an object of class \code{'stvar'}, created by, e.g., \code{fitSTVAR}.
#' @details The purpose of \code{iterate_more} is to provide a simple and convenient tool to finalize
#'   the estimation when the maximum number of iterations is reached when estimating a STVAR model
#'   with the main estimation function \code{fitSTVAR} or \code{fitSSTVAR}.
#' @return Returns an object of class \code{'stvar'} defining the estimated model
#' @seealso \code{\link{fitSTVAR}}, \code{\link{STVAR}}, \code{\link[stats]{optim}}
#' @inherit STVAR references
#' @examples
#' \donttest{
#' ## These are long running examples that use parallel computing!
#' ## Running the below examples takes approximately FILL IN HOW MANY minutes
#'
#' # p=1, M=2, d=2 relative_dens Gaussian STVAR, only 5 iterations of
#' # the variable matrix algorithm.
#' fit12 <- fitSTVAR(gdpdef, p=1, M=2, nrounds=1, seeds=1, maxit=5, use_parallel=FALSE)
#' fit12
#'
#' # Iterate more:
#' fit12 <- iterate_more(fit12)
#' fit12
#' }
#' @export

iterate_more <- function(stvar, maxit=100, calc_std_errors=TRUE) {
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
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function,
                                         weightfun_pars=stvar$model$weightfun_pars, cond_dist=cond_dist)
  minval <- get_minval(stvar$data)
  npars <- length(stvar$params)

  # Function to optimize
  loglik_fn <- function(params) {
    tryCatch(loglikelihood(data=data, p=p, M=M, params=params,
                           weight_function=weight_function, weightfun_pars=weightfun_pars,
                           cond_dist=cond_dist, parametrization=parametrization,
                           identification=identification, AR_constraints=AR_constraints,
                           mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                           B_constraints=B_constraints, to_return="loglik", check_params=TRUE, minval=minval),
             error=function(e) minval)
  }

  # Gradient of the log-likelihood function using central difference approximation
  h <- 6e-6
  I <- diag(rep(1, times=npars))
  loglik_grad <- function(params) {
    vapply(1:npars, function(i1) (loglik_fn(params + I[i1,]*h) - loglik_fn(params - I[i1,]*h))/(2*h), numeric(1))
  }

  res <- optim(par=stvar$params, fn=loglik_fn, gr=loglik_grad, method=c("BFGS"), control=list(fnscale=-1, maxit=maxit))
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
  warn_eigens(ret)
  ret
}



#' @title Maximum likelihood estimation of a structural STVAR model based on preliminary estimates from
#'   a reduced form model.
#'
#' @description \code{fitSSTVAR} uses a variable metric algorithm to estimate a structural STVAR model
#'   based on preliminary estimates from a reduced form model.
#'
#' @inheritParams fitSTVAR
#' @inheritParams STVAR
#' @param stvar a an object of class \code{'stvar'}, created by, e.g., \code{fitSTVAR},
#'   specifying a reduced form or a structural model
#' @param robust_method Should some robust estimation method be used in the estimation before switching
#'   to the gradient based variable metric algorithm? See details.
#'
#' @details When the structural model does not impose overidentifying constraints, it is directly
#'   obtained from the reduced form model, and estimation is not required. When overidentifying constraints
#'   are imposed, the model is estimated via ..
#'
#'   Structural models can be provided in the argument \code{stvar} if overidentifying constraints should be
#'   imposed.
#'
#'   Using a robust estimation method before switching to the variable metric can be useful if the intial
#'   estimates are not very close to the ML estimate of the structural model, as the variable metric algorithm
#'   (usually) converges to a nearby local maximum or saddle point. This is particularly the case when the imposed
#'   overidentifying restrictions are such that the unrestricted estimate is not close to satisfying them.
#'
#'   Employs the estimation function \code{optim} from the package \code{stats} that implements the optimization
#'   algorithms. See \code{?optim} for the documentation on the
#' @return Returns an object of class \code{'stvar'} defining the structural model
#' @seealso \code{\link{fitSTVAR}}, \code{\link{STVAR}}, \code{\link[stats]{optim}}
#' @inherit STVAR references
#' @examples
#' \donttest{
#' ## These are long running examples that use parallel computing!
#' ## Running the below examples takes approximately FILL IN HOW MANY minutes
#' }
#' @export

fitSSTVAR <- function(stvar, new_identification=c("recursive", "heteroskedasticity"), new_B_constraints=NULL,
                      maxit=100, robust_method=c("none", "Nelder-Mead", "SANN"), calc_std_errors=TRUE) {
  # Option for multiple estimation rounds? Useful when using robust methods? Or maybe lighter without it (ne need to parallelizes)
  # Option for two-phase estimation that initially uses Nelder-Mead or SANN (option to choose between them)?
  # HUOM! Jos tekee transition-function vertailupaperin, niin siinä tuskin edes bootstrappaillaan!
  check_stvar(stvar)
  stopifnot(maxit %% 1 == 0 & maxit >= 1)
  robust_method <- match.arg(robust_method)
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
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function,
                                         weightfun_pars=stvar$model$weightfun_pars, cond_dist=cond_dist)
  minval <- get_minval(stvar$data)
  old_npars <- length(stvar$params) # old npars, B_constraints not adjusted!

  # Check identification
  if(old_identification != "reduced_form") {
    if(old_identification == "recursive") {
      warning(paste("Recursively identified model supplied (not possible to impose overidentifying restrictions,",
                    "so the original model is returned"))
      return(stvar)
    } else if(old_identification != new_identification) {
      # Old model structural but identification changes
      stop(paste("The type of the identification of a structural models cannot be changed.",
                 "Consider using the initial reduced model to obtain the structural model with a specific identification."))
    }
  } else if(new_identification == "recursive") {
    # Old identification reduced form, recursively identified model should be returned.
    return(STVAR(data=data, p=p, M=M, d=d, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
                 parametrization=parametrization, identification=new_identification,
                 AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                 weight_constraints=weight_constraints, B_constraints=NULL))
  }


  # ret orig odel if old_ide

  identification <- new_identification # old in old_identification

  # Check the new B_constraints?
  B_constraints <- new_B_constraints # old in old_B_constraints
  check_constraints(p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    parametrization=parametrization, identification=identification,
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=B_constraints)
  npars <- n_params(p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    cond_dist=cond_dist, identification=identification,
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=B_constraints)

  # Set up initial values using the obtained parameters
  if(identification == "heteroskedasticity") {
    # Returns the model
  } else {
    # Not possible to arrive here currently, but will be updated later for new identification methods.
    stop("Unknown identification in fitSSTVAR")
  }

  ### Return the structural model without estimation if no estimation is required

  ### Estimation
  if(is.null(data)) {
    # No data, so returns the model with preliminary param values without estimation.
    message("There is no data for the model, so no estimation was conducted")
    # VAIKO! RETURNAA JO AIEMMINKIN JOS ESTIMOINTIA EI TARVITA ELI EI YLI-IDENTIFIOIVIA RAJOITTEITA
    # TÄLLÖIN TÄSSÄ KOHTAA VOISI VAIN STOPATA ERRORILLA.
  }


  # Function to optimize
  loglik_fn <- function(params) {
    tryCatch(loglikelihood(data=data, p=p, M=M, params=params,
                           weight_function=weight_function, weightfun_pars=weightfun_pars,
                           cond_dist=cond_dist, parametrization=parametrization,
                           identification=identification, AR_constraints=AR_constraints,
                           mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                           B_constraints=B_constraints, to_return="loglik", check_params=TRUE, minval=minval),
             error=function(e) minval)
  }

  # Estimation by robust method
  if(robust_method != "none") {

  }

  # Gradient of the log-likelihood function using central difference approximation
  h <- 6e-6
  I <- diag(rep(1, times=npars))
  loglik_grad <- function(params) {
    vapply(1:npars, function(i1) (loglik_fn(params + I[i1,]*h) - loglik_fn(params - I[i1,]*h))/(2*h), numeric(1))
  }

}


