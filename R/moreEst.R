#' @title Maximum likelihood estimation of a reduced form or structural STVAR model based on preliminary estimates
#'
#' @description \code{iterate_more} uses a variable metric algorithm to estimate a reduced form or structural STVAR model
#'  (object of class \code{'stvar'}) based on preliminary estimates.
#'
#' @inheritParams simulate.stvar
#' @inheritParams fitSTVAR
#' @inheritParams STVAR
#' @details The purpose of \code{iterate_more} is to provide a simple and convenient tool to finalize
#'   the estimation when the maximum number of iterations is reached when estimating a STVAR model
#'   with the main estimation function \code{fitSTVAR} or \code{fitSSTVAR}.
#' @return Returns an object of class \code{'stvar'} defining the estimated model
#' @seealso \code{\link{fitSTVAR}}, \code{\link{fitSSTVAR}}, \code{\link{STVAR}}, \code{\link[stats]{optim}}
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
  B_constraints <- stvar$model$B_constraints
  minval <- get_minval(stvar$data)
  npars <- length(stvar$params)

  # Function to optimize
  loglik_fn <- function(params) {
    tryCatch(loglikelihood(data=data, p=p, M=M, params=params, weight_function=weight_function,
                           cond_dist=cond_dist, parametrization=parametrization,
                           identification=identification, AR_constraints=AR_constraints,
                           mean_constraints=mean_constraints, B_constraints=B_constraints,
                           to_return="loglik", check_params=TRUE, minval=minval), error=function(e) minval)
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
               cond_dist=cond_dist,
               parametrization=parametrization,
               identification=identification,
               AR_constraints=AR_constraints,
               mean_constraints=mean_constraints,
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
