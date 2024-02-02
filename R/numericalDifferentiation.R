#' @title Calculate gradient or Hessian matrix
#'
#' @description \code{calc_gradient} or \code{calc_hessian} calculates the gradient or Hessian matrix
#'   of the given function at the given point using central difference numerical approximation.
#'   \code{get_gradient} or \code{get_hessian} calculates the gradient or Hessian matrix of the
#'   log-likelihood function at the parameter estimates of a class \code{'stvar'} object. \code{get_soc}
#'   returns eigenvalues of the Hessian matrix, and \code{get_foc} is the same as \code{get_gradient}
#'   but named conveniently.
#'
#' @inheritParams get_boldA_eigens
#' @param x a numeric vector specifying the point where the gradient or Hessian should be calculated.
#' @param fn a function that takes in argument \code{x} as the \strong{first} argument.
#' @param h difference used to approximate the derivatives.
#' @param ... other arguments passed to \code{fn}
#' @details In particular, the functions \code{get_foc} and \code{get_soc} can be used to check whether
#'   the found estimates denote a (local) maximum point, a saddle point, or something else. Note that
#'   profile log-likelihood functions can be conveniently plotted with the function \code{profile_logliks}.
#' @return Gradient functions return numerical approximation of the gradient and Hessian functions return
#'   numerical approximation of the Hessian. \code{get_soc} returns eigenvalues of the Hessian matrix.
#' @section Warning:
#'   No argument checks!
#' @examples
#' # Create a simple function:
#' foo <- function(x) x^2 + x
#'
#' # Calculate the gradient at x=1 and x=-0.5:
#' calc_gradient(x=1, fn=foo)
#' calc_gradient(x=-0.5, fn=foo)
#'
#' # Create a more complicated function
#' foo <- function(x, a, b) a*x[1]^2 - b*x[2]^2
#'
#' # Calculate the gradient at x=c(1, 2) with parameter values a=0.3 and b=0.1:
#' calc_gradient(x=c(1, 2), fn=foo, a=0.3, b=0.1)
#'
#' # Create a linear Gaussian VAR p=1 model:
#' theta_112 <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103,
#'  0.601786, -0.002945, 0.067224)
#' mod112 <- STVAR(data=gdpdef, p=1, M=1, params=theta_112)
#'
#' # Calculate the gradient of the log-likelihood function about the parameter values:
#' get_foc(mod112)
#'
#' # Calculate the eigenvalues of the Hessian matrix of the log-likelihood function
#' # about the parameter values:
#' get_soc(mod112)
#' @export

calc_gradient <- function(x, fn, h=6e-06, ...) {
  fn <- match.fun(fn)
  n <- length(x)
  I <- diag(1, nrow=n, ncol=n)
  h <- rep(h, times=n)
  vapply(1:n, function(i1) (fn(x + h[i1]*I[i1,], ...) - fn(x - h[i1]*I[i1,], ...))/(2*h[i1]), numeric(1))
}


#' @rdname calc_gradient
#' @export

calc_hessian <- function(x, fn, h=6e-06, ...) {
  fn <- match.fun(fn)
  n <- length(x)
  I <- diag(1, nrow=n, ncol=n)
  h <- rep(h, times=n)
  Hess <- matrix(ncol=n, nrow=n)
  for(i1 in 1:n) {
    for(i2 in i1:n) {
      dr1 <- (fn(x + h[i1]*I[i1,] + h[i2]*I[i2,], ...) - fn(x - h[i1]*I[i1,] + h[i2]*I[i2,], ...))/(2*h[i1])
      dr2 <- (fn(x + h[i1]*I[i1,] - h[i2]*I[i2,], ...) - fn(x - h[i1]*I[i1,] - h[i2]*I[i2,], ...))/(2*h[i1])
      Hess[i1, i2] <- (dr1 - dr2)/(2*h[i2])
      Hess[i2, i1] <- Hess[i1, i2] # Take use of symmetry
    }
  }
  Hess
}


#' @rdname calc_gradient
#' @export

get_gradient <- function(stvar) {
  check_stvar(stvar)
  foo <- function(x) {
    # Log-likelihood function as a function of the parameter
    loglikelihood(data=stvar$data, p=stvar$model$p, M=stvar$model$M, params=x,
                  weight_function=stvar$model$weight_function,
                  weightfun_pars=stvar$model$weightfun_pars,
                  cond_dist=stvar$model$cond_dist,
                  parametrization=stvar$model$parametrization,
                  identification=stvar$model$identification,
                  AR_constraints=stvar$model$AR_constraints,
                  mean_constraints=stvar$model$mean_constraints,
                  weight_constraints=stvar$model$weight_constraints,
                  B_constraints=stvar$model$B_constraints,
                  to_return="loglik", minval=NA)
  }
  calc_gradient(x=stvar$params, fn=foo)
}



#' @rdname calc_gradient
#' @export

get_hessian <- function(stvar) {
  check_stvar(stvar)
  foo <- function(x) {
    # Log-likelihood function as a function of the parameter
    loglikelihood(data=stvar$data, p=stvar$model$p, M=stvar$model$M, params=x,
                  weight_function=stvar$model$weight_function,
                  weightfun_pars=stvar$model$weightfun_pars,
                  cond_dist=stvar$model$cond_dist,
                  parametrization=stvar$model$parametrization,
                  identification=stvar$model$identification,
                  AR_constraints=stvar$model$AR_constraints,
                  mean_constraints=stvar$model$mean_constraints,
                  weight_constraints=stvar$model$weight_constraints,
                  B_constraints=stvar$model$B_constraints,
                  to_return="loglik", minval=NA)
  }
  calc_hessian(x=stvar$params, fn=foo)
}

#' @rdname calc_gradient
#' @export
get_foc <- function(stvar) {
  get_gradient(stvar)
}

#' @rdname calc_gradient
#' @export
get_soc <- function(stvar) {
  eigen(get_hessian(stvar))$values
}

