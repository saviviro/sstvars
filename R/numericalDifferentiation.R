#' @title Calculate gradient or Hessian matrix
#'
#' @description \code{calc_gradient} or \code{calc_hessian} calculates the gradient or Hessian matrix
#'   of the given function at the given point using central difference numerical approximation.
#'   \code{get_gradient} or \code{get_hessian} calculates the gradient or Hessian matrix of the
#'   log-likelihood function at the parameter estimates of a class \code{'FILL IN'} object. \code{get_soc}
#'   returns eigenvalues of the Hessian matrix, and \code{get_foc} is the same as \code{get_gradient}
#'   but named conveniently.
#'
#' @param x a numeric vector specifying the point where the gradient or Hessian should be calculated.
#' @param fn a function that takes in argument \code{x} as the \strong{first} argument.
#' @param h difference used to approximate the derivatives.
#' @param ... other arguments passed to \code{fn}
#' @details In particular, the functions \code{get_foc} and \code{get_soc} can be used to check whether
#'   the found estimates denote a (local) maximum point, a saddle point, or something else. Note that
#'   profile log-likelihood functions can be conveniently plotted with the function \code{profile_logliks}
#'   PROFILE_LOGLIKS IS NOT YET IMPLEMENTED.
#' @return Gradient functions return numerical approximation of the gradient and Hessian functions return
#'   numerical approximation of the Hessian. \code{get_soc} returns eigenvalues of the Hessian matrix.
#' @section Warning:
#'   No argument checks!
#' @examples
#'   # Simple function
#'   foo <- function(x) x^2 + x
#'   calc_gradient(x=1, fn=foo)
#'   calc_gradient(x=-0.5, fn=foo)
#'
#'   # More complicated function
#'   foo <- function(x, a, b) a*x[1]^2 - b*x[2]^2
#'   calc_gradient(x=c(1, 2), fn=foo, a=0.3, b=0.1)
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

get_gradient <- function(FILL_IN_OBJECT) {
  foo <- function(x) {
    # LOGLIKELIHOOD FUNCTION AS A FUNCTION OF THE PARAMETER VECTOR
  }
  # calc_gradient(x=FILL_IN_OBJECT$params, fn=foo)
}



#' @rdname calc_gradient
#' @export

get_hessian <- function(FILL_IN_OBJECT) {
  foo <- function(x) {
    # LOGLIKELIHOOD FUNCTION AS A FUNCTION OF THE PARAMETER VECTOR
  }
  # calc_hessian(x=FILL_IN_OBJECT$params, fn=foo)
}

#' @rdname calc_gradient
#' @export
get_foc <- function(FILL_IN_OBJECT) {
  # get_gradient(FILL_IN_OBJECT)
}

#' @rdname calc_gradient
#' @export
get_soc <- function(FILL_IN_OBJECT) {
  # eigen(get_hessian(FILL_IN_OBJECT))$value
}

