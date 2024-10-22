#' @title The density function of the univariate skewed t distribution
#'
#' @description \code{skewed_t_dens} calculates the density of the univariate skewed t distribution described in Hansen (1994).
#'
#' @param y a numeric vector containing the values at which the density is to be evaluated.
#' @param nu the degrees of freedom parameter value, a numeric scalar strictly larger than two.
#' @param lambda the skewness parameter value, a numeric scalar strictly between -1 and 1.
#' @details See Hansen (1994, Section 2.4) for the details of the skewed t distribution.
#' @return Returns a numeric vector of the same length as \code{y} containing the density values.
#' @references
#'  \itemize{
#'    \item Hansen B.E. 1994. Autoregressive Conditional Density estimation. \emph{Journal of Econometrics}, \strong{35}:3, 705-730.
#'  }
#' @keywords internal

skewed_t_dens <- function(y, nu, lambda) {
  logc_i <- lgamma(0.5*(1 + nu)) - 0.5*log(base::pi) - 0.5*log(nu - 2) - lgamma(0.5*nu) # (d x 1 )
  a_i <- 4*lambda*exp(logc_i)*(nu - 2)/(nu - 1) # (d x 1)
  logb_i <- 0.5*log(1 + 3*lambda^2 - a_i^2) # (d x 1)
  b_i <- exp(logb_i) # (d x 1)
  b_i*exp(logc_i)*(1 + 1/(nu - 2)*((b_i*y + a_i)/(1 + ifelse(y < -a_i/b_i, -lambda, lambda)))^2)^(-0.5*(1 + nu))
}


#' @title The density function of the univariate t distribution with zero mean and unit variance
#'
#' @description \code{stand_t_dens} calculates the density of the univariate t distribution with zero mean and unit variance,
#'   described, for example, in Virolainen (2024).
#'
#' @inheritParams skewed_t_dens
#' @details See Virolainen (2024) and the references therein, for example, for the details of the density function of
#'   t-distribution with zero mean and unit variance.
#' @return Returns a numeric vector of the same length as \code{y} containing the density values.
#' @references
#'  \itemize{
#'    \item Virolainen S. 2024. Identification by non-Gaussianity in structural threshold and
#'       smooth transition vector autoregressive models. Unpublished working
#'       paper, available as arXiv:2404.19707.
#'  }
#' @keywords internal

stand_t_dens <- function(y, nu) {
  logc_i <- lgamma(0.5*(1 + nu)) - 0.5*log(base::pi) - 0.5*log(nu - 2) - lgamma(0.5*nu) # (d x 1 )
  exp(logc_i)*(1 + y^2/(nu - 2))^(-0.5*(1 + nu))
}


#' @title Compute the bounding constant for acceptance-rejection sampling
#'
#' @description \code{bounding_const_M} calculates the bounding constant \eqn{M} used in the acceptance-rejection sampling algorithm
#'  for the univariate skewed \emph{t}-distribution described in Hansen (1994)
#'
#' @inheritParams skewed_t_dens
#' @details The function computes the bounding constant \eqn{M} required for the acceptance-rejection sampling method by evaluating
#'  the ratio of the skewed \emph{t}-density (\code{\link{skewed_t_dens}}) to the standard \emph{t}-density (\code{\link{stand_t_dens}})
#'  over a grid of \eqn{y} values ranging from \eqn{-10} to \eqn{10}. To improve the efficiency of the sampling algorithm, the degrees
#'  of freedom parameter for the proposal distribution is set to the minimum of \code{nu} and \eqn{3}, ensuring heavier tails in the
#'  proposal distribution when \code{nu} is large. A safety margin of 10\% is added to the maximum ratio to account for numerical
#'  inaccuracies and ensure that the inequality \eqn{f(y) \leq M \cdot q(y)} holds over the entire support.
#'
#' @return Returns a numeric scalar representing the estimated bounding constant \eqn{M} to be used in the acceptance-rejection
#'  sampling algorithm.
#' @inherit skewed_t_dens references
#' @keywords internal

bounding_const_M <- function(nu, lambda) {
  proposal_nu <- ifelse(nu > 3, 3, nu) # Proposal dist needs to have heavy tails to accommodate skewness in the target dist
  ratio_fun <- function(y) {
    skewed_t_dens(y, nu=nu, lambda=lambda)/stand_t_dens(y, nu=proposal_nu)
  }
  y_vals <- seq(from=-10, to=10, length.out=1000)
  ratios <- vapply(y_vals, FUN=ratio_fun, FUN.VALUE=numeric(1))
  max(ratios, na.rm=TRUE)*1.1 # Estimate of the bounding constant M with 10% safety margin
}
