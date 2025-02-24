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
#'   described, for example, in Virolainen (2025).
#'
#' @inheritParams skewed_t_dens
#' @details See Virolainen (2025) and the references therein, for example, for the details of the density function of
#'   t-distribution with zero mean and unit variance (assume the skewness parameter value is zero to obtain the
#'   non-skewed version of the t-distribution).
#' @return Returns a numeric vector of the same length as \code{y} containing the density values.
#' @references
#'  \itemize{
#'    \item Virolainen S. 2025. Identification by non-Gaussianity in structural threshold and
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


#' @title Generate random samples from the skewed t-distribution
#'
#' @description \code{generate_skewed_t} generates \code{n} random observations from the univariate skewed \emph{t}-distribution
#'  described in Hansen (1994) using the acceptance-rejection sampling method.
#'
#' @param n An integer specifying the number of random observations to generate. Must be a positive integer.
#' @param nu A numeric scalar specifying the degrees of freedom parameter for the skewed \emph{t}-distribution. Must be greater than 2.
#' @param lambda A numeric scalar specifying the skewness parameter for the skewed \emph{t}-distribution. Must be between \eqn{-1} and \eqn{1}.
#' @param bc_M An optional numeric scalar specifying the bounding constant \eqn{M} used in the acceptance-rejection algorithm.
#'   If not provided, it is computed using \code{\link{bounding_const_M}} with the given \code{nu} and \code{lambda}.
#' @details The function implements the acceptance-rejection algorithm to generate random samples from the skewed \emph{t}-distribution.
#'   The proposal distribution used is a standard \emph{t}-distribution with degrees of freedom \code{proposal_nu}, which is set to \eqn{3}
#'   when \code{nu > 3} to ensure heavier tails and accommodate the skewness of the target distribution.
#'
#'   If \code{bounding_const_M} is not provided, it is calculated using the \code{\link{bounding_const_M}} function. It is important that
#'   the same proposal distribution is used in both the computation of \code{bounding_const_M} and the acceptance-rejection sampling
#'   algorithm to ensure correctness.
#' @return A numeric vector of length \code{n} containing random observations from the skewed \emph{t}-distribution with
#'   parameters \code{nu} and \code{lambda}.
#' @inherit skewed_t_dens references
#' @keywords internal

generate_skewed_t <- function(n, nu, lambda, bc_M) {
  if(missing(bc_M)) {
    bc_M <- bounding_const_M(nu=nu, lambda=lambda)
  } # Note that the same proposal dist needs to be used in computing M as well as in the sampling algorithm
  proposal_nu <- ifelse(nu > 3, 3, nu) # Proposal dist needs to have heavy tails to accommodate skewness in the target dist
  y_vals <- numeric(n)

  # Iterate throught the algorithm
  for(i1 in 1:n) {
    accept <- FALSE
    while(!accept) {
      y_cand <- sqrt((proposal_nu - 2)/proposal_nu)*rt(1, df=proposal_nu)
      u <- runif(1)
      if(u < skewed_t_dens(y_cand, nu=nu, lambda=lambda)/(bc_M*stand_t_dens(y_cand, nu=proposal_nu))) {
        y_vals[i1] <- y_cand
        accept <- TRUE
      }
    }
  }
  y_vals # Sample of n values from the skewed t distribution
}
