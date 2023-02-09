#' @title Pick \eqn{\phi_{m,0}} or \eqn{\mu_{m}}, m=1,..,M vectors
#'
#' @description \code{pick_phi0} picks the intercept or mean parameters from the given parameter vector.
#'
#' @param M the number of regimes
#' @param d the number of time series in the system, i.e., the dimension
#' @param params a real valued vector specifying the parameter values.
#'   Should have the form \eqn{\theta = (\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M,\sigma,\alpha,\nu)},
#'   where
#'   \itemize{
#'     \item{\eqn{\phi_{m,0} = } the \eqn{(d \times 1)} intercept (or mean) vector of the \eqn{m}th regime.}
#'     \item{\eqn{\varphi_m = (vec(A_{m,1}),...,vec(A_{m,p}))} \eqn{(pd^2 \times 1)}.}
#'     \item{\eqn{\sigma = (vech(\Omega_1),...,vech(\Omega_M)} \eqn{(Md(d - 1)/2 \times 1)}.}
#'     \item{\eqn{\alpha} contains the transition weights parameters}
#'     \item{\eqn{\nu > 2} is the degrees of freedom parameter that is included only if \code{cond_dist="Student"}.}
#'   }
#'   \describe{
#'     \item{For models with \code{weight_function="relative_dens"}:}{\eqn{\alpha = (\alpha_1,...,\alpha_{M-1})}
#'           \eqn{(M - 1 \times 1)}, where \eqn{\alpha_m} \eqn{(1\times 1), m=1,...,M-1} are the transition weight parameters.}
#'     \item{For models with \code{weight_function="logit"}:}{\eqn{\alpha = (\gamma_1,...,\gamma_M)} \eqn{((M-1)k\times 1)},
#'           where \eqn{\gamma_m} \eqn{(k\times 1), m=1,...,M-1} contains the logit-regression coefficients of the \eqn{m}th regime.}
#'   }
#'   Above, \eqn{\phi_{m,0}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}th coefficient matrix of the \eqn{m}th
#'   mixture component, and \eqn{\Omega_{m}} denotes the error term covariance matrix of the \eqn{m}:th mixture component.
#'   If \code{parametrization=="mean"}, just replace each \eqn{\phi_{m,0}} with regimewise mean \eqn{\mu_{m}}.
#'   \eqn{vec()} is vectorization operator that stacks columns of a given matrix into a vector. \eqn{vech()} stacks columns
#'   of a given matrix from the principal diagonal downwards (including elements on the diagonal) into a vector.
#' @return Returns a \eqn{(dxM)} matrix containing \eqn{\phi_{m,0}} in the m:th column or
#'   \eqn{\mu_{m}} if the parameter vector is mean-parametrized, \eqn{, m=1,..,M}.
#' @details Does not support constrained parameter vectors.
#' @section Warning:
#'  No argument checks!
#' @references
#'  \itemize{
#'    \item TO BE FILLED IN
#'  }
#' @keywords internal

pick_phi0 <- function(M, d, params) {
  matrix(params[1:(d*M)], nrow=d, byrow=FALSE)
}


#' @title Pick coefficient matrix
#'
#' @description \code{pick_Ami} picks the coefficient matrix \eqn{A_{m,i}} from the given parameter vector.
#'
#' @inheritParams pick_phi0
#' @param p the autoregressive order of the model
#' @param m which regime?
#' @param i which lag in \eqn{1,...,p}?
#' @param unvec if \code{FALSE} then vectorized version of \eqn{A_{m,i}} will be returned instead of matrix.
#'   Default if \code{TRUE}.
#' @inherit pick_phi0 details references
#' @inheritSection pick_phi0 Warning
#' @return Returns the i:th lag coefficient matrix of m:th regime, \eqn{A_{m,i}}.
#' @keywords internal

pick_Ami <- function(p, M, d, params, m, i, unvec=TRUE) {
  qm1 <- d*M + d^2*p*(m - 1)
  Ami <- params[(qm1 + d^2*(i - 1) + 1):(qm1 + d^2*i)]
  if(unvec) {
    return(unvec(d=d, a=Ami))
  } else {
    return(Ami)
  }
}

