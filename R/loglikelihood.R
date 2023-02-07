
#' @title Log-likelihood function
#'
#' @description \code{loglikelihood} log-likelihood function of a smooth transition VAR model
#'
#' @param data data a matrix or class \code{'ts'} object with \code{d>1} columns. Each column is taken to represent
#'  a univariate time series. Missing values are not supported.
#' @param p a positive integer specifying the autoregressive order
#' @param M a positive integer specifying the number of regimes
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
#' @param weight_function what type of transition weights should be used?
#' @param cond_dist specifies the conditional distribution of the model as \code{"Gaussian"} or \eqn{"Student"}.
#' @param parametrization \code{"intercept"} or \code{"mean"} determining whether the model is parametrized with intercept
#'   parameters \eqn{\phi_{m,0}} or regime means \eqn{\mu_{m}}, m=1,...,M.
#' @param identification is it reduced form model or an identified structural model; if the latter, how is it identified?
#' @param to_return should the returned object be the log-likelihood, which is the default, or something else?
#'   See the section "Return" for all the options.
#' @param check_params should it be checked that the parameter vector satisfies the model assumptions? Can be skipped to save
#'   computation time if it does for sure.
#' @param minval the value that will be returned if the parameter vector does not lie in the parameter space
#'   (excluding the identification condition).
#' @param stat_tol numerical tolerance for stability of condition of the regimes: if the "bold A" matrix of any regime
#'   has eigenvalues larger that \code{1 - stat_tol} the parameter is considered to be outside the parameter space.
#'   Note that if tolerance is too small, numerical evaluation of the log-likelihood might fail and cause error.
#' @param posdef_tol numerical tolerance for positive definiteness of the error term covariance matrices: if
#'   the error term covariance matrix of any regime has eigenvalues smaller than this, the parameter is considered
#'   to be outside the parameter space. Note that if the tolerance is too small, numerical evaluation of the
#'   log-likelihood might fail and cause error.
#' @param df_tol the parameter vector is considered to be outside the parameter space if the degrees of
#'   freedom parameters is not larger than \code{2 + df_tol} (applies only if \code{cond_dist="Student"}).
#' @details FILL IN
#' @return
#'   \describe{
#'     \item{If \code{to_return="loglik"}:}{the log-likelihood of the specified model}
#'   }
#' @references
#'  \itemize{
#'    \item Lütkepohl H. 2005. New Introduction to Multiple Time Series Analysis,
#'          \emph{Springer}.
#'    \item McElroy T. 2017. Computation of vector ARMA autocovariances.
#'          \emph{Statistics and Probability Letters}, \strong{124}, 92-96.
#'  }
#' @keywords internal

loglikelihood <- function(data, p, M, params, weight_function=c("relative_dens", "logit"), cond_dist=c("Gaussian", "Student"),
                          parametrization=c("intercept", "mean"),
                          identification=c("reduced_form", "impact_responses", "heteroskedasticity"),
                          AR_constraints=NULL, mean_constraints=NULL, B_constraints=NULL,
                          to_return=c("loglik"), check_params=TRUE, minval=NULL,
                          stat_tol=1e-3, posdef_tol=1e-8, df_tol=1e-8) {
  # Match args
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  identification <- match.arg(identification)
  to_return <- match.arg(to_return)
  if(identification != "reduced form") stop("Only reduced form models are currently supported")
  if(!is.null(AR_constraints) || !is.null(mean_constraints) || !is.null(B_constraints)) stop("Constrained models are not
                                                                                             currently supported")

  # Compute some required statistics
  epsilon <- round(log(.Machine$double.xmin) + 10) # Logarithm of the smallest value that can be handled normally
  d <- ncol(data)
  n_obs <- nrow(data)
  T_obs <- n_obs - p

  # Collect the parameter values
  # First remove all constraints, if any, TO BE IMPLEMENTED

  # Check that the parameter vector lies in the parameter space
}
