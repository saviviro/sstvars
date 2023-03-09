#' @title Returns the default smallest allowed log-likelihood for given data.
#'
#' @description \code{get_minval} returns the default smallest allowed log-likelihood for given data.
#'
#' @inheritParams GAfit
#' @return Returns \code{-(10^(ceiling(log10(nrow(data)) + ncol(data))) - 1)}
#' @keywords internal

get_minval <- function(data) {
  -(10^(ceiling(log10(nrow(data)) + ncol(data))) - 1)
}


#' @title Calculate AIC, HQIC, and BIC
#'
#' @description \code{get_IC} calculates the information criteria values
#'   AIC, HQIC, and BIC divided by the number of observations.
#'
#' @param loglik log-likelihood value
#' @param npars number of (freely estimated) parameters in the model
#' @param T_obs numbers of observations with the \eqn{p} starting values excluded.
#' @return Returns a data frame containing the information criteria values
#'   divided by the number of observations.
#' @keywords internal

get_IC <- function(loglik, npars, T_obs) {
  AIC <- (-2*loglik + 2*npars)/T_obs
  HQIC <- (-2*loglik + 2*npars*log(log(T_obs)))/T_obs
  BIC <- (-2*loglik + npars*log(T_obs))/T_obs
  data.frame(AIC=AIC, HQIC=HQIC, BIC=BIC)
}
