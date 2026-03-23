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
#' @param loglik log-likelihood value, penalized or non-penalized one.
#' @param npars number of (freely estimated) parameters in the model
#' @param T_obs numbers of observations with the \eqn{p} starting values excluded.
#' @return Returns a data frame containing the information criteria values
#'   divided by the number of observations.
#' @internal

get_IC <- function(loglik, npars, T_obs) {
  AIC <- (-2*loglik + 2*npars)/T_obs
  HQIC <- (-2*loglik + 2*npars*log(log(T_obs)))/T_obs
  BIC <- (-2*loglik + npars*log(T_obs))/T_obs
  data.frame(AIC=AIC, HQIC=HQIC, BIC=BIC)
}

#' @title Calculate penalized AIC, HQIC, and BIC
#' @description \code{get_penalized_IC} calculates the penalized information criteria values
#'  AIC, HQIC, and BIC divided by the number of observations, i.e., using the penalized log-likelihood
#'  function if that was used for estimation.
#' @param sstvar an object of class \code{'sstvar'}.
#' @return Returns a data frame containing the penalized information criteria values.
#' @examples
#' # Logistic Student's t STVAR with p=1, M=2, and the first lag of the second variable
#' # as the switching variable, with penalized log-likelihood, with very strong penalization:
#' params12 <- c(0.62906848, 0.14245295, 2.41245785, 0.66719269, 0.3534745, 0.06041779, -0.34909745,
#'   0.61783824, 0.125769, -0.04094521, -0.99122586, 0.63805416, 0.371575, 0.00314754, 0.03440824,
#'   1.29072533, -0.06067807, 0.18737385, 1.21813844, 5.00884263, 7.70111672)
#' mod12 <- STVAR(data=gdpdef, p=1, M=2, params=params12, weight_function="logistic",
#'   weightfun_pars=c(2, 1), cond_dist="Student", penalized=TRUE, penalty_params=c(0.4, 0.4))
#'
#' mod12$IC # IC without penalization term in loglik
#' get_penalized_IC(mod12) # IC with penalization term in loglik
#' @export

get_penalized_IC <- function(stvar) {
  if(!stvar$penalized) {
    message("The model was not estimated with penalized log-likelihood, so nonpenalized information criteria are calculated.")
  }
  get_IC(loglik=stvar$loglik, npars=length(stvar$params), T_obs=nrow(stvar$data) - stvar$model$p)
}


#' @title Calculate "distance" between two (scaled) regimes
#'  \strong{\eqn{\upsilon_{m}}}\eqn{ = (\phi_{m},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m})}
#'
#' @description \code{regime_distance} calculates "distance" between two scaled regimes.
#'
#' @param regime_pars1 a length \eqn{pd^2+d+d(d+1)/2} vector
#'   \strong{\eqn{\upsilon_{m}}}\eqn{ = (\phi_{m},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m})}.
#' @param regime_pars2 a length \eqn{pd^2+d+d(d+1)/2} vector
#'   \strong{\eqn{\upsilon_{m}}}\eqn{ = (\phi_{m},}\strong{\eqn{\phi_{m}}}\eqn{,\sigma_{m})}.
#' @return Returns "distance" between \code{regime_pars1} and \code{regime_pars2}. Values are scaled
#'   before calculating the "distance". Read the source code for more details.
#' @section Warning:
#'  No argument checks!
#' @inherit in_paramspace references
#' @keywords internal

regime_distance <- function(regime_pars1, regime_pars2) {
  dist_fun <- function(x) {
    x <- abs(x)
    if(x < 1) {
      return(1)
    } else {
      return(10^ceiling(abs(log10(x))))
    }
  }
  scales1 <- vapply(regime_pars1, dist_fun, numeric(1))
  scales2 <- vapply(regime_pars2, dist_fun, numeric(1))
  c(sqrt(crossprod(regime_pars1/scales1 - regime_pars2/scales2)))
}


#' @title Get the new starting time of series that is forwarded some number of steps
#'
#' @description \code{get_new_start} calculates the new starting time of series
#'   that is forwarded some number of steps.
#'
#' @param y_start original starting time of the series
#' @param y_freq frequency of the series
#' @param steps_forward how many steps the series should be forwarded?
#' @return Returns a length two numeric vector with the "year" (or "major")
#'   time point in the first element the "quarter/month/week/day" (or "minor")
#'   time in the second element for a series that is forwarded from \code{y_start}
#'   \code{steps_forward} steps forward.
#' @keywords internal

get_new_start <- function(y_start, y_freq, steps_forward) {
  majors_forward <- steps_forward %/% y_freq
  minors_forward <- steps_forward %% y_freq
  new_start <- y_start + c(majors_forward, minors_forward)
  if(new_start[2] > y_freq) {
    new_start <- c(new_start[1] + 1, new_start[2] %% y_freq)
  }
  new_start
}
