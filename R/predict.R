#' @title Predict method for class 'stvar' objects
#'
#' @description \code{predict.stvar} is a predict method for class \code{'stvar'} objects.
#'
#' @inheritParams simulate.stvar
#' @param nsteps how many steps ahead should be predicted?
#' @param nsim to how many independent simulations should the forecast be based on?
#' @param pi a numeric vector specifying the confidence levels of the prediction intervals.
#' @param pred_type should the pointforecast be based on sample "median" or "mean"?
#' @param exo_weights if \code{weight_function="exogenous"}, provide a size \eqn{(nsteps x M)} matrix of exogenous
#'  transition weights for the regimes: \code{[step, m]} for \eqn{step} steps ahead and regime \eqn{m} weight. Ignored
#'  if \code{weight_function!="exogenous"}.
#' @details The forecasts are computed by simulating multiple sample paths of the future observations and
#'   using the sample medians or means as point forecasts and empirical quantiles as prediction intervals.
#' @return Returns a class '\code{stvarpred}' object containing, among the specifications,...
#'  \describe{
#'    \item{$pred}{Point forecasts}
#'    \item{$pred_ints}{Prediction intervals, as \code{[, , d]}.}
#'    \item{$trans_pred}{Point forecasts for the transition weights}
#'    \item{$trans_pred_ints}{Individual prediction intervals for transition weights, as \code{[, , m]}, m=1,..,M.}
#'  }
#' @seealso \code{\link{simulate.stvar}}
#' @inherit simulate.stvar references
#' @examples
#'  # p=2, M=2, d=2, Gaussian relative dens weights
#'  theta_222relg <- c(0.356914, 0.107436, 0.356386, 0.08633, 0.13996, 0.035172,
#'    -0.164575, 0.386816, 0.451675, 0.013086, 0.227882, 0.336084, 0.239257, 0.024173,
#'    -0.021209, 0.707502, 0.063322, 0.027287, 0.009182, 0.197066, 0.205831, 0.005157,
#'    0.025877, 1.092094, -0.009327, 0.116449, 0.592446)
#'  mod222relg <- STVAR(data=gdpdef, p=2, M=2, d=2, params=theta_222relg,
#'    weight_function="relative_dens")
#'
#'  # Predict 10 steps ahead, point forecast based on the conditional
#'  # mean and 90% prediction intervals; prediction based on 100 sample paths:
#'  pred1 <- predict(mod222relg, nsteps=10, nsim=100, pi=0.9, pred_type="mean")
#'  pred1
#'  plot(pred1)
#'
#'  # Predict 7 steps ahead, point forecast based on median and  90%, 80%,
#'  # and 70% prediction intervals; prediction based on 80 sample paths:
#'  pred2 <- predict(mod222relg, nsteps=7, nsim=80, pi=c(0.9, 0.8, 0.7),
#'   pred_type="median")
#'  pred2
#'  plot(pred2)
#' @export

predict.stvar <- function(object, ..., nsteps, nsim=1000, pi=c(0.95, 0.80), pred_type=c("mean", "median"),
                          exo_weights=NULL) {
  check_stvar(object, object_name="object")
  stvar <- object
  if(is.null(stvar$data)) stop("The model needs to contain data")
  stopifnot(all(pi > 0 & pi < 1))
  dat <- stvar$data
  pred_type <- match.arg(pred_type)
  if(!all_pos_ints(c(nsim, nsteps))) stop("nsim and n_ahaed should be positive integers")

  # Check the exogenous weights given for simulation
  if(stvar$model$weight_function == "exogenous") {
    check_exoweights(M=stvar$model$M, exo_weights=exo_weights, how_many_rows=nsteps, name_of_row_number="nsteps")
  }

  # Simulations
  simulations <- simulate.stvar(stvar, nsim=nsteps, init_values=dat, ntimes=nsim, drop=FALSE, exo_weights=exo_weights)
  sample <- simulations$sample
  alpha_mt <- simulations$transition_weights
  colnames(sample) <- colnames(dat)
  colnames(alpha_mt) <- vapply(1:stvar$model$M, function(m) paste("Regime", m), character(1))

  # Calculate quantiles from the third dimension of 3D simulation array
  dim3_quantiles <- function(x, q) {
    apply(x, c(1, 2), quantile, probs=q, names=FALSE)
  }

  # Calculate point forecasts
  if(pred_type == "mean") {
    pred <- rowMeans(sample, dims=2)
    trans_pred <- rowMeans(alpha_mt, dims=2)
  } else {
    pred <- dim3_quantiles(sample, q=0.5)
    trans_pred <- dim3_quantiles(alpha_mt, q=0.5)
  }
  if(is.null(colnames(pred))) colnames(pred) <- vapply(1:stvar$model$d, function(m) paste("Series", m), character(1))

  # Calculate prediction intervals
  lower <- (1 - pi)/2
  upper <- rev(1 - lower)
  q_tocalc <- c(lower, upper)
  pred_ints <- dim3_quantiles(sample, q_tocalc)
  trans_pred_ints <- dim3_quantiles(alpha_mt, q_tocalc)
  pred_ints <- aperm(pred_ints, perm=c(2, 1, 3))
  trans_pred_ints <- aperm(trans_pred_ints, perm=c(2, 1, 3))
  colnames(pred_ints) <- colnames(trans_pred_ints) <- q_tocalc

  # Return the results
  structure(list(stvar=stvar,
                 pred=pred,
                 pred_ints=pred_ints,
                 trans_pred=trans_pred,
                 trans_pred_ints=trans_pred_ints,
                 nsteps=nsteps,
                 nsim=nsim,
                 pi=pi,
                 pred_type=pred_type,
                 q=q_tocalc),
            class="stvarpred")
}
