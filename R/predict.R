#' @title Predict method for class 'stvar' objects
#'
#' @description \code{predict.stvar} is a predict method for class \code{'stvar'} objects. The forecasts
#'   are computed by simulating multiple sample paths of the future observations and using the
#'   sample medians or means as point forecasts and empirical quantiles as prediction intervals.
#'
#' @inheritParams simulate.stvar
#' @param nsteps how many steps ahead should be predicted?
#' @param nsim to how many independent simulations should the forecast be based on?
#' @param pi a numeric vector specifying the confidence levels of the prediction intervals.
#' @param pred_type should the pointforecast be based on sample "median" or "mean"?
#' @return Returns a class '\code{stvarpred}' object containing, among the specifications,...
#'  \describe{
#'    \item{$pred}{Point forecasts}
#'    \item{$pred_int}{Prediction intervals, as \code{[, , d]}.}
#'    \item{$trans_pred}{Point forecasts for the transition weights}
#'    \item{$trans_pred_int}{Individual prediction intervals for transition weights, as \code{[, , m]}, m=1,..,M.}
#'  }
#' @seealso \code{\link{simulate.stvar}}
#' @inherit simulate.stvar references
#' @examples
#' # FILL IN EXAMPLES
#' @export

predict.stvar <- function(object, ..., nsteps, nsim=1000, pi=c(0.95, 0.80), pred_type=c("mean", "median")) {
  stvar <- object
  if(is.null(stvar$data)) stop("The model needs to contain data")
  stopifnot(all(pi > 0 & pi < 1))
  dat <- stvar$data
  pred_type <- match.arg(pred_type)
  if(!all_pos_ints(nsim, nsteps)) stop("nsim and n_ahaed should be positive integers")

  # Simulations
  simulations <- simulate(stvar, nsim=nsteps, init_values=dat, ntimes=nsim)
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
  mix_pred_ints <- dim3_quantiles(alpha_mt, q_tocalc)

  # Return the results
  structure(list(stvar=stvar,
                 pred=pred,
                 pred_int=pred_ints,
                 trans_pred=trans_pred,
                 trans_pred_ints=trans_pred_ints,
                 nsteps=nsteps,
                 nsim=nsim,
                 pi=pi,
                 pred_type=pred_type,
                 q=q_tocalc),
            class="stvarpred")
}
