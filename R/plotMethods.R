#' @import graphics
#' @describeIn STVAR plot method for class 'stvar'
#' @param x object of class \code{'stvar'}
#' @param ... currently not used.
#' @details The plot displays the time series together with estimated transition weights.
#' @export

plot.stvar <- function(x, ...) {
  stvar <- x
  data <- as.ts(stvar$data)
  n_obs <- nrow(data)
  p <- stvar$model$p
  M <- stvar$model$M
  d <- ncol(data)
  params <- stvar$params
  ts_tw <- ts(rbind(matrix(NA, nrow=p, ncol=M), stvar$transition_weights),
              start=start(data), frequency=frequency(data)) # First p observations are starting values

  # Time series and transition
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))
  graphics::par(mfrow=c(2, 1), mar=c(2.5, 2.5, 2.1, 1))
  colpal_ts <- grDevices::colorRampPalette(c("darkgreen", "darkblue", "darkmagenta", "red3"))(d)
  colpal_tw <- grDevices::colorRampPalette(c("blue", "turquoise1", "green", "red"))(M)
  names_ts <- colnames(data)
  names_tw <- paste0("Regime", 1:M)
  draw_legend <- function(nams, cols) {
    legend("topleft", legend=nams, bty="n", col=cols, lty=1, lwd=2, text.font=2, cex=0.6, x.intersp=0.5, y.intersp=1)
  }

  ts.plot(data, gpars=list(main="Time series", col=colpal_ts, lty=1:d))
  draw_legend(names_ts, cols=colpal_ts)
  ts.plot(ts_tw, gpars=list(main="Transition weights", ylim=c(0, 1), col=colpal_tw, lty=2))
  draw_legend(names_tw, cols=colpal_tw)
}


#' @import graphics
#'
#' @describeIn predict.stvar predict method
#' @inheritParams print.stvarpred
#' @param nt a positive integer specifying the number of observations to be plotted
#'   along with the forecast.
#' @param trans_weights should forecasts for transition weights be plotted?
#' @export

plot.stvarpred <- function(x, ..., nt, trans_weights=TRUE) {
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))
  stvarpred <- x
  data <- as.ts(stvarpred$stvar$data)
  n_obs <- nrow(data)
  d <- ncol(data)
  q <- stvarpred$q
  M <- stvarpred$stvar$model$M
  transition_weights <- stvarpred$stvar$transition_weights
  n_trans <- nrow(transition_weights)
  trans_weights <- trans_weights & M > 1 # Don't plot transition weights if M == 1
  if(missing(nt)) {
    nt <- round(nrow(data)*0.15)
  } else {
    stopifnot(nt > 0 & nt %% 1 == 0)
    if(nt > nrow(data)) {
      warning("nt > nrow(data), using nt = nrow(data)")
      nt <- nrow(data)
    }
  }
  if(trans_weights) {
    par(mfrow=c(d + 1, 1), mar=c(2.5, 2.5, 2.1, 1))
  } else {
    par(mfrow=c(d, 1), mar=c(2.5, 2.5, 2.1, 1))
  }
  make_ts <- function(x, trans=FALSE) { # Make ts that has the first value the same as the last value of the observed series/estim. m.weights.
    if(trans) {
      last_obs <- transition_weights[n_trans,]
    } else {
      last_obs <- data[n_obs,]
    }
    ts(rbind(last_obs, x), start=time(data)[n_obs], frequency=frequency(data))
  }
  ts_pred <- make_ts(stvarpred$pred)
  ts_trans_pred <- make_ts(stvarpred$trans_pred, trans=TRUE)
  ts_dat <- ts(data[(n_obs - nt + 1):n_obs,], start=time(data)[n_obs - nt + 1], frequency=frequency(data))
  ts_trans <- ts(transition_weights[(n_trans - nt + 1):n_trans,], start=time(data)[n_obs - nt + 1], frequency=frequency(data))
  t0 <- time(ts_pred)
  ts_names <- attributes(data)$dimnames[[2]]
  reg_names <- attributes(stvarpred$trans_pred)$dimnames[[2]]
  pred_ints <- aperm(stvarpred$pred_ints, perm=c(1, 3, 2)) # [step, series, quantiles]
  trans_pred_ints <- aperm(stvarpred$trans_pred_ints, perm=c(1, 3, 2)) # [step, series, quantiles]
  all_val <- lapply(1:d, function(j) c(ts_dat[,j], ts_pred[,j], simplify2array(pred_ints, higher=TRUE)[, j, ])) # All values to get ylims

  # Prediction intervals, we lapply through quantiles [, , q]
  ts_fun_fact <- function(inds) function(pred_ints, trans=FALSE) lapply(inds, function(i1) make_ts(pred_ints[, , i1], trans))
  ts1_lapply <- ts_fun_fact(1:(length(q)/2)) # Lower bounds
  ts2_lapply <- ts_fun_fact((length(q)/2 + 1):length(q)) # Upper bounds
  ints1 <- pred_ints
  ints1_trans <- trans_pred_ints
  ts1 <- ts1_lapply(ints1)
  ts2 <- ts2_lapply(pred_ints)
  if(trans_weights) {
    ts1_trans <- ts1_lapply(ints1_trans, trans=TRUE)
    ts2_trans <- ts2_lapply(trans_pred_ints, trans=TRUE)
  }

  # Plot forecasts for the series
  draw_poly <- function(ts1_or_ts2, pred_ts, col) polygon(x=c(t0, rev(t0)), y=c(ts1_or_ts2, rev(pred_ts)), col=col, border=NA)
  col_pred <- grDevices::rgb(0, 0, 1, 0.2)
  for(i1 in 1:d) {
    ts.plot(ts_dat[,i1], ts_pred[,i1], gpars=list(col=c("black", "blue"),
                                                  lty=1:2,
                                                  ylim=c(floor(min(all_val[[i1]])),
                                                         ceiling(max(all_val[[i1]]))),
                                                  main=ts_names[i1]))
    for(i2 in 1:length(stvarpred$pi)) {
      draw_poly(ts1[[i2]][,i1], ts_pred[,i1], col=col_pred)
      draw_poly(ts2[[i2]][,i1], ts_pred[,i1], col=col_pred)
    }
  }

  # Plot forecasts for the transition weights
  if(trans_weights) {

    # Point forecasts
    colpal_mw <- grDevices::colorRampPalette(c("blue", "turquoise1", "green", "red"))(M)
    colpal_mw2 <- grDevices::adjustcolor(colpal_mw, alpha.f=0.5)
    ts.plot(ts_trans, ts_trans_pred, gpars=list(col=c(colpal_mw2, colpal_mw),
                                            ylim=c(0, 1), lty=c(rep(1, M), rep(2, M)),
                                            main="transition weights"))
    legend("topleft", legend=paste0("Regime ", 1:M), bty="n", col=colpal_mw, lty=1, lwd=2,
           text.font=2, cex=0.9, x.intersp=0.5, y.intersp=1)

    # Individual prediction intervals as for the transition weights
    colpal_mw3 <- grDevices::adjustcolor(colpal_mw, alpha.f=0.2)
    for(m in 1:M) { # Go through regimes
      for(i2 in 1:length(stvarpred$pi)) { # Go through the prediction intervals
        draw_poly(ts1_trans[[i2]][,m], ts_trans_pred[,m], col=colpal_mw3[m])
        draw_poly(ts2_trans[[i2]][,m], ts_trans_pred[,m], col=colpal_mw3[m])
      }
    }
  }
  invisible(stvarpred)
}

