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
