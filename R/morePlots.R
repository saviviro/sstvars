#' @import graphics
#'
#' @title Plot structural shock time series of a STVAR model
#'
#' @description \code{plot_struct_shocks} plots structural shock time series of a structural STVAR model.
#'   For reduced form models (not identified by non-Gaussianity), recursive identification is assumed.
#'
#' @inheritParams diagnostic_plot
#' @details Plot the time series of the structural shocks of a structural STVAR model.
#' @return No return value, called for its side effect of plotting the structural shock time series.
#' @inherit get_residuals references
#' @seealso \code{\link{diagnostic_plot}}, \code{\link{fitSTVAR}}, \code{\link{fitSSTVAR}}, \code{\link{STVAR}},
#' @examples
#' ## Gaussian STVAR p=1, M=2 model, with weighted relative stationary densities
#' # of the regimes as the transition weight function:
#' theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863,
#'   -0.3124, 0.505251, 0.298483, 0.030096, -0.176925, 0.838898, 0.310863, 0.007512,
#'   0.018244, 0.949533, -0.016941, 0.121403, 0.573269)
#' mod122 <- STVAR(data=gdpdef, p=1, M=2, params=theta_122relg)
#'
#' # Plot the times series structural shocks assuming recursive identification:
#' plot_struct_shocks(mod122)
#' @export

plot_struct_shocks <- function(stvar) {
  check_stvar(stvar)
  cond_dist <- stvar$model$cond_dist
  if(stvar$model$identification == "reduced_form" && !cond_dist %in% c("ind_Student", "ind_skewed_t")) {
    stvar$model$identification <- "recursive"
  }
  struct_shocks <- get_residuals(data=stvar$data, p=stvar$model$p, M=stvar$model$M, params=stvar$params,
                                 weight_function=stvar$model$weight_function, weightfun_pars=stvar$model$weightfun_pars,
                                 cond_dist=stvar$model$cond_dist, identification=stvar$model$identification,
                                 AR_constraints=stvar$model$AR_constraints, mean_constraints=stvar$model$mean_constraints,
                                 weight_constraints=stvar$model$weight_constraints, B_constraints=stvar$model$B_constraints,
                                 penalized=stvar$penalized, penalty_params=stvar$penalty_params, allow_unstab=stvar$allow_unstab,
                                 standardize=FALSE, structural_shocks=TRUE)

  d <- stvar$model$d
  p <- stvar$model$p
  names_ts <- colnames(as.ts(stvar$data))
  colnames(struct_shocks) <- names_ts
  struct_shocks <- ts(struct_shocks, start=get_new_start(y_start=start(stvar$data), y_freq=frequency(stvar$data), steps_forward=p),
                      frequency=frequency(stvar$data))
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))

  empty_plot <- function() plot(0, xaxt='n', yaxt='n', bty='n', pch='', ylab='', xlab='')
  par(mar=c(0, 3, 0, 0.5), las=1)
  layout(matrix(seq_len(d + 2), ncol=1), heights=c(0.35, rep(1, times=d), 0.35))
  empty_plot()
  text(x=1, y=0.0, "Time series of structural shocks", cex=1.5, font=2)
  for(d1 in 1:d) {
    xaxt <- "n"
    if(d1 == d) {
      xaxt <- "s"
    }
    yaxt1 <- round(min(struct_shocks[,d1]))
    yaxt2 <- round(max(struct_shocks[,d1]))
    #main <- ifelse(d1 == 1, "Structural shock time series", "")
    plot(struct_shocks[,d1], yaxt="n", xaxt=xaxt, type="l", col=grDevices::rgb(0, 0, 0, 1), ylab="", xlab="", main="")
    axis(2, at=yaxt1:yaxt2, labels=yaxt1:yaxt2)
    abline(h=0, col=grDevices::rgb(1, 0, 0, 0.3), lwd=2)
    legend("topleft", legend=names_ts[d1], bty="n", col="black", text.font=2, cex=0.65, x.intersp=0.5, y.intersp=1)
  }
  empty_plot()
  invisible(stvar)
}
