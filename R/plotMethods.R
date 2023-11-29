#' @import graphics
#' @describeIn STVAR plot method for class 'stvar'
#' @param x object of class \code{'stvar'}
#' @param ... currently not used.
#' @param plot_type should the series be plotted with the estimated transition weights
#'   or conditional means?
#' @details The plot displays the time series together with estimated transition weights.
#' @export

plot.stvar <- function(x, ..., plot_type=c("trans_weights", "cond_mean")) {
  stvar <- x
  stopifnot(!is.null(stvar$data))
  data <- as.ts(stvar$data)
  plot_type <- match.arg(plot_type)
  n_obs <- nrow(data)
  p <- stvar$model$p
  M <- stvar$model$M
  d <- ncol(data)
  params <- stvar$params
  ts_tw <- ts(rbind(matrix(NA, nrow=p, ncol=M), stvar$transition_weights),
              start=start(data), frequency=frequency(data)) # First p observations are starting values

  # Old graphical parameters, set back on exit.
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))

  # Time series and transition weights
  colpal_ts <- grDevices::colorRampPalette(c("darkgreen", "darkblue", "darkmagenta", "red3"))(d)
  colpal_tw <- grDevices::colorRampPalette(c("blue", "turquoise1", "green", "red"))(M)
  names_ts <- colnames(data)
  names_tw <- paste0("Regime", 1:M)

  if(plot_type == "trans_weights") {
    graphics::par(mfrow=c(2, 1), mar=c(2.5, 2.5, 2.1, 1))
    draw_legend <- function(nams, cols) {
      legend("topleft", legend=nams, bty="n", col=cols, lty=1, lwd=2, text.font=2, cex=0.6, x.intersp=0.5, y.intersp=1)
    }
    ts.plot(data, gpars=list(main="Time series", col=colpal_ts, lty=1:d))
    draw_legend(names_ts, cols=colpal_ts)
    ts.plot(ts_tw, gpars=list(main="Transition weights", ylim=c(0, 1), col=colpal_tw, lty=2))
    draw_legend(names_tw, cols=colpal_tw)
  } else { # plot_type == "cond_mean"
    if(is.null(stvar$regime_cmeans)) stop("Conditional means were not calculated when building this model")
    graphics::par(mfrow=c(d, 1), mar=c(0.5, 3, 2.1, 1), las=1)
    total_cmeans <- stvar$total_cmeans
    weight_x_reg <- lapply(1:d, function(d1) stvar$transition_weights*stvar$regime_cmeans[, d1, ])
    vals <- lapply(1:d, function(d1) c(total_cmeans[,d1], vec(weight_x_reg[[d1]]), data[,d1]))
    make_ts <- function(dat) ts(c(rep(NA, p), dat), start=start(data), frequency=frequency(data))

    for(d1 in 1:d) {
      xaxt <- "n"
      if(d1 == d) {
        xaxt <- "s"
        par(mar=c(2.5, 3, 0.5, 1))
      } else if(d1 > 1) {
        par(mar=c(0.5, 3, 0.5, 1))
      }
      ymin <- floor(min(vals[[d1]]))
      ymax <-  max(vals[[d1]]) # ceiling(max(vals[[d1]]))
      main <- ifelse(d1 == 1, "Conditional means", "")
      plot(data[,d1], ylim=c(ymin, ymax), xlab="", ylab="", xaxt=xaxt, main=main)
      lines(make_ts(total_cmeans[,d1]), col="grey", lty=2, lwd=2)
      for(m1 in 1:M) {
        lines(make_ts(weight_x_reg[[d1]][,m1]), col=colpal_tw[m1], lty=3)
      }
      legend("topleft", legend=names_ts[d1], bty="n", col="black", text.font=2, cex=0.65, x.intersp=0.5, y.intersp=1)
      if(d1 == 1) {
        legend("topright", legend=c("total", paste0("Regime ", 1:M)), bty="n", col=c("grey", colpal_tw),
               lty=c(2, rep(3, M)), lwd=2, text.font=2, cex=0.65, x.intersp=0.5, y.intersp=1)
      }
    }
  }
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


#' @describeIn GIRF plot method
#' @inheritParams print.girf
#' @param margs numeric vector of length four that adjusts the
#'  \code{[bottom_marginal, left_marginal, top_marginal, right_marginal]}
#'  as the relative sizes of the marginals to the figures of the responses.
#' @param ... graphical parameters passed to \code{plot} method plotting the GIRFs
#' @export

plot.girf <- function(x, margs, ...) {

  # Relevant statistics etc
  girf <- x
  girf_res <- girf$girf_res
  nresps <- ncol(girf_res[[1]]$point_est)
  resp_names <- colnames(girf_res[[1]]$point_est)
  ngirfs <- length(girf_res)

  # Graphical settings
  if(missing(margs)) {
    margs <- c(max(0.4, 0.4 + log(0.31 + log(nresps))/6),
               max(0.4, 0.4 + log(0.31 + log(ngirfs))),
               max(0.35, 0.35 + log(0.31 + log(nresps))/6),
               max(0.1, 0.1 + log(0.31 + log(ngirfs))/10))
    if(ngirfs == 1) margs[2] <- 0.3
    margs <- vapply(1:length(margs), function(i1) min(margs[i1], 1), numeric(1))
  } else {
    stopifnot(all(margs > 0))
  }
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))
  par(las=1, mar=c(0, 0, 0, 0))
  nrows <- nresps + 2 # + 2 for bottom and top marginals
  ncols <- 3*ngirfs # 3x for the left and right margin in each column of figures
  nfigs <- nrows*ncols
  layoutmat <- matrix(seq_len(nfigs), nrow=nrows, ncol=ncols, byrow=FALSE)
  layout(layoutmat, # Below -0.2 for not including the ylab and also adding to the right marginals
         widths=c(margs[2], 1, margs[4], rep(c(margs[2] - 0.2, 1, margs[4]), times=ngirfs - 1)),
         heights=c(margs[3], rep(1, times=nrows - 2), margs[1]))

  # Function to plot empty plots (for the marginals)
  empty_plot <- function() plot(0, xaxt='n', yaxt='n', bty='n', pch='', ylab='', xlab='')

  # Function to plot the GIRF for each response separately
  plot_girf <- function(resp_ind, main="", xaxt="n", ylab="", first_col=FALSE, ...) {

    # Plot point estimate
    point_est <- girf_i1$point_est[, resp_ind]
    conf_ints <- girf_i1$conf_ints[, , resp_ind]
    plot(x=0:(length(point_est) - 1), y=point_est, type="l", ylim=c(min(0, min(conf_ints)), max(0, max(conf_ints))),
         main="", ylab="", xlab="", xaxt=xaxt, lwd=2, col="blue", ...)
    if(first_col) { # Add yaxis label to the first column of responses
      mtext(resp_names[resp_ind], side=2, cex=0.8, font=2, las=0, padj=-4)
    }
    if(resp_ind == 1) mtext(main, padj=-0.5, cex=1, font=2)

    # Plot confidence intervals
    inds <- 0:girf$N
    draw_poly <- function(up_or_low) polygon(x=c(inds, rev(inds)), y=c(up_or_low, rev(point_est)),
                                             col=grDevices::rgb(0, 0, 1, 0.2), border=NA)

    for(i1 in 1:length(girf$ci)) {
      draw_poly(conf_ints[, i1]) # lower
      draw_poly(conf_ints[, ncol(conf_ints) + 1 - i1]) # upper
    }

    abline(h=0, lty=3, col="red")
  }

  # Loop through the shocks
  for(i1 in 1:ngirfs) {
    # Plot a column of empty plots as the left margins
    for(i2 in 1:(nrows + 1)) { # + 1 for the top margin of the first row of responses
      empty_plot()
    }

    # Plot the responses of each variable to shock i1
    girf_i1 <- girf_res[[i1]]

    # Plot the GIRF of shocks i1
    first_col <- i1 == 1
    plot_girf(resp_ind=1, main=paste("Shock", girf$shocks[i1]),
              ylab=resp_names[1], first_col=first_col)
    if(nresps > 2) {
      for(i2 in 2:(nresps - 1)) {
        plot_girf(resp_ind=i2, ylab=resp_names[i2], first_col=first_col)
      }
    }
    plot_girf(resp_ind=nresps, xaxt="s", ylab=resp_names[nresps], first_col=first_col)
    empty_plot() # To bottom margin of the last row of responses

    # Plot a column of empty plots as the right margins
    for(i2 in 1:nrows) {
      empty_plot()
    }
  }
}


#' @describeIn GFEVD plot method
#' @inheritParams print.gfevd
#' @param ... currently not used.
#' @export

plot.gfevd <- function(x, ...) {
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))

  gfevd <- x
  gfevd_res <- gfevd$gfevd_res
  n_gfevds <- dim(gfevd_res)[3]
  varnames <- dimnames(gfevd_res)[[3]]
  n_shocks <- dim(gfevd_res)[2]
  graphics::par(las=1, mfrow=c(1, 1), mar=c(2.6, 2.6, 2.6, 3.1))

  # Function to plot the GFEVD for each variable separately
  plot_gfevd <- function(var_ind, main) {

    one_gfevd <- gfevd_res[, , var_ind] # [horizon, shock]
    mycums <- as.matrix(1 - apply(one_gfevd[, 1:(ncol(one_gfevd) - 1), drop=FALSE], MARGIN=1, FUN=cumsum))
    if(ncol(mycums) > 1) mycums <- t(mycums)
    upper_ints <- cbind(rep(1, times=nrow(one_gfevd)), mycums)
    lower_ints <- cbind(mycums,rep(0, times=nrow(one_gfevd)))
    x_points <- seq(from=-0.5, to=gfevd$N + 0.5, by=1)

    # Plot the template
    colpal <- grDevices::adjustcolor(grDevices::topo.colors(ncol(one_gfevd)), alpha.f=0.3)
    xlim_adj <- ifelse(gfevd$N < 13, 0, 0.04*(gfevd$N - 12))
    plot(NA, ylim=c(0, 1), xlim=c(0 + xlim_adj, gfevd$N - xlim_adj),
         main=main, ylab="", xlab="")
    x_bars <- (0:(gfevd$N + 1) - 0.5)
    segments(x0=x_bars, y0=rep(0, times=length(x_bars)), x1=x_bars, y1=rep(1, times=length(x_bars)))

    # Go through the shocks
    for(i1 in 1:ncol(one_gfevd)) {
      for(i2 in 1:(length(x_points) - 1)) {
        polygon(x=c(x_points[c(i2, i2 + 1)], x_points[c(i2 + 1, i2)]),
                y=c(upper_ints[c(i2, i2), i1], lower_ints[c(i2, i2), i1]),
                col=colpal[i1], border=NA)
      }
    }

    # Add shock legengs
    ylim_ajd <- ifelse(n_shocks < 500, 0.02*n_shocks, 0.01*n_shocks)
    colpal2 <- grDevices::adjustcolor(colpal, alpha.f=3)
    xtext_adj <- ifelse(gfevd$N < 10, 0.75 - gfevd$N^(-1/3), 0.4 - gfevd$N/150)
    text(x=rep(gfevd$N + xtext_adj, n_shocks), y=seq(from=1, to=1 - 0.02*n_shocks, length.out=n_shocks),
         labels=paste("Shock", 1:n_shocks),
         col=colpal2, pos=4, font=2, cex=0.8, xpd=TRUE)

  }

  # Loop through the GFEVDs
  for(i1 in 1:n_gfevds) {
    if(i1 > 1) grDevices::devAskNewPage(TRUE)
    plot_gfevd(var_ind=i1, main=ifelse(i1 <= n_shocks,
                                       paste("GFEVD for ", varnames[i1]),
                                       paste("GFEVD for regime", i1 - n_shocks, "trans. weight")))

  }
}


