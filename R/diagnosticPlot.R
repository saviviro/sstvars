#' @import graphics
#'
#' @title Residual diagnostic plot for a STVAR model
#'
#' @description \code{diagnostic_plot} plots a multivariate residual diagnostic plot
#'   for either autocorrelation, conditional heteroskedasticity, or distribution,
#'   or simply draws the residual time series.
#'
#' @inheritParams get_boldA_eigens
#' @param type which type of diagnostic plot should be plotted?
#'   \itemize{
#'     \item{\code{"all"} all below sequentially.}
#'     \item{\code{"series"} the residual time series.}
#'     \item{\code{"ac"} the residual autocorrelation and cross-correlation functions.}
#'     \item{\code{"ch"} the squared residual autocorrelation and cross-correlation functions.}
#'     \item{\code{"dist"} the residual histogram with theoretical density (dashed line) and QQ-plots.}
#'   }
#' @param resid_type should standardized or raw residuals be used?
#' @param maxlag the maximum lag considered in types \code{"ac"} and \code{"ch"}.
#' @details Auto- and cross-correlations (types \code{"ac"} and \code{"ch"}) are calculated with the function
#'  \code{acf} from the package \code{stats} and the plot method for class \code{'acf'} objects is employed.
#'  If \code{cond_dist == "Student"}, the estimate of the degrees of freedom parameter is used in theoretical
#'  densities and quantiles.
#' @inherit get_residuals references
#' @seealso \code{\link{profile_logliks}}, \code{\link{fitSTVAR}}, \code{\link{STVAR}},
#'  \code{\link{LR_test}}, \code{\link{Wald_test}}, \code{\link{cond_moment_plot}}
#' @examples
#' # p=1, M=1, d=2, linear VAR
#' theta_112 <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103,
#'   0.601786, -0.002945, 0.067224)
#' mod112 <- STVAR(data=gdpdef, p=1, M=1, params=theta_112)
#' diagnostic_plot(mod112)
#'
#' # p=1, M=2, d=2, relative dens weight function
#' theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863,
#'   -0.3124, 0.505251, 0.298483, 0.030096, -0.176925, 0.838898, 0.310863, 0.007512,
#'   0.018244, 0.949533, -0.016941, 0.121403, 0.573269)
#' mod122 <- STVAR(data=gdpdef, p=1, M=2, params=theta_122relg)
#' diagnostic_plot(mod122, type="ac")
#' diagnostic_plot(mod122, type="ch")
#' @export

diagnostic_plot <- function(stvar, type=c("all", "series", "ac", "ch", "dist"), resid_type=c("standardized", "raw"),
                            maxlag=12) {

  if(is.null(stvar$data)) stop("The model needs to contain data!")
  type <- match.arg(type)
  resid_type <- match.arg(resid_type)
  if(resid_type == "standardized") {
    res <- stvar$residuals_std
  } else {
    res <- stvar$residuals_raw
  }
  d <- stvar$model$d
  names_ts <- colnames(as.ts(stvar$data))
  colnames(res) <- names_ts
  old_par <- par(no.readonly=TRUE)
  on.exit(par(old_par))
  if(type == "all") message("In total four figures of standardized residuals are plotted:
                            1) time series
                            2) autocorrelation function of residuals
                            3) autocorrelation function of squared residuals
                            4) histograms and QQ-plots")

  if(type == "series" || type == "all") {
    par(mfrow=c(d, 1), las=1)
    for(d1 in 1:d) {
      xaxt <- "n"
      if(d1 == 1) {
        par(mar=c(0.5, 2.5, 2.1, 1))
      } else if(d1 == d) {
        xaxt <- "s"
        par(mar=c(2.5, 2.5, 0.1, 1))
      } else {
        par(mar=c(0.5, 2.5, 0.1, 1))
      }
      yaxt1 <- round(min(res[,d1]))
      yaxt2 <- round(max(res[,d1]))
      main <- ifelse(all.equal(d1, 1), "Residual time series", "")
      plot(res[,d1], yaxt="n", xaxt=xaxt, type="l", col=grDevices::rgb(0, 0, 0, 1), ylab="", xlab="", main=main)
      axis(2, at=yaxt1:yaxt2, labels=yaxt1:yaxt2)
      abline(h=0, col=grDevices::rgb(1, 0, 0, 0.3), lwd=2)
      legend("topleft", legend=names_ts[d1], bty="n", col="black", text.font=2, cex=0.65, x.intersp=0.5, y.intersp=1)
    }
  }
  if(type == "ac" || type == "all") {
    par(mar=c(2.3, 2.8, 3.5, 1.0))
    acf(res, lag.max=maxlag, plot=TRUE, ylab="")
  }
  if(type == "ch" || type == "all") {
    par(mar=c(2.3, 2.8, 3.5, 1.0))
    acf(res^2, lag.max=maxlag, plot=TRUE, ylab="")
  }
  if(type == "dist" || type == "all") {
    par(mfrow=c(2, d), mar=c(2.5, 2.8, 2.1, 1.0))
    if(stvar$model$cond_dist == "Gaussian") {
      dens_fun <- function(y) dnorm(y)
      qqplot_fun <- function(y) qqnorm(y, main="", ylab="", xlab="")
      qqline_fun <- function(y) qqline(y, col="darkred")
    } else { # cond_dist == "Student"
      df <- stvar$params[length(stvar$params)] # The last param is always the df param
      dens_fun <- function(y) dt(y, df=df)
      qqplot_fun <- function(y){
        y <- sort(y, decreasing=FALSE) # Sorted sample quantiles
        T_obs <- length(y)
        p <- (1:T_obs - 0.5)/T_obs  # Probs; 0.5 substracted to avoid 0 and 1 that would result in -Inf and Inf
        t_quantiles <- qt(p, df=df)  # Theoretical quantiles; df taken from parent env
        plot(x=t_quantiles, y=y, main="", xlab="", ylab="") # Plot samples quantes agains theoretical quants
      }
      qqline_fun <- function(y) qqline(y, col="darkred", distribution=function(p) qt(p, df=df))
    }
    # Plot histograms with theoretical density
    for(i1 in 1:d) {
      hs <- hist(res[,i1], breaks="Scott", probability=TRUE, col="skyblue", plot=TRUE,
                 main=colnames(res)[i1], ylim=c(0, 0.5), ylab="", xlab="")
      x <- seq(from=min(hs$breaks), to=max(hs$breaks), length.out=1000)
      lines(x=x, y=dens_fun(x), lty=2, col="darkred", lwd=2)
    }
    # Plot QQ plots with theoretical quantiles
    for(i1 in 1:d) {
      qqplot_fun(res[,i1])
      qqline_fun(res[,i1])
    }
  }
}
