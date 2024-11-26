#' @title Estimate generalized impulse response function for
#'  structural STVAR models.
#'
#' @description \code{GIRF} estimates generalized impulse response function for
#'  structural STVAR models.
#'
#' @inheritParams iterate_more
#' @inheritParams simulate.stvar
#' @inheritParams fitSTVAR
#' @param which_shocks a numeric vector of length at most \eqn{d}
#'   (\code{=ncol(data)}) and elements in \eqn{1,...,d} specifying the
#'   structural shocks for which the GIRF should be estimated.
#' @param shock_size a non-zero scalar value specifying the common size for all scalar
#'   components of the structural shock. Note that the conditional covariance
#'   matrix of the structural shock is normalized to an identity matrix and that the
#'   (generalized) impulse responses may not be symmetric with respect to the sign
#'   and size of the shock.
#' @param N a positive integer specifying the horizon how far ahead should the
#'   generalized impulse responses be calculated.
#' @param R1 the number of repetitions used to estimate GIRF for each initial
#'   value.
#' @param R2 the number of initial values to use, i.e., to draw from \code{init_regime}
#'   if \code{init_values} are not specified. The confidence bounds
#'   will be sample quantiles of the GIRFs based on different initial values.
#'   Ignored if the argument \code{init_value} is specified.
#' @param which_cumulative a numeric vector with values in \eqn{1,...,d}
#'   (\code{d=ncol(data)}) specifying which the variables for which the impulse
#'   responses should be cumulative. Default is none.
#' @param scale should the GIRFs to some of the shocks be scaled so that they
#'   correspond to a specific magnitude of instantaneous or peak response
#'   of some specific variable (see the argument \code{scale_type})?
#'   Provide a length three vector where the shock of interest
#'   is given in the first element (an integer in \eqn{1,...,d}), the variable of
#'   interest is given in the second element (an integer in \eqn{1,...,d}), and
#'   the magnitude of its instantaneous or peak response in the third element
#'   (a non-zero real number). If the GIRFs of multiple shocks should be scaled, provide
#'   a matrix which has one column for each of the shocks with the columns being
#'   the length three vectors described above.
#' @param scale_type If argument \code{scale} is specified, should the GIRFs be
#'   scaled to match an instantaneous response (\code{"instant"}) or peak response
#'   (\code{"peak"}). If \code{"peak"}, the scale is based on the largest magnitude
#'   of peak response in absolute value. Ignored if \code{scale} is not specified.
#' @param scale_horizon If \code{scale_type == "peak"} what the maximum horizon up
#'   to which peak response is expected? Scaling won't based on values after this.
#' @param init_values a size \code{[p, d, R2]} array specifying the initial values in each slice
#'   for each Monte Carlo repetition, where d is the number of time series in the system and \code{R2}
#'   is an argument of this function. In each slice, the \strong{last} row will be used as initial values
#'   for the first lag, the second last row for second lag etc. If not specified, initial values will be
#'   drawn from the regime specified in \code{init_regimes}.
#' @param ci a numeric vector with elements in \eqn{(0, 1)} specifying the
#'   confidence levels of the confidence intervals.
#' @param ncores the number CPU cores to be used in parallel computing. Only
#'   single core computing is supported if an initial value is specified (and
#'   the GIRF won't thus be estimated multiple times).
#' @param exo_weights if \code{weight_function="exogenous"}, provide a size
#'  \eqn{(N+1 x M)} matrix of exogenous transition weights for the regimes: \code{[h, m]}
#'  for the (after-the-impact) period \eqn{h-1} and regime \eqn{m} weight (\code{[1, m]}
#'  is for the impact period). Ignored if \code{weight_function!="exogenous"}.
#' @param seeds a length \code{R2} vector containing the random number generator
#'   seed for estimation of each GIRF. A single number of an initial value is
#'   specified. or \code{NULL} for not initializing the seed. Exists for
#'   creating reproducible results.
#' @param use_parallel employ parallel computing? If \code{FALSE}, does not print
#'   anything.
#' @details The confidence bounds reflect uncertainty about the initial state (but
#'   not about the parameter estimates) if initial values are not
#'   specified. If initial values are specified, confidence intervals won't be
#'   estimated.
#'
#'   Note that if the argument \code{scale} is used, the scaled responses of
#'   the transition weights might be more than one in absolute value.
#'
#'   If \code{weight_function="exogenous"}, exogenous transition weights used in
#'   the Monte Carlo simulations for the future sample paths of the process must
#'   the given in the argument \code{exo_weights}. The same weights are used as
#'   the transition weights across the Monte Carlo repetitions.
#' @return Returns a class \code{'girf'} list with the GIRFs in the first
#'   element (\code{$girf_res}) and the used arguments the rest. The first
#'   element containing the GIRFs is a list with the \eqn{m}th element
#'   containing the point estimates for the GIRF in \code{$point_est} (the first
#'   element) and confidence intervals in \code{$conf_ints} (the second
#'   element). The first row is for the GIRF at impact \eqn{(n=0)}, the second
#'   for \eqn{n=1}, the third for \eqn{n=2}, and so on.
#'
#'   The element \code{$all_girfs} is a list containing results from all the individual GIRFs
#'   obtained from the MC repetitions. Each element is for one shock and results are in
#'   array of the form \code{[horizon, variables, MC-repetitions]}.
#' @seealso \code{\link{GFEVD}}, \code{\link{linear_IRF}}, \code{\link{fitSSTVAR}}
#'  \itemize{
#'    \item Kilian L., Lütkepohl H. 20017. Structural Vector Autoregressive Analysis. 1st edition.
#'      \emph{Cambridge University Press}, Cambridge.
#'  }

#' @examples
#'  \donttest{
#'  # These are long-running examples that use parallel computing.
#'  # It takes approximately 30 seconds to run all the below examples.
#'  # Note that larger R1 and R2 should be used for more reliable results;
#'  # small R1 and R2 are used here to shorten the estimation time.
#'
#'  # Recursively identified logistic Student's t STVAR(p=3, M=2) model with the first
#'  # lag of the second variable as the switching variable:
#'  params32logt <- c(0.5959, 0.0447, 2.6279, 0.2897, 0.2837, 0.0504, -0.2188, 0.4008,
#'   0.3128, 0.0271, -0.1194, 0.1559, -0.0972, 0.0082, -0.1118, 0.2391, 0.164, -0.0363,
#'   -1.073, 0.6759, 3e-04, 0.0069, 0.4271, 0.0533, -0.0498, 0.0355, -0.4686, 0.0812,
#'    0.3368, 0.0035, 0.0325, 1.2289, -0.047, 0.1666, 1.2067, 7.2392, 11.6091)
#'  mod32logt <- STVAR(gdpdef, p=3, M=2, params=params32logt, weight_function="logistic",
#'   weightfun_pars=c(2, 1), cond_dist="Student", identification="recursive")
#'
#'  # GIRF for one-standard-error positive structural shocks, N=30 steps ahead,
#'  # with the inital values drawn from the first regime.
#'  girf1 <- GIRF(mod32logt, which_shocks=1:2, shock_size=1, N=30, R1=50, R2=50,
#'   init_regime=2)
#'  print(girf1) # Print the results
#'  plot(girf1) # Plot the GIRFs
#'
#'  # GIRF for one-standard-error positive structural shocks, N=30 steps ahead,
#'  # with the inital values drawn from the second regime. The responses of the
#'  # GDP and GDP deflator growth rates are accumulated.
#'  girf2 <- GIRF(mod32logt, which_shocks=1:2, which_cumulative=1:2, shock_size=1,
#'   N=30, R1=50, R2=50, init_regime=2)
#'  plot(girf2) # Plot the GIRFs
#'
#'  # GIRF for two-standard-error negative structural shock - the first shock only.
#'  # N=50 steps ahead with the inital values drawn from the first regime. The responses
#'  # are scaled to correspond an instantanous increase of 0.5 of the first variable.
#'  girf3 <- GIRF(mod32logt, which_shocks=1, shock_size=-2, N=50, R1=50, R2=50,
#'   init_regime=1, scale_type="instant", scale=c(1, 1, 0.5))
#'  plot(girf3) # Plot the GIRFs
#'  }
#' @export

GIRF <- function(stvar, which_shocks, shock_size=1, N=30, R1=250, R2=250, init_regime=1, init_values=NULL,
                 which_cumulative=numeric(0), scale=NULL, scale_type=c("instant", "peak"), scale_horizon=N,
                 ci=c(0.95, 0.80), ncores=2, burn_in=1000, exo_weights=NULL, seeds=NULL, use_parallel=TRUE) {
  check_stvar(stvar)
  scale_type <- match.arg(scale_type)
  cond_dist <- stvar$model$cond_dist
  if(stvar$model$identification == "reduced_form" && cond_dist != "ind_Student" && cond_dist != "ind_skewed_t") {
    warning(paste("Reduced form model supplied, so using recursive identification"))
    stvar$model$identification <- "recursive"
  } else if(cond_dist == "ind_Student" || cond_dist == "ind_skewed_t") {
    stvar$model$identification <- "non-Gaussianity" # Readily identified by non-Gaussianity
  }

  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  identification <- stvar$model$identification
  stopifnot(N %% 1 == 0 && N > 0)
  stopifnot(scale_horizon %in% 0:N)
  if(!is.null(init_values)) {
    stopifnot(is.array(init_values) && all(dim(init_values) == c(p, d, R2)))
  }

  # Check the exogenous weights given for simulation
  if(stvar$model$weight_function == "exogenous") {
    check_exoweights(M=M, exo_weights=exo_weights, how_many_rows=N+1, name_of_row_number="N+1")
  }

  if(identification == "recursive") {
    B_constrs <- matrix(NA, nrow=d, ncol=d)
    B_constrs[upper.tri(B_constrs)] <- 0
    diag(B_constrs) <- 1 # Lower triangular Cholesky constraints
  } else if(identification == "heteroskedasticity" || identification == "non-Gaussianity") {
    if(is.null(stvar$model$B_constraints)) {
      B_constrs <- matrix(NA, nrow=d, ncol=d)
    } else {
      B_constrs <- stvar$model$B_constraints
    }
  }
  if(missing(which_shocks)) {
    which_shocks <- 1:d
  } else {
    stopifnot(all(which_shocks %in% 1:d))
    which_shocks <- unique(which_shocks)
  }
  if(!is.null(seeds) && length(seeds) != R2) stop("The argument 'seeds' needs be NULL or a vector of length 'R2'")
  stopifnot(init_regime %in% 1:M)
  stopifnot(length(ci) > 0 && all(ci > 0 & ci < 1))
  if(length(shock_size) != 1) {
    warning("The argument shock_size should be a numeric scalar. Using the first value.")
    shock_size <- shock_size[1]
  }
  stopifnot(shock_size != 0)
  if(!is.null(scale)) {
    scale <- as.matrix(scale)
    stopifnot(all(scale[1,] %in% 1:d)) # All shocks in 1,...,d
    stopifnot(length(unique(scale[1,])) == length(scale[1,])) # No duplicate scales for the same shock
    stopifnot(all(scale[2,] %in% 1:d)) # All variables in 1,...,d
    stopifnot(all(scale[3,] != 0)) # No zero initial magnitudes

    # For the considered shocks, check that there are not zero-constraints for the variable whose initial response is scaled.
    for(i1 in 1:ncol(scale)) {
      if(!is.na(B_constrs[scale[2, i1], scale[1, i1]]) && B_constrs[scale[2, i1], scale[1, i1]] == 0) {
        stop(paste("Instantaneous response of the variable that has a zero constraint for the considered shock cannot be scaled"))
      }
    }
  }
  if(length(which_cumulative) > 0) {
     which_cumulative <- unique(which_cumulative)
     stopifnot(all(which_cumulative %in% 1:d))
  }

  # Function that estimates GIRF
  get_one_girf <- function(shock_numb, shock_size, seed, rep_numb) {
    if(!is.null(init_values)) {
      init_v <- matrix(init_values[, , rep_numb], nrow=p, ncol=d)
    } else {
      init_v <- NULL
    }
    simulate.stvar(stvar, nsim=N+1, seed=seed, init_values=init_v, init_regime=init_regime,
                   ntimes=R1, burn_in=burn_in, exo_weights=exo_weights,
                   girf_pars=list(shock_numb=shock_numb, shock_size=shock_size))
  }

  GIRF_shocks <- vector("list", length=length(which_shocks)) # Storage for the GIRFs [[shock]]

  if(use_parallel) {
    if(ncores > parallel::detectCores()) {
      ncores <- parallel::detectCores()
      message("ncores was set to be larger than the number of cores detected")
    }
    if(is.null(init_values)) {
      message(paste("Using", ncores, "cores to estimate", R2,"GIRFs for", length(which_shocks), "structural shocks,",
                "each based on", R1, "Monte Carlo repetitions."))
    } else {
      message(paste("Using", ncores, "cores to estimate one GIRF for", length(which_shocks), "structural shocks, each based on",
              R1, "Monte Carlo repetitions."))
    }

    ### Calculate the GIRFs ###
    cl <- parallel::makeCluster(ncores)
    on.exit(try(parallel::stopCluster(cl), silent=TRUE)) # Close the cluster on exit, if not already closed.
    parallel::clusterExport(cl, ls(environment(GIRF)), envir=environment(GIRF)) # assign all variables from package:sstvars
    parallel::clusterEvalQ(cl, c(library(pbapply), library(Rcpp), library(RcppArmadillo), library(sstvars)))

    for(i1 in 1:length(which_shocks)) {
      message(paste0("Estimating GIRFs for structural shock ", which_shocks[i1], "..."))
      GIRF_shocks[[i1]] <- pbapply::pblapply(1:R2, function(i2) get_one_girf(shock_numb=which_shocks[i1],
                                                                             shock_size=shock_size, seed=seeds[i2],
                                                                             rep_numb=i2), cl=cl)
    }
    parallel::stopCluster(cl=cl)
  } else { # No parallel computing
    for(i1 in 1:length(which_shocks)) {
      GIRF_shocks[[i1]] <- lapply(1:R2, function(i2) get_one_girf(shock_numb=which_shocks[i1],
                                                                  shock_size=shock_size, seed=seeds[i2],
                                                                  rep_numb=i2))
    }
  }


  GIRF_results <- vector("list", length=length(which_shocks))
  all_GIRFS <- vector("list", length=length(which_shocks))
  names(all_GIRFS) <- paste0("shock", which_shocks)
  names(GIRF_results) <- paste0("shock", which_shocks)

  for(i1 in 1:length(which_shocks)) { # Go through shocks
    res_in_array <- array(unlist(GIRF_shocks[[i1]]), dim=c(N + 1, d + M, R2)) # [horz, vars, MCreps]
    if(length(which_cumulative) > 0) {
      for(i2 in which_cumulative) {
        res_in_array[, i2, ] <- apply(res_in_array[, i2, , drop=FALSE], MARGIN=3, FUN=cumsum) # Replace GIRF with cumulative GIRF
      }
    }

    # Scale the GIRFs if specified
    if(which_shocks[i1] %in% scale[1,]) { # GIRF of this shock should be scaled
      which_col <- which(which_shocks[i1] == scale[1,]) # which col of the scale-matrix contains the argument for this specific shock
      which_var <- scale[2, which_col] # According to initial/peak response of which variable the GIRFs should be scaled
      magnitude <- scale[3, which_col] # What should be the magnitude of the initial/peak response of this variable

      # To avoid potential problems with using == to compare numerical values
      my_comparison_fun <- function(vec1, scalar1) which(abs(vec1 - scalar1) < .Machine$double.eps)[1]
      for(i2 in 1:R2) { # Go through the MC repetitions
        # The scaling scalar is different for each MC repetition, because the instantaneous/peak movement is generally
        # different with different starting values.
        if(scale_type == "instant") {  # Scale by initial response
          one_scale <- magnitude/res_in_array[1, which_var, i2]
        } else {  # scale_type == "peak", "peak_max" or "peak_min", scale by peak response
          inds <- 1:(scale_horizon + 1) # +1 for period 0
          one_scale <- magnitude/res_in_array[my_comparison_fun(vec1=abs(res_in_array[inds, which_var, i2]),
                                                                scalar1=max(abs(res_in_array[inds, which_var, i2]))), which_var, i2]
          #one_scale <- magnitude/res_in_array[which(abs(res_in_array[, which_var, i2]) ==
          #             max(abs(res_in_array[, which_var, i2]))), which_var, i2]
        }
        res_in_array[, , i2] <- one_scale*res_in_array[, , i2]
      }
    }
    all_GIRFS[[i1]] <- res_in_array

    # Point estimates, confidence intervals
    colnames(res_in_array) <- colnames(GIRF_shocks[[1]][[1]])
    point_estimate <- apply(X=res_in_array, MARGIN=1:2, FUN=mean)
    lower <- (1 - ci)/2
    upper <- rev(1 - lower)
    q_tocalc <- c(lower, upper)
    q_tocalc <- sort(q_tocalc, decreasing=FALSE)
    conf_ints <- apply(res_in_array, MARGIN=1:2, FUN=quantile, probs=q_tocalc)
    conf_ints <- aperm(conf_ints, perm=c(2, 1, 3))
    rownames(point_estimate) <- rownames(conf_ints) <- 0:N
    GIRF_results[[i1]] <- list(point_est=point_estimate,
                               conf_ints=conf_ints)
  }

  if(use_parallel) message("Finished!")
  structure(list(girf_res=GIRF_results,
                 all_girfs=all_GIRFS,
                 shocks=which_shocks,
                 shock_size=shock_size,
                 N=N,
                 R1=R1,
                 R2=R2,
                 ci=ci,
                 which_cumulative=which_cumulative,
                 init_regime=init_regime,
                 init_values=init_values,
                 seeds=seeds,
                 burn_in=burn_in,
                 exo_weights=exo_weights,
                 stvar=stvar),
            class="girf")
}



#' @title Estimate generalized forecast error variance decomposition for structural
#'  STVAR models.
#'
#' @description \code{GFEVD} estimates generalized forecast error variance decomposition
#'  for structural STVAR models.
#'
#' @inheritParams GIRF
#' @param shock_size What sign and size should be used for all shocks? By the normalization, the conditional
#'   covariance matrix of the structural error is an identity matrix.
#' @param N a positive integer specifying the horizon how far ahead should the GFEVD be calculated.
#' @param initval_type What type initial values are used for estimating the GIRFs that the GFEVD is based on?
#'   \describe{
#'     \item{\code{"data"}:}{Estimate the GIRF for all the possible length \eqn{p} histories in the data.}
#'     \item{\code{"random"}:}{Estimate the GIRF for several random initial values generated from the a specific regime
#'        specified by the argument \code{init_regimes}. The number of initial values is set with the argument \code{R2}.}
#'     \item{\code{"fixed"}:}{Estimate the GIRF for the initial values specified with the argument \code{init_values}.}
#'   }
#' @param use_data_shocks \code{TRUE} for a special feature in which for every possible length \eqn{p} history in the data,
#'   the GFEVD is estimated for a shock that has the sign and size of the corresponding structural shock recovered from the data.
#'   See the details section.
#' @param R2 the number of initial values to be drawn/used if \code{initval_type="random"} or \code{"fixed"}.
#' @param seeds a numeric vector containing the random number generator seed for estimation
#'   of each GIRF. Should have the length...
#'   \itemize{
#'     \item ...\code{nrow(data) - p + 1} if \code{initval_type="data"}.
#'     \item ...\code{R2} if \code{initval_type="random"}.
#'     \item ...\code{1} if \code{initval_type="fixed."}.
#'   }
#'   Set to \code{NULL} for not initializing the seed. Exists for creating reproducible results.
#' @details  The GFEVD is a forecast error variance decomposition calculated with the generalized
#'   impulse response function (GIRF). See Lanne and Nyberg (2016) for details.
#'
#'   If \code{use_data_shocks == TRUE}, the GIRF is estimated for a shock that has the sign and size of the
#'   corresponding structural shock recovered from the fitted model. This is done for every possible length \eqn{p} history
#'   in the data. The GFEVD is then calculated as the average of the GFEVDs obtained from the GIRFs estimated for
#'   the data shocks. The plot and print methods can be used as usual for this GFEVD. However, this feature also
#'   obtain the contribution of each shock to the variance of the forecast errors at various horizons in specific
#'   historical points of time. This can be done by using the plot method with the argument \code{data_shock_pars}.
#'   Note that the arguments \code{shock_size}, \code{initval_type}, and \code{init_regime} are ignored if
#'   \code{use_data_shocks == TRUE}.
#' @return Returns and object of class 'gfevd' containing the GFEVD for all the variables and to
#'   the transition weights. Note that the decomposition does not exist at horizon zero for transition weights
#'   because the related GIRFs are always zero at impact.
#'   If \code{use_data_shocks=TRUE}, also contains the GFEVDs for each length \eqn{p} history in the data as
#'   4D array with dimensions \code{[horizon, variable, shock, time]}.
#' @seealso \code{\link{GIRF}}, \code{\link{linear_IRF}}, \code{\link{fitSSTVAR}}
#' @references
#'  \itemize{
#'    \item Lanne M. and Nyberg H. 2016. Generalized Forecast Error Variance Decomposition for Linear
#'      and Nonlineae Multivariate Models. \emph{Oxford Bulletin of Economics and Statistics}, \strong{78}, 4, 595-603.
#'  }
#' @examples
#'  \donttest{
#'  # These are long-running examples that use parallel computing.
#'  # It takes approximately 30 seconds to run all the below examples.
#'  # Note that larger R1 and R2 should be used for more reliable results;
#'  # small R1 and R2 are used here to shorten the estimation time.
#'
#'  # Recursively identifed logistic Student's t STVAR(p=3, M=2) model with the first
#'  # lag of the second variable as the switching variable:
#'  params32logt <- c(0.5959, 0.0447, 2.6279, 0.2897, 0.2837, 0.0504, -0.2188, 0.4008,
#'   0.3128, 0.0271, -0.1194, 0.1559, -0.0972, 0.0082, -0.1118, 0.2391, 0.164, -0.0363,
#'   -1.073, 0.6759, 3e-04, 0.0069, 0.4271, 0.0533, -0.0498, 0.0355, -0.4686, 0.0812,
#'    0.3368, 0.0035, 0.0325, 1.2289, -0.047, 0.1666, 1.2067, 7.2392, 11.6091)
#'  mod32logt <- STVAR(gdpdef, p=3, M=2, params=params32logt, weight_function="logistic",
#'   weightfun_pars=c(2, 1), cond_dist="Student", identification="recursive")
#'
#'  # GFEVD for one-standard-error positive structural shocks, N=30 steps ahead,
#'  # with fix initial values assuming all possible histories in the data.
#'  gfevd1 <- GFEVD(mod32logt, shock_size=1, N=30, initval_type="data", R1=10,
#'    seeds=1:(nrow(mod32logt$data)-2))
#'  print(gfevd1) # Print the results
#'  plot(gfevd1) # Plot the GFEVD
#'
#'  # GFEVD for one-standard-error positive structural shocks, N=30 steps ahead,
#'  # with fix initial values that are the last p observations of the data.
#'  gfevd2 <- GFEVD(mod32logt, shock_size=1, N=30, initval_type="fixed", R1=100, R2=1,
#'   init_values=array(mod32logt$data[(nrow(mod32logt$data) - 2):nrow(mod32logt$data),],
#'   dim=c(3, 2, 1)), seeds=1)
#'  plot(gfevd2) # Plot the GFEVD
#'
#'  # GFEVD for two-standard-error negative structural shocks, N=50 steps ahead
#'  # with the inital values drawn from the first regime. The responses of both
#'  # variables are accumulated.
#'  gfevd3 <- GFEVD(mod32logt, shock_size=-2, N=50, initval_type="random",
#'   R1=50, R2=50, init_regime=1)
#'  plot(gfevd3) # Plot the GFEVD
#'
#'  # GFEVD calculated for each lenght p history in the data in such a way that
#'  # for each history, the structural shock recoved from the fitted model is
#'  # used.
#'  gfevd4 <- GFEVD(mod32logt, N=20, use_data_shocks=TRUE, R1=10)
#'  plot(gfevd4) # Usual plot method
#'
#'  # Plot the contribution of the first to the variance of the forecast errors at
#'  # the historial points of time using the structural shocks recovered from the data:
#'  plot(gfevd4, data_shock_pars=c(1, 0)) # Contribution at impact
#'  plot(gfevd4, data_shock_pars=c(1, 2)) # Contribution after two periods
#'  plot(gfevd4, data_shock_pars=c(1, 4)) # Contribution after four periods
#'  }
#' @export

GFEVD <- function(stvar, shock_size=1, N=30, initval_type=c("data", "random", "fixed"), use_data_shocks=FALSE,
                  R1=250, R2=250, init_regime=1, init_values=NULL, which_cumulative=numeric(0), ncores=2,
                  burn_in=1000, exo_weights=NULL, seeds=NULL, use_parallel=TRUE) {
  check_stvar(stvar)
  initval_type <- match.arg(initval_type)
  cond_dist <- stvar$model$cond_dist
  if(stvar$model$identification == "reduced_form" && cond_dist != "ind_Student" && cond_dist != "ind_skewed_t") {
    warning(paste("Reduced form model supplied, so using recursive identification"))
    stvar$model$identification <- "recursive"
  } else if(cond_dist == "ind_Student" || cond_dist == "ind_skewed_t") {
    stvar$model$identification <- "non-Gaussianity" # Readily identified by non-Gaussianity
  }
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  if(!is.null(init_values)) {
    stopifnot(is.array(init_values) && all(dim(init_values) == c(p, d, R2)))
  }

  # Check the exogenous weights given for simulation
  if(stvar$model$weight_function == "exogenous") {
    check_exoweights(M=M, exo_weights=exo_weights, how_many_rows=N+1, name_of_row_number="N+1")
  }

  if(use_data_shocks) {
    initval_type <- "data" # Initval_type always data with data shocks
  }
  if(initval_type == "data") {
    if(is.null(stvar$data)) {
      stop("The model does not contain data! Add data with the function 'add_data' or select another 'initval_type'.")
    }
    stopifnot(nrow(stvar$data) >= p)
    # The number of length p histories in the data: the last history is not used if data shocks are used,
    # because the last history is used a by the time index T+1 for which the shock cannot be recovered.
    R2 <- nrow(stvar$data) - p + ifelse(use_data_shocks, 0, 1)
    all_initvals <- array(vapply(1:R2, function(i1) stvar$data[i1:(i1 + p - 1),], numeric(p*d)),
                          dim=c(p, d, R2)) # [, , i1] for i1 initval
  } else if(initval_type == "random") {
    stopifnot(init_regime%in% 1:M)
    all_initvals <- array(dim=c(1, 1, R2)) # This won't be used. NULL init_values will be used.
  } else if(initval_type == "fixed") {
    if(is.null(init_values)) stop(paste0("Initial values were not specified! Specify the initial values with the argument",
                                         "'init_values' or choose another 'initval_type'."))
    if(anyNA(init_values)) stop("init_values contains NA values")
    all_initvals <- init_values
  }
  if(!is.null(seeds) && length(seeds) != R2) stop("The argument 'seeds' has wrong length!")
  stopifnot(length(shock_size) == 1)
  if(length(which_cumulative) > 0) {
    which_cumulative <- unique(which_cumulative)
    stopifnot(all(which_cumulative %in% 1:d))
  }

  # Calculate the shock sizes for each history in the data
  if(use_data_shocks) {
    # Recover the structural shocks for each initial value in all_initvals:
    if(stvar$model$identification == "reduced_form" && cond_dist != "ind_Student" && cond_dist != "ind_skewed_t") {
      # Recover the structural shocks from the reduced form shocks using recursive identification:
      data_shocks <- get_residuals(data=stvar$data, p=p, M=M, params=stvar$params, weight_function=stvar$model$weight_function,
                                   weightfun_pars=stvar$model$weightfun_pars, cond_dist=stvar$model$cond_dist,
                                   parametrization=stvar$model$parametrization, identification="recursive",
                                   B_constraints=stvar$model$B_constraints, mean_constraints=stvar$model$mean_constraints,
                                   AR_constraints=stvar$model$AR_constraints, weight_constraints=stvar$model$weight_constraints,
                                   penalized=stvar$penalized, penalty_params=stvar$penalty_params,
                                   allow_unstab=stvar$allow_unstab, structural_shocks=TRUE)
    } else {
      data_shocks <- stvar$structural_shocks
    }
  }

  # Function that estimates GIRF
  get_one_girf <- function(shock_numb, shock_size, seed, init_values_for_1girf) {
    if(initval_type == "random") init_values_for_1girf <- NULL
    simulate.stvar(stvar, nsim=N+1, init_values=init_values_for_1girf, ntimes=R1, seed=seed, init_regime=init_regime,
                   burn_in=burn_in, exo_weights=exo_weights, girf_pars=list(shock_numb=shock_numb, shock_size=shock_size))
  }

  GIRF_shocks <- vector("list", length=d) # Storage for the GIRFs [[shock]]

  if(use_parallel) {
    if(ncores > parallel::detectCores()) {
      ncores <- parallel::detectCores()
      message("ncores was set to be larger than the number of cores detected")
    }
    if(initval_type != "fixed") {
      message(paste("Using", ncores, "cores to estimate", R2,"GIRFs for", d, "structural shocks,", "each based on",
                R1, "Monte Carlo repetitions."))
    } else {
      message(paste("Using", ncores, "cores to estimate one GIRF for", d, "structural shocks, each based on",
                R1, "Monte Carlo repetitions."))
    }

    ### Calculate the GIRFs ###
    cl <- parallel::makeCluster(ncores)
    on.exit(try(parallel::stopCluster(cl), silent=TRUE)) # Close the cluster on exit, if not already closed.
    parallel::clusterExport(cl, ls(environment(GFEVD)), envir = environment(GFEVD)) # assign all variables from package:gmvarkit
    parallel::clusterEvalQ(cl, c(library(pbapply), library(Rcpp), library(RcppArmadillo), library(sstvars)))

    for(i1 in 1:d) {
      message(paste0("Estimating GIRFs for structural shock ", i1, "..."))
      GIRF_shocks[[i1]] <- pbapply::pblapply(1:R2,
                                             function(i2) get_one_girf(shock_numb=i1,
                                                                       shock_size=ifelse(use_data_shocks, data_shocks[i2, i1], shock_size),
                                                                       seed=seeds[i2],
                                                                       init_values_for_1girf=matrix(all_initvals[, , i2],
                                                                                                    nrow=p, ncol=d)), cl=cl)
    }
    parallel::stopCluster(cl=cl)
  } else { # No parallel computing
    for(i1 in 1:d) {
      GIRF_shocks[[i1]] <- lapply(1:R2, function(i2) get_one_girf(shock_numb=i1,
                                                                  shock_size=ifelse(use_data_shocks, data_shocks[i2, i1], shock_size),
                                                                  seed=seeds[i2],
                                                                  init_values_for_1girf=matrix(all_initvals[, , i2], nrow=p, ncol=d)))
    }
  }

  if(is.null(colnames(stvar$data))) {
    varnames <- paste0("Variable", 1:d)
  } else {
    varnames <- colnames(stvar$data)
  }
  varnames <- c(varnames, paste0("tw", 1:M))
  shocknames <- paste0("Shock", 1:d)

  GIRF_square_cumsum <- array(dim=c(N + 1, length(varnames), d),
                              dimnames=list(0:N, varnames, shocknames)) # [horizon, variable, shock]
  for(i1 in 1:d) { # Go through the shocks
    res_in_array <- array(unlist(GIRF_shocks[[i1]]), dim=c(N + 1, length(varnames), R2))
    if(length(which_cumulative) > 0) {
      for(i2 in which_cumulative) {
        res_in_array[, i2, ] <- apply(res_in_array[, i2, , drop=FALSE], MARGIN=3, FUN=cumsum) # Replace GIRF with cumulative GIRF
      }
    }
    GIRF_square_cumsum[, , i1] <-  apply(apply(X=res_in_array, MARGIN=1:2, FUN=mean)^2, MARGIN=2, FUN=cumsum)
  }

  GFEVD_results <- array(dim=c(N + 1, d, length(varnames)),
                         dimnames=list(0:N, shocknames, varnames)) # [horizon, shock, variable]
  for(i1 in 1:ncol(GIRF_square_cumsum)) { # Go through the variables
    denominator <- rowSums(GIRF_square_cumsum[, i1, ]) # The denominators for h=0,1,...,N
    for(i2 in 1:d) { # Go through the shocks
      GFEVD_results[ , i2, i1] <- GIRF_square_cumsum[, i1, i2]/denominator
    }
  }

  # Is data shocks are used, calculate the GFEVD also separately for each length p history in the data
  if(use_data_shocks) {
    # Note that the dim is different to the other GFEVD_results below, because the interest is in the contribution of
    # each shock to the variance of the forecast errors at various horizons in specific historical points of time.
    data_GFEVD_results <- array(dim=c(N + 1, length(varnames), d, R2),
                                dimnames=list(0:N, varnames, shocknames, 1:R2)) # [horizon, variable, shock, time])
    for(i1 in 1:R2) { # Go through the length p histories
      GIRF_square_cumsum_i1 <- array(NA, dim=c(N+1, length(varnames), d)) # [horizon, variable, shock] for time period i1
      for(i2 in 1:d) { # Go through the shocks
        res_in_matrix <- GIRF_shocks[[i2]][[i1]] # [horizon, variable] GIRF of shock i2 at the time i1
        if(length(which_cumulative > 0)) { # Any GIRFs that should be accumulated?
          for(i3 in which_cumulative) {
            res_in_matrix[, i3] <- cumsum(res_in_matrix[, i3]) # Replace GIRF with cumulative GIRF
          }
        }
        # Squared cumulative GIRF of shock i2 at the time i1
        GIRF_square_cumsum_i1[, , i2] <- apply(res_in_matrix^2, MARGIN=2, FUN=cumsum) # [horizon, variable]
      }
      # Calculate the denonimators
      for(i2 in 1:length(varnames)) { # Go through the variables
        denominator <- rowSums(GIRF_square_cumsum_i1[, i2, ]) # The denominators for h=0,1,...,N (sums over socks)
        for(i3 in 1:d) { # Go through the shocks
          data_GFEVD_results[ , i2, i3, i1] <- GIRF_square_cumsum_i1[, i2, i3]/denominator
        }
      }
    }
  } else {
    data_GFEVD_results <- NULL
  }

  if(use_parallel) message("Finished!")
  structure(list(gfevd_res=GFEVD_results,
                 shock_size=shock_size,
                 N=N,
                 initval_type=initval_type,
                 R1=R1,
                 R2=R2,
                 which_cumulative=which_cumulative,
                 init_regime=init_regime,
                 init_values=init_values,
                 seeds=seeds,
                 burn_in=burn_in,
                 data_gfevd_res=data_GFEVD_results,
                 stvar=stvar),
            class="gfevd")
}
