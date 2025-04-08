#' @title Plot profile log-likelihood functions about the estimates
#'
#' @description \code{profile_logliks} plots profile log-likelihood functions about the estimates.
#'
#' @inheritParams loglikelihood
#' @inheritParams iterate_more
#' @param which_pars the profile log-likelihood function of which parameters should be plotted? An integer
#'  vector specifying the positions of the parameters in the parameter vector. The parameter vector has the
#'  form...
#' @param scale a numeric scalar specifying the interval plotted for each estimate:
#'  the estimate plus-minus \code{abs(scale*estimate)}.
#' @param nrows how many rows should be in the plot-matrix? The default is \code{max(ceiling(log2(length(which_pars)) - 1), 1)}.
#' @param ncols how many columns should be in the plot-matrix? The default is \code{ceiling(length(which_pars)/nrows)}.
#'   Note that \code{nrows*ncols} should not be smaller than the length of \code{which_pars}.
#' @param precision at how many points should each profile log-likelihood function be evaluated at?
#' @details When the number of parameters is large, it might be better to plot a smaller number of profile
#'  log-likelihood functions at a time using the argument \code{which_pars}.
#'
#' The red vertical line points the estimate.
#' @return  Only plots to a graphical device and doesn't return anything.
#' @inherit loglikelihood references
#' @seealso \code{\link{get_foc}}, \code{\link{get_soc}}, \code{\link{diagnostic_plot}}
#' @examples
#' # Threshold STVAR with p=1, M=2, the first lag of the second variable as switching variable:
#' pars <- c(0.5231, 0.1015, 1.9471, 0.3253, 0.3476, 0.0649, -0.035, 0.7513, 0.1651,
#'  -0.029, -0.7947, 0.7925, 0.4233, 5e-04, 0.0439, 1.2332, -0.0402, 0.1481, 1.2036)
#' mod12thres <- STVAR(data=gdpdef, p=1, M=2, params=pars, weight_function="threshold",
#'   weightfun_pars=c(2, 1))
#'
#' # Plot the profile log-likelihood functions of all parameters:
#' profile_logliks(mod12thres, precision=50) # Plots fast with precision=50
#'
#' # Plot only the profile log-likelihood function of the threshold parameter
#' # (which is the last parameter in the parameter vector):
#' profile_logliks(mod12thres, which_pars=length(pars), precision=100)
#'
#' # Plot only the profile log-likelihood functions of the intercept parameters
#' # (which are the first four parameters in the parameter vector, as d=2 and M=2):
#' profile_logliks(mod12thres, which_pars=1:4, precision=100)
#' @export

profile_logliks <- function(stvar, which_pars, scale=0.1, nrows, ncols, precision=50,
                            stab_tol=0.001, posdef_tol=1e-08, distpar_tol=1e-08, weightpar_tol=1e-08) {
  # Initial checks
  check_stvar(stvar)
  if(is.null(stvar$data)) stop("Cannot plot profile logliks if the model does not contain data.")

  # Model specs
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  params <- stvar$params
  weight_function <- stvar$model$weight_function
  weightfun_pars <- check_weightfun_pars(data=stvar$dat, p=p, M=M, d=d, weight_function=weight_function,
                                         weightfun_pars=stvar$model$weightfun_pars)
  cond_dist <- stvar$model$cond_dist
  parametrization <- stvar$model$parametrization
  identification <- stvar$model$identification
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  weight_constraints <- stvar$model$weight_constraints
  B_constraints <- stvar$model$B_constraints
  penalized <- stvar$penalized
  penalty_params <- stvar$penalty_params
  allow_unstab <- stvar$allow_unstab

  # Checks, default arguments etc
  if(missing(which_pars)) which_pars <- 1:length(params)
  if(!all_pos_ints(which_pars)) {
    stop("The argument 'which_pars' should contain strictly positive integers.")
  }
  if(any(which_pars > length(params))) {
    warning("Some elements in which_pars are larger than the number of parameters.")
    which_pars <- which_pars[which_pars <= length(params)]
    if(length(which_pars) == 0) {
      stop("All elements in which_pars are larger than the number of parameters.")
    }
  }
  if(anyDuplicated(which_pars) != 0) {
    warning("There are dublicates in which_pars")
    which_pars <- unique(which_pars)
  }
  npars <- length(which_pars)
  if(missing(nrows)) nrows <- max(ceiling(log2(npars) - 1), 1)
  if(missing(ncols)) ncols <- ceiling(npars/nrows)

  # Graphical settings: restore on exit.
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(mar=c(2.1, 2.1, 1.6, 1.1), mfrow=c(nrows, ncols))

  # Determine the numbers of each type of parameters:
  if(is.null(mean_constraints)) {
    n_mean_pars <- M*d
    g <- M # Number groups with the same mean parameters
  } else { # Means constrained
    n_mean_pars <- d*length(mean_constraints)
    g <- length(mean_constraints) # Number groups with the same mean parameters
  }
  if(is.null(AR_constraints)) {
    n_ar_pars <- M*p*d^2
  } else { # AR matrices constrained
    n_ar_pars <- ncol(AR_constraints)
  }
  if(cond_dist %in% c("ind_Student", "ind_skewed_t") || identification == "non-Gaussianity") {
    if(is.null(B_constraints)) {
      n_zeros <- 0
    } else {
      n_zeros <- M*sum(B_constraints == 0, na.rm=TRUE) # Zeros for all B_m
    }
    n_covmat_pars <- M*d^2 - n_zeros
  } else if(identification %in% c("reduced_form", "recursive")) {
    n_covmat_pars <- M*d*(d + 1)/2
  } else if(identification == "heteroskedasticity") {
    if(is.null(B_constraints)) {
      n_zeros <- 0
    } else {
      n_zeros <- sum(B_constraints == 0, na.rm=TRUE)
    }
    n_covmat_pars <- d^2 + d*(M - 1) - n_zeros
    W_row_ind <- rep(1, times=d) # row for each column
  }
  if(weight_function == "exogenous") {
    n_weight_pars <- 0
  } else {
    if(is.null(weight_constraints)) {
      if(weight_function == "relative_dens" || weight_function == "threshold") {
        n_weight_pars <- M - 1
      } else if(weight_function == "logistic" || weight_function == "exponential") {
        n_weight_pars <- 2
      } else if(weight_function == "mlogit") {
        n_weight_pars <- (M - 1)*(1 + length(weightfun_pars[[1]])*weightfun_pars[[2]])
      }
    } else { # Constraints on the weight parameters
      if(all(weight_constraints[[1]] == 0)) {
        n_weight_pars <- 0 # alpha = r, not in the parameter vector
      } else {
        n_weight_pars <- ncol(weight_constraints[[1]]) # The dimension of xi
      }
    }
  }

  if(cond_dist == "Gaussian") {
    n_dist_pars <- 0
  } else if(cond_dist == "Student") {
    n_dist_pars <- 1 # degrees of freedom param
  } else if(cond_dist == "ind_Student") {
    n_dist_pars <- d # df params
  } else { # cond_dist = "ind_skewed_t"
    n_dist_pars <- 2*d # df and skewness params
  }

  B_row_ind <- NULL # Initialize

  # Go though the parameters in which_pars
  for(i1 in which_pars) {
    pars <- params
    range <- abs(scale*pars[i1])
    vals <- seq(from=pars[i1] - range, to=pars[i1] + range, length.out=precision) # Loglik to be evaluated at these points
    logliks <- vapply(vals, function(val) { # Log-likelihoods about the estimate
      new_pars <- pars
      new_pars[i1] <- val # Change the single parameter value
      loglikelihood(data=stvar$data, p=p, M=M, params=new_pars,
                    weight_function=weight_function, weightfun_pars=weightfun_pars,
                    cond_dist=cond_dist, parametrization=parametrization,
                    identification=identification, AR_constraints=AR_constraints,
                    mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                    B_constraints=B_constraints, to_return="loglik", check_params=TRUE,
                    minval=NA, penalized=penalized, penalty_params=penalty_params,
                    allow_unstab=allow_unstab, stab_tol=stab_tol, posdef_tol=posdef_tol,
                    distpar_tol=distpar_tol, weightpar_tol=weightpar_tol)
    }, numeric(1))

    # Determine which type of parameter is i1 to determine the label for the individual plot
    if(i1 <= n_mean_pars) { # mean or intercept parameter
      cum_d <- (0:(g - 1))*d # The indeces after which regime changes
      m <- sum(i1 > cum_d) # Which regime?
      pos <- i1 - cum_d[m] # Which time series?
      if(parametrization == "intercept") {
        main <- substitute(phi[foo](foo2), list(foo=paste0(m, ",0"), foo2=pos))
      } else {
        main <- substitute(mu[foo](foo2), list(foo=m, foo2=pos))
      }
    } else if(i1 <= n_mean_pars + n_ar_pars) { # AR parameter
      if(is.null(AR_constraints)) {
        cum_q <- n_mean_pars + (0:(M - 1))*p*d^2 # The indeces after which regime changes
        m <- sum(i1 > cum_q) # To which regime the parameter is related to
        cum_lag <- cum_q[m] + (0:(p - 1))*d^2 # The indeces after which lag changes in the current regime
        j <- sum(i1 > cum_lag)  # To which lag the parameter is related to
        # Next, we want to obtain the row and column A_mj, the index is in vec(A_mj)
        cum_col <- cum_lag[j] + (0:(d - 1))*d # The index after which column changes
        col_ind <- sum(i1 > cum_col) # Column index of A_mj
        row_ind <- sum(i1 > cum_col[col_ind] + 0:(d - 1)) # Row index of A_mi
        main <- substitute(A[foo](foo2), list(foo=paste0(m, ",", j), foo2=paste0(row_ind, ",", col_ind)))
      } else {
        main <- substitute(psi(foo), list(foo=i1 - n_mean_pars))
      }
    } else if(i1 <= n_mean_pars + n_ar_pars + n_covmat_pars) { # Covariance matrix parameter
      if(identification %in% c("reduced_form", "recursive") && cond_dist != "ind_Student" && cond_dist != "ind_skewed_t") {
        cum_o <- n_mean_pars + n_ar_pars + (0:(M - 1))*d*(d + 1)/2 # The indeces after which regime changes
        m <- sum(i1 > cum_o) # Which regime
        cum_col <- cum_o[m] + c(0, cumsum(d - 0:(d - 1))) # The indeces after which column changes in vech(Omega_m)
        col_ind <- sum(i1 > cum_col) # Column in Omega_m
        row_inds <- unlist(lapply(1:d, function(i2) i2:d)) # At which row are we in for each pos in vech(Omega_m)
        row_ind <- row_inds[i1 - cum_o[m]] # Row in Omega_m
        main <- substitute(Omega[foo](foo2), list(foo=m, foo2=paste0(row_ind, ",", col_ind)))
      } else { # Structural model with different params than reduced form models, or cond_dist=ind_Student
        if(identification == "heteroskedasticity") {
          if(i1 <= n_mean_pars + n_ar_pars + d^2 - n_zeros) { # W params
            n_zeros_in_each_column <- vapply(1:d, function(i2) sum(B_constraints[,i2] == 0, na.rm=TRUE), numeric(1))
            zero_positions <- lapply(1:d, function(i2) (1:d)[B_constraints[,i2] == 0 & !is.na(B_constraints[,i2])]) # 0 constr pos in each col
            cum_wc <- c(0, cumsum(d - n_zeros_in_each_column)) # Index in W parameters after which a new column in W starts
            posw <- i1 - (n_mean_pars + n_ar_pars) # Index in W parameters
            col_ind <- sum(posw > cum_wc)
            while(TRUE) {
              if(W_row_ind[col_ind] %in% zero_positions[[col_ind]]) {
                W_row_ind[col_ind] <- W_row_ind[col_ind] + 1
              } else {
                break
              }
            }
            main <- substitute(W(foo), list(foo=paste0(W_row_ind[col_ind], ", ", col_ind)))
            W_row_ind[col_ind] <- W_row_ind[col_ind] + 1
          } else { # lambda parameters
            cum_lamb <- n_mean_pars + n_ar_pars + d^2 - n_zeros + c(0, cumsum(rep(d, times=M))) # Index after which the regime changes
            m <- sum(i1 > cum_lamb) + 1
            pos <- i1 - cum_lamb[m - 1] # which i=1,...,d in lambda_{mi}
            main <- substitute(lambda[foo](foo2), list(foo=m, foo2=pos))
          }
        } else { # identification by non-Gaussianity (also cond_dist="ind_Student" or "ind_skewed_t"): B_m params for m=1,...,M
          cum_b <- n_mean_pars + n_ar_pars + (0:(M - 1))*(d^2 - n_zeros/M) # Index after which the regime changes
          if(any(i1 == cum_b + 1)) { # Above n_zeros is the total number for all B_m, so divided by M
            B_row_ind <- rep(1, times=d) # row for each column, resets whenever a new regime starts
          }
          m <- sum(i1 > cum_b) # Which regime
          n_zeros_in_each_column <- vapply(1:d, function(i2) sum(B_constraints[,i2] == 0, na.rm=TRUE), numeric(1))
          zero_positions <- lapply(1:d, function(i2) (1:d)[B_constraints[,i2] == 0 & !is.na(B_constraints[,i2])]) # 0 constr pos in each col
          cum_bc <- c(0, cumsum(d - n_zeros_in_each_column)) # Index in B_m parameters after which a new column in B_m starts
          posb <- i1 - (n_mean_pars + n_ar_pars + (m - 1)*(d^2 - n_zeros/M)) # Index in B_m parameters
          col_ind <- sum(posb > cum_bc) # Above n_zeros is the total number for all B_m, so divided by M
          if(is.null(B_row_ind)) { # Which row of B_row_ind we are in, if started plotting mid B
            B_row_ind <- rep(1, times=d) # row for each column
            if(col_ind > 1) {
              zeros_in_prev_cols <- sum(n_zeros_in_each_column[1:(col_ind-1)])
            } else {
              zeros_in_prev_cols <- 0
            }
            row_ind_zeros_included <- posb - (col_ind*d - zeros_in_prev_cols)
            zero_positions_in_this_row <- zero_positions[[col_ind]]
            while(TRUE) {
              if(B_row_ind[col_ind] %in% zero_positions_in_this_row) {
                B_row_ind[col_ind] <- B_row_ind[col_ind] + 1
              } else {
                break
              }
            }
          } else { # Does not start plotting mid B, B_row_ind already exists
            while(TRUE) {
              if(B_row_ind[col_ind] %in% zero_positions[[col_ind]]) {
                B_row_ind[col_ind] <- B_row_ind[col_ind] + 1
              } else {
                break
              }
            }
          }
          main <- substitute(B[foo2](foo), list(foo=paste0(B_row_ind[col_ind], ", ", col_ind), foo2=m))
          B_row_ind[col_ind] <- B_row_ind[col_ind] + 1
        }
      }
    } else if(i1 <= n_mean_pars + n_ar_pars + n_covmat_pars + n_weight_pars) { # Transition weight parameter: never here for exo mods
      pos <- i1 - n_mean_pars - n_ar_pars - n_covmat_pars # Position in weightpars
      if(is.null(weight_constraints)) {
        if(weight_function == "relative_dens") {
          main <- substitute(alpha[foo], list(foo=pos))
        } else if(weight_function %in% c("logistic", "exponential")) {
          if(pos == 1) {
            main <- substitute(c) # Location parameter
          } else {
            main <- substitute(gamma) # Scale parameter
          }
        } else if(weight_function == "mlogit") {
          # The index after which the regime changes:
          cum_a <- n_mean_pars + n_ar_pars + n_covmat_pars + (0:(M - 2))*(1 + length(weightfun_pars[[1]])*weightfun_pars[[2]])
          m <- sum(i1 > cum_a) # Which regime
          pos0 <- sum(i1 > cum_a[m] + 0:(length(weightfun_pars[[1]])*weightfun_pars[[2]])) # Position in gamma_m
          main <- substitute(gamma[foo](foo2), list(foo=m, foo2=pos0))
        } else if(weight_function == "threshold") {
          main <- substitute(r[foo], list(foo=pos))
        }
      } else {
        if(n_weight_pars > 0) { # Zero weights pars if all are known constants
          main <- substitute(xi[foo], list(foo=pos))
        }
      }
    } else { # Must be a distribution parameter
      if(cond_dist == "Student") {
        main <- substitute(nu) # Must be the only df parameter
      } else if(cond_dist == "ind_Student") {
        which_df <- i1 - n_mean_pars - n_ar_pars - n_covmat_pars - n_weight_pars # Which df parameter (related to which time series)
        main <- substitute(nu[foo], list(foo=which_df))
      } else if(cond_dist == "ind_skewed_t") {
        which_distpar <- i1 - n_mean_pars - n_ar_pars - n_covmat_pars - n_weight_pars # which dist par
        if(which_distpar <= d) {
          main <- substitute(nu[foo], list(foo=which_distpar)) # df param
        } else {
          main <- substitute(lambda[foo], list(foo=which_distpar - d)) # skewness params
        }
      }
    }
    plot(x=vals, y=logliks, type="l", main=main, xlab="", ylab="")
    abline(v=pars[i1], col="red") # Points the estimate
  }
}
