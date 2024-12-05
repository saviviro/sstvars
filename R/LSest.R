#' @title Internal estimation function for estimating autoregressive and threshold parameters of
#'   TVAR models by the method of least squares.
#'
#' @description \code{estim_LS} estimates the autoregressive and threshold parameters of TVAR models
#'   by the method of least squares.
#'
#' @inheritParams loglikelihood
#' @param ncores the number CPU cores to be used in parallel computing.
#' @param min_obs_coef the smallest accepted number of observations (times variables) from each regime
#'  relative to the number of parameters in the regime.
#' @param use_parallel employ parallel computing? If \code{FALSE}, does not print anything.
#' @details Used internally in the multiple phase estimation procedure proposed by Koivisto,
#'  Luoto, and Virolainen (2025). Mean constraints are not supported. Only weight constraints that
#'  specify the threshold parameters as fixed values are supported. Only intercept parametrization is
#'  supported.
#' @return Returns the estimated parameters in a vector of the form
#'  \eqn{(\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M,\alpha}, where
#'  \itemize{
#'     \item{\eqn{\phi_{m,0} = } the \eqn{(d \times 1)} intercept vector of the \eqn{m}th regime.}
#'     \item{\eqn{\varphi_m = (vec(A_{m,1}),...,vec(A_{m,p}))} \eqn{(pd^2 \times 1)}.}
#'     \item{\eqn{\alpha = (r_1,...,r_{M-1})} the \eqn{(M-1\times 1)} vector of the threshold parameters.}
#'  }
#'  For models with...
#'   \describe{
#'     \item{AR_constraints:}{Replace \eqn{\varphi_1,...,\varphi_M} with \eqn{\psi} as described in the argument \code{AR_constraints}.}
#'     \item{weight_constraints:}{If linear constraints are imposed, replace \eqn{\alpha} with \eqn{\xi} as described in the
#'      argument \code{weigh_constraints}. If weight functions parameters are imposed to be fixed values, simply drop \eqn{\alpha}
#'      from the parameter vector.}
#'   }
#' @references
#'  \itemize{
#'    \item Hubrich K., Teräsvirta. T. 2013. Thresholds and Smooth Transitions in Vector Autoregressive Models.
#'      \emph{CREATES Research Paper 2013-18, Aarhus University.}
#'    \item Koivisto T., Luoto J., Virolainen S. 2025. Unpublished working paper.
#'    \item Tsay R. 1998. Testing and Modeling Multivariate Threshold Models.
#'      \emph{Journal of the American Statistical Association}, \strong{93}:443, 1188-1202.
#'    \item Virolainen S. 2024. Identification by non-Gaussianity in structural threshold and
#'       smooth transition vector autoregressive models. Unpublished working
#'       paper, available as arXiv:2404.19707.
#'  }
#' @keywords internal

estim_LS <- function(data, p, M, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                     weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t"),
                     parametrization=c("intercept", "mean"), AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL,
                     penalized=TRUE, penalty_params=c(0.05, 0.2), min_obs_coef=3, use_parallel=TRUE, ncores=2) {
  # Checks
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  if(!is.null(mean_constraints)) {
    stop("Mean constraints are not supported by the LS estimation")
  } else if(weight_function != "threshold") {
    stop("Only threshold weight function is supported by the LS estimation")
  } else if(parametrization != "intercept") {
    stop("Only the intercept parametrization is supported by the LS estimation")
  }
  stopifnot(is.numeric(penalty_params) && length(penalty_params) == 2 && all(penalty_params >= 0) && penalty_params[1] < 1)
  stab_tol <- penalty_params[1]
  tuning_par <- penalty_params[2]
  stopifnot(is.numeric(min_obs_coef) && length(min_obs_coef) == 1 && min_obs_coef > 1)

  # Check the weight constraints
  if(!is.null(weight_constraints)) {
    if(!all(weight_constraints[[1]] == 0)) {
      stop(paste("Only such weight_constraints that specify the threshold parameters some known fixed values",
                 "are supported in the least squares estimation."))
    }
  }
  data <- check_data(data, p=p)
  d <- ncol(data)
  weightfun_pars <- check_weightfun_pars(data=data, p=p, d=d, M=M, weight_function=weight_function,
                                         weightfun_pars=weightfun_pars, cond_dist=cond_dist)
  check_constraints(data=data, p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    parametrization=parametrization, identification="reduced_form",
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=NULL)

  # Obtain relevant statistics
  n_obs <- nrow(data)
  T_obs <- n_obs - p
  pars_per_regime <- d + ifelse(is.null(AR_constraints), p*d^2, ncol(AR_constraints)/M) + d^2
  T_min <- min_obs_coef/d*pars_per_regime # Minimum number of obs in each regime
  if(T_obs/M < T_min) { # Try smaller T_min
    stop(paste("The number of observations is too small for reasonable estimation (according to the argument 'min_obs_coef').",
               "Decrease the order p or the number of regimes M."))
  }

  ################################
  ## Create estim etc functions ##
  ################################

  ## Least squares estimation function given thresholds for models without AR constraints
  LS_without_AR_constraints <- function(thresholds) {
    # threshold = length M-1 vector of the thresholds r_1,...,r_{M-1}; if M=1 anything is ok
    # Other arguments are taken from the parent environment.

    # Storage for the estimates
    all_intercepts <- matrix(NA, nrow=d, ncol=M) # (d x M), [, m] for regime m
    all_AR_mats <- array(NA, dim=c(d, d*p, M)) # [, , m] for A_{m,1} : ... : A_{m,p}

    # Storage for the sums of squares of residuals
    all_rss <- numeric(M) # The sum of squares of residuals for each regime

    # In Y, i:th row denotes the vector \bold{y_{i-1}} = (y_{i-1},...,y_{i-p}) (dpx1),
    # assuming the observed data is y_{-p+1},...,y_0,y_1,...,y_{T}. The last row is for
    # the vector (y_{T},...,y_{T-p}).
    Y <- reform_data(data, p) # (T_obs + 1 x dp)
    Y2 <- Y[1:T_obs, , drop=FALSE] # Last row removed; not needed when only lagged observations used

    # The transition weights: (T x M) matrix, the t:th row is for the time point t and the m:th column is for the regime m.
    alpha_mt <- get_alpha_mt(data, Y2=Y2, p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                             weightpars=thresholds) # Transition weights (T x M), t:th row

    # Go through the regimes
    for(m in 1:M) {
      # Obtain the time periods during which regime m prevails
      m_periods <- which(abs(alpha_mt[, m] - 1) < 1e-6) # The time periods t when regime m prevails

      # Create the matrix Y for regime m; the first p values in data are the initial valus
      Y_m <- data[m_periods + p, , drop=FALSE] # (T_m x d)]

      # Create the matrix X for regime m
      X_m <- cbind(1, Y2[m_periods, , drop=FALSE]) # (T_m x d(p+1))

      # Calculate the lest squares estimates
      C_m <- tryCatch(solve(crossprod(X_m, X_m), crossprod(X_m, Y_m)), # (d x d(p+1)), [\phi_{m,0} : A_{m,1} : ... : A_{m,p}]
                      error=function(e) matrix(0, nrow=d*p + 1, ncol=d)) # zero estimates are legal but bad, dummy estimates
      tC_m <- t(C_m)

      # Store the estimates
      all_intercepts[, m] <- tC_m[, 1]
      all_AR_mats[, , m] <- tC_m[, -1]

      # Calculate the sum of squares of residuals
      U_m <- Y_m - X_m%*%C_m # (T_m x d), residuals
      all_rss[m] <- sum(diag((crossprod(U_m, U_m)))) # The sum of squares of residuals
    }

    # Obtain the estimates in the vector (\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M)
    estims <- c(all_intercepts, all_AR_mats)

    # Return the estimates and the sum of squares of residuals, the last element is the for rss
    c(estims, sum(all_rss))
  }

  ## Least squares estimation function given thresholds for models with AR constraints
  LS_with_AR_constraints <- function(thresholds) {
    # threshold = length M-1 vector of the thresholds r_1,...,r_{M-1}; if M=1 anything is ok
    # Other arguments are taken from the parent environment.

    # Create the constraint matrix \tilde{C} with dummy constraints for the intercepts as well
    q <- ncol(AR_constraints) # The number of intercept and AR parameters under constraints
    C_tilde <- rbind(cbind(diag(rep(1, times=M*d)), matrix(0, nrow=M*d, ncol=q)),
                     cbind(matrix(0, nrow=M*p*d^2, ncol=M*d), AR_constraints))

    # Storage for the estimates
    estims <- array(NA, dim=c(d, d*p + 1, M)) # [, , m] for [\phi_{m,0} : A_{m,1} : ... : A_{m,p}]

    # In Y, i:th row denotes the vector \bold{y_{i-1}} = (y_{i-1},...,y_{i-p}) (dpx1),
    # assuming the observed data is y_{-p+1},...,y_0,y_1,...,y_{T}. The last row is for
    # the vector (y_{T},...,y_{T-p}).
    Y <- reform_data(data, p) # (T_obs + 1 x dp)
    Y2 <- Y[1:T_obs, , drop=FALSE] # Last row removed; not needed when only lagged observations used

    # The transition weights: (T x M) matrix, the t:th row is for the time point t and the m:th column is for the regime m.
    alpha_mt <- get_alpha_mt(data, Y2=Y2, p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                             weightpars=thresholds) # Transition weights (T x M), t:th row

    # Now the observations are stored to a (Td x 1) vector defining the matrix Y.
    Y <- as.matrix(vec(t(data[(p+1):nrow(data),]))) # (Td x 1), we don't take the first p observations (init values)

    # Initialize the matrix of regressors X.
    X <- matrix(0, nrow=d*T_obs, ncol=M*d + M*p*d^2) # (Td x Md + Mpd^2)

    # We construct the (Td x Md + Mpd^2) matrix X consisting of (d x Md + Mpd^2) blocks X_t for each time period t.
    # Each block X_t can be partioned to the intercept part X_t_phi (d x Md) and the AR part X_t_A (d x Mpd^2).
    I_d <- diag(rep(1, times=d)) # The (d x d) identity matrix
    for(t in 1:T_obs) { # We need to loop through all time periods
      m <- which(abs(alpha_mt[t, ] - 1) < 1e-6) # The regime m that prevails at time t
      X_t_phi <- cbind(matrix(0, nrow=d, ncol=(m - 1)*d), I_d, matrix(0, nrow=d, ncol=(M - m)*d)) # (d x Md)
      X_t_A <- matrix(0, nrow=d, ncol=M*p*d^2) # (d x Mpd^2), initialize the AR part
      for(j in 1:p) { # Go through the lags and fill in the blocks
        # The Block_{A_{m,j}} is the columns: (dp*(m - 1)+ d*(j - 1) + 1):(dp*(m - 1) + d*j)
        X_t_A[, ((m - 1)*p*d^2 + (j - 1)*d^2 + 1):((m - 1)*p*d^2 + j*d^2)] <- kronecker(t(data[t + p - j,]), I_d) # (d x d^2) block
      }
      # Will in the X_t block to the X matrix
      X[((t - 1)*d + 1):((t - 1)*d + d), ] <- cbind(X_t_phi, X_t_A) # (d x Md + Mpd^2)
    }

    # Integrate the constraints into the design matrix
    X_bold <- X%*%C_tilde # (Td x q)

    estims <- tryCatch(solve(crossprod(X_bold, X_bold), crossprod(X_bold, Y)), # (M*d + q x 1)
                       error=function(e) matrix(0, nrow=M*d + q, ncol=1)) # zero estimates are legal but bad, dummy estimates

    # Calculate the sum of squares of residuals
    U <- Y - X_bold%*%estims # (Td x 1), residuals
    rss <- crossprod(U, U) # The sum of squares of residuals

    # Return the estimates and the sum of squares of residuals, the last element is the for rss
    c(estims, rss)
  }

  ## A function to check whether the stability condition is satisfied for the AR matrices, and
  ## if not, to what extend it is not satisfied.
  stab_exceeded <- function(estims) {
    # Estims should be a vector of the form (\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M)
    if(!is.null(AR_constraints)) { # Expand the AR constraints
      pars_to_check <- c(estims[1:(M*d)], AR_constraints%*%estims[(M*d + 1):(M*d + ncol(AR_constraints))])
    } else {
      pars_to_check <- estims[1:(M*d + M*p*d^2)]
    }
    all_phi0 <- pick_phi0(M=M, d=d, params=pars_to_check)
    all_A <- pick_allA(p=p, M=M, d=d, params=pars_to_check)
    all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
    all_stab_exceeds <- matrix(nrow=nrow(all_boldA[, , 1]), ncol=M) # Square of how much modulus of eigenvalues exceed 1 - stab_tol
    for(m in 1:M) { # Check stability condition for each regime
      abs_eigs <- abs(eigen(all_boldA[, , m], symmetric=FALSE, only.values=TRUE)$values)
      all_stab_exceeds[, m] <- pmax(0, abs_eigs - (1 - stab_tol))^2 # How much abs eigens exceed 1 - stab_tol, squared
    }
    sum(all_stab_exceeds) # Sum of the squared exceeded values of stab cond
  }

  ###########################################
  ## Create threshold vectors and estimate ##
  ###########################################

  ## Create the set of threshold vectors for the optimization; M=1 will use numeric(0) and run the LS only once
  if(is.null(weight_constraints)) {
    switch_var_series <- data[,weightfun_pars[1]] # The switching variable time series
    sv_sorted_full <- sort(switch_var_series, decreasing=FALSE) # The sorted switch variable series

    # Remove the values from the sorted switch variable series that would leave less than T_min observations
    # smaller or larger than the threshold.
    switch_var_sorted <- sv_sorted_full[T_min:(length(switch_var_series) - T_min)]
    min_switchvar <- min(switch_var_sorted)
    max_switchvar <- max(switch_var_sorted)

    # The maximum number of grid points is calculated so that the number M-1 dimensional of multisets
    # is at most 20000 for M <= 4.
    if(M >= 2) { # M=1 case is separately handled
      max_thresholds <- 1000 # The maximum number of threshold values to considered
      if(M == 2) {
        if(length(switch_var_sorted) < max_thresholds) {
          grid_points <- switch_var_sorted
        } else {
          grid_points <- seq(from=min_switchvar, to=max_switchvar, length.out=max_thresholds)
        }
      } else if(M == 3) {
        grid_points <- seq(from=min_switchvar, to=max_switchvar, length.out=min(200, max_thresholds))
      } else if(M == 4) {
        grid_points <- seq(from=min_switchvar, to=max_switchvar, length.out=50)
      } else {
        grid_points <- seq(from=min_switchvar, to=max_switchvar, length.out=30)
      }
      thresholds <- t(utils::combn(x=grid_points, m=M - 1, simplify=TRUE)) # M-1 dim multisets of lexically ordered grid points
    }
    if(M == 2) {
      thresvecs <- thresholds # Each row for each threshold vector (scalar in this case)
    } else if(M > 2) {
      obs_between_thresholds <-  matrix(NA, nrow=nrow(thresholds), ncol=M - 2) # The number of observations between the thresholds
      for(m in 1:(M - 2)) {
        n_at_most_upper <- findInterval(x=thresholds[, m + 1], vec=sv_sorted_full, left.open=TRUE, rightmost.closed=TRUE)
        n_at_most_lower <- findInterval(x=thresholds[, m], vec=sv_sorted_full, left.open=TRUE, rightmost.closed=TRUE)
        obs_between_thresholds[,m] <- n_at_most_upper - n_at_most_lower # Number of observations between lower and upper threshold
      }
      # The threshold vectors with enough observations in all regimes
      #print(which(rowSums(obs_between_thresholds >= T_min) == ncol(obs_between_thresholds)))
      thresvecs <- thresholds[which(rowSums(obs_between_thresholds >= T_min) == ncol(obs_between_thresholds)), , drop=FALSE]
    }
  } else { # thresholds fixed to known numbers
    thresvecs <- matrix(weight_constraints[[2]], nrow=1, ncol=M - 1)
  }
  # Each row in thresvecs of a threshold vector

  ## Estimate the model for all thresholds vectors in thresvecs
  estim_fun <- if(is.null(AR_constraints)) LS_without_AR_constraints else LS_with_AR_constraints
  estim_length <- if(is.null(AR_constraints)) M*d + M*p*d^2 + 1 else M*d + ncol(AR_constraints) + 1

  if(M == 1) {
    if(use_parallel) message(paste("PHASE 1: Estimating the AR and weight parameters by least squares..."))
    estims <- as.matrix(estim_fun(numeric(0)))
    all_stab_ex <- stab_exceeded(estims[,1])
  } else {
    if(use_parallel) {
      if(ncores > parallel::detectCores()) {
        ncores <- parallel::detectCores()
      }
      n_thresvecs <- ifelse(M == 1 || !is.null(weight_constraints), 1, nrow(thresvecs))
      message(paste0("PHASE 1: Estimating the AR and weight parameters by least squares for ", n_thresvecs,
                     " vectors of thresholds...")) # "PHASE 1" print i related to the multiple-phase estimation procedure
      cl <- parallel::makeCluster(ncores)
      on.exit(try(parallel::stopCluster(cl), silent=TRUE)) # Close the cluster on exit, if not already closed.
      parallel::clusterExport(cl, ls(environment(estim_LS)), envir=environment(estim_LS)) # assign all variables from package:sstvars
      parallel::clusterEvalQ(cl, c(library(pbapply), library(sstvars)))
      estims <- as.matrix(simplify2array(pbapply::pblapply(1:nrow(thresvecs), FUN=function(i1) estim_fun(thresvecs[i1,]), cl=cl)))

      if(penalized) {
        if(M > 2) {
          message(paste0("Checking the stability condition for all the LS estimates..."))
          all_stab_ex <- simplify2array(pbapply::pblapply(1:nrow(thresvecs), FUN=function(i1) stab_exceeded(estims[,i1]), cl=cl))
        } else { # Less prints, since the calculations are fast enough
          all_stab_ex <- simplify2array(pbapply::pblapply(1:nrow(thresvecs), FUN=function(i1) stab_exceeded(estims[,i1]), cl=cl))
        }
      }
      parallel::stopCluster(cl=cl)
    } else { # No parallel computing
      estims <- as.matrix(vapply(1:nrow(thresvecs), FUN=function(i1) estim_fun(thresvecs[i1,]),
                                 FUN.VALUE=numeric(estim_length)))

      if(penalized) {
        all_stab_ex <- vapply(1:nrow(thresvecs), FUN=function(i1) stab_exceeded(estims[,i1]), FUN.VALUE=numeric(1))
      }
    }
  }
  # Each column in estims corresponds to each vector of thresholds

  ## Obtain the LS estimates, possibly among stable estimates
  if(penalized) {
    # Determine the tuning parameter value that controls the extend of the penalization
    all_rss <- estims[nrow(estims),]
    min_rss <- min(all_rss) # The smallest residual sum of squares
    penalty_coef <- tuning_par*min_rss # The penalty coefficient

    # Obtain the index with the smallest penalized sum of squares
    penalized_stab_ex <- penalty_coef*all_stab_ex # Penalization for non-stable estimates
    all_pen_rss <- all_rss + penalized_stab_ex # Penalized sum of squares of residuals
    min_rss_index <- which.min(all_pen_rss)[1] # The index for which the penalized sum of squares is the smallest

  } else {
    # Find the index for which the sum of squares of residuals is the smallest (regardless of stability)
    min_rss_index <- which.min(estims[nrow(estims),])[1]
  }


  ## Obtain and return the estimates corresponding the smallest sum of squares of residuals
  int_and_ar_estims <- estims[1:(nrow(estims) - 1), min_rss_index]
  if(M == 1 || !is.null(weight_constraints)) {
    threshold_estims <- numeric(0) # No threshold estimates
  } else {
    threshold_estims <- thresvecs[min_rss_index,]
  }
  c(int_and_ar_estims, threshold_estims) # Return the estimates
}


#' @title Internal estimation function for estimating autoregressive and weight parameters of
#'   STVAR models by the method of nonlinear least squares.
#'
#' @description \code{estim_NLS} estimates the autoregressive and weight parameters of STVAR models
#'   by the method of least squares (\code{relative_dens} weight function is not supported).
#'
#' @inheritParams estim_LS
#' @details Used internally in the multiple phase estimation procedure proposed by Virolainen (2025).
#'  The weight function \code{relative_dens} is not supported.  Mean constraints are not supported.
#'  Only weight constraints that specify the weight parameters as fixed values are supported.
#'  Only intercept parametrization is supported.
#' @return Returns the estimated parameters in a vector of the form
#'  \eqn{(\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M,\alpha}, where
#'  \itemize{
#'    \item{\eqn{\phi_{m,0} = } the \eqn{(d \times 1)} intercept vector of the \eqn{m}th regime.}
#'    \item{\eqn{\varphi_m = (vec(A_{m,1}),...,vec(A_{m,p}))} \eqn{(pd^2 \times 1)}.}
#'    \item{\eqn{\alpha}} is the vector of the weight parameters: \describe{
#'      \item{\code{weight_function="relative_dens"}:}{\eqn{\alpha = (\alpha_1,...,\alpha_{M-1})}
#'        \eqn{(M - 1 \times 1)}, where \eqn{\alpha_m} \eqn{(1\times 1), m=1,...,M-1} are the transition weight parameters.}
#'      \item{\code{weight_function="logistic"}:}{\eqn{\alpha = (c,\gamma)}
#'        \eqn{(2 \times 1)}, where \eqn{c\in\mathbb{R}} is the location parameter and \eqn{\gamma >0} is the scale parameter.}
#'      \item{\code{weight_function="mlogit"}:}{\eqn{\alpha = (\gamma_1,...,\gamma_M)} \eqn{((M-1)k\times 1)},
#'        where \eqn{\gamma_m} \eqn{(k\times 1)}, \eqn{m=1,...,M-1} contains the multinomial logit-regression coefficients
#'        of the \eqn{m}th regime. Specifically, for switching variables with indices in \eqn{I\subset\lbrace 1,...,d\rbrace}, and with
#'        \eqn{\tilde{p}\in\lbrace 1,...,p\rbrace} lags included, \eqn{\gamma_m} contains the coefficients for the vector
#'        \eqn{z_{t-1} = (1,\tilde{z}_{\min\lbrace I\rbrace},...,\tilde{z}_{\max\lbrace I\rbrace})}, where
#'        \eqn{\tilde{z}_{i} =(y_{it-1},...,y_{it-\tilde{p}})}, \eqn{i\in I}. So \eqn{k=1+|I|\tilde{p}}
#'        where \eqn{|I|} denotes the number of elements in \eqn{I}.}
#'      \item{\code{weight_function="exponential"}:}{\eqn{\alpha = (c,\gamma)}
#'        \eqn{(2 \times 1)}, where \eqn{c\in\mathbb{R}} is the location parameter and \eqn{\gamma >0} is the scale parameter.}
#'      \item{\code{weight_function="threshold"}:}{\eqn{\alpha = (r_1,...,r_{M-1})}
#'         \eqn{(M-1 \times 1)}, where \eqn{r_1,...,r_{M-1}} are the thresholds.}
#'      \item{\code{weight_function="exogenous"}:}{Omit \eqn{\alpha} from the parameter vector.}
#'    }
#'  }
#'  For models with...
#'   \describe{
#'     \item{AR_constraints:}{Replace \eqn{\varphi_1,...,\varphi_M} with \eqn{\psi} as described in the argument \code{AR_constraints}.}
#'     \item{weight_constraints:}{If linear constraints are imposed, replace \eqn{\alpha} with \eqn{\xi} as described in the
#'      argument \code{weigh_constraints}. If weight functions parameters are imposed to be fixed values, simply drop \eqn{\alpha}
#'      from the parameter vector.}
#'   }
#' @references
#'  \itemize{
#'    \item Hubrich K., Teräsvirta. T. 2013. Thresholds and Smooth Transitions in Vector Autoregressive Models.
#'      \emph{CREATES Research Paper 2013-18, Aarhus University.}
#'    \item Virolainen S. 2024. Identification by non-Gaussianity in structural threshold and
#'       smooth transition vector autoregressive models. Unpublished working
#'       paper, available as arXiv:2404.19707.
#'    \item Virolainen S. 2025. Unpublished working paper.
#'  }
#' @keywords internal

estim_NLS <- function(data, p, M, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                      weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t"),
                      parametrization=c("intercept", "mean"), AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL,
                      penalized=TRUE, penalty_params=c(0.05, 0.2), min_obs_coef=3, use_parallel=TRUE, ncores=2) {
  # Checks
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  if(!is.null(mean_constraints)) {
    stop("Mean constraints are not supported by the NLS estimation")
  } else if(weight_function == "relative_dens") {
    stop("The relative_dens weight function is not supported by the NLS estimation")
  } else if(parametrization != "intercept") {
    stop("Only the intercept parametrization is supported by the LS estimation")
  }
  # Check the weight constraints
  if(!is.null(weight_constraints)) {
    if(!all(weight_constraints[[1]] == 0)) {
      stop(paste("Only such weight_constraints that specify the threshold parameters some known fixed values",
                 "are supported in the least squares estimation."))
    }
  }
  data <- check_data(data, p=p)
  d <- ncol(data)
  weightfun_pars <- check_weightfun_pars(data=data, p=p, d=d, M=M, weight_function=weight_function,
                                         weightfun_pars=weightfun_pars, cond_dist=cond_dist)
  check_constraints(data=data, p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    parametrization=parametrization, identification="reduced_form",
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=NULL)

  stopifnot(is.numeric(penalty_params) && length(penalty_params) == 2 && all(penalty_params >= 0) && penalty_params[1] < 1)
  stab_tol <- penalty_params[1]
  tuning_par <- penalty_params[2]
  stopifnot(is.numeric(min_obs_coef) && length(min_obs_coef) == 1 && min_obs_coef > 1)

  # Obtain relevant statistics
  n_obs <- nrow(data)
  T_obs <- n_obs - p
  if(weight_function != "exogenous") {
    T_min <- min_obs_coef/d*ifelse(is.null(AR_constraints), (p + 1)*d^2 + d,
                        ncol(AR_constraints)/M + d + d^2) # Minimum number of obs in each regime
    if(T_obs/M < T_min) { # Try smaller T_min
       stop(paste("The number of observations is too small for reasonable estimation (according to the argument 'min_obs_coef').",
                  "Decrease the order p or the number of regimes M."))
    }
  }

  ################################
  ## Create estim etc functions ##
  ################################

  # Logarithm of the smallest value that can be handled normally (used in get_alpha_mt)
  epsilon <- round(log(.Machine$double.xmin) + 10)

  ## Nonlinear least squares estimation function given weight parameters
  NLS_est <- function(weightpars, AR_constraints) {
    # weightpars = vector of the weight parameters, AR_constraints = as usual
    # Other arguments are taken from the parent env

    if(!is.null(AR_constraints)) {
      # Create the matrix C_tilde that sets dummy constraints for the intercepts:
      q <- ncol(AR_constraints) # The number of intercept and AR parameters under constraints
      C_tilde <- rbind(cbind(diag(rep(1, times=M*d)), matrix(0, nrow=M*d, ncol=q)),
                       cbind(matrix(0, nrow=M*p*d^2, ncol=M*d), AR_constraints))

      # Create the permutation matrix that expand \tilde{C}\tilde{\psi} to \beta for AR constraints:
      P <- matrix(0, nrow=M*d + M*p*d^2, ncol=M*d + M*p*d^2)
      I_d <- diag(nrow=d)
      I_pd2 <- diag(nrow=p*d^2)
      for(m in 1:M) { # Go through the regimes
        # Insert the identity matrix for the intercepts
        P[((m - 1)*d + (m - 1)*p*d^2 + 1):(m*d + (m - 1)*p*d^2), ((m - 1)*d + 1):(m*d)] <- I_d

        # Insert the identity matrix for the AR matrices
        P[(m*d + (m - 1)*p*d^2 + 1):(m*d + m*p*d^2), (M*d + (m - 1)*p*d^2 + 1):(M*d + m*p*d^2)] <- I_pd2
      }

      # Calculate the multiplication of P and C_tilde that expands psi_tilde to the usual parameter vector:
      PC <- P%*%C_tilde
    }

    # In Y, i:th row denotes the vector \bold{y_{i-1}} = (y_{i-1},...,y_{i-p}) (dpx1),
    # assuming the observed data is y_{-p+1},...,y_0,y_1,...,y_{T}. The last row is for
    # the vector (y_{T},...,y_{T-p}).
    Y <- reform_data(data, p) # (T_obs + 1 x dp)
    Y2 <- Y[1:T_obs, , drop=FALSE] # Last row removed; not needed when only lagged observations used

    # The transition weights: (T x M) matrix, the t:th row is for the time point t and the m:th column is for the regime m.
    alpha_mt <- get_alpha_mt(data, Y2=Y2, p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                             weightpars=weightpars, epsilon=epsilon) # Transition weights (T x M), t:th row

    # Check whether there are enough observations from each regime
    if(weight_function != "exogenous") {
      if(any(colSums(alpha_mt) < T_min)) { # Enough observations from each regime?
        all_Psi <- array(0, dim=c(d, M*d + M*p*d^2, T_obs)) # [, , t] for \Psi_t, dummy Psi
        if(is.null(AR_constraints)) {
          estims <- matrix(0, nrow=M*d + M*p*d^2, ncol=1) # Dummy estimates, legal but very bad
        } else {
          estims <- matrix(0, nrow=ncol(C_tilde), ncol=1) # Dummy estimates, legal but very bad
        }
        calc_estims <- FALSE
      } else {
        calc_estims <- TRUE
      }
    } else {
      calc_estims <- TRUE
    }

    if(calc_estims) { # Exo weights or enough observations, no dummy estims at least yet
      # Storage for the data matrices \Psi_t in the form y_t = \Psi_t\beta + u_t,
      # \beta = (vec(\Phi_1),...,vec(\Phi_M)), \Phi_m = [\phi_{m,0} : A_{m,1} : ... : A_{m,p}].
      all_Psi <- array(NA, dim=c(d, M*d + M*p*d^2, T_obs)) # [, , t] for \Psi_t
      all_cPsi <- array(NA, dim=c(M*d + M*p*d^2, M*d + M*p*d^2, T_obs)) # [, , t] for crossprod(\Psi_t)
      all_tPsiy <- array(NA, dim=c(M*d + M*p*d^2, 1, T_obs)) # [, , t] for \Psi_t'y_t

      # Storages for models with AR constraints
      if(!is.null(AR_constraints)) {
        all_PsiPC <- array(NA, dim=c(d, ncol(C_tilde), T_obs)) # [, , t] for \Psi_tPC
        all_cPsiPC <- array(NA, dim=c(ncol(C_tilde), ncol(C_tilde), T_obs)) # [, , t] for crossprod(\Psi_tPC)
        all_tPsiPCy <- array(NA, dim=c(ncol(C_tilde), 1, T_obs)) # [,  ,t] for (\Psi_tPC)'y_t
      }

      # Go through the time periods
      for(t in 1:T_obs) {

        # Create the vector Z_t = (1, y_{t-1},...,y_{t-p}) (1 + dp x 1), and calculate Knonecker product between Z_t' and I_d
        Z_kron <- kronecker(matrix(c(1, Y[t,]), nrow=1), diag(d)) # (d x d(1 + dp)), the Kronecker product

        # Go through the regimes
        for(m in 1:M) {
          all_Psi[, ((m - 1)*(d + p*d^2) + 1):(m*(d + p*d^2)), t] <- alpha_mt[t, m]*Z_kron # Fill in Psi_t, also for AR_constraints
        }

        if(!is.null(AR_constraints)) {
          all_PsiPC[, , t] <- all_Psi[, , t]%*%PC # Multiply Psi_t with PC
          all_cPsiPC[, , t] <- crossprod(all_PsiPC[, , t]) # Fill in Psi_t'Psi_t
          all_tPsiPCy[, , t] <- crossprod(all_PsiPC[, , t], data[p + t,]) # Fill in \Psi_t'y_t
        } else {
          all_cPsi[, , t] <- crossprod(all_Psi[, , t]) # Fill in Psi_t'Psi_t
          all_tPsiy[, , t] <- crossprod(all_Psi[, , t], data[p + t,]) # Fill in \Psi_t'y_t
        }
      }

      # Sum over the time periods in cPsi and tPsiy
      if(!is.null(AR_constraints)) {
        sum_cPsi <- apply(all_cPsiPC, MARGIN=c(1, 2), FUN=sum)
        sum_tPsiy <- apply(all_tPsiPCy, MARGIN=c(1, 2), FUN=sum)
      } else {
        sum_cPsi <- apply(all_cPsi, MARGIN=c(1, 2), FUN=sum)
        sum_tPsiy <- apply(all_tPsiy, MARGIN=c(1, 2), FUN=sum)
      }
      ## Obtain the estimates in the vector \beta = (vec(\Phi_1),...,vec(\Phi_M)), where \Phi_m = [\phi_{m,0} : A_{m,1} : ... : A_{m,p}]
      # (with AR_constraints \beta = (\phi_{1,0},...,\phi_{M,0},\psi)
      #estims <- solve(sum_cPsi, sum_tPsiy) # (M*d + M*p*d^2 x 1)
      estims <- tryCatch(solve(sum_cPsi, sum_tPsiy), # (M*d + M*p*d^2 x 1), fails if the system is singular
                         error=function(e) {
                           if(is.null(AR_constraints)) {
                             return(matrix(0, nrow=M*d + M*p*d^2, ncol=1))
                           } else {
                             return(matrix(0, nrow=ncol(C_tilde), ncol=1))
                           }
                         }) # zero estimates are legal but bad, dummy estimates
    }

    ## Calculate the residual sums of squares
    if(is.null(AR_constraints)) {
      est_to_use <- estims
    } else {
      est_to_use <- PC%*%estims # Estimates in the standard form
    }
    all_rss <- numeric(T_obs) # Storage for all residual sums of squares
    for(t in 1:T_obs) { # Go through the time periods again
      u_t <- data[p + t,] - all_Psi[, , t]%*%est_to_use # The residual for time period t
      all_rss[t] <- crossprod(u_t, u_t) # The residual sum of squares for time period t
    }

    ## Obtain the estimates in the vector (\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M)
    # (for models with AR constraints (\phi_{1,0},...,\phi_{M,0},\psi), already in the correct form):
    if(is.null(AR_constraints)) {
      estims <- matrix(estims, nrow=d) # estims in the form [Phi_1,...,Phi_M]
      int_cols <- 0:(M - 1)*(d*p + 1) + 1 # which columns have the intercept parameters
      estims <- c(estims[, int_cols], estims[, -int_cols]) # Estimates in the form (\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M)
    }

    ## Return the estimates and the sum of squares of residuals, the last element is the for rss
    c(estims, sum(all_rss))
  }

  ## A function to check whether the stability condition is satisfied for the AR matrices, and
  ## if not, to what extend it is not satisfied.
  stab_exceeded <- function(estims) {
    # Estims should be a vector of the form (\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M)
    if(!is.null(AR_constraints)) { # Expand the AR constraints
      pars_to_check <- c(estims[1:(M*d)], AR_constraints%*%estims[(M*d + 1):(M*d + ncol(AR_constraints))])
    } else {
      pars_to_check <- estims[1:(M*d + M*p*d^2)]
    }
    all_phi0 <- pick_phi0(M=M, d=d, params=pars_to_check)
    all_A <- pick_allA(p=p, M=M, d=d, params=pars_to_check)
    all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
    all_stab_exceeds <- matrix(nrow=nrow(all_boldA[, , 1]), ncol=M) # Square of how much modulus of eigenvalues exceed 1 - stab_tol
    for(m in 1:M) { # Check stability condition for each regime
      abs_eigs <- abs(eigen(all_boldA[, , m], symmetric=FALSE, only.values=TRUE)$values)
      all_stab_exceeds[, m] <- pmax(0, abs_eigs - (1 - stab_tol))^2 # How much abs eigens exceed 1 - stab_tol, squared
    }
    sum(all_stab_exceeds) # Sum of the squared exceeded values of stab cond
  }

  ############################################
  ## Create weight par vectors and estimate ##
  ############################################

  ## Create the set of weight parameters for the optimization; M=1 will use numeric(0) and run the NLS only once
  if(is.null(weight_constraints) && weight_function != "exogenous") {
    if(weight_function != "mlogit") {
      switch_var_series <- data[,weightfun_pars[1]] # The switching variable time series
      sv_sorted_full <- sort(switch_var_series, decreasing=FALSE) # The sorted switch variable series
      switch_var_sorted <- sv_sorted_full[T_min:(length(switch_var_series) - T_min)] # Switch var vals >=T_min vals below and above
      min_switchvar <- min(switch_var_sorted)
      max_switchvar <- max(switch_var_sorted)
    }
    if(weight_function %in% c("logistic", "exponential")) {
      # Here always M=2, the first weight parameter is location parameter, and the second one is strictly positive scale parameter
      c_grid <- seq(from=min_switchvar, to=max_switchvar, length.out=100) # The grid for the location parameter
      gamma_grid <- seq(from=0.1, to=100, length.out=100) # The grid for the scale parameter
      weightparvecs <- unname(simplify2array(expand.grid(c_grid, gamma_grid))) # Each row for a vector of weight parameters
    } else if(weight_function == "mlogit") {
      n_weight_pars <- (M - 1)*(1 + length(weightfun_pars[[1]])*weightfun_pars[[2]])
      weightparvecs <- unname(simplify2array(expand.grid(replicate(n=n_weight_pars,
                                                                   expr=seq(from=-20, to=20,
                                                                            length.out=max(2, ceiling(10000^(1/n_weight_pars)))),
                                                                   simplify=FALSE)))) # Roughly 10000-100000 grid points
    } else if(weight_function == "threshold") {
      # Remove the values from the sorted switch variable series that would leave less than T_min observations
      # smaller or larger than the threshold.
      switch_var_sorted <- sv_sorted_full[T_min:(length(switch_var_series) - T_min)]
      min_switchvar <- min(switch_var_sorted)
      max_switchvar <- max(switch_var_sorted)

      # The maximum number of grid points is calculated so that the number M-1 dimensional of multisets
      # is at most 20000 for M <= 4.
      if(M >= 2) { # M=1 case is separately handled
        max_thresholds <- 1000 # The maximum number of threshold values to considered
        if(M == 2) {
          if(length(switch_var_sorted) < max_thresholds) {
            grid_points <- switch_var_sorted
          } else {
            grid_points <- seq(from=min_switchvar, to=max_switchvar, length.out=max_thresholds)
          }
        } else if(M == 3) {
          grid_points <- seq(from=min_switchvar, to=max_switchvar, length.out=min(200, max_thresholds))
        } else if(M == 4) {
          grid_points <- seq(from=min_switchvar, to=max_switchvar, length.out=50)
        } else {
          grid_points <- seq(from=min_switchvar, to=max_switchvar, length.out=30)
        }
        thresholds <- t(utils::combn(x=grid_points, m=M - 1, simplify=TRUE)) # M-1 dim multisets of lexically ordered grid points
      }
      if(M == 2) {
        weightparvecs <- thresholds # Each row for each threshold vector (scalar in this case)
      } else if(M > 2) {
        obs_between_thresholds <-  matrix(NA, nrow=nrow(thresholds), ncol=M - 2) # The number of observations between the thresholds
        for(m in 1:(M - 2)) {
          n_at_most_upper <- findInterval(x=thresholds[, m + 1], vec=sv_sorted_full, left.open=TRUE, rightmost.closed=TRUE)
          n_at_most_lower <- findInterval(x=thresholds[, m], vec=sv_sorted_full, left.open=TRUE, rightmost.closed=TRUE)
          obs_between_thresholds[,m] <- n_at_most_upper - n_at_most_lower # Number of observations between lower and upper threshold
        }
        # The threshold vectors with enough observations in all regimes
        weightparvecs <- thresholds[which(rowSums(obs_between_thresholds >= T_min) == ncol(obs_between_thresholds)), , drop=FALSE]
      }
    }

  } else if(weight_function == "exogenous") { # Exogenous weights
    weightparvecs <- weightfun_pars
  } else { # weight parameters fixed to known numbers
    weightparvecs <- matrix(weight_constraints[[2]], nrow=1)
  }
  # Each row in weightparvecs correspond to one vector of weight parameters

  ## Estimate the model for all weight pars in weight_pars
  estim_length <- if(is.null(AR_constraints)) M*d + M*p*d^2 + 1 else M*d + ncol(AR_constraints) + 1

  if(M == 1) {
    if(use_parallel) message(paste("PHASE 1: Estimating the AR and weight parameters by nonlinear least squares..."))
    estims <- as.matrix(NLS_est(numeric(0), AR_constraints=AR_constraints))
    all_stab_ex <- stab_exceeded(estims[,1])
  } else {
    if(use_parallel) {
      if(ncores > parallel::detectCores()) {
        ncores <- parallel::detectCores()
      }
      n_weightvecs <- ifelse(M == 1 || !is.null(weight_constraints), 1, nrow(weightparvecs))
      message(paste0("PHASE 1: Estimating the AR and weight parameters by least squares for ", n_weightvecs,
                     " vectors of weight parameters...")) # "PHASE 1" print i related to the multiple-phase estimation procedure
      cl <- parallel::makeCluster(ncores)
      on.exit(try(parallel::stopCluster(cl), silent=TRUE)) # Close the cluster on exit, if not already closed.
      parallel::clusterExport(cl, ls(environment(estim_LS)), envir=environment(estim_LS)) # assign all variables from package:sstvars
      parallel::clusterEvalQ(cl, c(library(pbapply), library(sstvars)))
      estims <- as.matrix(simplify2array(pbapply::pblapply(1:nrow(weightparvecs),
                                                           FUN=function(i1) NLS_est(weightparvecs[i1,],
                                                                                    AR_constraints=AR_constraints), cl=cl)))

      if(penalized) {
        if(M > 2) {
          message(paste0("Checking the stability condition for all the LS estimates..."))
          all_stab_ex <- simplify2array(pbapply::pblapply(1:nrow(weightparvecs), FUN=function(i1) stab_exceeded(estims[,i1]), cl=cl))
        } else { # Less prints, since the calculations are fast enough
          all_stab_ex <- simplify2array(pbapply::pblapply(1:nrow(weightparvecs), FUN=function(i1) stab_exceeded(estims[,i1]), cl=cl))
        }
      }
      parallel::stopCluster(cl=cl)
    } else { # No parallel computing
      estims <- as.matrix(vapply(1:nrow(weightparvecs), FUN=function(i1) NLS_est(weightparvecs[i1,], AR_constraints=AR_constraints),
                                 FUN.VALUE=numeric(estim_length)))

      if(penalized) {
        all_stab_ex <- vapply(1:nrow(weightparvecs), FUN=function(i1) stab_exceeded(estims[,i1]), FUN.VALUE=numeric(1))
      }
    }
  }
  # Each column in estims corresponds to each vector of thresholds

  ## Obtain the LS estimates, possibly among stable estimates
  if(penalized) {
    # Determine the tuning parameter value that controls the extend of the penalization
    all_rss <- estims[nrow(estims),]
    min_rss <- min(all_rss) # The smallest residual sum of squares
    penalty_coef <- tuning_par*min_rss # The penalty coefficient

    # Obtain the index with the smallest penalized sum of squares
    penalized_stab_ex <- penalty_coef*all_stab_ex # Penalization for non-stable estimates
    all_pen_rss <- all_rss + penalized_stab_ex # Penalized sum of squares of residuals
    min_rss_index <- which.min(all_pen_rss)[1] # The index for which the penalized sum of squares is the smallest

  } else {
    # Find the index for which the sum of squares of residuals is the smallest (regardless of stability)
    min_rss_index <- which.min(estims[nrow(estims),])[1]
  }


  ## Obtain and return the estimates corresponding the smallest sum of squares of residuals
  int_and_ar_estims <- estims[1:(nrow(estims) - 1), min_rss_index]
  if(M == 1 || !is.null(weight_constraints) || weight_function == "exogenous") {
    weightpar_estims <- numeric(0) # No weightpar estimates
  } else {
    weightpar_estims <- weightparvecs[min_rss_index,]
  }
  c(int_and_ar_estims, weightpar_estims) # Return the estimates
}
