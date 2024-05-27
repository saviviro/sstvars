#' @title Estimate linear impulse response function based on a single regime of a structural STVAR model.
#'
#' @description \code{linear_IRF} estimates linear impulse response function based on a single regime
#'   of a structural STVAR model.
#'
#' @inheritParams fitbsSSTVAR
#' @param stvar an object of class \code{'stvar'} defining a structural or reduced form
#'   STVAR model. For a reduced form model, the shocks are automatically identified by
#'   the lower triangular Cholesky decomposition.
#' @param N a positive integer specifying the horizon how far ahead should the
#'   linear impulse responses be calculated.
#' @param regime Based on which regime the linear IRF should be calculated?
#'   An integer in \eqn{1,...,M}.
#' @param which_cumulative a numeric vector with values in \eqn{1,...,d}
#'   (\code{d=ncol(data)}) specifying which the variables for which the linear impulse
#'   responses should be cumulative. Default is none.
#' @param scale should the linear IRFs to some of the shocks be scaled so that they
#'   correspond to a specific instantaneous response of some specific
#'   variable? Provide a length three vector where the shock of interest
#'   is given in the first element (an integer in \eqn{1,...,d}), the variable of
#'   interest is given in the second element (an integer in \eqn{1,...,d}), and
#'   its instantaneous response in the third element (a non-zero real number).
#'   If the linear IRFs of multiple shocks should be scaled, provide a matrix which has one
#'   column for each of the shocks with the columns being the length three vectors described above.
#' @param ci a real number in \eqn{(0, 1)} specifying the confidence level of the
#'   confidence intervals calculated via a fixed-design wild residual bootstrap method.
#'   Available only for models that impose linear autoregressive dynamics
#'   (excluding changes in the volatility regime).
#' @param bootstrap_reps the number of bootstrap repetitions for estimating confidence bounds.
#' @param ncores the number of CPU cores to be used in parallel computing when bootstrapping confidence bounds.
#' @param seed a real number initializing the seed for the random generator.
#' @param ... parameters passed to the plot method \code{plot.irf} that plots
#'   the results.
#' @details If the autoregressive dynamics of the model are linear (i.e., either M == 1 or mean and AR parameters
#'   are constrained identical across the regimes), confidence bounds can be calculated based on a fixed-design
#'   wild residual bootstrap method. We employ the method described in Herwartz and Lütkepohl (2014); see also
#'   the relevant chapters in Kilian and Lütkepohl (2017).
#'
#'   Employs the estimation function \code{optim} from the package \code{stats} that implements the optimization
#'   algorithms. The robust optimization method Nelder-Mead is much faster than SANN but can get stuck at a local
#'   solution. See \code{?optim} and the references therein for further details.
#'
#'   For model identified by non-Gaussianity, the signs and ordering of the shocks are normalized by assuming
#'   that the first non-zero element of each column of the impact matrix of Regime 1 is strictly positive and they are
#'   in a decreasing order. Use the argument \code{scale} to obtain IRFs scaled for specific impact responses.
#' @return Returns a class \code{'irf'} list with  with the following elements:
#'   \describe{
#'     \item{\code{$point_est}:}{a 3D array \code{[variables, shock, horizon]} containing the point estimates of the IRFs.
#'        Note that the first slice is for the impact responses and the slice i+1 for the period i. The response of the
#'        variable 'i1' to the shock 'i2' is subsetted as \code{$point_est[i1, i2, ]}.}
#'     \item{\code{$conf_ints}:}{bootstrapped confidence intervals for the IRFs in a \code{[variables, shock, horizon, bound]}
#'        4D array. The lower bound is obtained as \code{$conf_ints[, , , 1]}, and similarly the upper bound as
#'         \code{$conf_ints[, , , 2]}. The subsetted 3D array is then the bound in a form similar to \code{$point_est}.}
#'     \item{\code{$all_bootstrap_reps}:}{IRFs from all of the bootstrap replications in a \code{[variables, shock, horizon, rep]}.
#'        4D array. The IRF from replication i1 is obtained as \code{$all_bootstrap_reps[, , , i1]}, and the subsetted 3D array
#'        is then the in a form similar to \code{$point_est}.}
#'     \item{Other elements:}{contains some of the arguments the \code{linear_IRF} was called with.}
#'   }
#' @seealso \code{\link{GIRF}}, \code{\link{GFEVD}}, \code{\link{fitSTVAR}}, \code{\link{STVAR}},
#'   \code{\link{reorder_B_columns}}, \code{\link{swap_B_signs}}
#' @references
#'  \itemize{
#'    \item Herwartz H. and Lütkepohl H. 2014. Structural vector autoregressions with Markov switching:
#'      Combining conventional with statistical identification of shocks. \emph{Journal of Econometrics},
#'      183, pp. 104-116.
#'    \item Kilian L. and Lütkepohl H. 2017. Structural Vectors Autoregressive Analysis.
#'          \emph{Cambridge University Press}, Cambridge.
#'  }
#' @examples
#' \donttest{
#' ## These are long running examples that take approximately 10 seconds to run.
#' ## A small number of bootstrap replications is used below to shorten the
#' ## running time (in practice, a larger number of replications should be used).
#'
#' # p=1, M=1, d=2, linear VAR model with independent Student's t shocks identified
#' # by non-Gaussianity (arbitrary weight function applied here):
#' theta_112it <- c(0.644, 0.065, 0.291, 0.021, -0.124, 0.884, 0.717, 0.105, 0.322,
#'   -0.25, 4.413, 3.912)
#' mod112 <- STVAR(data=gdpdef, p=1, M=1, params=theta_112it, cond_dist="ind_Student",
#'  identification="non-Gaussianity", weight_function="threshold", weightfun_pars=c(1, 1))
#' mod112 <- swap_B_signs(mod112, which_to_swap=1:2)
#'
#' # Estimate IRFs 20 periods ahead, bootstrapped 90% confidence bounds based on
#' # 10 bootstrap replications. Linear model so robust estimation methods are
#' # not required.
#' irf1 <- linear_IRF(stvar=mod112, N=20, regime=1, ci=0.90, bootstrap_reps=1,
#'  robust_method="none", seed=1, ncores=1)
#' plot(irf1)
#' print(irf1, digits=3)
#'
#' # p=1, M=2, d=2, Gaussian STVAR with relative dens weight function,
#' # shocks identified recursively.
#' theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863,
#'   -0.3124, 0.505251, 0.298483, 0.030096, -0.176925, 0.838898, 0.310863, 0.007512,
#'   0.018244, 0.949533, -0.016941, 0.121403, 0.573269)
#' mod122 <- STVAR(data=gdpdef, p=1, M=2, params=theta_122relg, identification="recursive")
#'
#' # Estimate IRF based on the first regime 30 period ahead. Scale IRFs so that
#' # the instantaneous response of the first variable to the first shock is 0.3,
#' # and the response of the second variable to the second shock is 0.5.
#' # response of the Confidence bounds
#' # are not available since the autoregressive dynamics are nonlinear.
#' irf2 <- linear_IRF(stvar=mod122, N=30, regime=1, scale=cbind(c(1, 1, 0.3), c(2, 2, 0.5)))
#' plot(irf2)
#'
#'  # Estimate IRF based on the second regime without scaling the IRFs:
#' irf3 <- linear_IRF(stvar=mod122, N=30, regime=2)
#' plot(irf3)
#'
#' # p=3, M=2, d=3, Students't logistic STVAR model with the first lag of the second
#' # variable as the switching variable. Autoregressive dynamics restricted linear,
#' # but the volatility regime varies in time, allowing the shocks to be identified
#' # by conditional heteroskedasticity.
#' theta_322 <- c(0.7575, 0.6675, 0.2634, 0.031, -0.007, 0.5468, 0.2508, 0.0217, -0.0356,
#'  0.171, -0.083, 0.0111, -0.1089, 0.1987, 0.2181, -0.1685, 0.5486, 0.0774, 5.9398, 3.6945,
#'  1.2216, 8.0716, 8.9718)
#' mod322 <- STVAR(data=gdpdef, p=3, M=2, params=theta_322, weight_function="logistic",
#'   weightfun_pars=c(2, 1), cond_dist="Student", mean_constraints=list(1:2),
#'   AR_constraints=rbind(diag(3*2^2), diag(3*2^2)), identification="heteroskedasticity",
#'   parametrization="mean")
#'
#' ## Estimate IRFs 30 periods ahead, bootstrapped 90% confidence bounds based on
#' # 10 bootstrap replications. Responses of the second variable are accumulated.
#' irf4 <- linear_IRF(stvar=mod322, N=30, regime=1, ci=0.90, bootstrap_reps=10,
#'  which_cumulative=2, seed=1)
#' plot(irf4)
#' }
#' @export

linear_IRF <- function(stvar, N=30, regime=1, which_cumulative=numeric(0), scale=NULL, ci=NULL,
                       bootstrap_reps=100, ncores=2, robust_method=c("Nelder-Mead", "SANN", "none"),
                       maxit_robust=1000, seed=NULL, ...) {
  # Get the parameter values etc
  stopifnot(all_pos_ints(c(N, regime, ncores, bootstrap_reps)))
  stopifnot(regime <= stvar$model$M)
  stopifnot(!is.null(stvar$data))
  if(!is.null(seed)) stopifnot(is.numeric(seed) && length(seed) == 1)
  if(stvar$model$identification == "reduced_form" && stvar$model$cond_dist == "ind_Student") {
    stvar$model$identification <- "non-Gaussianity" # Readily identified by non-Gaussianity
  }
  data <- stvar$data
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  params <- stvar$params
  weight_function <- stvar$model$weight_function
  weightfun_pars <- check_weightfun_pars(data=data, p=p, M=M, d=d, weight_function=weight_function,
                                         weightfun_pars=stvar$model$weightfun_pars)
  cond_dist <- stvar$model$cond_dist
  parametrization <- stvar$model$parametrization
  identification <- stvar$model$identification
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  weight_constraints <- stvar$model$weight_constraints
  B_constraints <- stvar$model$B_constraints
  stopifnot(regime <= M)
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params,
                                    weight_function=weight_function, weightfun_pars=weightfun_pars,
                                    cond_dist=cond_dist, identification=identification,
                                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                    weight_constraints=weight_constraints, B_constraints=B_constraints)
  if(stvar$model$parametrization == "mean") {
    params <- change_parametrization(p=p, M=M, d=d, params=params,
                                     weight_function=weight_function, weightfun_pars=weightfun_pars,
                                     identification=identification, cond_dist=cond_dist,
                                     AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL,
                                     B_constraints=NULL, change_to="intercept")
  }
  all_mu <- get_regime_means(p=p, M=M, d=d, params=params,
                             weight_function=weight_function, weightfun_pars=weightfun_pars,
                             cond_dist=cond_dist, parametrization="intercept",
                             identification=identification,
                             AR_constraints=NULL, mean_constraints=NULL,
                             weight_constraints=NULL, B_constraints=NULL)
  all_phi0 <- pick_phi0(M=M, d=d, params=params)
  all_A <- pick_allA(p=p, M=M, d=d, params=params)
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, cond_dist=cond_dist,
                           identification=identification)
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  weightpars <- pick_weightpars(p=p, M=M, d=d, params=params, weight_function=weight_function,
                                weightfun_pars=weightfun_pars, cond_dist=cond_dist)
  distpars <- pick_distpars(d=d, params=params, cond_dist=cond_dist)

  # Check whether it is possible to calculate confidence intervals by wild residual bootstrap
  AR_mats_identical <- all(apply(all_boldA, MARGIN=3, FUN=function(x) identical(x, all_boldA[,,1])))
  means_identical <- !is.null(mean_constraints) && length(mean_constraints) == 1 && all(mean_constraints[[1]] == 1:M)
  ci_possible <- (means_identical && AR_mats_identical) || M == 1

  # Check the argument scale and which_cumulative
  if(identification == "heteroskedasticity" || identification == "non-Gaussianity") { # ind_Students mods ident set to non-Gaus
    if(is.null(B_constraints)) {
      B_constrs <- matrix(NA, nrow=d, ncol=d)
    } else {
      B_constrs <- B_constraints
    }
  } else { # Reduced form or recursive identification
    B_constrs <- matrix(NA, nrow=d, ncol=d)
    B_constrs[upper.tri(B_constrs)] <- 0
    diag(B_constrs) <- 1 # Lower triangular Cholesky constraints
  }
  if(!is.null(scale)) {
    scale <- as.matrix(scale)
    stopifnot(all(scale[1,] %in% 1:d)) # All shocks in 1,...,d
    stopifnot(length(unique(scale[1,])) == length(scale[1,])) # No duplicate scales for the same shock
    stopifnot(all(scale[2,] %in% 1:d)) # All variables in 1,...,d
    stopifnot(all(scale[3,] != 0)) # No zero initial magnitudes

    # For the considered shocks, check that there are no zero-constraints for the variable
    # whose initial response is scaled.
    for(i1 in 1:ncol(scale)) {
      if(!is.na(B_constrs[scale[2, i1], scale[1, i1]]) && B_constrs[scale[2, i1], scale[1, i1]] == 0) {
        if(identification == "heteroskedasticity") {
          stop(paste("Instantaneous response of the variable that has a zero constraint for",
                     "the considered shock cannot be scaled"))
        } else { # Reduced form or recursive identification
          stop(paste("Instantaneous response of the variable that has a zero constraint for the considered",
                     "shock cannot be scaled",
                     "(lower triangular recursive identification is assumed for reduced form models)"))
        }
      }
    }
  }
  if(length(which_cumulative) > 0) {
    which_cumulative <- unique(which_cumulative)
    stopifnot(all(which_cumulative %in% 1:d))
  }

  # Obtain the impact matrix of the regime the IRF is to calculated for
  if(identification == "heteroskedasticity") { # Shocks identified by heteroskedasticity
    W <- pick_W(p=p, M=M, d=d, params=params, identification=identification)
    if(regime > 1) { # Include lambdas
      lambdas <- matrix(pick_lambdas(p=p, M=M, d=d, params=params, identification=identification), nrow=d, byrow=FALSE)
      Lambda_m <- diag(lambdas[, regime - 1])
      B_matrix <- W%*%sqrt(Lambda_m)
    } else { # regime == 1
      B_matrix <- W
    }
  } else if(identification == "non-Gaussianity") { # ind_Student mods ident set to non-Gaus
    if(is.null(ci) || !ci_possible) {
      B_matrix <- all_Omega[, , regime] # impact matrix readily parametrized
    } else {
      # ci calculated, so impact matrices normalized so that the first nonzero element
      # in each column of B_1 is strictly positive and they are in a decreasing order.
      # Determine which columns should swap signs:
      first_non_zero_entries <- vapply(1:d, function(i1) all_Omega[, i1, 1][all_Omega[, i1, 1] != 0][1], numeric(1))
      cols_to_swap <- which(first_non_zero_entries < 0) # Swap signs of the columns with negative first non-zero element

      # Swap the signs of the corresponding columns of the impact matrices
      new_stvar <- swap_B_signs(stvar, which_to_swap=cols_to_swap, calc_std_errors=FALSE)
      tmp_params <- reform_constrained_pars(p=p, M=M, d=d, params=new_stvar$params,
                                            weight_function=weight_function, weightfun_pars=weightfun_pars,
                                            cond_dist=cond_dist, identification=identification,
                                            AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                            weight_constraints=weight_constraints, B_constraints=B_constraints)
      B_matrix <- pick_Omegas(p=p, M=M, d=d, params=tmp_params, cond_dist=cond_dist,
                              identification=identification)[, , regime]
      # No need to swap the ordering of the df params, as they do not affect the linear IRF.
    }

  } else { # identification == "reduced_form" or "recursive"
    if(identification == "reduced_form") {
      message("Reduced form model supplied, using lower triangular recursive identification.")
    }
    # Shocks identified by lower-triangular Cholesky decomposition
    B_matrix <- t(chol(all_Omega[, , regime]))
  }

  ## Function to calculate IRF
  get_IRF <- function(p, d, N, boldA, B_matrix) {
    J_matrix <- create_J_matrix(d=d, p=p)
    all_boldA_powers <- array(NA, dim=c(d*p, d*p, N+1)) # The first [, , i1] is for the impact period, i1+1 for period i1
    all_Phi_i <- array(NA, dim=c(d, d, N+1)) # JA^iJ' matrices; [, , 1] is for the impact period, i1+1 for period i1
    all_Theta_i <- array(NA, dim=c(d, d, N+1)) # IR-matrices [, , 1] is for the impact period, i1+1 for period i1
    for(i1 in 1:(N + 1)) { # Go through the periods, i1=1 for the impact period, i1+1 for the period i1 after the impact
      if(i1 == 1) {
        all_boldA_powers[, , i1] <- diag(d*p) # Diagonal matrix for the power 0
      } else {
        all_boldA_powers[, , i1] <- all_boldA_powers[, , i1 - 1]%*%boldA # boldA^{i1-1} because i1=1 is for the zero period
      }
      all_Phi_i[, , i1] <- tcrossprod(J_matrix%*%all_boldA_powers[, , i1], J_matrix)
      all_Theta_i[, , i1] <- all_Phi_i[, , i1]%*%B_matrix
    }
    all_Theta_i # all_Theta_i[variable, shock, horizon] -> all_Theta_i[variable, shock, ] subsets the IRF!
  }

  ## Calculate the impulse response functions:
  point_est <- get_IRF(p=p, d=d, N=N,
                       boldA=all_boldA[, , regime], # boldA= Companion form AR matrix of the selected regime
                       B_matrix=B_matrix)

  dimnames(point_est)[[1]] <- colnames(stvar$data)
  dimnames(point_est)[[2]] <- paste("Shock", 1:d)
  # all_Theta_i[variable, shock, horizon] -> all_Theta_i[variable, shock, ] subsets the IRF!

  ## Fixed design wild residual bootstrap for calculating confidence bounds
  if(!is.null(ci) && !ci_possible) {
    warning("Confidence bounds are not available as the autoregressive dynamics are not linear")
    all_bootstrap_IRF <- NULL
  } else if(!is.null(ci) && ci_possible) { # Bootstrap confidence bounds
    ## Create initial values for the two-phase estimation algorithm: does not vary across the bootstrap reps
    new_params <- stvar$params

    # For all models, bootstrapping conditions on the estimated transition weight parameters,
    # so they need to be removed:
    if(is.null(weight_constraints) && M > 1 && weight_function != "exogenous") { # No weight constraints
      n_weightpars <- length(weightpars) - ifelse(weight_function == "relative_dens", 1, 0)
      new_params <- c(new_params[1:(length(new_params) - n_weightpars - length(distpars))],
                      distpars) # Removes weight params
    } else if(M > 1 && weight_function != "exogenous") { # Weight constraints employed
      if(weight_constraints$R != 0) { # Linear weight constraints (no changes if R == 0)
        new_params <- c(new_params[1:(length(new_params) - nrow(weight_constraints$R) - length(distpars))],
                        distpars) # Removes weight params
      }
    } # If M==1 or weight_function=="exogenous", no changes

    # Set the new fixed weight constraints to be the originally fitted weights params
    if(weight_function == "exogenous") {
      new_weight_constraints <- NULL
    } else if(weight_function == "relative_dens") {
      new_weight_constraints <- list(R=0, r=weightpars[-length(weightpars)]) # Removes alpha_M
    } else {
      new_weight_constraints <- list(R=0, r=weightpars)
    }

    # For models identified by heteroskedasticity, bootstrapping also conditions on the estimated lambda parameters
    # to keep the shocks in a fixed ordering (which is given for recursively identified models). Also make sure that
    # each column of W has a strict sign constraint: if not normalize diagonal elements to positive.
    if(identification == "heteroskedasticity") {
      W <- pick_W(p=p, M=M, d=d, params=params, identification=identification)
      all_lambdas <- pick_lambdas(p=p, M=M, d=d, params=params, identification=identification)
      new_fixed_lambdas <- all_lambdas # If fixed_lambdas already used, they don't change
      if(!is.null(B_constraints)) {
        new_B_constraints <- B_constraints
        for(i1 in 1:ncol(new_B_constraints)) { # Iterate through each column of W
          col_vec <- B_constraints[, i1]
          if(all(is.na(col_vec) | col_vec == 0, na.rm=TRUE)) { # Are all elements in col_vec NA or zero?
            if(is.na(new_B_constraints[i1, i1])) { # Check if the diagonal element is NA
              new_B_constraints[i1, i1] <- 1 # Impose a positive sign constraints to the diagonal
            } else { # Zero constraint in the diagonal elements
              # Impose positive sign constraint on the first non-zero element
              new_B_constraints[which(is.na(col_vec))[1], i1] <- 1
            }
          }
        }
      } else { # No B_constraints
        new_B_constraints <- matrix(NA, nrow=d, ncol=d)
        diag(new_B_constraints) <- 1 # Impose positive sign constraints to the diagonal
      }
      other_constraints <- list(fixed_lambdas=new_fixed_lambdas) # other_constraints if used internally only

      # Finally, we need to make sure that the W params in new_params are in line with the constraints new_B_constraints.
      # This amounts checking the strict sign constraints and swapping the signs of the columns that don't
      # match the sign constraints.
      if(sum(W == 0, na.rm=TRUE) != sum(B_constraints == 0, na.rm=TRUE)) {
        # Throws an error since Wvec wont work properly if W contains exact zeros that are not constrained to zeros.
        stop(paste("A parameter value in W exactly zero but not constrained to zero.",
                   "Please adjust B_constraints so that the exact zeros match the constraints"))
      }
      # Determine which columns to swap: compare the first non-NA and non-zero element of the column of new_B_constraints
      # to the corresponding element of the corresponding column of W, and swap the signs of the column if the signs don't match.
      for(i1 in 1:ncol(new_B_constraints)) { # Loop through the columns
        col_new_W <- new_B_constraints[,i1]
        col_old_W <- W[,i1]
        which_to_compare <- which(!is.na(col_new_W) & col_new_W != 0)[1] # The first element that imposes a sign constraint
        if(sign(col_new_W[which_to_compare]) != sign(col_old_W[which_to_compare])) { # Different sign than the constrained one
          W[,i1] <- -W[,i1] # Swap the signs of the column
        }
      }
      # New params with the new W that corresponds to new_B_constraints. Note that AR parameters are assumed
      # identical across the regimes here.
      new_params <- c(stvar$params[1:(d + p*d^2)], Wvec(W), distpars) # No lambdas or weightpars
    } else if(identification == "non-Gaussianity") {
      # For models identified by non-Gaussianity, we keep the shocks in a fixed ordering by assuming that the
      # first non-zero element of each column of the impact matrix of Regime 1 is strictly positive and they are
      # in a decreasing order. This already done in the "new_stvar" object above.

      # Mark that the param space should be constrained to the fixed signs and ordering of the columns of B_1:
      other_constraints <- list(B1_constraints="fixed_sign_and_order")

      # Obtain the new parameters and B_constraints:
      new_params <- new_stvar$params
      new_B_constraints <- new_stvar$model$B_constraints
      # Note that we don't impose sign constraints via B_constraints here, which would impose them to
      # all B_m. Instead, we impose the sign constraints to the first regime's B matrix and restrict
      # the parameter space by other_constraints.
    } else { # No changes
      new_B_constraints <- B_constraints
      other_constraints <- NULL
    }

    ## Obtain residuals
    # Each y_t fixed, so the initial values y_{-p+1},...,y_0 are fixed in any case.
    # For y_1,...,y_T, new residuals are drawn at each bootstrap rep.
    # First, obtain the original residuals:
    mu_mt <- stvar$regime_cmeans[, , regime]
    u_t <- stvar$residuals_raw

    # Function to get one boostrapped IRF
    get_one_bootstrap_IRF <- function(seed) { # Take rest of the arguments from parent environment
      set.seed(seed) # Set seed for data generation

      # Create new data
      eta_t <- sample(c(-1, 1), size=nrow(u_t), replace=TRUE, prob=c(0.5, 0.5))
      new_resid <- eta_t*u_t # each row of u_t multiplied by -1 or 1 based on eta_t
      new_data <- rbind(data[1:p,], # Fixed initial values
                        mu_mt + new_resid) # Bootstrapped data

      # Estimate the model to the new data
      bs_params <- fitbsSSTVAR(data=new_data, p=p, M=M, params=new_params,
                               weight_function=weight_function, weightfun_pars=weightfun_pars,
                               cond_dist=cond_dist, parametrization=parametrization,
                               identification=identification, AR_constraints=AR_constraints,
                               mean_constraints=mean_constraints, weight_constraints=new_weight_constraints,
                               B_constraints=new_B_constraints, other_constraints=other_constraints,
                               seed=seed)

      # Get the IRF from the bootstrap replication
      tmp_params <- reform_constrained_pars(p=p, M=M, d=d, params=bs_params,
                                            weight_function=weight_function, weightfun_pars=weightfun_pars,
                                            cond_dist=cond_dist, identification=identification,
                                            AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                            weight_constraints=new_weight_constraints, B_constraints=new_B_constraints,
                                            other_constraints=other_constraints)
      tmp_all_A <- pick_allA(p=p, M=M, d=d, params=tmp_params)
      tmp_all_boldA <- form_boldA(p=p, M=M, d=d, all_A=tmp_all_A)
      tmp_all_Omega <- pick_Omegas(p=p, M=M, d=d, params=tmp_params, cond_dist=cond_dist, identification=identification)

      # Obtain the impact matrix of the regime the IRF is to calculated for
      if(identification == "heteroskedasticity") { # Shocks identified by heteroskedasticity
        tmp_W <- pick_W(p=p, M=M, d=d, params=tmp_params, identification=identification)
        if(regime > 1) { # Include lambdas
          tmp_lambdas <- matrix(pick_lambdas(p=p, M=M, d=d, params=tmp_params, identification=identification), nrow=d, byrow=FALSE)
          tmp_Lambda_m <- diag(tmp_lambdas[, regime - 1])
          tmp_B_matrix <- W%*%sqrt(tmp_Lambda_m)
        } else { # regime == 1
          tmp_B_matrix <- tmp_W
        }
      } else if(identification == "non-Gaussianity") { # ind_Student models always here
        tmp_B_matrix <- tmp_all_Omega[, , regime] # impact matrix readily parametrized
      } else { # identification == "reduced_form" or "recursive"
        tmp_B_matrix <- t(chol(tmp_all_Omega[, , regime])) # Shocks identified by lower-triangular Cholesky decomposition
      }

      # Calculate and return the IRF
      get_IRF(p=p, d=d, N=N, boldA=tmp_all_boldA[, , regime], B_matrix=tmp_B_matrix)
    }

    ## Calculate the bootstrap replications using parallel computing
    if(ncores > parallel::detectCores()) {
      ncores <- parallel::detectCores()
      message("ncores was set to be larger than the number of cores detected")
    }
    message(paste("Using", ncores, "cores for", bootstrap_reps, "bootstrap replications..."))
    cl <- parallel::makeCluster(ncores)
    on.exit(try(parallel::stopCluster(cl), silent=TRUE)) # Close the cluster on exit, if not already closed.
    parallel::clusterExport(cl, ls(environment(fitSTVAR)), envir=environment(fitSTVAR)) # assign all variables from package:sstvars
    parallel::clusterEvalQ(cl, c(library(pbapply), library(Rcpp), library(RcppArmadillo), library(sstvars)))
    set.seed(seed); seeds <- sample.int(1e+6, size=bootstrap_reps, replace=TRUE) # Seeds for the bootstrap replications
    all_bootstrap_IRF <- pbapply::pblapply(1:bootstrap_reps, function(i1) get_one_bootstrap_IRF(seed=seeds[i1]), cl=cl)
    parallel::stopCluster(cl=cl)
  } else {
    all_bootstrap_IRF <- NULL
  }

  ## Accumulate IRF based on which_cumulative
  if(length(which_cumulative) > 0) {
    for(which_var in which_cumulative) {
      # Accumulate the impulse responses of the variables in which_cumulative
      point_est[which_var, , ] <- t(apply(point_est[which_var, , , drop=FALSE], MARGIN=2, FUN=cumsum))
      if(!is.null(all_bootstrap_IRF)) { # Do the same accumulation for each bootstrap replication:
        for(i2 in 1:length(all_bootstrap_IRF)) {
          all_bootstrap_IRF[[i2]][which_var, , ] <- t(apply(all_bootstrap_IRF[[i2]][which_var, , , drop=FALSE],
                                                            MARGIN=2, FUN=cumsum))
        }
      }
    }
  }

  ## Scale the IRFs
  if(!is.null(scale)) {
    for(i1 in 1:ncol(scale)) {
      which_shock <- scale[1, i1]
      which_var <- scale[2, i1]
      scale_size <- scale[3, i1]
      # Scale the IRFs of which_shock to correspond scale_size impact response of the variable which_var:
      multiplier <- scale_size/point_est[which_var, which_shock, 1]
      point_est[, which_shock, ] <- multiplier*point_est[, which_shock, ] # Impact response to scale_size
      if(!is.null(all_bootstrap_IRF)) { # Do the same scaling for each bootstrap replication:
        for(i2 in 1:length(all_bootstrap_IRF)) {
          multiplier <- scale_size/all_bootstrap_IRF[[i2]][which_var, which_shock, 1]
          all_bootstrap_IRF[[i2]][, which_shock, ] <- multiplier*all_bootstrap_IRF[[i2]][, which_shock, ]
        }
      }
    }
  }

  ## Calculate the confidence bounds
  if(!is.null(all_bootstrap_IRF)) {
    # First we convert the list into a 4D array in order to use apply to calculate empirical qunatiles
    all_bootstrap_IRF_4Darray <- array(NA, dim = c(dim(all_bootstrap_IRF[[1]]), length(all_bootstrap_IRF)))
    for (i1 in 1:length(all_bootstrap_IRF)) { # Fill the arrays
      all_bootstrap_IRF_4Darray[, , , i1] <- all_bootstrap_IRF[[i1]]
    }

    # Calculate empirical quantiles to obtain ci
    quantile_fun <- function(x) quantile(x, probs=c((1 - ci)/2, 1 - (1 - ci)/2), na.rm=TRUE)
    conf_ints <- apply(all_bootstrap_IRF_4Darray, MARGIN=1:3, FUN=quantile_fun)
    conf_ints <- aperm(conf_ints, perm=c(2, 3, 4, 1))
    dimnames(conf_ints)[[1]] <- colnames(data)
    dimnames(conf_ints)[[2]] <- paste("Shock", 1:d)
  } else {
    conf_ints <- NULL
    all_bootstrap_IRF_4Darray <- NULL
  }

  # Return the results
  structure(list(point_est=point_est,
                 conf_ints=conf_ints,
                 all_bootstrap_reps=all_bootstrap_IRF_4Darray,
                 N=N,
                 ci=ci,
                 scale=scale,
                 which_cumulative=which_cumulative,
                 seed=seed,
                 stvar=stvar),
            class="irf")
}
