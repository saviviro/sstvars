#' @title Estimate linear impulse response function based on a single regime of a structural STVAR model.
#'
#' @description \code{linear_IRF} estimates linear impulse response function based on a single regime
#'   of a structural STVAR model.
#'
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
#' @param nrounds on how many estimation rounds should each bootstrap estimation be based on?
#'   Does not have to be very large since initial estimates used are based on the initially fitted model.
#'   Larger number of rounds gives more reliable results but is computationally more demanding.
#' @param seeds a numeric vector of length \code{bootstrap_reps} initializing the seed for the random
#'   generator for each bootstrap replication.
#' @param ... parameters passed to the plot method \code{plot.irf} that plots
#'   the results.
#' @details If the autoregressive dynamics of the model are linear (i.e., either M == 1 or mean and AR parameters
#'   are constrained identical across the regimes), confidence bounds can be calculated based on a fixed-design
#'   wild residual bootstrap method. We employ the method described in Herwartz and Lütkepohl (2014); see also
#'   the relevant chapters in Kilian and Lütkepohl (2017).
#' @return Returns a class \code{'irf'} list with the linear IRFs in ... FILL IN!
#' @seealso \code{\link{GIRF}}, \code{\link{GFEVD}}, \code{\link{fitSTVAR}}, \code{\link{STVAR}},
#'   \code{\link{reorder_W_columns}}, \code{\link{swap_W_signs}}
#' @references
#'  \itemize{
#'    \item Herwartz H. and Lütkepohl H. 2014. Structural vector autoregressions with Markov switching:
#'      Combining conventional with statistical identification of shocks. \emph{Journal of Econometrics},
#'      183, pp. 104-116.
#'    \item Kilian L. and Lütkepohl H. 2017. Structural Vectors Autoregressive Analysis.
#'          \emph{Cambridge University Press}, Cambridge.
#'  }
#' @examples
#'  # FILL IN
#' @export

linear_IRF <- function(stvar, N=30, regime=1, which_cumulative=numeric(0), scale=NULL, ci=NULL,
                       bootstrap_reps=100, ncores=2, nrounds=1, seed=NULL, ...) {
  # Get the parameter values etc
  stopifnot(all_pos_ints(c(N, regime)))
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  params <- stvar$params
  weight_function <- stvar$model$weight_function
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function,
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
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, identification=identification)
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  weightpars <- pick_weightpars(p=p, M=M, d=d, params=params, weight_function=weight_function,
                            weightfun_pars=weightufun_pars, cond_dist=cond_dist)
  distpars <- pick_distpars(params=params, cond_dist=cond_dist)

  # Obtain the impact matrix of the regime the IRF is to calculated for
  if(identification == "heteroskedasticity") { # Shocks identified by heteroskedasticity
    W <- pick_W(p=p, M=M, d=d, params=params, structural_pars=structural_pars)
    if(regime > 1) { # Include lambdas
      lambdas <- matrix(pick_lambdas(p=p, M=M, d=d, params=params, structural_pars=structural_pars), nrow=d, byrow=FALSE)
      Lambda_m <- diag(lambdas[, regime - 1])
      B_matrix <- W%*%sqrt(Lambda_m)
    } else { # regime == 1
      B_matrix <- W
    }
  } else { # identification == "reduced_form" or "recursive"
    if(identification == "reduced_form") {
      message("Reduced form model supplied, using lower triangular recursive identification.")
    }
    # Shocks identified by lower-triangular Cholesky decomposition
    B_matrix <- t(chol(all_Omega[, , regime]))
  }

  # Calculate the impulse response functions
  # Function to calculate IRF
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

  ## Confidence bounds by fixed design wild residual bootstrap
  AR_mats_identical <- all(apply(all_boldA, MARGIN=3, FUN=function(x) identical(x, all_boldA[,,1])))
  means_identical <- !is.null(mean_constraints) && length(mean_constraints) == 1 && all(mean_constraints[[1]] == 1:M)
  ci_possible <- (means_identical && AR_mats_identical) || M == 1


  if(!is.null(ci) && !ci_possible) {
    warning("Confidence bounds are not available as the autoregressive dynamics are not linear")
  } else if(!is.null(ci) && ci_possible) { # Bootstrap confidence bounds
    ## Create initial values for the two-phase estimation algorithm: does not vary across the bootstrap reps
    new_params <- stvar$params

    # For all models, bootstrapping conditions on the estimated transition weight parameters,
    # so they need to be removed:
    if(is.null(weight_constraints) && M > 1) { # No weight constraints
      n_weightpars <- length(weightpars) - ifelse(weight_function=="relative_dens", 1, 0)
      new_params <- c(new_params[1:(length(new_params) - n_weightpars - length(distpars))],
                      distpars) # Removes weight params
    } else { # Weight constraints employed
      if(weight_constraints$R != 0) { # Linear weight constraints (no changes if R == 0)
        new_params <- c(new_params[1:(length(new_params) - nrow(weight_constraints$R) - length(distpars))],
                        distpars) # Removes weight params
      }
    }
    # Set the new fixed weight constraints to be the originally fitted weights params
    if(weight_function == "relative_dens") {
      new_weight_constraints <- weightpars[-length(weightpars)] # Removes alpha_M
    } else {
      new_weight_constraints <- weightpars
    }
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
  } else { # No changes
    new_B_constraints <- B_constraints
    other_constraints <- NULL
  }

  ## Obtain residuals
  # Each y_t fixed, so the initial values y_{-p+1},...,y_0 are fixed in any case.
  # For y_1,...,y_T, new residuals are drawn at each bootstrap rep.
  # First, obtain the original residuals:
  u_t <- stvar$residuals_raw

  get_one_bootstrap_IRF <- function(seed) { # Take rest of the arguments from parent environment
    set.seed(seed) # Set seed for data generation
    estim_seeds <- sample.int(n=1e+6, size=ncalls) # Seeds for estimation

    # Create new data
    eta_t <- sample(c(-1, 1), size=nrow(u_t), replace=TRUE, prob=c(0.5, 0.5))
    new_resid <- eta_t*u_t # each row of u_t multiplied by -1 or 1 based on eta_t
    new_data <- rbind(data[1:p,], # Fixed initial values
                      mu_mt + new_resid) # Bootstrapped data

    ## Estimate the model to the new data

    #### NEED TO FILL IN THE ESTIMATION PART HERE ####
    #### AFTER CREATING THE ESTIMATION BOOTSTRAP FUNCTION ####
  }

  # Return the results
  structure(list(point_est=point_est,
                 N=N,
                 ci=ci,
                 which_cumulative=which_cumulative,
                 seed=seed,
                 stvar=stvar),
            class="irf")
}
