#' @title Estimate linear impulse response function based on a single regime of a structural STVAR model.
#'
#' @description \code{linear_IRF} estimates linear impulse response function based on a single regime
#'   of a structural STVAR model.
#'
#' @param stvar an object of class \code{'stvar'} defining a structural or reduced form
#'   stvar model. For a reduced form model, the shocks are automatically identified by
#'   the lower triangular Cholesky decomposition.
#' @param N a positive integer specifying the horizon how far ahead should the
#'   linear impulse responses be calculated.
#' @param regime Based on which regime the linear IRF should be calculated?
#'   An integer in \eqn{1,...,M}.
#' @param which_cumulative a numeric vector with values in \eqn{1,...,d}
#'   (\code{d=ncol(data)}) specifying which the variables for which the linear impulse
#'   responses should be cumulative. Default is none.
#' @param scale should the linear IRFs to some of the shocks be scaled so that they
#'   correspond to a specific magnitude of instantaneous or peak response
#'   of some specific variable (see the argument \code{scale_type})?
#'   Provide a length three vector where the shock of interest
#'   is given in the first element (an integer in \eqn{1,...,d}), the variable of
#'   interest is given in the second element (an integer in \eqn{1,...,d}), and
#'   the magnitude of its instantaneous or peak response in the third element
#'   (a non-zero real number). If the linear IRFs of multiple shocks should be scaled,
#'   provide a matrix which has one column for each of the shocks with the columns being
#'   the length three vectors described above.
#' @param scale_type If argument \code{scale} is specified, should the linear IRFs be
#'   scaled to match an instantaneous response (\code{"instant"}) or peak response
#'   (\code{"peak"}). If \code{"peak"}, the scale is based on the largest magnitude
#'   of peak response in absolute value. Ignored if \code{scale} is not specified.
#' @param scale_horizon If \code{scale_type == "peak"} what the maximum horizon up
#'   to which peak response is expected? Scaling won't based on values after this.
#' @param ci a numeric vector with elements in \eqn{(0, 1)} specifying the
#'   confidence levels of the confidence intervals calculated via a bootstrap
#'   method, see the details section. Available only for models that impose linear
#'   autoregressive dynamics (excluding changes in the volatility regime).
#' @param seed a length one numeric vector initializing the seed for the random generator.
#' @param ... parameters passed to the plot method \code{plot.irf} that plots
#'   the results.
#' @details  FILL IN DETAILS! NOTE THAT THIS FUNCTION IS STILL UNDER CONSTRUCTION!
#' @return Returns a class \code{'irf'} list with the linear IRFs in ... FILL IN!
#' @seealso \code{\link{GIRF}}, \code{\link{GFEVD}}, \code{\link{fitSTVAR}}, \code{\link{STVAR}},
#'   \code{\link{reorder_W_columns}}, \code{\link{swap_W_signs}}
#' @references
#'  \itemize{
#'    \item Kilian L. and Lütkepohl H. 2017. Structural Vectors Autoregressive Analysis.
#'          \emph{Cambridge University Press}, Cambridge.
#'  }
#' @examples
#'  # FILL IN
#' @export

linear_IRF <- function(stvar, N=30, regime=1, which_cumulative=numeric(0),
                       scale=NULL, scale_type=c("instant", "peak"), scale_horizon=N,
                       ci=NULL, seed=NULL, ...) {
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
  boldA <- all_boldA[, , regime] # Companion form AR matrix for the selected regime
  J_matrix <- create_J_matrix(d=d, p=p)
  all_boldA_powers <- array(NA, dim=c(d*p, d*p, N+1)) # The first [, , i1] is for the impact period, i1+1 for period i1
  all_phi_i <- array(NA, dim=c(d, d, N+1)) # JA^iJ' matrices; [, , 1] is for the impact period, i1+1 for period i1
  all_Phi_i <- array(NA, dim=c(d*p, d*p, N+1)) # IR-matrices [, , 1] is for the impact period, i1+1 for period i1
  for(i1 in 1:(N + 1)) { # Go through the periods, i1=1 for the impact period, i1+1 for the period i1 after the impact
    if(i1 == 1) {
      all_boldA_powers[, , i1] <- diag(d*p) # Diagonal matrix for the power 0
    } else {
      all_boldA_powers[, , i1] <- all_boldA_powers[, , i1 - 1]%*%boldA # boldA^{i1-1} because i1=1 is for the zero period
    }
    all_phi_i[, , i1] <- J_matrix%*%all_boldA_powers[, , i1]%*%t(J_matrix)
    all_Phi_i[, , i1] <- all_phi_i[, , i1]%*%B_matrix
  }
  # all_Phi_i[variable, shock, horizon] -> all_Phi_i[variable, shock, ] subsets the IRF!

  # NOTE which cumulative does not do anything yet! Maybe original + cumulative for all?
  # Depends on how the confidence bounds are created!

  # Return the results
  structure(list(all_irfs=all_Phi_i,
                 N=N,
                 ci=ci,
                 which_cumulative=which_cumulative,
                 seed=seed,
                 stvar=stvar),
            class="irf")
}
