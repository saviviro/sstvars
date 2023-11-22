#' @title Maximum likelihood estimation of a reduced form or structural STVAR model based on preliminary estimates
#'
#' @description \code{iterate_more} uses a variable metric algorithm to estimate a reduced form or structural STVAR model
#'  (object of class \code{'stvar'}) based on preliminary estimates.
#'
#' @inheritParams fitSTVAR
#' @inheritParams STVAR
#' @param stvar an object of class \code{'stvar'}, created by, e.g., \code{fitSTVAR}.
#' @details The purpose of \code{iterate_more} is to provide a simple and convenient tool to finalize
#'   the estimation when the maximum number of iterations is reached when estimating a STVAR model
#'   with the main estimation function \code{fitSTVAR} or \code{fitSSTVAR}.
#' @return Returns an object of class \code{'stvar'} defining the estimated model
#' @seealso \code{\link{fitSTVAR}}, \code{\link{STVAR}}, \code{\link[stats]{optim}}
#' @inherit STVAR references
#' @examples
#' \donttest{
#' ## These are long running examples that use parallel computing!
#' ## Running the below examples takes approximately FILL IN HOW MANY minutes
#'
#' # p=1, M=2, d=2 relative_dens Gaussian STVAR, only 5 iterations of
#' # the variable matrix algorithm.
#' fit12 <- fitSTVAR(gdpdef, p=1, M=2, nrounds=1, seeds=1, maxit=5, use_parallel=FALSE)
#' fit12
#'
#' # Iterate more:
#' fit12 <- iterate_more(fit12)
#' fit12
#' }
#' @export

iterate_more <- function(stvar, maxit=100, calc_std_errors=TRUE) {
  check_stvar(stvar)
  stopifnot(maxit %% 1 == 0 & maxit >= 1)
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  data <- stvar$data
  weight_function <- stvar$model$weight_function
  cond_dist <- stvar$model$cond_dist
  parametrization <- stvar$model$parametrization
  identification <- stvar$model$identification
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  weight_constraints <- stvar$model$weight_constraints
  B_constraints <- stvar$model$B_constraints
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function,
                                         weightfun_pars=stvar$model$weightfun_pars, cond_dist=cond_dist)
  minval <- get_minval(stvar$data)
  npars <- length(stvar$params)

  # Function to optimize
  loglik_fn <- function(params) {
    tryCatch(loglikelihood(data=data, p=p, M=M, params=params,
                           weight_function=weight_function, weightfun_pars=weightfun_pars,
                           cond_dist=cond_dist, parametrization=parametrization,
                           identification=identification, AR_constraints=AR_constraints,
                           mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                           B_constraints=B_constraints, to_return="loglik", check_params=TRUE, minval=minval),
             error=function(e) minval)
  }

  # Gradient of the log-likelihood function using central difference approximation
  h <- 6e-6
  I <- diag(rep(1, times=npars))
  loglik_grad <- function(params) {
    vapply(1:npars, function(i1) (loglik_fn(params + I[i1,]*h) - loglik_fn(params - I[i1,]*h))/(2*h), numeric(1))
  }

  res <- optim(par=stvar$params, fn=loglik_fn, gr=loglik_grad, method=c("BFGS"), control=list(fnscale=-1, maxit=maxit))
  if(res$convergence == 1) message("The maximum number of iterations was reached! Consired iterating more.")

  ret <- STVAR(data=data, p=p, M=M, d=d, params=res$par,
               weight_function=weight_function,
               weightfun_pars=weightfun_pars,
               cond_dist=cond_dist,
               parametrization=parametrization,
               identification=identification,
               AR_constraints=AR_constraints,
               mean_constraints=mean_constraints,
               weight_constraints=weight_constraints,
               B_constraints=B_constraints,
               calc_std_errors=calc_std_errors)

  ret$all_estimates <- stvar$all_estimates
  ret$all_logliks <- stvar$all_logliks
  ret$which_converged <- stvar$which_converged
  if(!is.null(stvar$which_round)) {
    ret$which_round <- stvar$which_round
    ret$all_estimates[[stvar$which_round]] <- ret$params
    ret$all_logliks[stvar$which_round] <- ret$loglik
    ret$which_converged[stvar$which_round] <- res$convergence == 0
  }
  warn_eigens(ret)
  ret
}



#' @title Maximum likelihood estimation of a structural STVAR model based on preliminary estimates from
#'   a reduced form model.
#'
#' @description \code{fitSSTVAR} uses a variable metric algorithm to estimate a structural STVAR model
#'   based on preliminary estimates from a reduced form model.
#'
#' @inheritParams fitSTVAR
#' @inheritParams STVAR
#' @param maxit_robust the maximum number of iterations on the first phase robust estimation, if employed.
#' @param stvar a an object of class \code{'stvar'}, created by, e.g., \code{fitSTVAR},
#'   specifying a reduced form or a structural model
#' @param robust_method Should some robust estimation method be used in the estimation before switching
#'   to the gradient based variable metric algorithm? See details.
#'
#' @details When the structural model does not impose overidentifying constraints, it is directly
#'   obtained from the reduced form model, and estimation is not required. When overidentifying constraints
#'   are imposed, the model is estimated via ..
#'
#'   Structural models can be provided in the argument \code{stvar} if overidentifying constraints should be
#'   imposed.
#'
#'   Using a robust estimation method before switching to the variable metric can be useful if the initial
#'   estimates are not very close to the ML estimate of the structural model, as the variable metric algorithm
#'   (usually) converges to a nearby local maximum or saddle point. This is particularly the case when the imposed
#'   overidentifying restrictions are such that the unrestricted estimate is not close to satisfying them.
#'
#'   Employs the estimation function \code{optim} from the package \code{stats} that implements the optimization
#'   algorithms. See \code{?optim} for the documentation on the
#' @return Returns an object of class \code{'stvar'} defining the structural model
#' @seealso \code{\link{fitSTVAR}}, \code{\link{STVAR}}, \code{\link[stats]{optim}}
#' @inherit STVAR references
#' @examples
#' \donttest{
#' ## These are long running examples that use parallel computing!
#' ## Running the below examples takes approximately FILL IN HOW MANY minutes
#' }
#' @export

fitSSTVAR <- function(stvar, new_identification=c("recursive", "heteroskedasticity"), new_B_constraints=NULL,
                      maxit=1000, maxit_robust=1000, robust_method=c("none", "Nelder-Mead", "SANN"), calc_std_errors=TRUE) {
  check_stvar(stvar)
  stopifnot(maxit %% 1 == 0 & maxit >= 1)
  robust_method <- match.arg(robust_method)
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  data <- stvar$data
  params <- stvar$params
  weight_function <- stvar$model$weight_function
  cond_dist <- stvar$model$cond_dist
  parametrization <- stvar$model$parametrization
  old_identification <- stvar$model$identification
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  weight_constraints <- stvar$model$weight_constraints
  old_B_constraints <- stvar$model$B_constraints
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function,
                                         weightfun_pars=stvar$model$weightfun_pars, cond_dist=cond_dist)
  minval <- get_minval(stvar$data)
  old_npars <- length(stvar$params) # old npars, B_constraints not adjusted!

  # Check identification
  if(old_identification != "reduced_form") {
    if(old_identification == "recursive") {
      warning(paste("Recursively identified model supplied (not possible to impose overidentifying restrictions,",
                    "so the original model is returned"))
      return(stvar)
    } else if(old_identification != new_identification) {
      # Old model structural but identification changes
      stop(paste("The type of the identification of a structural models cannot be changed.",
                 "Consider using the initial reduced model to obtain the structural model with a specific identification."))
    }
  } else if(new_identification == "recursive") {
    # Old identification reduced form, recursively identified model should be returned.
    return(STVAR(data=data, p=p, M=M, d=d, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
                 parametrization=parametrization, identification=new_identification,
                 AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                 weight_constraints=weight_constraints, B_constraints=NULL, calc_std_errors=calc_std_errors))
  }

  identification <- new_identification # old in old_identification
  B_constraints <- new_B_constraints # old in old_B_constraints

  # Check the new B_constraints
  check_constraints(p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    parametrization=parametrization, identification=identification,
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=B_constraints)

  ## The numbers of parameters (needed to create the new initial param vector)
  if(is.null(mean_constraints)) {
    n_mean_pars <- M*d
  } else { # Means constrained
    n_mean_pars <- d*length(mean_constraints)
  }
  if(is.null(AR_constraints)) {
    n_ar_pars <- M*p*d^2
  } else { # AR matrices constrained
    n_ar_pars <- ncol(AR_constraints)
  }
  # Covmatpars
  if(old_identification %in% c("reduced_form", "recursive")) {
    n_covmat_pars <- M*d*(d + 1)/2 # No B_constraints available here
    n_zeros_in_W <- 0
  } else { # identification == "heteroskedasticity
    if(is.null(old_B_constraints)) {
      n_zeros_in_W <- 0
    } else {
      n_zeros_in_W <- sum(old_B_constraints == 0, na.rm=TRUE)
    }
    n_covmat_pars <- d^2 + d*(M - 1) - n_zeros_in_W
  }
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
  if(cond_dist == "Gaussian") {
    n_dist_pars <- 0
  } else { # cond_dist == "Student
    n_dist_pars
  }
  n_pars <- n_mean_pars + n_ar_pars + n_covmat_pars + weight_pars + n_dist_pars


  # Set up initial values using the obtained parameters
  if(identification == "heteroskedasticity") {
    if(M == 1) {
      stop("Identification by heteroskedasticity requires at least two regimes")
    } else if(M == 2 && is.null(B_constraints) && old_identification == "reduced_form") {
      # Nothing to estimate, decompose the error term covariance matrices and create new params
      # and return the new model
      return(get_hetsked_sstvar(stvar, calc_std_errors=calc_std_errors))
    } else if(old_B_constraints == B_constraints && old_identification == "heteroskedasticity") {
      # The model does not change, return the original model
      return(stvar)
    }
    ## If arrived here, overidentifying constraints are imposed, either via M>2 or B_constraints.
    if(old_identification == "reduced_form") {
      if(M == 2) {
        params <- get_hetsked_sstvar(stvar, calc_std_errors=FALSE)$params # Change params to hetsked params
      } else {
        # M > 2, so needs to obtain initial estimates some other way. We decompose the first two Omegas
        # to obtain W and the first lambdas; then we obtain the lambda_m by decomposing the Omega_1 and Omega_m.
        # Maybe not the optimal solution, but easily scales to any number of regimes and gives the most weight
        # to the first regime, and the second most to the third regime.

        # First, pick reduced form error covariance matrices:
        params_std <- reform_constrained_pars(p=p, M=M, params=params, cond_dist=cond_dist,
                                              weight_function=weight_function, weightfun_pars=weightfun_pars,
                                              parametrization=parametrization, identification=old_identification,
                                              AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                              weight_constraints=weight_constraints, B_constraints=old_B_constraints)
        all_Omegas <- pick_Omegas(p=p, M=M, params=params_std, identification=old_identification)

        # W and lambda_2:
        W_and_lambdas <- numeric(d^2 + (M - 1)*d)
        W_and_lambdas[1:(d^2 + d)] <- diag_Omegas(Omega1=all_Omegas[, , 1], Omega2=all_Omegas[, , 2])

        # lambda_3,...,lambda_M:
        for(i1 in 3:M) {
          tmp <- diag_Omegas(Omega1=all_Omegas[, , 1], Omega2=all_Omegas[, , i1])
          W_and_lambdas[(d^2 + d*(i1 - 2) + 1):(d^2 + d*(i1 - 1))] <- tmp[(d^2 + 1):(d^2 + d)]
        }
      }
    }

    # Obtain the W pars
    if(old_identification == "reduced_form" && M > 2) {
      oldWpars <- W_and_lambdas[1:d^2]
    } else {
      oldWpars <- params[(n_mean_pars + n_ar_pars + 1):(n_mean_pars + n_ar_pars + d^2 - n_zeros_in_W)]
    }

    n_sign_changes <- 0
    oldWpars_inW <- matrix(0, nrow=d, ncol=d) # Fill in the parameters including non-parametrized zeros:
    if(is.null(old_B_constraints)) {
      old_B_constraints <- matrix(NA, nrow=d, ncol=d) # No constraints imposed here but avoids null problems below
    }
    oldWpars_inW[old_B_constraints != 0 | is.na(old_B_constraints)] <- oldWpars

    # Go through the matrices old_B_constraints, new_B_constraints, and changes the oldWpars accordingly
    newWpars <- matrix(NA, nrow=d, ncol=d)
    for(i1 in 1:d) { # i1 marks the row
      for(i2 in 1:d) { # i2 marks the column
        so <- sign(old_B_constraints[i1, i2])
        sn <- sign(new_B_constraints[i1, i2])
        if(is.na(so)) { # If no constraint, just handle it as if there was sign constraint of the estimate's sign
          so <- sign(oldWpars_inW[i1, i2])
        }
        if(is.na(sn) && so != 0) { # No new constraint, old constraint is not zero
          newWpars[i1, i2] <- oldWpars_inW[i1, i2] # Insert old param
        } else if(is.na(sn) && so == 0) { # No new constraints, old constraint is zero
          newWpars[i1, i2] <- 1e+6  # Insert close to zero value: not exactly zero to avoid the parameter disappearing in Wvec
        } else if(so == sn) { # Same constraint in new and old W
          newWpars[i1, i2] <- oldWpars_inW[i1, i2] # Insert old param
        } else if(sn == 0) { # Zero constraint in new W
          newWpars[i1, i2] <- 0 # Insert zero (which will be removed as it is not parametrized)
        } else if(so > sn) { # sn must be negative, so could be zero or positive
          newWpars[i1, i2] <- -0.01 # Insert small negative number
          n_sign_changes <- n_sign_changes + 1
        } else { # It must be that so < sn, which implies sn > 0, while so could be zero or negative
          newWpars[i1, i2] <- 0.01 # Insert s mall positive number
          n_sign_changes <- n_sign_changes + 1
        }
      }
    }
    newWpars <- Wvec(newWpars) # New initial parameters for W

    # New initial params
    if(old_identification == "reduced_form" && M > 2) {
      new_params <- c(params[1:(n_mean_pars + n_ar_pars)], newWpars, W_and_lambdas[(d^2 + 1):(d^2 + (M - 1)*d)],
                      params[(n_mean_pars + n_ar_pars + n_covmat_pars + 1):length(params)])
    } else { # Het.sked model already originally or M=2
      new_params <- c(params[1:(n_mean_pars + n_ar_pars)], newWpars,
                      params[(n_mean_pars + n_ar_pars + length(oldWpars) + 1):length(params)])
    }
    if(n_sign_changes > 0) {
      cat(paste0("There was ", n_sign_changes, "sign changes in W when creating preliminary estimates for the new model.",
                  "The sign changes make the results more unrealiable. To obtain more reliable results, consider using",
                  "the function 'swap_W_signs' to create a model that has the sign constraints readily satisfied and",
                  "then applying this function.\n"))
    }
  } else {
    # Not possible to arrive here currently, but will be updated later for new identification methods.
    stop("Unknown identification in fitSSTVAR")
  }

  ### Check ups for the initial estimates of the new model
  new_loglik <- loglikelihood(data=data, p=p, M=M, params=new_params, cond_dist=cond_dist,
                              weight_function=weight_function, weightfun_pars=weightfun_pars,
                              parametrization=parametrization, identification=identification,
                              AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                              weight_constraints=weight_constraints, B_constraints=B_constraints,
                              to_return="loglik", check_params=TRUE, minval=NULL)

  if(is.null(new_loglik)) {
    cat("Problem with the new parameter vector - try a different new_W?\n See:\n")
    check_params(p=p, M=M, params=new_params, cond_dist=cond_dist,
                 weight_function=weight_function, weightfun_pars=weightfun_pars,
                 parametrization=parametrization, identification=identification,
                 AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                 weight_constraints=weight_constraints, B_constraints=B_constraints,)
  }

  cat("The log-likelihood of the supplied model:   ", round(c(stvar$loglik), 3),
      "\nConstrained log-likelihood prior estimation:", round(new_loglik, 3), "\n\n")


  ### Estimation
  if(is.null(data)) {
    # No data, so returns the model with preliminary param values without estimation.
    message("There is no data for the model, so no estimation was done")
    return(STVAR(data=data, p=p, M=M, d=d, params=new_params, cond_dist=cond_dist,
                 weight_function=weight_function, weightfun_pars=weightfun_pars,
                 parametrization=parametrization, identification=identification,
                 AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                 weight_constraints=weight_constraints, B_constraints=B_constraints,
                 calc_std_errors=FALSE))
  }


  # Function to optimize
  loglik_fn <- function(params) {
    tryCatch(loglikelihood(data=data, p=p, M=M, params=new_params,
                           weight_function=weight_function, weightfun_pars=weightfun_pars,
                           cond_dist=cond_dist, parametrization=parametrization,
                           identification=identification, AR_constraints=AR_constraints,
                           mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                           B_constraints=B_constraints, to_return="loglik", check_params=TRUE, minval=minval),
             error=function(e) minval)
  }

  ## Gradient of the log-likelihood function using central difference approximation
  h <- 6e-6
  I <- diag(rep(1, times=npars))
  loglik_grad <- function(params) {
    vapply(1:npars, function(i1) (loglik_fn(params + I[i1,]*h) - loglik_fn(params - I[i1,]*h))/(2*h), numeric(1))
  }

  ## Estimation by robust method
  if(robust_method != "none") {
    #cat(paste0("Estimation with the", robust_method, "algorithm..."), "\n")
    robust_results <- stats::optim(par=new_params, fn=loglik_fn, gr=loglik_grad, method=robust_method,
                                   control=list(fnscale=-1, maxit=maxit_robust))
    new_params <- robust_results$par
    cat(paste0("The log-likelihood after robust estimation: ", round(robust_results$value, 3)), "\n")
  }

  ## Estimation by the variable metric algorithm...
  #cat(paste0("Estimation with the variable metric algorithm..."), "\n")
  final_results <- stats::optim(par=new_params, fn=loglik_fn, gr=loglik_grad, method="BFGS",
                             control=list(fnscale=-1, maxit=maxit))
  new_params <- final_results$par
  cat(paste0("The log-likelihood after the estimation: ", round(final_results$value), 3), "\n")

  #cat("Finished!")
  STVAR(data=data, p=p, M=M, d=d, params=new_params, cond_dist=cond_dist,
        weight_function=weight_function, weightfun_pars=weightfun_pars,
        parametrization=parametrization, identification=identification,
        AR_constraints=AR_constraints, mean_constraints=mean_constraints,
        weight_constraints=weight_constraints, B_constraints=B_constraints,
        calc_std_errors=calc_std_errors)
}


