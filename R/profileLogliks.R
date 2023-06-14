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
#' \donttest{
#' # Running all the below examples takes approximately FILL IN HOW MANY MINUTES
#'
#' # FILL IN
#' }
#' @export

profile_logliks <- function(stvar, which_pars, scale=0.02, nrows, ncols, precision=200,
                            stab_tol=0.001, posdef_tol=1e-08, distpar_tol=1e-08, weightpar_tol=1e-08) {
  # Initial checks
  stopifnot(class(stvar) == "stvar")
  if(is.null(stvar$data)) stop("Cannot plot profile logliks if the model does not contain data.")

  # Model specs
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  params <- stvar$params
  weight_function <- stvar$model$weight_function
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function, weightfun_pars=stvar$model$weightfun_pars)
  cond_dist <- stvar$model$cond_dist
  parametrization <- stvar$model$parametrization
  identification <- stvar$model$identification
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  weight_constraints <- stvar$model$weight_constraints
  B_constraints <- stvar$model$B_constraints

  # Checks, default arguments etc
  if(identification != "reduced_form") stop("Structural models are not yet implemented to profile_logliks")
  if(!is.null(B_constraints)) stop("B_constained models are not yet implemented to profile_logliks")
  if(missing(which_pars)) which_pars <- 1:length(params)
  if(!all_pos_ints(which_pars) || any(which_pars > length(params))) {
    stop("The argument 'which_pars' should contain strictly positive integers not larger than length of the parameter vector.")
  } else if(anyDuplicated(which_pars) != 0) {
    stop("There are dublicates in which_pars")
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
  } else { # Means constrained
    n_mean_pars <- d*length(mean_constraints)
  }
  if(is.null(AR_constraints)) {
    n_ar_pars <- M*p*d^2
  } else { # AR matrices constrained
    n_ar_pars <- ncol(AR_constraints)
  }
  if(is.null(B_constraints)) {
    n_covmat_pars <- M*d*(d + 1)/2
  } else { # Constraints on the impact matrix
    stop("B_constraints not yet implemented to n_parms")
    n_covmat_pars <- NULL
  }
  if(is.null(weight_constraints)) {
    if(weight_function == "relative_dens" || weight_function == "threshold") {
      n_weight_pars <- M - 1
    } else if(weight_function == "logistic" || weight_function == "exponential") {
      n_weight_pars <- 2
    } else if(weight_function == "mlogit") {
      n_weight_pars <- (M - 1)*(1 + length(weightfun_pars[[1]])*weightfun_pars[[2]])
    } else {
      stop("Unknown weightfunction in n_params")
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
  } else { # cond_dist == "Student"
    n_dist_pars <- 1 # degrees of freedom param
  }

  for(i1 in which_pars) { # Go though the parameters
    pars <- params
    range <- abs(scale*pars[i1])
    vals <- seq(from=pars[i1] - range, to=pars[i1] + range, length.out=precision) # Loglik to be evaluated at these points
    logliks <- vapply(vals, function(val) { # Log-likelihoods about the estimate
      new_pars <- pars
      new_pars[i1] <- val # Change the single parameter value
      loglikelihood(data=stvar$data, p=p, M=M, params=new_pars, weight_function=weight_function,
                    weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                    parametrization=parametrization, identification=identification,
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=B_constraints,
                    to_return="loglik", check_params=TRUE, minval=NA,
                    stab_tol=stab_tol, posdef_tol=posdef_tol, distpar_tol=distpar_tol, weightpar_tol=weightpar_tol)
    }, numeric(1))

    # Determine which type of parameter is i1 to determine the label for the individual plot
  }
}
