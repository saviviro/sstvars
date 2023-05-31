#' @title Two-phase maximum likelihood estimation of a reduced form smooth transition VAR model
#'
#' @description \code{fitSTVAR} estimates a reduced form smooth transition VAR model in two phases:
#'   in the first phase, it uses a genetic algorithm to find starting values for a gradient based
#'   variable metric algorithm, which it then uses to finalize the estimation in the second phase.
#'   Parallel computing is utilized to perform multiple rounds of estimations in parallel.
#'
#' @inheritParams GAfit
#' @param nrounds the number of estimation rounds that should be performed.
#' @param ncores the number CPU cores to be used in parallel computing.
#' @param maxit the maximum number of iterations in the variable metric algorithm.
#' @param seeds a length \code{nrounds} vector containing the random number generator seed
#'  for each call to the genetic algorithm, or \code{NULL} for not initializing the seed.
#' @param print_res should summaries of estimation results be printed?
#' @param use_parallel employ parallel computing?
#' @param filter_estimates should the likely inappropriate estimates be filtered? See details.
#' @param ... additional settings passed to the function \code{GAfit} employing the genetic algorithm.
#' @details
#'  If you wish to estimate a structural model, estimate first the reduced form model and then use the
#'  use the function \code{fitSSTVAR} to estimate the structural model based on the estimated reduced form model.
#'
#'  Because of complexity and high multimodality of the log-likelihood function, it is \strong{not certain}
#'  that the estimation algorithm will end up in the global maximum point. It is expected that most of the
#'  estimation rounds will end up in some local maximum or saddle point instead. Therefore, a (sometimes large)
#'  number of estimation rounds is required for reliable results. Because of the nature of the model,
#'  the estimation may fail especially in the cases where the number of regimes is chosen too large.
#'
#'  The estimation process is computationally heavy and it might take considerably long time for large models with
#'  large number of observations. If the iteration limit \code{maxit} in the variable metric algorithm is reached,
#'  one can continue the estimation by iterating more with the function \code{iterate_more}. Alternatively, one may
#'  use the found estimates as starting values for the genetic algorithm and and employ another round of estimation
#'  (see \code{?GAfit} how to set up an initial population with the dot parameters).
#'
#'  \strong{If the estimation algorithm fails to create an initial population for the genetic algorithm},
#'  it usually helps to scale the individual series so that the AR coefficients (of a VAR model) will be
#'  relative small, preferably less than one. Even if one is able to create an initial population, it should
#'  be preferred to scale the series so that most of the AR coefficients will not be very large, as the
#'  estimation algorithm works better with relatively small AR coefficients. If needed, another package can be used
#'  to fit linear VARs to the series to see which scaling of the series results in relatively small AR coefficients.
#'
#'  \strong{weight_constraints:} If you are using weight constraints other than restricting some of the weight parameters to known
#'  constants, make sure the constraints are sensible. Otherwise, the estimation may fail due to the estimation algorithm not being
#'  able to generate reasonable random guesses for the values of the constrained weight parameters.
#'
#'  \strong{Filtering inappropriate estimates:} If \code{filter_estimates == TRUE}, the code will automatically filter
#'  through estimates that it deems "inappropriate". That is, estimates that are not likely solutions of interest.
#'  Specifically, solutions that incorporate a near-singular error term covariance matrix (any eigenvalue less than \eqn{0.002})
#'  or transition weights such that they are close to zero for almost all \eqn{t} for at least one regime.
#'  You are welcome to set \code{filter_estimates=FALSE} and find the solutions of interest yourself by using the
#'  function \code{alt_stvar}.
#'
#'  \strong{weight_constraints}
#'
#' @return Returns an object of class \code{'stvar'} defining the estimated reduced form smooth transition VAR model.
#' @section S3 methods:
#'   The following S3 methods are supported for class \code{'stvar'}: \code{logLik}, \code{residuals}, \code{print}, \code{summary},
#'    \code{predict}, \code{simulate}, and \code{plot}. NONE OF THESE IS IMPLEMENTED YET!
#' @section About the warning "chol(): given matrix is not symmetric":
#'  Sometimes the function produces warnings "chol(): given matrix is not symmetric" because sometimes the underlying C++ code
#'  calculating conditional covariance matrices notices that due to numerical error in the float-point representation, the matrix
#'  is not symmetric. The warning is produces even if symmetricity is forced with arma::symmatu() (and .eval()), so there seems to be some
#'  issue with the LAPACK library that arma::chol uses to check that the matrix is symmetric. It is necessary to take the decomposition
#'  in C++ for computational efficiency in order to calculate its inverse for each t in each evaluation of the log-likelihood function.
#'  These warnings can be ignored.
#' @seealso \code{\link{GAfit}}
#' @references
#'  \itemize{
#'    \item FILL IN REFEENCES
#'  }
#' @examples
#' \donttest{
#' ## These are long running examples that use parallel computing!
#' # Running all the below examples will take approximately FILL IN HOW MANY minutes.
#' # SET SEEDS/ROUNDS
#'
#' # p=3, M=2, d=2 relative_dens Gaussian STVAR
#' fit32 <- fitSTVAR(gdpdef, p=3, M=2, nrounds=20, seeds=1:20)
#' fit32
#' summary(fit32)
#' plot(fit32)
#'
#' # p=1, M=2, d=2 relative_dens Gaussian STVAR with the AR matrices
#' # constrained to be identical in both regimes
#' C_122 <- rbind(diag(1*2^2), diag(1*2^2))
#' fit12c <- fitSTVAR(gdpdef, p=1, M=2, AR_constraints=C_122, nrounds=1, seeds=1, use_parallel=FALSE)
#'
#' # p=1, M=2, d=2 relative_dens Gaussian STVAR with the means constrained
#' # to be identical in both regimes
#' fit12m <- fitSTVAR(gdpdef, p=1, M=2, mean_constraints=list(1:2), nrounds=1, seeds=1, use_parallel=FALSE)
#'
#' # p=1, M=2, d=2 relative_dens Gaussian STVAR with the means and AR matrices constrained
#' # to be identical in both regimes
#' fit12cm <- fitSTVAR(gdpdef, p=1, M=2, AR_constraints=C_122, mean_constraints=list(1:2), nrounds=1, seeds=1, use_parallel=FALSE)
#'
#' # p=2, M=2, d=2, relative_dens Gaussian STVAR; constrain AR-parameters to be the same for all regimes
#' # and constrain the of-diagonal elements of AR-matrices to be zero.
#' mat0 <- matrix(c(1, rep(0, 10), 1, rep(0, 8), 1, rep(0, 10), 1), nrow=2*2^2, byrow=FALSE)
#' C_222 <- rbind(mat0, mat0)
#' fit22c <- fitSTVAR(gdpdef, p=2, M=2, AR_constraints=C_222, nrounds=1, seeds=1, use_parallel=FALSE)
#'
#' # p=2, M=2, d=2, constrain AR-parameters to be the same for all regimes
#' # and constrain the of-diagonal elements of AR-matrices to be zero. Moreover,
#' # constrain the means identical in both regimes.
#' fit22cm <- fitSTVAR(gdpdef, p=2, M=2, AR_constraints=C_222, mean_constraints=list(1:2),
#'   nrounds=1, seeds=1, use_parallel=FALSE)
#'
#' # p=3, M=2, d=3
#' f32 <- fitSTVAR(usamone, p=3, M=2, nrounds=20, seeds=61:80)
#' summary(f32)
#' plot(f32)
#'
#' # mlogit STVAR, p=3, M=2, d=2, the second variable is the switching variable with one lag.
#' fitmlogit32 <- fitSTVAR(gdpdef, p=3, M=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1),
#'  nrounds=20, seeds=1:20)
#'
#' # mlogit STVAR, p=3, M=2, d=2, the second variable is the switching variable with one lag,
#' # constrain AR matrices equal across the regimes.
#' C_322 <- rbind(diag(3*2^2), diag(3*2^2))
#' fitmlogit32c <- fitSTVAR(gdpdef, p=3, M=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1),
#'  AR_constraints=C_322, nrounds=10, seeds=1:10)
#'
#' # relative_dens STVAR, p=1, M=2, d=2, with the weight parameter fixed to the constant r=0.8.
#' fit12w <- fitSTVAR(gdpdef, p=1, M=2, weight_constraints=list(R=0, r=0.8), nrounds=1, seeds=10, use_parallel=FALSE)
#'
#' # mlogit STVAR, p=1, M=2, d=2, second variable as switching variable with one lag such that the constant term is -5.
#' fitmlogit12w <- fitSTVAR(gdpdef, p=1, M=2, weight_function="mlogit", weightfun_pars=list(vars=2, lags=1),
#'   weight_constraints=list(R=matrix(c(0, 1), nrow=2), r=c(-5, 0)), nrounds=20, seeds=1:20, use_parallel=FALSE)
#' }
#' @export

fitSTVAR <- function(data, p, M, weight_function=c("relative_dens", "mlogit"), weightfun_pars=NULL,
                     cond_dist=c("Gaussian", "Student"), parametrization=c("intercept", "mean"),
                     AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL,
                     nrounds=(M + 1)^5, ncores=2, maxit=1000,
                     seeds=NULL, print_res=TRUE, use_parallel=TRUE, filter_estimates=TRUE, ...) {
  # Initial checks etc
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  check_pMd(p=p, M=M)
  if(!all_pos_ints(c(nrounds, ncores, maxit))) stop("Arguments nrounds, ncores, and maxit must be positive integers")
  stopifnot(length(nrounds) == 1)
  if(!is.null(seeds) && length(seeds) != nrounds) stop("The argument 'seeds' should be NULL or a vector of length 'nrounds'")
  data <- check_data(data=data, p=p)
  d <- ncol(data)
  n_obs <- nrow(data)
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars)
  check_constraints(p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=NULL)
  if(!is.null(mean_constraints) && parametrization == "intercept") {
    cat("mean_constraints can be applied for mean-parametrized models only. Switching to parametrization = 'mean'.\n")
    parametrization <- "mean"
  }
  npars <- n_params(p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                    B_constraints=NULL, identification="reduced_form")
  if(npars >= d*nrow(data)) stop("There are at least as many parameters in the model as there are observations in the data,
                                 so you should use smaller p or M.")
  dot_params <- list(...)
  minval <- ifelse(is.null(dot_params$minval), get_minval(data), dot_params$minval)
  red_criteria <- ifelse(rep(is.null(dot_params$red_criteria), 2), c(0.05, 0.01), dot_params$red_criteria)

  if(use_parallel) {
    if(ncores > parallel::detectCores()) {
      ncores <- parallel::detectCores()
      message("ncores was set to be larger than the number of cores detected")
    }
    if(ncores > nrounds) {
      ncores <- nrounds
      message("ncores was set to be larger than the number of estimation rounds")
    }
    cat(paste("Using", ncores, "cores for", nrounds, "estimations rounds..."), "\n")

    ### Optimization with the genetic algorithm ###
    cl <- parallel::makeCluster(ncores)
    on.exit(try(parallel::stopCluster(cl), silent=TRUE)) # Close the cluster on exit, if not already closed.
    parallel::clusterExport(cl, ls(environment(fitSTVAR)), envir=environment(fitSTVAR)) # assign all variables from package:fitSTVAR
    parallel::clusterEvalQ(cl, c(library(pbapply), library(Rcpp), library(RcppArmadillo), library(sstvars)))

    doParallel::registerDoParallel(cl)

    cat("Optimizing with a genetic algorithm...\n")
    GAresults <- pbapply::pblapply(1:nrounds, function(i1) GAfit(data=data, p=p, M=M,
                                                                 weight_function=weight_function,
                                                                 weightfun_pars=weightfun_pars,
                                                                 cond_dist=cond_dist,
                                                                 parametrization=parametrization,
                                                                 AR_constraints=AR_constraints,
                                                                 mean_constraints=mean_constraints,
                                                                 weight_constraints=weight_constraints,
                                                                 seed=seeds[i1], ...), cl=cl)

  } else {
    cat("Optimizing with a genetic algorithm...\n")

    tmpfunGA <- function(i1, ...) {
      cat(i1, "/", nrounds, "\r")
      GAfit(data=data, p=p, M=M,
            weight_function=weight_function,
            weightfun_pars=weightfun_pars,
            cond_dist=cond_dist,
            parametrization=parametrization,
            AR_constraints=AR_constraints,
            mean_constraints=mean_constraints,
            weight_constraints=weight_constraints,
            seed=seeds[i1], ...)
    }
    GAresults <- lapply(1:nrounds, function(i1) tmpfunGA(i1, ...))
  }


  loks <- vapply(1:nrounds, function(i1) loglikelihood(data=data, p=p, M=M,
                                                       params=GAresults[[i1]],
                                                       weight_function=weight_function,
                                                       weightfun_pars=weightfun_pars,
                                                       cond_dist=cond_dist,
                                                       parametrization=parametrization,
                                                       identification="reduced_form",
                                                       AR_constraints=AR_constraints,
                                                       mean_constraints=mean_constraints,
                                                       weight_constraints=weight_constraints,
                                                       B_constraints=NULL,
                                                       to_return="loglik",
                                                       check_params=TRUE,
                                                       minval=minval), numeric(1))

  if(print_res) {
    print_loks <- function() {
      printfun <- function(txt, FUN) cat(paste(txt, round(FUN(loks), 3)), "\n")
      printfun("The lowest loglik: ", min)
      printfun("The largest loglik:", max)
    }
    cat("Results from the genetic algorithm:\n")
    print_loks()
  }

  ### Optimization with the variable metric algorithm ###
  loglik_fn <- function(params) {
    tryCatch(loglikelihood(data=data, p=p, M=M, params=params,
                           weight_function=weight_function, weightfun_pars=weightfun_pars,
                           cond_dist=cond_dist, parametrization=parametrization,
                           identification="reduced_form", AR_constraints=AR_constraints,
                           mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                           B_constraints=NULL, to_return="loglik", check_params=TRUE, minval=minval),
             error=function(e) minval)
  }

  # Gradient of the log-likelihood function using central difference approximation
  h <- 6e-6
  I <- diag(rep(1, npars))
  loglik_grad <- function(params) {
    vapply(1:npars, function(i1) (loglik_fn(params + I[i1,]*h) - loglik_fn(params - I[i1,]*h))/(2*h), numeric(1))
  }

  cat("Optimizing with a variable metric algorithm...\n")
  if(use_parallel) {
    NEWTONresults <- pbapply::pblapply(1:nrounds, function(i1) optim(par=GAresults[[i1]], fn=loglik_fn, gr=loglik_grad, method="BFGS",
                                                                     control=list(fnscale=-1, maxit=maxit)), cl=cl)
    parallel::stopCluster(cl=cl)
  } else {
    tmpfunNE <- function(i1) {
      cat(i1, "/", nrounds, "\r")
      optim(par=GAresults[[i1]], fn=loglik_fn, gr=loglik_grad, method="BFGS",
            control=list(fnscale=-1, maxit=maxit))
    }
    NEWTONresults <- lapply(1:nrounds, function(i1) tmpfunNE(i1))
  }

  loks <- vapply(1:nrounds, function(i1) NEWTONresults[[i1]]$value, numeric(1)) # Log-likelihoods
  converged <- vapply(1:nrounds, function(i1) NEWTONresults[[i1]]$convergence == 0, logical(1)) # Which coverged

  if(print_res) {
    cat("Results from the variable metric algorithm:\n")
    print_loks()
  }

  ### Obtain estimates and filter the inapproriate estimates
  all_estimates <- lapply(NEWTONresults, function(x) x$par)

  if(filter_estimates) {
    cat("Filtering inappropriate estimates...\n")
    ord_by_loks <- order(loks, decreasing=TRUE) # Ordering from largest loglik to smaller

    # Go through estimates, take the estimate that yield the higher likelihood
    # among estimates that are do not include wasted regimes or near-singular
    # error term covariance matrices.
    for(i1 in 1:length(all_estimates)) {
       which_round <- ord_by_loks[i1] # Est round with i1:th largest loglik
       pars <- all_estimates[[which_round]]
       pars_std <- reform_constrained_pars(p=p, M=M, d=d, params=pars,
                                           weight_function=weight_function, weightfun_pars=weightfun_pars,
                                           cond_dist=cond_dist, identification="reduced_form",
                                           AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                           weight_constraints=weight_constraints, B_constraints=NULL) # Pars in standard form for pick pars fns
       Omega_eigens <- get_omega_eigens_par(p=p, M=M, d=d, params=pars_std,
                                            weight_function=weight_function, weightfun_pars=weightfun_pars,
                                            cond_dist=cond_dist, identification="reduced_form",
                                            AR_constraints=NULL, mean_constraints=NULL,
                                            weight_constraints=NULL, B_constraints=NULL)
       Omegas_ok <- !any(Omega_eigens < 0.002)
       if(weight_function == "relative_dens") {
         alphas <- pick_weightpars(p=p, M=M, d=d, params=pars_std, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                   cond_dist=cond_dist)
         weightpars_ok <- !any(alphas < 0.01)
       } else {
         weightpars_ok <- TRUE
       }
       tweights <- loglikelihood(data=data, p=p, M=M, params=pars_std,
                                 weight_function=weight_function, weightfun_pars=weightfun_pars,
                                 cond_dist=cond_dist, parametrization=parametrization,
                                 identification="reduced_form", AR_constraints=NULL,
                                 mean_constraints=NULL, B_constraints=NULL, weight_constraints=NULL,
                                 to_return="tw", check_params=TRUE, minval=matrix(0, nrow=n_obs-p, ncol=M))
       tweights_ok <- !any(vapply(1:M, function(m) sum(tweights[,m] > red_criteria[1]) < red_criteria[2]*n_obs, logical(1)))
       if(Omegas_ok && tweights_ok && weightpars_ok) {
         which_best_fit <- which_round # The estimation round of the appropriate estimate with the largest loglik
         break
       }
       if(i1 == length(all_estimates)) {
         message("No 'appropriate' estimates were found!
                 Check that all the variables are scaled to vary in similar magninutes, also not very small or large magnitudes.
                 Consider running more estimation rounds or study the obtained estimates one-by-one with the function alt_stvar.")
         if(M > 2) {
           message("Consider also using smaller M. Too large M leads to identification problems.")
         }
         which_best_fit <- which(loks == max(loks))[1]
       }
    }
  } else {
    which_best_fit <- which(loks == max(loks))[1]
  }
  params <- all_estimates[[which_best_fit]] # The params to return

  ### Obtain standard errors, calculate IC ###
  # Sort regimes if no constraints are employed
  if(is.null(AR_constraints) && is.null(mean_constraints) && is.null(weight_constraints)) {
    params <- sort_regimes(p=p, M=M, d=d, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
                           cond_dist=cond_dist, identification="reduced_form")
    all_estimates <- lapply(all_estimates, function(pars) sort_regimes(p=p, M=M, d=d, params=pars,
                                                                       weight_function=weight_function,
                                                                       weightfun_pars=weightfun_pars,
                                                                       cond_dist=cond_dist,
                                                                       identification="reduced_form"))
  }

  if(NEWTONresults[[which_best_fit]]$convergence == 1) {
    message("Iteration limit was reached when estimating the best fitting individual!
            Consider further estimation with the function 'iterate_more'")
  }

  transition_weights <- loglikelihood(data=data, p=p, M=M, params=params,
                                      weight_function=weight_function, weightfun_pars=weightfun_pars,
                                      cond_dist=cond_dist, parametrization=parametrization,
                                      identification="reduced_form", AR_constraints=AR_constraints,
                                      mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                      B_constraints=NULL, to_return="tw", check_params=TRUE, minval=minval)
  if(any(vapply(1:M, function(i1) sum(transition_weights[,i1] > red_criteria[1]) < red_criteria[2]*n_obs, logical(1)))) {
    message("At least one of the regimes in the estimated model seems to be wasted in the best fitting individual!")
  }

  ### Wrap up ###
  cat("Calculating approximate standard errors...\n")
  ret <- STVAR(data=data, p=p, M=M, d=d, params=params,
               weight_function=weight_function,
               weightfun_pars=weightfun_pars,
               cond_dist=cond_dist,
               parametrization=parametrization,
               identification="reduced_form",
               AR_constraints=AR_constraints,
               mean_constraints=mean_constraints,
               weight_constraints=weight_constraints,
               B_constraints=NULL,
               calc_std_errors=TRUE)
  ret$all_estimates <- all_estimates
  ret$all_logliks <- loks
  ret$which_converged <- converged
  ret$which_round <- which_best_fit # Which estimation round induced the largest log-likelihood?
  warn_eigens(ret)
  cat("Finished!\n")
  ret
}


