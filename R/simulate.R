#' @title Simulate method for class 'stvar' objects
#'
#' @description \code{simulate.stvar} is a simulate method for class 'stvar' objects.
#'
#' @param object an object of class \code{'stvar'}.
#' @param nsim number of observations to be simulated.
#' @param seed set seed for the random number generator?
#' @param ... currently not in use.
#' @param init_values a size \eqn{(pxd)} matrix specifying the initial values, where d is the number
#'   of time series in the system. The \strong{last} row will be used as initial values for the first lag,
#'   the second last row for second lag etc. If not specified, initial values will be drawn from
#'   the regime specified in \code{init_regimes} (for Gaussian models only).
#' @param init_regime an integer in \eqn{1,...,M} specifying the regime from which
#'   the initial values should be generated from. The initial values will be generated
#'   from the stationary distribution of the specific regime. Due to the (lack of)
#'   knowledge of the stationary distribution, only model with Gaussian conditional distribution
#'   are supported. For Student's t models, specify \code{init_values}.
#' @param ntimes how many sets of simulations should be performed?
#' @param drop if \code{TRUE} (default) then the components of the returned list are coerced to lower dimension if \code{ntimes==1}, i.e.,
#'   \code{$sample} and \code{$transition_weights} will be matrices, and \code{$component} will be vector.
#' @details The argument \code{ntimes} is intended for forecasting, which is used by the predict method (see \code{?predict.stvar}).
#' @return If \code{drop==TRUE} and \code{ntimes==1} (default): \code{$sample}, \code{$component}, and \code{$transition_weights} are matrices.
#'   Otherwise, returns a list with...
#'   \describe{
#'     \item{\code{$sample}}{a size (\code{nsim}\eqn{ x d x }\code{ntimes}) array containing the samples: the dimension \code{[t, , ]} is
#'      the time index, the dimension \code{[, d, ]} indicates the marginal time series, and the dimension \code{[, , i]} indicates
#'      the i:th set of simulations.}
#'     \item{\code{$transition_weights}}{a size (\code{nsim}\eqn{ x M x }\code{ntimes}) array containing the transition weights
#'      corresponding to the sample: the dimension \code{[t, , ]} is the time index, the dimension \code{[, m, ]} indicates the regime,
#'      and the dimension \code{[, , i]} indicates the i:th set of simulations.}
#'   }
#' @seealso \code{\link{fitSTVAR}}, \code{\link{STVAR}}, \code{\link{predict.stvar}}
#' @inherit loglikelihood references
#' @examples
#'  # FILL IN EXAMPLES
#' @export

simulate.stvar <- function(object, nsim=1, seed=NULL, ..., init_values=NULL, init_regime,
                            ntimes=1, drop=TRUE) {
  # Checks etc
  if(!is.null(seed)) set.seed(seed)
  epsilon <- round(log(.Machine$double.xmin) + 10)
  stvar <- object
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  weight_function <- stvar$model$weight_function
  cond_dist <- stvar$model$cond_dist
  identification <- stvar$model$identification
  if(identification != "reduced_form") stop("Structural models are not yet implemented to simulate.stvar")
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  B_constraints <- stvar$model$B_constraints
  if(!is.null(AR_constraints) || !is.null(mean_constraints) || !is.null(B_constraints)) {
    stop("Constained models are not yet implemented to simulate.stvar")
  }
  if(cond_dist != "Gaussian") stop("Other that Gaussian models are not yet implemented to simulate.stvar")
  if(is.null(init_values) & missing(init_regime)) {
    stop("Either init_values or init_regime needs to be specified")
  }
  if(!missing(init_regime)) {
    if(cond_dist != "Gaussian") {
      stop("init_regime is currently implemented for Gaussian models only.
           Please specify init_values instead.")
    }
    stopifnot(init_regime %in% 1:M)
  }
  if(!all_pos_ints(c(nsim, ntimes))) stop("Arguments nsim and ntimes must be positive integers")
  if(!is.null(init_values)) {
    if(!is.matrix(init_values)) stop("init_values must be a numeric matrix")
    if(anyNA(init_values)) stop("init_values contains NA values")
    if(ncol(init_values) != d | nrow(init_values) < p) stop("init_values must contain d columns and at least p rows")
  }

  # Collect parameter values
  params <- stvar$params
  ## REFORM CONSTRAINED PARS
  if(stvar$model$parametrization == "mean") {
    params <- change_parametrization(p=p, M=M, d=d, param=params, AR_onstraints=NULL,
                                     mean_constraints=NULL, change_to="intercept")
  }
  all_mu <- get_regime_means(stvar)
  all_phi0 <- pick_phi0(M=M, d=d, params=params)
  all_A <- pick_allA(p=p, M=M, d=d, params=params)
  all_Omegas <- pick_Omegas(p=p, M=M, d=d, params=params) # Note that structural models not implemented here
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  weightpars <- pick_weightpars(p=p, M=M, d=d, params=params, weight_function=weight_function, cond_dist=cond_dist)
  # pick_distpars

  # Calculate statistics that remain constant through the iterations
  if(weight_function == "relative_dens") {
     Sigmas <- get_Sigmas(p=p, M=M, d=d, all_A=all_A, all_boldA=all_boldA, all_Omegas=all_Omegas)
     inv_Sigmas <- array(NA, dim=c(d*p, d*p, M)) # Store inverses of the (dpxdp) covariance matrices
     det_Sigmas <- numeric(M) # Store determinants of the (dpxdp) covariance matrices
     chol_Sigmas <- array(dim=c(d*p, d*p, M)) # Cholesky decompositions of the  (dpxdp) covariance matrices
     for(m in 1:M) {
       chol_Sigmas[, , m] <- chol(Sigmas[, , m]) # Upper triangle
       inv_Sigmas[, , m] <- chol2inv(chol_Sigmas[, , m]) # Faster inverse
       det_Sigmas[m] <- prod(diag(chol_Sigmas[, , m]))^2 # Faster determinant
     }
  } else {
    stop("Other than relative_dens weight functions are not yet implemented to simulate.stvar")
  }

  # Set/generate initial values
  if(is.null(init_values)) {
    # Generate the initial values from the stationary distribution of init_regime; Gaussian dist
    mu <- rep(all_mu[, init_regime], p)
    L <- t(chol_Sigmas[, , init_regime]) # Lower triangle
    init_values <- matrix(mu + L%*%rnorm(d*p), nrow=p, ncol=d, byrow=TRUE) # i:th row for the i:th length d random vector
  } else {
    # Initial values given as argument
    init_values <- init_values[(nrow(init_values) - p + 1):nrow(init_values), , drop=FALSE]
  }


  # Container for the simulated values and initial values. First row row initial values vector, and t:th row for (y_{i-1},...,y_{i-p})
  Y <- matrix(nrow=nsim + 1, ncol=d*p)
  Y[1,] <- reform_data(init_values, p=p)

  # Initialize data structures
  sample <- array(dim=c(nsim, d, ntimes))
  component <- matrix(nrow=nsim, ncol=ntimes)
  transition_weights <- array(dim=c(nsim, M, ntimes))

  # Some functions to be used
  if(weight_function == "relative_dens") {
    # Get log multivariate normal densities for calculating the transition weights
    get_logmvdvalues <- function(Y, i1) {
      vapply(1:M,
             function(m) -0.5*d*p*log(2*pi) - 0.5*log(det_Sigmas[m]) - 0.5*(crossprod(Y[i1,] - rep(all_mu[, m], p),
                                                                                      inv_Sigmas[, , m])%*%(Y[i1,] - rep(all_mu[, m], p))),
             numeric(1))
    } # Returns M x 1 vector; transformed into a matrix in get_alpha_mt
  }

  # Run through the time periods and repetitions
  for(j1 in 1:seq_len(ntimes)) {
    for(i1 in seq_len(nsim)) {
      # Calculate transition weights
      if(weight_function == "relative_dens") {
        log_mvdvalues <- get_logmvdvalues(Y=Y, i1=i1)
        alpha_mt <- get_alpha_mt(M=M, log_mvdvalues=log_mvdvalues, epsilon=epsilon)
      } else {
        stop("Only relative_dens weight function is implemented to simulate.stvar")
      }
      transition_weights[i1, , j1] <- alpha_mt

      # Calculate conditional mean
      mu_yt <- get_mu_yt_Cpp(obs=Y[,i1], all_phi0=all_phi0, all_A=all_A, alpha_mt=alpha_mt)

      # Calculate conditional covariance matrix
      Omega_yt <- matrix(rowSums(vapply(1:M, function(m) alpha_mt[, m]*as.vector(all_Omegas[, , m]),
                                        numeric(d*d))), nrow=d, ncol=d)
      # Calculate B_t
      if(identification == "reduced_form") {
        B_t <- get_symmetric_sqrt(Omega_yt)
      } else {
        stop("Structural models are not yet implemented to simulate.stvar")
      }

      # Draw the structural error
      if(cond_dist == "Gaussian") {
        e_t <- rnorm(d)
      } else {
        stop("Other than Gaussian models are not yet implemented to simulat.stvar")
      }

      # Calculate the current observation
      sample[i1, , j1] <- mu_yt + B_t%*%e_t

      # Update storage matrix Y (overwrites when ntimes > 1)
      if(p == 1) {
        Y[i1 + 1,] <- sample[i1, , j1]
      } else {
        Y[i1 + 1,] <- c(sample[i1, , j1], Y[i1, 1:(d*p - d)])
      }
    }
  }

  if(ntimes == 1 & drop) {
    sample <- matrix(sample, nrow=nsim, ncol=d)
    transition_weights <- matrix(transition_weights, nrow=nsim, ncol=M)
  }

  list(sample=sample,
       transition_weights=transition_weights)
}
