#' @title Genetic algorithm for preliminary estimation of a GMVAR, StMVAR, or G-StMVAR model
#'
#' @description \code{GAfit} estimates the specified GMVAR, StMVAR, or G-StMVAR model using a genetic algorithm.
#'   It's designed to find starting values for gradient based methods.
#'
#' @inheritParams loglikelihood
#' @inheritParams random_covmat
#' @param ngen a positive integer specifying the number of generations to be ran through in
#'   the genetic algorithm.
#' @param popsize a positive even integer specifying the population size in the genetic algorithm.
#'   Default is \code{10*n_params}.
#' @param smart_mu a positive integer specifying the generation after which the random mutations
#'   in the genetic algorithm are "smart". This means that mutating individuals will mostly mutate fairly
#'   close (or partially close) to the best fitting individual (which has the least regimes with time varying
#'   mixing weights practically at zero) so far.
#' @param initpop a list of parameter vectors from which the initial population of the genetic algorithm
#'   will be generated from. The parameter vectors should have the form
#'   \eqn{\theta = (\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M,\sigma,\alpha,\nu)},
#'   where
#'   \itemize{
#'     \item{\eqn{\phi_{m,0} = } the \eqn{(d \times 1)} intercept (or mean) vector of the \eqn{m}th regime.}
#'     \item{\eqn{\varphi_m = (vec(A_{m,1}),...,vec(A_{m,p}))} \eqn{(pd^2 \times 1)}.}
#'     \item{\eqn{\sigma = (vech(\Omega_1),...,vech(\Omega_M)} \eqn{(Md(d + 1)/2 \times 1)}.}
#'     \item{\eqn{\alpha} contains the transition weights parameters}
#'     \item{\eqn{\nu > 2} is the degrees of freedom parameter that is included only if \code{cond_dist="Student"}.}
#'   }
#'   \describe{
#'     \item{For models with \code{weight_function="relative_dens"}:}{\eqn{\alpha = (\alpha_1,...,\alpha_{M-1})}
#'           \eqn{(M - 1 \times 1)}, where \eqn{\alpha_m} \eqn{(1\times 1), m=1,...,M-1} are the transition weight parameters.}
#'     \item{For models with \code{weight_function="logit"}:}{\eqn{\alpha = (\gamma_1,...,\gamma_M)} \eqn{((M-1)k\times 1)},
#'           where \eqn{\gamma_m} \eqn{(k\times 1), m=1,...,M-1} contains the logit-regression coefficients of the \eqn{m}th
#'            regime.}
#'   }
#'   Above, \eqn{\phi_{m,0}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}th coefficient matrix of the \eqn{m}th
#'   mixture component, and \eqn{\Omega_{m}} denotes the error term covariance matrix of the \eqn{m}:th mixture component.
#'   If \code{parametrization=="mean"}, just replace each \eqn{\phi_{m,0}} with regimewise mean \eqn{\mu_{m}}.
#'   \eqn{vec()} is vectorization operator that stacks columns of a given matrix into a vector. \eqn{vech()} stacks columns
#'   of a given matrix from the principal diagonal downwards (including elements on the diagonal) into a vector.
#' @param conditional a logical argument specifying whether the conditional or exact log-likelihood function
#' @param mu_scale a size \eqn{(dx1)} vector defining \strong{means} of the normal distributions from which each
#'   mean parameter \eqn{\mu_{m}} is drawn from in random mutations. Default is \code{colMeans(data)}. Note that
#'   mean-parametrization is always used for optimization in \code{GAfit} - even when \code{parametrization=="intercept"}.
#'   However, input (in \code{initpop}) and output (return value) parameter vectors can be intercept-parametrized.
#' @param mu_scale2 a size \eqn{(dx1)} strictly positive vector defining \strong{standard deviations} of the normal
#'   distributions from which each mean parameter \eqn{\mu_{m}} is drawn from in random mutations.
#'   Default is \code{2*sd(data[,i]), i=1,..,d}.
#' @param omega_scale a size \eqn{(dx1)} strictly positive vector specifying the scale and variability of the
#'   random covariance matrices in random mutations. The covariance matrices are drawn from (scaled) Wishart
#'   distribution. Expected values of the random covariance matrices are \code{diag(omega_scale)}. Standard
#'   deviations of the diagonal elements are \code{sqrt(2/d)*omega_scale[i]}
#'   and for non-diagonal elements they are \code{sqrt(1/d*omega_scale[i]*omega_scale[j])}.
#'   Note that for \code{d>4} this scale may need to be chosen carefully. Default in \code{GAfit} is
#'   \code{var(stats::ar(data[,i], order.max=10)$resid, na.rm=TRUE), i=1,...,d}. This argument is ignored if
#'   structural model is considered.
#' @param ar_scale a positive real number adjusting how large AR parameter values are typically proposed in construction
#'   of the initial population: larger value implies larger coefficients (in absolute value). After construction of the
#'   initial population, a new scale is drawn from \code{(0, 0.)} uniform distribution in each iteration.
#' @param upper_ar_scale the upper bound for \code{ar_scale} parameter (see above) in the random mutations. Setting
#'  this too high might lead to failure in proposing new parameters that are well enough inside the parameter space,
#'  and especially with large \code{p} one might want to try smaller upper bound (e.g., 0.5).
#' @param ar_scale2 a positive real number adjusting how large AR parameter values are typically proposed in some
#'   random mutations (if AR constraints are employed, in all random mutations): larger value implies \strong{smaller}
#'   coefficients (in absolute value). \strong{Values larger than 1 can be used if the AR coefficients are expected to
#'   be very small. If set smaller than 1, be careful as it might lead to failure in the creation of parameter candidates
#'   that satisfy the stability condition.}
#' @param regime_force_scale a non-negative real number specifying how much should natural selection favor individuals
#'   with less regimes that have almost all mixing weights (practically) at zero. Set to zero for no favoring or large
#'   number for heavy favoring. Without any favoring the genetic algorithm gets more often stuck in an area of the
#'   parameter space where some regimes are wasted, but with too much favouring the best genes might never mix into
#'   the population and the algorithm might converge poorly. Default is \code{1} and it gives \eqn{2x} larger surviving
#'   probability weights for individuals with no wasted regimes compared to individuals with one wasted regime.
#'   Number \code{2} would give \eqn{3x} larger probability weights etc.
#' @param red_criteria a length 2 numeric vector specifying the criteria that is used to determine whether a regime is
#'   redundant (or "wasted") or not.
#'   Any regime \code{m} which satisfies \code{sum(transitionWeights[,m] > red_criteria[1]) < red_criteria[2]*n_obs} will
#'   be considered "redundant". One should be careful when adjusting this argument (set \code{c(0, 0)} to fully disable
#'   the 'redundant regime' features from the algorithm).
#' @param pre_smart_mu_prob A number in \eqn{[0,1]} giving a probability of a "smart mutation" occuring randomly in each
#'   iteration before the iteration given by the argument \code{smart_mu}.
#' @param to_return should the genetic algorithm return the best fitting individual which has "positive enough" mixing
#'   weights for as many regimes as possible (\code{"alt_ind"}) or the individual which has the highest log-likelihood
#'   in general (\code{"best_ind"}) but might have more wasted regimes?
#' @param minval a real number defining the minimum value of the log-likelihood function that will be considered.
#'   Values smaller than this will be treated as they were \code{minval} and the corresponding individuals will
#'   never survive. The default is \code{-(10^(ceiling(log10(n_obs)) + d) - 1)}.
#' @param seed a single value, interpreted as an integer, or NULL, that sets seed for the random number generator in
#'   the beginning of the function call. If calling \code{GAfit} from \code{fitGSMVAR}, use the argument \code{seeds}
#'   instead of passing the argument \code{seed}.
#' @details
#'  Only reduced form models are supported!
#'
#'  The core of the genetic algorithm is mostly based on the description by \emph{Dorsey and Mayer (1995)}.
#'  It utilizes a slightly modified version of the individually adaptive crossover and mutation rates described
#'  by \emph{Patnaik and Srinivas (1994)} and employs (50\%) fitness inheritance discussed by
#'  \emph{Smith, Dike and Stegmann (1995)}.
#'
#'  By "redundant" or "wasted" regimes we mean regimes that have the time varying mixing weights practically at
#'  zero for almost all t. A model including redundant regimes would have about the same log-likelihood value without
#'  the redundant regimes and there is no purpose to have redundant regimes in a model.
#'
#'  Structural models are not supported here, as they are best estimated based on reduced form parameter estimates
#'  using the function FILL IN THE FUNCTION WHEN IT HAS BEEN CREATED.
#' @return Returns the estimated parameter vector which has the form described in \code{initpop}.
#' @references
#'  \itemize{
#'    \item Ansley C.F., Kohn R. 1986. A note on reparameterizing a vector autoregressive
#'          moving average model to enforce stationarity. \emph{Journal of statistical computation
#'          and simulation}, \strong{24}:2,  99-106.
#'    \item Dorsey R. E. and Mayer W. J. 1995. Genetic algorithms for estimation problems with multiple optima,
#'          nondifferentiability, and other irregular features. \emph{Journal of Business & Economic Statistics},
#'         \strong{13}, 53-66.
#'    \item Patnaik L.M. and Srinivas M. 1994. Adaptive Probabilities of Crossover and Mutation in Genetic Algorithms.
#'          \emph{Transactions on Systems, Man and Cybernetics} \strong{24}, 656-667.
#'    \item Smith R.E., Dike B.A., Stegmann S.A. 1995. Fitness inheritance in genetic algorithms.
#'          \emph{Proceedings of the 1995 ACM Symposium on Applied Computing}, 345-350.
#'  }

GAfit <- function(data, p, M, weight_function=c("relative_dens", "logit"), cond_dist=c("Gaussian", "Student"),
                  parametrization=c("intercept", "mean"), AR_constraints=NULL, mean_constraints=NULL,
                  ngen=200, popsize, smart_mu=min(100, ceiling(0.5*ngen)), initpop=NULL, mu_scale, mu_scale2, omega_scale,
                  ar_scale=0.2, upper_ar_scale=1, ar_scale2=1, regime_force_scale=1, red_criteria=c(0.05, 0.01),
                  pre_smart_mu_prob=0, to_return=c("alt_ind", "best_ind"), minval, seed=NULL) {

  # Required values and preliminary checks
  set.seed(seed)
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  if(!is.null(AR_constraints) || !is.null(mean_constraints) || !is.null(B_constraints)) {
    stop("Constrained models are not currently supported")
  }
  check_pMd(p=p, M=M)
  data <- check_data(data=data, p=p)
  d <- ncol(data)
  n_obs <- nrow(data)
  n_pars <- n_params(p=p, M=M, d=d, weight_function=weight_function, cond_dist=cond_dist,
                     identification=identification, AR_constraints=AR_constraints,
                     mean_constraints=mean_constraints, B_constraints=B_constraints)
  # FILL IN ARGUMENT CHECKS FOR CONSTRAINTS

  # Defaults and checks
  if(!all_pos_ints(c(ngen, smart_mu))) stop("Arguments ngen and smart_mu have to be positive integers")
  if(missing(popsize)) {
    popsize <- 50*ceiling(sqrt(npars))
  } else if(popsize < 2 | popsize %% 2 != 0) {
    stop("The population size popsize must be even positive integer")
  }
  if(missing(minval)) {
    minval <- get_minval(data)
  } else if(!is.numeric(minval)) {
    stop("Argument minval must be numeric")
  }
  if(missing(mu_scale)) {
    mu_scale <- colMeans(data)
  } else if(length(mu_scale) != d) {
    stop("Argument mu_scale must be numeric vector with length d")
  }
  if(missing(mu_scale2)) {
    mu_scale2 <- vapply(1:d, function(i1) 2*sd(data[,i1]), numeric(1))
  } else if(length(mu_scale2) != d | any(mu_scale2 <= 0)) {
    stop("Argument mu_scale2 must be strictly positive vector with length d")
  }
  if(missing(omega_scale)) {
    omega_scale <- vapply(1:d, function(i1) var(stats::ar(data[,i1], order.max=10)$resid, na.rm=TRUE), numeric(1))
  } else if(!(length(omega_scale) == d & all(omega_scale > 0))) {
    stop("omega_scale must be numeric vector with length d and positive elements")
  }
  stopifnot(pre_smart_mu_prob >= 0 && pre_smart_mu_prob <= 1)
  if(length(ar_scale) != 1 | ar_scale <= 0) {
    stop("ar_scale must be positive and have length one")
  }
  if(length(ar_scale2) != 1 | ar_scale2 <= 0) {
    stop("ar_scale2 must be positive and have length one")
  } else if(ar_scale2 > 1.5) {
    warning("Large ar_scale2 might lead to failure of the estimation process")
  }
  if(length(regime_force_scale) != 1 | regime_force_scale < 0) {
    stop("regime_force_scale should be non-negative real number")
  }

  # The initial population
  if(is.null(initpop)) {
    # CONSTRUCT THE INITIAL POPULATION
    n_attempts <- 20
    G <- numeric(0)
    for(i1 in 1:n_attempts) {
      inds <- replicate(popsize, random_ind(p=p, M=M, d=d, weight_function=weight_function,
                                            cond_dist=cond_dist, AR_constraints=AR_constraints,
                                            mean_constraints=mean_constraints,
                                            force_stability=is.null(AR_constraints),
                                            mu_scale=mu_scale, mu_scale2=mu_scale2,
                                            omega_scale=omega_scale, ar_scale=ar_scale,
                                            ar_scale2=ar_scale2))
      ind_loks <- vapply(1:popsize, function(i2) loglikelihood(data=data, p=p, M=M, params=ind[,i2],
                                                               weight_function=weight_function, cond_dist=cond_dist,
                                                               parametrization="mean", identification="reduced_form",
                                                               AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                                               B_constraints=NULL, to_return="loglik", check_params=TRUE,
                                                               minval=minval), numeric(1))
      G <- cbind(G, inds[, ind_loks > minval]) # Take good enough individuals
      if(ncol(G) >= popsize) {
        G <- G[, 1:popsize]
        break
      } else if(i1 == nattempts) {
        if(length(G) == 0) {
          stop("Failed to create initial population with good enough individuals. Scaling the individual series
               so that the AR coefficients (of a VAR model) will not be very large (preferably less than one)
               may solve the problem. If needed, another package may be used to fit linear VARs so see which
               scalings produce relatively small AR coefficient estimates.")
        } else {
          G <- G[, sample.int(ncol(G), size=popsize, replace=TRUE)]
        }
      }
    }
  } else {  # Initial population set by the user
    stopifnot(is.list(initpop))
    for(i1 in 1:length(initpop)) {
      ind <- initpop[[i1]]
      tryCatch(check_params(p=p, M=M, d=d, params=ind,  weight_function=weight_function, cond_dist=cond_dist,
                            parametrization=parametrization, identification="reduced_form", AR_constraints=AR_constraints,
                            mean_constraints=mean_constraints, B_constraints=NULL),
               error=function(e) stop(paste("Problem with individual", i1, "in the initial population: "), e))
      if(parametrization == "intercept") {
        ind <- change_parametrization(p=p, M=M, d=d, params=ind, AR_constraints=AR_constraints,
                                      mean_constraints=mean_constaints, change_to="mean")
      }
      if(is.null(AR_constraints) && is.null(mean_constraints)) {
        initpop[[i1]] <- sort_components(p=p, M=M, d=d, params=ind, weight_function=weight_function,
                                         cond_dist=cond_dist, identification="reduced_form")
      } else {
        initpop[[i1]] <- ind
      }
    }
    G <- replicate(popsize, initpop[[sample.int(length(initpop), size=1)]])
  }

  # Calculates the number of redundant regimes
  n_redundants <- function(M, tw) {
    sum(vapply(1:M, function(m) sum(tw[,m] > red_criteria[1]) < red_criteria[2]*n_obs, logical(1)))
  }

  # Initial setup
  generations <- array(NA, dim=c(npars, popsize, ngen))
  logliks <- matrix(minval, nrow=ngen, ncol=popsize)
  redundants <- matrix(M, nrow=ngen, ncol=popsize) # Store the number of redundant regimes of each individual
  which_redundant_alt <- 1:M

  fill_lok_and_red <- function(i1, i2, lok_and_tw) {
    if(!is.list(lok_and_tw)) {
      logliks[i1, i2] <<- minval
      redundants[i1, i2] <<- M
    } else {
      logliks[i1, i2] <<- lok_and_tw$loglik
      redundants[i1, i2] <<- n_redundants(M=M, tw=lok_and_tw$tw) # Number of redundant regimes
    }
  }

  # Run through the generations
  for(i1 in 1:ngen) {

  }
}
