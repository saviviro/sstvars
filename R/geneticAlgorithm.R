#' @import stats
#'
#' @title Genetic algorithm for preliminary estimation of a STVAR models
#'
#' @description \code{GAfit} estimates the specified STVAR model using a genetic algorithm.
#'   It is designed to find starting values for gradient based methods and NOT to obtain
#'   final estimates constituting a local maximum.
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
#'     \item{\describe{
#'       \item{if \code{cond_dist="Gaussian"} or \code{"Student"}:}{\eqn{\sigma = (vech(\Omega_1),...,vech(\Omega_M))}
#'         \eqn{(Md(d + 1)/2 \times 1)}.}
#'       \item{if \code{cond_dist="ind_Student"}:}{\eqn{\sigma = (vec(B_1),...,vec(B_M)} \eqn{(Md^2 \times 1)}.}
#'       }
#'     }
#'     \item{\eqn{\alpha} contains the transition weights parameters (see below)}
#'     \item{\describe{
#'       \item{if \code{cond_dist = "Gaussian")}:}{Omit \eqn{\nu} from the parameter vector.}
#'       \item{if \code{cond_dist="Student"}:}{\eqn{\nu > 2} is the single degrees of freedom parameter.}
#'       \item{if \code{cond_dist="ind_Student"}:}{\eqn{\nu = (\nu_1,...,\nu_M)} \eqn{(M \times 1)}, \eqn{nu_m > 2}.}
#'       }
#'     }
#'   }
#'   \describe{
#'     \item{\code{weight_function="relative_dens"}:}{\eqn{\alpha = (\alpha_1,...,\alpha_{M-1})}
#'           \eqn{(M - 1 \times 1)}, where \eqn{\alpha_m} \eqn{(1\times 1), m=1,...,M-1} are the transition weight parameters.}
#'    \item{\code{weight_function="logistic"}:}{\eqn{\alpha = (c,\gamma)}
#'           \eqn{(2 \times 1)}, where \eqn{c\in\mathbb{R}} is the location parameter and \eqn{\gamma >0} is the scale parameter.}
#'     \item{\code{weight_function="mlogit"}:}{\eqn{\alpha = (\gamma_1,...,\gamma_M)} \eqn{((M-1)k\times 1)},
#'           where \eqn{\gamma_m} \eqn{(k\times 1)}, \eqn{m=1,...,M-1} contains the multinomial logit-regression coefficients
#'           of the \eqn{m}th regime. Specifically, for switching variables with indices in \eqn{I\subset\lbrace 1,...,d\rbrace}, and with
#'          \eqn{\tilde{p}\in\lbrace 1,...,p\rbrace} lags included, \eqn{\gamma_m} contains the coefficients for the vector
#'          \eqn{z_{t-1} = (1,\tilde{z}_{\min\lbrace I\rbrace},...,\tilde{z}_{\max\lbrace I\rbrace})}, where
#'          \eqn{\tilde{z}_{i} =(y_{it-1},...,y_{it-\tilde{p}})}, \eqn{i\in I}. So \eqn{k=1+|I|\tilde{p}}
#'          where \eqn{|I|} denotes the number of elements in \eqn{I}.}
#'     \item{\code{weight_function="exponential"}:}{\eqn{\alpha = (c,\gamma)}
#'           \eqn{(2 \times 1)}, where \eqn{c\in\mathbb{R}} is the location parameter and \eqn{\gamma >0} is the scale parameter.}
#'     \item{\code{weight_function="threshold"}:}{\eqn{\alpha = (r_1,...,r_{M-1})}
#'           \eqn{(M-1 \times 1)}, where \eqn{r_1,...,r_{M-1}} are the threshold values.}
#'     \item{\code{weight_function="exogenous"}:}{Omit \eqn{\alpha} from the parameter vector.}
#'     \item{AR_constraints:}{Replace \eqn{\varphi_1,...,\varphi_M} with \eqn{\psi} as described in the argument \code{AR_constraints}.}
#'     \item{mean_constraints:}{Replace \eqn{\phi_{1,0},...,\phi_{M,0}} with \eqn{(\mu_{1},...,\mu_{g})} where
#'           \eqn{\mu_i, \ (d\times 1)} is the mean parameter for group \eqn{i} and \eqn{g} is the number of groups.}
#'     \item{weight_constraints:}{If linear constraints are imposed, replace \eqn{\alpha} with \eqn{\xi} as described in the
#'      argument \code{weigh_constraints}. If weight functions parameters are imposed to be fixed values, simply drop \eqn{\alpha}
#'      from the parameter vector.}
#'   }
#'   Above, \eqn{\phi_{m,0}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}th coefficient matrix of the \eqn{m}th
#'   regime, \eqn{\Omega_{m}} denotes the positive definite error term covariance matrix of the \eqn{m}th regime, and \eqn{B_m}
#'   is the invertible \eqn{(d\times d)} impact matrix of the \eqn{m}th regime. \eqn{\nu_m} is the degrees of freedom parameter
#'   of the \eqn{m}th regime.
#'   If \code{parametrization=="mean"}, just replace each \eqn{\phi_{m,0}} with regimewise mean \eqn{\mu_{m}}.
#'   \eqn{vec()} is vectorization operator that stacks columns of a given matrix into a vector. \eqn{vech()} stacks columns
#'   of a given matrix from the principal diagonal downwards (including elements on the diagonal) into a vector. \eqn{Bvec()}
#'   is a vectorization operator that stacks the columns of a given impact matrix \eqn{B_m} into a vector so that the elements
#'   that are constrained to zero by the argument \code{B_constraints} are excluded.
#' @param mu_scale a size \eqn{(dx1)} vector defining \strong{means} of the normal distributions from which each
#'   mean parameter \eqn{\mu_{m}} is drawn from in random mutations. Default is \code{colMeans(data)}. Note that
#'   mean-parametrization is always used for optimization in \code{GAfit} - even when \code{parametrization=="intercept"}.
#'   However, input (in \code{initpop}) and output (return value) parameter vectors can be intercept-parametrized.
#' @param mu_scale2 a size \eqn{(dx1)} strictly positive vector defining \strong{standard deviations} of the normal
#'   distributions from which each mean parameter \eqn{\mu_{m}} is drawn from in random mutations.
#'   Default is \code{vapply(1:d, function(i1) sd(data[,i1]), numeric(1))}.
#' @param omega_scale a size \eqn{(dx1)} strictly positive vector specifying the scale and variability of the
#'   random covariance matrices in random mutations. The covariance matrices are drawn from (scaled) Wishart
#'   distribution. Expected values of the random covariance matrices are \code{diag(omega_scale)}. Standard
#'   deviations of the diagonal elements are \code{sqrt(2/d)*omega_scale[i]}
#'   and for non-diagonal elements they are \code{sqrt(1/d*omega_scale[i]*omega_scale[j])}.
#'   Note that for \code{d>4} this scale may need to be chosen carefully. Default in \code{GAfit} is
#'   \code{var(stats::ar(data[,i], order.max=10)$resid, na.rm=TRUE), i=1,...,d}. This argument is ignored if
#'   \code{cond_dist == "ind_Student"}.
#' @param B_scale a size \eqn{(d \times 1)} strictly positive vector specifying the mean and variability of the
#'   random impact matrices in random mutations. In Regime 1, the mean of the error term covariance matrix
#'   implied by the random impact matrix will be \code{0.95*diag(B_scale)} and in the rest of the regimes \code{diag(B_scale)},
#'   whereas the variability increases with \code{B_scale}.
#'   Default in \code{GAfit} is \code{var(stats::ar(data[,i], order.max=10)$resid, na.rm=TRUE), i=1,...,d}.
#'   This argument is ignored if \code{cond_dist != "ind_Student"}.
#' @param ar_scale a positive real number between zero and one adjusting how large AR parameter values are typically
#'   proposed in construction of the initial population: larger value implies larger coefficients (in absolute value).
#'   After construction of the initial population, a new scale is drawn from \code{(0, upper_ar_scale)} uniform
#'   distribution in each iteration.
#' @param weight_scale For...
#'   \describe{
#'     \item{\code{weight_function \%in\% c("relative_dens", "exogenous")}:}{not used.}
#'     \item{\code{weight_function \%in\% c("logistic", "exponential")}:}{length three vector with the mean (in the first element)
#'        and standard deviation (in the second element) of the normal distribution the location parameter is drawn from
#'        in random mutations. The third element is the standard deviation of the normal distribution from whose absolute value
#'        the location parameter is drawn from.}
#'     \item{\code{weight_function == "mlogit"}:}{length two vector with the mean (in the first element)
#'        and standard deviation (in the second element) of the normal distribution the coefficients of the logit sub model's
#'        constant terms are drawn from in random mutations. The third element is the standard deviation of the normal distribution
#'        from which the non-constant regressors' coefficients are drawn from.}
#'     \item{\code{weight_function == "threshold"}:}{a lenght two vector with the lower bound, in the first element
#'        and the upper bound, in the second element, of the uniform distribution threshold parameters are drawn from
#'        in random mutations.}
#'   }
#' @param upper_ar_scale the upper bound for \code{ar_scale} parameter (see above) in the random mutations. Setting
#'  this too high might lead to failure in proposing new parameters that are well enough inside the parameter space,
#'  and especially with large \code{p} one might want to try smaller upper bound (e.g., 0.5). With large \code{p} or
#'  \code{d}, \code{upper_ar_scale} is restricted from above, see the details section.
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
#'   the beginning of the function call. If calling \code{GAfit} from \code{fitSTVAR}, use the argument \code{seeds}
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
#'  Some of the AR coefficients are drawn with the algorithm by Ansley and Kohn (1986). However,
#'  when using large \code{ar_scale} with large \code{p} or \code{d}, numerical inaccuracies caused
#'  by the imprecision of the float-point presentation may result in errors or nonstationary AR-matrices.
#'  Using smaller \code{ar_scale} facilitates the usage of larger \code{p} or \code{d}. Therefore, we bound
#'  \code{upper_ar_scale} from above by \eqn{1-pd/150} when \code{p*d>40} and by \eqn{1} otherwise.
#'
#'  Structural models are not supported here, as they are best estimated based on reduced form parameter estimates
#'  using the function \code{fitSSTVAR}.
#' @return Returns the estimated parameter vector which has the form described in \code{initpop},
#'  \strong{with the exception} that for models with \code{cond_dist == "ind_Student"} or
#'  \code{identification="non-Gaussianity"}, the parameter vector is parametrized with \eqn{B_1,B_2^*,...,B_M^*}
#'  instead of  \eqn{B_1,B_2,...,B_M}, where \eqn{B_m^* = B_m - B_1}. Use the function \code{change_parametrization}
#'  to change back to the original parametrization if desired.
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

GAfit <- function(data, p, M, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                  weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student"), parametrization=c("intercept", "mean"),
                  AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL,
                  ngen=200, popsize, smart_mu=min(100, ceiling(0.5*ngen)), initpop=NULL,
                  mu_scale, mu_scale2, omega_scale, B_scale, weight_scale,
                  ar_scale=0.2, upper_ar_scale=1, ar_scale2=1, regime_force_scale=1, red_criteria=c(0.05, 0.01),
                  pre_smart_mu_prob=0, to_return=c("alt_ind", "best_ind"), minval, seed=NULL) {

  # Required values and preliminary checks
  set.seed(seed)
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  to_return <- match.arg(to_return)
  check_pMd(p=p, M=M, weight_function=weight_function, identification="reduced_form")
  data <- check_data(data=data, p=p)
  d <- ncol(data)
  n_obs <- nrow(data)
  weightfun_pars <- check_weightfun_pars(data=data, p=p, d=d, M=M, weight_function=weight_function,
                                         weightfun_pars=weightfun_pars, cond_dist=cond_dist)
  check_constraints(data=data, p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    parametrization=parametrization, identification="reduced_form",
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=NULL)
  npars <- n_params(p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                    B_constraints=NULL, identification="reduced_form")

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
    mu_scale2 <- vapply(1:d, function(i1) sd(data[,i1]), numeric(1))
  } else if(length(mu_scale2) != d | any(mu_scale2 <= 0)) {
    stop("Argument mu_scale2 must be strictly positive vector with length d")
  }
  if(missing(omega_scale)) {
    omega_scale <- vapply(1:d, function(i1) var(stats::ar(data[,i1], order.max=10)$resid, na.rm=TRUE), numeric(1))
  } else if(!(length(omega_scale) == d & all(omega_scale > 0))) {
    stop("omega_scale must be numeric vector with length d and positive elements")
  }
  if(missing(B_scale)) {
    B_scale <- vapply(1:d, function(i1) var(stats::ar(data[,i1], order.max=10)$resid, na.rm=TRUE), numeric(1))
  } else if(!(length(B_scale) == d && all(B_scale > 0))) {
    stop("B_scale must be numeric vector with length d and strictly positive elements")
  }
  if(missing(weight_scale)) {
    if(weight_function == "logistic" || weight_function == "exponential") {
      weight_scale <- c(mean(data[,weightfun_pars[1]]), 3*sd(data[,weightfun_pars[1]]), 8*sd(data[,weightfun_pars[1]]))
    } else if(weight_function == "mlogit") {
      weight_scale <- c(mean(data[,weightfun_pars[[1]]]), 8*sd(data[,weightfun_pars[[1]]]), 8*sd(data[,weightfun_pars[[1]]]))
    } else if(weight_function == "threshold") {
      quants <- unname(quantile(data[,weightfun_pars[1]], probs=c(0.2, 0.8)))
      weight_scale <- c(quants[1], quants[2])
    } else {
      weight_scale <- NULL # Not used
    }
  } else {
    if(weight_function %in% c("logistic", "exponential")) {
      if(length(weight_scale) != 3 || !is.vector(weight_scale) || !is.numeric(weight_scale) ||
         weight_scale[2] <= 0 || weight_scale[3] <= 0) {
        stop("For logistic and exponential weight functions, the argument weight_scale should be a length three numeric vector
              with strictly positive second and third elements")
      }
    } else if(weight_function == "mlogit") {
      if(length(weight_scale) != 2 || !is.vector(weight_scale) || !is.numeric(weight_scale) || weight_scale[2] <= 0) {
        stop("For mlogit weight function, the argument weight_scale should be a length two numeric vector
              with strictly positive second element")
      }
    } else if(weight_function == "threshold") {
      if(length(weight_scale) != 2 || !is.numeric(weight_scale)) {
        stop("For threshold weight function, the argument weight_scale should be a length two numeric vecor")
      }
      if(weight_scale[1] < min(data[,weightfun_pars[[1]]])) {
        warning("Lower bound for thresholds was set smaller than the smallest observation of the switching variable.
                Using lower bound set to the smallest observation.")
        weight_scale[1] <- min(data[,weightfun_pars[[1]]])
      }
      if(weight_scale[2] > max(data[,weightfun_pars[[1]]])) {
        warning("Upper bound for thresholds was set larger than the largest observation of the switching variable.
                Using upper bound set to the largest observation.")
        weight_scale[2] <- max(data[,weightfun_pars[[1]]])
      }
    }
  }

  stopifnot(pre_smart_mu_prob >= 0 && pre_smart_mu_prob <= 1)
  if(length(ar_scale) != 1 | ar_scale <= 0) {
    stop("ar_scale must be positive and have length one")
  }
  if(upper_ar_scale > 1) {
    stop("upper_ar_scale should be at most one")
  } else if(p*d > 40) {
    if(upper_ar_scale > 1 - p*d/150) {
      upper_ar_scale <- max(1 - p*d/150, 0.05)
    }
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
    n_attempts <- 20
    G <- numeric(0)
    for(i1 in 1:n_attempts) {
      inds <- replicate(popsize, random_ind(p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                            cond_dist=cond_dist, AR_constraints=AR_constraints,
                                            mean_constraints=mean_constraints,
                                            weight_constraints=weight_constraints,
                                            force_stability=is.null(AR_constraints),
                                            mu_scale=mu_scale, mu_scale2=mu_scale2,
                                            omega_scale=omega_scale, B_scale=B_scale, ar_scale=ar_scale,
                                            weight_scale=weight_scale, ar_scale2=ar_scale2))

      ind_loks <- vapply(1:popsize, function(i2) loglikelihood(data=data, p=p, M=M, params=inds[,i2],
                                                               weight_function=weight_function, weightfun_pars=weightfun_pars,
                                                               cond_dist=cond_dist, parametrization="mean",
                                                               identification="reduced_form",
                                                               AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                                               weight_constraints=weight_constraints, B_constraints=NULL,
                                                               to_return="loglik", check_params=TRUE,
                                                               minval=minval, alt_par=TRUE), numeric(1))
      G <- cbind(G, inds[, ind_loks > minval]) # Take good enough individuals
      if(ncol(G) >= popsize) {
        G <- G[, 1:popsize]
        break
      } else if(i1 == n_attempts) {
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
      tryCatch(check_params(data=data, p=p, M=M, d=d, params=ind, weight_function=weight_function,
                            weightfun_pars=weightfun_pars, cond_dist=cond_dist, parametrization=parametrization,
                            identification="reduced_form", AR_constraints=AR_constraints,
                            mean_constraints=mean_constraints, weight_constraints=weight_constraints, B_constraints=NULL),
               error=function(e) stop(paste("Problem with individual", i1, "in the initial population: "), e))
      if(parametrization == "intercept") {
        ind <- change_parametrization(p=p, M=M, d=d, params=ind, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                      cond_dist=cond_dist, identification="reduced_form", AR_constraints=AR_constraints,
                                      mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                      B_constraints=NULL, change_to="mean")
      }
      if(cond_dist == "ind_Student") { # Sort and sign change columns so that the first row of B_1 is positive and in a decreasing order.
        ind <- sort_impactmats(p=p, M=M, d=d, params=ind, weight_function=weight_function, weightfun_pars=weightfun_pars,
                              cond_dist=cond_dist, AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                              weight_constraints=weight_constraints)
        ind <- change_parametrization(p=p, M=M, d=d, params=ind, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                      cond_dist=cond_dist, identification="reduced_form", AR_constraints=AR_constraints,
                                      mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                      B_constraints=NULL, change_to="alt") # Change B_1^*,..,B_M^*
      }
      if(is.null(AR_constraints) && is.null(mean_constraints) && is.null(weight_constraints)) {
        initpop[[i1]] <- sort_regimes(p=p, M=M, d=d, params=ind, weight_function=weight_function, weightfun_pars=weightfun_pars,
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
    generations[, , i1] <- G

    # Calculate log-likelihoods and fitness inheritance
    if(i1 == 1) {
      # No fitness inheritance
      for(i2 in 1:popsize) {
        loks_and_tw <- loglikelihood(data=data, p=p, M=M, params=G[,i2], weight_function=weight_function,
                                     weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                                     parametrization="mean", identification="reduced_form",
                                     AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                     weight_constraints=weight_constraints, B_constraints=NULL,
                                     to_return="loglik_and_tw", check_params=TRUE, minval=minval, alt_par=TRUE)
        fill_lok_and_red(i1, i2, loks_and_tw)
      }
    } else {
      # Proportional fitness inheritance: individual has 50% change to inherit fitness if it's a result of crossover.
      # Variable "I" tells the proportions of parent material.
      I2 <- rep(I, each=2)
      which_did_co <- which(1 - which_not_co == 1)
      if(length(which_did_co) > 0) {
        which_inherit <- sample(x=which_did_co, size=round(0.5*length(which_did_co)), replace=FALSE)
      } else {
        which_inherit <- numeric(0)
      }
      # survivor_liks holds the parent loglikelihood values: for odd number they are (index, index+1) and for even (index-1, index)
      if(length(which_inherit) > 0 & abs(max_lik - mean_lik) > abs(0.03*mean_lik)) { # No inheritance if massive mutations
        for(i2 in which_inherit) {
          if(i2 %% 2 == 0) {
            logliks[i1, i2] <- ((npars - I2[i2])/npars)*survivor_liks[i2-1] + (I2[i2]/npars)*survivor_liks[i2]
            redundants[i1, i2] <- max(c(survivor_redundants[i2-1], survivor_redundants[i2]))
          } else {
            logliks[i1, i2] <- (I2[i2]/npars)*survivor_liks[i2] + ((npars - I2[i2])/npars)*survivor_liks[i2+1]
            redundants[i1, i2] <- max(c(survivor_redundants[i2], survivor_redundants[i2+1]))
          }
        }
        which_no_inherit <- (1:popsize)[-which_inherit]
      } else {
        which_no_inherit <- 1:popsize
      }
      for(i2 in which_no_inherit) { # Calculate the rest log-likelihoods
        if((mutate[i2] == 0 & which_not_co[i2] == 1) | all(H[,i2] == G[,i2])) { # If nothing changed
          logliks[i1, i2] <- survivor_liks[i2]
          redundants[i1, i2] <- survivor_redundants[i2]
        } else {
          if(stat_mu == TRUE & mutate[i2] == 1) { # Stability condition satisfied
            loks_and_tw <- tryCatch(loglikelihood(data=data, p=p, M=M, params=G[,i2], weight_function=weight_function,
                                                  weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                                                  parametrization="mean", identification="reduced_form",
                                                  AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                                  weight_constraints=weight_constraints, B_constraints=NULL,
                                                  to_return="loglik_and_tw", check_params=FALSE, minval=minval, alt_par=TRUE),
                                    error=function(e) minval)
          } else {
            loks_and_tw <- tryCatch(loglikelihood(data=data, p=p, M=M, params=G[,i2], weight_function=weight_function,
                                                  weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                                                  parametrization="mean", identification="reduced_form",
                                                  AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                                  weight_constraints=weight_constraints, B_constraints=NULL,
                                                  to_return="loglik_and_tw", check_params=TRUE, minval=minval, alt_par=TRUE),
                                    error=function(e) minval)
          }
          fill_lok_and_red(i1, i2, loks_and_tw)
        }
      }
    }

    # Take care of individuals that are not good enough + calculate the numbers redundant regimes
    if(anyNA(logliks[i1,])) {
      which_inds_na <- which(is.na(logliks[i1,]))
      logliks[i1, which_inds_na] <- minval
      warning(paste0("Removed NaNs from the logliks of inds ", which_inds_na, " in the GA round ", i1,
                     " with the seed ", seed, "."))
    }
    logliks[i1, which(logliks[i1,] < minval)] <- minval
    redundants[i1, which(logliks[i1,] <= minval)] <- M

    ## Selection and the reproduction pool ##
    if(length(unique(logliks[i1,])) == 1) {
      choosing_probs <- rep(1, popsize) # If all individuals are the same, the surviving probability weight is 1.
    } else {
      T_values <- logliks[i1,] + abs(min(logliks[i1,])) # Function T as described by Dorsey R. E. ja Mayer W. J., 1995
      T_values <- T_values/(1 + regime_force_scale*redundants[i1,]) # Favor individuals with least number of redundant regimes
      choosing_probs <- T_values/sum(T_values) # The surviving probability weights
    }

    # Draw popsize individuals with replacement and form the reproduction pool H.
    survivors <- sample(1:popsize, size=popsize, replace=TRUE, prob=choosing_probs)
    H <- G[,survivors]

    # Calculate mean and max log-likelihood of the survivors
    survivor_liks <- logliks[i1, survivors]
    survivor_redundants <- redundants[i1, survivors]
    max_lik <- max(survivor_liks)
    mean_lik <- mean(survivor_liks)
    if(max_lik == mean_lik) mean_lik <- mean_lik + 0.1 # +0.1 to avoid dividing by zero when all the individuals are the same

    ## Cross-overs ##
    # Individually adaptive cross-over rates as described by Patnaik and Srinivas (1994) with the modification of
    # setting the crossover rate to be at least 0.4 for all individuals (so that the best genes mix in the population too).
    indeces <- seq(from=1, to=popsize - 1, by=2)
    parent_max <- vapply(indeces, function(i2) max(survivor_liks[i2], survivor_liks[i2+1]), numeric(1))
    co_rates <- vapply(1:length(indeces), function(i2) max(min((max_lik - parent_max[i2])/(max_lik - mean_lik), 1), 0.4), numeric(1))

    # Do the crossovers
    which_co <- rbinom(n=popsize/2, size=1, prob=co_rates)
    I <- round(runif(n=popsize/2, min=0.5 + 1e-16, max=npars - 0.5 - 1e-16)) # Break points
    H2 <- vapply(1:(popsize/2), function(i2) {
      i3 <- indeces[i2]
      if(which_co[i2] == 1) {
        c(c(H[1:I[i2], i3], H[(I[i2]+1):npars, i3+1]), c(H[1:I[i2], i3+1], H[(I[i2]+1):npars, i3]))
      } else {
        c(H[,i3], H[,i3+1])
      }
    }, numeric(2*npars))
    H2 <- matrix(H2, nrow=npars, byrow=FALSE)

    # Get the best individual so far and check for reduntant regimes
    best_index0 <- which(logliks == max(logliks), arr.ind=TRUE)
    best_index <- best_index0[order(best_index0[,1], decreasing=FALSE)[1],] # First generation when the best loglik occurred
    best_ind <- generations[, best_index[2], best_index[1]]
    best_mw <- loglikelihood(data=data, p=p, M=M, params=best_ind, weight_function=weight_function, weightfun_pars=weightfun_pars,
                             cond_dist=cond_dist, parametrization="mean", identification="reduced_form", AR_constraints=AR_constraints,
                             mean_constraints=mean_constraints, weight_constraints=weight_constraints, B_constraints=NULL,
                             to_return="tw", check_params=FALSE, minval=minval, alt_par=TRUE)

    # Which regimes are wasted:
    which_redundant <- which(vapply(1:M, function(i2) sum(best_mw[,i2] > red_criteria[1]) < red_criteria[2]*n_obs, logical(1)))

    # Keep track of "the alternative best individual" that has (weakly) less reduntant regimes than the current best one.
    if(length(which_redundant) <= length(which_redundant_alt)) {
      alt_ind <- best_ind
      which_redundant_alt <- which_redundant
    }

    ## Mutations ##
    which_not_co <- rep(1 - which_co, each=2)
    if(abs(max_lik - mean_lik) <= abs(0.03*mean_lik)) {
      mu_rates <- rep(0.7, popsize) # Massive mutations if converging
    } else {
      # Individually adaptive mutation rates, Patnaik and Srinivas (1994); we only mutate those who did not crossover.
      mu_rates <- 0.5*vapply(1:popsize, function(i2) min(which_not_co[i2], (max_lik - survivor_liks[i2])/(max_lik - mean_lik)),
                             numeric(1))
    }

    mutate <- rbinom(n=popsize, size=1, prob=mu_rates)
    which_mutate <- which(mutate == 1)
    pre_smart_mu <- runif(1, min=1e-6, max=1-1e-6) < pre_smart_mu_prob
    ar_scale <- runif(1, min=1e-6, max=upper_ar_scale - 1e-6) # Random AR scale
    if(i1 <= smart_mu & length(which_mutate) >= 1 & !pre_smart_mu) { # Random mutations
      if(!is.null(AR_constraints) || runif(1, min=1e-6, max=1 - 1e-6) > 0.5) { # Does not always satisfy the stability conditions
        stat_mu <- FALSE
      } else { # Force stability condition with an algorithm (slower but can skip stability check), Ansley and Kohn (1986)
        stat_mu <- TRUE
      }
      H2[,which_mutate] <- vapply(1:length(which_mutate), function(x) random_ind(p=p, M=M, d=d,
                                                                                 weight_function=weight_function,
                                                                                 weightfun_pars=weightfun_pars,
                                                                                 cond_dist=cond_dist,
                                                                                 AR_constraints=AR_constraints,
                                                                                 mean_constraints=mean_constraints,
                                                                                 weight_constraints=weight_constraints,
                                                                                 force_stability=stat_mu,
                                                                                 mu_scale=mu_scale,
                                                                                 mu_scale2=mu_scale2,
                                                                                 omega_scale=omega_scale,
                                                                                 B_scale=B_scale,
                                                                                 ar_scale=ar_scale,
                                                                                 weight_scale=weight_scale,
                                                                                 ar_scale2=ar_scale2), numeric(npars))

    } else if(length(which_mutate) >= 1) { # Smart mutations
      stat_mu <- FALSE

      # If redundant regimes - smart mutate more
      if(length(which_redundant) >= 1) {
        mu_rates <- vapply(1:popsize, function(i2) which_not_co[i2]*max(0.1, mu_rates[i2]), numeric(1))
        mutate0 <- rbinom(n=popsize, size=1, prob=mu_rates)
        which_mutate0 <- which(mutate0 == 1)
        if(length(which_mutate0) > length(which_mutate)) {
          mutate <- mutate0
          which_mutate <- which_mutate0
        }
      }

      # Mutating accuracy
      accuracy <- abs(rnorm(length(which_mutate), mean=10, sd=15))

      ## 'Smart mutation': mutate close to a well fitting individual. We obviously don't mutate close to
      # redundant regimes but draw them at random ('rand_to_use' in what follows).
      if(!is.null(AR_constraints) | !is.null(mean_constraints) | !is.null(weight_constraints) |
         length(which_redundant) <= length(which_redundant_alt) | runif(1) > 0.5) {
        # The first option for smart mutations: smart mutate to 'alt_ind' which is the best fitting individual
        # with the least redundant regimes.
        # Note that best_ind == alt_ind when length(which_redundant) <= length(which_redundant_alt).
        ind_to_use <- alt_ind
        rand_to_use <- which_redundant_alt
      } else {
        # Alternatively, if there exists an alternative individual with strictly less redundant regimes
        # than in the best_ind, a "regime combining procedure" might take place: take a redundant regime
        # of the best_ind and replace it with a nonredundant regime taken from alt_ind. Then, do smart
        # mutation close to this new individual. For simplicity, regime combining is not considered for
        # models imposing constraints.

        # We want to take such nonredundant regime from alt_ind that is not similar to the nonredundant
        # regimes of best_ind. In order to choose such regime, we compare all nonredundant regimes of
        # best_ind to all nonredundant regimes of alt_ind. Then, we choose the nonredundant regime of
        # alt_ind which has the largest "distance" to the closest regime of best_ind. Note that weight
        # nor distribution parameters are included the regime combination procedure.

        # Choose regime of best_ind to be changed
        which_to_change <- which_redundant[1]

        # Pick the nonredundant regimes of best_ind and alt_ind
        n_regpars <- p*d^2 + d + ifelse(cond_dist == "ind_Student", d^2, d*(d+1)/2)
        non_red_regs_best <- vapply((1:M)[-which_redundant], function(i2) pick_regime(p=p, M=M, d=d, params=best_ind, cond_dist=cond_dist,
                                                                                      m=i2), numeric(n_regpars))

        alt_regs_to_pick  <- 1:M
        if(length(which_redundant_alt) != 0) { # Needed because (1:M)[-numeric(0)] = integer(0)
          alt_regs_to_pick <- alt_regs_to_pick[-which_redundant_alt]
        }
        non_red_regs_alt <- vapply(alt_regs_to_pick, function(i2) pick_regime(p=p, M=M, d=d, params=alt_ind, cond_dist=cond_dist,
                                                                              m=i2), numeric(n_regpars))

        # Calculate the "distances" between the nonredundant regimes
        #  Row for each non-red-reg-best and column for each non-red-reg-al:
        dist_to_regime <- matrix(nrow=ncol(non_red_regs_best), ncol=ncol(non_red_regs_alt))
        for(i2 in 1:nrow(dist_to_regime)) {
          dist_to_regime[i2,] <- vapply(1:ncol(non_red_regs_alt), function(i3) regime_distance(regime_pars1=non_red_regs_best[,i2],
                                                                                               regime_pars2=non_red_regs_alt[,i3]),
                                        numeric(1))
        }


        # Which alt_ind regime, i.e. column should be used? Choose the one that with largest distance to the closest
        # regime avoid duplicating similar regimes
        which_reg_to_use <- which(apply(dist_to_regime, 2, min) == max(apply(dist_to_regime, 2, min)))[1]

        # The obtain the regime of alt_ind that is used to replace a redundant regime in best_ind
        reg_to_use <- non_red_regs_alt[,which_reg_to_use]

        # Change the chosen regime of best_ind to be the one chosen from alt_ind
        ind_to_use <- change_regime(p=p, M=M, d=d, params=best_ind, cond_dist=cond_dist, m=which_to_change, regime_pars=reg_to_use)

        # Should some regimes still be random?
        rand_to_use <- which_redundant[which_redundant != which_to_change]
      }
      # Do the smart mutations
      H2[,which_mutate] <- vapply(1:length(which_mutate), function(i2) smart_ind(p=p, M=M, d=d,
                                                                                 params=ind_to_use,
                                                                                 weight_function=weight_function,
                                                                                 weightfun_pars=weightfun_pars,
                                                                                 cond_dist=cond_dist,
                                                                                 AR_constraints=AR_constraints,
                                                                                 mean_constraints=mean_constraints,
                                                                                 weight_constraints=weight_constraints,
                                                                                 accuracy=accuracy[i2],
                                                                                 which_random=rand_to_use,
                                                                                 mu_scale=mu_scale,
                                                                                 mu_scale2=mu_scale2,
                                                                                 omega_scale=omega_scale,
                                                                                 B_scale=B_scale,
                                                                                 ar_scale=ar_scale,
                                                                                 ar_scale2=ar_scale2), numeric(npars))
    }

    # Sort components according to the transition weight parameters (for some weight functions). No sorting if constraints are employed.
    if(is.null(AR_constraints) && is.null(mean_constraints) && is.null(weight_constraints)) {
      H2 <- vapply(1:popsize, function(i2) sort_regimes(p=p, M=M, d=d, params=H2[,i2], weight_function=weight_function,
                                                        weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                                                        identification="reduced_form"), numeric(npars))
    }

    # Sort and sign change the columns of the impact matrices so that the first row of B_1 is positive and in a decreasing order.
    if(cond_dist == "ind_Student") {
      H2 <- vapply(1:popsize, function(i2) sort_impactmats(p=p, M=M, d=d, params=H2[,i2], weight_function=weight_function,
                                                           weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                                                           AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                                           weight_constraints=weight_constraints), numeric(npars))
    }

    # Save the results and set up new generation
    G <- H2
  }

  if(to_return == "best_ind") {
    ret <- best_ind
  } else {
    ret <- alt_ind
  }
  # # GA always optimizes with mean parametrization, and cond_dist="ind_Student" models by parametrizing impact matrices
  # # instead of B_1,...,B_M, with B_1,B_2*,...,B_M*, where B_m* = B_m - B_1 for m=2,...,M.
  # # Switch parametriaztion back to the original one:
  # if(cond_dist == "ind_Student") {
  #   ret <- change_parametrization(p=p, M=M, d=d, params=ret, weight_function=weight_function, weightfun_pars=weightfun_pars,
  #                                 cond_dist=cond_dist, identification="reduced_form", AR_constraints=AR_constraints,
  #                                 mean_constraints=mean_constraints, weight_constraints=weight_constraints,
  #                                 B_constraints=NULL, change_to="alt")
  # }

  # # GA always optimizes with mean parametrization
  # Return intercept parametrized estimate if parametrization=="intercept".
  if(parametrization == "mean") { # This is always the case with mean_constraints
    return(ret)
  } else {
    return(change_parametrization(p=p, M=M, d=d, params=ret, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                  cond_dist=cond_dist, identification="reduced_form", AR_constraints=AR_constraints,
                                  mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                  B_constraints=NULL, change_to="intercept"))
  }
}
