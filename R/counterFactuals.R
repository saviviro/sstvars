#' @title Compute historical counterfactual for structural STVAR models.
#'
#' @description \code{cfact_hist} computes historical counterfactual for structural STVAR models.
#'
#' @inheritParams simulate.stvar
#' @param type a character string indicating the type of counterfactual to be computed: should the path of the policy
#'  variable be fixed to some hypothetical path (\code{type="fixed_path"}) in given points of time or should the responses
#'  of the policy variable to lagged and contemporaneous movements of some given variable be muted (\code{type="muted_response"})?
#'  See details for more information.
#' @param policy_var a positive integer between \eqn{1} and \eqn{d} indicating the index of the policy variable considered in the
#'  counterfactual scenario.
#' @param mute_var a positive integer between \eqn{1} and \eqn{d} indicating the index of the variable to whose movements the policy variable
#'  specified in the argument \code{policy_var} should not react to in the counterfactual scenario. This indicates also the index of the shock
#'  to which the policy variable should not react to. It is assumed that \code{mute_var != policy_var}. This argument is only used when
#'  \code{type="muted_response"}.
#' @param cfact_start a positive integer between \eqn{1} and \eqn{T} indicating the starting period for the counterfactual behavior
#'  of the specified policy variable. It is assumed that the observed time series in indexed as \eqn{y_{t-p+1},...,y_{0},y_1,...,y_T},
#'  i.e., that the first \eqn{p} observations are the initial values, and the "time period one" observation is the \eqn{p+1}th row in
#'  the data matrix.
#' @param cfact_end a positive integer between \code{cfact_start} and \eqn{T} indicating the ending period for the counterfactual
#'  behavior of the specified policy variable.
#' @param cfact_path a numeric vector of length \code{cfact_end-cfact_start+1} indicating the hypothetical path of the policy variable
#'  specified in the argument \code{policy_var}. This argument is only used when \code{type="fixed_path"}.
#' @details Two types of historical counterfactuals are accommodated where in given historical points of time
#'  either (1) the policy variable of interest takes some hypothetical path (\code{type="fixed_path"}), or (2)
#'  its responses to lagged and contemporaneous movements of some given variable are shut off (\code{type="muted_response"}).
#'  In both cases, the counterfactual scenarios are simulated by creating hypothetical shocks to the policy variable of interest
#'  that yield the counterfactual outcome. This approach has the appealing feature that the counterfactual deviations from the
#'  policy reaction function are treated as policy surprises, allowing them to propagate normally, so that the dynamics of the model
#'  are not, per se, tampered but just the policy surprises are.
#'
#'  \strong{Important:} This function assumes that when the policy variable of interest is the \eqn{i_1}th variable, the shock
#'  to it that is manipulated is the \eqn{i_1}th shock. This should be automatically satisfied for recursively identified models,
#'  whereas for model identified by heteroskedasticity or non-Gaussianity, the ordering of the shocks can be generally changed
#'  without loss of generality with the function \code{reorder_B_columns}. In Type (2) counterfactuals it is additionally assumed
#'  that, if the variable to whose movements the policy variable should not react to is the \eqn{i_2}:th variable, the shock to it
#'  is the \eqn{i_2}th shock. If it is not clear whether the $i_2$th shock of interest can be interpreted as a shock to a variable
#'  (but has a broader definition such as "a demand shock"), the Type (2)counterfactual scenario is interpreted as follows: the $i_1$th
#'  variable does not react to lagged movements of the $i_2$th variable nor to the $i_2$th shock.
#'
#'  See the seminal paper of Bernanke et al (1997) for discussing about the "Type (1)" counterfactuals and
#'  Kilian and Lewis (2011) for discussion about the "type (2)" counterfactuals. See Kilian and Lütkepohl (2017), Section 4.3
#'  for further discussion about the historical counterfactuals. The literature cited about considers linear models, but it is
#'  explained in the vignette of this package how this function computes the historical counterfactuals for the STVAR models in
#'  a way that accommodates nonlinear time-varying dynamics.
#' @return Returns a class \code{'histdecomp'} list with the following elements:
#'   \describe{
#'     \item{FILL IN}{FILL IN}
#'     \item{stvar}{The original STVAR model object.}
#'  }
#' @seealso \code{\link{GIRF}}, \code{\link{GFEVD}}, \code{\link{linear_IRF}}, \code{\link{hist_decomp}}, \code{\link{fitSSTVAR}}
#' @references
#'  \itemize{
#'    \item Bernanke B., Gertler M., Watson M. 1997. Systematic monetary policy and the effects of oilprice shocks.
#'      \emph{Brookings Papers on Economic Activity}, 1, 91—142.
#'    \item Kilian L., Lewis L. 2011. Does the fed respond to oil price shocks? \emph{The Economic Journal}, \strong{121}:555.
#'    \item Kilian L., Lütkepohl H. 2017. Structural Vector Autoregressive Analysis. 1st edition.
#'      \emph{Cambridge University Press}, Cambridge.
#'  }
#' @examples
#' # Recursively identified logistic Student's t STVAR(p=3, M=2) model with the first
#' # lag of the second variable as the switching variable:
#' params32logt <- c(0.5959, 0.0447, 2.6279, 0.2897, 0.2837, 0.0504, -0.2188, 0.4008,
#'   0.3128, 0.0271, -0.1194, 0.1559, -0.0972, 0.0082, -0.1118, 0.2391, 0.164, -0.0363,
#'   -1.073, 0.6759, 3e-04, 0.0069, 0.4271, 0.0533, -0.0498, 0.0355, -0.4686, 0.0812,
#'   0.3368, 0.0035, 0.0325, 1.2289, -0.047, 0.1666, 1.2067, 7.2392, 11.6091)
#' mod32logt <- STVAR(gdpdef, p=3, M=2, params=params32logt, weight_function="logistic",
#'   weightfun_pars=c(2, 1), cond_dist="Student", identification="recursive")
#'
#' # FILL IN
#' @export

cfact_hist <- function(stvar, type=c("fixed_path", "muted_response"), policy_var=1, mute_var, cfact_start=1, cfact_end=1, cfact_path) {
  type <- match.arg(type)
  data <- stvar$data
  p <- stvar$model$p
  M <- stvar$model$M
  d <- ncol(data)
  T_obs <- nrow(data) - p
  cond_dist <- stvar$model$cond_dist
  parametrization <- stvar$model$parametrization
  identification <- stvar$model$identification
  weight_function <- stvar$model$weight_function
  weightfun_pars <- stvar$model$weightfun_pars

  ## Check the arguments
  if(!is.numeric(policy_var) || length(policy_var) != 1 || policy_var < 1 || policy_var > d || policy_var%%1 != 0) {
    stop("The argument policy_var should be a positive integer between 1 and d")
  } else if(!is.numeric(cfact_start) || length(cfact_start) != 1 || cfact_start < 1 || cfact_start > T_obs || cfact_start%%1 != 0) {
    stop("The argument cfact_start should be a positive integer between 1 and T")
  } else if(!is.numeric(cfact_end) || length(cfact_end) != 1 || cfact_end < cfact_start || cfact_end > T_obs || cfact_end%%1 != 0) {
    stop("The argument cfact_end should be a positive integer between cfact_start and T")
  }
  if(type == "fixed_path") {
    if(!is.numeric(cfact_path) || length(cfact_path) != cfact_end - cfact_start + 1) {
      stop("The argument cfact_path should be a numeric vector of length cfact_end-cfact_start+1")
    } else if(!is.numeric(cfact_path) || any(is.na(cfact_path))) {
      stop("The argument cfact_path should not contain NA values")
    }
  } else { # type == "muted_response"
    if(missing(mute_var)) {
      stop("The argument mute_var is missing with no default")
    } else if(!is.numeric(mute_var) || length(mute_var) != 1 || mute_var < 1 || mute_var > d || mute_var%%1 != 0) {
      stop("The argument mute_var should be a positive integer between 1 and d")
    } else if(policy_var == mute_var) {
      stop("The arguments policy_var and mute_var should not be equal")
    }
  }

  ## Obtain the parameter values in the "non constrained form", and also switch to reduced_form parameter vector:
  params <- reform_constrained_pars(p=p, M=M, d=d, params=stvar$params, weight_function=weight_function, cond_dist=cond_dist,
                                    identification=identification, AR_constraints=stvar$model$AR_constraints,
                                    mean_constraints=stvar$model$mean_constraints, weight_constraints=stvar$model$weight_constraints,
                                    B_constraints=stvar$model$B_constraints, other_constraints=NULL, weightfun_pars=weightfun_pars)

  ## Pick params
  if(stvar$model$parametrization == "mean") { # Change to intercapt parametrization
    params <- change_parametrization(p=p, M=M, d=d, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                     identification=identification, cond_dist=cond_dist, AR_constraints=NULL, mean_constraints=NULL,
                                     weight_constraints=NULL, B_constraints=NULL, change_to="intercept")
  }
  all_mu <- get_regime_means(p=p, M=M, d=d, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
                             cond_dist=cond_dist, parametrization="intercept", identification=identification,
                             AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL) # Not necessarily valid if allow_unstab
  all_phi0 <- pick_phi0(M=M, d=d, params=params)
  all_A <- pick_allA(p=p, M=M, d=d, params=params) # [d, d, p, M]
  all_Omegas <- pick_Omegas(p=p, M=M, d=d, params=params, cond_dist=cond_dist, identification=identification) # [d, d, M], B_m for ind_Stud and ind_skewed_t
  weightpars <- pick_weightpars(p=p, M=M, d=d, params=params, weight_function=weight_function, cond_dist=cond_dist, weightfun_pars=weightfun_pars)
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A) # The bold A matrices of the regimes to be used later
  distpars <- pick_distpars(d=d, params=params, cond_dist=cond_dist)

  # Structural pars (recursive just takes Cholesky and model identified by non-Gaussianity have B_m in all_Omegas)
  if(identification == "heteroskedasticity") {
    W <- pick_W(p=p, M=M, d=d, params=params, identification=identification)
    lambdas <- matrix(pick_lambdas(p=p, M=M, d=d, params=params, identification=identification), nrow=d, ncol=M-1)
  }

  ## Calculate required statistics that remain constant through the iterations
  if(weight_function == "relative_dens") {  # relative_dens weight function uses this
    Sigmas <- get_Sigmas(p=p, M=M, d=d, all_A=all_A, all_boldA=all_boldA, all_Omegas=all_Omegas)
    inv_Sigmas <- array(NA, dim=c(d*p, d*p, M)) # Store inverses of the (dpxdp) covariance matrices
    det_Sigmas <- numeric(M) # Store determinants of the (dpxdp) covariance matrices
    chol_Sigmas <- array(dim=c(d*p, d*p, M)) # Cholesky decompositions of the (dpxdp) covariance matrices
    for(m in 1:M) {
      chol_Sigmas[, , m] <- chol(Sigmas[, , m]) # Upper triangle
      inv_Sigmas[, , m] <- chol2inv(chol_Sigmas[, , m]) # Faster inverse
      det_Sigmas[m] <- prod(diag(chol_Sigmas[, , m]))^2 # Faster determinant
    }
  } else if(weight_function == "mlogit") {
    all_gamma_m <- matrix(weightpars, ncol=M-1)
    vars <- weightfun_pars[[1]]
    lags <- weightfun_pars[[2]]
    lowers <- (1:lags - 1)*d # We want add vars to each of these
    # Indices of switching variables in cbind(1, Y): we add +1 to the indices since the column of ones on the left,
    # and then the index is added to always account for the constant term.
    inds_of_switching_vars <- c(1, as.vector(matrix(lowers, nrow=length(vars), ncol=length(lowers), byrow=TRUE) + vars + 1))
  }

  ## Obtain the structural shocks recovered from the model
  if(cond_dist %in% c("ind_Student", "ind_skewed_t") && identification == "reduced_form") {
    identification_to_use <- "non-Gaussianity"
  } else if(cond_dist %in% c("Gaussian", "Student") && identification == "reduced_form") {
    identification_to_use <- "recursive"
  } else {
    identification_to_use <- identification
  }
  all_e_t <- get_residuals(data=stvar$data, p=p, M=M, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
                           cond_dist=cond_dist, parametrization=parametrization, identification=identification_to_use, AR_constraints=NULL,
                           mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL, standardize=TRUE, structural_shocks=TRUE,
                           penalized=stvar$penalized, penalty_params=stvar$penalty_params, allow_unstab=stvar$allow_unstab)

  ## Obtain the transition weights
  alpha_mt <- stvar$transition_weights # [T_obs, M]

  # Some functions to be used to obtain the transition weights during the counterfactual simulation
  if(weight_function == "relative_dens") {
    # Get log multivariate normal densities for calculating the transition weights
    get_logmvdvalues <- function(Y, i1) {
      vapply(1:M, function(m) -0.5*d*p*log(2*base::pi) - 0.5*log(det_Sigmas[m]) - 0.5*(crossprod(Y[i1,] - rep(all_mu[, m], p),
                                                                                                 inv_Sigmas[, , m])%*%(Y[i1,] - rep(all_mu[, m], p))),
             numeric(1))
    } # Returns M x 1 vector; transformed into a matrix in get_alpha_mt
  } else if(weight_function %in% c("logistic", "exponential", "threshold")) {
    # Nothing to do here, can use get_alpha_mt directly
  } else if(weight_function == "mlogit") {
    # Calculate the sub model regressions for calculating the transition weights
    get_regression_values <- function(Y, i1) {
      regressions_mt <- numeric(M) # Vector of zeros
      for(i2 in 1:ncol(all_gamma_m)) {
        regressions_mt[i2] <- crossprod(all_gamma_m[,i2], cbind(1, Y)[i1, inds_of_switching_vars])
      }
      regressions_mt
    }
  } else if(weight_function == "exogenous") {
    # Nothing to do here, uses exo_weights directly
  }

  ## Obtain the data in a convenient form:
  # i:th row denotes the vector \bold{y_{i-1}} = (y_{i-1},...,y_{i-p}) (dpx1), assuming the observed data is y_{-p+1},...,y_0,y_1,...,y_{T}.
  Y <- reform_data(stvar$data, p) # (T+1 x dp)
  #Y <- Y[1:T_obs, , drop=FALSE] # Last row removed; not needed when calculating something based on lagged observations

  ## Create a container for the counterfactual observations, including also the original data prior to the counterfactual period:
  cfact_data <- matrix(NA, nrow=nrow(data), ncol=d) # Note that the first p periods are for the initial values
  cfact_data[1:(cfact_start + p - 1),] <- data[1:(cfact_start + p - 1),]

  ## Create a container for the counterfactual transition weights, including also the original ones prior to the counterfactual period:
  cfact_alpha_mt <- matrix(NA, nrow=T_obs, ncol=M) # Note that there are no initial values here, so the indexing is different to cfact_data
  cfact_alpha_mt[1:(cfact_start - 1),] <- alpha_mt[1:(cfact_start - 1),] # Original transition weights

  ## Create a container for the counterfactual shocks, including also the original ones prior to the counterfactual period:
  cfact_e_t <- matrix(NA, nrow=T_obs, ncol=d) # Note that there are no initial values, so the indexing is different to cfact_data
  cfact_e_t[1:(cfact_start - 1),] <- all_e_t[1:(cfact_start - 1),] # Original transition weights

  ## Simulate the counterfactual observations
  for(t in cfact_start:T_obs) { # Loop through the time periods starting from the beginning of the counterfactual
    t_row_in_data <- t + p # The row in the data matrix (but not in the matrix Y)

    # Calculate the time period t transition weights
  }
}


#' @title Compute the intercept \eqn{\phi_{y,t}=\sum_{m=1}^M \alpha_{mt} \phi_{m}} parameter value for a single time period
#'
#' @description \code{get_next_phi_yt} computes the intercept \eqn{\phi_{y,t}=\sum_{m=1}^M \alpha_{mt} \phi_{m}} paramater
#'  value for a single time period based on the regime intercepts and transition weights.
#'
#' @param all_phi0 a \eqn{(d \times M)} matrix such that the \eqn{m}th column contains the intercept parameters of the \eqn{m}th regime.
#' @param alpha_mt an \eqn{(M \times 1)} vector containing the time period \eqn{t} transition weights.
#' @details This is used in simulation of the counterfactual scenarios.
#' @return Returns a \eqn{(d \times 1)} vector of the intercept parameter values for the time period \eqn{t}.
#' @keywords internal

get_phi_yt <- function(all_phi0, alpha_mt) {
  # Calculate the intercept parameter value for the time period t
  as.numeric(all_phi0%*%alpha_mt) # [d, M] x [M, 1] = [d, 1]
}
