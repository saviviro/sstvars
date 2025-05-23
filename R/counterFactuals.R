
#' @title Compute the intercept \eqn{\phi_{y,t}=\sum_{m=1}^M \alpha_{mt} \phi_{m}} parameter value for a single time period
#'
#' @description \code{get_phi_yt} computes the intercept \eqn{\phi_{y,t}=\sum_{m=1}^M \alpha_{mt} \phi_{m}} parameter
#'  value for a single time period based on the regime intercepts and transition weights.
#'
#' @param all_phi0 a \eqn{(d \times M)} matrix such that the \eqn{m}th column contains the intercept parameters of the \eqn{m}th regime.
#' @param alpha_mt an \eqn{(M \times 1)} vector containing the time period \eqn{t} transition weights.
#' @details This is used in simulation of the counterfactual scenarios.
#' @return Returns the \eqn{(d \times 1)} vector of the intercept parameter values for the time period \eqn{t}.
#' @keywords internal

get_phi_yt <- function(all_phi0, alpha_mt) {
  # Calculate the intercept parameter value for the time period t
  as.numeric(all_phi0%*%alpha_mt) # [d, M] x [M, 1] = [d, 1]
}


#' @title Compute the autoregression matrices \eqn{A_{y,t,i}\equiv \sum_{m=1}^M\alpha_{m,t}A_{m,i}} for all lags \eqn{i=1,...,p}
#'  for a single time period
#'
#' @description \code{get_allA_yti} computes the autoregression matrices \eqn{A_{y,t,i}\equiv \sum_{m=1}^M\alpha_{m,t}A_{m,i}}, for all lags \eqn{i=1,...,p}
#'  for a single time period, based on the regime autoregression matrices and transition weights.
#'
#' @param all_A  4D array containing the coefficient matrices of all regimes so that coefficient matrix
#'  \eqn{A_{m,i}} can be obtained by choosing \code{[, , i, m]} (as obtained from \code{pick_allA}).
#' @param alpha_mt an \eqn{(M \times 1)} vector containing the time period \eqn{t} transition weights.
#' @details This is used in simulation of the counterfactual scenarios.
#' @return Returns the 3D array containing the coefficient matrices for the given time period so that the lag \eqn{i} coefficient matrix \eqn{A_{y,t,i}}
#' can be obtained by choosing \code{[, , i]}.
#' @keywords internal

get_allA_yti <- function(all_A, alpha_mt) {
  # Calculate the autoregression matrices A_{y,t,i} for all lags i=1,...,p
  # all_A: [d, d, p, M], alpha_mt: length M
  d <- dim(all_A)[1]
  p <- dim(all_A)[3]
  M <- dim(all_A)[4]
  all_A_yti <- array(0, dim=c(d, d, p)) # Pre-allocate [d, d, p]

  for(i1 in seq_len(p)) {
    slice_i1 <- all_A[, ,i1 , , drop=FALSE] # Extract [d, d, M] slice for lag i
    all_A_yti[, , i1] <- apply(slice_i1, MARGIN=c(1, 2),  FUN=function(x) sum(x*alpha_mt)) # At each (row, col), sum over regimes m with weights
  }

  all_A_yti
}


#' @title Compute the conditional mean \eqn{\mu_{y,t}=\phi_{y,t} + \sum_{i=1}^pA_{y,t,i}y_{t-i}} for a single time period
#'
#' @description \code{get_mu_yt} computes the conditional mean \eqn{\mu_{y,t}=\phi_{y,t} + \sum_{i=1}^pA_{y,t,i}y_{t-i}} for a single time period
#'  based on the intercepts, AR matrices, and the vector of lagged observations.
#'
#' @param phi_yt a \eqn{(d \times M)} matrix such that the \eqn{m}th column contains the intercept parameters of the \eqn{m}th regime.
#' @param all_A_yti a 3D array containing the coefficient matrices for the given time period so that the lag \eqn{i} coefficient matrix
#'  \eqn{A_{y,t,i}} can be obtained by choosing \code{[, , i]}.
#' @param bold_y_t_minus_1 a \eqn{(dp \times 1)} vector \eqn{\boldsymbol{y}_{t-1}=(y_{t-1},...,y_{t-p})} containing the lagged observations
#'  for the time period \eqn{t}.
#' @details This is used in simulation of the counterfactual scenarios.
#' @return Returns the \eqn{(d \times 1)} vector of the conditional mean for the time period \eqn{t}.
#' @keywords internal

get_mu_yt <- function(phi_yt, all_A_yti, bold_y_t_minus_1) {
  d <- dim(all_A_yti)[1]
  p <- dim(all_A_yti)[3]

  # Create AR matrix for all lags
  big_A <- matrix(0, nrow=d, ncol=d*p) # [d, dp]
  for(i1 in seq_len(p)) {
    big_A[, ((i1 - 1)*d + 1):(i1*d)] <- all_A_yti[, , i1] # Fill the i-th block of big_A
  }

  as.vector(phi_yt + big_A%*%bold_y_t_minus_1) # [d, 1] + [d, dp] x [dp, 1] = [d, 1]
}


#' @title Compute the impact matrix \eqn{B_{y,t}} for a single time period
#'
#' @description \code{get_B_yt} computes the impact matrix \eqn{B_{y,t}}. For \code{"ind_Student"} and \code{"ind_skewed_t"} models
#'  \eqn{B_{y,t}=\sum_{m=1}^M\alpha_{m,t}B_m}. For models identified by heteroskedasticity \eqn{B_{y,t}=W\sqrt{\sum_{m=1}^M\alpha_{m,t}\Lambda_m}}.
#'  For recursive identification \eqn{B_{y,t}} is obtained from the Cholesky decomposition of the conditional covariance matrix.
#'
#' @inheritParams loglikelihood
#' @param all_Omegas a 3D array such that the covariance matrix (or impact matrix \eqn{B_m}) of the \eqn{m}th regime is obtained from \code{all_Omegas[, , m]}.
#' @param alpha_mt an \eqn{(M \times 1)} vector containing the time period \eqn{t} transition weights.
#' @param W a \eqn{(d \times d)} matrix containing the matrix \eqn{W} for models identified by heteroskedasticity (as returned by \code{pick_W}).
#' @param lambdas a \eqn{(d(M-1)\times 1)} vector \eqn{\lambda_2,...,\lambda_M} for models identified by heteroskedasticity (as returned by \code{pick_lambdas}).
#' @details This is used in simulation of the counterfactual scenarios.
#' @return Returns the \eqn{(d \times d)} impact matrix for the time period \eqn{t}.
#' @keywords internal

get_B_yt <- function(all_Omegas, alpha_mt, W, lambdas, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t"),
                     identification=c("reduced_form", "recursive", "non-Gaussianity", "heteroskedasticity")) {
  cond_dist <- match.arg(cond_dist)
  identification <- match.arg(identification)

  if(cond_dist %in% c("ind_Student", "ind_skewed_t")) {
    #d <- dim(all_Omegas)[1]
    #return(matrix(rowSums(vapply(1:dim(all_Omegas)[3], function(m) alpha_mt[m]*as.vector(all_Omegas[, , m]),
    #                             numeric(d*d))), nrow=d, ncol=d)) # Faster, but not as readable
    return(apply(sweep(all_Omegas, MARGIN=3, STATS=alpha_mt, FUN="*"),
                 MARGIN=c(1, 2), FUN=sum)) # Multiply each slice B_m by its weight and sum over the regimes
  } else if(identification == "heteroskedasticity") {
    d <- dim(all_Omegas)[1]
    M <- dim(all_Omegas)[3]
    lambdas <- matrix(lambdas, ncol=M - 1)
    tmp <- array(dim=c(d, d, M)) # Store alpha_mt[m]*Lambda_m
    tmp[, , 1] <- alpha_mt[1]*diag(d) # m=1, Lambda = I_d
    for(m in 2:M) {
      tmp[, , m] <- alpha_mt[m]*diag(lambdas[, m - 1])
    }
    return(W%*%sqrt(apply(tmp, MARGIN=1:2, FUN=sum))) # Calculate B_yt as in Virolainen 2025 (JBES)
  } else { # Recursive identification, B_yt calculated from the conditional covariance matrix
    Omega_yt <- apply(sweep(all_Omegas, MARGIN=3, STATS=alpha_mt, FUN="*"),
                      MARGIN=c(1, 2), FUN=sum) # Multiply each slice Omega_m by its weight and sum over the regimes
    return(t(chol(Omega_yt))) # Cholesky decomposition of the conditional covariance matrix, zeros in the upper triangle
  }
}


#' @title Compute the observation \eqn{y_t=\mu_{y,t} + B_{y,t}e_t} for a single time period
#'
#' @description \code{get_y_t} computes the observation \eqn{y_t=\mu_{y,t} + B_{y,t}e_t} for a single time period
#'  based on the conditional mean, impact matrix, and shock vector.
#'
#' @param mu_yt a \eqn{(d \times 1)} vector of the conditional mean for the time period \eqn{t}.
#' @param B_yt a \eqn{(d \times d)} impact matrix for the time period \eqn{t}.
#' @param e_t a \eqn{(d \times 1)} vector of the structural shocks for the time period \eqn{t}.
#' @details This is used in simulation of the counterfactual scenarios.
#' @return Returns the \eqn{(d \times 1)} vector of observations for the time period \eqn{t}.
#' @keywords internal

get_y_t <- function(mu_yt, B_yt, e_t) {
  # Compute the corresponding observation for the time period t based on the obtained shock vector:
  as.vector(mu_yt + B_yt%*%e_t) # [d, 1] + [d, d] x [d, 1] = [d, 1]
}


#' @title Simulate historical counterfactual for structural STVAR models.
#'
#' @description \code{cfact_hist} simulates historical counterfactual for structural STVAR models.
#'
#' @inheritParams simulate.stvar
#' @inheritParams linear_IRF
#' @param cfact_type a character string indicating the type of counterfactual to be computed: should the path of the policy
#'  variable be fixed to some hypothetical path (\code{cfact_type="fixed_path"}) in given points of time or should the responses
#'  of the policy variable to lagged and contemporaneous movements of some given variable be muted (\code{cfact_type="muted_response"})?
#'  See details for more information.
#' @param policy_var a positive integer between \eqn{1} and \eqn{d} indicating the index of the policy variable considered in the
#'  counterfactual scenario.
#' @param mute_var a positive integer between \eqn{1} and \eqn{d} indicating the index of the variable to whose movements the policy variable
#'  specified in the argument \code{policy_var} should not react to in the counterfactual scenario. This indicates also the index of the shock
#'  to which the policy variable should not react to. It is assumed that \code{mute_var != policy_var}. This argument is only used when
#'  \code{cfact_type="muted_response"}.
#' @param cfact_start a positive integer between \eqn{1} and \eqn{T} indicating the starting period for the counterfactual behavior
#'  of the specified policy variable. It is assumed that the observed time series in indexed as \eqn{y_{t-p+1},...,y_{0},y_1,...,y_T},
#'  i.e., that the first \eqn{p} observations are the initial values, and the "time period one" observation is the \eqn{p+1}th row in
#'  the data matrix.
#' @param cfact_end a positive integer between \code{cfact_start} and \eqn{T} indicating the ending period for the counterfactual
#'  behavior of the specified policy variable.
#' @param cfact_path a numeric vector of length \code{cfact_end-cfact_start+1} indicating the hypothetical path of the policy variable
#'  specified in the argument \code{policy_var}. This argument is only used when \code{cfact_type="fixed_path"}.
#' @details Two types of historical counterfactuals are accommodated where in given historical points of time
#'  either (1) the policy variable of interest takes some hypothetical path (\code{cfact_type="fixed_path"}), or (2)
#'  its responses to lagged and contemporaneous movements of some given variable are shut off (\code{cfact_type="muted_response"}).
#'  In both cases, the counterfactual scenarios are simulated by creating hypothetical shocks to the policy variable of interest
#'  that yield the counterfactual outcome. This approach has the appealing feature that the counterfactual deviations from the
#'  policy reaction function are treated as policy surprises, allowing them to propagate normally, so that the dynamics of the model
#'  are not, per se, tampered but just the policy surprises are.
#'
#'  \strong{Important:} This function assumes that when the policy variable of interest is the \eqn{i_1}th variable, the shock
#'  to it that is manipulated is the \eqn{i_1}th shock. This should be automatically satisfied for recursively identified models,
#'  whereas for model identified by heteroskedasticity or non-Gaussianity, the ordering of the shocks can be generally changed
#'  without loss of generality with the function \code{reorder_B_columns}. In Type (2) counterfactuals it is additionally assumed
#'  that, if the variable to whose movements the policy variable should not react to is the \eqn{i_2}th variable, the shock to it
#'  is the \eqn{i_2}th shock. If it is not clear whether the \eqn{i_2}th shock can be interpreted as a shock to a variable
#'  (but has a broader definition such as "a demand shock"), the Type (2) counterfactual scenario is interpreted as follows: the \eqn{i_1}th
#'  variable does not react to lagged movements of the \eqn{i_2}th variable nor to the \eqn{i_2}th shock.
#'
#'  See the seminal paper of Bernanke et al (1997) for discussing about the "Type (1)" counterfactuals and
#'  Kilian and Lewis (2011) for discussion about the "Type (2)" counterfactuals. See Kilian and Lütkepohl (2017), Section 4.5
#'  for further discussion about the historical counterfactuals. The literature cited about considers linear models, but it is
#'  explained in the vignette of this package how this function computes the historical counterfactuals for the STVAR models in
#'  a way that accommodates nonlinear time-varying dynamics.
#' @return Returns a class \code{'cfacthist'} list with the following elements:
#'   \describe{
#'     \item{cfact_data}{A matrix of size \eqn{(T+p \times d)} containing the counterfactual time series. Note that the first \eqn{p} rows
#'      are for the initial values prior the time period \eqn{t=1}.}
#'     \item{cfact_shocks}{A matrix of size \eqn{(T \times d)} containing the counterfactual shocks.}
#'     \item{cfact_weights}{A matrix of size \eqn{(T \times M)} containing the counterfactual transition weights.}
#'     \item{stvar}{The original STVAR model object.}
#'     \item{input}{A list containing the arguments used to calculate the counterfactual.}
#'  }
#' @seealso \code{\link{GIRF}}, \code{\link{GFEVD}}, \code{\link{linear_IRF}}, \code{\link{hist_decomp}}, \code{\link{cfact_fore}},
#'  \code{\link{cfact_girf}}, \code{\link{fitSSTVAR}}
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
#' # Simulate historical counterfactual where the first variable takes the values 5 and -5
#' # in the first and second time periods, respectively.
#' cfact1 <- cfact_hist(mod32logt, cfact_type="fixed_path", policy_var=1, cfact_start=1,
#'   cfact_end=2, cfact_path=c(5, -5))
#' print(cfact1, start=c(1959, 1), end=c(1960, 4)) # Print cfact data from 1959Q1 to 1960Q4
#' plot(cfact1) # Plot the observed and counterfactual data
#'
#' # Simulate historical counterfactual where the first variable does not respond to lagged
#' # movements of the second variable nor to the second shock in time periods from 10 to 100.
#' cfact2 <- cfact_hist(mod32logt, cfact_type="muted_response", policy_var=1, mute_var=2,
#'  cfact_start=10, cfact_end=100)
#' print(cfact2, start=c(1960, 4), end=c(1963, 4)) # Print cfact data from 1960Q4 to 1963Q4
#' plot(cfact2) # plot the observed and counterfactual data
#' @export

cfact_hist <- function(stvar, cfact_type=c("fixed_path", "muted_response"), policy_var=1, mute_var=NULL, cfact_start=1, cfact_end=1, cfact_path=NULL) {
  check_stvar(stvar, object_name="stvar")
  epsilon <- round(log(.Machine$double.xmin) + 10)
  cfact_type <- match.arg(cfact_type)
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
  if(cfact_type == "fixed_path") {
    if(is.null(cfact_path)) {
      stop("The argument cfact_path needs to be specified")
    } else if(!is.numeric(cfact_path) || length(cfact_path) != cfact_end - cfact_start + 1) {
      stop("The argument cfact_path should be a numeric vector of length cfact_end-cfact_start+1")
    } else if(!is.numeric(cfact_path) || any(is.na(cfact_path))) {
      stop("The argument cfact_path should not contain NA values")
    }
  } else { # cfact_type == "muted_response"
    if(is.null(mute_var)) {
      stop("The argument mute_var needs to be specified")
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
  W <- pick_W(p=p, M=M, d=d, params=params, identification=identification) # NULL for non het.sked ident models
  lambdas <- pick_lambdas(p=p, M=M, d=d, params=params, identification=identification) # NULL for non het.sked ident models

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
  e_t <- get_residuals(data=stvar$data, p=p, M=M, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
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
  # t:th row denotes the vector \bold{y_{i-1}} = (y_{t-1},...,y_{t-p}) (dpx1), assuming the observed data is y_{-p+1},...,y_0,y_1,...,y_{T}.
  Y <- reform_data(stvar$data, p=p) # (T+1 x dp)

  # Create a similar convenient form container for the counterfactual values, including also the original data prior to the countefactual period.
  # First row row initial values vector, and t:th row for (y_{t-1},...,y_{t-p})
  cfact_Y <- matrix(nrow=T_obs + 1, ncol=d*p)
  cfact_Y[1:cfact_start,] <- Y[1:cfact_start,] # Original data, the first row is for the initial value vector

  ## Create a container for the counterfactual observations, including also the original data prior to the counterfactual period:
  cfact_data <- matrix(NA, nrow=nrow(data), ncol=d) # Note that the first p periods are for the initial values
  cfact_data[1:(cfact_start + p - 1),] <- data[1:(cfact_start + p - 1),]

  ## Create a container for the counterfactual transition weights, including also the original ones prior to the counterfactual period:
  cfact_alpha_mt <- matrix(NA, nrow=T_obs, ncol=M) # Note that there are no initial values here, so the indexing is different to cfact_data
  cfact_alpha_mt[1:(cfact_start - 1),] <- alpha_mt[1:(cfact_start - 1),] # Original transition weights

  ## Create a container for the counterfactual shocks, including also the original ones prior to the counterfactual period:
  cfact_e_t <- matrix(NA, nrow=T_obs, ncol=d) # Note that there are no initial values, so the indexing is different to cfact_data
  cfact_e_t[1:(cfact_start - 1),] <- e_t[1:(cfact_start - 1),] # Original transition weights

  ## Simulate the counterfactual observations
  for(t in cfact_start:T_obs) { # Loop through the time periods starting from the beginning of the counterfactual
    t_row_in_data <- t + p # The row in the data matrix (but not in the matrix Y)

    # Calculate the time period t transition weights
    if(M == 1) {
      alpha_mt_t <- c(1)
    } else {
      if(weight_function == "relative_dens") {
        log_mvdvalues <- get_logmvdvalues(Y=cfact_Y, i1=t)
        alpha_mt_t <- get_alpha_mt(M=M, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                   weightpars=weightpars, log_mvdvalues=log_mvdvalues, epsilon=epsilon)
      } else if(weight_function %in% c("logistic", "exponential", "threshold")) {
        alpha_mt_t <- get_alpha_mt(M=M, d=d, Y2=cfact_Y[t, , drop=FALSE], weight_function=weight_function,
                                   weightfun_pars=weightfun_pars, weightpars=weightpars, epsilon=epsilon)
      } else if(weight_function == "mlogit") {
        regression_values <- get_regression_values(Y=cfact_Y, i1=t) # Uses regressions as logmvd values in relative_dens fun
        alpha_mt_t <- get_alpha_mt(M=M, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                   weightpars=rep(1, times=M), log_mvdvalues=regression_values, epsilon=epsilon)
      } else if(weight_function == "exogenous") {
        alpha_mt_t <- alpha_mt[t, , drop=FALSE] # The original transition weights at the time t
      }
    }
    alpha_mt_t <- as.vector(alpha_mt_t)
    cfact_alpha_mt[t, ] <- alpha_mt_t # Store the transition weights

    # Calculate the intercept parameter value for the time period t
    phi_yt <- get_phi_yt(all_phi0=all_phi0, alpha_mt=alpha_mt_t) # [d, 1]

    # Calculate the autoregression matrices A_{y,t,i} for all lags i=1,...,p, for the time period t:
    all_A_yti <- get_allA_yti(all_A=all_A, alpha_mt=alpha_mt_t) # [d, d, p], lag i is obtained from [, , i]

    # Calculate the impact matrix B_{y,t} for the time period t:
    B_yt <- get_B_yt(all_Omegas=all_Omegas, alpha_mt=alpha_mt_t, W=W, lambdas=lambdas, cond_dist=cond_dist,
                     identification=identification) # [d, d]
    if(B_yt[policy_var, policy_var] == 0) {
      if(t %in% cfact_start:cfact_end) {
        stop(paste("The obtained impact matrix B_yt implies that the shock policy_var has zero impact on the policy variable.",
                   "Thus, it is not possible to manipulate shocks to the policy variable to obtain any countefactual scenarios."))
      }
    } else if(abs(B_yt[policy_var, policy_var]) < 1e-6) {
      if(t %in% cfact_start:cfact_end) {
        warning(paste("The shock to the policy variable seems to have a very small effect to the policy variable.",
                      "This can create weird results in the counterfactual scenario."))
      }
    }

    # Compute the conditional mean of the model, is not affected by the period t shocks:
    mu_yt <- get_mu_yt(phi_yt=phi_yt, all_A_yti=all_A_yti, bold_y_t_minus_1=cfact_Y[t,]) # cfact_t t:th row is one lagged

    if(t %in% cfact_start:cfact_end) { # Counterfactual period
      # Compute the counterfactual shock
      if(cfact_type == "fixed_path") {
        e_t_to_use <- e_t[t,] # Recovered shocks for time period t
        cfact_path_t <- cfact_path[t - cfact_start + 1] # The hypothetical path of the policy variable
        effect_of_other_shocks <- as.numeric(crossprod(B_yt[policy_var, -policy_var], e_t_to_use[-policy_var])) # The effect of the other shocks to the policy var
        cfact_policy_e_t <- (cfact_path_t - mu_yt[policy_var] -
                               effect_of_other_shocks)/B_yt[policy_var, policy_var] # The counterfactual shock to the policy var, yielding the cfactual path
        e_t_to_use[policy_var] <- cfact_policy_e_t # Insert the counterfactual shock to the policy variable
        cfact_e_t[t,] <- e_t_to_use # Store the counterfactual shock vector
      } else if(cfact_type == "muted_response") {
        # Lagged effects of mute_var on policy_var:
        lagged_effects <- as.numeric(crossprod(all_A_yti[policy_var, mute_var, ], # (p x 1) vector of i_1i_2:th elements of A_yt so that lag i is in the i:th element.
                                               matrix(cfact_Y[t,], nrow=d)[mute_var,])) # i:th column is the vector y_{t-i}, and the mute_var:th row is the mute_var
        # Contemporaneous effects of mute_var:th shock on policy_var:
        cont_effects <- B_yt[policy_var, mute_var]*e_t[t, mute_var]
        # Calculate the counterfactual shock to the policy variable that mutes the response to mute_var:
        e_t_to_use <- e_t[t,] # Recovered shocks for time period t
        cfact_policy_e_t <- e_t_to_use[policy_var] - (1/B_yt[policy_var, policy_var])*(lagged_effects + cont_effects) # The cfact shock to the policy var
        e_t_to_use[policy_var] <- cfact_policy_e_t # Insert the counterfactual shock to the policy variable
        cfact_e_t[t,] <- e_t_to_use # Store the counterfactual shock vector
      }
    } else { # Outside the counterfactual period
      cfact_e_t[t,] <- e_t[t,] # Use the structural shocks recovered from the fitted model
    }

    # Compute the corresponding observation for the time period t based on the obtained shock vector:
    cfact_data[t_row_in_data, ] <- get_y_t(mu_yt=mu_yt, B_yt=B_yt, e_t=cfact_e_t[t,])

    # Update cfact_Y for the next iteration:
    cfact_Y[t + 1,] <- reform_data(cfact_data[1:t_row_in_data,], p=p)[t + 1,] # The t:th row of cfact_Y is the vector y_{t-1},...,y_{t-p}
  }

  colnames(cfact_data) <- colnames(data) # Set the column names

  # Return the results
  structure(list(cfact_data=ts(cfact_data, start=start(stvar$data), frequency=frequency(stvar$data)),
                 cfact_e_t=cfact_e_t,
                 cfact_alpha_mt=cfact_alpha_mt,
                 stvar=stvar,
                 input=list(policy_var=policy_var, mute_var=mute_var, cfact_type=cfact_type,
                            cfact_start=cfact_start, cfact_end=cfact_end, cfact_path=cfact_path)), class="cfacthist")
}



#' @title Simulate counterfactual forecast scenarios for structural STVAR models.
#'
#' @description \code{cfact_fore} simulates counterfactual forecast scenarios for structural STVAR models.
#'
#' @inheritParams predict.stvar
#' @inheritParams linear_IRF
#' @param nsteps how many steps ahead should be predicted, i.e., the forecast horizon?
#' @param cfact_type a character string indicating the type of counterfactual to be computed: should the path of the policy
#'  variable be fixed to some hypothetical path (\code{cfact_type="fixed_path"}) in given forecast horizons or should the responses
#'  of the policy variable to lagged and contemporaneous movements of some given variable be muted (\code{cfact_type="muted_response"})?
#'  See details for more information.
#' @param policy_var a positive integer between \eqn{1} and \eqn{d} indicating the index of the policy variable considered in the
#'  counterfactual scenario.
#' @param mute_var a positive integer between \eqn{1} and \eqn{d} indicating the index of the variable to whose movements the policy variable
#'  specified in the argument \code{policy_var} should not react to in the counterfactual scenario. This indicates also the index of the shock
#'  to which the policy variable should not react to. It is assumed that \code{mute_var != policy_var}. This argument is only used when
#'  \code{cfact_type="muted_response"}.
#' @param cfact_start a positive integer between \eqn{1} and \code{nsteps} indicating the starting forecast horizon period for the counterfactual
#'  behavior of the specified policy variable.
#' @param cfact_end a positive integer between \code{cfact_start} and \code{nsteps} indicating the ending period for the counterfactual
#'  behavior of the specified policy variable.
#' @param cfact_path a numeric vector of length \code{cfact_end-cfact_start+1} indicating the hypothetical path of the policy variable
#'  specified in the argument \code{policy_var}. This argument is only used when \code{cfact_type="fixed_path"}.
#' @details Two types of counterfactual forecast scenarios are accommodated where in given forecast horizons
#'  either (1) the policy variable of interest takes some hypothetical path (\code{cfact_type="fixed_path"}), or (2)
#'  its responses to lagged and contemporaneous movements of some given variable are shut off (\code{cfact_type="muted_response"}).
#'  In both cases, the counterfactual scenarios are simulated by creating hypothetical shocks to the policy variable of interest
#'  that yield the counterfactual outcome. This approach has the appealing feature that the counterfactual deviations from the
#'  policy reaction function are treated as policy surprises, allowing them to propagate normally, so that the dynamics of the model
#'  are not, per se, tampered but just the policy surprises are.
#'
#'  \strong{Important:} This function assumes that when the policy variable of interest is the \eqn{i_1}th variable, the shock
#'  to it that is manipulated is the \eqn{i_1}th shock. This should be automatically satisfied for recursively identified models,
#'  whereas for model identified by heteroskedasticity or non-Gaussianity, the ordering of the shocks can be generally changed
#'  without loss of generality with the function \code{reorder_B_columns}. In Type (2) counterfactuals it is additionally assumed
#'  that, if the variable to whose movements the policy variable should not react to is the \eqn{i_2}th variable, the shock to it
#'  is the \eqn{i_2}th shock. If it is not clear whether the \eqn{i_2}th shock can be interpreted as a shock to a variable
#'  (but has a broader definition such as "a demand shock"), the Type (2) counterfactual scenario is interpreted as follows: the \eqn{i_1}th
#'  variable does not react to lagged movements of the \eqn{i_2}th variable nor to the \eqn{i_2}th shock.
#'
#'  See the seminal paper of Bernanke et al (1997) for discussing about the "Type (1)" counterfactuals and
#'  Kilian and Lewis (2011) for discussion about the "Type (2)" counterfactuals. See Kilian and Lütkepohl (2017), Section 4.4
#'  for further discussion about counterfactual forecast scenarios. The literature cited about considers linear models, but it is
#'  explained in the vignette of this package how this function computes the counterfactual forecast scenarios for the STVAR models in
#'  a way that accommodates nonlinear time-varying dynamics.
#' @return Returns a class \code{'cfactfore'} list with the following elements:
#'   \describe{
#'     \item{$cfact_pred}{Counterfactual forecast in a class '\code{stvarpred}' object (see \code{?predict.stvar}).}
#'     \item{$pred}{Forecast that does not impose the counterfactual scenario, in a class '\code{stvarpred}' object.}
#'     \item{stvar}{The original STVAR model object.}
#'     \item{input}{A list containing the arguments used to calculate the counterfactual.}
#'  }
#' @seealso \code{\link{GIRF}}, \code{\link{GFEVD}}, \code{\link{linear_IRF}}, \code{\link{hist_decomp}}, \code{\link{cfact_hist}},
#'  \code{\link{cfact_girf}},  \code{\link{fitSSTVAR}}
#' @inherit cfact_hist references
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
#' # Counteractual forecast scenario 5 steps ahead (using only 100 Monte Carlo repetitions
#' # to save computation time), where the first variable takes values 1, -2, and 3 in the
#' # horizons 1, 2, and 3, respectively:
#' set.seed(1)
#' cfact1 <- cfact_fore(mod32logt, nsteps=5, nsim=100, cfact_type="fixed_path", policy_var=1,
#'  cfact_start=1, cfact_end=3, cfact_path=c(1, -2, 3))
#' cfact1 # Print the results
#' plot(cfact1) # Plot the factual and counterfactual forecasts
#'
#' # Counteractual forecast scenario 5 steps ahead (using only 100 Monte Carlo repetitions
#' # to save computation time), where the first variable does not respond to lagged
#' # movements of the second variable nor to the second shock in time periods from 1 to 3:
#' set.seed(1)
#' cfact2 <- cfact_fore(mod32logt, nsteps=5, nsim=100, cfact_type="muted_response", policy_var=1,
#'  mute_var=2, cfact_start=1, cfact_end=3)
#' cfact2 # Print the results
#' plot(cfact2) # Plot the factual and counterfactual forecasts
#' @export

cfact_fore <- function(stvar, nsteps, nsim=1000, pi=0.95, pred_type=c("mean", "median"), exo_weights=NULL,
                       cfact_type=c("fixed_path", "muted_response"), policy_var=1, mute_var=NULL, cfact_start=1, cfact_end=1, cfact_path=NULL) {
  check_stvar(stvar, object_name="stvar")
  cfact_type <- match.arg(cfact_type)
  data <- stvar$data
  d <- ncol(data)

  ## Check the arguments
  if(!is.numeric(policy_var) || length(policy_var) != 1 || policy_var < 1 || policy_var > d || policy_var%%1 != 0) {
    stop("The argument policy_var should be a positive integer between 1 and d")
  } else if(!is.numeric(cfact_start) || length(cfact_start) != 1 || cfact_start < 1 || cfact_start > nsteps || cfact_start%%1 != 0) {
    stop("The argument cfact_start should be a positive integer between 1 and nsteps")
  } else if(!is.numeric(cfact_end) || length(cfact_end) != 1 || cfact_end < cfact_start || cfact_end > nsteps || cfact_end%%1 != 0) {
    stop("The argument cfact_end should be a positive integer between cfact_start and nsteps")
  }
  if(cfact_type == "fixed_path") {
    if(is.null(cfact_path)) {
      stop("The argument cfact_path needs to be specified")
    } else if(!is.numeric(cfact_path) || length(cfact_path) != cfact_end - cfact_start + 1) {
      stop("The argument cfact_path should be a numeric vector of length cfact_end-cfact_start+1")
    } else if(!is.numeric(cfact_path) || any(is.na(cfact_path))) {
      stop("The argument cfact_path should not contain NA values")
    }
  } else { # cfact_type == "muted_response"
    if(is.null(mute_var)) {
      stop("The argument mute_var needs to be specified")
    } else if(!is.numeric(mute_var) || length(mute_var) != 1 || mute_var < 1 || mute_var > d || mute_var%%1 != 0) {
      stop("The argument mute_var should be a positive integer between 1 and d")
    } else if(policy_var == mute_var) {
      stop("The arguments policy_var and mute_var should not be equal")
    }
  }

  ## Compute the counterfactual forecast:
  all_cfact_pred <- predict.stvar(stvar, nsteps=nsteps, nsim=nsim, pi=pi, pred_type=pred_type,
                                  exo_weights=exo_weights, girf_pars=list(cfact_pars=list(cfact_metatype="counterfactual_fore",
                                                                                          cfact_type=cfact_type,
                                                                                          policy_var=policy_var,
                                                                                          mute_var=mute_var,
                                                                                          cfact_start=cfact_start,
                                                                                          cfact_end=cfact_end,
                                                                                          cfact_path=cfact_path)))
  ## Compute the non-countefactual forecast:
  all_pred <- predict.stvar(stvar, nsteps=nsteps, nsim=nsim, pi=pi, pred_type=pred_type, exo_weights=exo_weights)

  # Return the results
  structure(list(cfact_pred=all_cfact_pred,
                 pred=all_pred,
                 stvar=stvar,
                 input=list(policy_var=policy_var,
                            mute_var=mute_var,
                            cfact_type=cfact_type,
                            cfact_start=cfact_start,
                            cfact_end=cfact_end,
                            cfact_path=cfact_path,
                            nsteps=nsteps,
                            nsim=nsim,
                            pi=pi,
                            pred_type=pred_type,
                            exo_weights=exo_weights)), class="cfactfore")
}



#' @title Simulate counterfactual generalized impulse response functions for structural STVAR models.
#'
#' @description \code{cfact_girf} simulates counterfactual generalized impulse response functions for structural STVAR models.
#'
#' @inheritParams GIRF
#' @param cfact_type a character string indicating the type of counterfactual to be computed: should the path of the policy
#'  variable be fixed to some hypothetical path (\code{cfact_type="fixed_path"}) in given impulse response horizons or should the responses
#'  of the policy variable to lagged and contemporaneous movements of some given variable be muted (\code{cfact_type="muted_response"})?
#'  See details for more information.
#' @param policy_var a positive integer between \eqn{1} and \eqn{d} indicating the index of the policy variable considered in the
#'  counterfactual scenario. Note that \code{policy_var} is assumed to satisfy \code{!(policy_var \%in\% which_shocks)}.
#' @param mute_var a positive integer between \eqn{1} and \eqn{d} indicating the index of the variable to whose movements the policy variable
#'  specified in the argument \code{policy_var} should not react to in the counterfactual scenario. This indicates also the index of the shock
#'  to which the policy variable should not react to. It is assumed that \code{mute_var != policy_var}. This argument is only used when
#'  \code{cfact_type="muted_response"}.
#' @param cfact_start a positive integer between \eqn{0} and \code{N} indicating the starting impulse response horizon period for the
#'  counterfactual behavior of the specified policy variable.
#' @param cfact_end a positive integer between \code{cfact_start} and \code{N} indicating the ending period for the counterfactual
#'  behavior of the specified policy variable.
#' @param cfact_path a numeric vector of length \code{cfact_end-cfact_start+1} indicating the hypothetical path of the policy variable
#'  specified in the argument \code{policy_var}. This argument is only used when \code{cfact_type="fixed_path"}.
#' @details Two types of counterfactual generalized impulse response functions (GIRFs) are accommodated where in given impulse response
#'  horizons either (1) the policy variable of interest takes some hypothetical path (\code{cfact_type="fixed_path"}), or (2)
#'  its responses to lagged and contemporaneous movements of some given variable are shut off (\code{cfact_type="muted_response"}).
#'  In both cases, the counterfactual scenarios are simulated by creating hypothetical shocks to the policy variable of interest
#'  that yield the counterfactual outcome. This approach has the appealing feature that the counterfactual deviations from the
#'  policy reaction function are treated as policy surprises, allowing them to propagate normally, so that the dynamics of the model
#'  are not, per se, tampered but just the policy surprises are.
#'
#'  \strong{Important:} This function assumes that when the policy variable of interest is the \eqn{i_1}th variable, the shock
#'  to it that is manipulated is the \eqn{i_1}th shock. This should be automatically satisfied for recursively identified models,
#'  whereas for model identified by heteroskedasticity or non-Gaussianity, the ordering of the shocks can be generally changed
#'  without loss of generality with the function \code{reorder_B_columns}. In Type (2) counterfactuals it is additionally assumed
#'  that, if the variable to whose movements the policy variable should not react to is the \eqn{i_2}th variable, the shock to it
#'  is the \eqn{i_2}th shock. If it is not clear whether the \eqn{i_2}th shock can be interpreted as a shock to a variable
#'  (but has a broader definition such as "a demand shock"), the Type (2) counterfactual scenario is interpreted as follows: the \eqn{i_1}th
#'  variable does not react to lagged movements of the \eqn{i_2}th variable nor to the \eqn{i_2}th shock.
#'
#'  See the seminal paper of Bernanke et al (1997) for discussing about the "Type (1)" counterfactuals and
#'  Kilian and Lewis (2011) for discussion about the "Type (2)" counterfactuals. See Kilian and Lütkepohl (2017), Section 4.5
#'  for further discussion about counterfactuals. The literature cited about considers linear models, but it is
#'  explained in the vignette of this package how this function computes the historical counterfactuals for the STVAR models in
#'  a way that accommodates nonlinear time-varying dynamics.
#' @return Returns a class \code{'cfactgirf'} list with the following elements:
#'   \describe{
#'     \item{\code{$girf}}{An object of class \code{'girf'} containing the counterfactual GIRFs (see \code{?GIRF}).}
#'     \item{\code{$stvar}}{The original STVAR model object.}
#'     \item{\code{$input}}{A list containing the arguments used to calculate the counterfactual.}
#'  }
#' @seealso \code{\link{GIRF}}, \code{\link{GFEVD}}, \code{\link{linear_IRF}}, \code{\link{hist_decomp}}, \code{\link{cfact_hist}},
#'  \code{\link{cfact_fore}},  \code{\link{fitSSTVAR}}
#' @inherit cfact_hist references
#' @examples
#' \donttest{
#' # Recursively identified logistic Student's t STVAR(p=3, M=2) model with the first
#' # lag of the second variable as the switching variable:
#' params32logt <- c(0.5959, 0.0447, 2.6279, 0.2897, 0.2837, 0.0504, -0.2188, 0.4008,
#'   0.3128, 0.0271, -0.1194, 0.1559, -0.0972, 0.0082, -0.1118, 0.2391, 0.164, -0.0363,
#'   -1.073, 0.6759, 3e-04, 0.0069, 0.4271, 0.0533, -0.0498, 0.0355, -0.4686, 0.0812,
#'   0.3368, 0.0035, 0.0325, 1.2289, -0.047, 0.1666, 1.2067, 7.2392, 11.6091)
#' mod32logt <- STVAR(gdpdef, p=3, M=2, params=params32logt, weight_function="logistic",
#'   weightfun_pars=c(2, 1), cond_dist="Student", identification="recursive")
#'
#' # Counterfactual GIRFs for Shock 2 with horizon N=5 (using only R1=R2=10 Monte Carlo repetitions
#' # to save computation time), where the first variable takes values 1, -2, and 3 in the
#' # horizons 1, 2, and 3, respectively:
#' cfact1 <- cfact_girf(mod32logt, which_shocks=2, N=5, R1=10, R2=10, init_regime=1, seeds=1:10,
#'  cfact_type="fixed_path", policy_var=1, cfact_start=1, cfact_end=3, cfact_path=c(1, -2, 3))
#' cfact1 # Print the results
#' plot(cfact1) # Plot the counterfactual GIRF
#'
#' # Counterfactual GIRFs for Shock 2 with horizon N=5 (using only R1=R2=10 Monte Carlo repetitions
#' # to save computation time), where the first variable does not respond to lagged movements
#' # of the second variable nor to the second shock in time periods from 1 to 3:
#' cfact2 <- cfact_girf(mod32logt, which_shocks=2, N=5, R1=10, R2=10, init_regime=1, seeds=1:20,
#'  cfact_type="muted_response", policy_var=1, mute_var=2, cfact_start=1, cfact_end=3)
#' cfact2 # Print the results
#' plot(cfact2) # Plot the counterfactual GIRF
#' }
#' @export

cfact_girf <- function(stvar, which_shocks, shock_size=1, N=30, R1=200, R2=250, init_regime=1, init_values=NULL,
                       which_cumulative=numeric(0), scale=NULL, scale_type=c("instant", "peak"), scale_horizon=N,
                       ci=c(0.95, 0.80), use_data_shocks=FALSE, data_girf_pars=c(0, 0.75, 0, 0, 1.5), ncores=2,
                       burn_in=1000, exo_weights=NULL, seeds=NULL, use_parallel=TRUE,
                       cfact_type=c("fixed_path", "muted_response"), policy_var=1, mute_var=NULL, cfact_start=1, cfact_end=1, cfact_path=NULL) {
  check_stvar(stvar, object_name="stvar")
  cfact_type <- match.arg(cfact_type)
  data <- stvar$data
  d <- ncol(data)

  ## Check the arguments
  if(!is.numeric(policy_var) || length(policy_var) != 1 || policy_var < 1 || policy_var > d || policy_var%%1 != 0) {
    stop("The argument policy_var should be a positive integer between 1 and d")
  } else if(!is.numeric(cfact_start) || length(cfact_start) != 1 || cfact_start < 0 || cfact_start > N || cfact_start%%1 != 0) {
    stop("The argument cfact_start should be a positive integer between 0 and N")
  } else if(!is.numeric(cfact_end) || length(cfact_end) != 1 || cfact_end < cfact_start || cfact_end > N || cfact_end%%1 != 0) {
    stop("The argument cfact_end should be a positive integer between cfact_start and N")
  }
  if(cfact_type == "fixed_path") {
    if(is.null(cfact_path)) {
      stop("The argument cfact_path needs to be specified")
    } else if(!is.numeric(cfact_path) || length(cfact_path) != cfact_end - cfact_start + 1) {
      stop("The argument cfact_path should be a numeric vector of length cfact_end-cfact_start+1")
    } else if(!is.numeric(cfact_path) || any(is.na(cfact_path))) {
      stop("The argument cfact_path should not contain NA values")
    }
  } else { # cfact_type == "muted_response"
    if(is.null(mute_var)) {
      stop("The argument mute_var needs to be specified")
    } else if(!is.numeric(mute_var) || length(mute_var) != 1 || mute_var < 1 || mute_var > d || mute_var%%1 != 0) {
      stop("The argument mute_var should be a positive integer between 1 and d")
    } else if(policy_var == mute_var) {
      stop("The arguments policy_var and mute_var should not be equal")
    }
  }
  if(policy_var %in% which_shocks) {
    stop("The argument policy_var should not be in the argument which_shocks")
  }

  ## Compute the counterfactual GIRF:
  mygirf <- GIRF_int(stvar, which_shocks=which_shocks, shock_size=shock_size, N=N, R1=R1, R2=R2,
                     init_regime=init_regime, init_values=init_values, which_cumulative=which_cumulative,
                     scale=scale, scale_type=scale_type, scale_horizon=scale_horizon,
                     ci=ci, use_data_shocks=use_data_shocks, data_girf_pars=data_girf_pars,
                     ncores=ncores, burn_in=burn_in, exo_weights=exo_weights,
                     seeds=seeds, use_parallel=use_parallel,
                     cfact_pars=list(cfact_metatype="counterfactual_girf",
                                     cfact_type=cfact_type,
                                     policy_var=policy_var,
                                     mute_var=mute_var,
                                     cfact_start=cfact_start,
                                     cfact_end=cfact_end,
                                     cfact_path=cfact_path))

  # Return the results
  structure(list(girf=mygirf,
                 stvar=stvar,
                 input=list(which_shocks=which_shocks,
                            shock_size=shock_size,
                            N=N,
                            R1=R1,
                            R2=R2,
                            init_regime=init_regime,
                            init_values=init_values,
                            which_cumulative=which_cumulative,
                            scale=scale,
                            scale_type=scale_type,
                            scale_horizon=scale_horizon,
                            ci=ci,
                            use_data_shocks=use_data_shocks,
                            data_girf_pars=data_girf_pars,
                            ncores=ncores,
                            burn_in=burn_in,
                            exo_weights=exo_weights,
                            seeds=seeds,
                            use_parallel=use_parallel,
                            policy_var=policy_var,
                            mute_var=mute_var,
                            cfact_type=cfact_type,
                            cfact_start=cfact_start,
                            cfact_end=cfact_end,
                            cfact_path=cfact_path)), class="cfactgirf")
}
