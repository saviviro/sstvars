
#' @title Log-likelihood function
#'
#' @description \code{loglikelihood} log-likelihood function of a smooth transition VAR model
#'
#' @param data data a matrix or class \code{'ts'} object with \code{d>1} columns. Each column is taken to represent
#'  a univariate time series. Missing values are not supported.
#' @param p a positive integer specifying the autoregressive order
#' @param M a positive integer specifying the number of regimes
#' @param params a real valued vector specifying the parameter values.
#'   Should have the form \eqn{\theta = (\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M,\sigma,\alpha,\nu)},
#'   where
#'   \itemize{
#'     \item{\eqn{\phi_{m,0} = } the \eqn{(d \times 1)} intercept (or mean) vector of the \eqn{m}th regime.}
#'     \item{\eqn{\varphi_m = (vec(A_{m,1}),...,vec(A_{m,p}))} \eqn{(pd^2 \times 1)}.}
#'     \item{\eqn{\sigma = (vech(\Omega_1),...,vech(\Omega_M)} \eqn{(Md(d + 1)/2 \times 1)}.}
#'     \item{\eqn{\alpha = } the \eqn{(a\times 1)} vector containing the transition weight parameters.}
#'     \item{\eqn{\nu > 2} is the degrees of freedom parameter that is included only if \code{cond_dist="Student"}.}
#'   }
#'   For models with...
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
#'     \item{AR_constraints:}{Replace \eqn{\varphi_1,...,\varphi_M} with \eqn{\psi} as described in the argument \code{AR_constraints}.}
#'     \item{mean_constraints:}{Replace \eqn{\phi_{1,0},...,\phi_{M,0}} with \eqn{(\mu_{1},...,\mu_{g})} where
#'           \eqn{\mu_i, \ (d\times 1)} is the mean parameter for group \eqn{i} and \eqn{g} is the number of groups.}
#'     \item{weight_constraints:}{Replace \eqn{\alpha} with \eqn{\xi} as described in the argument \code{weigh_constraints}.}
#'   }
#'   Above, \eqn{\phi_{m,0}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}th coefficient matrix of the \eqn{m}th
#'   mixture component, and \eqn{\Omega_{m}} denotes the error term covariance matrix of the \eqn{m}:th mixture component.
#'   If \code{parametrization=="mean"}, just replace each \eqn{\phi_{m,0}} with regimewise mean \eqn{\mu_{m}}.
#'   \eqn{vec()} is vectorization operator that stacks columns of a given matrix into a vector. \eqn{vech()} stacks columns
#'   of a given matrix from the principal diagonal downwards (including elements on the diagonal) into a vector.
#' @param weight_function What type of transition weights \eqn{\alpha_{m,t}} should be used?
#'  \describe{
#'    \item{\code{"relative_dens"}:}{\eqn{\alpha_{m,t}=
#'      \frac{\alpha_mf_{m,dp}(y_{t-1},...,y_{t-p+1})}{\sum_{n=1}^M\alpha_nf_{n,dp}(y_{t-1},...,y_{t-p+1})}}, where
#'      \eqn{\alpha_m\in (0,1)} are weight parameters that satisfy \eqn{\sum_{m=1}^M\alpha_m=1} and
#'      \eqn{f_{m,dp}(\cdot)} is the \eqn{dp}-dimensional stationary density of the \eqn{m}th regime corresponding to \eqn{p}
#'      consecutive observations. Available for Gaussian conditional distribution only.}
#'    \item{\code{"logistic"}:}{\eqn{M=2}, \eqn{\alpha_{1,t}=1-\alpha_{2,t}},
#'      and \eqn{\alpha_{2,t}=[1+\exp\lbrace -\gamma(y_{it-j}-c) \rbrace]^{-1}}, where \eqn{y_{it-j}} is the lag \eqn{j}
#'      observation of the \eqn{i}th variable, \eqn{c} is a location parameter, and \eqn{\gamma > 0} is a scale parameter.}
#'    \item{\code{"mlogit"}:}{\eqn{\alpha_{m,t}=\frac{\exp\lbrace \gamma_m'z_{t-1} \rbrace}
#'      {\sum_{n=1}^M\exp\lbrace \gamma_n'z_{t-1} \rbrace}}, where \eqn{\gamma_m} are coefficient vectors, \eqn{\gamma_M=0},
#'      and \eqn{z_{t-1}} \eqn{(k\times 1)} is the \eqn{\mathcal{F}_{t-1}}-measurable vector containing a constant and
#'      the (lagged) switching variables.}
#'    \item{\code{"exponential"}:}{\eqn{M=2}, \eqn{\alpha_{1,t}=1-\alpha_{2,t}},
#'      and \eqn{\alpha_{2,t}=1-\exp\lbrace -\gamma(y_{it-j}-c) \rbrace}, where \eqn{y_{it-j}} is the lag \eqn{j}
#'      observation of the \eqn{i}th variable, \eqn{c} is a location parameter, and \eqn{\gamma > 0} is a scale parameter.}
#'    \item{\code{"threshold"}:}{\eqn{\alpha_{m,t} = 1} if \eqn{r_{m-1}<y_{it-j}\leq r_{m}} and \eqn{0} otherwise, where
#'       \eqn{-\infty\equiv r_0<r_1<\cdots <r_{M-1}<r_M\equiv\infty} are thresholds \eqn{y_{it-j}} is the lag \eqn{j}
#'       observation of the \eqn{i}th variable.}
#'  }
#'  See the vignette for more details about the weight functions.
#' @param weightfun_pars \describe{
#'   \item{If \code{weight_function == "relative_dens"}:}{Not used.}
#'   \item{If \code{weight_function \%in\% c("logistic", "exponential", "threshold")}:}{a numeric vector with the switching variable
#'     \eqn{i\in\lbrace 1,...,d \rbrace} in the first and the lag \eqn{j\in\lbrace 1,...,p \rbrace} in the second element.}
#'   \item{If \code{weight_function == "mlogit"}:}{a list of two elements:
#'     \describe{
#'       \item{The first element \code{$vars}:}{a numeric vector containing the variables that should used as switching variables
#'         in the weight function in an increasing order, i.e., a vector with unique elements in \eqn{\lbrace 1,...,d \rbrace}.}
#'       \item{The second element \code{$lags}:}{an integer in \eqn{\lbrace 1,...,p \rbrace} specifying the number of lags to be
#'         used in the weight function.}
#'     }
#'   }
#' }
#' @param cond_dist specifies the conditional distribution of the model as \code{"Gaussian"} or \code{"Student"}.
#' @param parametrization \code{"intercept"} or \code{"mean"} determining whether the model is parametrized with intercept
#'   parameters \eqn{\phi_{m,0}} or regime means \eqn{\mu_{m}}, m=1,...,M.
#' @param identification is it reduced form model or an identified structural model; if the latter, how is it identified
#'   (see the vignette or the references for details)?
#' @param AR_constraints a size \eqn{(Mpd^2 x q)} constraint matrix \eqn{C} specifying linear constraints
#'   to the autoregressive parameters. The constraints are of the form
#'   \eqn{(\varphi_{1},...,\varphi_{M}) = C\psi}, where \eqn{\varphi_{m} = (vec(A_{m,1}),...,vec(A_{m,p})) \ (pd^2 x 1),\ m=1,...,M},
#'   contains the coefficient matrices and \eqn{\psi} \eqn{(q x 1)} contains the related parameters.
#'   For example, to restrict the AR-parameters to be the identical across the regimes, set \eqn{C =}
#'   [\code{I:...:I}]' \eqn{(Mpd^2 x pd^2)} where \code{I = diag(p*d^2)}.
#' @param mean_constraints Restrict the mean parameters of some regimes to be identical? Provide a list of numeric vectors
#'   such that each numeric vector contains the regimes that should share the common mean parameters. For instance, if
#'   \code{M=3}, the argument \code{list(1, 2:3)} restricts the mean parameters of the second and third regime to be
#'   identical but the first regime has freely estimated (unconditional) mean. Ignore or set to \code{NULL} if mean parameters
#'   should not be restricted to be the same among any regimes. \strong{This constraint is available only for mean parametrized models;
#'   that is, when \code{parametrization="mean"}.}
#' @param weight_constraints a list of two elements, \eqn{R} in the first element and \eqn{r} in the second element,
#'   specifying linear constraints on the transition weight parameters \eqn{\alpha}.
#'   The constraints are of the form \eqn{\alpha = R\xi + r}, where \eqn{R} is a known \eqn{(a\times l)}
#'   constraint matrix of full column rank (\eqn{a} is the dimension of \eqn{\alpha}), \eqn{r} is a known \eqn{(a\times 1)} constant,
#'   and \eqn{\xi} is an unknown \eqn{(l\times 1)} parameter. \strong{Alternatively}, set \eqn{R=0} in order to constrain the
#'   the weight parameter to the constant \eqn{r} (in this case, \eqn{\alpha} is dropped from the constrained parameter vector).
#' @param B_constraints NOT YET IMPLEMENTED!
#' @param to_return should the returned object be the log-likelihood, which is the default, or something else?
#'   See the section "Return" for all the options.
#' @param check_params should it be checked that the parameter vector satisfies the model assumptions? Can be skipped to save
#'   computation time if it does for sure.
#' @param minval the value that will be returned if the parameter vector does not lie in the parameter space
#'   (excluding the identification condition).
#' @param stab_tol numerical tolerance for stability of condition of the regimes: if the "bold A" matrix of any regime
#'   has eigenvalues larger that \code{1 - stat_tol} the parameter is considered to be outside the parameter space.
#'   Note that if tolerance is too small, numerical evaluation of the log-likelihood might fail and cause error.
#' @param posdef_tol numerical tolerance for positive definiteness of the error term covariance matrices: if
#'   the error term covariance matrix of any regime has eigenvalues smaller than this, the parameter is considered
#'   to be outside the parameter space. Note that if the tolerance is too small, numerical evaluation of the
#'   log-likelihood might fail and cause error.
#' @param distpar_tol the parameter vector is considered to be outside the parameter space if the degrees of
#'   freedom parameters is not larger than \code{2 + distpar_tol} (applies only if \code{cond_dist="Student"}).
#' @param weightpar_tol numerical tolerance for weight parameters being in the parameter space. Values closer to
#'   to the border of the parameter space than this are considered to be "outside" the parameter space.
#' @details Calculates the log-likelihood of the specified model.
#' @return
#'   \describe{
#'     \item{If \code{to_return="loglik"}:}{the log-likelihood of the specified model.}
#'     \item{If \code{to_return=="tw"}:}{a size \eqn{((n_obs-p)\times M)} matrix containing the transition weights: for m:th component
#'       in m:th column.}
#'     \item{If \code{to_return=="loglik_and_tw"}:}{a list of two elements. The first element (\code{$loglik}) contains the
#'       log-likelihood and the second element (\code{$tw}) contains the transition weights.}
#'     \item{If \code{to_return=="terms"}:}{a size \eqn{((n_obs-p)\times 1)} numeric vector containing the terms \eqn{l_{t}}.}
#'     \item{If \code{to_return=="regime_cmeans"}:}{an \code{[n_obs-p, d, M]} array containing the regimewise conditional means.}
#'     \item{If \code{to_return=="total_cmeans"}:}{a \code{[n_obs-p, d]} matrix containing the conditional means of the process.}
#'     \item{If \code{to_return=="total_ccovs"}:}{an \code{[d, d, n_obs-p]} array containing the conditional covariance matrices of
#'       the process.}
#'   }
#' @references
#'  \itemize{
#'    \item Lütkepohl H. 2005. New Introduction to Multiple Time Series Analysis,
#'          \emph{Springer}.
#'    \item McElroy T. 2017. Computation of vector ARMA autocovariances.
#'          \emph{Statistics and Probability Letters}, \strong{124}, 92-96.
#'    \item References to the STVAR models TO BE INCLUDED.
#'  }
#' @keywords internal

loglikelihood <- function(data, p, M, params, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold"),
                          weightfun_pars=NULL, cond_dist=c("Gaussian", "Student"), parametrization=c("intercept", "mean"),
                          identification=c("reduced_form", "impact_responses", "heteroskedasticity", "other"),
                          AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL,
                          to_return=c("loglik", "tw", "loglik_and_tw", "terms", "regime_cmeans", "total_cmeans", "total_ccovs"),
                          check_params=TRUE, minval=NULL, stab_tol=1e-3, posdef_tol=1e-8, distpar_tol=1e-8, weightpar_tol=1e-8) {
  # Match args
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  identification <- match.arg(identification)
  to_return <- match.arg(to_return)
  if(identification != "reduced_form") stop("Structural models are not yet implemented to loglikelihood")
  if(!is.null(B_constraints)) stop("B_constrained models are not yet implemented to loglikelihood")

  # Compute some required statistics
  epsilon <- round(log(.Machine$double.xmin) + 10) # Logarithm of the smallest value that can be handled normally
  d <- ncol(data)
  n_obs <- nrow(data)
  T_obs <- n_obs - p

  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars)


  # Collect the parameter values
  # First remove all constraints, if any; also switch to reduced form parameter vector;
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, weight_function=weight_function,
                                    cond_dist=cond_dist, identification=identification,
                                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                    weight_constraints=weight_constraints, B_constraints=B_constraints,
                                    weightfun_pars=weightfun_pars)

  # Pick params
  if(parametrization == "intercept") { # [d, M]
    all_phi0 <- pick_phi0(M=M, d=d, params=params)
  } else {
    all_mu <- pick_phi0(M=M, d=d, params=params) # mean parameters instead of intercepts
  }
  all_A <- pick_allA(p=p, M=M, d=d, params=params) # [d, d, p, M]
  all_Omegas <- pick_Omegas(p=p, M=M, d=d, params=params) # [d, d, M]
  weightpars <- pick_weightpars(p=p, M=M, d=d, params=params, weight_function=weight_function,
                                cond_dist=cond_dist, weightfun_pars=weightfun_pars)
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  distpars <- pick_distpars(params=params, cond_dist=cond_dist)

  # Check that the parameter vector lies in the parameter space
  if(check_params) {
    if(!in_paramspace(p=p, M=M, d=d, weight_function=weight_function, cond_dist=cond_dist,
                      all_boldA=all_boldA, all_Omegas=all_Omegas, weightpars=weightpars,
                      distpars=distpars, weightfun_pars=weightfun_pars,
                      stab_tol=stab_tol, posdef_tol=posdef_tol,
                      distpar_tol=distpar_tol, weightpar_tol=weightpar_tol)) {
      return(minval)
    }
  }

  # i:th row denotes the vector \bold{y_{i-1}} = (y_{i-1},...,y_{i-p}) (dpx1),
  # assuming the observed data is y_{-p+1},...,y_0,y_1,...,y_{T}
  Y <- reform_data(data, p)
  Y2 <- Y[1:T_obs,] # Last row removed; not needed when calculating something based on lagged observations

  # Calculate unconditional regime-specific expected values (column per component) or phi0-parameters if using mean-parametrization
  Id <- diag(nrow=d)
  if(parametrization == "intercept") {
    all_mu <- vapply(1:M, function(m) solve(Id - rowSums(all_A[, , , m, drop=FALSE], dims=2), all_phi0[,m]), numeric(d)) # sum over dims+1=3
  } else {
    all_phi0 <- vapply(1:M, function(m) (Id - rowSums(all_A[, , , m, drop=FALSE], dims=2))%*%all_mu[,m], numeric(d))
  }

  # Calculate the transition weights [T_obs, M] with [t,m] indexing (nothing for the initial values here)
  alpha_mt <- get_alpha_mt(data=data, Y2=Y2, p=p, M=M, d=d, weight_function=weight_function, all_A=all_A, all_boldA=all_boldA,
                           all_Omegas=all_Omegas, weightpars=weightpars, weightfun_pars=weightfun_pars, all_mu=all_mu, epsilon=epsilon)
  if(to_return == "tw") {
    return(alpha_mt)
  }

  # Calculate the conditional means mu_{m,t}
  # The dimensions of mu_mt will be: [t, p, m]
  all_A2 <- array(all_A, dim=c(d, d*p, M)) # cbind coefficient matrices of each component: m:th component is obtained at [, , m]

  mu_yt <- get_mu_yt_Cpp(obs=Y2, all_phi0=all_phi0, all_A=all_A2, alpha_mt=alpha_mt)

  # R implementation saved below for speed comparisons
  #mu_mt <- array(vapply(1:M, function(m) t(all_phi0[, m] + tcrossprod(all_A2[, , m], Y2)), numeric(d*T_obs)), dim=c(T_obs, d, M)) # [, , m]
  #mu_yt <- vapply(1:d, function(i1) rowSums(alpha_mt*mu_mt[,i1,]), numeric(T_obs)) # [T_obs, d]

  # Return conditional moments if those were to be returned (R implementation used, as computation speed is no issue here)
  if(to_return == "regime_cmeans") { # Regime-specific conditional menas
    return(array(vapply(1:M, function(m) t(all_phi0[, m] + tcrossprod(all_A2[, , m], Y2)), numeric(d*T_obs)), dim=c(T_obs, d, M))) # [, , m]
  } else if(to_return == "total_cmeans") { # Cond means of the process: weighted sum of regime-specific conditional means
    #mu_mt <- array(vapply(1:M, function(m) t(all_phi0[, m] + tcrossprod(all_A2[, , m], Y2)), numeric(d*T_obs)), dim=c(T_obs, d, M)) # [, , m]
    #return(matrix(rowSums(vapply(1:M, function(m) alpha_mt[,m]*mu_mt[, , m], numeric(d*T_obs))), nrow=T_obs, ncol=d, byrow=FALSE))
    return(mu_yt)
  } else if(to_return == "total_ccovs") { # Cond covariance matrices of the process: weighted sum of regime-specific cond cov mats
    all_covmats <- array(rowSums(vapply(1:M, function(m) rep(alpha_mt[, m], each=d*d)*as.vector(all_Omegas[, , m]),
                                        numeric(d*d*T_obs))), dim=c(d, d, T_obs))
    return(all_covmats)
  }

  # Calculate the conditional log-likelihood; the initial values are not used here
  if(cond_dist == "Gaussian") { # Gaussian conditional distribution
    all_lt <- -0.5*d*log(2*pi) + Gaussian_densities_Cpp(obs=data[(p+1):nrow(data),], means=mu_yt, covmats=all_Omegas, alpha_mt=alpha_mt)

    # BELOW IS THE R IMPLEMENTATION FOR SPEED COMPARISONS
    #all_covmats <- array(rowSums(vapply(1:M, function(m) rep(alpha_mt[, m], each=d*d)*as.vector(all_Omegas[, , m]),
    # numeric(d*d*T_obs))), dim=c(d, d, T_obs))
    #obs_minus_cmean <- data[(p+1):nrow(data),] - mu_yt
    #all_lt <- numeric(T_obs)
    #tmp0 <- -0.5*d*log(2*pi)
    # for(i1 in 1:T_obs) {
    #    # Calculate the l_t log multinormal density for each observation
    #   cond_covmat <- matrix(0, nrow=d, ncol=d)
    #   for(i2 in 1:M) {
    #     cond_covmats[, , i2] <- alpha_mt[i1, i2]*all_Omegas[, , i2]
    #   }
    #   cond_covmat <- apply(cond_covmats, MARGIN=1:2, sum)
    #   all_lt[i1] <- tmp0 - 0.5*log(det(cond_covmat)) - 0.5*crossprod(obs_minus_cmean[i1,],
    #                                                                          chol2inv(chol(cond_covmat))%*%(obs_minus_cmean[i1,]))
    # }
  } else if(cond_dist == "Student") {
    logCd <- lgamma(0.5*(d + distpars)) - 0.5*d*log(base::pi) - 0.5*d*log(distpars - 2) - lgamma(0.5*distpars)
    all_lt <- logCd + Student_densities_Cpp(obs=data[(p+1):nrow(data),], means=mu_yt, covmats=all_Omegas, alpha_mt=alpha_mt, df=distpars)

    # # BELOW AN R IMPLEMENTATION FOR SPEED COMPARISONS
    # all_covmats <- array(rowSums(vapply(1:M, function(m) rep(alpha_mt[, m], each=d*d)*as.vector(all_Omegas[, , m]),
    #                                     numeric(d*d*T_obs))), dim=c(d, d, T_obs))
    # obs_minus_cmean <- data[(p+1):nrow(data),] - mu_yt
    # all_lt <- numeric(T_obs)
    # for(i1 in 1:T_obs) {
    #    # Calculate the l_t log multistudent density for each observation
    #   all_lt[i1] <- logCd - 0.5*log(det(all_covmats[, , i1])) -
    #     0.5*(d + distpars)*log(1 + crossprod(obs_minus_cmean[i1,], chol2inv(chol(all_covmats[, , i1]))%*%(obs_minus_cmean[i1,]))/(distpars - 2))
    # }
  }
  if(to_return == "terms") {
    return(all_lt)
  } else if(to_return == "loglik_and_tw") {
    return(list(loglik=sum(all_lt),
                tw=alpha_mt))
  }
  sum(all_lt)
}




#' @title Get the transition weights alpha_mt
#'
#' @description \code{get_alpha_mt} computes the transition weights.
#'
#' @inheritParams loglikelihood
#' @inheritParams in_paramspace
#' @inheritParams get_Sigmas
#' @param Y2 the data arranged as obtained from \code{reform_data(data, p)} but excluding the last row
#' @param all_mu an \eqn{(d \times M)} matrix containing the unconditional regime-specific means
#' @param epsilon the smallest number such that its exponent is wont classified as numerically zero
#'   (around \code{-698} is used).
#' @param log_mvdvalues a \eqn{T x M} matrix containing log multivariate normal densities (can be used with
#'   relative dens weight function only)
#' @details Note that we index the time series as \eqn{-p+1,...,0,1,...,T}.
#' @return Returns the mixing weights a \eqn{(T x M)} matrix, so that the t:th row is for the time point t
#'   and m:th column is for the regime m.
#' @inherit in_paramspace references
#' @keywords internal

get_alpha_mt <- function(data, Y2, p, M, d, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold"),
                         weightfun_pars=NULL, all_A, all_boldA, all_Omegas, weightpars, all_mu, epsilon, log_mvdvalues=NULL) {
  weight_function <- match.arg(weight_function)
  if(is.null(log_mvdvalues)) {
    T_obs <- ifelse(missing(data), 1, nrow(data) - p) # simulate.stvar uses without data, needs to return 1 if M=1.
    if(M == 1) {
      return(as.matrix(rep(1, times=T_obs)))
    }
  } else {
    if(M == 1) {
      if(!is.matrix(log_mvdvalues)) { # Only one observation
        return(as.matrix(1)) # Only one observation and one regime
      } else {
        return(as.matrix(rep(1, times=nrow(log_mvdvalues)))) # Multiple observations but only one regime
      }
    }
  }
  if(weight_function == "logistic") {
    # According to "lag" in weightfun_pars[2], only the column d(lag-1) + variable are used where "variable i is the switching variable"
    subY2 <- Y2[,d*(weightfun_pars[2] - 1) + weightfun_pars[1]] # Returns a vector
    in_exp <- -weightpars[2]*(subY2 - weightpars[1])
    in_exp[in_exp > -epsilon] <- -epsilon # Values larger than that would produce Inf and "break" the loglikelihood function; epsilon -698
    alpha_2t <- (1 + exp(in_exp))^(-1) # Weights of the second regime
    return(unname(cbind(1 - alpha_2t, alpha_2t)))
  } else if(weight_function == "exponential") {
    # According to "lag" in weightfun_pars[2], only the column d(lag-1) + variable are used where "variable i is the switching variable"
    subY2 <- Y2[,d*(weightfun_pars[2] - 1) + weightfun_pars[1]] # Returns a vector
    in_exp <- -weightpars[2]*(subY2 - weightpars[1])^2
    in_exp[in_exp > -epsilon] <- -epsilon # Values larger than that would produce Inf and "break" the loglikelihood function; epsilon -698
    alpha_2t <- 1 - exp(in_exp) # Weights of the second regime
    return(unname(cbind(1 - alpha_2t, alpha_2t)))
  } else if(weight_function == "threshold") {
    # According to "lag" in weightfun_pars[2], only the column d(lag-1) + variable are used where "variable i is the switching variable"
    subY2 <- Y2[,d*(weightfun_pars[2] - 1) + weightfun_pars[1]] # Returns a vector
    alpha_mt <- matrix(0, nrow=length(subY2), ncol=M) # [t, m] fill with ones according to the thresholds; weightpars in increasing order
    for(m in 1:M) {
      if(m == 1) {
        alpha_mt[subY2 <= weightpars[m], m] <- 1
      } else if(m < M) { # m > 1 && m < M; m>1 known if we end up here
        alpha_mt[subY2 > weightpars[m-1] & subY2 <= weightpars[m], m] <- 1
      } else { # m == M; m < M false if we end up here so we know this
        alpha_mt[subY2 > weightpars[m-1], m] <- 1
      }
    }
    return(alpha_mt) # Each row has one 1 and the rest are zero
  }

  if(weight_function == "mlogit" && is.null(log_mvdvalues)) {

    # M-1 vectors gamma_m, since gamma_M = 0.
    all_gamma_m <- matrix(weightpars, ncol=M-1) # Column per gamma_m, m=1,...,M-1
    #all_gamma_m <- cbind(all_gamma_m, 0) # Add gamma_M = 0 as the M:th column

    # To get the regressor matrix, we need matrix such that i:th row=(1,z_{min(J)},...,z_{max(J)})

    # Subset columns of Y2 according to vars and lags in weightfun_pars
    vars <- weightfun_pars[[1]]
    lags <- weightfun_pars[[2]]

    # According to "lags" in weightfun_pars[[2]], only the columns 1,...,d*lags are used
    subY2 <- Y2[,1:(d*lags)]

    # Then, take only the columns related to the chosen switching-variables
    # y_{i-j}, j=1,...,lags, has the length d, and each y_{i-j} only uses the cols in weightfun_pars[[1]]

    #   vars <- 2
    #    1:(d*weightfun_pars[[2]])
    # in each repetition in 1...lags, we take the elements given by vars
    # For instance, if vars=2, lags=2 and d=2, we take the elements 2,4

    # If vars=1,3, lags=2, d=3, we take the elements 1,3,4,6
    #    vars <- c(1,3); lags <- 2; d <- 3
    # Define the starting indices - 1 for each lag in 1,...,lags, add vars to obtain the indices.
    # e.g., starts=0,3, then 0 + c(1,3) = 1,3 and 3 + c(1,3) = 4,6 --> 1,3,4,6
    # Use matrices to calculate all additions without loops at the same time

    lowers <- (1:lags - 1)*d # We want add vars to each of these
    tmp <- matrix(lowers, nrow=length(vars), ncol=length(lowers), byrow=TRUE) # rep lowers as the rows
    subY2 <- subY2[,as.vector(tmp + vars)]  # add vars to each column, and obtain the columns to subset
    all_z_tilde <- cbind(1, subY2)  # i:th row=(1,z_{min(J)},...,z_{max(J)})

    # Calculate the regressions gamma_m'z_{t-1} based on all_gamma_m and all_z_tilde
    regressions_mt <- matrix(0, nrow=nrow(Y2), ncol=M) # The last column is a column of zeros
    for(i1 in 1:ncol(all_gamma_m)) { # i1=1,...,M-1
      #for(i2 in 1:nrow(Y2)) {
      #  regressions_mt[i2, i1] <- crossprod(all_gamma_m[,i1], all_z_tilde[i2,])
      #}
      regressions_mt[,i1] <- t(all_gamma_m[,i1])%*%t(all_z_tilde)
    }

    # Note that if abs(regressions_mt) > epsilon, there will be infs when taking exponent. Therefore, similar
    # procedure to the relative_dens weights that handle too large values correctly but computationally fast is required.
    # We can directly make use of the code for relative_dens weights, because the rest of the the calculations are identical
    # to relative_dens weight functions with log_mvdvalues <- regressions_mt; weightpars <- rep(1, times=M)
    # i.e., we regressions instead of log_mvdvalues, and we do not weight the exponents of the regressions.
    log_mvdvalues <- regressions_mt
    weightpars <- rep(1, times=M) # Overwrites weightpars; the original ones are not needed here anymore
  }

  if(weight_function == "relative_dens" || weight_function == "mlogit") { # mlogit defines log_mvdvalues and weightpars above
    if(is.null(log_mvdvalues)) {
      # Calculate the covariance matrices Sigma_{m,p} (Lutkepohl 2005, eq. (2.1.39) or the algorithm proposed by McElroy 2017)
      Sigmas <- get_Sigmas(p=p, M=M, d=d, all_A=all_A, all_boldA=all_boldA, all_Omegas=all_Omegas) # Store the (dpxdp) covariance matrices

      #obs_minus_mean <- Y2 - rep(all_mu[,m], times=p)
      #const_term <- -0.5*d*log(2*pi) - 0.5*log(det(Sigmas))
      #inv_cholcovmat = solve(chol(Sigmas[, , m]))
      #log_mvdvalues <- vapply(1:M, function(m) -0.5*d*log(2*pi)
      # - 0.5*log(det(Sigmas[, , m])) - 0.5*solve(chol(Sigmas[, , m])), numeric(T_obs))
      # log_mvdvalues_test0 <- matrix(nrow=T_obs, ncol=M)
      # for(m in 1:M) {
      #   obs_minus_mean <- t(t(Y2) - rep(all_mu[,m], times=p))
      #   for(i1 in 1:T_obs) {
      #     log_mvdvalues_test0[i1, m] <-  -0.5*d*log(2*pi) - 0.5*log(det(Sigmas[, , m])) -
      #       0.5*t(obs_minus_mean[i1,])%*%solve(Sigmas[, , m])%*%obs_minus_mean[i1,]
      #   }
      # }
      # Calculate the dp-dimensional multinormal densities in logarithm with the package mvnfast:
      # i:th row for index i-1 etc, m:th column for m:th component.
      # We calculate in logarithm because the non-log values may be too close to zero for machine accuracy (if they are too close to zero
      # for all regimes and computer handles them as zero, we would divide by zero when calculating the transition weights).
      # Cholesky decomposition is taken in R to avoid unnecessary warnings that caused by numerical error
      # that makes the matrices to be not exactly symmetric but only up to numerical tolerance.
      log_mvdvalues <- -0.5*d*log(2*pi) + vapply(1:M, function(m) Gaussian_densities_const_Cpp(obs=Y2,
                                                                                               mean=matrix(rep(all_mu[,m], times=p), nrow=1),
                                                                                               cholcovmat=chol(Sigmas[, , m])),
                                                 numeric(T_obs)) # [T_obs, M] removes the period T+1 weights
    }

    # Calculate the transition weights based on the log-multivariate density values
    if(!is.matrix(log_mvdvalues)) log_mvdvalues <- t(as.matrix(log_mvdvalues)) # Only one time point but multiple regimes
    log_mvdvalues_orig <- log_mvdvalues
    small_logmvns <- log_mvdvalues < epsilon
    if(any(small_logmvns)) {
      # If too small or large non-log-density values are present (i.e., that would yield -Inf or Inf),
      # we replace them with ones that are not too small or large but imply the same mixing weights
      # up to negligible numerical tolerance (tested in gmvarkit).
      which_change <- rowSums(small_logmvns) > 0 # Which rows contain too small values
      to_change <- log_mvdvalues[which_change, , drop=FALSE]
      largest_vals <- do.call(pmax, split(to_change, f=rep(1:ncol(to_change), each=nrow(to_change)))) # The largest values of those rows
      diff_to_largest <- to_change - largest_vals # Differences to the largest value of the row

      # For each element in each row, check the (negative) distance from the largest value of the row. If the difference
      # is smaller than epsilon, replace the with epsilon. The results are then the new log_mvn values.
      diff_to_largest[diff_to_largest < epsilon] <- epsilon

      # Replace the old log_mvdvalues with the new ones
      log_mvdvalues[which_change,] <- diff_to_largest
    }

    mvnvalues <- exp(log_mvdvalues)
    denominator <- as.vector(mvnvalues%*%weightpars)
    alpha_mt <- (mvnvalues/denominator)%*%diag(weightpars)
  } else {
    stop("Only relative_dens and mlogit weight functions are currently implemented to get_alpha_mt")
  }
  alpha_mt
}
