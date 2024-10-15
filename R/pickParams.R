#' @title Pick \eqn{\phi_{m,0}} or \eqn{\mu_{m}}, m=1,..,M vectors
#'
#' @description \code{pick_phi0} picks the intercept or mean parameters from the given parameter vector.
#'
#' @param M the number of regimes
#' @param d the number of time series in the system, i.e., the dimension
#' @param params a real valued vector specifying the parameter values.
#'   Should have the form \eqn{\theta = (\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M,\sigma,\alpha,\nu)},
#'   where (see exceptions below):
#'   \itemize{
#'     \item{\eqn{\phi_{m,0} = } the \eqn{(d \times 1)} intercept (or mean) vector of the \eqn{m}th regime.}
#'     \item{\eqn{\varphi_m = (vec(A_{m,1}),...,vec(A_{m,p}))} \eqn{(pd^2 \times 1)}.}
#'     \item{\describe{
#'       \item{if \code{cond_dist="Gaussian"} or \code{"Student"}:}{\eqn{\sigma = (vech(\Omega_1),...,vech(\Omega_M))}
#'         \eqn{(Md(d + 1)/2 \times 1)}.}
#'       \item{if \code{cond_dist="ind_Student"} or \code{"ind_skewed_t"}:}{\eqn{\sigma = (vec(B_1),...,vec(B_M)} \eqn{(Md^2 \times 1)}.}
#'       }
#'     }
#'     \item{\eqn{\alpha = } the \eqn{(a\times 1)} vector containing the transition weight parameters (see below).}
#'     \item{\describe{
#'       \item{if \code{cond_dist = "Gaussian")}:}{Omit \eqn{\nu} from the parameter vector.}
#'       \item{if \code{cond_dist="Student"}:}{\eqn{\nu > 2} is the single degrees of freedom parameter.}
#'       \item{if \code{cond_dist="ind_Student"}:}{\eqn{\nu = (\nu_1,...,\nu_d)} \eqn{(d \times 1)}, \eqn{\nu_i > 2}.}
#'       \item{if \code{cond_dist="ind_skewed_t"}:}{\eqn{\nu = (\nu_1,...,\nu_d,\lambda_1,...,\lambda_d)} \eqn{(2d \times 1)},
#'        \eqn{\nu_i > 2} and \eqn{\lambda_i \in (0, 1)}.}
#'       }
#'     }
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
#'     \item{\code{weight_function="exogenous"}:}{Omit \eqn{\alpha} from the parameter vector.}
#'     \item{\code{identification="heteroskedasticity"}:}{\eqn{\sigma = (vec(W),\lambda_2,...,\lambda_M)}, where
#'           \eqn{W} \eqn{(d\times d)} and \eqn{\lambda_m} \eqn{(d\times 1)}, \eqn{m=2,...,M}, satisfy
#'           \eqn{\Omega_1=WW'} and \eqn{\Omega_m=W\Lambda_mW'}, \eqn{\Lambda_m=diag(\lambda_{m1},...,\lambda_{md})},
#'           \eqn{\lambda_{mi}>0}, \eqn{m=2,...,M}, \eqn{i=1,...,d}.}
#'   }
#'   Above, \eqn{\phi_{m,0}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}th coefficient matrix of the \eqn{m}th
#'   regime, \eqn{\Omega_{m}} denotes the positive definite error term covariance matrix of the \eqn{m}th regime, and \eqn{B_m}
#'   is the invertible \eqn{(d\times d)} impact matrix of the \eqn{m}th regime. \eqn{\nu_m} is the degrees of freedom parameter
#'   of the \eqn{m}th regime. If \code{parametrization=="mean"}, just replace each \eqn{\phi_{m,0}} with regimewise mean
#'   \eqn{\mu_{m}}.
#' @return Returns a \eqn{(d\times M)} matrix containing \eqn{\phi_{m,0}} in the m:th column or
#'   \eqn{\mu_{m}} if the parameter vector is mean-parametrized, \eqn{, m=1,..,M}.
#' @details Does not support constrained parameter vectors.
#' @section Warning:
#'  No argument checks!
#' @keywords internal

pick_phi0 <- function(M, d, params) {
  matrix(params[1:(d*M)], nrow=d, byrow=FALSE)
}


#' @title Pick coefficient matrix
#'
#' @description \code{pick_Ami} picks the coefficient matrix \eqn{A_{m,i}} from the given parameter vector.
#'
#' @inheritParams pick_phi0
#' @param p the autoregressive order of the model
#' @param m which regime?
#' @param i which lag in \eqn{1,...,p}?
#' @param unvec if \code{FALSE} then vectorized version of \eqn{A_{m,i}} will be returned instead of matrix.
#'   Default if \code{TRUE}.
#' @inherit pick_phi0 details
#' @inheritSection pick_phi0 Warning
#' @return Returns the i:th lag coefficient matrix of m:th regime, \eqn{A_{m,i}}.
#' @keywords internal

pick_Ami <- function(p, M, d, params, m, i, unvec=TRUE) {
  qm1 <- d*M + d^2*p*(m - 1)
  Ami <- params[(qm1 + d^2*(i - 1) + 1):(qm1 + d^2*i)]
  if(unvec) {
    return(unvec(d=d, a=Ami))
  } else {
    return(Ami)
  }
}


#' @title Pick coefficient matrices
#'
#' @description \code{pick_Am} picks the coefficient matrices \eqn{A_{m,i} (i=1,..,p)}
#'   from the given parameter vector for a given regime, so that they are arranged in
#'   a 3D array with the third dimension indicating each lag.
#'
#' @inheritParams pick_Ami
#' @return Returns a 3D array containing the coefficient matrices of the given regime.
#'  The coefficient matrix \eqn{A_{m,i}} can be obtained by choosing \code{[, , i]}.
#' @inherit pick_Ami details
#' @inheritSection pick_Ami Warning
#' @keywords internal

pick_Am <- function(p, M, d, params, m, structural_pars=NULL) {
  array(params[(d*M + d^2*p*(m - 1) + 1):(d*M + d^2*p*m)], dim=c(d, d, p))
}


#' @title Pick all coefficient matrices
#'
#' @description \code{pick_allA} picks all coefficient matrices \eqn{A_{m,i} (i=1,..,p, m=1,..,M)}
#'   from the given parameter vector so that they are arranged in a 4D array with the fourth dimension
#'   indicating each regime and third dimension indicating each lag.
#'
#' @inheritParams pick_Am
#' @return Returns a 4D array containing the coefficient matrices of the all components. Coefficient matrix
#'  \eqn{A_{m,i}} can be obtained by choosing \code{[, , i, m]}.
#' @inherit pick_Ami details
#' @inheritSection pick_Ami Warning
#' @keywords internal

pick_allA <- function(p, M, d, params) {
  array(params[(d*M + 1):(d*M + d^2*p*M)], dim=c(d, d, p, M))
}


#' @title Pick covariance matrices
#'
#' @description \code{pick_Omegas} picks the covariance matrices \eqn{\Omega_{m} (m=1,..,M)}
#'  from the given parameter vector so that they are arranged in a 3D array with the third
#'  dimension indicating each component.
#'
#' @inherit pick_Am
#' @return Returns a 3D array containing...
#'   \describe{
#'     \item{If \code{identification == "non-Gaussianity"} or \code{cond_dist \%in\% c("ind_Student", "ind_skewed_t")}:}{the impact
#'           matrices of the regimes, \eqn{B_m} in \code{[, , m]}.}
#'     \item{If otherwise:}{the covariance matrices of the given model, \eqn{\Omega_m} in \code{[, , m]}.}
#'   }
#' @details Constrained parameter vectors are not supported.
#' @inheritSection pick_Ami Warning
#' @keywords internal

pick_Omegas <- function(p, M, d, params, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t"),
                        identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity")) {
  cond_dist <- match.arg(cond_dist)
  identification <- match.arg(identification)
  Omegas <- array(dim=c(d, d, M))
  if(identification == "non-Gaussianity" || cond_dist == "ind_Student" || cond_dist == "ind_skewed_t") {
    for(m in 1:M) {
      Omegas[, , m] <- unvec(d=d, a=params[(M*d + d^2*p*M + (m - 1)*d^2 + 1):(M*d + d^2*p*M + m*d^2)])
    }
  } else if(identification == "heteroskedasticity") {
    W <- unvec(d=d, a=params[(M*d + d^2*p*M + 1):(M*d + d^2*p*M + d^2)])
    Omegas[, , 1] <- tcrossprod(W)
    for(m in 2:M) {
      lambdas <- params[(d*M*(1 + d*p) + d^2 + d*(m - 2) + 1):(d*M*(1 + d*p) + d^2 + d*(m - 1))]
      Omegas[, , m] <- W%*%tcrossprod(diag(lambdas), W)
    }
  } else { # Regular parameter vector
    qm1 <- d*M*(1 + p*d) + (1:M - 1)*d*(d + 1)/2
    for(m in 1:M) {
      Omegas[, , m] <- unvech(d=d, a=params[(qm1[m] + 1):(qm1[m] + d*(d + 1)/2)])
    }
  }
  Omegas
}


#' @title Pick transition weight parameters
#'
#' @description \code{pick_weightpars} picks the transition weight parameters from the given parameter vector.
#'
#' @inheritParams pick_Ami
#' @inheritParams loglikelihood
#' @return
#'   \describe{
#'     \item{If \code{weight_function = "relative_dens"}:}{Returns a length \eqn{M} vector containing the transition weight
#'           parameters \eqn{\alpha_{m}, m=1,...,M}, including the non-parametrized \eqn{\alpha_{M}}.}
#'    \item{\code{weight_function="logistic"}:}{Returns a length two vector \eqn{(c,\gamma)}, where
#'          \eqn{c\in\mathbb{R}} is the location parameter and \eqn{\gamma >0} is the scale parameter.}
#'     \item{If \code{weight_function = "mlogit"}:}{Returns a length \eqn{(M-1)k} vector \eqn{(\gamma_1,...,\gamma_M)},
#'           where \eqn{\gamma_m} \eqn{(k\times 1)}, \eqn{m=1,...,M-1} (\eqn{\gamma_M=0}) contains the mlogit-regression
#'           coefficients of the \eqn{m}th regime. Specifically, for switching variables with indices in
#'           \eqn{J\subset\lbrace 1,...,d\rbrace}, and with \eqn{\tilde{p}\in\lbrace 1,...,p\rbrace} lags included,
#'           \eqn{\gamma_m} contains the coefficients for the vector
#'           \eqn{z_{t-1} = (1,\tilde{z}_{\min\lbrace I\rbrace},...,\tilde{z}_{\max\lbrace I\rbrace})}, where
#'           \eqn{\tilde{z}_{i} =(y_{j,t-1},...,y_{j,t-\tilde{p}})}, \eqn{i\in I}. So \eqn{k=1+|I|\tilde{p}}
#'           where \eqn{|I|} denotes the number of elements in \eqn{I}.}
#'     \item{\code{weight_function="exponential"}:}{Returns a length two vector \eqn{(c,\gamma)}, where
#'           \eqn{c\in\mathbb{R}} is the location parameter and \eqn{\gamma >0} is the scale parameter.}
#'     \item{\code{weight_function="threshold"}:}{Returns a length \eqn{M-1} vector \eqn{(r_1,...,r_{M-1})}, where
#'           \eqn{r_1,...,r_{M-1}} are the threshold values.}
#'     \item{\code{weight_function="exogenous"}:}{Returns \code{numeric(0)}.}
#'   }
#' @inheritSection pick_Ami Warning
#' @keywords internal

pick_weightpars <- function(p, M, d, params,
                            weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                            weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t")) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  if(cond_dist == "Gaussian") {
    n_distpars <- 0
  } else if(cond_dist == "Student") {
    n_distpars <- 1
  } else if (cond_dist == "ind_Student") {
    n_distpars <- d
  } else { # Cond_dist == "ind_skewed_t"
    n_distpars <- 2*d
  }
  if(M == 1) {
    if(weight_function == "relative_dens") {
      return(1)
    } else {
      return(numeric(0))
    }
  }
  if(weight_function == "relative_dens") {
    alphas <- params[(length(params) - M - n_distpars + 2):(length(params) - n_distpars)]
    return(c(alphas, 1 - sum(alphas)))
  } else if(weight_function == "logistic" || weight_function == "exponential") {
    return(params[(length(params) - n_distpars - 1):(length(params) - n_distpars)]) # two params: c and gamma
  } else if(weight_function == "mlogit") {
    # (M-1)*k = (M-1)*|I|\tilde{p} pars to return
    return(params[(length(params) - (M - 1)*(1 + length(weightfun_pars[[1]])*weightfun_pars[[2]]) -
                     n_distpars + 1):(length(params) - n_distpars)])
  } else if(weight_function == "threshold") {
    return(params[(length(params) - M - n_distpars + 2):(length(params) - n_distpars)])
  } else if(weight_function == "exogenous") {
    return(numeric(0))
  }
}


#' @title Pick distribution parameters
#'
#' @description \code{pick_distpars} picks all the distribution parameters from
#'   the parameter vector
#'
#' @inheritParams pick_weightpars
#' @return Returns...
#'   \describe{
#'     \item{If \code{cond_dist == "Gaussian"}:}{a numeric vector of length zero.}
#'     \item{If \code{cond_dist == "Student"}:}{the degrees of freedom parameter.}
#'     \item{If \code{cond_dist == "ind_Student"}:}{a numeric vector of length \eqn{d} containing the degrees of freedom parameters.}
#'     \item{If \code{cond_dist == "ind_skewed_t"}:}{a numeric vector \eqn{(\nu_1,...,\nu_d,\lambda_1,...,\lambda_d)} of length \eqn{2d}
#'           containing the degrees of freedom and skewness parameters.}
#'   }
#' @keywords internal

pick_distpars <- function(d, params, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t")) {
  cond_dist <- match.arg(cond_dist)
  if(cond_dist == "Gaussian") {
    return(numeric(0))
  } else if(cond_dist == "Student") {
    return(params[length(params)])
  } else if(cond_dist == "ind_Student") {
    return(params[(length(params) - d + 1):length(params)])
  } else { # Cond_dist == "ind_skewed_t"
    return(params[(length(params) - 2*d + 1):length(params)])
  }
}


#' @title Pick regime parameters
#'
#' @description \code{pick_regime} picks the regime parameters
#'   \eqn{(\phi_{m,0},vec(A_{m,1}),...,vec(A_{m,p}),vech(\Omega_m))}
#' @inheritParams pick_Am
#' @details Constrained models nor structural models are supported.
#' @return Returns the vector...
#'   \describe{
#'     \item{If \code{identification == "non-Gaussianity"} or \code{cond_dist \%in\% c("ind_Student", "ind_skewed_t")}:}{
#'           \eqn{(\phi_{m,0},vec(A_{m,1}),...,vec(A_{m,p}),vec(B_m))}.}
#'     \item{If otherwise:}{\eqn{(\phi_{m,0},vec(A_{m,1}),...,vec(A_{m,p}),vech(\Omega_m))}.}
#'   }
#'   Note that neither weight parameters or distribution parameters are picked.
#' @keywords internal

pick_regime <- function(p, M, d, params, m, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t"),
                        identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity")) {
  identification <- match.arg(identification)
  cond_dist <- match.arg(cond_dist)
  if(identification == "non-Gaussianity" || cond_dist == "ind_Student" || cond_dist == "ind_skewed_t") {
    return(c(params[((m - 1)*d + 1):(m*d)],
             params[(M*d + (m - 1)*p*d^2 + 1):(M*d + m*p*d^2)],
             params[(M*d + M*p*d^2 + (m - 1)*d^2 + 1):(M*d + M*p*d^2 + m*d^2)]))
  } else {
    return(c(params[((m - 1)*d + 1):(m*d)],
             params[(M*d + (m - 1)*p*d^2 + 1):(M*d + m*p*d^2)],
             params[(M*d + M*p*d^2 + (m - 1)*d*(d + 1)/2 + 1):(M*d + M*p*d^2 + m*d*(d + 1)/2)]))
  }
}


#' @title Pick the structural parameter matrix W
#'
#' @description \code{pick_W} picks the structural parameter matrix W from the parameter vector
#'   of a structural model identified by heteroskedasticity.
#' @inheritParams loglikelihood
#' @details Constrained parameter vectors are not supported. Not even constraints in \eqn{W}!
#' @return Returns a \eqn{(d x d)} matrix \eqn{W} for structural models identified by heteroskedasticity
#'   and \code{NULL} for other models.
#' @references
#' \itemize{
#'    \item Lütkepohl H., Netšunajev A. 2017. Structural vector autoregressions with smooth transition in variances.
#'      \emph{Journal of Economic Dynamics & Control}, \strong{84}, 43-57.
#'  }
#' @keywords internal

pick_W <- function(p, M, d, params, identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity")) {
  identification <- match.arg(identification)
  if(identification != "heteroskedasticity") return(NULL)
  unvec(d=d, a=params[(M*d + d^2*p*M + 1):(M*d + d^2*p*M + d^2)])
}


#' @title Pick the structural parameter eigenvalues 'lambdas'
#'
#' @description \code{pick_lambdas} picks the structural parameters eigenvalue 'lambdas' from the parameter vector
#'   of a structural model identified by heteroskedasticity.
#' @inheritParams loglikelihood
#' @return Returns the length \code{(d*(M - 1))} vector \eqn{(\lambda_{2},...,\lambda_{M})}
#'  (see the argument \code{params}) for structural models identified by heteroskedasticity,
#'  \code{numeric(0)} if \eqn{M=1}, and \code{NULL} for other models.
#' @inherit pick_W details references
#' @keywords internal

pick_lambdas <- function(p, M, d, params, identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity")) {
  identification <- match.arg(identification)
  if(identification != "heteroskedasticity") {
    return(NULL)
  } else if(M == 1) {
    return(numeric(0))
  }
  params[(M*d + d^2*p*M + d^2 + 1):((M*d + d^2*p*M + d^2 + d*(M - 1)))]
}

