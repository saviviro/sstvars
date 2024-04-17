#' @title Create random VAR model \eqn{(dxd)} coefficient matrices \eqn{A}.
#'
#' @description \code{random_coefmats} generates random VAR model coefficient matrices.
#'
#' @inheritParams loglikelihood
#' @param how_many how many \eqn{(dxd)} coefficient matrices \eqn{A} should be drawn?
#' @param scale non-diagonal elements will be drawn from mean zero normal distribution
#'   with \code{sd=0.3/scale} and diagonal elements from one with \code{sd=0.6/scale}.
#'   Larger scale will hence more likely result stationary coefficient matrices, but
#'   will explore smaller area of the parameter space. Can be for example
#'   \code{1 + log(2*mean(c((p-0.2)^(1.25), d)))}.
#' @return Returns \eqn{((how_many*d^2)x1)} vector containing vectorized coefficient
#'  matrices \eqn{(vec(A_{1}),...,vec(A_{how_many}))}. Note that if \code{how_many==p},
#'  then the returned vector equals \strong{\eqn{\phi_{m}}}.
#' @keywords internal

random_coefmats <- function(d, how_many, scale) {
  as.vector(vapply(1:how_many, function(i1) {
    x <- rnorm(d*d, mean=0, sd=0.3/scale)
    x[1 + 0:(d - 1) * (d + 1)] <- rnorm(d, mean=0, sd=0.6/scale)
    x
  }, numeric(d*d)))
}


#' @title Create random stationary VAR model \eqn{(dxd)} coefficient matrices \eqn{A}.
#'
#' @description \code{random_coefmats2} generates random VAR model coefficient matrices.
#'
#' @inheritParams loglikelihood
#' @param ar_scale a positive real number. Larger values will typically result larger AR coefficients.
#' @details The coefficient matrices are generated using the algorithm proposed by Ansley
#'   and Kohn (1986) which forces stationarity. It's not clear in detail how \code{ar_scale}
#'   affects the coefficient matrices. Read the cited article by Ansley and Kohn (1986) and
#'   the source code for more information.
#'
#'   Note that when using large \code{ar_scale} with large \code{p} or \code{d}, numerical
#'   inaccuracies caused by the imprecision of the float-point presentation may result in errors
#'   or nonstationary AR-matrices. Using smaller \code{ar_scale} facilitates the usage of larger
#'   \code{p} or \code{d}.
#' @return Returns \eqn{((pd^2)x1)} vector containing stationary vectorized coefficient
#'  matrices \eqn{(vec(A_{1}),...,vec(A_{p})}.
#' @references
#'  \itemize{
#'    \item Ansley C.F., Kohn R. 1986. A note on reparameterizing a vector autoregressive
#'       moving average model to enforce stationarity.
#'       \emph{Journal of statistical computation and simulation}, \strong{24}:2, 99-106.
#'  }
#' @keywords internal

random_coefmats2 <- function(p, d, ar_scale=1) {
  # First generate matrices P_1,..,P_p with singular values less than one
  stopifnot(ar_scale > 0 && ar_scale <= 1)
  Id <- diag(nrow=d)
  all_P <- array(dim=c(d, d, p))
  for(i1 in 1:p) {
    A <- matrix(rnorm(d*d, sd=ar_scale), nrow=d)
    B <- t(chol(Id + tcrossprod(A, A)))
    all_P[, , i1] <- solve(B, A)
  }

  all_phi <- array(dim=c(d, d, p, p)) # [ , , i, j] for phi_{i, j}
  all_phi_star <- array(dim=c(d, d, p, p)) # [ , , i, j] for phi_{i, j}*

  # Set initial values
  L <- L_star <- Sigma <- Sigma_star <- Gamma <- Id

  # Recursion algorithm (Ansley and Kohn 1986, lemma 2.1)
  for(s in 0:(p - 1)) {
    all_phi[, , s+1, s+1] <- L%*%all_P[, , s+1]%*%solve(L_star)
    all_phi_star[, , s+1, s+1] <- tcrossprod(L_star, all_P[, , s+1])%*%solve(L)

    if(s >= 1) {
      for(k in 1:s) {
        all_phi[, , s+1, k] <- all_phi[, , s, k] - all_phi[, , s+1, s+1]%*%all_phi_star[, , s, s-k+1]
        all_phi_star[, , s+1, k] <- all_phi_star[, , s, k] - all_phi_star[, , s+1, s+1]%*%all_phi[, , s, s-k+1]
      }
    }

    if(s < p - 1) { # These are not needed in the last round because only coefficient matrices will be returned.
      Sigma_next <- Sigma - all_phi[, , s+1, s+1]%*%tcrossprod(Sigma_star, all_phi[, , s+1, s+1])
      if(s < p + 1) {
        Sigma_star <- Sigma_star - all_phi_star[, , s+1, s+1]%*%tcrossprod(Sigma, all_phi_star[, , s+1, s+1])
        L_star <- t(chol(Sigma_star))
      }
      Sigma <- Sigma_next
      L <- t(chol(Sigma))
    }
  }
  as.vector(all_phi[, , p, 1:p])
}


#' @title Create random VAR model error term covariance matrix
#'
#' @description \code{random_covmat} generates random VAR model \eqn{(dxd)} error term covariance matrix \eqn{\Omega}
#'   from (scaled) Wishart distribution.
#'
#' @inheritParams loglikelihood
#' @inheritParams GAfit
#' @return Returns a \eqn{(d(d+1)/2x1)} vector containing vech-vectorized covariance matrix
#'       \eqn{\Omega}.
#' @keywords internal

random_covmat <- function(d, omega_scale) {
  smart_covmat(d=d, Omega=diag(x=omega_scale), accuracy=1)
}



#' @title Create random VAR model \eqn{(dxd)} error term covariance matrix \eqn{\Omega}
#'   fairly close to the given \strong{positive definite} covariance matrix using (scaled)
#'   Wishart distribution
#'
#' @description \code{smart_covmat} generates random VAR model \eqn{(dxd)} error term covariance matrix \eqn{\Omega}
#'   from (scaled) Wishart distribution that is fairly close to the given matrix.
#'
#' @inheritParams loglikelihood
#' @param Omega a symmetric positive definite \eqn{(dxd)} covariance matrix specifying
#'   expected value of the matrix to be generated.
#' @param accuracy a positive real number adjusting how close to the given covariance matrix
#'   the returned individual should be.
#'
#'   The standard deviation of each diagonal element is...
#'   \itemize{
#'     \item \eqn{\omega_{i,i}/}\code{accuracy} when \code{accuracy > d/2}
#'     \item and \code{sqrt(2/d)*}\eqn{\omega_{i,i}} when \code{accuracy <= d/2}.
#'   }
#'   Wishart distribution is used for reduced form models, but for more details read the source code.
#' @inherit random_covmat return
#' @keywords internal

smart_covmat <- function(d, Omega, accuracy) {
  if(accuracy <= d/2) {
    covmat <- rWishart(n=1, df=d, Sigma=Omega/d)[, , 1]
  } else {
    covmat <- (1 - d/(2*accuracy))*Omega + rWishart(n=1, df=(d^2)/2, Sigma=Omega/(d*accuracy))[, , 1]
  }
  vech(covmat)
}



#' @title Create random VAR model impact matrix
#'
#' @description \code{random_impactmat} generates random VAR model \eqn{(dxd)} impact matrix \eqn{B}
#'   with its elements drawn from specific normal distributions (see the source code). If not the first
#'   regime, will create the matrix \eqn{B_m*}.
#'
#' @inheritParams loglikelihood
#' @inheritParams GAfit
#' @param is_regime1 is the impact matrix for Regime 1? Regime 1 impact matrix is constrained so the elements
#'   in its first row are in a decreasing ordering and the diagonal elements are strictly positive.
#' @details If the impact matrix is not for Regime 1, will create the matrix \eqn{B_m*}, which is related
#'   to the impact matrix \eqn{B_m} of Regime m as \eqn{B_m* = B_m - B_1}.
#' @return Returns a \eqn{(d^2 \times 1)} vector containing the vectorized impact matrix \eqn{B}.
#' @keywords internal

random_impactmat <- function(d, B_scale, is_regime1=TRUE) {
  new_B <- matrix(nrow=d, ncol=d)
  if(is_regime1) {
    for(i1 in 1:d) {
      new_B[i1,] <- rnorm(d, sd=sqrt(0.95*B_scale[i1]/d)) # Random elements to each row
    }
    return(order_B(new_B)) # First element in each column normalized to positive and columns ordered to a decreasing order
  } else {
    for(i1 in 1:d) {
      new_B[i1,] <- rnorm(d, sd=sqrt(0.05*B_scale[i1]/d)) # Random elements to each row
    }
    return(new_B) # No constraints since not the first impact matrix
  }
}


#' @title Create a random VAR model \eqn{(dxd)} error impact matrix \eqn{B}
#'   fairly close to the given \strong{invertible} impact matrix.
#'
#' @description \code{smart_impactmat} generates a random VAR model \eqn{(dxd)} error impact matrix \eqn{B}
#'   fairly close to the given \strong{invertible} impact matrix.
#'
#' @inheritParams random_impactmat
#' @param B an invertible \eqn{(dxd)} impact matrix specifying
#'   expected value of the matrix to be generated. Should have strictly positive diagonal entries in
#'   a decreasing order.
#' @param accuracy a positive real number adjusting how close to the given impact matrix
#'   the returned individual should be. Larger value implies higher accuracy.
#' @inherit random_impactmat return
#' @keywords internal

smart_impactmat <- function(d, B, accuracy, is_regime1=TRUE) {
  new_B <- matrix(rnorm(d^2, mean=B, sd=pmax(0.2, abs(B))/accuracy), nrow=d)
  if(is_regime1) {
    return(order_B(new_B)) # First element in each column normalized to positive and columns ordered to a decreasing order
  } else {
    return(new_B) # No constraints since not the first impact matrix
  }
}


#' @title Create random distribution parameter values
#'
#' @description \code{random_distpars} generates random distribution parameter values
#'
#' @inheritParams loglikelihood
#' @return Returns a numeric vector ...
#'   \describe{
#'     \item{If \code{cond_dist == "Gaussian"}:}{of length zero.}
#'     \item{If \code{cond_dist == "Student"}:}{of length one containing a df param strictly larger than two.}
#'     \item{If \code{cond_dist == "ind_Student"}:}{of length d containing a df params strictly larger than two.}
#'   }
#' @keywords internal

random_distpars <- function(d, cond_dist) {
  if(cond_dist == "Gaussian") {
    return(numeric(0))
  } else if(cond_dist == "Student") {
    return(2.000001 + rgamma(1, shape=0.3, rate=0.007))
  } else if(cond_dist == "ind_Student") {
    return(2.000001 + rgamma(d, shape=0.3, rate=0.007))
  }
}


#' @title Create random distribution parameter values close to given values
#'
#' @description \code{smart_distpars} generates random degrees of freedom parameter values
#'   close to given values.
#'
#' @param distpars the old distribution parameters (of all regimes)
#' @param accuracy a positive real number adjusting how close to the given distribution parameters
#'   the returned values should be.
#' @return Returns a numeric vector ...
#'   \describe{
#'     \item{If \code{cond_dist == "Gaussian"}:}{of length zero.}
#'     \item{If \code{cond_dist == "Student"}:}{of length one containing a degrees of freedom parameter value
#'      (strictly larger than two).}
#'     \item{If \code{cond_dist == "ind_Student"}:}{of length d containing a degrees of freedom parameter values
#'      (strictly larger than two).}
#'   }
#' @keywords internal

smart_distpars <- function(distpars, accuracy, cond_dist) {
  if(cond_dist == "Gaussian") {
    return(numeric(0))
  } else if(cond_dist == "Student" || cond_dist == "ind_Student") {
    new_distpars <- rnorm(length(distpars), mean=distpars, sd=pmax(0.2, abs(distpars))/accuracy) # smart df
    return(pmax(2.01, new_distpars)) # Make sure all df are above the strict lower bound 2
  }
}


#' @title Create random transition weight parameter values
#'
#' @description \code{random_weightpars} generates random transition weight parameter values
#'
#' @inheritParams random_ind
#' @return Returns a numeric vector ...
#'   \describe{
#'     \item{If \code{weight_function == "relative_dens"}:}{a length \code{M-1} vector \eqn{(\alpha_1,...,\alpha_{M-1})}.}
#'     \item{If \code{weight_function == "logistic"}:}{a length two vector \eqn{(c,\gamma)},
#'           where \eqn{c\in\mathbb{R}} is the location parameter and \eqn{\gamma >0} is the scale parameter.}
#'     \item{If \code{weight_function == "mlogit"}:}{a length \eqn{((M-1)k\times 1)} vector \eqn{(\gamma_1,...,\gamma_{M-1})},
#'           where \eqn{\gamma_m} \eqn{(k\times 1)}, \eqn{m=1,...,M-1} contains the mlogit-regression coefficients of the \eqn{m}th
#'           regime. Specifically, for switching variables with indices in \eqn{I\subset\lbrace 1,...,d\rbrace}, and with
#'           \eqn{\tilde{p}\in\lbrace 1,...,p\rbrace} lags included, \eqn{\gamma_m} contains the coefficients for the vector
#'           \eqn{z_{t-1} = (1,\tilde{z}_{\min\lbrace I\rbrace},...,\tilde{z}_{\max\lbrace I\rbrace})}, where
#'           \eqn{\tilde{z}_{i} =(y_{it-1},...,y_{it-\tilde{p}})}, \eqn{i\in I}. So \eqn{k=1+|I|\tilde{p}}
#'           where \eqn{|I|} denotes the number of elements in \eqn{I}.}
#'     \item{If \code{weight_function == "exponential"}:}{a length two vector \eqn{(c,\gamma)},
#'           where \eqn{c\in\mathbb{R}} is the location parameter and \eqn{\gamma >0} is the scale parameter.}
#'     \item{If \code{weight_function == "threshold"}:}{a length \eqn{M-1} vector \eqn{(r_1,...,r_{M-1})},
#'           where \eqn{r_1,...,r_{M-1}} are the threshold values in an increasing order.}
#'     \item{If \code{weight_function == "exogenous"}:}{of length zero.}
#'   }
#' @keywords internal

random_weightpars <- function(M, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                              weightfun_pars=NULL, AR_constraints=NULL, mean_constraints=NULL,
                              weight_constraints=NULL, weight_scale) {
  weight_function <- match.arg(weight_function)
  if(M == 1 || weight_function == "exogenous") {
    return(numeric(0))
  }
  if(weight_function == "relative_dens") {
    if(is.null(weight_constraints)) {
      alphas <- runif(n=M, min=0.00, max=1)
      # Sort and standardize alphas; don't sort if AR_constraints or mean_constraints are used
      if(is.null(AR_constraints) && is.null(mean_constraints)) {
        alphas <- sort(alphas, decreasing=TRUE)
        alphas <- alphas/sum(alphas)
      }
      ret <- alphas[-M]
    } else {
      if(all(weight_constraints[[1]] == 0)) {
        ret <- numeric(0) # alpha = r constant, so it is not parametrized
      } else {
        ret <- sort(runif(n=ncol(weight_constraints[[1]]), min=0, max=0.8), decreasing=TRUE)
      }
    }
  } else if(weight_function == "logistic" || weight_function == "exponential") {
    if(is.null(weight_constraints)) {
      ret <- c(rnorm(n=1, mean=weight_scale[1], sd=weight_scale[2]), abs(rnorm(n=1, mean=0, sd=weight_scale[3])) + 1e-8)
    } else {
      if(all(weight_constraints[[1]] == 0)) {
        ret <- numeric(0) # alpha = r constant, so it is not parametrized
      } else {
        ret <- rnorm(n=ncol(weight_constraints[[1]]), mean=0, sd=weight_scale[3])
      }
    }
  } else if(weight_function == "mlogit") {
    if(is.null(weight_constraints)) {
      k <- 1 + length(weightfun_pars[[1]])*weightfun_pars[[2]] # how many pars in each gamma_m
      ret <- rnorm(n=(M - 1)*k, mean=0, sd=weight_scale[3])
      # Replace different coefficients for constant terms:
      ret[(1:(M - 1) - 1)*k + 1] <- rnorm(n=(M - 1), mean=weight_scale[1], sd=weight_scale[2])
    } else {
      if(all(weight_constraints[[1]] == 0)) {
        ret <- numeric(0) # alpha = r constant, so it is not parametrized
      } else {
        ret <- rnorm(n=ncol(weight_constraints[[1]]), mean=0, sd=weight_scale[3])
      }
    }
  } else if(weight_function == "threshold") {
    if(is.null(weight_constraints)) {
      ret <- sort(runif(n=M-1, min=weight_scale[1], max=weight_scale[2]), decreasing=FALSE)
    } else {
      if(all(weight_constraints[[1]] == 0)) {
        ret <- numeric(0) # alpha = r constant, so it is not parametrized
      } else {
        # We don't sort here so that we don't accidentally always be outside the parameter space with some weird weight constraints
        ret <- runif(n=ncol(weight_constraints[[1]]), min=weight_scale[1], max=weight_scale[2])
      }
    }
  }
  ret
}

#' @title Create random transition weight parameter values
#'
#' @description \code{smart_weightpars} generates random transition weight parameter values
#'  relatively close to the ones given in \code{weight_pars}
#'
#' @inheritParams smart_ind
#' @param weight_pars a vector containing transition weight parameter values.
#'   \describe{
#'     \item{If \code{weight_function == "relative_dens"}:}{a length \code{M-1} vector \eqn{(\alpha_1,...,\alpha_{M-1})}.}
#'     \item{If \code{weight_function == "logistic"}:}{a length two vector \eqn{(c,\gamma)},
#'           where \eqn{c\in\mathbb{R}} is the location parameter and \eqn{\gamma >0} is the scale parameter.}
#'     \item{If \code{weight_function == "mlogit"}:}{a length \eqn{((M-1)k\times 1)} vector \eqn{(\gamma_1,...,\gamma_{M-1})},
#'           where \eqn{\gamma_m} \eqn{(k\times 1)}, \eqn{m=1,...,M-1} contains the mlogit-regression coefficients of the \eqn{m}th
#'           regime. Specifically, for switching variables with indices in \eqn{I\subset\lbrace 1,...,d\rbrace}, and with
#'           \eqn{\tilde{p}\in\lbrace 1,...,p\rbrace} lags included, \eqn{\gamma_m} contains the coefficients for the vector
#'           \eqn{z_{t-1} = (1,\tilde{z}_{\min\lbrace I\rbrace},...,\tilde{z}_{\max\lbrace I\rbrace})}, where
#'           \eqn{\tilde{z}_{i} =(y_{it-1},...,y_{it-\tilde{p}})}, \eqn{i\in I}. So \eqn{k=1+|I|\tilde{p}}
#'           where \eqn{|I|} denotes the number of elements in \eqn{I}.}
#'     \item{If \code{weight_function == "exponential"}:}{a length two vector \eqn{(c,\gamma)},
#'           where \eqn{c\in\mathbb{R}} is the location parameter and \eqn{\gamma >0} is the scale parameter.}
#'     \item{If \code{weight_function == "threshold"}:}{a length \eqn{M-1} vector \eqn{(r_1,...,r_{M-1})},
#'           where \eqn{r_1,...,r_{M-1}} are the threshold values in an increasing order.}
#'     \item{If \code{weight_function == "exogenous"}:}{of length zero.}
#'   }
#' @inherit random_weightpars return
#' @keywords internal

smart_weightpars <- function(M, weight_pars,
                             weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                             weight_constraints=NULL, accuracy) {
  weight_function <- match.arg(weight_function)
  if(M == 1 || weight_function == "exogenous") {
    return(numeric(0))
  }
  if(!is.null(weight_constraints)) {
    if(all(weight_constraints[[1]] == 0)) {
      return(numeric(0)) # alpha = r constant, so it is not parametrized
    }
  }
  if(weight_function == "relative_dens") {
    # Note: this assumes that the unparametrized alpha is not in weight_pars!
    weight_pars <- abs(rnorm(n=length(weight_pars) + 1, mean=c(weight_pars, 1 - sum(weight_pars)),
                         sd=pmax(0.2, c(weight_pars, 1 - sum(weight_pars)))/accuracy))
    ret <- (weight_pars/sum(weight_pars))[-M] # Standardize alphas; sorts in GA in sort_regimes if no constraints are imposed
  } else if(weight_function == "mlogit" || weight_function == "logistic" || weight_function == "exponential") {
    ret <- rnorm(n=length(weight_pars), mean=weight_pars, sd=pmax(0.2, weight_pars)/accuracy)
  } else if(weight_function == "threshold") {
    # Threshold pars always need to be sorted;
    ret <- sort(rnorm(n=length(weight_pars), mean=weight_pars, sd=pmax(0.2, weight_pars)/accuracy), decreasing=FALSE)
  }
  ret
}



#' @title Create random mean parametrized parameter vector
#'
#' @description \code{random_ind} generates random mean parametrized parameter vector.
#'
#' @inheritParams GAfit
#' @inheritParams loglikelihood
#' @param force_stability Should the algorithm proposed by Ansley and Kohn be used to generate
#'   AR matrices that always satisfy the stability condition? Not supported if AR constraints are
#'   employed.
#' @details Structural models are not supported!
#' @return Returns random mean parametrized parameter vector that has the same form as the argument \code{params}
#'   in the other functions, for instance, in the function \code{loglikelihood}.
#' @references
#'  \itemize{
#'    \item Ansley C.F., Kohn R. 1986. A note on reparameterizing a vector autoregressive
#'       moving average model to enforce stationarity.
#'       \emph{Journal of statistical computation and simulation}, \strong{24}:2, 99-106.
#'  }
#' @keywords internal

random_ind <- function(p, M, d, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                       weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student"), AR_constraints=NULL, mean_constraints=NULL,
                       weight_constraints=NULL, force_stability=is.null(AR_constraints),
                       mu_scale, mu_scale2, omega_scale, B_scale, weight_scale, ar_scale=1, ar_scale2=1) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  g <- ifelse(is.null(mean_constraints), M, length(mean_constraints)) # Number of groups of regimes with the same mean parameters

  # Generate mean params
  if(is.null(mean_constraints)) {
    mean_pars <- as.vector(replicate(n=M, expr=rnorm(d, mean=mu_scale, sd=mu_scale2)))
  } else { # mean_constraints used
    mean_pars <- rnorm(d*g, mean=mu_scale, sd=mu_scale2)
  }

  # Generate AR params
  if(is.null(AR_constraints) && force_stability) {
    coefmat_pars <- as.vector(replicate(n=M, random_coefmats2(p=p, d=d, ar_scale=ar_scale)))
  } else {
    scale_A <- ar_scale2*ifelse(is.null(AR_constraints),
                                1 + log(2*mean(c((p - 0.2)^(1.25), d))),
                                1 + (sum(AR_constraints)/(M*d^2))^0.85)
    if(is.null(AR_constraints)) {
      coefmat_pars <- as.vector(replicate(n=M, expr=random_coefmats(d=d, how_many=p, scale=scale_A)))
    } else { # AR_constraints employed
      coefmat_pars <- rnorm(ncol(AR_constraints), mean=0, sd=0.5/scale_A) # random psi
    }
  }

  # Generate covmat params
  if(cond_dist == "ind_Student") { # Covmat pars impact matrix params
    if(M == 1) {
      covmat_pars <- as.vector(random_impactmat(d=d, B_scale=B_scale, is_regime1=TRUE))
    } else {
      covmat_pars <- c(random_impactmat(d=d, B_scale=B_scale, is_regime1=TRUE), # Regime 1 impact matrix is constrained
                               replicate(n=M-1, expr=random_impactmat(d=d, B_scale=B_scale, is_regime1=FALSE))) # Regime 2,...,M impact matrices
    }
  } else { # cond_dist == "Gaussian" or "Student"
    covmat_pars <- as.vector(replicate(n=M, expr=random_covmat(d=d, omega_scale=omega_scale)))
  }

  # Generate weight params
  if(M > 1) {
    weight_pars <- random_weightpars(M=M, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                     AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                     weight_constraints=weight_constraints, weight_scale=weight_scale)
  } else {
    weight_pars <- numeric(0) # Replicate returns a list if the value is numeric(0)
  }

  # Generate distribution params (df etc)
  dist_pars <- random_distpars(d=d, cond_dist=cond_dist)

  # Return the parameter vector
  c(mean_pars, coefmat_pars, covmat_pars, weight_pars, dist_pars)
}


#' @title Create random parameter vector that is fairly close to a given parameter vector
#'
#' @description \code{smart_ind} creates random mean parametrized parameter vector that is
#'   model fairly close to a given parameter vector. The result may not be satisfy the stability
#'   condition.
#'
#' @inheritParams loglikelihood
#' @inheritParams GAfit
#' @inheritParams random_coefmats2
#' @param accuracy a positive real number adjusting how close to the given parameter vector the returned individual should be.
#'   Larger number means larger accuracy. Read the source code for details.
#' @param which_random a vector with length between 1 and M specifying the mixture components that should be random instead of
#'   close to the given parameter vector. This does not consider constrained AR or lambda parameters.
#' @details Structural models are not supported!
#' @inherit random_ind return references
#' @keywords internal

smart_ind <- function(p, M, d, params, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                      weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student"), AR_constraints=NULL, mean_constraints=NULL,
                      weight_constraints=NULL,  accuracy=1, which_random=numeric(0),
                      mu_scale, mu_scale2, omega_scale, B_scale, ar_scale=1, ar_scale2=1) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  scale_A <- ar_scale2*(1 + log(2*mean(c((p - 0.2)^(1.25), d))))
  params_std <- reform_constrained_pars(p=p, M=M, d=d, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                        cond_dist=cond_dist, identification="reduced_form", AR_constraints=AR_constraints,
                                        mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                        B_constraints=NULL) # Used so that pick_pars-functions works
  dist_pars <- pick_distpars(d=d, params=params_std, cond_dist=cond_dist)
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params_std, cond_dist=cond_dist,
                           identification="reduced_form") # Impact matrices for cond_dist == "ind_Student"
  new_pars <- numeric(length(params))
  if(is.null(AR_constraints) && is.null(mean_constraints) && is.null(weight_constraints)) {
    all_means_and_A <- params[1:(d*M + M*p*d^2)] # all mu + A if called from GAfit
    new_pars[1:(d*M + M*p*d^2)] <- rnorm(n=length(all_means_and_A), # Fills in smart means and A, later replaced with random if random regime
                                         mean=all_means_and_A, sd=pmax(0.2, abs(all_means_and_A))/accuracy)
    for(m in 1:M) {
      if(any(which_random) == m) { # Not a smart regime
        new_pars[((m - 1)*d + 1):(m*d)] <- rnorm(d, mean=mu_scale, sd=mu_scale2) # Random mean
        if(runif(1) > 0.5) {
          new_pars[(M*d + (m - 1)*p*d^2 + 1):(M*d + m*p*d^2)] <- random_coefmats(d=d, how_many=p, scale=scale_A) # Without algorithm
        } else {
          new_pars[(M*d + (m - 1)*p*d^2 + 1):(M*d + m*p*d^2)] <- random_coefmats2(p=p, d=d, ar_scale=ar_scale) # Use the algorithm
        }
        if(cond_dist == "ind_Student") { # Impact matrix parameters
          new_pars[(M*d + M*p*d^2 + (m - 1)*d^2 + 1):(M*d + M*p*d^2 + m*d^2)] <- random_impactmat(d=d, B_scale=B_scale, is_regime1=(m == 1))
        } else { # cond_dist == "Gaussian" or "Student", cov mat parameters
          new_pars[(M*d + M*p*d^2 + (m - 1)*d*(d + 1)/2 + 1):(M*d + M*p*d^2 + m*d*(d + 1)/2)] <- random_covmat(d=d, omega_scale=omega_scale)
        }
      } else { # Smart_regime, smart mean and AR parameters already generated
        if(cond_dist == "ind_Student") {  # Impact matrix parameters
          new_pars[(M*d + M*p*d^2 + (m - 1)*d^2 + 1):(M*d + M*p*d^2 + m*d^2)] <- smart_impactmat(d=d, B=all_Omega[, , m],
                                                                                                 accuracy=accuracy,
                                                                                                 is_regime1=(m == 1))
        } else { # cond_dist == "Gaussian" or "Student", cov mat parameters
          new_pars[(M*d + M*p*d^2 + (m - 1)*d*(d + 1)/2 + 1):(M*d + M*p*d^2 + m*d*(d + 1)/2)] <- smart_covmat(d=d, Omega=all_Omega[, , m],
                                                                                                              accuracy=accuracy)
        }
      }
    }
    if(M > 1 && weight_function != "exogenous") {
      weight_pars <- pick_weightpars(p=p, M=M, d=d, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                     cond_dist=cond_dist)
      if(weight_function == "relative_dens") {
        weight_pars <- weight_pars[-length(weight_pars)] # Remove alpha_M
      }
      new_pars[(M*d + M*p*d^2 + M*d*(d + 1)/2 + 1):
                 (M*d + M*p*d^2 + M*d*(d + 1)/2 + length(weight_pars))] <- smart_weightpars(M=M,
                                                                                            weight_pars=weight_pars,
                                                                                            weight_function=weight_function,
                                                                                            weight_constraints=weight_constraints,
                                                                                            accuracy=accuracy)
    }
  } else { # AR, mean, or weight constraints employed
    # mean parameters
    g <- ifelse(is.null(mean_constraints), M, length(mean_constraints)) # Number of groups of regimes with the same mean parameters
    if(length(which_random) == 0) {
      smart_regs <- 1:M
    } else {
      smart_regs <- (1:M)[-which_random]
    }
    all_means <- matrix(params[1:(d*g)], nrow=d, ncol=g, byrow=FALSE)
    mean_pars <- as.vector(vapply(1:g, function(m) {
      which_reg <- ifelse(is.null(mean_constraints), m, mean_constraints[[m]]) # Can be many if mean_constraints used
      if(any(which_reg %in% smart_regs)) { # Smart parameters
        rnorm(d, mean=all_means[,m], sd=abs(all_means[,m]/accuracy))
      } else { # Random parameters
        rnorm(d, mean=mu_scale, sd=mu_scale2)
      }
    }, numeric(d)))
    # AR parameters
    if(is.null(AR_constraints)) {
      q <- M*p*d^2
      all_A <- pick_allA(p=p, M=M, d=d, params=params_std)
      AR_pars <- vapply(1:M, function(m) {
        if(any(which_random) == m) { # Random AR matrix
          if(runif(1) > 0.5) { # Use algorithm to force stability of coefficient matrices
            random_coefmats2(p=p, d=d, ar_scale=ar_scale)
          } else {
            random_coefmats(d=d, how_many=p, scale=scale_A)
          }
        } else { # Smart AR matrix
          all_Am <- as.vector(all_A[, , , m])
          rnorm(n=length(all_Am), mean=all_Am, sd=pmax(0.2, abs(all_Am))/accuracy)
        }
      }, numeric(p*d^2))
    } else {
      q <- ncol(AR_constraints)
      psi <- params[(g*d + 1):(g*d + q)]
      AR_pars <- rnorm(q, mean=psi, sd=pmax(0.2, abs(psi))/accuracy)
    }

    # Covariance matrix / impact matrix parameters
    if(cond_dist == "ind_Student") { # impact mat params
      covmat_pars <- as.vector(vapply(1:M, function(m) {
        if(any(which_random == m)) {
          random_impactmat(d=d, B_scale=B_scale, is_regime1=(m == 1))
        } else {
          smart_impactmat(d=d, B=all_Omega[, , m], accuracy=accuracy, is_regime1=(m == 1))
        }
      }, numeric(d^2)))
    } else { # cond_dist == "Gaussian" or "Student", cov mat params
      covmat_pars <- as.vector(vapply(1:M, function(m) {
        if(any(which_random == m)) {
          random_covmat(d=d, omega_scale=omega_scale)
        } else {
          smart_covmat(d=d, Omega=all_Omega[, , m], accuracy=accuracy)
        }
      }, numeric(d*(d + 1)/2)))
    }

    # Weight function parameters
    if(M > 1 && weight_function != "exogenous") {
      if(is.null(weight_constraints)) {
        weight_pars <- pick_weightpars(p=p, M=M, d=d, params=params_std, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                       cond_dist=cond_dist)
        if(weight_function == "relative_dens") {
          weight_pars <- weight_pars[-length(weight_pars)]
        }
        weight_pars <- smart_weightpars(M=M, weight_pars=weight_pars, weight_function=weight_function,
                                        weight_constraints=weight_constraints, accuracy=accuracy)
      } else { # Weight constraints used
        if(all(weight_constraints[[1]] == 0)) {
          weight_pars <- numeric(0) # alpha = r known constant so it is not parametrized here
        } else {
          weight_pars <- smart_weightpars(M=M,
                                          weight_pars=params[(M*g + q + M*d*(d + 1)/2 + 1):(M*g + q + M*d*(d + 1)/2 +
                                                                                              ncol(weight_constraints[[1]]))],
                                          weight_function=weight_function,
                                          weight_constraints=weight_constraints, accuracy=accuracy)
        }
      }
    } else {
      weight_pars <- numeric(0)
    }
    new_pars[1:(d*g + q + length(covmat_pars) + length(weight_pars))] <- c(mean_pars, AR_pars, covmat_pars, weight_pars)
  }
  if(cond_dist == "Student") {
    # Add degrees of freedom parameter
    new_pars[length(new_pars)] <- smart_distpars(distpars=dist_pars, accuracy=accuracy, cond_dist=cond_dist)
  } else if(cond_dist == "ind_Student") {
    # Add degrees of freedom parameters
    new_pars[(length(new_pars) - d + 1):length(new_pars)] <- smart_distpars(distpars=dist_pars, accuracy=accuracy, cond_dist=cond_dist)
  }
  new_pars
}
