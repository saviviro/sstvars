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
  stopifnot(ar_scale > 0)
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
#'   fairly close to a given \strong{positive definite} covariance matrix using (scaled)
#'   Wishart distribution
#'
#' @description \code{random_covmat} generates random VAR model \eqn{(dxd)} error term covariance matrix \eqn{\Omega}
#'   from (scaled) Wishart distribution that is fairly close to the given matrix.
#'
#' @inheritParams loglikelihood
#' @param Omega a symmetric positive definite \eqn{(dxd)} covariance matrix specifying
#'   expected value of the matrix to be generated.
#' @param accuracy a positive real number adjusting how close to the given covariance matrix
#'   the returned individual should be.
#'
#'   For \strong{reduced form models} standard deviation of each diagonal element is for reduced form
#'   models
#'   \itemize{
#'     \item \eqn{\omega_{i,i}/}\code{accuracy} when \code{accuracy > d/2}
#'     \item and \code{sqrt(2/d)*}\eqn{\omega_{i,i}} when \code{accuracy <= d/2}.
#'   }
#'   Wishart distribution is used for reduced form models, but for more details read the source code.
#' @inherit random_covmat return
#' @keywords internal

smart_covmat <- function(d, M, Omega, accuracy) {
  if(accuracy <= d/2) {
    covmat <- rWishart(n=1, df=d, Sigma=Omega/d)[, , 1]
  } else {
    covmat <- (1 - d/(2*accuracy))*Omega + rWishart(n=1, df=(d^2)/2, Sigma=Omega/(d*accuracy))[, , 1]
  }
  vech(covmat)
}


#' @title Create random distribution parameter values
#'
#' @description \code{random_distpars} generates random distribution parameter values
#'
#' @inheritParams loglikelihood
#' @return Returns a numeric vector ...
#'   \decribe{
#'     \item{If \code{cond_dist == "Gaussian"}:}{of length zero.}
#'     \item{If \code{cond_dist == "Student"}:}{of length one containing a df param strictly larger than two.}
#'   }
#' @keywords internal

random_distpars <- function(cond_dist) {
  if(weight_function == "Gaussian") {
    return(numeric(0))
  } else if(weight_function == "Student") {
    return(2.000001 + rgamma(1, shape=0.3, rate=0.007))
  }
}


#' @title Create random degrees of freedom parameter values close to given values
#'
#' @description \code{smart_df} generates random degrees of freedom parameter values
#'   close to given values.
#'
#' @param df the old degrees of freedom parameters (of all regimes)
#' @param accuracy a positive real number adjusting how close to the given degrees of freedom parameters
#'   the returned df should be.
#' @return a numeric vector containing the degrees of freedom parameter values
#' @keywords internal

smart_df <- function(df, accuracy) {
  if(length(df) == 0) {
    return(numeric(0))
  }
  new_df <- rnorm(length(df), mean=df, sd=pmax(0.2, abs(df))/accuracy) # smart df
  pmax(2.01, new_df) # Make sure all df are above the strict lower bound 2
}



#' @title Create random transition weight parameter values
#'
#' @description \code{random_weightpars} generates random transition weight parameter values
#'
#' @inheritParams random_ind
#' @return Returns a numeric vector ...
#'   \describe{
#'     \item{If \code{weight_function == "relative_dens"}:}{a length \code{M-1} vector \eqn{(\alpha_1,...,\alpha_{M-1})}.}
#'     \item{If \code{weight_function == "logit"}:}{NOT YET IMPLEMENTED}
#'   }
#' @keywords internal

random_weightpars <- function(M, weight_function, AR_constraints, mean_constraints) {
  if(weight_function == "relative_dens") {
    alphas <- runif(n=M)
    # Sort and standardize alphas; don't sort if AR_constraints or mean_constraints are used
    if(is.null(AR_constraints) && is.null(mean_constraints)) {
      alphas <- sort(alphas, decreasing=TRUE)
      alphas <- alphas/sum(alphas)
    }
    ret <- alphas[-M]
  } else if(weight_function == "logit") {
    stop("Logit weights not yet implemented in random_weightpars")
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
#'     \item{If \code{weight_function == "logit"}:}{NOT YET IMPLEMENTED}
#'   }
#' @return Returns a numeric vector ...
#'   \describe{
#'     \item{If \code{weight_function == "relative_dens"}:}{a length \code{M-1} vector \eqn{(\alpha_1,...,\alpha_{M-1})}.}
#'     \item{If \code{weight_function == "logit"}:}{NOT YET IMPLEMENTED}
#'   }
#' @keywords internal

smart_weightpars <- function(M, weight_pars, weight_function, accuracy) {
  weight_pars <- rnorm(n=length(weight_pars) + 1, mean=c(weight_pars, 1-sum(weight_pars)),
                       sd=pmax(0.2, c(weight_pars, 1-sum(weight_pars)))/accuracy)
  if(weight_function == "relative_dens") {
    ret <- (weight_pars/sum(weight_pars))[-M]
    # Sort and standardize alphas; don't sort if AR_constraints or mean_constraints are used
  } else if(weight_function == "logit") {
    stop("Logit weights not yet implemented in random_weightpars")
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

random_ind <- function(p, M, d, weight_function=c("relative_dens", "logit"), cond_dist=c("Gaussian", "Student"),
                       AR_constraints=NULL, mean_constraints=NULL, force_stability=is.null(AR_constraints),
                       mu_scale, mu_scale2, omega_scale, ar_scale=1, ar_scale2=1) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  #g <- ifelse(is.null(same_means), M, length(same_means)) # Number of groups of regimes with the same mean parameters

  # Generate mean params
  if(is.null(mean_constraints)) {
    mean_pars <- as.vector(replicate(n=M, expr=rnorm(d, mean=mu_scale, sd=mu_scale2)))
  } else {
    stop("mean_constraints not yet supported in random_ind")
  }

  # Generate AR params
  if(is.null(AR_constraints) && force_stability) {
    coefmat_pars <- as.vector(replicate(n=M, random_coefmats2(p=p, d=d, ar_scale=ar_scale)))
  } else {
    if(!is.null(AR_constraints)) stop("AR_constraints not yet supported in random_ind")
    scale_A <- ar_scale2*ifelse(is.null(constraints),
                                1 + log(2*mean(c((p - 0.2)^(1.25), d))),
                                1 + (sum(constraints)/(M*d^2))^0.85)
    coefmat_pars <- as.vector(replicate(n=M, expr=random_coefmats(d=d, how_many=p, scale=scale_A)))
  }

  # Generate covmat params
  covmat_pars <- as.vector(replicate(n=M, expr=random_covmat(d=d, omega_scale=omega_scale)))

  # Generate weight params
  weight_pars <- as.vector(replice(n=M, expr=random_weightpars(M=M, weight_function=weight_function,
                                                               AR_constraints=AR_constraints,
                                                               mean_constraints=mean_constraints)))

  # Generate distribution params (df etc)
  dist_pars <- random_distpars(cond_dist=cond_dist)

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

smart_ind <- function(p, M, d, params, weight_function=c("relative_dens", "logit"), cond_dist=c("Gaussian", "Student"),
                      AR_constraints=NULL, mean_constraints=NULL, accuracy=1, which_random=numeric(0), mu_scale, mu_scale2,
                      omega_scale, ar_scale=1, ar_scale2=1) {
  weight_function <- match.arg(weight_functions)
  cond_dist <- match.arg(cond_dist)
  if(!is.null(AR_constraints)) stop("AR_constraints not yet implemented to smart_ind!")
  if(!is.null(mean_constraints)) stop("mean_constraints not yet implemented to smart_ind!")
  if(cond_dist != "Gaussian") stop("Only Gaussian cond_dist is currently implemented to smart_ind!")
  scale_A <- ar_scale2*(1 + log(2*mean(c((p - 0.2)^(1.25), d))))
  # ? reform_constrained_pars here?
  weight_pars <- pick_weightpars(p=p, M=M, d=d, params=params, weight_function=weight_function, cond_dist=cond_dist)
  # ? dist_pars <- pick_distpars
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params)
  new_pars <- numeric(length(params))
  if(is.null(AR_constraints) || is.null(mean_constraints)) {
    all_phi0_and_A <- params[1:(d*M + M*p*d^2)] # all mu + A if called from GAfit
    new_pars[1:(d*M + M*p*d^2)] <- rnorm(n=length(all_phi0_and_A), mean=all_phi0_and_A, sd=pmax(0.2, abs(all_phi0_A)))
    for(m in 1:M) {
      if(any(which_random) == m) { # Not a smart regime
        new_pars[((m - 1)*d + 1):(m*d)] <- rnorm(d, mean=mu_scale, sd=mu_scale2) # Random mean
        if(runif(1) > 0.5) {
          new_pars[(M*d + (m - 1)*p*d^2 + 1):(M*d + m*p*d^2)] <- random_coefmats(d=d, how_many=p, scale=scale_A) # Without algorithm
        } else {
          new_pars[(M*d + (m - 1)*p*d^2 + 1):(M*d + m*p*d^2)] <- random_coefmats2(p=p, d=d, ar_scale=ar_scale) # Use the algorithm
        }
        new_pars[(M*d + M*p*d^2 + (m - 1)*d*(d + 1)/2 + 1):(M*d + M*p*d^2 + m*d*(d + 1)/2)] <- random_covmat(d=d, omega_scale=omega_scale)
      } else { # Smart_regime
        new_pars[(M*d + M*p*d^2 + (m - 1)*d*(d + 1)/2 + 1):(M*d + M*p*d^2 + m*d*(d + 1)/2)] <- smart_covmat(d=d, omega_scale=omega_scale,
                                                                                                            M=M, Omega=all_Omega[,m],
                                                                                                            accuracy=accuracy)
      }
    }
    new_pars[(M*d + M*p*d^2 + M*d*(d + 1)/2 + 1):(M*d + M*p*d^2 + M*d*(d + 1)/2 + 1)] <- smart_weightpars(M=M, weight_pars=weight_pars,
                                                                                                          weight_function=weight_function,
                                                                                                          accuracy=accuracy)
  } else {
    # AR or mean constraints employed

  }
  if(cond_dist == "Student") {
    # ADD DF PARAMETERS HERE
  }
  new_pars
}
