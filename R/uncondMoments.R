#' @title Calculate the dp-dimensional covariance matrix of p consecutive
#'  observations of a VAR process
#'
#' @description \code{VAR_pcovmat} calculate the dp-dimensional covariance matrix of p consecutive
#'  observations of a VAR process with the algorithm proposed by McElroy (2017).
#'
#' @inheritParams in_paramspace
#' @param all_Am \code{[d, d, p]} array containing the AR coefficient matrices
#' @param Omega_m the \eqn{(d\times d)} positive definite error term covariance matrix
#' @details
#'  Most of the code in this function is adapted from the one provided in the
#'  supplementary material of McElroy (2017). Reproduced under GNU General
#'  Public License, Copyright (2015) Tucker McElroy.
#' @return Returns the \eqn{(dp \times dp)} covariance matrix.
#' @references
#'  \itemize{
#'    \item McElroy T. 2017. Computation of vector ARMA autocovariances.
#'          \emph{Statistics and Probability Letters}, \strong{124}, 92-96.
#'  }
#' @keywords internal

VAR_pcovmat <- function(p, d, all_Am, Omega_m) {
  # all_Am = all_A[, , , m]
  # Omega_m = all_Omega[, , m]

  # The K commutation matrix
  Kcommut <- function(vect) matrix(t(matrix(vect, nrow=d, ncol=d)), ncol=1)
  Kmat <- apply(diag(d^2), MARGIN=1, FUN=Kcommut)

  # Step 1: vectorized error term covariance matrix for lag zero
  gamMAvec <- matrix(Omega_m, nrow=d^2, ncol=1)

  # Step 2: error term versus y_t covariance matrix for lag zero
  Amat <- array(0, dim=c(d^2, p + 1, d^2, 2*p + 1))
  Arow <-  matrix(nrow=d^2, ncol=d^2*(p + 1))
  start_inds <- seq(from=ncol(Arow) - d^2 + 1, to=1, by=-d^2)
  end_inds <- seq(from=ncol(Arow), to=d^2, by=-d^2)
  Arow[,start_inds[1]:end_inds[1]] <- diag(d^2)
  for(i1 in 1:p) {
    Arow[,start_inds[i1 + 1]:end_inds[i1 + 1]] <- -1*diag(d)%x%all_Am[, , i1]
  }
  for(i1 in 1:(p + 1)) {
    Amat[, i1, , i1:(i1 + p)] <- Arow
  }
  newA <- array(Amat[, 1:(p + 1), , 1:p], dim=c(d^2, p + 1, d^2, p))
  for(i1 in 1:(p + 1)) {
    for(i2 in 1:p) {
      newA[, i1, , i2] <- newA[, i1, , i2]%*%Kmat
    }
  }
  Amat <- cbind(matrix(Amat[, , , p + 1], nrow=d^2*(p + 1), ncol=d^2),
                matrix(Amat[, , , (p + 2):(2*p + 1)], nrow=d^2*(p + 1), ncol=d^2*(p)) + matrix(newA[, , , p:1], nrow=d^2*(p + 1), ncol=d^2*p))

  Bmat <- array(0, dim=c(d^2, 1, d^2, p + 1))
  Brow <-  matrix(nrow=d^2, ncol=d^2*(p + 1))
  start_inds <- rev(start_inds)
  end_inds <- rev(end_inds)
  Brow[,start_inds[1]:end_inds[1]] <- diag(d^2)
  for(i1 in 1:p) {
    Brow[,start_inds[i1 + 1]:end_inds[i1 + 1]] <- -1*diag(d)%x%all_Am[, , i1]
  }
  Bmat[, 1, , 1:(1 + p)] <- Brow
  gamMix <- solve(matrix(Bmat[, , , 1], nrow=d^2, ncol=d^2), gamMAvec)

  # Step 3: we pad out Gamma_WX(0) with zeros
  gamMixTemp <- c(gamMix, rep(0, times=p*d^2))

  # Step 4: compute the Gamma_XX(h) autocovariances for lags h=0,...,p
  gamVAR <- array(array(solve(Amat, gamMixTemp), dim=c(d, d, p + 1))[, , 1:p], dim=c(d, d, p))

  # After obtaining the autocovariances for lags h=0,...,p-1,
  # we construct the dp-dimensional covariance matrix for p consecutive
  # observations of a VAR process.

  # gamVAR contains lag=0 autocovariance in [, , 1], and lag=i in [, , i + 1].
  # Moreover, we use Gamma_Y(-h) = t(Gamma_Y(h)) and store the transposes
  # (as taking transpose multiple times uses more computation time):
  tgamVAR <- array(dim=c(d, d, p))
  for(i1 in 1:p) {
    tgamVAR[, , i1] <- t(gamVAR[, , i1])
  }

  # Finally, we fill in the covariance matrix for p consecutive observations
  # of the VAR process:
  Sigma_m <- matrix(nrow=d*p, ncol=d*p)
  start_inds <- seq(from=1, to=d*(p - 1) + 1, by=d)
  end_inds <- seq(from=d, to=d*p, by=d)
  for(i1 in 1:p) { # Go through row blocks
    for(i2 in 1:p) { # Go through column blocks
      # If i1 end_ind is larger than i2 end_ind, we consider blocks below block
      # diagonal (use transpose), and if i1 end_ind is smaller than i2 end_ind,
      # we consider blocks above block diagonal. If i1 end_ind == i2 end_ind,
      # we are at the block diagonal.
      if(end_inds[i1] > end_inds[i2]) { # Below diagonal, use transpose
        Sigma_m[start_inds[i1]:end_inds[i1], start_inds[i2]:end_inds[i2]] <- tgamVAR[, , abs(i1 - i2) + 1]
      } else {
        Sigma_m[start_inds[i1]:end_inds[i1], start_inds[i2]:end_inds[i2]] <- gamVAR[, , abs(i1 - i2) + 1]
      }
    }
  }
  Sigma_m
}


#' @title Calculate the dp-dimensional covariance matrices \eqn{\Sigma_{m,p}} in the transition weights
#'  with \code{weight_function="relative_dens"}
#'
#' @description \code{get_Sigmas} calculatesthe dp-dimensional covariance matrices \eqn{\Sigma_{m,p}} in
#'  the transition weights with \code{weight_function="relative_dens"} so that the algorithm proposed
#'  by McElroy (2017) employed whenever it reduces the computation time.
#'
#' @inheritParams in_paramspace
#' @inheritParams form_boldA
#' @param all_Omegas a \code{[d, d, M]} array containing the covariance matrix Omegas
#' @details
#'  Calculates the dp-dimensional covariance matrix using the formula (2.1.39) in Lütkepohl (2005) when
#'  \code{d*p < 12} and using the algorithm proposed by McElroy (2017) otherwise.
#'
#'  The code in the implementation of the McElroy's (2017) algorithm (in the function \code{VAR_pcovmat}) is
#'  adapted from the one provided in the supplementary material of McElroy (2017). Reproduced under GNU General
#'  Public License, Copyright (2015) Tucker McElroy.
#' @return Returns a \code{[dp, dp, M]} array containing the dp-dimensional covariance matrices for each regime.
#' @references
#'  \itemize{
#'    \item Lütkepohl H. 2005. New Introduction to Multiple Time Series Analysis,
#'            \emph{Springer}.
#'    \item McElroy T. 2017. Computation of vector ARMA autocovariances.
#'          \emph{Statistics and Probability Letters}, \strong{124}, 92-96.
#'  }
#' @keywords internal

get_Sigmas <- function(p, M, d, all_A, all_boldA, all_Omegas) {
  Sigmas <- array(NA, dim=c(d*p, d*p, M)) # Store the (dpxdp) covariance matrices
  if(d*p < 12) { # d*p < 12
    # Calculate the covariance matrices Sigma_{m,p} using the equation (2.1.39) in Lütkepohl (2005).
    I_dp2 <- diag(nrow=(d*p)^2)
    ZER_lower <- matrix(0, nrow=d*(p - 1), ncol=d*p)
    ZER_right <- matrix(0, nrow=d, ncol=d*(p - 1))
    for(m in 1:M) {
      kronmat <- I_dp2 - kronecker(all_boldA[, , m], all_boldA[, , m])
      sigma_epsm <- rbind(cbind(all_Omegas[, , m], ZER_right), ZER_lower)
      Sigma_m <- solve(kronmat, vec(sigma_epsm))
      Sigmas[, , m] <- Sigma_m
    }
  } else { # d*p >= 12
    # Calculate the covariance matrices Sigma_{m,p} using the algorithm proposed my McElroy (2017).
    for(m in 1:M) {
      Sigmas[, , m] <- VAR_pcovmat(p=p, d=d, all_Am=all_A[, , , m], Omega_m=all_Omegas[, , m])
    }
  }
  Sigmas
}


#' @title Calculate regime means \eqn{\mu_{m}}
#'
#' @description \code{get_regime_means} calculates regime means \eqn{\mu_{m} = (I - \sum A)^(-1))}
#'   from the given parameter vector.
#'
#' @inheritParams loglikelihood
#' @inheritParams stab_conds_satisfied
#' @return Returns a \eqn{(d\times M)} matrix containing regime mean \eqn{\mu_{m}} in the m:th column, \eqn{m=1,..,M}.
#' @section Warning:
#'  No argument checks!
#' @inherit stab_conds_satisfied references
#' @keywords internal

get_regime_means <- function(p, M, d, params,
                             weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                             weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t"),
                             parametrization=c("intercept", "mean"),
                             identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity"),
                             AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  identification <- match.arg(identification)
  weightfun_pars <- check_weightfun_pars(p=p, M=M, d=d, weight_function=weight_function,
                                         weightfun_pars=weightfun_pars, cond_dist=cond_dist)
  check_constraints(p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    parametrization=parametrization, identification=identification,
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=B_constraints)
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, weight_function=weight_function,
                                    weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                                    identification=identification, AR_constraints=AR_constraints,
                                    mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                    B_constraints=B_constraints)

  if(parametrization == "intercept") {
    params <- change_parametrization(p=p, M=M, d=d, params=params, weight_function=weight_function,
                                     weightfun_pars=weightfun_pars, identification=identification,
                                     AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL,
                                     B_constraints=NULL, change_to="mean")
  }
  pick_phi0(M=M, d=d, params=params)
}



#' @title Calculate regimewise autocovariance matrices
#'
#' @description \code{get_regime_autocovs} calculates the regimewise autocovariance matrices \eqn{\Gamma_{m}(j)}
#'  \eqn{j=0,1,...,p}
#'
#' @inheritParams loglikelihood
#' @inheritParams reform_constrained_pars
#' @return Returns an \eqn{(d \times d \times p+1 \times M)} array containing the first p regimewise autocovariance matrices.
#'   The subset \code{[, , j, m]} contains the j-1:th lag autocovariance matrix of the m:th regime.
#' @references
#'   \itemize{
#'    \item Lütkepohl H. 2005. New Introduction to Multiple Time Series Analysis,
#'            \emph{Springer}.
#'   }
#' @keywords internal

get_regime_autocovs <- function(p, M, d, params,
                                weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold", "exogenous"),
                                weightfun_pars=NULL, cond_dist=c("Gaussian", "Student", "ind_Student", "ind_skewed_t"),
                                identification=c("reduced_form", "recursive", "heteroskedasticity", "non-Gaussianity"),
                                AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL) {
  # Match args
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  identification <- match.arg(identification)

  weightfun_pars <- check_weightfun_pars(p=p, M=M, d=d, weight_function=weight_function,
                                         weightfun_pars=weightfun_pars, cond_dist=cond_dist)

  # Collect the required parameter values
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, weight_function=weight_function,
                                    cond_dist=cond_dist, identification=identification,
                                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                                    weight_constraints=weight_constraints, B_constraints=B_constraints,
                                    weightfun_pars=weightfun_pars)

  all_A <- pick_allA(p=p, M=M, d=d, params=params) # [d, d, p, M]
  all_Omegas <- pick_Omegas(p=p, M=M, d=d, params=params, cond_dist=cond_dist, identification=identification) # [d, d, M]
  if(cond_dist == "ind_Student" || cond_dist == "ind_skewed_t" || identification == "non-Gaussiniaty") {
    # Calculate the covariance matrices of the regimes from the impact matrices
    for(m in 1:M) {
      all_Omegas[, , m] <- tcrossprod(all_Omegas[, , m])
    }
  }

  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)

  # Calculate teh regimewise autocovarinces
  I_dp2 <- diag(nrow=(d*p)^2)
  ZER_lower <- matrix(0, nrow=d*(p-1), ncol=d*p)
  ZER_right <- matrix(0, nrow=d, ncol=d*(p - 1))
  # For each m=1,..,M, store the (dxd) covariance matrices Gamma_{y,m}(0),...,Gamma{y,m}(p-1),,Gamma{y,m}(p):
  all_Gammas <- array(NA, dim=c(d, d, p + 1, M))
  for(m in 1:M) {
    # Calculate the (dpxdp) Gamma_{Y,m}(0) covariance matrix (Lütkepohl 2005, eq. (2.1.39))
    kronmat <- I_dp2 - kronecker(all_boldA[, , m], all_boldA[, , m])
    sigma_epsm <- rbind(cbind(all_Omegas[, , m], ZER_right), ZER_lower)
    Gamma_m <- matrix(solve(kronmat, vec(sigma_epsm)), nrow=d*p, ncol=d*p, byrow=FALSE)

    # Obtain the Gamma_{y,m}(0),...,Gamma_{y,m}(p-1) covariance matrices from Gamma_{Y,m}(0)
    all_Gammas[, , , m] <- c(as.vector(Gamma_m[1:d,]), rep(NA, d*d))

    # Calculate the Gamma{y,m}(p) recursively from Gamma_{y,m}(0),...,Gamma_{y,m}(p-1) (Lütkepohl 2005, eq. (2.1.37))
    all_Gammas[, , p + 1, m] <- rowSums(vapply(1:p, function(i1) all_A[, ,i1 , m]%*%all_Gammas[, , p + 1 - i1, m], numeric(d*d)))
  }
  all_Gammas
}



#' @title Calculate the unconditional means, variances, the first p autocovariances, and the first p autocorrelations
#'  of the regimes of the model.
#'
#' @description \code{uncond_moments} calculates the unconditional means, variances, the first p autocovariances,
#'  and the first p autocorrelations of the regimes of the model.
#'
#' @inheritParams get_boldA_eigens
#' @return Returns a list with three components:
#'   \describe{
#'     \item{\code{$regime_means}}{a \eqn{M \times d} matrix vector containing the unconditional mean of the regime
#'           \eqn{m} in the \eqn{m}th column.}
#'     \item{\code{$regime_vars}}{a \eqn{M \times d} matrix vector containing the unconditional marginal variances
#'           of the regime \eqn{m} in the \eqn{m}th column.}
#'     \item{\code{$regime_autocovs}}{an \eqn{(d x d x p+1, M)} array containing the lag 0,1,...,p autocovariances of the process.
#'           The subset \code{[, , j, m]} contains the lag \code{j-1} autocovariance matrix (lag zero for the variance) for
#'           the regime \eqn{m}.}
#'     \item{\code{$regime_autocors}}{the autocovariance matrices scaled to autocorrelation matrices.}
#'   }
#' @inherit get_regime_autocovs references
#' @examples
#' # Two-variate Gaussian STVAR p=1, M=2 model with the weighted relative stationary
#' # densities of the regimes as the transition weight function:
#' theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863,
#' -0.3124, 0.505251, 0.298483, 0.030096, -0.176925, 0.838898, 0.310863, 0.007512,
#' 0.018244, 0.949533, -0.016941, 0.121403, 0.573269)
#' mod122 <- STVAR(data=gdpdef, p=1, M=2, params=theta_122relg, weight_function="relative_dens")
#'
#' # Calculate the unconditional moments of model:
#' tmp122 <- uncond_moments(mod122)
#'
#' # Print the various unconditional moments calculated:
#' tmp122$regime_means[,1] # Unconditional means of the first regime
#' tmp122$regime_means[,2] # Unconditional means of the second regime
#' tmp122$regime_vars[,1] # Unconditional variances of the first regime
#' tmp122$regime_vars[,2] # Unconditional variances of the second regime
#' tmp122$regime_autocovs[, , , 1] # a.cov. matrices of the first regime
#' tmp122$regime_autocovs[, , , 2] # a.cov. matrices of the second regime
#' tmp122$regime_autocors[, , , 1] # a.cor. matrices of the first regime
#' tmp122$regime_autocors[, , , 2] # a.cor. matrices of the second regime
#'
#' # A two-variate linear Gaussian VAR p=1 model:
#' theta_112 <- c(0.649526, 0.066507, 0.288526, 0.021767, -0.144024, 0.897103,
#'  0.601786, -0.002945, 0.067224)
#' mod112 <- STVAR(data=gdpdef, p=1, M=1, params=theta_112)
#'
#' # Calculate the unconditional moments of model:
#' tmp112 <- uncond_moments(mod112)
#'
#' # Print the various unconditional moments calculated:
#' tmp112$regime_means # Unconditional means
#' tmp112$regime_vars # Unconditional variances
#' tmp112$regime_autocovs # Unconditional autocovariance matrices
#' tmp112$regime_autocovs[, , 1, 1] # a.cov. matrix of lag zero (of the first regime)
#' tmp112$regime_autocovs[, , 2, 1] # a.cov. matrix of lag one (of the first regime)
#' tmp112$regime_autocors # Unconditional autocorrelation matrices
#' @export

uncond_moments <- function(stvar) {
  check_stvar(stvar)
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  params <- stvar$params
  weigth_function <- stvar$model$weight_function
  weightfun_pars <- stvar$model$weightfun_pars
  cond_dist <- stvar$model$cond_dist
  parametrization <- stvar$model$parametrization
  identification <- stvar$model$identification
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  weight_constraints <- stvar$model$weight_constraints
  B_constraints <- stvar$model$B_constraints

  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, weight_function=weigth_function,
                                    weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                                    identification=identification, AR_constraints=AR_constraints,
                                    mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                    B_constraints=B_constraints)

  if(stvar$allow_unstab) { # Check that the stability condition is satisfied
    if(!stab_conds_satisfied(p=p, M=M, d=d, params=params)) {
      warnings("Cannot calculate unconditonal moments: the stability condition is not satisfied for all regimes.")
      return(NULL)
    }
  }

  reg_means <- get_regime_means(p=p, M=M, d=d, params=params, weight_function=weigth_function,
                                weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                                parametrization=parametrization, identification=identification,
                                AR_constraints=NULL, mean_constraints=NULL,
                                weight_constraints=NULL, B_constraints=NULL)

  reg_autocovs <- get_regime_autocovs(p=p, M=M, d=d, params=params, weight_function=weigth_function,
                                      weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                                      identification=identification, AR_constraints=NULL,
                                      mean_constraints=NULL, weight_constraints=NULL,
                                      B_constraints=NULL)

  # Obtain the unconditional variances from the diagonal of the first slice of reg_autocovs[, , , m]
  reg_vars <- vapply(1:M, function(m) diag(reg_autocovs[, , 1, m]), numeric(d))

  # Obtain the regimewise autocorrelation matrices from the regimewise autocovariance matrices
  reg_autocors <- array(NA, dim=c(d, d, p + 1, M))
  for(m in 1:M) {
    for(i1 in 1:(p + 1)) {
      for(i2 in 1:d) {
        for(i3 in 1:d) {
          reg_autocors[i2, i3, i1, m] <- reg_autocovs[i2, i3, i1, m]/sqrt(reg_vars[i2, m]*reg_vars[i3, m])
        }
      }
    }
  }

  # Return the results
  list(regime_means=reg_means,
       regime_vars=reg_vars,
       regime_autocovs=reg_autocovs,
       regime_autocors=reg_autocors)
}
