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
#' @param all_Omega a \code{[d, d, M]} array containing the covariance matrix Omegas
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
#' @return Returns a \eqn{(dxM)} matrix containing regime mean \eqn{\mu_{m}} in the m:th column, \eqn{m=1,..,M}.
#' @section Warning:
#'  No argument checks!
#' @inherit stab_conds_satisfied references
#' @keywords internal

get_regime_means <- function(p, M, d, params, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold"),
                             weightfun_pars=NULL, cond_dist=c("Gaussian", "Student"), parametrization=c("intercept", "mean"),
                             identification=c("reduced_form", "recursive", "heteroskedasticity"),
                             AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  identification <- match.arg(identification)
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                         cond_dist=cond_dist)
  check_constraints(p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    parametrization=parametrization, identification=identification,
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=B_constraints)
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params, weight_function=weight_function,
                                    weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                                    identification=identification, AR_constraints=AR_constraints,
                                    mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                    B_constraints=B_constraints)
  if(!is.null(B_constraints)) {
    stop("B_constraints are not yet implemented to get_regime_means")
  }
  if(identification != "reduced_form") stop("Structural models are not yet implemented to get_regime_means")

  if(parametrization == "intercept") {
    params <- change_parametrization(p=p, M=M, d=d, params=params, weight_function=weight_function,
                                     weightfun_pars=weightfun_pars, AR_constraints=NULL,
                                     mean_constraints=NULL, weight_constraints=NULL,
                                     B_constraints=NULL, change_to="mean")
  }
  pick_phi0(M=M, d=d, params=params)
}
