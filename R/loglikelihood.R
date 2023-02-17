
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
#'     \item{\eqn{\alpha} contains the transition weights parameters}
#'     \item{\eqn{\nu > 2} is the degrees of freedom parameter that is included only if \code{cond_dist="Student"}.}
#'   }
#'   \describe{
#'     \item{For models with \code{weight_function="relative_dens"}:}{\eqn{\alpha = (\alpha_1,...,\alpha_{M-1})}
#'           \eqn{(M - 1 \times 1)}, where \eqn{\alpha_m} \eqn{(1\times 1), m=1,...,M-1} are the transition weight parameters.}
#'     \item{For models with \code{weight_function="logit"}:}{\eqn{\alpha = (\gamma_1,...,\gamma_M)} \eqn{((M-1)k\times 1)},
#'           where \eqn{\gamma_m} \eqn{(k\times 1), m=1,...,M-1} contains the logit-regression coefficients of the \eqn{m}th regime.}
#'   }
#'   Above, \eqn{\phi_{m,0}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}th coefficient matrix of the \eqn{m}th
#'   mixture component, and \eqn{\Omega_{m}} denotes the error term covariance matrix of the \eqn{m}:th mixture component.
#'   If \code{parametrization=="mean"}, just replace each \eqn{\phi_{m,0}} with regimewise mean \eqn{\mu_{m}}.
#'   \eqn{vec()} is vectorization operator that stacks columns of a given matrix into a vector. \eqn{vech()} stacks columns
#'   of a given matrix from the principal diagonal downwards (including elements on the diagonal) into a vector.
#' @param weight_function what type of transition weights should be used?
#' @param cond_dist specifies the conditional distribution of the model as \code{"Gaussian"} or \eqn{"Student"}.
#' @param parametrization \code{"intercept"} or \code{"mean"} determining whether the model is parametrized with intercept
#'   parameters \eqn{\phi_{m,0}} or regime means \eqn{\mu_{m}}, m=1,...,M.
#' @param identification is it reduced form model or an identified structural model; if the latter, how is it identified?
#' @param to_return should the returned object be the log-likelihood, which is the default, or something else?
#'   See the section "Return" for all the options.
#' @param check_params should it be checked that the parameter vector satisfies the model assumptions? Can be skipped to save
#'   computation time if it does for sure.
#' @param minval the value that will be returned if the parameter vector does not lie in the parameter space
#'   (excluding the identification condition).
#' @param stat_tol numerical tolerance for stability of condition of the regimes: if the "bold A" matrix of any regime
#'   has eigenvalues larger that \code{1 - stat_tol} the parameter is considered to be outside the parameter space.
#'   Note that if tolerance is too small, numerical evaluation of the log-likelihood might fail and cause error.
#' @param posdef_tol numerical tolerance for positive definiteness of the error term covariance matrices: if
#'   the error term covariance matrix of any regime has eigenvalues smaller than this, the parameter is considered
#'   to be outside the parameter space. Note that if the tolerance is too small, numerical evaluation of the
#'   log-likelihood might fail and cause error.
#' @param df_tol the parameter vector is considered to be outside the parameter space if the degrees of
#'   freedom parameters is not larger than \code{2 + df_tol} (applies only if \code{cond_dist="Student"}).
#' @details FILL IN
#' @return
#'   \describe{
#'     \item{If \code{to_return="loglik"}:}{the log-likelihood of the specified model.}
#'     \item{If \code{to_return=="tw"}:}{a size \eqn{((n_obs-p)\times M)} matrix containing the transition weights: for m:th component
#'       in m:th column.}
#'     \item{If \code{to_return=="terms"}:}{a size \eqn{((n_obs-p)\times 1)} numeric vector containing the terms \eqn{l_{t}}.}
#'     \item{If \code{to_return=="regime_cmeans"}:}{an \code{[n_obs-p, d, M]} array containing the regimewise conditional means.}
#'     \item{If \code{to_return=="total_cmeans"}:}{a \code{[n_obs-p, d]} matrix containing the conditional means of the process.}
#'     \item{If \code{to_return=="total_ccovs"}:}{an \code{[d, d, n_obs-p]} array containing the conditional covariance matrices of the process.}
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

loglikelihood <- function(data, p, M, params, weight_function=c("relative_dens", "logit"), cond_dist=c("Gaussian", "Student"),
                          parametrization=c("intercept", "mean"),
                          identification=c("reduced_form", "recursive", "heteroskedasticity"),
                          AR_constraints=NULL, mean_constraints=NULL, B_constraints=NULL,
                          to_return=c("loglik", "tw", "terms", "regime_cmeans", "total_cmeans", "total_ccovs"),
                          check_params=TRUE, minval=NULL, stab_tol=1e-3, posdef_tol=1e-8, df_tol=1e-8) {
  # Match args
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  identification <- match.arg(identification)
  to_return <- match.arg(to_return)
  if(identification != "reduced_form") stop("Only reduced form models are currently supported")
  if(!is.null(AR_constraints) || !is.null(mean_constraints) || !is.null(B_constraints)) stop("Constrained models are not
                                                                                             currently supported")

  # Compute some required statistics
  epsilon <- round(log(.Machine$double.xmin) + 10) # Logarithm of the smallest value that can be handled normally
  d <- ncol(data)
  n_obs <- nrow(data)
  T_obs <- n_obs - p

  # Collect the parameter values
  # First remove all constraints, if any; also switch to reduced form parameter vector; TO BE IMPLEMENTED

  # Pick params
  if(parametrization == "intercept") { # [d, M]
    all_phi0 <- pick_phi0(M=M, d=d, params=params)
  } else {
    all_mu <- pick_phi0(M=M, d=d, params=params) # mean parameters instead of intercepts
  }
  all_A <- pick_allA(p=p, M=M, d=d, params=params) # [d, d, p, M]
  all_Omegas <- pick_Omegas(p=p, M=M, d=d, params=params) # [d, d, M]
  weightpars <- pick_weightpars(p=p, M=M, d=d, params=params, weight_function=weight_function, cond_dist=cond_dist)
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  df <- numeric(0) # FILL IN WHEN STUDENT IS IMPLEMENTED

  # Check that the parameter vector lies in the parameter space
  if(check_params) {
    if(!in_paramspace(p=p, M=M, d=d, weight_function=weight_function, cond_dist=cond_dist,
                      all_boldA=all_boldA, all_Omegas=all_Omegas, weightpars=weightpars, df=df,
                      stab_tol=stab_tol, posdef_tol=posdef_tol, df_tol=df_tol)) {
      return(minval)
    }
  }

  # i:th row denotes the vector \bold{y_{i-1}} = (y_{i-1},...,y_{i-p}) (dpx1),
  # assuming the observed data is y_{-p+1},...,y_0,y_1,...,y_{T}
  Y <- reform_data(data, p)
  Y2 <- Y[1:T_obs,] # Last row removed; not needed when calculating something based on lagged observations

  # NOTE! IF parametrization == "intercept" && weight_function == "logit" WE DONT NEED TO CALCULATE all_mu
  # THAT MIGHT SAVE COMPUTATION TIME IN GRADIENT BASED ESTIMATION!
  # Calculate unconditional regime-specific expected values (column per component) or phi0-parameters if using mean-parametrization
  Id <- diag(nrow=d)
  if(parametrization == "intercept") {
    all_mu <- vapply(1:M, function(m) solve(Id - rowSums(all_A[, , , m, drop=FALSE], dims=2), all_phi0[,m]), numeric(d)) # rowSums: sum over dims+1=3
  } else {
    all_phi0 <- vapply(1:M, function(m) (Id - rowSums(all_A[, , , m, drop=FALSE], dims=2))%*%all_mu[,m], numeric(d))
  }

  # Calculate the transition weights [T_obs, M] with [t,m] indexing (nothing for the initial values here)
  alpha_mt <- get_alpha_mt(data=data, Y2=Y2, p=p, M=M, d=d, weight_function=weight_function, all_A=all_A, all_boldA=all_boldA,
                           all_Omegas=all_Omegas, weightpars=weightpars, all_mu=all_mu, epsilon=epsilon)
  if(to_return == "tw") {
    return(alpha_mt)
  }

  # Calculate the conditional means mu_{m,t}
  # The dimensions of mu_mt will be: [t, p, m]
  all_A2 <- array(all_A, dim=c(d, d*p, M)) # cbind coefficient matrices of each component: m:th component is obtained at [, , m]
  mu_mt <- array(vapply(1:M, function(m) t(all_phi0[, m] + tcrossprod(all_A2[, , m], Y2)), numeric(d*T_obs)), dim=c(T_obs, d, M)) # [, , m]

  # mu_yt VAIHDA NOPEAMPAAN KUN SAADAAN JOTAIN LUKUJA!!!
  mu_yt <- matrix(nrow=T_obs, ncol=d) # [t, d]
  for(i1 in 1:T_obs) {
    for(i2 in 1:d) {
      mu_yt[i1, i2] <- sum(alpha_mt[i1,]*mu_mt[i1,i2,])
    }
  }

  # Return conditional moments if those were to be returned
  if(to_return == "regime_cmeans") { # Regime-specific conditional menas
    return(mu_mt)
  } else if(to_return == "total_cmeans") { # Cond means of the process: weighted sum of regime-specific conditional means
    return(matrix(rowSums(vapply(1:M, function(m) alpha_mt[,m]*mu_mt[, , m], numeric(d*T_obs))), nrow=T_obs, ncol=d, byrow=FALSE))
  } else { # Cond covariance matrices of the process: weighted sum of regime-specific cond cov mats
    all_covmats <- array(rowSums(vapply(1:M, function(m) rep(alpha_mt[, m], each=d*d)*as.vector(all_Omegas[, , m]),
                                        numeric(d*d*T_obs))), dim=c(d, d, T_obs))
    if(to_return == "total_ccovs") {
      return(all_covmats)
    }
  }

  # Calculate the conditional log-likelihood
  dat <- data[(p + 1):n_obs,] # Initial values are not used here (conditional means and variances are already calculated)
  mvd_vals <- matrix(nrow=T_obs, ncol=M)
  if(cond_dist == "Gaussian") { # Gaussian conditiona distribution
    all_lt <- numeric(T_obs)
    for(i1 in 1:T_obs) {
      # Calculate the l_t multinormal density for each observation
      # # +p in data because the first row is for the initial values
      all_lt[i1] <- -0.5*d*log(2*pi) - 0.5*log(det(all_covmats[, , i1])) - 0.5*crossprod(data[i1+p,] - mu_yt[i1,],
                                                                                         solve(all_covmats[, , i1],
                                                                                               data[i1+p,] - mu_yt[i1,]))
    }

    # TÄÄLLÄ:
    # https://gallery.rcpp.org/articles/dmvnorm_arma/
    # on tämän: https://gallery.rcpp.org/articles/dmvnorm_arma/
    # koodit saatavilla ja käytettävissä. Ja siellä sanotaan:
    # "For instance, the inverse Cholesky decomposition can be put inside the main loop, if varying covariance matrices are necessary."
    # Eli tuo varmaan on syytä ottaa käyttöön sopivalla modifikaatiolla! Sen voi varmaan ottaa käyttöön myös alpha_mt:ssä?
    # Kun compile-codea löytyy joka tapauksessa käyttämään.
    # MUTTA ensin hidas versio, jolla saadaan jotkut luvut, jotta voidaan testata koodin rikkoitumista yksikkötesteillä.
    # Jätä hidas tapa laskea kommentteihin, jotta voi myöhemmin CRAN-tiimille perustellessa laskea eroja estimointiajoissa
    # ja sillä perustella compiled koodin käyttöä.

    # mvnfast taitaa käyttää samaa kovarianssimatriisia kaikilla t. Eli joutuunee laskemaan käsin?
    # Notable computational effort is required to invert all the covariance matrcies for all the
    # total_ccovs antaa ne kovarianssimatriisit
    # RccpArmadillo? https://gallery.rcpp.org/articles/dmvnorm_arma/
    #backsolve (chol2inv nopeampi)

    # tmp <- matrix(round(rnorm(16), 2), nrow=4)
    # testcovmat <- crossprod(tmp)
    # tmp2 <- backsolve(chol(testcovmat), x=diag(4)) # noin 10 mikrosekunttia; tmp2%*%t(tmp2) == testcovmat (jälkimmäinen noin 500 nanosekunttia)
    # sum(log(diag(tmp2))) # = -0.5*log(det); noin 3 mikrosekunttia - diag on jostain syystä hidas; siksi hidas koska ei voi laskea kaikille t saman
    # samanaikaisesti; mutta voi testata miten kun tuon quadratic formin saa laskettua nopeammin niin toimiiko


    #  microbenchmark::microbenchmark(chol2inv(chol(testcovmat))) # Noin 5.5 mikrosekunttia eli tämä nopeampi
    #  microbenchmark::microbenchmark(det(testcovmat)) # Noin 5.5 mikrosekunttia eli chol2inv + det hieman nopeampi?; menee noin 12 microsekunntia yht
    # --> toista 300 havainnolle niin tulee 1.65 millisekunttia; mikähän on koko loglikin evualuoinnin aika?
    # Mutta tuo ei laske determinanttia, joten backsolvella menee varmaan yht nopeammin?
    # HUOM hidastuu huomattavasti jos aikasarjan pituus kasvaa! Looppi rccp:llä? Koska joutuu mennä kaikki t:t läpi.. Miksi sum vie myös sekunnin?

    # Voi myös koittaa loopata mwnfastia? Mutta se on varmaan hitain tapa..

    # 12*300*0.001 = 3.6 millisekunttia eli ehkä 5 millisekunttia tms kuluu kokonaisuudessaan?
    # Seuraavaksi varmaan tee loppuun että saat jotkut arvot ja ala testailemaan aikoja ja nopeutuksia? Muista testata eri d ja T arvoilla

    # gmvarkitissä GMVAR(3,2) d=4 loglik 300 obsilla vie noin 2.5 millisekunttia (läppärillä).

    # mvnfast::dmvn käyttää .Call("dmvtCpp", X_ = X, mu_ = mu, sigma_ = sigma, df = -1,  log_ = log, ncores_ = ncores, isChol_ = isChol)
    # Voiko tuota dmvtCpp:tä käyttää suoraan myös sstvarssista? Tämä on mvnfast/src paketin oma funktio
    # https://github.com/mfasiolo/mvnfast/blob/master/src/dmvtCpp.cpp
    # Tuosta voi GNU lisenssin perusteella muokata omaan tarpeeseen ja testata nopeuttaako.

    # KUN OLET TESTANNUT JONKIN MENETELMÄT, ETTÄ SAA LOGLIKARVOJA; KOKEILE VOIKO MYÖS ALPHA_MT NOPEUTTAA KÄYTÄTMÄLLÄ:
    #rooti <- backsolve(chol(Sigmas[, , M]), x=diag(d))
    #quads <- colSums((crossprod(rooti,(t(Y2) - all_mu[,m])))^2) # Saako rep(all_mu[,m], times=p) implementoilla nopeammaksi?
    #exp(-(d/2)*log(2*pi) + sum(log(diag(rooti))) - .5*quads) # =  sum(log(diag(rooti))) = det

    #mvd_vals <- vapply(1:M, function(m) mvnfast::dmvn(X=dat - mu_mt[, , m], mu=rep(0, times=d), sigma=all_Omega[, , m], log=FALSE,
    #                   ncores=1, isChol=FALSE), numeric(T_obs))
  } else if(cond_dist == "Student") {
    # Entä RCCP toimiiko Studentille? Pitää vain koodata
    stop("Student's t cond_dist is not implemented yet!")
  }
  if(to_return == "terms") {
    return(all_lt)
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
#' @details Note that we index the time series as \eqn{-p+1,...,0,1,...,T}.
#' @return Returns the mixing weights a \eqn{(T x M)} matrix, so that the t:th row is for the time point t
#'   and m:th column is for the regime m.
#' @inherit in_paramspace references
#' @keywords internal

get_alpha_mt <- function(data, Y2, p, M, d, weight_function, all_A, all_boldA, all_Omegas, weightpars, all_mu, epsilon) {
  T_obs <- nrow(data) - p
  if(M == 1) {
    return(as.matrix(rep(1, times=T_obs)))
  }
  if(weight_function == "relative_dens") {
    # Calculate the covariance matrices Sigma_{m,p} (Lutkepohl 2005, eq. (2.1.39) or the algorithm proposed by McElroy 2017)
    Sigmas <- get_Sigmas(p=p, M=M, d=d, all_A=all_A, all_boldA=all_boldA, all_Omegas=all_Omegas) # Store the (dpxdp) covariance matrices
    chol_Sigmas <- array(dim=c(d*p, d*p, M))
    for(m in 1:M) {
      chol_Sigmas[, , m] <- chol(Sigmas[, , m]) # Take Cholesky here to avoid unnecessary warnings from mvnfast::dmvn
    }

    # Calculate the dp-dimensional multinormal densities in logarithm with the package mvnfast:
    # i:th row for index i-1 etc, m:th column for m:th component.
    # We calculate in logarithm because the non-log values may be too close to zero for machine accuracy (if they are too close to zero
    # for all regimes and computer handles them as zero, we would divide by zero when calculating the transition weights).
    log_mvdvalues <- vapply(1:M, function(m) mvnfast::dmvn(X=Y2, mu=rep(all_mu[,m], p), sigma=chol_Sigmas[, , m],
                                                                   log=TRUE, ncores=1, isChol=TRUE),
                            numeric(T_obs)) # [T_obs, M] removes the period T+1 weights



    # Calculate the transition weights based on the log-multivariate density values
    if(!is.matrix(log_mvdvalues)) log_mvdvalues <- t(as.matrix(log_mvdvalues)) # Only one time point but multiple regimes
    log_mvdvalues_orig <- log_mvdvalues
    small_logmvns <- log_mvdvalues < epsilon
    if(any(small_logmvns)) {
      # If too small or large non-log-density values are present (i.e., that would yield -Inf or Inf),
      # we replace them with ones that are not too small or large but imply the same mixing weights
      # up to negligible numerical tolerance (tested in gmvarkit).
      which_change <- rowSums(small_logmvns) > 0 # Which rows contain too small  values
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
  } else if(weight_function == "logit") {
    stop("Logit weight function is not yet implemented to get_alpha_mt")
  }
  alpha_mt
}
