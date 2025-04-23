#' @title Compute historical decompositions for structural STVAR models.
#'
#' @description \code{hist_decomp} compute historical decompositions for structural STVAR models.
#'
#' @inheritParams simulate.stvar
#' @details The historical decomposition quantifies the cumulative effects the shocks to the movements of
#'   the variables (see, e.g., Kilian and Lütkepohl, 2017, Section~4.3) The historical decompositions are
#'   computed as described in Wong (2018). Note that due to the effect of the "initial conditions" and the
#'   "steady state component", which are not attributed to the shocks, the cumulative effects of the shocks
#'   do not sum to the observed time series.
#' @return Returns a class \code{'histdecomp'} list with the following elements:
#'   \describe{
#'     \item{init_cond_comp}{A matrix of size \eqn{(T \times d)} containing the contributions of the initial conditions to the movements of
#'      the variables at each time point; the element \code{t, i} giving the contribution at the time \code{t} on the variable \code{i}.}
#'     \item{steady_state_comp}{A matrix of size \eqn{(T \times d)} containing the contributions of the steady state component to the movements of
#'      the variables at each time point; the element \code{t, i} giving the contribution at the time \code{t} on the variable \code{i}.}
#'     \item{shock_comp}{A matrix of size \eqn{(T \times d)} containing the contributions of the shocks to the movements of the variables at each
#'      time point; the element \code{t, i} giving the contribution at the time \code{t} on the variable \code{i}.}
#'     \item{contributions_of_shocks}{A 3D array of size \eqn{(T \times d \times d)} containing the cumulative contributions of the shocks to the
#'      movements of the variables at each time point; the element \code{t, i1, i2} giving the contribution of the shock \code{i1} to the variable
#'      \code{i2} at the time \code{t}.}
#'     \item{stvar}{The original STVAR model object.}
#'  }
#' @seealso \code{\link{GIRF}}, \code{\link{GFEVD}}, \code{\link{linear_IRF}}, \code{\link{fitSSTVAR}}
#'  \itemize{
#'    \item Kilian L., Lütkepohl H. 2017. Structural Vector Autoregressive Analysis. 1st edition.
#'      \emph{Cambridge University Press}, Cambridge.
#'    \item Wong H. 2018. Historical decomposition for nonlinear vector autoregressive models.
#'      \emph{CAMA Working Paper No. 62/2017, available as SSRN:3057759}
#'  }
#' @examples
#'  ## FILL IN
#' @export

hist_decomp <- function(stvar) {
  check_stvar(stvar)
  if(is.null(stvar$data)) stop("The model needs to contain data")
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  T_obs <- nrow(stvar$data) - p
  cond_dist <- stvar$model$cond_dist
  parametrization <- stvar$model$parametrization
  identification <- stvar$model$identification
  weight_function <- stvar$model$weight_function
  weightfun_pars <- stvar$model$weightfun_pars

  ## Obtain the parameter values in the "non constrained form", and also switch to reduced_form parameter vector:
  params <- reform_constrained_pars(p=p, M=M, d=d, params=stvar$params, weight_function=weight_function, cond_dist=cond_dist,
                                    identification=identification, AR_constraints=stvar$model$AR_constraints,
                                    mean_constraints=stvar$model$mean_constraints, weight_constraints=stvar$model$weight_constraints,
                                    B_constraints=stvar$model$B_constraints, other_constraints=NULL, weightfun_pars=weightfun_pars)

  ## Pick params
  all_A <- pick_allA(p=p, M=M, d=d, params=params) # [d, d, p, M]
  all_Omegas <- pick_Omegas(p=p, M=M, d=d, params=params, cond_dist=cond_dist, identification=identification) # [d, d, M], B_m for ind_Stud and ind_skewed_t
  weightpars <- pick_weightpars(p=p, M=M, d=d, params=params, weight_function=weight_function, cond_dist=cond_dist, weightfun_pars=weightfun_pars)
  all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A) # The bold A matrices of the regimes to be used later
  distpars <- pick_distpars(d=d, params=params, cond_dist=cond_dist)
  if(parametrization == "intercept") { # [d, M]
    all_phi0 <- pick_phi0(M=M, d=d, params=params)
  } else {
    all_phi0 <- vapply(1:M, function(m) (Id - rowSums(all_A[, , , m, drop=FALSE], dims=2))%*%all_mu[,m], numeric(d))
  }

  ## Obtain the structural shocks
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

  ## Create the (dp x d) matrix H:
  H <- rbind(diag(d), matrix(0, nrow=(d*(p-1)), ncol=d))

  ## Create the [d, d, T] array of impact matrices B_{y,t}:
  if(cond_dist %in% c("ind_Student", "ind_skewed_t")) { # Parametrization directly via the impact matrix
    B_yt <- get_Bt_Cpp(all_Omegas=all_Omegas, alpha_mt=alpha_mt)
  } else if(identification == "heteroskedasticity") {
    W <- pick_W(p=p, M=M, d=d, params=params, identification=identification) # [d, d, M])
    lambdas <- matrix(pick_lambdas(p=p, M=M, d=d, params=params, identification=identification), nrow=d, ncol=M - 1) # [d, M - 1]
    if(M > 1) lambdas <- cbind(1, matrix(lambdas, nrow=d, ncol=M - 1)) # First column is column of ones for the first regime
    B_yt <- array(0, dim=c(d, d, T_obs))
    for(i1 in 1:T_obs) {
      if(M == 1) {
        B_t <- W
      } else {
        tmp <- array(dim=c(d, d, M))
        for(m in 1:M) {
          tmp[, , m] <- alpha_mt[i1, m]*diag(lambdas[, m])
        }
        B_yt[, , i1] <- W%*%sqrt(apply(tmp, MARGIN=1:2, FUN=sum))
      }
    }
  } else { # Non-identified reduced form or recursive identification
    B_yt <- array(0, dim=c(d, d, T_obs))
    for(i1 in 1:T_obs) {
      B_yt[, , i1] <- t(chol(stvar$total_ccovs[, , i1])) # Lower triangular Cholesky decomposition of Omega_t
    }
  }

  ## Create the [dp, dp, T] array of time-varying "bold A" (companion form) matrices bold A_{y,t},
  ## AND
  ## create also the [T_obs, d] matrix of the time-varying intercepts of the process \phi_{y,t}:
  boldA_yt <- array(0, dim=c(d*p, d*p, T_obs))
  phi_yt <- matrix(0, nrow=T_obs, ncol=d) # [T_obs, d]
  for(i1 in 1:T_obs) {
    tmp <- array(dim=c(d*p, d*p, M))
    tmp2 <- matrix(nrow=d, ncol=M)
    for(m in 1:M) {
      tmp[, , m] <- alpha_mt[i1, m]*all_boldA[, , m] # bold_A matrices weighted by transition weights
      tmp2[, m] <- alpha_mt[i1, m]*all_phi0[, m] # Intercept terms weighted by transition weights
    }
    boldA_yt[, , i1] <- apply(tmp, MARGIN=1:2, FUN=sum) # bold A matrix of the process obtained by weighting regime bold A
    phi_yt[i1, ] <- apply(tmp2, MARGIN=1, FUN=sum) # Intercept term of the process obtained by weighting regime intercepts
  }

  ## Calculate cumulative products of the bold A_yt matrices to a [dp, dp, T] array so that
  # the i:th slice contains the product of the bold A matrices from t=1 to t=i:
  cum_boldA_yt <- array(0, dim=c(d*p, d*p, T_obs))
  cum_boldA_yt[, , 1] <- boldA_yt[, , 1]
  for(i1 in 2:T_obs) {
    cum_boldA_yt[, , i1] <- boldA_yt[, , i1]%*%cum_boldA_yt[, , i1 - 1]
  }

  ## Obtain the data in a convenient form:
  # i:th row denotes the vector \bold{y_{i-1}} = (y_{i-1},...,y_{i-p}) (dpx1),
  # assuming the observed data is y_{-p+1},...,y_0,y_1,...,y_{T}
  Y <- reform_data(stvar$data, p) # (T+1 x dp)
  #Y <- Y[1:T_obs, , drop=FALSE] # Last row removed; not needed when calculating something based on lagged observations

  ## Calculate the component quantifying the contribution of the initial conditions (for all t=1,...,T):
  init_cond_comp <- matrix(NA, nrow=T_obs, ncol=d*p) # [t, i] = contrbtn of init conds to i:th variable at time t; the rest cols are for the companion form stuff
  for(i1 in 1:T_obs) {
    init_cond_comp[i1, ] <- cum_boldA_yt[, , i1]%*%Y[1,] # Initial conditions are the first row of Y
  }
  init_cond_comp <- init_cond_comp[, 1:d] # [t, i] = contribution of init conds to i:th variable at time t; extra cols removed here

  ## Calculate the so-called "steady state component" (for all t=1,...,T)
  ## AND
  ## the component quantifying the cumulative contributions of the shocks, as well as the cumulative contributions of the shocks.
  steady_state_comp <- matrix(NA, nrow=T_obs, ncol=d) # [t, i] = contrbution of steady state comp to i:th var at time t
  shock_comp <- matrix(NA, nrow=T_obs, ncol=d) # [t, i] = contribution of shocks to i:th variable at time t
  contributions_of_shocks <- array(NA, dim=c(T_obs, d, d)) # [t, shock, var], [t, i1, i2] = cum contrbtn of i1:th shock to i2:th variable at time t
  # Time t=1:
  steady_state_comp[1,] <- (H%*%phi_yt[1,])[1:d] # Removes the companion form extra stuff
  C_t0 <- H%*%B_yt[, , 1]
  shock_comp[1,] <- (C_t0%*%all_e_t[1,])[1:d] # Removes the companion form extra stuff
  for(k in 1:d) { # shock = k, variable = s
    for(s in 1:d) {
      contributions_of_shocks[1, k, s] <- C_t0[s, k]*all_e_t[1, k] # [t=1, shock=k, var=s]
    }
  }
  for(t in 2:T_obs) { # Time t=2,...,T
    # Calculate the coefficient matrices C_{t,j}, C_{t,0}=HB_{y,t}, C_{t,j} = (\prod_{i=0}^{j-1}A_{y,t-i})HB_{y,t}, for j=0,...,t-1:
    # j=0:
    C_tj <- array(0, dim=c(d*p, d, t)) # C_{t,j} in [dp, d, j+1]
    C_tj[, , 1] <- H%*%B_yt[, , t] # [dp, d] matrix C_{t,0}
    for(j in 1:(t-1)) { # j=1,...,t-1
      if(j == 1) {
        tmp_cum_A_yt <- boldA_yt[, , t] # [dp, dp] matrix \prod_{i=0}^{j-1}A_{y,t-i}
      } else {
        tmp_cum_A_yt <- tmp_cum_A_yt%*%boldA_yt[, , t - j + 1] # [dp, dp] matrix \prod_{i=0}^{j-1}A_{y,t-i}
      }
      C_tj[, , j + 1] <- tmp_cum_A_yt%*%C_tj[, , 1] # [dp, d] matrix C_{t,j}
    }
    # Calculate the steady state component:
    tmp <- matrix(nrow=d*p, ncol=t)
    inv_B_yt <- solve(B_yt[, , t]) # [d, d] matrix
    for(j in 0:(t-1)) {
      tmp[, j + 1] <- C_tj[, , j + 1]%*%inv_B_yt%*%phi_yt[t-j,] # [dp, 1] matrix
    }
    steady_state_comp[t, ] <- rowSums(tmp)[1:d] # (d x 1), removing the extra companion form stuff

    # Calculate the shock component:
    tmp <- matrix(nrow=d*p, ncol=t)
    for(j in 0:(t-1)) {
      tmp[, j + 1] <- C_tj[, , j + 1]%*%all_e_t[t-j,] # [dp, 1] matrix
    }
    shock_comp[t, ] <- rowSums(tmp)[1:d] # (d x 1), removing the extra companion form stuff

    # Calculate the contributions of the shocks:
    for(k in 1:d) { # shock = k
      for(s in 1:d) { # variable = s
        # C_tj[s, k, ] = c_sk elements in the order j=0,1,...,t-j and we need to match j:th C_tj with t-j:th e_t.
        # So, e.g., the first slice in C_tj (j=0) is multiplied with e_t (t-j=0), and the second slice in C_tj (j=1) is multiplied with e_{t-1} (t-j=1), etc.
        contributions_of_shocks[t, k, s] <- crossprod(C_tj[s, k, ], rev(all_e_t[1:t, k])) # [t, shock=k, var=s]
      }
    }
  }
  dimnames(contributions_of_shocks) <- list(1:T_obs, paste0("Shock ", 1:d), colnames(stvar$data))

  # Return the results
  structure(list(stvar=stvar,
                 init_cond_comp=init_cond_comp,
                 steady_state_comp=steady_state_comp,
                 shock_comp=shock_comp,
                 contributions_of_shocks=contributions_of_shocks),
            class="histdecomp")
}
