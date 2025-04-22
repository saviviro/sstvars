#' @title Compute historical decompositions for structural STVAR models.
#'
#' @description \code{hist_decomp} compute historical decompositions for structural STVAR models.
#'
#' @inheritParams simulate.stvar
#' @param tmp
#' @details The historical decomposition quantifies the cumulative effects the shocks to the movements of
#'   the variables (see, e.g., Kilian and Lütkepohl, 2017, Section~4.3) The historical decompositions are
#'   computed as described in Wong (2018). Note that due to the effect of the "initial conditions" and the
#'   "steady state component", which are not attributed to the shocks, the cumulative effects of the shocks
#'   do not sum to the observed time series.
#' @return Returns a class \code{'histdecomp'} list with ...
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

  # Obtain the parameter values in the "non constrained form".
  params <- reform_constrained_pars(p=p, M=M, d=d, params=stvar$params, weight_function=weight_function, cond_dist=cond_dist,
                                    identification=identification, AR_constraints=stvar$model$AR_constraints,
                                    mean_constraints=stvar$model$mean_constraints, weight_constraints=stvar$model$weight_constraints,
                                    B_constraints=stvar$model$B_constraints, other_constraints=NULL, weightfun_pars=weightfun_pars)

  # Pick params
  all_A <- pick_allA(p=p, M=M, d=d, params=params) # [d, d, p, M]
  all_Omegas <- pick_Omegas(p=p, M=M, d=d, params=params, cond_dist=cond_dist, identification=identification) # [d, d, M], B_m for ind_Stud and ind_skewed_t
  weightpars <- pick_weightpars(p=p, M=M, d=d, params=params, weight_function=weight_function, cond_dist=cond_dist, weightfun_pars=weightfun_pars)
  # all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A) # Use similar for for each time period to get bold A for A_t mats.
  distpars <- pick_distpars(d=d, params=params, cond_dist=cond_dist)
  if(parametrization == "intercept") { # [d, M]
    all_phi0 <- pick_phi0(M=M, d=d, params=params)
  } else {
    all_phi0 <- vapply(1:M, function(m) (Id - rowSums(all_A[, , , m, drop=FALSE], dims=2))%*%all_mu[,m], numeric(d))
  }

  ## Create the (dp x d) matrix H:
  H <- rbind(diag(d), matrix(0, nrow=(d*(p-1)), ncol=d))

  ## Create the [d, d, T] array of impact matrices B_{y,t}:
  B_yt <- array(0, dim=c(d, d, T_obs))

  ## Obtain the data in a convenient form:
  # i:th row denotes the vector \bold{y_{i-1}} = (y_{i-1},...,y_{i-p}) (dpx1),
  # assuming the observed data is y_{-p+1},...,y_0,y_1,...,y_{T}
  Y <- reform_data(stvar$data, p) # (T+1 x dp)
  #Y <- Y[1:T_obs, , drop=FALSE] # Last row removed; not needed when calculating something based on lagged observations

  # Huomaa, että illustroinnin voi tehdä joko suoraan total contributionina tai
  # Kilian and Leen tapaa percent deviationina from the mean. Itse en ehkä pidä
  # jälkimmäisestä, koska jos mean on esim melkein nolla, tai tosi suuri pienellä varianssilla,
  # niin se vääristää niitä arvoja.
}
