% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/moreEst.R
\name{fitSSTVAR}
\alias{fitSSTVAR}
\title{Maximum likelihood estimation of a structural STVAR model based on preliminary estimates from
  a reduced form model.}
\usage{
fitSSTVAR(
  stvar,
  identification = c("recursive", "heteroskedasticity", "non-Gaussianity"),
  B_constraints = NULL,
  B_pm_reg = NULL,
  B_perm = NULL,
  B_signs = NULL,
  maxit = 1000,
  maxit_robust = 1000,
  h = 0.001,
  robust_method = c("Nelder-Mead", "SANN", "none"),
  print_res = TRUE,
  calc_std_errors = TRUE
)
}
\arguments{
\item{stvar}{a an object of class \code{'stvar'}, created by, e.g., \code{fitSTVAR},
specifying a reduced form or a structural model}

\item{identification}{Which identification should the structural model use?
(see the vignette or the references for details)
 \describe{
   \item{\code{"recursive"}:}{The usual lower-triangular recursive identification of the shocks via their impact responses.}
   \item{\code{"heteroskedasticity"}:}{Identification by conditional heteroskedasticity, which imposes constant relative
     impact responses for each shock.}
 }}

\item{B_constraints}{Employ further constraints on the impact matrix?
A \eqn{(d \times d)} matrix with its entries imposing constraints on the impact matrix \eqn{B_t}:
\code{NA} indicating that the element is unconstrained, a positive value indicating strict positive sign constraint,
a negative value indicating strict negative sign constraint, and zero indicating that the element is constrained to zero.
Currently only available for models with \code{identification="heteroskedasticity"} due to the (in)availability of appropriate
parametrizations that allow such constraints to be imposed.}

\item{B_pm_reg}{an integer between \eqn{1} and \eqn{M} specifying the regime the permutations and sign changes of \eqn{B_m}
specified in the arguments \code{B_perm} and \code{B_signs} are applied to.}

\item{B_perm}{a numeric vector of length \eqn{d} specifying the permutation of the columns of the impact matrix \eqn{B_m}
of a single regime specified in the argument \code{B_pm_reg} prior to re-estimating the model. Applicable only for models
with \code{cond_dist = "ind_Student"} or \code{"ind_skewed_t"}.}

\item{B_signs}{a numeric vector specifying the columns of the impact matrix of a single regime specified in the argument
\code{B_pm_reg} that should be multiplied by -1 \strong{prior} to reordering them according to \code{B_perm} (if specified).
 Applicable only for models with \code{cond_dist = "ind_Student"} or \code{"ind_skewed_t"}.}

\item{maxit}{the maximum number of iterations in the variable metric algorithm.}

\item{maxit_robust}{the maximum number of iterations on the first phase robust estimation, if employed.}

\item{h}{the strictly positive difference used in the finite difference approximation of the gradient used in numerical optimization.}

\item{robust_method}{Should some robust estimation method be used in the estimation before switching
to the gradient based variable metric algorithm? See details.}

\item{print_res}{should summaries of estimation results be printed?}

\item{calc_std_errors}{Calculate approximate standard errors (based on standard asymptotics)?}
}
\value{
Returns an S3 object of class \code{'stvar'} defining a smooth transition VAR model. The returned list
 contains the following components (some of which may be \code{NULL} depending on the use case):
   \item{data}{The input time series data.}
   \item{model}{A list describing the model structure.}
   \item{params}{The parameters of the model.}
   \item{std_errors}{Approximate standard errors of the parameters, if calculated.}
   \item{transition_weights}{The transition weights of the model.}
   \item{regime_cmeans}{Conditional means of the regimes, if data is provided.}
   \item{total_cmeans}{Total conditional means of the model, if data is provided.}
   \item{total_ccovs}{Total conditional covariances of the model, if data is provided.}
   \item{uncond_moments}{A list of unconditional moments including regime autocovariances, variances, and means.}
   \item{residuals_raw}{Raw residuals, if data is provided.}
   \item{residuals_std}{Standardized residuals, if data is provided.}
   \item{structural_shocks}{Recovered structural shocks, if applicable.}
   \item{loglik}{Log-likelihood of the model, if data is provided.}
   \item{IC}{The values of the information criteria (AIC, HQIC, BIC) for the model, if data is provided.}
   \item{all_estimates}{The parameter estimates from all estimation rounds, if applicable.}
   \item{all_logliks}{The log-likelihood of the estimates from all estimation rounds, if applicable.}
   \item{which_converged}{Indicators of which estimation rounds converged, if applicable.}
   \item{which_round}{Indicators of which round of optimization each estimate belongs to, if applicable.}
   \item{seeds}{The seeds used in the estimation in \code{fitSTVAR}, if applicable.}
   \item{LS_estimates}{The least squares estimates of the parameters in the form
     \eqn{(\phi_{1},...,\phi_{M},\varphi_1,...,\varphi_M,\alpha} (intercepts replaced by unconditional means
     if mean parametrization is used), if applicable.}
}
\description{
\code{fitSSTVAR} uses a robust method and a variable metric algorithm to estimate
  a structural STVAR model based on preliminary estimates from a reduced form model.
}
\details{
When the structural model does not impose overidentifying constraints, it is directly
  obtained from the reduced form model, and estimation is not required. When overidentifying constraints
  are imposed, the model is estimated subject to the constraints.

  Using the robust estimation method before switching to the variable metric can be useful if the initial
  estimates are not very close to the ML estimate of the structural model, as the variable metric algorithm
  (usually) converges to a nearby local maximum or saddle point. However, if the initial estimates are far from
  the ML estimate, the resulting solution is likely local only due to the complexity of the model. Note that
  Nelder-Mead algorithm is much faster than SANN but can get stuck at a local solution.
  This is particularly the case when the imposed overidentifying restrictions are such that the unrestricted
  estimate is not close to satisfying them. Nevertheless, in most practical cases, the model is just identified
  and estimation is not required, and often reasonable overidentifying constraints are close to the unrestricted estimate.

  Employs the estimation function \code{optim} from the package \code{stats} that implements the optimization
  algorithms. See \code{?optim} for the documentation on the optimization methods.

  The arguments \code{B_pm_reg}, \code{B_perm}, and \code{B_signs} can be used to explore estimates based various orderings
  and sign changes of the columns of the impact matrices \eqn{B_m} of specific regimes. This can be useful in the presence
  of weak identification with respect to the ordering or signs of the columns \eqn{B_2,...,B_M} (see Virolainen 2025).
}
\examples{
\donttest{
## These are long running examples that take approximately 1 minute to run.

## Estimate first a reduced form Gaussian STVAR p=3, M=2 model with the weighted relative
# stationary densities of the regimes as the transition weight function, and the means and
# AR matrices constrained to be identical across the regimes:
fit32cm <- fitSTVAR(gdpdef, p=3, M=2, AR_constraints=rbind(diag(3*2^2), diag(3*2^2)),
  weight_function="relative_dens", mean_constraints=list(1:2), parametrization="mean",
  nrounds=1, seeds=1, ncores=1)

# Then, we estimate/create various structural models based on the reduced form model.
# Create a structural model with the shocks identified recursively:
fit32cms_rec <- fitSSTVAR(fit32cm, identification="recursive")

# Create a structural model with the shocks identified by conditional heteroskedasticity:
fit32cms_hetsked <- fitSSTVAR(fit32cm, identification="heteroskedasticity")
fit32cms_hetsked # Print the estimates

# Estimate a structural model with the shocks identified by conditional heteroskedasticity
# and overidentifying constraints imposed on the impact matrix: positive diagonal element
# and zero upper right element:
fit32cms_hs2 <- fitSSTVAR(fit32cm, identification="heteroskedasticity",
 B_constraints=matrix(c(1, NA, 0, 1), nrow=2))

# Estimate a structural model with the shocks identified by conditional heteroskedasticity
# and overidentifying constraints imposed on the impact matrix: positive diagonal element
# and zero off-diagonal elements:
fit32cms_hs3 <- fitSSTVAR(fit32cms_hs2, identification="heteroskedasticity",
 B_constraints=matrix(c(1, 0, 0, 1), nrow=2))

# Estimate first a reduced form two-regime Threshold VAR p=1 model with
# with independent skewed t shocks, and the first lag of the second variable
# as the switching variable, and AR matrices constrained to be identical
# across the regimes:
fit12c <- fitSTVAR(gdpdef, p=1, M=2, cond_dist="ind_skewed_t",
 AR_constraints=rbind(diag(1*2^2), diag(1*2^2)), weight_function="threshold",
 weightfun_pars=c(2, 1), nrounds=1, seeds=1, ncores=1)

# Due to the independent non-Gaussian shocks, the structural shocks are readily
# identified. The following returns the same model but marked as structural
# with the shocks identified by non-Gaussianity:
fit12c <- fitSSTVAR(fit12c)

# Estimate a model based on a reversed ordering of the columns of the impact matrix B_2:
fit12c2 <- fitSSTVAR(fit12c, B_pm_reg=2, B_perm=c(2, 1))

# Estimate a model based on reversed signs of the second column of B_2 and reversed
# ordering of the columns of B_2:
fit12c3 <- fitSSTVAR(fit12c, B_pm_reg=2, B_perm=c(2, 1), B_signs=2)
}
}
\references{
\itemize{
   \item Kilian L., Lütkepohl H. 20017. Structural Vector Autoregressive Analysis. 1st edition.
   \emph{Cambridge University Press}, Cambridge.
   \item Lütkepohl H., Netšunajev A. 2017. Structural vector autoregressions with smooth transition in variances.
     \emph{Journal of Economic Dynamics & Control}, \strong{84}, 43-57.
   \item Virolainen S. 2025. Identification by non-Gaussianity in structural
      smooth transition vector autoregressive models. Unpublished working
      paper, available as arXiv:2404.19707.
 }
}
\seealso{
\code{\link{fitSTVAR}}, \code{\link{STVAR}}, \code{\link[stats]{optim}}
}
