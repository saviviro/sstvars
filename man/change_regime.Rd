% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameterReforms.R
\name{change_regime}
\alias{change_regime}
\title{Change the parameters of a specific regime of the given parameter vector}
\usage{
change_regime(
  p,
  M,
  d,
  params,
  m,
  regime_pars,
  cond_dist = c("Gaussian", "Student", "ind_Student", "ind_skewed_t")
)
}
\arguments{
\item{p}{the autoregressive order of the model}

\item{M}{the number of regimes}

\item{d}{the number of time series in the system, i.e., the dimension}

\item{params}{a real valued vector specifying the parameter values.
Should have the form \eqn{\theta = (\phi_{1},...,\phi_{M},\varphi_1,...,\varphi_M,\sigma,\alpha,\nu)},
where (see exceptions below):
\itemize{
  \item{\eqn{\phi_{m} = } the \eqn{(d \times 1)} intercept (or mean) vector of the \eqn{m}th regime.}
  \item{\eqn{\varphi_m = (vec(A_{m,1}),...,vec(A_{m,p}))} \eqn{(pd^2 \times 1)}.}
  \item{\describe{
    \item{if \code{cond_dist="Gaussian"} or \code{"Student"}:}{\eqn{\sigma = (vech(\Omega_1),...,vech(\Omega_M))}
      \eqn{(Md(d + 1)/2 \times 1)}.}
    \item{if \code{cond_dist="ind_Student"} or \code{"ind_skewed_t"}:}{\eqn{\sigma = (vec(B_1),...,vec(B_M)} \eqn{(Md^2 \times 1)}.}
    }
  }
  \item{\eqn{\alpha = } the \eqn{(a\times 1)} vector containing the transition weight parameters (see below).}
  \item{\describe{
    \item{if \code{cond_dist = "Gaussian")}:}{Omit \eqn{\nu} from the parameter vector.}
    \item{if \code{cond_dist="Student"}:}{\eqn{\nu > 2} is the single degrees of freedom parameter.}
    \item{if \code{cond_dist="ind_Student"}:}{\eqn{\nu = (\nu_1,...,\nu_d)} \eqn{(d \times 1)}, \eqn{\nu_i > 2}.}
    \item{if \code{cond_dist="ind_skewed_t"}:}{\eqn{\nu = (\nu_1,...,\nu_d,\lambda_1,...,\lambda_d)} \eqn{(2d \times 1)},
     \eqn{\nu_i > 2} and \eqn{\lambda_i \in (0, 1)}.}
    }
  }
}
For models with...
\describe{
  \item{\code{weight_function="relative_dens"}:}{\eqn{\alpha = (\alpha_1,...,\alpha_{M-1})}
        \eqn{(M - 1 \times 1)}, where \eqn{\alpha_m} \eqn{(1\times 1), m=1,...,M-1} are the transition weight parameters.}
 \item{\code{weight_function="logistic"}:}{\eqn{\alpha = (c,\gamma)}
        \eqn{(2 \times 1)}, where \eqn{c\in\mathbb{R}} is the location parameter and \eqn{\gamma >0} is the scale parameter.}
  \item{\code{weight_function="mlogit"}:}{\eqn{\alpha = (\gamma_1,...,\gamma_M)} \eqn{((M-1)k\times 1)},
        where \eqn{\gamma_m} \eqn{(k\times 1)}, \eqn{m=1,...,M-1} contains the multinomial logit-regression coefficients
        of the \eqn{m}th regime. Specifically, for switching variables with indices in \eqn{I\subset\lbrace 1,...,d\rbrace}, and with
       \eqn{\tilde{p}\in\lbrace 1,...,p\rbrace} lags included, \eqn{\gamma_m} contains the coefficients for the vector
       \eqn{z_{t-1} = (1,\tilde{z}_{\min\lbrace I\rbrace},...,\tilde{z}_{\max\lbrace I\rbrace})}, where
       \eqn{\tilde{z}_{i} =(y_{it-1},...,y_{it-\tilde{p}})}, \eqn{i\in I}. So \eqn{k=1+|I|\tilde{p}}
       where \eqn{|I|} denotes the number of elements in \eqn{I}.}
  \item{\code{weight_function="exponential"}:}{\eqn{\alpha = (c,\gamma)}
        \eqn{(2 \times 1)}, where \eqn{c\in\mathbb{R}} is the location parameter and \eqn{\gamma >0} is the scale parameter.}
  \item{\code{weight_function="threshold"}:}{\eqn{\alpha = (r_1,...,r_{M-1})}
        \eqn{(M-1 \times 1)}, where \eqn{r_1,...,r_{M-1}} are the threshold values.}
  \item{\code{weight_function="exogenous"}:}{Omit \eqn{\alpha} from the parameter vector.}
  \item{\code{identification="heteroskedasticity"}:}{\eqn{\sigma = (vec(W),\lambda_2,...,\lambda_M)}, where
        \eqn{W} \eqn{(d\times d)} and \eqn{\lambda_m} \eqn{(d\times 1)}, \eqn{m=2,...,M}, satisfy
        \eqn{\Omega_1=WW'} and \eqn{\Omega_m=W\Lambda_mW'}, \eqn{\Lambda_m=diag(\lambda_{m1},...,\lambda_{md})},
        \eqn{\lambda_{mi}>0}, \eqn{m=2,...,M}, \eqn{i=1,...,d}.}
}
Above, \eqn{\phi_{m}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}th coefficient matrix of the \eqn{m}th
regime, \eqn{\Omega_{m}} denotes the positive definite error term covariance matrix of the \eqn{m}th regime, and \eqn{B_m}
is the invertible \eqn{(d\times d)} impact matrix of the \eqn{m}th regime. \eqn{\nu_m} is the degrees of freedom parameter
of the \eqn{m}th regime. If \code{parametrization=="mean"}, just replace each \eqn{\phi_{m}} with regimewise mean
\eqn{\mu_{m}}.}

\item{m}{which regime?}

\item{regime_pars}{\describe{
  \item{If \code{cond_dist="Gaussian"} or \code{cond_dist="Student"}:}{rhe \eqn{(dp + pd^2 + d(d+1)/2)} vector
   \eqn{(\phi_{m},vec(A_{m,1}),...,\vec(A_{m,p}),vech(\Omega_m))}.}
  \item{If \code{cond_dist="ind_Student"} or \code{"ind_skewed_t"}:}{the \eqn{(dp + pd^2 + d^2 + 1)} vector
   \eqn{(\phi_{m},vec(A_{m,1}),...,\vec(A_{m,p}),vec(B_m))}.}
}}
}
\value{
Returns parameter vector with \code{m}:th regime changed to \code{regime_pars}.
}
\description{
\code{change_regime} changes the regime parameters
  \eqn{(\phi_{m},vec(A_{m,1}),...,\vec(A_{m,p}),vech(\Omega_m))} (replace \eqn{vech(\Omega_m)}
  by \eqn{vec(B_m)} for \code{cond_dist="ind_Student"})
  of the given regime to the new given parameters.
}
\details{
Does not support constrained models or structural models. Weight parameters and distribution
  parameters are not changed.
}
\references{
\itemize{
   \item Kheifets I.L., Saikkonen P.J. 2020. Stationarity and ergodicity of Vector STAR models.
     \emph{Econometric Reviews}, \strong{39}:4, 407-414.
   \item Lütkepohl H. 2005. New Introduction to Multiple Time Series Analysis, \emph{Springer}.
   \item Lanne M., Virolainen S. 2025. A Gaussian smooth transition vector autoregressive model:
      An application to the macroeconomic effects of severe weather shocks. Unpublished working
      paper, available as arXiv:2403.14216.
   \item Virolainen S. 2025. Identification by non-Gaussianity in structural
      smooth transition vector autoregressive models. Unpublished working
      paper, available as arXiv:2404.19707.
 }
 @keywords internal
}
\keyword{internal}
