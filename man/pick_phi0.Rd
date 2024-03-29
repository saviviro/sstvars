% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pickParams.R
\name{pick_phi0}
\alias{pick_phi0}
\title{Pick \eqn{\phi_{m,0}} or \eqn{\mu_{m}}, m=1,..,M vectors}
\usage{
pick_phi0(M, d, params)
}
\arguments{
\item{M}{the number of regimes}

\item{d}{the number of time series in the system, i.e., the dimension}

\item{params}{a real valued vector specifying the parameter values.
Should have the form \eqn{\theta = (\phi_{1,0},...,\phi_{M,0},\varphi_1,...,\varphi_M,\sigma,\alpha,\nu)},
where
\itemize{
  \item{\eqn{\phi_{m,0} = } the \eqn{(d \times 1)} intercept (or mean) vector of the \eqn{m}th regime.}
  \item{\eqn{\varphi_m = (vec(A_{m,1}),...,vec(A_{m,p}))} \eqn{(pd^2 \times 1)}.}
  \item{\eqn{\sigma = (vech(\Omega_1),...,vech(\Omega_M)} \eqn{(Md(d + 1)/2 \times 1)}.}
  \item{\eqn{\alpha} contains the transition weights parameters.}
  \item{\eqn{\nu > 2} is the degrees of freedom parameter that is included only if \code{cond_dist="Student"}.}
}
\describe{
  \item{\code{weight_function="relative_dens"}:}{\eqn{\alpha = (\alpha_1,...,\alpha_{M-1})}
        \eqn{(M - 1 \times 1)}, where \eqn{\alpha_m} \eqn{(1\times 1), m=1,...,M-1} are the transition weight parameters.}
 \item{\code{weight_function="logistic"}:}{\eqn{\alpha = (c,\gamma)}
        \eqn{(2 \times 1)}, where \eqn{c\in\mathbb{R}} is the location parameter and \eqn{\gamma >0} is the scale parameter.}
  \item{\code{weight_function="mlogit"}:}{\eqn{\alpha = (\gamma_1,...,\gamma_M)} \eqn{((M-1)k\times 1)},
        where \eqn{\gamma_m} \eqn{(k\times 1)}, \eqn{m=1,...,M-1} contains the multinomial logit-regression coefficients
        of the \eqn{m}th regime. Specifically, for switching variables with indices in \eqn{I\subset\lbrace 1,...,d\rbrace},
        and with \eqn{\tilde{p}\in\lbrace 1,...,p\rbrace} lags included, \eqn{\gamma_m} contains the coefficients for the vector
       \eqn{z_{t-1} = (1,\tilde{z}_{\min\lbrace I\rbrace},...,\tilde{z}_{\max\lbrace I\rbrace})}, where
       \eqn{\tilde{z}_{i} =(y_{it-1},...,y_{it-\tilde{p}})}, \eqn{i\in I}. So \eqn{k=1+|I|\tilde{p}}
       where \eqn{|I|} denotes the number of elements in \eqn{I}.}
  \item{\code{weight_function="exponential"}:}{\eqn{\alpha = (c,\gamma)}
        \eqn{(2 \times 1)}, where \eqn{c\in\mathbb{R}} is the location parameter and \eqn{\gamma >0} is the scale parameter.}
  \item{\code{weight_function="threshold"}:}{\eqn{\alpha = (r_1,...,r_{M-1})}
        \eqn{(M-1 \times 1)}, where \eqn{r_1,...,r_{M-1}} are the threshold values.}
}
Above, \eqn{\phi_{m,0}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}th coefficient matrix of the \eqn{m}th
mixture component, and \eqn{\Omega_{m}} denotes the error term covariance matrix of the \eqn{m}:th mixture component.
If \code{parametrization=="mean"}, just replace each \eqn{\phi_{m,0}} with regimewise mean \eqn{\mu_{m}}.
\eqn{vec()} is vectorization operator that stacks columns of a given matrix into a vector. \eqn{vech()} stacks columns
of a given matrix from the principal diagonal downwards (including elements on the diagonal) into a vector.}
}
\value{
Returns a \eqn{(dxM)} matrix containing \eqn{\phi_{m,0}} in the m:th column or
  \eqn{\mu_{m}} if the parameter vector is mean-parametrized, \eqn{, m=1,..,M}.
}
\description{
\code{pick_phi0} picks the intercept or mean parameters from the given parameter vector.
}
\details{
Does not support constrained parameter vectors.
}
\section{Warning}{

 No argument checks!
}

\keyword{internal}
