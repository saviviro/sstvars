% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulate.R
\name{simulate_from_regime}
\alias{simulate_from_regime}
\title{Simulate observations from a regime of a STVAR model}
\usage{
simulate_from_regime(
  stvar,
  regime = 1,
  nsim = 1,
  init_values = NULL,
  use_transweights = TRUE
)
}
\arguments{
\item{stvar}{an object of class \code{'stvar'}.}

\item{regime}{an integer in \eqn{1,...,M} determining the regime from which to simulate observations from}

\item{nsim}{number of observations to be simulated.}

\item{init_values}{a size \eqn{(p\times d)} matrix specifying the initial values, where d is the number
of time series in the system. The \strong{last} row will be used as initial values for the first lag,
the second last row for second lag etc. If not specified, initial values are set to the unconditional
mean of the regime.}

\item{use_transweights}{if \code{TRUE} will calculate the transition weights of the provided model, simulate
\eqn{p + 100} observations more, calculate the transition weights for the last \eqn{100} observations, and
return the the consecutive \eqn{p} observations have the highest transition weight for the specified regime.}
}
\value{
\describe{
    \item{If \code{use_transweights=FALSE}:}{Returns a \eqn{(nsim \times d)} matrix such that the \eqn{t}th row
      contains the \eqn{t}th simulated observation.}
    \item{If \code{use_transweights=TRUE}:}{Returns a \eqn{(p \times d)} such that the \eqn{t}th row constrains
      the \eqn{t}th observations.}
  }
}
\description{
\code{simulate_from_regime} allows to simulate observations from a single
  regime of a STVAR model
}
\details{
Does not take random number generator seed as an argument to avoid unwanted behavior,
   because \code{simulate_from_regime} is mostly called from \code{simulate.stvar}
   that takes a seed as its argument, and \code{simulate_from_regime} calls \code{simulate.stvar} to simulate the observations.
   Specifically, \code{simulate_from_regime} generates a STVAR model from the given regime, sets up the initial values to the
   (if not specified), and then calls \code{simulate.stvar} accordingly.
}
\references{
\itemize{
   \item Anderson H., Vahid F. 1998. Testing multiple equation systems for common nonlinear components.
     \emph{Journal of Econometrics}, \strong{84}:1, 1-36.
   \item Hansen B.E. 1994. Autoregressive Conditional Density estimation.
     \emph{Journal of Econometrics}, \strong{35}:3, 705-730.
   \item Kheifets I.L., Saikkonen P.J. 2020. Stationarity and ergodicity of Vector STAR models.
     \emph{International Economic Review}, \strong{35}:3, 407-414.
   \item Lanne M., Virolainen S. 2025. A Gaussian smooth transition vector autoregressive model:
      An application to the macroeconomic effects of severe weather shocks. Unpublished working
      paper, available as arXiv:2403.14216.
   \item Lütkepohl H. 2005. New Introduction to Multiple Time Series Analysis,
         \emph{Springer}.
   \item McElroy T. 2017. Computation of vector ARMA autocovariances.
         \emph{Statistics and Probability Letters}, \strong{124}, 92-96.
   \item Kilian L., Lütkepohl H. 20017. Structural Vector Autoregressive Analysis. 1st edition.
     \emph{Cambridge University Press}, Cambridge.
   \item Tsay R. 1998. Testing and Modeling Multivariate Threshold Models.
     \emph{Journal of the American Statistical Association}, \strong{93}:443, 1188-1202.
   \item Virolainen S. 2025. Identification by non-Gaussianity in structural
      smooth transition vector autoregressive models. Unpublished working
      paper, available as arXiv:2404.19707.
 }
}
\seealso{
\code{\link{simulate.stvar}}
}
\keyword{internal}
