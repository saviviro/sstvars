% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateParams.R
\name{random_weightpars}
\alias{random_weightpars}
\title{Create random transition weight parameter values}
\usage{
random_weightpars(M, weight_function, AR_constraints, mean_constraints)
}
\arguments{
\item{M}{a positive integer specifying the number of regimes}

\item{weight_function}{what type of transition weights should be used?}
}
\value{
Returns a numeric vector ...
  \describe{
    \item{If \code{weight_function == "relative_dens"}:}{a length \code{M-1} vector \eqn{(\alpha_1,...,\alpha_{M-1})}.}
    \item{If \code{weight_function == "logit"}:}{NOT YET IMPLEMENTED}
  }
}
\description{
\code{random_weightpars} generates random transition weight parameter values
}
\keyword{internal}