% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/counterFactuals.R
\name{get_phi_yt}
\alias{get_phi_yt}
\title{Compute the intercept \eqn{\phi_{y,t}=\sum_{m=1}^M \alpha_{mt} \phi_{m}} parameter value for a single time period}
\usage{
get_phi_yt(all_phi0, alpha_mt)
}
\arguments{
\item{all_phi0}{a \eqn{(d \times M)} matrix such that the \eqn{m}th column contains the intercept parameters of the \eqn{m}th regime.}

\item{alpha_mt}{an \eqn{(M \times 1)} vector containing the time period \eqn{t} transition weights.}
}
\value{
Returns the \eqn{(d \times 1)} vector of the intercept parameter values for the time period \eqn{t}.
}
\description{
\code{get_phi_yt} computes the intercept \eqn{\phi_{y,t}=\sum_{m=1}^M \alpha_{mt} \phi_{m}} parameter
 value for a single time period based on the regime intercepts and transition weights.
}
\details{
This is used in simulation of the counterfactual scenarios.
}
\keyword{internal}
