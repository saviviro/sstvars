% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generateParams.R
\name{random_ind}
\alias{random_ind}
\title{Create random mean parametrized parameter vector}
\usage{
random_ind(
  p,
  M,
  d,
  weight_function = c("relative_dens", "logit"),
  cond_dist = c("Gaussian", "Student"),
  AR_constraints = NULL,
  mean_constraints = NULL,
  force_stability = is.null(AR_constraints),
  mu_scale,
  mu_scale2,
  omega_scale,
  ar_scale = 1,
  ar_scale2 = 1
)
}
\arguments{
\item{p}{a positive integer specifying the autoregressive order}

\item{M}{a positive integer specifying the number of regimes}

\item{weight_function}{what type of transition weights should be used?}

\item{cond_dist}{specifies the conditional distribution of the model as \code{"Gaussian"} or \eqn{"Student"}.}

\item{force_stability}{Should the algorithm proposed by Ansley and Kohn be used to generate
AR matrices that always satisfy the stability condition? Not supported if AR constraints are
employed.}

\item{mu_scale}{a size \eqn{(dx1)} vector defining \strong{means} of the normal distributions from which each
mean parameter \eqn{\mu_{m}} is drawn from in random mutations. Default is \code{colMeans(data)}. Note that
mean-parametrization is always used for optimization in \code{GAfit} - even when \code{parametrization=="intercept"}.
However, input (in \code{initpop}) and output (return value) parameter vectors can be intercept-parametrized.}

\item{mu_scale2}{a size \eqn{(dx1)} strictly positive vector defining \strong{standard deviations} of the normal
distributions from which each mean parameter \eqn{\mu_{m}} is drawn from in random mutations.
Default is \code{2*sd(data[,i]), i=1,..,d}.}

\item{omega_scale}{a size \eqn{(dx1)} strictly positive vector specifying the scale and variability of the
random covariance matrices in random mutations. The covariance matrices are drawn from (scaled) Wishart
distribution. Expected values of the random covariance matrices are \code{diag(omega_scale)}. Standard
deviations of the diagonal elements are \code{sqrt(2/d)*omega_scale[i]}
and for non-diagonal elements they are \code{sqrt(1/d*omega_scale[i]*omega_scale[j])}.
Note that for \code{d>4} this scale may need to be chosen carefully. Default in \code{GAfit} is
\code{var(stats::ar(data[,i], order.max=10)$resid, na.rm=TRUE), i=1,...,d}. This argument is ignored if
structural model is considered.}

\item{ar_scale}{a positive real number adjusting how large AR parameter values are typically proposed in construction
of the initial population: larger value implies larger coefficients (in absolute value). After construction of the
initial population, a new scale is drawn from \code{(0, 0.)} uniform distribution in each iteration.}

\item{ar_scale2}{a positive real number adjusting how large AR parameter values are typically proposed in some
random mutations (if AR constraints are employed, in all random mutations): larger value implies \strong{smaller}
coefficients (in absolute value). \strong{Values larger than 1 can be used if the AR coefficients are expected to
be very small. If set smaller than 1, be careful as it might lead to failure in the creation of parameter candidates
that satisfy the stability condition.}}
}
\value{
Returns random mean parametrized parameter vector that has the same form as the argument \code{params}
  in the other functions, for instance, in the function \code{loglikelihood}.
}
\description{
\code{random_ind} generates random mean parametrized parameter vector.
}
\details{
Structural models are not supported!
}
\references{
\itemize{
   \item Ansley C.F., Kohn R. 1986. A note on reparameterizing a vector autoregressive
      moving average model to enforce stationarity.
      \emph{Journal of statistical computation and simulation}, \strong{24}:2, 99-106.
 }
}
\keyword{internal}
