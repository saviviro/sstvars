% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hypoTests.R
\name{Wald_test}
\alias{Wald_test}
\title{Perform Wald test for a STVAR model}
\usage{
Wald_test(stvar, A, c)
}
\arguments{
\item{stvar}{an object of class \code{'stvar'} generated by \code{fitSTVAR} or \code{STVAR}, containing
the model specified by the alternative hypothesis (i.e., \strong{the unconstrained model}).}

\item{A}{a size \eqn{(k\times}\code{n_params}\eqn{)} matrix with full row rank specifying a part of the null hypothesis
where \code{n_params} is the number of parameters in the (unconstrained) model.
See details for more information.}

\item{c}{a length \eqn{k} vector specifying a part of the null hypothesis. See details for more information.}
}
\value{
A list with class "hypotest" containing the test results and arguments used to calculate the test.
}
\description{
\code{Wald_test} performs a Wald test for a STVAR model
}
\details{
Denoting the true parameter value by \eqn{\theta_{0}}, we test the null hypothesis \eqn{A\theta_{0}=c}.
  Under the null, the test statistic is asymptotically \eqn{\chi^2}-distributed with \eqn{k}
  (\code{=nrow(A)}) degrees of freedom. The parameter \eqn{\theta_{0}} is assumed to have the same form as in
  the model supplied in the argument \code{stvar} and it is presented in the documentation of the argument
  \code{params} in the function \code{STVAR} (see \code{?STVAR}).

  \strong{The test is based on the assumption of the standard result of asymptotic normality!}
  Also note that this function does \strong{not} check whether the model assumptions hold under the null.
}
\examples{
 # Logistic Student's t STVAR with p=1, M=2, and the first lag of the second variable
 # as the switching variable (parameter values were obtained by maximum likelihood estimation;
 # fitSTVAR is not used here because the estimation is computationally demanding).
 params12 <- c(0.62906848, 0.14245295, 2.41245785, 0.66719269, 0.3534745, 0.06041779, -0.34909745,
  0.61783824, 0.125769, -0.04094521, -0.99122586, 0.63805416, 0.371575, 0.00314754, 0.03440824,
  1.29072533, -0.06067807, 0.18737385, 1.21813844, 5.00884263, 7.70111672)
 fit12 <- STVAR(data=gdpdef, p=1, M=2, params=params12, weight_function="logistic",
  weightfun_pars=c(2, 1), cond_dist="Student")
 fit12

 # Test whether the location parameter equals 1.
 # For this model, the parameter vector has the length 21 and
 # location parameter is in the 19th element:
 A <- matrix(c(rep(0, times=18), 1, 0, 0), nrow=1, ncol=21)
 c <- 1
 Wald_test(fit12, A=A, c=c)

 # Test whether the intercepts and autoregressive matrices are identical across the regimes:
 # fit12 has parameter vector of length 21. In the first regime, the intercepts are in the
 # elements 1,2 and the AR parameters in the elements 5,...,8. In the second regime,
 # the intercepts are in the elements 3,4, and the AR parameters the elements 9,...,12.
 A <- rbind(cbind(diag(2), -diag(2), matrix(0, nrow=2, ncol=17)), # intercepts
            cbind(matrix(0, nrow=4, ncol=4), diag(4), -diag(4), matrix(0, nrow=4, ncol=9))) # AR
 c <- rep(0, times=6)
 Wald_test(fit12, A=A, c=c)
}
\references{
\itemize{
   \item Buse A. (1982). The Likelihood Ratio, Wald, and Lagrange Multiplier Tests: An Expository Note.
     \emph{The American Statistician}, 36(3a), 153-157.
 }
}
\seealso{
\code{\link{LR_test}}, \code{\link{Rao_test}}, \code{\link{fitSTVAR}}, \code{\link{STVAR}},
  \code{\link{diagnostic_plot}}, \code{\link{profile_logliks}}, \code{\link{Portmanteau_test}}
}
