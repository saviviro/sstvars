% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/residuals.R
\name{get_residuals}
\alias{get_residuals}
\title{Calculate residuals of a smooth transition VAR}
\usage{
get_residuals(
  data,
  p,
  M,
  params,
  weight_function = c("relative_dens", "logistic", "mlogit", "exponential", "threshold",
    "exogenous"),
  weightfun_pars = NULL,
  cond_dist = c("Gaussian", "Student", "ind_Student", "ind_skewed_t"),
  parametrization = c("intercept", "mean"),
  identification = c("reduced_form", "recursive", "heteroskedasticity",
    "non-Gaussianity"),
  AR_constraints = NULL,
  mean_constraints = NULL,
  weight_constraints = NULL,
  B_constraints = NULL,
  standardize = TRUE,
  structural_shocks = FALSE,
  penalized = FALSE,
  penalty_params = c(0.05, 0.5),
  allow_unstab = FALSE
)
}
\arguments{
\item{data}{a matrix or class \code{'ts'} object with \code{d>1} columns. Each column is taken to represent
a univariate time series. Missing values are not supported.}

\item{p}{a positive integer specifying the autoregressive order}

\item{M}{a positive integer specifying the number of regimes}

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
        \eqn{(M-1 \times 1)}, where \eqn{r_1,...,r_{M-1}} are the thresholds.}
  \item{\code{weight_function="exogenous"}:}{Omit \eqn{\alpha} from the parameter vector.}
  \item{AR_constraints:}{Replace \eqn{\varphi_1,...,\varphi_M} with \eqn{\psi} as described in the argument \code{AR_constraints}.}
  \item{mean_constraints:}{Replace \eqn{\phi_{1},...,\phi_{M}} with \eqn{(\mu_{1},...,\mu_{g})} where
        \eqn{\mu_i, \ (d\times 1)} is the mean parameter for group \eqn{i} and \eqn{g} is the number of groups.}
  \item{weight_constraints:}{If linear constraints are imposed, replace \eqn{\alpha} with \eqn{\xi} as described in the
   argument \code{weigh_constraints}. If weight functions parameters are imposed to be fixed values, simply drop \eqn{\alpha}
   from the parameter vector.}
  \item{\code{identification="heteroskedasticity"}:}{\eqn{\sigma = (vec(W),\lambda_2,...,\lambda_M)}, where
        \eqn{W} \eqn{(d\times d)} and \eqn{\lambda_m} \eqn{(d\times 1)}, \eqn{m=2,...,M}, satisfy
        \eqn{\Omega_1=WW'} and \eqn{\Omega_m=W\Lambda_mW'}, \eqn{\Lambda_m=diag(\lambda_{m1},...,\lambda_{md})},
        \eqn{\lambda_{mi}>0}, \eqn{m=2,...,M}, \eqn{i=1,...,d}.}
  \item{B_constraints:}{For models identified by heteroskedasticity, replace \eqn{vec(W)} with \eqn{\tilde{vec}(W)}
        that stacks the columns of the matrix \eqn{W} in to vector so that the elements that are constrained to zero
        are not included. For models identified by non-Gaussianity, replace \eqn{vec(B_1),...,vec(B_M)} with
        similarly with vectorized versions \eqn{B_m} so that the elements that are constrained to zero are not included.}
}
Above, \eqn{\phi_{m}} is the intercept parameter, \eqn{A_{m,i}} denotes the \eqn{i}th coefficient matrix of the \eqn{m}th
regime, \eqn{\Omega_{m}} denotes the positive definite error term covariance matrix of the \eqn{m}th regime, and \eqn{B_m}
is the invertible \eqn{(d\times d)} impact matrix of the \eqn{m}th regime. \eqn{\nu_m} is the degrees of freedom parameter
of the \eqn{m}th regime.
If \code{parametrization=="mean"}, just replace each \eqn{\phi_{m}} with regimewise mean \eqn{\mu_{m}}.
\eqn{vec()} is vectorization operator that stacks columns of a given matrix into a vector. \eqn{vech()} stacks columns
of a given matrix from the principal diagonal downwards (including elements on the diagonal) into a vector. \eqn{Bvec()}
is a vectorization operator that stacks the columns of a given impact matrix \eqn{B_m} into a vector so that the elements
that are constrained to zero by the argument \code{B_constraints} are excluded.}

\item{weight_function}{What type of transition weights \eqn{\alpha_{m,t}} should be used?
\describe{
  \item{\code{"relative_dens"}:}{\eqn{\alpha_{m,t}=
    \frac{\alpha_mf_{m,dp}(y_{t-1},...,y_{t-p+1})}{\sum_{n=1}^M\alpha_nf_{n,dp}(y_{t-1},...,y_{t-p+1})}}, where
    \eqn{\alpha_m\in (0,1)} are weight parameters that satisfy \eqn{\sum_{m=1}^M\alpha_m=1} and
    \eqn{f_{m,dp}(\cdot)} is the \eqn{dp}-dimensional stationary density of the \eqn{m}th regime corresponding to \eqn{p}
    consecutive observations. Available for Gaussian conditional distribution only.}
  \item{\code{"logistic"}:}{\eqn{M=2}, \eqn{\alpha_{1,t}=1-\alpha_{2,t}},
    and \eqn{\alpha_{2,t}=[1+\exp\lbrace -\gamma(y_{it-j}-c) \rbrace]^{-1}}, where \eqn{y_{it-j}} is the lag \eqn{j}
    observation of the \eqn{i}th variable, \eqn{c} is a location parameter, and \eqn{\gamma > 0} is a scale parameter.}
  \item{\code{"mlogit"}:}{\eqn{\alpha_{m,t}=\frac{\exp\lbrace \gamma_m'z_{t-1} \rbrace}
    {\sum_{n=1}^M\exp\lbrace \gamma_n'z_{t-1} \rbrace}}, where \eqn{\gamma_m} are coefficient vectors, \eqn{\gamma_M=0},
    and \eqn{z_{t-1}} \eqn{(k\times 1)} is the vector containing a constant and the (lagged) switching variables.}
  \item{\code{"exponential"}:}{\eqn{M=2}, \eqn{\alpha_{1,t}=1-\alpha_{2,t}},
    and \eqn{\alpha_{2,t}=1-\exp\lbrace -\gamma(y_{it-j}-c) \rbrace}, where \eqn{y_{it-j}} is the lag \eqn{j}
    observation of the \eqn{i}th variable, \eqn{c} is a location parameter, and \eqn{\gamma > 0} is a scale parameter.}
  \item{\code{"threshold"}:}{\eqn{\alpha_{m,t} = 1} if \eqn{r_{m-1}<y_{it-j}\leq r_{m}} and \eqn{0} otherwise, where
     \eqn{-\infty\equiv r_0<r_1<\cdots <r_{M-1}<r_M\equiv\infty} are thresholds \eqn{y_{it-j}} is the lag \eqn{j}
     observation of the \eqn{i}th variable.}
  \item{\code{"exogenous"}:}{Exogenous nonrandom transition weights, specify the weight series in \code{weightfun_pars}.}
}
See the vignette for more details about the weight functions.}

\item{weightfun_pars}{\describe{
  \item{If \code{weight_function == "relative_dens"}:}{Not used.}
  \item{If \code{weight_function \%in\% c("logistic", "exponential", "threshold")}:}{a numeric vector with the switching variable
    \eqn{i\in\lbrace 1,...,d \rbrace} in the first and the lag \eqn{j\in\lbrace 1,...,p \rbrace} in the second element.}
  \item{If \code{weight_function == "mlogit"}:}{a list of two elements:
    \describe{
      \item{The first element \code{$vars}:}{a numeric vector containing the variables that should used as switching variables
        in the weight function in an increasing order, i.e., a vector with unique elements in \eqn{\lbrace 1,...,d \rbrace}.}
      \item{The second element \code{$lags}:}{an integer in \eqn{\lbrace 1,...,p \rbrace} specifying the number of lags to be
        used in the weight function.}
    }
  }
  \item{If \code{weight_function == "exogenous"}:}{a size (\code{nrow(data) - p} x \code{M}) matrix containing the exogenous
    transition weights as \code{[t, m]} for time \eqn{t} and regime \eqn{m}. Each row needs to sum to one and only weakly positive
    values are allowed.}
}}

\item{cond_dist}{specifies the conditional distribution of the model as \code{"Gaussian"}, \code{"Student"}, \code{"ind_Student"},
or \code{"ind_skewed_t"}, where \code{"ind_Student"} the Student's \eqn{t} distribution with independent components, and
\code{"ind_skewed_t"} is the skewed \eqn{t} distribution with independent components (see Hansen, 1994).}

\item{parametrization}{\code{"intercept"} or \code{"mean"} determining whether the model is parametrized with intercept
parameters \eqn{\phi_{m}} or regime means \eqn{\mu_{m}}, m=1,...,M.}

\item{identification}{is it reduced form model or an identified structural model; if the latter, how is it identified
(see the vignette or the references for details)?
\describe{
  \item{\code{"reduced_form"}:}{Reduced form model.}
  \item{\code{"recursive"}:}{The usual lower-triangular recursive identification of the shocks via their impact responses.}
  \item{\code{"heteroskedasticity"}:}{Identification by conditional heteroskedasticity, which imposes constant relative
    impact responses for each shock.}
  \item{\code{"non-Gaussianity"}:}{Identification by non-Gaussianity; requires mutually independent non-Gaussian shocks, thus,
    currently available only with the conditional distribution \code{"ind_Student"}.}
}}

\item{AR_constraints}{a size \eqn{(Mpd^2 \times q)} constraint matrix \eqn{C} specifying linear constraints
to the autoregressive parameters. The constraints are of the form
\eqn{(\varphi_{1},...,\varphi_{M}) = C\psi}, where \eqn{\varphi_{m} = (vec(A_{m,1}),...,vec(A_{m,p})) \ (pd^2 \times 1),\ m=1,...,M},
contains the coefficient matrices and \eqn{\psi} \eqn{(q \times 1)} contains the related parameters.
For example, to restrict the AR-parameters to be the identical across the regimes, set \eqn{C =}
[\code{I:...:I}]' \eqn{(Mpd^2 \times pd^2)} where \code{I = diag(p*d^2)}.}

\item{mean_constraints}{Restrict the mean parameters of some regimes to be identical? Provide a list of numeric vectors
such that each numeric vector contains the regimes that should share the common mean parameters. For instance, if
\code{M=3}, the argument \code{list(1, 2:3)} restricts the mean parameters of the second and third regime to be
identical but the first regime has freely estimated (unconditional) mean. Ignore or set to \code{NULL} if mean parameters
should not be restricted to be the same among any regimes. This constraint is available only for mean parametrized models;
that is, when \code{parametrization="mean"}.}

\item{weight_constraints}{a list of two elements, \eqn{R} in the first element and \eqn{r} in the second element,
specifying linear constraints on the transition weight parameters \eqn{\alpha}.
The constraints are of the form \eqn{\alpha = R\xi + r}, where \eqn{R} is a known \eqn{(a\times l)}
constraint matrix of full column rank (\eqn{a} is the dimension of \eqn{\alpha}), \eqn{r} is a known \eqn{(a\times 1)} constant,
and \eqn{\xi} is an unknown \eqn{(l\times 1)} parameter. \strong{Alternatively}, set \eqn{R=0} to constrain the
weight parameters to the constant \eqn{r} (in this case, \eqn{\alpha} is dropped from the constrained parameter vector).}

\item{B_constraints}{a \eqn{(d \times d)} matrix with its entries imposing constraints on the impact matrix \eqn{B_t}:
\code{NA} indicating that the element is unconstrained, a positive value indicating strict positive sign constraint,
a negative value indicating strict negative sign constraint, and zero indicating that the element is constrained to zero.
Currently only available for models with \code{identification="heteroskedasticity"} or \code{"non-Gaussianity"} due to the
(in)availability of appropriate parametrizations that allow such constraints to be imposed.}

\item{standardize}{standardize the residuals to identity matrix covariance matrix?}

\item{structural_shocks}{If \code{TRUE}, returns structural shocks instead of residuals
(not available if \code{identification == "reduced_form"}, argument \code{standardize}
is if structural shocks are to be returned).}

\item{penalized}{Perform penalized LS estimation that minimizes penalized RSS in which estimates close to breaking or not satisfying the
usual stability condition are penalized? If \code{TRUE}, the tuning parameter is set by the argument \code{penalty_params[2]},
and the penalization starts when the eigenvalues of the companion form AR matrix are larger than \code{1 - penalty_params[1]}.}

\item{penalty_params}{a numeric vector with two positive elements specifying the penalization parameters:
the first element determined how far from the boundary of the stability region the penalization starts
(a number between zero and one, smaller number starts penalization closer to the boundary) and the second element
is a tuning parameter for the penalization (a positive real number, a higher value penalizes non-stability more).}

\item{allow_unstab}{If \code{TRUE}, estimates not satisfying the stability condition are allowed. Always \code{FALSE} if
\code{weight_function="relative_dens"}.}
}
\value{
Returns a \eqn{(T \times d)} matrix containing...
   \describe{
     \item{If \code{standardize == TRUE}:}{the standardized Pearson residuals.}
     \item{If \code{standardize == FALSE}:}{the nonstandardized residuals.}
     \item{If \code{structural_shocks == TRUE}:}{the structural shocks.}
   }
  Note that the starting time is \eqn{p + 1} counted from the beginning of the starting time of the data,
  as the first \eqn{p} observations are used as initial values.
}
\description{
\code{get_residuals} calculates residuals of a smooth transition VAR
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
\keyword{internal}
