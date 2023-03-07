#' @title Two-phase maximum likelihood estimation of a reduced form smooth transition VAR model
#'
#' @description \code{fitSTVAR} estimates a reduced form smooth transition VAR model in two phases:
#'   in the first phase, it uses a genetic algorithm to find starting values for a gradient based
#'   variable metric algorithm, which it then uses to finalize the estimation in the second phase.
#'   Parallel computing is utilized to perform multiple rounds of estimations in parallel.
#'
#' @inheritParams GAfit
#' @param nrounds the number of estimation rounds that should be performed.
#' @param ncores the number CPU cores to be used in parallel computing.
#' @param maxit the maximum number of iterations in the variable metric algorithm.
#' @param seeds a length \code{nrounds} vector containing the random number generator seed
#'  for each call to the genetic algorithm, or \code{NULL} for not initializing the seed.
#' @param print_res should summaries of estimation results be printed?
#' @param ... additional settings passed to the function \code{GAfit} employing the genetic algorithm.
#' @details
#'  If you wish to estimate a structural model, estimate first the reduced form model and then use the
#'  use the function FILL IN to estimate the structural model based on the estimated reduced form model.
#'
#'  Because of complexity and high multimodality of the log-likelihood function, it is \strong{not certain}
#'  that the estimation algorithm will end up in the global maximum point. It is expected that most of the
#'  estimation rounds will end up in some local maximum or saddle point instead. Therefore, a (sometimes large)
#'  number of estimation rounds is required for reliable results. Because of the nature of the model,
#'  the estimation may fail especially in the cases where the number of regimes is chosen too large.
#'
#'  The estimation process is computationally heavy and it might take considerably long time for large models with
#'  large number of observations. If the iteration limit \code{maxit} in the variable metric algorithm is reached,
#'  one can continue the estimation by iterating more with the function \code{iterate_more}. Alternatively, one may
#'  use the found estimates as starting values for the genetic algorithm and and employ another round of estimation
#'  (see \code{?GAfit} how to set up an initial population with the dot parameters).
#'
#'  \strong{If the estimation algorithm fails to create an initial population for the genetic algorithm},
#'  it usually helps to scale the individual series so that the AR coefficients (of a VAR model) will be
#'  relative small, preferably less than one. Even if one is able to create an initial population, it should
#'  be preferred to scale the series so that most of the AR coefficients will not be very large, as the
#'  estimation algorithm works better with relatively small AR coefficients. If needed, another package can be used
#'  to fit linear VARs to the series to see which scaling of the series results in relatively small AR coefficients.
#'
#' @return Returns an object of class \code{'stvar'} defining the estimated reduced form smooth transition VAR model.
#' @section S3 methods:
#'   The following S3 methods are supported for class \code{'stvar'}: \code{logLik}, \code{residuals}, \code{print}, \code{summary},
#'    \code{predict}, \code{simulate}, and \code{plot}. NONE OF THESE IS IMPLEMENTED YET!
#' @seealso \code{\link{GAfit}}
#' @references
#'  \itemize{
#'    \item FILL IN REFEENCES
#'  }
#' @examples
#' \donttest{
#' ## These are long running examples that use parallel computing!
#' # Running all the below examples will take approximately FILL IN HOW MANY minutes.
#'
#' # FILL IN EXAMPLES
#' }
#' @export

fitSTVAR <- function(data, p, M, weight_function=c("relative_dens", "logit"), cond_dist=c("Gaussian", "Student"),
                     parametrization=c("intercept", "mean"), AR_constraints=NULL, mean_constraints=NULL,
                     ncalls=(M + 1)^5, ncores=2, maxit=1000, seeds=NULL, print_res=TRUE, ...) {
}
