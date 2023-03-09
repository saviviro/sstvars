#' @title Create a class 'stvar' object defining a reduced form or structural smooth transition VAR model
#'
#' @description \code{STVAR} creates a class \code{'stvar'} object that defines
#'  a reduced form or structural smooth transition VAR model
#'
#' @inheritParams loglikelihood
#' @param data a matrix or class \code{'ts'} object with \code{d>1} columns. Each column is taken to represent
#'  a single times series. \code{NA} values are not supported. Ignore if defining a model without data is desired.
#' @param d number of times series in the system, i.e. \code{ncol(data)}. This can be
#'   used to define STVAR models without data and can be ignored if \code{data} is provided.
#' @param calc_std_errors should approximate standard errors be calculated?
#' @details If data is provided, then also residuals are computed and included in the returned object.
#' @return Returns an object of class \code{'stvar'} defining the specified reduced form or structural
#'  smooth transition VAR model. Can be used to work with other functions provided in \code{sstvars}.
#' @section About S3 methods:
#'   If data is not provided, only the \code{print} and \code{simulate} methods are available.
#'   If data is provided, then in addition to the ones listed above, \code{predict} method is also available.
#'   See \code{?simulate.stvar} and \code{?predict.stvar} for details about the usage.
#' @seealso \code{\link{fitSTVAR}}, \code{\link{add_data}}, \code{\link{swap_parametrization}}, \code{\link{GIRF}}
#' @references
#'  \itemize{
#'    \item TO BE FILLED IN
#'  }
#' @examples
#' # FILL IN
#' @export

STVAR <- function(data, p, M, d, params, weight_function=c("relative_dens", "logit"), cond_dist=c("Gaussian", "Student"),
                  parametrization=c("intercept", "mean"), identification=c("reduced_form", "recursive", "heteroskedasticity"),
                  AR_constraints=NULL, mean_constraints=NULL, B_constraints=NULL, calc_std_erros=FALSE) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  identification <- match.arg(identification)
  if(missing(data) & missing(d)) stop("data or d must be provided")
  if(missing(data) || is.null(data)) {
    data <- NULL
  } else {
    data <- check_data(data=data, p=p)
    if(missing(d)) {
      d <- ncol(data)
    } else if(ncol(data) != d) {
      warning("ncol(data) does not equal d. Using d = ncol(data)")
      d <- ncol(data)
    }
  }
  check_pMd(p=p, M=M, d=d)
  # CHECK CONSTRAINTS HERE
  check_params(p=p, M=M, d=d, params=params, weight_function=weight_function, cond_dist=cond_dist,
               parametrization=parametrization, identification=identification, AR_constraints=AR_constraints,
               mean_constraints=mean_constraints, B_constraints=B_constraints)
  npars <- n_params(p=p, M=M, d=d, params=params, weight_function=weight_function, cond_dist=cond_dist,
                    parametrization=parametrization, identification=identification, AR_constraints=AR_constraints,
                    mean_constraints=mean_constraints, B_constraints=B_constraints)

  if(is.null(data)) {
    lok_and_tw <- list(loglik=NA, mw=NA)
    IC <- data.frame(AIC=NA, HQIC=NA, BIC=NA)
    residuals <- NA
  } else {
    if(npars >= d*nrow(data)) warning("There are at least as many parameters in the model as there are observations in the data")
    lok_and_ts <- loglikelihood(p=p, M=M, d=d, params=params, weight_function=weight_function, cond_dist=cond_dist,
                                parametrization=parametrization, identification=identification, AR_constraints=AR_constraints,
                                mean_constraints=mean_constraints, B_constraints=B_constraints, to_return="loglik_and_tw",
                                minval=NA)
    residuals <- NA # implement function get_residuals
    IC <- NA # implement function get_IC (divide by the number of observatios)
  }
}
