#' @describeIn STVAR Log-likelihood method
#' @param object object of class \code{'stvar'}.
#' @param ... currently not used.
#' @export
logLik.stvar <- function(object, ...) object$loglik


#' @describeIn STVAR residuals method to extract Pearson residuals
#' @inheritParams logLik.stvar
#' @export
residuals.stvar <- function(object, ...) {
  res <- object$residuals_std
  colnames(res) <- colnames(object$data)
  res
}


#' @describeIn STVAR summary method
#' @inheritParams logLik.stvar
#' @inheritParams print.stvar
#' @export
summary.stvar <- function(object, ..., digits=2, standard_error_print=FALSE) {
  stvar <- object
  structure(list(stvar=stvar,
                 abs_boldA_eigens=get_boldA_eigens(stvar),
                 omega_eigens=get_omega_eigens(stvar),
                 digits=digits,
                 standard_error_print=standard_error_print),
            class="stvarsum")
}
