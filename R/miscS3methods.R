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
summary.stvar <- function(object, ..., digits=2) {
  stvar <- object
  structure(list(stvar=stvar,
                 abs_boldA_eigens=get_boldA_eigens(stvar),
                 omega_eigens=get_omega_eigens(stvar),
                 regime_means=get_regime_means(p=stvar$model$p, M=stvar$model$M, d=stvar$model$d,
                                               weight_function=stvar$model$weight_function,
                                               cond_dist=stvar$model$cond_dist,
                                               parametrization=stvar$model$parametrization,
                                               identification=stvar$model$identification,
                                               AR_constraints=stvar$model$AR_constraints,
                                               mean_constraints=stvar$model$mean_constraints,
                                               B_constraints=stvar$model$B_constraints),
                 digits=digits),
            class="stvarsum")
}
