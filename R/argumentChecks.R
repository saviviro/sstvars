#' @title Check the stability condition for each of the regimes
#'
#' @description \code{stab_conds_satisfied} checks whether the stability condition is satisfied
#'   for each of the regimes.
#'
#' @inheritParams pick_allA
#' @param all_boldA 3D array containing the \eqn{((dp)x(dp))} "bold A" (companion form) matrices of each regime,
#'   obtained from \code{form_boldA}. Will be computed if not given.
#' @param tolerance Returns \code{FALSE} if modulus of any eigenvalue of "bold A" is larger or equal to \code{1-tolerance}.
#' @inherit pick_phi0 details
#' @return Returns \code{TRUE} if the stability condition is satisfied for all regimes and \code{FALSE} if not.
#'   According to the argument \code{tolerance}, \code{stab_conds_satisfied} may return \code{FALSE} when the parameter
#'   vector satisfies the stability conditions but is very close to the boundary (this is used to ensure numerical stability
#'   in the estimation of the model parameters).
#' @section Warning:
#'  No argument checks!
#' @inherit pick_phi0 details
#' @inherit form_boldA references
#' @keywords internal

stab_conds_satisfied <- function(p, M, d, params, all_boldA=NULL, tolerance=1e-3) {
  if(is.null(all_boldA)) {
    all_A <- pick_allA(p=p, M=M, d=d, params=params)
    all_boldA <- form_boldA(p=p, M=M, d=d, all_A=all_A)
  }
  for(m in 1:M) {
    if(any(abs(eigen(all_boldA[, , m], symmetric=FALSE, only.values=TRUE)$'values') >= 1 - tolerance)) {
      return(FALSE)
    }
  }
  TRUE
}
