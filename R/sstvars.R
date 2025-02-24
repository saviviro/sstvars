#' @title sstvars: toolkit for reduced form and structural smooth transition vector autoregressive models
#'
#' @description \code{sstvars} is a package for reduced form and structural smooth transition vector
#'   autoregressive models. The package implements various transition weight functions, conditional distributions,
#'   identification methods, and parameter restrictions. The model parameters are estimated with the method of maximum
#'   likelihood or penalized maximum likelihood by running multiple rounds of either a two-phase estimation procedure
#'   or a three-phase procedure. In the former, a genetic algorithm is used to find starting values for a gradient based
#'   variable metric algorithm. In the latter, nonlinear least squares (NLS) first used obtain initial estimates for some
#'   of the parameters, then a genetic algorithm is used to find starting values for the rest of the parameters conditional
#'   on the NLS estimates, and finally a gradient based variable metric algorithm is initialized from the estimates obtained
#'   from the previous two steps. For evaluating the adequacy of the estimated models, \code{sstvars} utilizes residuals based
#'   diagnostics and provides functions for graphical diagnostics as well as for calculating formal diagnostic tests.
#'   \code{sstvars} also accommodates the estimation of linear impulse response functions, nonlinear generalized impulse response
#'   functions, and generalized forecast error variance decompositions. Further functionality includes hypothesis testing,
#'   plotting the profile log-likelihood functions about the estimate, simulation from STVAR processes, and forecasting, for example.
#'
#'   The vignette is a good place to start, and see also the readme file.
#'
#' @docType package
#' @author you <savi.virolainen@helsinki.fi>
#' @import Rcpp RcppArmadillo parallel pbapply
#' @importFrom Rcpp evalCpp
#' @useDynLib sstvars
#' @name sstvars-package
"_PACKAGE"
