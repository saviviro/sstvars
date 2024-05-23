#' @title sstvars: toolkit for reduced form and structural smooth transition vector autoregressive models
#'
#' @description \code{sstvars} is a package for reduced form and structural smooth transition vector
#'   autoregressive models. The package implements various transition weight functions, conditional distributions,
#'   identification methods, and parameter restrictions. The model parameters are estimated with the method of maximum
#'   likelihood by running multiple rounds of a two-phase estimation procedure in which a genetic algorithm is used
#'   to find starting values for a gradient based method. For evaluating the adequacy of the estimated models,
#'   \code{sstvars} utilizes residuals based diagnostics and provides functions for graphical diagnostics and for calculating
#'   formal diagnostic tests. \code{sstvars} also accommodates the estimation of linear impulse response functions, nonlinear
#'   generalized impulse response functions, and generalized forecast error variance decompositions. Further functionality includes
#'   hypothesis testing, plotting the profile log-likelihood functions about the estimate, simulation from STVAR processes,
#'   and forecasting, for example.
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
