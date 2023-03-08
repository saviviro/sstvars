#' @title sstvars: toolkit for structural smooth transition vector autoregressive models
#'
#' @description \code{sstvars} is a package for reduced form and structural smooth transition vector
#'   autoregressive models with focus on the structural analysis. It provides functions for ...
#'
#'   Some of the code is taken from the R package gmvarkit (also created by the author of this package)
#'   but is copied rather than using the code in gmvarkit directly in order to allow independent development.
#'
#'   The readme file is a good place to start and the vignette might be useful too.
#'
#' @docType package
#' @author you <savi.virolainen@helsinki.fi>
#' @import Rcpp RcppArmadillo parallel pbapply
#' @importFrom Rcpp evalCpp
#' @useDynLib sstvars
#' @name sstvars
NULL
