
#' @title Vectorization operator
#'
#' @description \code{vec} stacks columns of the given matrix to form a vector.
#'
#' @param A a size \eqn{(dxd)} square matrix to be vectorized.
#' @return a vector of size \eqn{(d^2x1)}.
#' @section Warning:
#'  No argument checks!
#' @keywords internal

vec <- function(A) {
  as.vector(A)
}

#' @title Reverse vectorization operator
#'
#' @description \code{unvec} forms a square matrix from a vector of
#'  stacked columns, stacked by \code{vec}.
#'
#' @param a a size \eqn{(d^2x1)} vector to be unvectorized into a \eqn{(dxd)} matrix.
#' @param d the number of rows in the square matrix to be formed.
#' @return a matrix of size \eqn{(dxd)}.
#' @section Warning:
#'  No argument checks!
#' @keywords internal

unvec <- function(d, a) {
  matrix(a, nrow=d, byrow=FALSE)
}


#' @title Vectorization operator that removes zeros
#'
#' @description \code{Bvec} stacks columns of the given matrix to form a vector
#'   and removes elements that are zeros.
#'
#' @param Bmat a size \eqn{(dxd)} square matrix to be vectorized.
#' @return a vector of length \eqn{d^2 - n_zeros} where \eqn{n_zeros} is the
#'   number of zero entries in the matrix \code{B}.
#' @section Warning:
#'  No argument checks!
#' @keywords internal

Bvec <- function(Bmat) {
  Bmat[Bmat != 0]
}


#' @title Reverse vectorization operator that restores zeros
#'
#' @description \code{unBvec} forms a square matrix from a vector of
#'   stacked columns where zeros are removed according to impact
#'   matrix constraints.
#'
#' @inheritParams unvec
#' @param Bvector a length \eqn{d^2 - n_zeros} vector where \eqn{n_zeros} is the
#'   number of zero entries in the impact matrix.
#'
#' @return a \eqn{(d x d)} impact matrix \eqn{B}.
#' @section Warning:
#'  No argument checks!
#' @keywords internal

unBvec <- function(Bvector, d, B_constraints) {
  if(is.null(B_constraints)) stop("B_constraints are needed")
  new_B <- numeric(d^2)
  new_B[B_constraints != 0] <- Bvector
  matrix(new_B, nrow=d, byrow=FALSE)
}




#' @title Parsimonious vectorization operator for symmetric matrices
#'
#' @description \code{vech} stacks columns of the given matrix from
#'   the principal diagonal downwards (including elements on the diagonal) to form a vector.
#'
#' @param A a size \eqn{(dxd)} symmetric matrix to be vectorized parsimoniously.
#' @return a vector of size \eqn{(d(d+1)/2x1)}.
#' @section Warning:
#'  No argument checks!
#' @keywords internal

vech <- function(A) {
  A[lower.tri(x=A, diag=TRUE), drop=TRUE]
}


#' @title Reverse operator of the parsimonious vectorization operator \code{vech}
#'
#' @description \code{unvech} creates a symmetric matrix from the given vector by
#'   copying the lower triangular part to be the upper triangular part as well.
#'
#' @param a a size \eqn{(d(d+1)/2x1)} vector to be unvectorized into a symmetric \eqn{(dxd)} matrix.
#' @param d number of rows the square matrix to be formed.
#' @return a symmetric matrix of size \eqn{(dxd)}.
#' @section Warning:
#'  No argument checks!
#' @keywords internal

unvech <- function(d, a) {
  A <- matrix(nrow=d, ncol=d)
  upA <- upper.tri(A)
  A[!upA] <- a
  A[upA] <- t(A)[upA]
  A
}


#' @title Simultaneously diagonalize two covariance matrices
#'
#' @description \code{diag_Omegas} Simultaneously diagonalizes two covariance matrices using
#'   eigenvalue decomposition.
#'
#' @param Omega1 a positive definite \eqn{(dxd)} covariance matrix \eqn{(d>1)}
#' @param Omega2 another positive definite \eqn{(dxd)} covariance matrix
#' @details See the return value and Muirhead (1982), Theorem A9.9 for details.
#' @return Returns a length \eqn{d^2 + d} vector where the first \eqn{d^2} elements
#'   are \eqn{vec(W)} with the columns of \eqn{W} being (specific) eigenvectors of
#'   the matrix \eqn{\Omega_2\Omega_1^{-1}} and the rest \eqn{d} elements are the
#'   corresponding eigenvalues "lambdas". The result satisfies \eqn{WW' = Omega1} and
#'   \eqn{Wdiag(lambdas)W' = Omega2}.
#'
#'   If \code{Omega2} is not supplied, returns a vectorized symmetric (and pos. def.)
#'   square root matrix of \code{Omega1}.
#' @section Warning:
#'  No argument checks! Does not work with dimension \eqn{d=1}!
#' @references
#' \itemize{
#'   \item Muirhead R.J. 1982. Aspects of Multivariate Statistical Theory, \emph{Wiley}.
#' }
#' @examples
#' d <- 2
#' W0 <- matrix(1:(d^2), nrow=2)
#' lambdas0 <- 1:d
#' (Omg1 <- W0%*%t(W0))
#' (Omg2 <- W0%*%diag(lambdas0)%*%t(W0))
#' res <- diag_Omegas(Omg1, Omg2)
#' W <- matrix(res[1:(d^2)], nrow=d, byrow=FALSE)
#' tcrossprod(W) # == Omg1
#' lambdas <- res[(d^2 + 1):(d^2 + d)]
#' W%*%diag(lambdas)%*%t(W) # == Omg2
#' @export

diag_Omegas <- function(Omega1, Omega2) {
  eig1 <- eigen(Omega1, symmetric=TRUE)
  D <- diag(eig1$values) # Pos. def.
  H <- eig1$vectors # Orthogonal
  sqrt_omg1 <- H%*%sqrt(D)%*%t(H) # Symmetric and pos. def.
  if(missing(Omega2)) return(vec(sqrt_omg1))
  inv_sqrt_omg1 <- solve(sqrt_omg1)
  eig2 <- eigen(inv_sqrt_omg1%*%Omega2%*%inv_sqrt_omg1, symmetric=TRUE)
  lambdas <- eig2$values
  V <- eig2$vectors # Orthogonal
  W <- sqrt_omg1%*%V # all.equal(W%*%t(W), Omega1); all.equal(W%*%diag(lambdas)%*%t(W), Omega2)
  c(vec(W), lambdas)
}


#' @title In the decomposition of the covariance matrices (Muirhead, 1982, Theorem A9.9), change
#'   the order of the covariance matrices.
#'
#' @description \code{redecompose_Omegas} exchanges the order of the covariance matrices in
#'   the decomposition of Muirhead (1982, Theorem A9.9) and returns the new decomposition.
#'
#' @inheritParams is_stationary
#' @param W a length \code{d^2} vector containing the vectorized W matrix.
#' @param lambdas a length \code{d*(M-1)} vector of the form \strong{\eqn{\lambda_{2}}}\eqn{,...,}\strong{\eqn{\lambda_{M}}}
#'   where \strong{\eqn{\lambda_{m}}}\eqn{=(\lambda_{m1},...,\lambda_{md})}
#' @param perm a vector of length \code{M} giving the new order of the covarince matrices
#'   (relative to the current order)
#' @details We consider the following decomposition of positive definite covariannce matrices:
#'  \eqn{\Omega_1 = WW'}, \eqn{\Omega_m = W\Lambda_{m}W'}, \eqn{m=2,..,M} where \eqn{\Lambda_{m} = diag(\lambda_{m1},...,\lambda_{md})}
#'  contains the strictly postive eigenvalues of \eqn{\Omega_m\Omega_1^{-1}} and the column of the invertible \eqn{W} are the
#'  corresponding eigenvectors. Note that this decomposition does not necessarily exists for \eqn{M > 2}.
#'
#'  See Muirhead (1982), Theorem A9.9 for more details on the decomposition and the source code for more details on the reparametrization.
#' @return Returns a \eqn{d^2 + (M - 1)*d x 1} vector of the form \code{c(vec(new_W), new_lambdas)}
#'   where the lambdas parameters are in the regimewise order (first regime 2, then 3, etc) and the
#'   "new W" and "new lambdas" are constitute the new decomposition with the order of the covariance
#'   matrices given by the argument \code{perm}. Notice that if the first element of \code{perm}
#'   is one, the W matrix will be the same and the lambdas are just re-ordered.
#'
#'   \strong{Note that unparametrized zero elements ARE present in the returned W!}
#' @section Warning:
#'  No argument checks! Does not work with dimension \eqn{d=1} or with only
#'  one mixture component \eqn{M=1}.
#' @inherit diag_Omegas references
#' @examples
#'  d <- 2
#'  M <- 2
#'  Omega1 <- matrix(c(2, 0.5, 0.5, 2), nrow=d)
#'  Omega2 <- matrix(c(1, -0.2, -0.2, 1), nrow=d)
#'
#'  # Decomposition with Omega1 as the first covariance matrix:
#'  decomp1 <- diag_Omegas(Omega1, Omega2)
#'  W <- matrix(decomp1[1:d^2], nrow=d, ncol=d)
#'  lambdas <- decomp1[(d^2 + 1):length(decomp1)]
#'  tcrossprod(W) # = Omega1
#'  W%*%tcrossprod(diag(lambdas), W) # = Omega2
#'
#'  # Reorder the covariance matrices in the decomposition so that now
#'  # the first covariance matrix is Omega2:
#'  decomp2 <- redecompose_Omegas(M=M, d=d, W=as.vector(W), lambdas=lambdas,
#'                                perm=2:1)
#'  new_W <- matrix(decomp2[1:d^2], nrow=d, ncol=d)
#'  new_lambdas <- decomp2[(d^2 + 1):length(decomp2)]
#'  tcrossprod(new_W) # = Omega2
#'  new_W%*%tcrossprod(diag(new_lambdas), new_W) # = Omega1
#' @export

redecompose_Omegas <- function(M, d, W, lambdas, perm=1:sum(M)) {
  M <- sum(M)
  if(all(perm == 1:M)) {
    return(c(W, lambdas))
  } else if(perm[1] == 1) {
    new_lambdas <- matrix(lambdas, nrow=d, ncol=M - 1, byrow=FALSE)[, perm[2:M] - 1]
    return(c(W, new_lambdas))
  }
  # If the first covariance matrix changes, the decomposition needs to be reparametrized.
  lambdas <- cbind(rep(1, d), # Lambda matrix of the first regime is always identity matrix
                   matrix(lambdas, nrow=d, ncol=M - 1, byrow=FALSE)) # [, m], m=1, ..., M

  new_W <- matrix(W, nrow=d, ncol=d, byrow=FALSE)%*%diag(sqrt(lambdas[, perm[1]]))
  new_lambdas <- matrix(nrow=d, ncol=M - 1)
  for(i1 in 1:(M - 1)) {
    new_lambdas[, i1] <- diag(diag(1/lambdas[, perm[1]])%*%diag(lambdas[, perm[i1 + 1]]))
  }
  c(vec(new_W), vec(new_lambdas))
}

