
#' @title Vectorization operator
#'
#' @description \code{vec} stacks columns of the given matrix to form a vector.
#'
#' @param A a size \eqn{(dxd)} square matrix to be vectorized.
#' @return a vector of size \eqn{(d^2x1)}.
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
#' @keywords internal

unvec <- function(d, a) {
  matrix(a, nrow=d, byrow=FALSE)
}


#' @title Vectorization operator that removes zeros
#'
#' @description \code{Wvec} stacks columns of the given matrix to form a vector
#'   and removes elements that are zeros.
#'
#' @param W a size \eqn{(dxd)} square matrix to be vectorized.
#' @return a vector of length \eqn{d^2 - n_zeros} where \eqn{n_zeros} is the
#'   number of zero entries in the matrix \code{W}.
#' @keywords internal

Wvec <- function(W) {
  W[W != 0]
}


#' @title Reverse vectorization operator that restores zeros
#'
#' @description \code{unWvec} forms a square matrix from a vector of
#'   stacked columns where zeros are removed according to impact
#'   matrix constraints.
#'
#' @inheritParams unvec
#' @param Wvector a length \eqn{d^2 - n_zeros} vector where \eqn{n_zeros} is the
#'   number of zero entries in the impact matrix.
#'
#' @return a \eqn{(d x d)} impact matrix \eqn{W}.
#' @keywords internal

unWvec <- function(Wvector, d, B_constraints) {
  new_W <- numeric(d^2)
  new_W[B_constraints != 0] <- Wvector
  matrix(new_W, nrow=d, byrow=FALSE)
}




#' @title Parsimonious vectorization operator for symmetric matrices
#'
#' @description \code{vech} stacks columns of the given matrix from
#'   the principal diagonal downwards (including elements on the diagonal) to form a vector.
#'
#' @param A a size \eqn{(dxd)} symmetric matrix to be vectorized parsimoniously.
#' @return a vector of size \eqn{(d(d+1)/2x1)}.
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
#' # Create two (2x2) coviance matrices using the parameters W and lambdas:
#' d <- 2 # The dimension
#' W0 <- matrix(1:(d^2), nrow=2) # W
#' lambdas0 <- 1:d # The eigenvalues
#' (Omg1 <- W0%*%t(W0)) # The first covariance matrix
#' (Omg2 <- W0%*%diag(lambdas0)%*%t(W0)) # The second covariance matrix
#'
#' # Then simultaneously diagonalize the covariance matrices:
#' res <- diag_Omegas(Omg1, Omg2)
#'
#' # Recover W:
#' W <- matrix(res[1:(d^2)], nrow=d, byrow=FALSE)
#' tcrossprod(W) # == Omg1, the first covariance matrix
#'
#' # Recover lambdas:
#' lambdas <- res[(d^2 + 1):(d^2 + d)]
#' W%*%diag(lambdas)%*%t(W) # == Omg2, the second covariance matrix
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
#'   the ordering of the covariance matrices.
#'
#' @description \code{redecompose_Omegas} exchanges the order of the covariance matrices in
#'   the decomposition of Muirhead (1982, Theorem A9.9) and returns the new decomposition.
#'
#' @param M the number of regimes in the model
#' @param d the number of time series in the system
#' @param W a length \code{d^2} vector containing the vectorized W matrix.
#' @param lambdas a length \code{d*(M-1)} vector of the form \strong{\eqn{\lambda_{2}}}\eqn{,...,}\strong{\eqn{\lambda_{M}}}
#'   where \strong{\eqn{\lambda_{m}}}\eqn{=(\lambda_{m1},...,\lambda_{md})}
#' @param perm a vector of length \code{M} giving the new order of the covariance matrices
#'   (relative to the current order)
#' @details We consider the following decomposition of positive definite covariannce matrices:
#'  \eqn{\Omega_1 = WW'}, \eqn{\Omega_m = W\Lambda_{m}W'}, \eqn{m=2,..,M} where
#'  \eqn{\Lambda_{m} = diag(\lambda_{m1},...,\lambda_{md})} contains the strictly postive eigenvalues of
#'  \eqn{\Omega_m\Omega_1^{-1}} and the column of the invertible \eqn{W} are the corresponding eigenvectors.
#'  Note that this decomposition does not necessarily exists for \eqn{M > 2}.
#'
#'  See Muirhead (1982), Theorem A9.9 for more details on the decomposition and the source code for more details on the
#'  reparametrization.
#' @return Returns a \eqn{d^2 + (M - 1)d \times 1} vector of the form \code{c(vec(new_W), new_lambdas)}
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
#' # Create two (2x2) coviance matrices:
#' d <- 2 # The dimension
#' M <- 2 # The number of covariance matrices
#' Omega1 <- matrix(c(2, 0.5, 0.5, 2), nrow=d)
#' Omega2 <- matrix(c(1, -0.2, -0.2, 1), nrow=d)
#'
#' # The decomposition with Omega1 as the first covariance matrix:
#' decomp1 <- diag_Omegas(Omega1, Omega2)
#' W <- matrix(decomp1[1:d^2], nrow=d, ncol=d) # Recover W
#' lambdas <- decomp1[(d^2 + 1):length(decomp1)] # Recover lambdas
#' tcrossprod(W) # = Omega1
#' W%*%tcrossprod(diag(lambdas), W) # = Omega2
#'
#' # Reorder the covariance matrices in the decomposition so that now
#' # the first covariance matrix is Omega2:
#' decomp2 <- redecompose_Omegas(M=M, d=d, W=as.vector(W), lambdas=lambdas,
#'                               perm=2:1)
#' new_W <- matrix(decomp2[1:d^2], nrow=d, ncol=d) # Recover W
#' new_lambdas <- decomp2[(d^2 + 1):length(decomp2)] # Recover lambdas
#' tcrossprod(new_W) # = Omega2
#' new_W%*%tcrossprod(diag(new_lambdas), new_W) # = Omega1
#' @export

redecompose_Omegas <- function(M, d, W, lambdas, perm=1:M) {
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


#' @title Calculate symmetric square root matrix of a positive definite covariance
#'   matrix
#'
#' @description \code{get_symmetric_sqrt} calculates symmetric square root matrix
#'  of a positive definite covariance matrix
#'
#' @param Omega a positive definite covariance matrix
#' @return a vectorized symmetric square root matrix of the matrix \code{Omega}.
#' @section Warning:
#'  No argument checks!
#' @keywords internal

get_symmetric_sqrt <- function(Omega) {
  diag_Omegas(Omega1=Omega)
}


#' @title Compute the j:th power of a square matrix A
#'
#' @description \code{mat_power} computes the j:th power of a square matrix A using
#' exponentiation by squaring.
#'
#' @param A A square numeric matrix.
#' @param j A natural number representing the power to which the matrix will be raised.
#'
#' @return A matrix which is A raised to the power j.
#' @keywords internal

mat_power <- function(A, j) {
  if(j == 0) {
    return(diag(nrow(A)))
  }
  res <- diag(nrow(A))
  tmp <- A
  while(j > 0) {
    if(j%%2 == 1) {
      res <- res%*%tmp
    }
    tmp <- tmp%*%tmp
    j <- j%/%2
  }
  res
}


#' @title Create a special matrix J
#'
#' @description \code{create_J_matrix} generates a d x dp matrix J, where the first d x d block is the identity matrix I_d,
#' and the rest is filled with zeros.
#'
#' @param d An integer representing the dimension of the identity matrix.
#' @param p An integer representing the factor by which to extend the matrix with zeros.
#'
#' @return A \eqn{d x dp} matrix \eqn{J} where the first \eqn{d x d} block is the identity matrix \eqn{I_d},
#'   and the rest is filled with zeros.
#' @keywords internal

create_J_matrix <- function(d, p) {
  J <- matrix(0, nrow=d, ncol=d*p) # Initialize a matrix with zeros
  diag(J[ , 1:d]) <- 1 # Fill the first d x d block with the identity matrix
  J
}


#' @title Create Matrix F_i
#'
#' @description \code{create_Fi_matrix} function generates a \eqn{(T_obs \times T_obs)} matrix
#' with 1's on the \code{(i+1)}-th sub-diagonal and 0's elsewhere.
#'
#' @param i Integer, the lag for which the matrix is constructed, corresponding to the sub-diagonal filled with 1's.
#' @param T_obs Integer, the number of time periods or observations in the dataset, corresponding to the dimensions of the matrix.
#' @details Used in \code{Portmanteau_test}. The matrix \eqn{F_i} in Lütkepohl (2005), p.158.
#' @return A \eqn{(T_obs \times T_obs)} matrix with 1's on the \code{(i+1)}-th sub-diagonal and 0's elsewhere.
#' @references
#'  \itemize{
#'    \item Lütkepohl, H. (2005). New Introduction to Multiple Time Series Analysis. Springer.
#'  }
#' @keywords internal

create_Fi_matrix <- function(i, T_obs) {
  Fi <- matrix(0, nrow=T_obs, ncol=T_obs)  # Create a T_obs x T_obs matrix of zeros
  if(i == 0) {
    return(diag(T_obs))
  } else if(i < T_obs) {
    Fi[cbind((i+1):T_obs, 1:(T_obs-i))] <- 1 # Place 1s on the (i+1)-th sub-diagonal
  } # Matrix of zeros if i=T_obs.
  Fi
}


#' @title Reorder columns of a square matrix so that the first nonzero elements are
#'  in decreasing order
#'
#' @description \code{order_B} takes a square matrix \code{B} as input and reorders its columns
#' so that the diagonal entries of \code{B} are in decreasing order. The function
#' is designed to be computationally efficient and ensures that the input is a square matrix.
#'
#' @param B A square numeric matrix.
#' @return A square matrix \code{B} with columns reordered so that its
#' diagonal entries are in a decreasing order.
#' @keywords internal

order_B <- function(B) {
  # Change signs of columns so that the first element of each column is positive:
  for(i1 in 1:ncol(B)) {
    if(B[1, i1] < 0) {
      B[, i1] <- -B[, i1]
    }
  }

  # Reorder the columns of B so that the diagonal elements are in decreasing order:
  B[, order(B[1, ], decreasing=TRUE), drop=FALSE]
}
