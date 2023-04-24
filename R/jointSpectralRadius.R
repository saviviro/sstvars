
#' @title Calculate the dth induced matrix of a square matrix
#'
#' @description \code{d_lift} the dth induced matrix of a square matrix
#'  by so-called d-lifting as described in Parrilo and Jadbabaie (2008).
#'
#' @param A a real square matrix
#' @param d a strictly positive integer specifying the order of the lift.
#' @param multisets all \eqn{d}-element multisets in \eqn{\lbrace 1,...,nrow(A) \rbrace} in
#'   a lexicographic order (obtained with get_multisets_Cpp).
#' @details Calculated the d-lifting using the formula presented in Parrilo and Jadbabaie (2008),
#'   Equation~(8).
#' @return Returns the dth induced matrix of A: a square matrix with \code{choose(nrow(A) + d - 1, d)}
#'   rows and columns.
#' @references
#'  \itemize{
#'    \item P.A. Parrilo, A. Jadbabaie. 2008. Approximation of the joint spectral
#'       radius using sum of squares. \emph{Linear Algebra and its Applications},
#'       \strong{428}, 2385-2402.
#'  }
#' @keywords internal

d_lift <- function(A, d) {
  stopifnot(all_pos_ints(d))
  stopifnot(nrow(A) == ncol(A))
  n <- nrow(A)
  N <- choose(n + d - 1, d)
  if(N > 1e+6) {
    stop("Aborted since the large dimension required in the calculations may cause memory issues and likely take very long!")
  }
  all_multisets <- get_multisets_Cpp(n=n, d=d, N=N) # Each row is one multiset

  get_mu_multiset <- function(multiset) { # Calculate \mu(alpha) as described in Parrilo and Jadbabaie (2008), p. 14.
    # Product of the factorials of the multiplicies of the multiset
    prod(vapply(unique(multiset), function(elem) factorial(sum(multiset == elem)), numeric(1)))
  }
  all_mu_multisets <- apply(all_multisets, MARGIN=1, FUN=get_mu_multiset) # mu(multiset) calculated for all multisets

  A_lifted <- matrix(nrow=N, ncol=N) # Initialize A_lifted

  get_elements <- Vectorize(function(i, j) { # Used in obtaining "submatrices" A(alpha, beta) with alpha and beta multisets
    A[alpha[j], beta[i]] # Notice the indexing (j and i "swapped")
  })

  # Calculate the elements of A_lifted using Equation (8) of Parrilo and Jadbabaie (2008)
  for(i1 in 1:N) { # i1 = row index
    for(i2 in 1:N) { # i2 = column index
       alpha <- all_multisets[i1,] # get_element uses this
       beta <- all_multisets[i2,] # get_elements uses this
       A_submat <- outer(seq_len(d), seq_len(d), get_elements) # "Submatrix" A(alpha,beta)
       A_lifted[i1, i2] <-  get_permanent_Cpp(A_submat)/sqrt(all_mu_multisets[i1]*all_mu_multisets[i2]) # Equation (8)
    }
  }

  # # Below works also for creating the "submatrix" A_submat but is slow as uses loops
  # A_submat <- matrix(0, nrow = length(alpha), ncol = length(beta))
  # # Loop through the elements of a and b and populate matrix B
  # for (i in 1:length(alpha)) {
  #   for (j in 1:length(beta)) {
  #     A_submat[i, j] <- A[alpha[j], beta[i]]
  #   }
  # }

  A_lifted
}



#' @title Calculate upper bound for the joint spectral radius of the "bold A" matrices
#'
#' @description \code{bound_jsr} calculates an upper bound for the joint spectral radius of the "bold A" matrices
#' as described in Parrilo and Jadbabaie (2008), Theorems 4.2 and 4.3 and Equations (8) and (11).
#'
#' @param all_boldA all \eqn{((dp)x(dp))} "bold A" (companion form) matrices in a 3D array,
#'  so that \code{[, , m]} gives the matrix the regime \code{m}.
#' @param accuracy what should the relative accuracy of the bounds be? Note only upper bound is used
#'  and that the bound holds with any accuracy (but the bounds get tighter with increased accuracy).
#' @details Upper bound calculated the formula presented in Parrilo and Jadbabaie (2008), Equation~(11),
#'  whereas the accuracy is based on  Parrilo and Jadbabaie (2008), Table~1. Kheifets and Saikkonen (2020)
#'  show that if the joint spectral radius is smaller than one, the STVAR process is ergodic stationary.
#'  Therefore, if the upper bound is smaller than one, the process is stationary ergodic. However, as the condition
#'  is not necessary but sufficient and also because the bound might be too conservative, upper bound larger than
#'  one does not imply that the process is not ergodic stationary. You can try larger accuracy, and if the bound
#'  is still larger than one, the result does not tell whether the process is ergodic stationary or not.
#'
#'  Note that due to the extremely large dimenions required in the calculation of the bound, using large accuracy
#'  might not be feasible for other than very small models (small p*d^2). You can also try other implementations
#'  for bounding the joint spectral radius, for instance, the JSR toolbox in Matlab (Jungers 2023).
#' @return Returns an upper bound for the joint spectral radius of the "companion form AR matrices" of the regimes.
#' @references
#'  \itemize{
#'    \item I.L. Kheifets, P.J. Saikkonen. 2020. Stationarity and ergodicity of Vector STAR models.
#'      \emph{Econometric Reviews}, \strong{39}:4, 407-414.
#'    \item P.A. Parrilo, A. Jadbabaie. 2008. Approximation of the joint spectral
#'       radius using sum of squares. \emph{Linear Algebra and its Applications},
#'       \strong{428}, 2385-2402.
#'    \item R. Jungers (2023). The JSR toolbox (https://www.mathworks.com/matlabcentral/fileexchange/33202-the-jsr-toolbox),
#'       MATLAB Central File Exchange.
#'  }
#' @keywords internal

bound_jsr <- function(all_boldA, accuracy=c("0.707", "0.840", "0.917", "0.957", "0.978")) {
  # Determine "2d" (here d) based on the givas accuracy as in Parrilo & Jadbabaie (2008), Table 1.
  # Note that the upper bound holds at all accuracies; it just becomes tighter when accuracy is increased.
  accuracy <- match.arg(accuracy)
  if(accuracy == "0.707") {
    d <- 2
  } else if(accuracy == "0.840") {
    d <- 4
  } else if(accuracy == "0.917") {
    d <- 8
  } else if(accuracy == "0.957") {
    d <- 16
  } else { # accuracy == "0.978"
    d <- 32
  }
  n <- nrow(all_boldA[, , 1])
  N <- choose(n + d - 1, d)
  if(N > 1500) {
    cat("The calculations likely take relatively long due to the large dimension required!\n")
    cat("If it takes too long, you can try first a smaller accuracy.\n")
    cat("You may also try other implementations, e.g., the JSR toolbox in Matlab (Jungers 2023).")
  } else if(N > 3000) {
    cat("The large dimension required in the calculations may cause memory issues and likely take very long!\n")
    cat("You should probably try smaller accuracy or other implementations, e.g., the JSR toolbox in Matlab (Jungers 2023).")
  }
  all_multisets <- get_multisets_Cpp(n=n, d=d, N=N)

  M <- dim(all_boldA)[3]
  sum_of_d_lifted <- matrix(0, nrow=N, ncol=N)
  for(i1 in 1:M) {
    sum_of_d_lifted <- sum_of_d_lifted + d_lift(A=all_boldA[, , i1], d=d)
  }
  max(eigen(sum_of_d_lifted)$values)^(1/d) # Upper bound for the JSR as in Parrilo & Jadbabaie (2008), Equation (11)
}


#' @title Calculate upper bound for the joint spectral radius of the "companion form AR matrices" of the regimes
#'
#' @description \code{bound_joint_spectral_radius} calculates an upper bound for the joint spectral radius of the
#'  "companion form AR matrices" matrices of the regimes as described in Parrilo and Jadbabaie (2008),
#'   Theorems 4.2 and 4.3 and Equations (8) and (11).
#'
#' @inheritParams diagnostic_plot
#' @inheritParams bound_jsr
#' @inherit bound_jsr details references return
#' @examples
#' FILL IN EXAMPLES
#' @export


bound_joint_spectral_radius <- function(stvar, accuracy=c("0.707", "0.840", "0.917", "0.957", "0.978")) {
  accuracy <- match.arg(accuracy)
}
