#
wrap_multisets <- function(n, d) {
  # Laita kaksi accuracy vaihtoehtoa ja selvitä järkevä max dim?
  N <- choose(n + d - 1, d)
  if(N > 1e+8) {
    stop("Aborted since the large dimension required in the calculations may cause memory issues and likely take very long!")
  }
  get_multisets_Cpp(n=n, d=d, N=N)
}

#' @title Calculate the dth induced matrix of a square matrix
#'
#' @description \code{d_lift} the dth induced matrix of a square matrix
#'  by so-called d-lifting as described in Parrilo and Jadbabaie (2008).
#'
#' @param A a real square matrix
#' @param d a striclty positive integer specifying the order of the lift.
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
  # MODIFY THIS SO THAT THE MULTISETS ARE CREATED OUTSIDE, SINCE THE SAME MULTISETS ARE USED FOR MANY MATRICES?

  stopifnot(all_pos_ints(d))
  stopifnot(nrow(A) == ncol(A))
  n <- nrow(A)
  N <- choose(n + d - 1, d)
  if(N > 1e+8) {
    stop("Aborted since the large dimension required in the calculations may cause memory issues and likely take very long!")
  }
  all_multisets <- get_multisets_Cpp(n=n, d=d, N=N) # Each row is one multiset

  get_mu_multiset <- function(multiset) { # Calculate \mu(alpha) as described in Parrilo and Jadbabaie (2008), p. 14.
    # Product of the factorials of the multiplicies of the multiset
    prod(vapply(unique(multiset), function(elem) factorial(sum(multiset == elem)), numeric(1)))
  }
  all_mu_multisets <- apply(all_multisets, MARGIN=1, FUN=get_mu_multiset) # mu(multiset) calculated for all multisets

  A_lifted <- matrix(nrow=N, ncol=N)

  get_elements <- Vectorize(function(i, j) { # Used in obtaining subsets A(alpha, beta)
    A[alpha[j], beta[i]]
  })

  for(i1 in 1:N) { # i1 = row index
    for(i2 in 1:N) { # i2 = column index
       alpha <- all_multisets[i1,]
       beta <- all_multisets[i2,]
       A_submat <- outer(seq_len(d), seq_len(d), get_elements)
       print(A_submat)
    }
  }
}

A <- matrix(c(11, 21, 31, 12, 22, 32, 13, 23, 33), nrow=3)
A <- matrix(c(0.5, 0.2, -0.1, 0.6), nrow=2)
A <- matrix(c(11, 21, 12, 22), nrow=2)

A <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3)

a <- c(1, 1, 2)
b <- c(1, 2, 2)

A <- matrix(c(1, 2, 3, 4, 5, 6, 7, 8, 9), nrow = 3, ncol = 3)
a <- c(1, 3)
b <- c(2, 3)
# Initialize B with the same dimensions as the length of a and b



# Create a function to extract elements from A using a and b
get_elements <- Vectorize(function(i, j) {
  A[a[j], b[i]]
})

# Create matrix B using outer() function and the get_elements() function
B <- outer(1:length(a), 1:length(b), Vectorize(get_elements))

# # Below works for creating the submatrix A_submat but is slow
# A_submat <- matrix(0, nrow = length(alpha), ncol = length(beta))
# # Loop through the elements of a and b and populate matrix B
# for (i in 1:length(alpha)) {
#   for (j in 1:length(beta)) {
#     A_submat[i, j] <- A[alpha[j], beta[i]]
#   }
# }

