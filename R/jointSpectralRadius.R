
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
  if(N > 10000) {
    stop("Aborted since the large dimension required in the calculations may cause memory issues and they likely take very long!")
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
       A_lifted[i1, i2] <- get_permanent_Cpp(A_submat)/sqrt(all_mu_multisets[i1]*all_mu_multisets[i2]) # Equation (8)
    }
  }

  # # Below works also for creating the "submatrix" A_submat but is slow as uses loops
  # A_submat <- matrix(0, nrow=length(alpha), ncol=length(beta))
  # # Loop through the elements of a and b and populate matrix B
  # for(i in 1:length(alpha)) {
  #   for(j in 1:length(beta)) {
  #     A_submat[i, j] <- A[alpha[j], beta[i]]
  #   }
  # }

  A_lifted
}



#' @title Calculate upper bound for the joint spectral radius of the "bold A" matrices
#'
#' @description \code{bound_jsr_JP} calculates an upper bound for the joint spectral radius of the "bold A" matrices
#' as described in Parrilo and Jadbabaie (2008), Theorems 4.2 and 4.3 and Equations (8) and (11).
#'
#' @param all_boldA all \eqn{((dp)x(dp))} "bold A" (companion form) matrices in a 3D array,
#'  so that \code{[, , m]} gives the matrix the regime \code{m}.
#' @param accuracy what should the relative accuracy of the bounds be? Note only upper bound is used
#'  and that the bound holds with any accuracy (but the bounds get tighter with increased accuracy).
#' @details Upper bound calculated the formula presented in Parrilo and Jadbabaie (2008), Equation~(11),
#'  whereas the accuracy is based on  Parrilo and Jadbabaie (2008), Table 1. Note that the accuracy is
#'  for the bounds including the lower bound which we do not calculate, as only upper bound used in practice.
#'  Specifically, Kheifets and Saikkonen (2020) show that if the joint spectral radius is smaller than one,
#'  the STVAR process is ergodic stationary. Therefore, if the upper bound is smaller than one,
#'  the process is stationary ergodic. However, as the condition is not necessary but sufficient and also because
#'  the bound might be too conservative, upper bound larger than one does not imply that the process is not ergodic
#'  stationary. You can try higher accuracy, and if the bound is still larger than one, the result does not tell
#'  whether the process is ergodic stationary or not.
#'
#'  Note that due to the extremely large dimenions required in the calculation of the bound, using high accuracy
#'  might not be feasible for other than very small models (small p*d^2), as one may run out of memory or the
#'  computations just take very long. \strong{That is, it is advisable to start with low accuracy,
#'  and only increase it if the sufficient condition for ergodic stationarity is not satisfied in the upper bound,
#'  i.e., if the bound is larger than one.}
#'  You can also try other implementations for bounding the joint spectral radius, for instance,
#'  the JSR toolbox in Matlab (Jungers 2023).
#' @return Returns an upper bound for the joint spectral radius of the "companion form AR matrices" of the regimes.
#' @references
#'  \itemize{
#'    \item I.L. Kheifets, P.J. Saikkonen. 2020. Stationarity and ergodicity of Vector STAR models.
#'      \emph{Econometric Reviews}, \strong{39}:4, 407-414.
#'    \item P.A. Parrilo, A. Jadbabaie. 2008. Approximation of the joint spectral
#'       radius using sum of squares. \emph{Linear Algebra and its Applications}, \strong{428}, 2385-2402.
#'    \item R. Jungers (2023). The JSR toolbox (https://www.mathworks.com/matlabcentral/fileexchange/33202-the-jsr-toolbox),
#'       MATLAB Central File Exchange.
#'  }
#' @keywords internal

bound_jsr_JP <- function(all_boldA, accuracy=c("0.707", "0.840", "0.917", "0.957", "0.978")) {
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
  if(N > 28) {
    cat("The calculations likely take relatively long due to the large dimension required!\n")
    cat("If it takes too long, you can try first a lower accuracy.\n")
    cat("You may also try other implementations, e.g., the JSR toolbox in Matlab (Jungers 2023).\n")
  } else if(N > 34) {
    cat("The large dimension required in the calculations takes very long and migh cause memory issues!\n")
    cat("You should probably try lower accuracy or other implementations, e.g., the JSR toolbox in Matlab (Jungers 2023).\n")
  }
  all_multisets <- get_multisets_Cpp(n=n, d=d, N=N)

  M <- dim(all_boldA)[3]
  all_d_lifted <- array(NA, dim=c(N, N, M))
  for(i1 in 1:M) {
    all_d_lifted[, , i1] <- d_lift(A=all_boldA[, , i1], d=d)
  }
  sum_of_d_lifted <- apply(all_d_lifted, MARGIN=1:2, FUN=sum)
  max(abs(eigen(sum_of_d_lifted)$values))^(1/d) # Upper bound for the JSR as in Parrilo & Jadbabaie (2008), Equation (11)
}


#' @title Calculate upper bound for the joint spectral radius of a set of matrices
#'
#' @description \code{bound_jsr_G} calculates lowr and upper bounds for the joint spectral radious of a set of square matrices,
#'  typically the "bold A" matrices, using the algorithm by Gripenberg (1996)
#'
#' @param S the set of matrices the bounds should be calculated for in an array, in VAR applications,
#'  all \eqn{((dp)x(dp))} "bold A" (companion form) matrices in a 3D array, so that \code{[, , m]} gives the matrix
#'  the regime \code{m}.
#' @param epsilon a strictly positive real number specifying the absolute accuracy of the bounds (which will be approximate
#'  in practice). Smaller number gives better accuracy but requires more substantial computational effort.
#' @details The bounds are calculated using the Gripenberg's (1996) branch-and-bound method, which is also discussed
#'  in Chand and Blondel (2013). Specifically, Kheifets and Saikkonen (2020) show that if the joint spectral radius
#'  of the companion form AR matrices of the regimes is smaller than one, the STVAR process is ergodic stationary. Therefore,
#'  if the upper bound is smaller than one, the process is stationary ergodic. However, as the condition is not
#'  necessary but sufficient and also because the bound might be too conservative, upper bound larger than one
#'  does not imply that the process is not ergodic stationary. You can try higher accuracy, and if the bound is
#'  still larger than one, the result does not tell whether the process is ergodic stationary or not.
#'
#'  Note that with high precision (small \code{epsilon}), the computational effort required are substantial and
#'  the estimation may take very long.
#'
#'  You can also try other implementations for bounding the joint spectral radius, for instance,
#'  the JSR toolbox in Matlab (Jungers 2023).
#' @return Returns an upper bound for the joint spectral radius of the "companion form AR matrices" of the regimes.
#' @references
#'  \itemize{
#'  \item C-T Chang and V.D. Blondel. 2013 . An experimental study of approximation algorithms for the joint spectral radius.
#'      \emph{Numerical algorithms}, \strong{64}, 181-202.
#'    \item Gripenberg, G. 1996. Computing the joint spectral radius. \emph{Linear Algebra and its Applications},
#'      234, 43–60.
#'    \item I.L. Kheifets, P.J. Saikkonen. 2020. Stationarity and ergodicity of Vector STAR models.
#'      \emph{Econometric Reviews}, \strong{39}:4, 407-414.
#'    \item R. Jungers (2023). The JSR toolbox (https://www.mathworks.com/matlabcentral/fileexchange/33202-the-jsr-toolbox),
#'       MATLAB Central File Exchange.
#'  }
#' @keywords internal

bound_jsr_G <- function(S, epsilon=0.01, print_progress=TRUE) {
  n <- dim(S)[1] # The dimension of the n x n matrices
  m <- dim(S)[3] # The number of matrices
  maxit <- 1000 # Maximum number of iterations (just some large number)
  all_alpha <- array(NA, dim=maxit) # Storage for the lower bound in each iteration
  all_beta <- array(NA, dim=maxit)  # Storage for the upper bound in each iteration

  # Function to calculate the Fobelius norm, which we employ as our norm in the algorithm:
  Fobelius_norm <- function(A) sqrt(sum(A^2))

  ### Step 1, initialize

  # Calculate the initial lower and upper bounds for the joint spectral radius
  all_alpha[1] <- max(vapply(1:m, function(i1) max(abs(eigen(S[, , i1])$values)),
                             numeric(1))) # Max of the spectral radiusses of the matrices in S
  all_beta[1] <- max(vapply(1:m, function(i1) Fobelius_norm(S[, , i1]), numeric(1))) # Max of the Fobelius norms of the matrices in S

  ### Step 2, iteration process

  # Function to calculate the mu(MP): takes in the set of matrices that are involved in the concerned matrix product.
  # in the order they appear in the matrix product
  mu <- function(MP, k) { # Calculate the mu(MP) for the set of matrices MP in an [n, n, k] array, for iteration k>2
    all_prods <- array(NA, dim=c(n, n, k))
    res <- numeric(k)
    all_prods[, , 1] <- MP[, , 1]
    res[1] <- Fobelius_norm(MP[, , 1])
    for(i1 in 2:k) {
      all_prods[, , i1] <- all_prods[, , i1-1]%*%MP[, , i1]
      res[i1] <- Fobelius_norm(all_prods[, , i1])^(1/i1)
    }
    min(res)
  }

  mat_prod <- function(P) { # Calculate the product of a set of matrices in P in an [n, n, k] array
    all_prods <- array(NA, dim=c(n, n, dim(P)[3]))
    all_prods[, , 1] <- P[, , 1]
    for(i1 in 2:dim(P)[3]) {
      all_prods[, , i1] <- all_prods[, , i1-1]%*%P[, , i1]
    }
    all_prods[, , dim(P)[3]]
  }


  # First iteration, k=1:
  # The set MP is just the set S and the "products" are the single matrices in S.
  all_mu <- vapply(1:m, function(i1) Fobelius_norm(S[, , i1]), numeric(1)) # Calculate the mu(S)
  all_n_candidates <- array(NA, dim=maxit)  # The number of candidate matrix products in each iteration
  all_n_candidates[1] <- m # The number of candidate matrix products in the first iteration is the number of matrices in S
  which_new_candidates <- 1:m # The indices of the new candidate matrix products in the first iteration
  all_matprod_inds_old <- matrix(1:2, ncol=2) # Initialize matrix for storing the indices of the matrices involved in the matrix products
  # all_matprods_old <- array(S, dim=c(n, n, m, 1)) # All matrix products used in the previous iteration are just the matrices in S

  for(k in 2:maxit) {
    # For each iteration, calculate the set MP, i.e., the sets of matrices that are involved in the concerned matrix products.
    # The set MP is the set of all possible matrix products obtained by pre-multiplying a matrix product from the previous set of
    # candidate solutions by a matrix from the initial set S. The matrix products involved will consist of k matrices each.

    # Initialize matrix for storing the indices of the matrices involved in the matrix products:
    # Each column contains a matrix product given by the indices of the matrices in the set S
    n_new_matprods <- m*all_n_candidates[k-1]
    all_matprod_inds <- matrix(nrow=k, ncol=n_new_matprods)  # [matrix_indices, index_of_the_product]
    all_matprod_inds[2:k,] <- all_matprod_inds_old[,which_new_candidates] # Copy the indices from the previous iteration: repeats m times
    all_matprod_inds[1,] <- c(vapply(1:m, function(i1) rep(i1, times=all_n_candidates[k-1]),
                                     numeric(all_n_candidates[k-1]))) # Add the indices of the matrices from the initial set S
    # Note that due to premultiplication, the first row is for the new matrices in the product.

    # # Matrix products from the previous iteration and then we add the new matrix products by adding one pre-multiplication
    # # by a matrix from the initial set S to a matrix product from the previous iteration.
    # all_matprods <- array(dim=c(n, n, n_new_matprods, k)) # Initialize array for storing the matrix products:
    # # [n, n, index_of_the_product, product_up_to_this_row_in_all_mat_prod_inds_from_the_bottom_to_the_top]
    # # (above, multiplication is still done in the order from the first row to the last row)
    # for(i1 in 2:k) { # Fill in the old matrix products
    #   # i1=1 is not included because it is for the new matrix products
    #   # all_matprods_old[, , , 1] corresponds to the new matrix matrix product of the previous iteration, so it
    #   # will be filled to the second super-slice, etc.
    #   all_matprods[, , , i1] <- all_matprods_old[, , , i1-1] # Repeats the array m times to fill the bigger sub-array
    # }
    # # Calculate and fill in the new matrix products:
    # for(i1 in 1:n_new_matprods) {
    #   all_matprods[, , i1, 1] <- S[, , all_matprod_inds[1, i1]]%*%all_matprods[, , i1, 2]
    # }

    # Calculate the mu(MP) for the set of matrices MP in an [n, n, k] array, for iteration k>2
    all_mu <- numeric(n_new_matprods) # Initialize storage for the mu(MP) values

    for(i1 in 1:n_new_matprods) { # The set MP is the matrix products in all_matprods[, , , 1]
      all_mu[i1] <- mu(S[, , all_matprod_inds[,i1]], k=k) # S[, , all_matprod_inds[,i1]] obtains the matrices in the matrix product
    }

    # Construct the new set of candidates MP, i.e., the products for which mu(MP) > alpha_{k-1} + epsilon
    which_new_candidates <- which(all_mu > all_alpha[k-1] + epsilon)
    if(length(which_new_candidates) == 0) {
      # If there are no new candidates, return the best bounds so far
      no_new_cands_break <- TRUE
      if(print_progress) cat("\nFinnished!                                         \n")
      break
    }

    new_candidates <- all_matprod_inds[, which_new_candidates, drop=FALSE] # New candidate products
    all_n_candidates[k] <- length(which_new_candidates) # The number of new candidate products

    # Calculate the new lower bound alpha_k
    all_alpha[k] <- max(all_alpha[k - 1],
                        max(vapply(1:ncol(new_candidates), function(i1) max(abs(eigen(mat_prod(S[, , new_candidates]))$values))^(1/k),
                                   numeric(1))))

    # Calculate the new upper bound beta_k
    all_beta[k] <- min(all_beta[k - 1], max(all_alpha[k - 1] + epsilon, max(all_mu[which_new_candidates])))

    # Update "old stuff" for the next round
    all_matprod_inds_old <- all_matprod_inds

    if(print_progress) {
      cat(paste0("Iteration: ", k, ", current bounds: ", round(all_alpha[k], 4), ", ", round(all_beta[k], 4)), "\r")
    }

    # Stop iteration if the lower and upper bounds are close enough
    if(all_beta[k] - all_alpha[k] <= epsilon) {
      no_new_cands_break <- FALSE
      if(print_progress) cat("\nFinnished!                                         \n")
      break
    }
  }

  # Return lower and upper bounds for the joint spectral radius
  c(max(all_alpha, na.rm=TRUE), min(all_beta, na.rm=TRUE))
}



#' @title Calculate upper bound for the joint spectral radius of the "companion form AR matrices" of the regimes
#'
#' @description \code{bound_JSR} calculates an upper bound for the joint spectral radius of the
#'  "companion form AR matrices" matrices of the regimes as described in Parrilo and Jadbabaie (2008),
#'   Theorems 4.2 and 4.3 and Equations (8) and (11).
#'
#' @inheritParams diagnostic_plot
#' @inheritParams bound_jsr_JP
#' @inherit bound_jsr_JP details references return
#' @examples
#' # p=1, M=2, d=2, relative dens weight function
#' theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863,
#'   -0.3124, 0.505251, 0.298483, 0.030096, -0.176925, 0.838898, 0.310863, 0.007512,
#'   0.018244, 0.949533, -0.016941, 0.121403, 0.573269)
#' mod122 <- STVAR(data=gdpdef, p=1, M=2, params=theta_122relg)
#' # Absolute values of the eigenvalues of the "companion form AR matrices":
#' summary(mod122)$abs_boldA_eigens
#' # It is a necessary (but not sufficient!) condition for ergodic stationary that
#' # the spectral radius of the "companion form AR matrices" are smaller than one
#' # for all of the regimes. A sufficient (but not necessary) condition for
#' # ergodic stationary is that the joint spectral radius of the #companion form
#' # AR matrices" of the regimes is smaller than one. Therefore, we calculate
#' # upper bounds for the joint spectral radius.
#'
#' # Upper bound for the joint spectral radius of the "companion form AR matrices":
#' bound_JSR(mod122, accuracy="0.707")
#' # accuracy "0.707" gives a rough approximation, but the upper bound is already
#' # smaller than one, so the sufficient condition for ergodic stationarity is satisfied.
#'
#' # Higher accuracy gives tighter upper bound:
#' bound_JSR(mod122, accuracy="0.840")
#' bound_JSR(mod122, accuracy="0.917")
#' # bound_JSR(mod122, accuracy="0.957") # Takes a while to compute!
#' # bound_JSR(mod122, accuracy="0.978") # Takes much longer to compute!
#' @export

bound_JSR <- function(stvar, accuracy=c("0.707", "0.840", "0.917", "0.957", "0.978")) {
  accuracy <- match.arg(accuracy)
  all_A <- pick_allA(p=stvar$model$p, M=stvar$model$M, d=stvar$model$d, params=stvar$params)
  all_boldA <- form_boldA(p=stvar$model$p, M=stvar$model$M, d=stvar$model$d, all_A=all_A)
  bound_jsr_JP(all_boldA=all_boldA, accuracy=accuracy)
}


