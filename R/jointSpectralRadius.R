
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
#'  typically the "bold A" matrices, using the algorithm by Gripenberg (1996).
#'
#' @param S the set of matrices the bounds should be calculated for in an array, in VAR applications,
#'  all \eqn{((dp)x(dp))} "bold A" (companion form) matrices in a 3D array, so that \code{[, , m]} gives the matrix
#'  the regime \code{m}.
#' @param epsilon a strictly positive real number that approximately defines the length of the interval between the lower
#'   and upper bounds in Gripenberg's method. A smaller epsilon value results in a narrower interval, thus providing better
#'   accuracy for the bounds, but at the cost of increased computational effort.
#' @param ncores the number of cores to be used in parallel computing in the Gripenberg's algorithm.
#' @param print_progress logical: should the progress of the Gripenberg's algorithm be printed?
#' @details The bounds are calculated using the Gripenberg's (1996) branch-and-bound method, which is also discussed
#'  in Chand and Blondel (2013). Specifically, Kheifets and Saikkonen (2020) show that if the joint spectral radius
#'  of the companion form AR matrices of the regimes is smaller than one, the STVAR process is ergodic stationary. Therefore,
#'  if the upper bound is smaller than one, the process is stationary ergodic. However, as the condition is not
#'  necessary but sufficient and also because the bound might be too conservative, upper bound larger than one
#'  does not imply that the process is not ergodic stationary. You can try higher accuracy, and if the bound is
#'  still larger than one, the result does not tell whether the process is ergodic stationary or not.
#'
#'  Note that with high precision (small \code{epsilon}), the computational effort required are substantial and
#'  the estimation may take very long, even though the function takes use of parallel computing. This is because
#'  with small epsilon the the number of candidate solutions in each iteration may grow exponentially and a large
#'  number of iterations may be required. For this reason, the algorithm starts with a large epsilon, and then
#'  decreases it when new candidate solutions are not found, until the desired epsilon is reached.
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

bound_jsr_G <- function(S, epsilon=0.01, ncores=2, print_progress=TRUE) {
  n <- dim(S)[1] # The dimension of the n x n matrices
  m <- dim(S)[3] # The number of matrices
  maxit <- 1000 # Maximum number of iterations (just some large number)
  all_alpha <- array(NA, dim=maxit) # Storage for the lower bound in each iteration
  all_beta <- array(NA, dim=maxit)  # Storage for the upper bound in each iteration
  stopifnot(length(ncores) == 1 && ncores %% 1 == 0 && ncores > 0)
  if(ncores > parallel::detectCores()) {
    ncores <- parallel::detectCores()
    cat("ncores was set larger than the number of cores available in the system. Using", ncores, "cores.\n")
  }
  epsilon_goal <- epsilon # The epsilon value that we want to achieve
  if(epsilon < 0.2) {
    epsilon <- 0.2 # The epsilon value that we use in the algorithm (can be larger than epsilon_goal)
  }

  ### Step 0: create functions and storages used in the algorithm

  # Function to calculate the norm that we employ as our norm in the algorithm:
  #matrix_norm <- function(A) sqrt(sum(A^2)) # Frobenius norm
  matrix_norm <- function(A) sqrt(max(eigen(crossprod(A), symmetric=TRUE)$values)) # Spectral norm

  # An environment for storing the matrix products calculated in the algorithm similarly to a hash map
  matprod_env <- new.env(hash=TRUE, parent=emptyenv())

  # Store the matrices in S to the hash map:
  for(i1 in 1:m) {
    matprod_env[[paste(i1, collapse="")]] <- S[, , i1]
  }

  # An environment for storing the (Frobenius) norms of the matrix products calculated in the algorithm:
  norm_env <- new.env(hash=TRUE, parent=emptyenv())

  # Store the (Frobenius) norms of the matrices S to the hash map:
  S_norms <- numeric(m) # Storage for the Frobenius norms of S
  for(i1 in 1:m) {
    S_norms[i1] <- matrix_norm(S[, , i1])
    norm_env[[paste(i1, collapse="")]] <- S_norms[i1]
  }

  # A function to calculate the product of a set of matrices in P in an [n, n, k] array
  matprod <- function(P) {
    all_prods <- array(NA, dim=c(n, n, dim(P)[3]))
    all_prods[, , 1] <- P[, , 1]
    for(i1 in 2:dim(P)[3]) {
      all_prods[, , i1] <- all_prods[, , i1-1]%*%P[, , i1]
    }
    all_prods[, , dim(P)[3]]
  }

  # A function to calculate the mu(MP): takes in a vector containing the indices of the matrices in S involved
  # in the matrix product in the order they appear in the matrix product.
  mu <- function(inds, k, matprod_env) { # k>2 assumed
    res <- numeric(k)
    res[1] <- S_norms[inds[1]] # First matrix of the product is a matrix in S.
    for(i1 in 2:k) { # Go through the rest of the products M_1M_2...M_i in the matrix product
      if(!exists(paste(inds[1:i1], collapse=""), envir=matprod_env)) { # Check whether the product is calculated
        # If not and, check whether the product M_2...M_{i} is calculated
        # (M_1 is omitted since due to pre-multiplication, M_2,..,M_{i} should be calculated)
        if(exists(paste(inds[2:i1], collapse=""), envir=matprod_env)) {
          # If yes, calculate the product M_1(M_2...M_{i}) and store it to the hash map:
          matprod_env[[paste(inds[1:i1], collapse="")]] <- S[, , inds[1]]%*%matprod_env[[paste(inds[2:i1], collapse="")]]
        } else {
          # If not (should never happen), calculate the product M_1M_2...M_{i} from scratch and store it to the hash map:
          matprod_env[[paste(inds[1:i1], collapse="")]] <- S[, , inds[1]]%*%matprod(S[, , inds[2:i1]])
        }
      }
      # Calculate the Frobenius norm of the product:
      if(exists(paste(inds[1:i1], collapse=""), envir=norm_env)) { # Check whether the norm is calculated
        # If yes, use the stored norm
        res[i1] <- norm_env[[paste(inds[1:i1], collapse="")]]
      } else {
        # If not, calculate the norm and store it to the hash map:
        res[i1] <- matrix_norm(matprod_env[[paste(inds[1:i1], collapse="")]])^(1/i1)
        norm_env[[paste(inds[1:i1], collapse="")]] <- res[i1]
      }
    }
    min(res) # Return the minimum of the Frobenius norms of the products
  }

  # A function to calculate a product of matrices: takes in a vector containing the indices of
  # the matrices in S involved in the matrix product in the order they appear in the matrix product.
  matprod_hash <- function(inds, matprod_env) {
    if(exists(paste(inds, collapse=""), envir=matprod_env)) { # Check if already calculated
      return(matprod_env[[paste(inds, collapse="")]]) # If so, return it.
    } else {
      return(matprod(S[, , inds])) # If not, calculate from scratch and return it (should never happen)
    }
  }

  ### Step 1, initialize

  # Set up the cluster for parallel computing
  cl <- parallel::makeCluster(ncores)
  on.exit(try(parallel::stopCluster(cl), silent=TRUE)) # Close the cluster on exit, if not already closed.
  parallel::clusterExport(cl, varlist=c("matprod_env", "S", "S_norms", "norm_env", "matrix_norm", "mu", "matprod_hash"),
                          envir=environment())

  # Calculate the initial lower and upper bounds for the joint spectral radius
  all_alpha[1] <- max(vapply(1:m, function(i1) max(abs(eigen(S[, , i1])$values)),
                             numeric(1))) # Max of the spectral radiusses of the matrices in S
  all_beta[1] <- max(vapply(1:m, function(i1) matrix_norm(S[, , i1]), numeric(1))) # Max of the Frobenius norms of the matrices in S

  S_norms <- vapply(1:m, function(i1) matrix_norm(S[, , i1]), numeric(1)) # The Frobenius norms of the matrices in S

  ### Step 2, iteration process

  # First iteration, k=1:
  # The set MP is just the set S and the "products" are the single matrices in S.
  all_mu <- vapply(1:m, function(i1) matrix_norm(S[, , i1]), numeric(1)) # Calculate the mu(S)
  all_n_candidates <- array(NA, dim=maxit)  # The number of candidate matrix products in each iteration
  all_n_candidates[1] <- m # The number of candidate matrix products in the first iteration is the number of matrices in S
  which_new_candidates <- 1:m # The indices of the new candidate matrix products in the first iteration
  all_matprod_inds_old <- matrix(1:2, ncol=2) # Initialize matrix for storing the indices of the matrices involved in the matrix prod

  if(print_progress) {
    cat(paste0("Iteration: ", 1, ", current bounds: ", round(all_alpha[1], 4), ", ", round(all_beta[1], 4)), "\r")
  }

  for(k in 2:maxit) {
    # For each iteration, calculate the sets of matrices that are involved in the concerned matrix products.
    # This set is the set of all possible matrix products obtained by pre-multiplying a matrix product from the previous set of
    # candidate solutions by a matrix from the initial set S. The matrix products involved will consist of k matrices each.

    ## Initialize matrix for storing the indices of the matrices involved in the matrix products:
    # Each column contains a matrix product given by the indices of the matrices in the set S
    n_new_matprods <- m*all_n_candidates[k-1]
    all_matprod_inds <- matrix(nrow=k, ncol=n_new_matprods)  # [matrix_indices, index_of_the_product]
    all_matprod_inds[2:k,] <- all_matprod_inds_old[,which_new_candidates] # Copy the indices from the previous iter: repeats m times
    all_matprod_inds[1,] <- c(vapply(1:m, function(i1) rep(i1, times=all_n_candidates[k-1]),
                                     numeric(all_n_candidates[k-1]))) # Add the indices of the matrices from the initial set S
    # Note that due to pre-multiplication, the first row is for the new matrices in the product.

    ## Calculate the mu(MP) for the set of matrices MP in an [n, n, k] array, for iteration k>2
    all_mu <- numeric(n_new_matprods) # Initialize storage for the mu(MP) values
    all_mu <- parallel::parSapply(cl=cl, X=1:n_new_matprods,
                                  FUN=function(i1) mu(inds=all_matprod_inds[,i1], k=k, matprod_env=matprod_env))

    ## Update the matprod_env to the cluster for calculation of alphas and for the next iteration.
    parallel::clusterExport(cl=cl, varlist="matprod_env", envir=environment())

    ## Construct the new set of candidates MP, i.e., the products for which mu(MP) > alpha_{k-1} + epsilon
    which_new_candidates <- which(all_mu > all_alpha[k-1] + epsilon)
    if(length(which_new_candidates) == 0) {
      if(epsilon == epsilon_goal) {
        # If there are no new candidates with epsilon_goal, return the best bounds so far
        if(print_progress) cat("\nFinnished!                                         \n")
        break
      } else {
        # If no new candidates with the larger epsilon, switch to smaller epsilon
        while(TRUE) { # Decrease epsilon until epsilon_goal is obtained or new candidates are found
          if(epsilon/2 > epsilon_goal) {
            epsilon <- epsilon/2 # Decrease epsilon
          } else {
            epsilon <- epsilon_goal
          }
          which_new_candidates <- which(all_mu > all_alpha[k-1] + epsilon)
          if(length(which_new_candidates) > 0) {
            break
          } else if(length(which_new_candidates == 0) && epsilon == epsilon_goal) {
            # If there are no new candidates with epsilon_goal, return the best bounds so far
            if(print_progress) cat("\nFinnished!                                         \n")
            break
          }
        }
      }
    }

    new_candidates <- all_matprod_inds[, which_new_candidates, drop=FALSE] # New candidate products
    all_n_candidates[k] <- length(which_new_candidates) # The number of new candidate products

    ## Calculate the new lower bound alpha_k
    all_alpha[k] <- max(all_alpha[k - 1],
                        max(parallel::parSapply(cl=cl, X=1:ncol(new_candidates),
                                                FUN=function(i1) max(abs(eigen(matprod_hash(inds=new_candidates[,i1],
                                                                                            matprod_env=matprod_env))$values))^(1/k))))

    ## Calculate the new upper bound beta_k
    all_beta[k] <- min(all_beta[k - 1], max(all_alpha[k - 1] + epsilon, max(all_mu[which_new_candidates])))

    ## Update "old stuff" for the next round
    all_matprod_inds_old <- all_matprod_inds

    if(print_progress) {
      cat(paste0("Iteration: ", k, ", current bounds: ", round(all_alpha[k], 4), ", ", round(all_beta[k], 4)), "\r")
    }

    ## Stop iteration if the lower and upper bounds are close enough
    if(all_beta[k] - all_alpha[k] <= epsilon) {
      if(print_progress) cat("\nFinnished!                                         \n")
      break
    }
    if(k == maxit) {
      cat("\nThe maximum number of iterations reached!                                         \n")
    }
  }

  # Return lower and upper bounds for the joint spectral radius
  parallel::stopCluster(cl=cl) # Stop the cluster
  c(max(all_alpha, na.rm=TRUE), min(all_beta, na.rm=TRUE))
}


#' @title Calculate upper bound for the joint spectral radius of the "companion form AR matrices" of the regimes
#'
#' @description \code{bound_JSR} calculates an bounds for the joint spectral radius of the
#'  "companion form AR matrices" matrices of the regimes to assess the validity of the stationarity condition.
#'
#' @inheritParams diagnostic_plot
#' @inheritParams bound_jsr_G
#' @param method which method should be used for bounding the joint spectral radius? Gripenberg's method is highly recommended
#'   over JP (see the details section).
#' @param JP_accuracy what should the relative accuracy of the upper bound of Jadabaie and Parrilo be (for \code{method="JP"})?
#' @details A sufficient condition for ergodic stationarity of the STVAR processes implemented in \code{sstvars} is that the joint
#'  spectral radius of the "companion form AR matrices" of the regimes is smaller than one (Kheifets and Saikkonen, 2020).
#'  This function calculates an upper (and lower) bound for the JSR and is implemented to assess the validity of this condition
#'  in practice. If the bound is smaller than one, the model is deemed ergodic stationary.
#'
#'  Currently, two methods are implemented: the branch-and-bound method by Gripenberg (1996) and the upper bound by
#'  Jadbabaie and Parrilo (2008). We highly recommend using the Gripenberg's method, as it is mainly much faster and less prone to
#'  memory issues than the upper bound by Jadbabaie and Parrilo (2008). Calculation of the latter with good enough accuracy may not be
#'  feasible for other than very small models. However, for large models also Gripenberg's method may take very long if tight bounds
#'  are required. When \code{print_progress == TRUE}, the tightest bounds found so-far are printed in each iteration of Gripenberg's
#'  algorithm, so you can also just terminate the algorithm when the bounds are tight enough for your purposes. Consider also
#'  adjusting the argument \code{epsilon}, as larger epsilon does not just make the bounds less tight but also speeds up the algorithm
#'  significantly.
#'
#'  Various methods for bounding the JSR are discussed and compared in Chang and Blondel (2013).
#' @return Returns lower and upper bounds for the joint spectral radius of the "companion form AR matrices" of the regimes.
#'   If the upper bound by Jadbabaie and Parrilo (2008) is calculated, only the upper bound is returned.
#' @references
#'  \itemize{
#'  \item Chang C-T, Blondel, V.D. 2013. An experimental study of approximation algorithms for the joint spectral radius.
#'      \emph{Numerical algorithms}, \strong{64}, 181-202.
#'    \item Gripenberg, G. 1996. Computing the joint spectral radius. \emph{Linear Algebra and its Applications},
#'      234, 43–60.
#'    \item Kheifets I.L., Saikkonen P.J. 2020. Stationarity and ergodicity of Vector STAR models.
#'      \emph{Econometric Reviews}, \strong{39}:4, 407-414.
#'    \item Parrilo P.A., Jadbabaie A. 2008. Approximation of the joint spectral
#'       radius using sum of squares. \emph{Linear Algebra and its Applications}, \strong{428}, 2385-2402.
#'  }
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
#' # ergodic stationary is that the joint spectral radius of the companion form
#' # AR matrices" of the regimes is smaller than one. Therefore, we calculate
#' # an upper (and lower) bound for the joint spectral radius.
#'
#' ## Bounds by Gripenberg's method (default and recommonded).
#' # Since the largest modulus of the companion form AR matrices is not very close
#' # to one, we likely won't need very thight bounds to verify the JSR is smaller
#' # than one. Thus, we set epsilon=0.01 so that the interval between the lower
#' # and upper bound is roughly 0.01:
#' bound_JSR(mod122, epsilon=0.01, method="Gripenberg")
#' # The upper bound is smaller than one, so the model is ergodic stationary.
#'
#' # If we want tighter bounds, we can set smaller epsilon, e.g., epsilon=0.001:
#' bound_JSR(mod122, epsilon=0.001, method="Gripenberg")
#'
#' ## Upper bound by the method of Jadbabaie and Parrilo (2008):
#' # We start with the lowest relative accuracy to see if it is already enough
#' # to verify that the JSR is smaller than one:
#' bound_JSR(mod122, method="JP", JP_accuracy="0.707")
#' # accuracy "0.707" gives a rough approximation, but the upper bound is already
#' # smaller than one, so the sufficient condition for ergodic stationarity is
#' # satisfied.
#'
#' # Higher accuracy gives tighter upper bound:
#' bound_JSR(mod122, method="JP", JP_accuracy="0.840")
#' bound_JSR(mod122, method="JP", JP_accuracy="0.917")
#' # bound_JSR(mod122, method="JP", JP_accuracy="0.957") # Takes a while to compute!
#' @export

bound_JSR <- function(stvar, epsilon=0.01, ncores=2, print_progress=TRUE, method=c("Gripenberg", "JP"),
                      JP_accuracy=c("0.707", "0.840", "0.917", "0.957", "0.978")) {
  check_stvar(stvar)
  stopifnot(is.numeric(epsilon) && epsilon > 0)
  JP_accuracy <- match.arg(JP_accuracy)
  method <- match.arg(method)
  all_A <- pick_allA(p=stvar$model$p, M=stvar$model$M, d=stvar$model$d, params=stvar$params)
  all_boldA <- form_boldA(p=stvar$model$p, M=stvar$model$M, d=stvar$model$d, all_A=all_A)
  if(method == "Gripenberg") {
    return(bound_jsr_G(S=all_boldA, epsilon=epsilon, ncores=ncores, print_progress=print_progress))
  } else { # method == "JP"
    return(bound_jsr_JP(all_boldA=all_boldA, accuracy=JP_accuracy))
  }
}


