#' @title Calculate upper bound for the joint spectral radius of a set of matrices
#'
#' @description \code{bound_jsr_G} calculates lower and upper bounds for the joint spectral radious of a set of square matrices,
#'  typically the "bold A" matrices, using the algorithm by Gripenberg (1996).
#'
#' @param S the set of matrices the bounds should be calculated for in an array, in STVAR applications,
#'  all \eqn{((dp)x(dp))} "bold A" (companion form) matrices in a 3D array, so that \code{[, , m]} gives the matrix
#'  the regime \code{m}.
#' @param epsilon a strictly positive real number that approximately defines the goal of length of the interval between the lower
#'   and upper bounds. A smaller epsilon value results in a narrower interval, thus providing better accuracy for the bounds,
#'   but at the cost of increased computational effort. Note that the bounds are always wider than \code{epsilon} and it is not
#'   obvious what \code{epsilon} should be chosen obtain bounds of specific tightness.
#' @param adaptive_eps logical: if \code{TRUE}, starts with a large epsilon and then decreases it gradually whenever the progress
#'   of the algorithm requires, until the value given in the argument \code{epsilon} is reached. Usually speeds up the algorithm
#'   substantially but is an unconventional approach, and there is no guarantee that the algorithm converges appropriately towards
#'   bounds with the tightness given by the argument \code{epsilon}.
#' @param ncores the number of cores to be used in parallel computing.
#' @param print_progress logical: should the progress of the algorithm be printed?
#' @details The upper and lower bounds are calculated using the Gripenberg's (1996) branch-and-bound method, which is also discussed
#'  in Chang and Blondel (2013). This function can be generally used for approximating the JSR of a set of square matrices, but the
#'  main intention is STVAR applications (for models created with \code{sstvars}, the function \code{bound_JSR} should be preferred).
#'  Specifically, Kheifets and Saikkonen (2020) show that if the joint spectral radius of the companion form AR matrices of the regimes
#'  is smaller than one, the STVAR process is ergodic stationary. Virolainen (2024) shows the same result for his parametrization of
#'  of threshold and smooth transition vector autoregressive models. Therefore, if the upper bound is smaller than one, the process is
#'  stationary ergodic. However, as the condition is not necessary but sufficient and also because the bound might be too conservative,
#'  upper bound larger than one does not imply that the process is not ergodic stationary. You can try higher accuracy, and if the bound
#'  is still larger than one, the result does not tell whether the process is ergodic stationary or not.
#'
#'  Note that with high precision (small \code{epsilon}), the computational effort required are substantial and
#'  the estimation may take long, even though the function takes use of parallel computing. This is because
#'  with small epsilon the the number of candidate solutions in each iteration may grow exponentially and a large
#'  number of iterations may be required. For this reason, \code{adaptive_eps=TRUE} can be considered for large matrices,
#'  in which case the algorithm starts with a large epsilon, and then decreases it when new candidate solutions are
#'  not found, until the epsilon given by the argument \code{epsilon} is reached.
#' @return Returns an upper bound for the joint spectral radius of the "companion form AR matrices" of the regimes.
#' @seealso \code{\link{bound_JSR}}
#' @references
#'  \itemize{
#'  \item C-T Chang and V.D. Blondel. 2013 . An experimental study of approximation algorithms for the joint spectral radius.
#'      \emph{Numerical algorithms}, \strong{64}, 181-202.
#'    \item Gripenberg, G. 1996. Computing the joint spectral radius. \emph{Linear Algebra and its Applications},
#'      234, 43–60.
#'    \item I.L. Kheifets, P.J. Saikkonen. 2020. Stationarity and ergodicity of Vector STAR models.
#'      \emph{Econometric Reviews}, \strong{39}:4, 407-414.
#'    \item Virolainen S. 2024. Identification by non-Gaussianity in structural threshold and
#'       smooth transition vector autoregressive models. Unpublished working
#'       paper, available in ArXiv.
#'  }
#' @examples
#' \donttest{
#' ## Below examples take approximately 5 seconds to run.
#'
#' # A set of two (5x5) square matrices:
#' set.seed(1); S1 <- array(rnorm(20*20*2), dim=c(5, 5, 2))
#'
#' # Bound the joint spectral radius of the set of matrices S1, with the
#' # approximate tightness epsilon=0.01:
#' bound_jsr_G(S1, epsilon=0.01, adaptive_eps=FALSE)
#'
#' # Obtain bounds faster with adaptive_eps=TRUE:
#' bound_jsr_G(S1, epsilon=0.01, adaptive_eps=TRUE)
#' # Note that the upper bound is not the same as with adaptive_eps=FALSE.
#'
#' # A set of three (3x3) square matrices:
#' set.seed(2); S2 <- array(rnorm(3*3*3), dim=c(3, 3, 3))
#'
#' # Bound the joint spectral radius of the set of matrices S2:
#' bound_jsr_G(S2, epsilon=0.01, adaptive_eps=FALSE)
#'
#' # Larger epsilon terminates the iteration earlier and results in wider bounds:
#' bound_jsr_G(S2, epsilon=0.05, adaptive_eps=FALSE)
#'
#' # A set of eight (2x2) square matrices:
#' set.seed(3); S3 <- array(rnorm(2*2*8), dim=c(2, 2, 8))
#'
#' # Bound the joint spectral radius of the set of matrices S3:
#' bound_jsr_G(S3, epsilon=0.01, adaptive_eps=FALSE)
#' }
#' @export

bound_jsr_G <- function(S, epsilon=0.01, adaptive_eps=FALSE, ncores=2, print_progress=TRUE) {
  n <- dim(S)[1] # The dimension of the n x n matrices
  m <- dim(S)[3] # The number of matrices
  maxit <- 1000 # Maximum number of iterations (just some large number)
  all_alpha <- array(NA, dim=maxit) # Storage for the lower bound in each iteration
  all_beta <- array(NA, dim=maxit)  # Storage for the upper bound in each iteration
  stopifnot(length(ncores) == 1 && ncores %% 1 == 0 && ncores > 0)
  if(ncores > parallel::detectCores()) {
    ncores <- parallel::detectCores()
    message("ncores was set larger than the number of cores available in the system. Using", ncores, "cores.\n")
  }
  epsilon_goal <- epsilon # The epsilon value that we want to achieve
  if(adaptive_eps && epsilon < 0.2) {
    epsilon <- 0.2 # The initial epsilon value that we use in the algorithm (can be larger than epsilon_goal)
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
  all_matprod_inds_old <- matrix(1:m, ncol=m) # Initialize matrix for storing the indices of the matrices involved in the matrix prod

  if(print_progress) {
    message(paste0("Iteration: ", 1, ", current bounds: ", round(all_alpha[1], 4), ", ", round(all_beta[1], 4)), "\r")
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
      if(isTRUE(all.equal(epsilon, epsilon_goal))) {
        # If there are no new candidates with epsilon_goal, return the best bounds so far
        if(print_progress) {
          if(k %% 10 != 0) {
            message(paste0("Iteration: ", k, ", current bounds: ", round(max(all_alpha, na.rm=TRUE), 4), ", ",
                           round(min(all_beta, na.rm=TRUE), 4)), "\r")
          }
          message("\nFinnished!                                         ")
        }
        break
      } else {
        # If no new candidates with the larger epsilon, switch to smaller epsilon
        # NOTE: The bounds obtained with the initial epsilon are valid by Gripenberg's (1996) Theorem 3.
        # When the epsilon decreases, we kind of restart the algorithm with the difference that instead of
        # using the matrices in S as the initial set of candidate "products", we use the candidate products found
        # with the larger epsilon. In the proof of Theorem 3, Gripenberg (1996) shows that the upper bound
        # is valid for arbitrary set of candidate products comprised of the matrices in S. Smaller epsilon discards
        # less candidate products, and candidate products that were not discarded by a larger epsilon would anyway be
        # in the set of candidate products when running the algorithm with the smaller epsilon. Taking max
        # (or sup) over the mu-values obtained with the decreased set of candidate products that is a subset of the larger set
        # cannot be larger than the max over the larger set. Therefore, the bounds obtained with the decreased epsilon should
        # be valid. However, the asymptotic converge to bounds of desired tightness is not guaranteed, since some of the candidate
        # products that potentially contribute to tighter bounds may have been discarded with the larger epsilon.
        # ALSO NOTE: Sometimes the algorithm takes longer and converges closer to the desired tightness because the lower bound
        # increases slower, leading to more iterations that tighten the upper bound.
        while(TRUE) { # Decrease epsilon until epsilon_goal is obtained or new candidates are found
          if(epsilon/2 > epsilon_goal) {
            epsilon <- epsilon/2 # Decrease epsilon until candidates found or epsilon goal is met
          } else {
            epsilon <- epsilon_goal
          }
          which_new_candidates <- which(all_mu > all_alpha[k-1] + epsilon)
          if(length(which_new_candidates) > 0) {
            break # Breaks the inner while loop if new candidates are found
          } else if(length(which_new_candidates) == 0 && isTRUE(all.equal(epsilon, epsilon_goal))) {
            break # Breaks the inner while loop if no new candidates are found and epsilon_goal is met
          }
        }
        if(length(which_new_candidates) == 0) {
          # If there are no new candidates with epsilon_goal, return the best bounds so far
          if(print_progress) {
            if(k %% 10 != 0) {
              message(paste0("Iteration: ", k, ", current bounds: ", round(max(all_alpha, na.rm=TRUE), 4), ", ",
                             round(min(all_beta, na.rm=TRUE), 4)), "\r")
            }
            message("\nFinnished!                                         ")
          }
          break
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

    if(print_progress && k %% 10 == 0) {
      message(paste0("Iteration: ", k, ", current bounds: ", round(all_alpha[k], 4), ", ", round(all_beta[k], 4)), "\r")
    }

    ## Stop iteration if the lower and upper bounds are close enough
    if(all_beta[k] - all_alpha[k] <= epsilon) {
      if(print_progress) {
        if(k %% 10 != 0) {
          message(paste0("Iteration: ", k, ", current bounds: ", round(max(all_alpha, na.rm=TRUE), 4), ", ",
                         round(min(all_beta, na.rm=TRUE), 4)), "\r")
        }
        message("\nFinnished!                                         ")
      }
      break
    }
    if(k == maxit) {
      if(print_progress) {
        if(k %% 10 != 0) {
          message(paste0("Iteration: ", k, ", current bounds: ", round(max(all_alpha, na.rm=TRUE), 4), ", ",
                         round(min(all_beta, na.rm=TRUE), 4)), "\r")
        }
      }
      message("\nThe maximum number of iterations reached!                                         ")
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
#' @details A sufficient condition for ergodic stationarity of the STVAR processes implemented in \code{sstvars} is that the joint
#'  spectral radius of the "companion form AR matrices" of the regimes is smaller than one (Kheifets and Saikkonen, 2020).
#'  This function calculates an upper (and lower) bound for the JSR and is implemented to assess the validity of this condition
#'  in practice. If the bound is smaller than one, the model is deemed ergodic stationary.
#'
#'  Implements the branch-and-bound method by Gripenberg (1996) in the conventional form (\code{adaptive_eps=FALSE}) and in a form
#'  incorporating "adaptive tightness" (\code{adaptive_eps=FALSE}). The latter approach is unconventional and does not guarantee
#'  appropriate convergence of the bounds close to the desired tightness given in the argument \code{epsilon}, but it usually
#'  substantially speeds up the algorithm. When \code{print_progress==TRUE}, the tightest bounds found so-far are printed in each
#'  iteration of the algorithm, so you can also just terminate the algorithm when the bounds are tight enough for your purposes.
#'  Consider also adjusting the argument \code{epsilon}, in particular when \code{adaptive_eps=FALSE}, as larger epsilon does not
#'  just make the bounds less tight but also speeds up the algorithm significantly. See Chang and Blondel (2013) for a discussion
#'  on variuous methods for bounding the JSR.
#' @return Returns lower and upper bounds for the joint spectral radius of the "companion form AR matrices" of the regimes.
#' @seealso \code{\link{bound_jsr_G}}
#' @inherit bound_jsr_G references
#' @examples
#' \donttest{
#' ## Below examples take approximately 5 seconds to run.
#'
#' # Gaussian STVAR p=1, M=2 model with weighted relative stationary densities
#' # of the regimes as the transition weight function:
#' theta_122relg <- c(0.734054, 0.225598, 0.705744, 0.187897, 0.259626, -0.000863,
#'   -0.3124, 0.505251, 0.298483, 0.030096, -0.176925, 0.838898, 0.310863, 0.007512,
#'   0.018244, 0.949533, -0.016941, 0.121403, 0.573269)
#' mod122 <- STVAR(data=gdpdef, p=1, M=2, params=theta_122relg)
#'
#' # Absolute values of the eigenvalues of the "companion form AR matrices":
#' summary(mod122)$abs_boldA_eigens
#' # It is a necessary (but not sufficient!) condition for ergodic stationary that
#' # the spectral radius of the "companion form AR matrices" are smaller than one
#' # for all of the regimes. A sufficient (but not necessary) condition for
#' # ergodic stationary is that the joint spectral radius of the companion form
#' # AR matrices" of the regimes is smaller than one. Therefore, we calculate
#' # bounds for the joint spectral radius.
#'
#' ## Bounds by Gripenberg's (1996) branch-and-bound method:
#' # Since the largest modulus of the companion form AR matrices is not very close
#' # to one, we likely won't need very thight bounds to verify the JSR is smaller
#' # than one. Thus, using a small epsilon would make the algorithm unnecessarily slow,
#' # so we use the (still quite small) epsilon=0.01:
#' bound_JSR(mod122, epsilon=0.01, adaptive_eps=FALSE)
#' # The upper bound is smaller than one, so the model is ergodic stationary.
#'
#' # If we want tighter bounds, we can set smaller epsilon, e.g., epsilon=0.001:
#' bound_JSR(mod122, epsilon=0.001, adaptive_eps=FALSE)
#'
#' # Using adaptive_eps=TRUE usually speeds up the algorithm when the model
#' # is large, but with the small model here, the speed-difference is small:
#' bound_JSR(mod122, epsilon=0.001, adaptive_eps=TRUE)
#' }
#' @export

bound_JSR <- function(stvar, epsilon=0.01, adaptive_eps=FALSE, ncores=2, print_progress=TRUE) {
  check_stvar(stvar)
  stopifnot(is.numeric(epsilon) && epsilon > 0)
  all_A <- pick_allA(p=stvar$model$p, M=stvar$model$M, d=stvar$model$d, params=stvar$params)
  all_boldA <- form_boldA(p=stvar$model$p, M=stvar$model$M, d=stvar$model$d, all_A=all_A)
  bound_jsr_G(S=all_boldA, epsilon=epsilon, ncores=ncores, print_progress=print_progress)
}


