#' @title Calculate standard errors for estimates of a smooth transition VAR model
#'
#' @description \code{standard_errors} calculates approximate standard errors for the smooth transition
#'   VAR model using square roots of the diagonal of inverse of observed information matrix
#'   and central-difference approximation for the differentiation.
#'
#' @inheritParams loglikelihood
#' @details This function assumes the standard asymptotic distribution of the estimator
#' @return A vector containing the approximate standard errors of the estimates.
#' @inherit in_paramspace references
#' @keywords internal

standard_errors <- function(data, p, M, params, weight_function=c("relative_dens", "logit"), weightfun_pars=NULL,
                            cond_dist=c("Gaussian", "Student"), parametrization=c("intercept", "mean"),
                            identification=c("reduced_form", "recursive", "heteroskedasticity"),
                            AR_constraints=NULL, mean_constraints=NULL, B_constraints=NULL, minval) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  identification <- match.arg(identification)
  d <- ncol(data)
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars)
  if(missing(minval)) {
    minval <- get_minval(data)
  }

  # The log-likelihood function to differentiate
  loglik_fn <- function(params) {
    tryCatch(loglikelihood(data=data, p=p, M=M, params=params,
                           weight_function=weight_function, weightfun_pars=weightfun_pars,
                           cond_dist=cond_dist, parametrization=parametrization,
                           identification=identification, AR_constraints=AR_constraints,
                           mean_constraints=mean_constraints, B_constraints=B_constraints,
                           check_params=TRUE, to_return="loglik", minval=minval),
             error=function(e) NA)
  }

  # Calculate Hessian
  Hess <- calc_hessian(x=params, fn=loglik_fn, h=6e-6)

  # Inverse of the observed information matrix
  inv_obs_inf <- tryCatch(solve(-Hess), error=function(e) matrix(NA, nrow=length(params), ncol=length(params)))

  # Calculate the standard errors
  unlist(lapply(diag(inv_obs_inf), function(x) ifelse(is.na(x) | x < 0, NA, sqrt(x))))
}



#' @title Print standard errors of a smooth transition VAR model in the same form as the model estimates are printed
#'
#' @description \code{print_std_errors} prints the approximate standard errors of a smooth transition VAR model in the
#'   same form as the parameters of objects of class \code{'stvar'} are printed.
#'
#' @inheritParams get_boldA_eigens
#' @inheritParams print.stvar
#' @details \strong{Note that the approximate standard errors are based on the unverified assumption of
#'   asymptotic normality of the estimator!}
#'
#'   The main purpose of \code{print_std_errors} is to provide a convenient tool to match the standard
#'   errors to certain parameter estimates. Note that if the model is intercept parametrized, there won't
#'   be standard errors for the unconditional means, and vice versa.
#'
#'   Note that if linear constraints are imposed and they involve summations or multiplications, then the AR
#'   parameter standard errors are printed separately as they don't correspond one-to-one to the model parameter
#'   standard errors.
#' @seealso \code{\link{fitSTVAR}}, \code{\link{STVAR}}, \code{\link{print.stvar}},
#'  \code{\link{swap_parametrization}}
#' @inherit STVAR references
#' @examples
#' \donttest{
#' # STVAR p=1, M=2 models
#' mod12 <- fitSTVAR(gdpdef, p=1, M=2, nrounds=1, seeds=4)
#' mod12
#' print_std_errors(mod12)
#' }
#' @export

print_std_errors <- function(stvar, digits=3) {
  if(!all_pos_ints(digits)) stop("Argument digits must be positive integer")
  format_value <- format_valuef(digits)
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  var_names <- colnames(stvar$data)
  if(is.null(var_names)) var_names <- paste0("Var.", 1:d)
  weight_function <- stvar$model$weight_function
  cond_dist <- stvar$model$cond_dist
  parametrization <- stvar$model$parametrization
  identification <- stvar$model$identification
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  B_constraints <- stvar$model$B_constraints
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function, weightfun_pars=stvar$model$weightfun_pars)
  npars <- length(stvar$params)
  pars <- stvar$std_errors
  pars <- reform_constrained_pars(p=p, M=M, d=d, params=pars, weight_function=weight_function,
                                  weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                                  identification=identification, AR_constraints=AR_constraints,
                                  mean_constraints=mean_constraints, B_constraints=B_constraints)
  all_phi0_or_mu <- pick_phi0(M=M, d=d, params=pars)
  all_A <- pick_allA(p=p, M=M, d=d, params=pars)
  if(identification == "reduced_form") {
    all_Omega <- pick_Omegas(p=p, M=M, d=d, params=pars)
  } else {
    # No standard errors for cov. mats. as the model is parametrized with W and lambdas
    all_Omega <- array(" ", dim=c(d, d, M))
    stop("Structural models are not yet implemented to print_std_errors")
  }
  weightpars <- pick_weightpars(p=p, M=M, d=d, params=pars, weight_function=weight_function, weightfun_pars=weightfun_pars)
  if(cond_dist != "Gaussian") stop("Only Gaussian models are implemented to print_std_errors")
  if(weight_function == "relative_dens") {
    weightpars[M] <- NA # No standard error for the last alpha
  } else if(weight_function == "logit") {
    all_gamma_m <- matrix(weightpars, ncol=M-1) # Column per gamma_m, m=1,...,M-1, gamma_M=0.
  } else {
    stop("only relative dens and logit weigghtfunctions are implemented to print_std_errors")
  }

  if(parametrization == "mean") {
    all_mu <- all_phi0_or_mu
    all_phi0 <- matrix(" ", nrow=d, ncol=M)
  } else {
    all_mu <- matrix(NA, nrow=d, ncol=M)
    all_phi0 <- all_phi0_or_mu
  }
  if(!is.null(AR_constraints)) {
    # The constrained AR parameter standard errors multiplied open in 'pars' are valid iff
    # the constraint matrix contains zeros and ones only, and there is at most one one in
    # each row (no multiplications or summations).
    if(any(AR_constraints != 1 & AR_constraints != 0) | any(rowSums(AR_constraints) > 1)) {
      sep_AR <- TRUE # The AR parameter std errors must be printed separately
      all_A <- array(" ", dim=c(d, d, p, M))
      AR_stds <- stvar$std_errors[(M*d + 1):(M*d + ncol(AR_constraints))] # Constrained AR param std errors
    } else {
      sep_AR <- FALSE
    }
  } else {
    sep_AR <- FALSE # No constraints imposed
  }

  cat(weight_function, cond_dist, "STVAR model,", paste0(identification, ","),
      ifelse(is.null(AR_constraints), "no AR_constraints,", "AR_constraints used,"),
      ifelse(is.null(mean_constraints), paste0("no mean_constraints", ifelse(is.null(B_constraints), "", ",")),
             paste0("mean_constraints used", ifelse(is.null(B_constraints), "", ","))),
      ifelse(identification == "reduced_form", "", ifelse(is.null(B_constraints), "no B_constraints", "B_constraints used")))
  cat("\n", paste0(" p = ", p, ", "))
  cat(paste0("M = ", M, ", "))
  cat(paste0("d = ", d, ", #parameters = " , npars, ","))
  if(weight_function == "logit") {
    cat("\n ", paste0("Switching variables: ", paste0(var_names[weightfun_pars[[1]]], collapse=", "), " with ",
                      weightfun_pars[[2]], ifelse(weightfun_pars[[2]] == 1, " lag.", " lags.")))
  }
  cat("\n\n")
  cat("APPROXIMATE STANDARD ERRORS\n\n")

  left_brackets <- rep("[", times=d)
  right_brackets <- rep("]", times=d)
  plus <- c("+", rep(" ", d - 1))
  round_lbrackets <- rep("(", times=d)
  round_rbrackets <- rep(")", times=d)
  Y <- paste0("Y", 1:d)
  tmp_names <- paste0("tmp", 1:(p*(d + 2) + d + 2))

  # PRINT DISTRIBUTION PARAM SOMEWHERE FOR STUDENT COND DIST

  for(m in seq_len(M)) {
    count <- 1
    cat(paste("Regime", m))
    cat("\n")
    if(weight_function == "relative_dens") {
      if(m < M) cat(paste("Weight param:", format_value(weightpars[m])), "\n")
    } else if(weight_function == "logit") {
      if(m < M) cat(paste("Weight params:", paste0(format_value(all_gamma_m[,m]), collapse=", ")), "\n")
    }
    if(parametrization == "mean") cat("Regime means:", paste0(format_value(all_mu[,m]), collapse=", "), "\n")
    cat("\n")
    df <- data.frame(Y=Y,
                     eq=c("=", rep(" ", d - 1)),
                     eq=left_brackets,
                     phi0=format_value(all_phi0[, m, drop=FALSE]),
                     eq=right_brackets,
                     plus)
    for(i1 in seq_len(p)) {
      Amp_colnames <- c(paste0("A", i1), tmp_names[count:(count + d - 1 - 1)]); count <- count + d - 1
      df[, tmp_names[count]] <- left_brackets; count <- count + 1
      df[, Amp_colnames] <- format_value(all_A[, ,i1 , m])
      df[, tmp_names[count]] <- right_brackets; count <- count + 1
      df[, tmp_names[count]] <- paste0(Y, ".", i1); count <- count + 1
      df <- cbind(df, plus)
    }
    df[, tmp_names[p*(d + 2) + 1]] <- left_brackets
    df[, c("Omega", tmp_names[(p*(d + 2) + 2):(p*(d + 2) + d)])] <- format_value(all_Omega[, , m])
    df[, tmp_names[p*(d + 2) + d + 1]] <- right_brackets
    df[, "1/2"] <- rep(" ", d)
    df[, tmp_names[p*(d + 2) + d + 2]] <- paste0("eps", 1:d)
    names_to_omit <- unlist(lapply(c("plus", "eq", "round_lbrackets", "round_rbrackets", tmp_names),
                                   function(nam) grep(nam, colnames(df))))
    colnames(df)[names_to_omit] <- " "
    print(df)
    cat("\n")
  }
  if(sep_AR) cat(paste0("AR parameters: ", paste0(format_value(AR_stds), collapse=", ")), "\n\n")

  if(identification != "reduced_form") {
    stop("Structural models are not yet implemented to print_std_errors")
    cat("Structural parameters:\n")
    W <- format_value(pick_W(p=p, M=M, d=d, params=pars, structural_pars=structural_pars))

    if(M > 1) {
      lambdas <- format_value(pick_lambdas(p=p, M=M, d=d, params=pars, structural_pars=structural_pars))
      lambdas <- matrix(lambdas, nrow=d, ncol=M - 1, byrow=FALSE) # Column for each regime
    }

    tmp <- c(rep(" ", times=d - 1), ",")
    df2 <- data.frame(left_brackets, W=W[,1])
    for(i1 in 2:d) {
      df2 <- cbind(df2, W[, i1])
      colnames(df2)[1 + i1] <- "tmp"
    }
    df2 <- cbind(df2, right_brackets)
    if(M > 1) {
      tmp <- c(rep(" ", times=d - 1), ",")
      for(i1 in 1:(M - 1)) {
        if(sep_lambda) {
          lmb <- rep(NA, times=d)
        } else {
          lmb <- lambdas[,i1]
        }
        df2 <- cbind(df2, tmp, left_brackets, lmb, right_brackets)
        colnames(df2)[grep("lmb", colnames(df2))] <- paste0("lamb", i1 + 1)
      }
    }
    names_to_omit <- unlist(lapply(c("left_brackets", "right_brackets", "tmp"), function(nam) grep(nam, colnames(df2))))
    colnames(df2)[names_to_omit] <- " "
    print(df2)
    cat("\n")
    W_orig <- gsmvar$model$structural_pars$W
    n_zero <- sum(W_orig == 0, na.rm=TRUE)
    n_free <- sum(is.na(W_orig))
    n_sign <- d^2 - n_zero - n_free
    if(sep_lambda) cat(paste0("lambda parameters: ", paste0(format_value(lambda_stds), collapse=", ")), "\n\n")
    cat("The B-matrix is subject to", n_zero, "zero constraints and", n_sign, "sign constraints.\n")
    cat("\n")
  }
  invisible(stvar)
}
