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

standard_errors <- function(data, p, M, params, weight_function=c("relative_dens", "logistic", "mlogit", "exponential", "threshold"),
                            weightfun_pars=NULL, cond_dist=c("Gaussian", "Student"), parametrization=c("intercept", "mean"),
                            identification=c("reduced_form", "recursive", "heteroskedasticity"),
                            AR_constraints=NULL, mean_constraints=NULL, weight_constraints=NULL, B_constraints=NULL,
                            minval) {
  weight_function <- match.arg(weight_function)
  cond_dist <- match.arg(cond_dist)
  parametrization <- match.arg(parametrization)
  identification <- match.arg(identification)
  d <- ncol(data)
  check_pMd(p=p, M=M, d=d, weight_function=weight_function, identification=identification)
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                         cond_dist=cond_dist)
  check_constraints(p=p, M=M, d=d, weight_function=weight_function, weightfun_pars=weightfun_pars,
                    parametrization=parametrization, identification=identification,
                    AR_constraints=AR_constraints, mean_constraints=mean_constraints,
                    weight_constraints=weight_constraints, B_constraints=B_constraints)
  if(missing(minval)) {
    minval <- get_minval(data)
  }

  # The log-likelihood function to differentiate
  loglik_fn <- function(params) {
    tryCatch(loglikelihood(data=data, p=p, M=M, params=params,
                           weight_function=weight_function, weightfun_pars=weightfun_pars,
                           cond_dist=cond_dist, parametrization=parametrization,
                           identification=identification, AR_constraints=AR_constraints,
                           mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                           B_constraints=B_constraints, check_params=TRUE, to_return="loglik",
                           minval=minval),
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
#' ## These are long-running examples that take approximately 30 seconds to run.
#'
#' # Gaussian STVAR p=1, M=2, model with weighted relative stationary densities
#' # of the regimes as the transition weight function:
#' mod12 <- fitSTVAR(gdpdef, p=1, M=2, nrounds=1, seeds=4) # Estimate the model
#' mod12 # Print the estimates
#' print_std_errors(mod12) # Print approximate standard errors of the estimates
#' # Note that the standard errors are based on the assumption of the standard
#' # asymptotic limiting Gaussian distribution of the estimator.
#' }
#' @export

print_std_errors <- function(stvar, digits=3) {
  check_stvar(stvar)
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
  weight_constraints <- stvar$model$weight_constraints
  B_constraints <- stvar$model$B_constraints
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function, weightfun_pars=stvar$model$weightfun_pars)
  npars <- length(stvar$params)
  pars <- stvar$std_errors
  pars <- reform_constrained_pars(p=p, M=M, d=d, params=pars, weight_function=weight_function,
                                  weightfun_pars=weightfun_pars, cond_dist=cond_dist,
                                  identification=identification, AR_constraints=AR_constraints,
                                  mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                  B_constraints=B_constraints)
  all_phi0_or_mu <- pick_phi0(M=M, d=d, params=pars)
  all_A <- pick_allA(p=p, M=M, d=d, params=pars)
  distpars <- pick_distpars(params=pars, cond_dist=cond_dist)
  if(identification == "reduced_form") {
    all_Omega <- pick_Omegas(p=p, M=M, d=d, params=pars, identification=identification)
  } else {
    # No standard errors for cov. mats. as the model is parametrized with W and lambdas
    all_Omega <- array(" ", dim=c(d, d, M))
  }
  weightpars <- pick_weightpars(p=p, M=M, d=d, params=pars, weight_function=weight_function, weightfun_pars=weightfun_pars)
  if(weight_function == "relative_dens") {
    weightpars[M] <- NA # No standard error for the last alpha
  } else if(weight_function == "logistic") {
    # Weightpars are ok as is.
  } else if(weight_function == "mlogit") {
    all_gamma_m <- matrix(weightpars, ncol=M-1) # Column per gamma_m, m=1,...,M-1, gamma_M=0.
  } else {
    stop("Unkown weight function in print_std_errors")
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

  cat(weight_function, cond_dist, "STVAR model,",
      ifelse(identification == "reduced_form", "reduced form model",
             ifelse(identification == "recursive", "recursive identification,", paste0("identified by ", identification, ","))),
      ifelse(is.null(AR_constraints), "no AR_constraints,", "AR_constraints used,"),
      ifelse(is.null(mean_constraints), paste0("no mean_constraints,", ifelse(is.null(B_constraints), "", ",")),
             paste0("mean_constraints used,", ifelse(is.null(B_constraints), "", ","))),
      ifelse(identification %in% c("reduced_form", "recursive"), "", ifelse(is.null(B_constraints),
                                                                            "no B_constraints,", "B_constraints used,")))
  cat("\n", paste0(" p = ", p, ", "))
  cat(paste0("M = ", M, ", "))
  cat(paste0("d = ", d, ", #parameters = " , npars, ","))
  if(weight_function == "mlogit") {
    cat("\n ", paste0("Switching variables: ", paste0(var_names[weightfun_pars[[1]]], collapse=", "), " with ",
                      weightfun_pars[[2]], ifelse(weightfun_pars[[2]] == 1, " lag.", " lags.")))
  } else if(weight_function %in% c("logistic", "exponential", "threshold")) {
    cat("\n ", paste0("Switching variable: ", paste0(var_names[weightfun_pars[1]], collapse=", "), " with lag ",
                      weightfun_pars[2], "."))
  }
  cat("\n\n")
  cat("APPROXIMATE STANDARD ERRORS ASSUMING ASYMPTOTIC NORMALITY\n\n")

  left_brackets <- rep("[", times=d)
  right_brackets <- rep("]", times=d)
  plus <- c("+", rep(" ", d - 1))
  round_lbrackets <- rep("(", times=d)
  round_rbrackets <- rep(")", times=d)
  Y <- paste0("Y", 1:d)
  tmp_names <- paste0("tmp", 1:(p*(d + 2) + d + 2))

  for(m in seq_len(M)) {
    count <- 1
    cat(paste("Regime", m))
    cat("\n")
    if(cond_dist == "Student") {
      if(m == 1) {
        cat(paste0("Degrees of freedom: ", format_value(distpars), " (for all regimes)"), "\n")
      }
    }
    if(weight_function == "relative_dens") {
      if(m < M) cat(paste("Weight param:", format_value(weightpars[m])), "\n")
    } else if(weight_function == "mlogit") {
      if(m < M) cat(paste("Weight params:", paste0(format_value(all_gamma_m[,m]), collapse=", ")), "\n")
    } else if(weight_function %in% c("logistic", "exponential")) {
      if(m == M) {
        cat(paste("Weight params:", paste0(format_value(weightpars[1]), " (location), ",
                                           format_value(weightpars[2]), " (scale)")), "\n")
      }
    } else if(weight_function == "threshold") {
      if(m < M) {
        cat(paste0("Upper threshold: ", format_value(weightpars[m])), "\n")
      }
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

  if(!identification %in% c("reduced_form", "recursive")) {
    cat("Structural parameters:\n")
    if(identification == "heteroskedasticity") {
      W <- format_value(pick_W(p=p, M=M, d=d, params=pars, identification=identification))

      tmp <- c(rep(" ", times=d - 1), ",")
      df2 <- data.frame(left_brackets, W=W[,1])
      for(i1 in 2:d) {
        df2 <- cbind(df2, W[, i1])
        colnames(df2)[1 + i1] <- "tmp"
      }
      df2 <- cbind(df2, right_brackets)
      if(M > 1) {
        lambdas <- format_value(pick_lambdas(p=p, M=M, d=d, params=pars, identification=identification))
        tmp <- c(rep(" ", times=d - 1), ",")
        lambdas <- matrix(lambdas, nrow=d, ncol=M - 1, byrow=FALSE) # Column for each regime
        for(i1 in 1:(sum(M) - 1)) {
          lmb <- lambdas[,i1]
          df2 <- cbind(df2, tmp, left_brackets, lmb, right_brackets)
          colnames(df2)[grep("lmb", colnames(df2))] <- paste0("lamb", i1 + 1)
        }
      }
      names_to_omit <- unlist(lapply(c("left_brackets", "right_brackets", "tmp"), function(nam) grep(nam, colnames(df2))))
      colnames(df2)[names_to_omit] <- " "
      print(df2)
      cat("\n")
      n_zero <- sum(B_constraints == 0, na.rm=TRUE)
      n_free <- sum(is.na(B_constraints))
      n_sign <- d^2 - n_zero - n_free
      cat("The impact matrix is subject to", n_zero, "zero constraints and", n_sign, "sign constraints.\n")
      cat("\n")
    }
  }
  invisible(stvar)
}
