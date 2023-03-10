#' @title Function factory for value formatting
#'
#' @description \code{format_valuef} is a function factory for
#'   formatting values with certain number of digits.
#'
#' @param digits the number of decimals to print
#' @return Returns a function that takes an atomic vector as argument
#'   and returns it formatted to character with \code{digits} decimals.
#' @keywords internal

format_valuef <- function(digits) {
  function(x) tryCatch(format(round(x, digits), nsmall=digits), error=function(e) x)
}


#' @describeIn STVAR print method
#' @param x an object of class \code{'stvar'}.
#' @param ... currently not used.
#' @param digits number of digits to be printed.
#' @param summary_print if set to \code{TRUE} then the print
#'   will include log-likelihood and information criteria values.
#' @export

print.stvar <- function(x, ..., digits=2, summary_print=FALSE) {
  stvar <- x
  stopifnot(digits >= 0 & digits%%1 == 0)
  format_value <- format_valuef(digits)
  p <- stvar$model$p
  M <- stvar$model$M
  d <- stvar$model$d
  params <- stvar$params
  weight_function <- stvar$model$weight_function
  cond_dist <- stvar$model$cond_dist
  parametrization <- stvar$model$parametrization
  identification <- stvar$model$identification
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  B_constraints <- stvar$model$B_constraints
  IC <- stvar$IC
  all_mu <- round(get_regime_means(p=p, M=M, d=d, params=params, weight_function=weight_function,
                                   cond_dist=cond_dsit, parametrization=parametrization,
                                   identification=identification, AR_constraints=AR_constraints,
                                   mean_constraints=mean_constrains, B_constraints=B_constraints), digits)
  npars <- length(params)
  T_obs <- ifelse(is.null(stvar$data), NA, nrow(stvar$data))
  # REFORM CONSTRAINED PARS HERE
  if(!is.null(AR_constraints) || !is.null(mean_constraints) || !is.null(B_constraints)) {
    print("Constrained models are not yet implemented to get_regime_means")
    return(invisible(stvar))
  } else if(identification != "reduced_form") {
    print("Structural models are not yet implemented to get_regime_means")
    return(invisible(stvar))
  }
  if(stvar$model$parametrization == "mean") {
    params <- change_parametrization(p=p, M=M, d=d, params=params, AR_constraints=NULL,
                                     mean_constraints=NULL, change_to="intercept")
  }
  all_phi0 <- pick_phi0(M=M, d=d, params=params)
  all_A <- pick_allA(p=p, M=M, d=d, params=params)
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params)
  weightpars <- pick_weightpars(p=p, M=M, d=d, params=params, weight_function=weight_function,
                                cond_dist=cond_dist)
  # pick dist_pars
  cat(weight_function, cond_dist, "STVAR model, identification:", identification,"\n")
  cat(paste0(" p = ", p, ", "))
  cat(paste0("M = ", M, ", "))
  cat(paste0("d = ", d, ", #parameters = " , npars, ","),
      ifelse(is.na(T_obs), "\n", paste0("#observations = ", T_obs, " x ", d, ",\n")),
      "log-likelihood,", parametrization, "parametrization,",
      ifelse(is.null(AR_constraints), "no AR param constraints,", "AR params constrained,"),
      ifelse(is.null(mean_constraints), "no mean param constraints,", "mean params constrained,"),
      ifelse(identification == "reduced_form", "", ifelse(is.null(B_constraints),
                                                          "no B_constraints", "B_constraints employed")), "\n")
  cat("\n")

  # IMPLEMENT GET BOLDA EIGEN AND GET OMEGA EIGENS (with para args?)
  if(summary_print) {
    all_boldA_eigens <- get_boldA_eigens(stvar)
    all_omega_eigens <- get_omega_eigens(stvar)
    form_val2 <- function(txt, val) paste(txt, format_value(val))
    cat(paste(form_val2(" log-likelihood:", stvar$loglik),
              form_val2("AIC:", IC$AIC),
              form_val2("HQIC:", IC$HQIC),
              form_val2("BIC:", IC$BIC),
              sep=", "), "\n\n")
  }

  plus <- c("+", rep(" ", times=d-1))
  arch_scalar <- c(rep(" ", times=d-1), "ARCH_mt")
  round_lbrackets <- rep("(", times=d)
  round_rbrackets <- rep(")", times=d)
  Y <- paste0("y", 1:d)
  tmp_names <- paste0("tmp", 1:(p*(d + 2) + d + 2))

  for(m in seq_len(sum(M))) {
    count <- 1
    if(model == "GMVAR") {
      regime_type <- "GMVAR"
    } else if(model == "StMVAR") {
      regime_type <- "StMVAR"
      M1 <- 0
    } else {
      M1 <- M[1]
      regime_type <- ifelse(m <= M1, "GMVAR", "StMVAR")
    }
    cat(paste("Regime", m))
    if(model == "G-StMVAR") cat(paste0(" (", regime_type, " type)"))
    cat("\n")
    if(summary_print) {
      cat(paste("Moduli of 'bold A' eigenvalues: ", paste0(format_value(all_boldA_eigens[,m]), collapse=", ")),"\n")
      cat(paste("Cov. matrix 'Omega' eigenvalues:", paste0(format_value(all_omega_eigens[,m]), collapse=", ")),"\n")
    }
    cat(paste("Mixing weight:", format_value(alphas[m])), "\n")
    cat("Regime means:", paste0(format_value(all_mu[,m]), collapse=", "))
    if(regime_type == "StMVAR") { # Print degrees of freedom parameter for StMVAR type regimes
      cat("\nDf parameter: ", format_value(all_df[m - M1]))
    }
    cat("\n\n")

    left_brackets <- rep("[", times=d)
    right_brackets <- rep("]", times=d)
    df <- data.frame(Y=Y,
                     eq=c("=", rep(" ", d - 1)),
                     eq=left_brackets,
                     phi0=format_value(all_phi0[, m, drop=FALSE]),
                     eq=rep("]", times=d),
                     plus)
    for(i1 in seq_len(p)) {
      Amp_colnames <- c(paste0("A", i1), tmp_names[count:(count + d - 1 - 1)]); count <- count + d - 1
      df[, tmp_names[count]] <- left_brackets; count <- count + 1
      df[, Amp_colnames] <- format_value(all_A[, ,i1 , m])
      df[, tmp_names[count]] <- rep("]", times=d); count <- count + 1
      df[, tmp_names[count]] <- paste0(Y, ".", i1); count <- count + 1
      df <- cbind(df, plus)
    }
    if(regime_type == "StMVAR") { # Time varying ARCH scalar multiplying the constant part of error term covariance matrix
      df <- cbind(df, round_lbrackets, arch_scalar)
    }
    df[, tmp_names[p*(d + 2) + 1]] <- left_brackets
    df[, c("Omega", tmp_names[(p*(d + 2) + 2):(p*(d + 2) + d)])] <- format_value(all_Omega[, , m])
    df[, tmp_names[p*(d + 2) + d + 1]] <- rep("]", times=d)
    if(regime_type == "StMVAR") {
      df <- cbind(df, round_rbrackets)
    }
    df[, "1/2"] <- rep(" ", d)
    df[, tmp_names[p*(d + 2) + d + 2]] <- paste0("eps", 1:d)
    names_to_omit <- unlist(lapply(c("plus", "eq", "arch_scalar", "round_lbrackets", "round_rbrackets", tmp_names),
                                   function(nam) grep(nam, colnames(df))))
    colnames(df)[names_to_omit] <- " "
    print(df)
    cat("\n")
    if(summary_print) {
      cat("Error term correlation matrix:\n")
      print(cov2cor(all_Omega[, , m]), digits=digits)
      cat("\n")
    }
  }
  if(!is.null(structural_pars)) {
    cat("Structural parameters:\n")
    W <- format_value(pick_W(p=p, M=M, d=d, params=params, structural_pars=structural_pars))

    tmp <- c(rep(" ", times=d - 1), ",")
    df2 <- data.frame(left_brackets, W=W[,1])
    for(i1 in 2:d) {
      df2 <- cbind(df2, W[, i1])
      colnames(df2)[1 + i1] <- "tmp"
    }
    df2 <- cbind(df2, right_brackets)
    if(sum(M) > 1) {
      lambdas <- format_value(pick_lambdas(p=p, M=M, d=d, params=params, structural_pars=structural_pars))
      tmp <- c(rep(" ", times=d - 1), ",")
      lambdas <- matrix(lambdas, nrow=d, ncol=sum(M) - 1, byrow=FALSE) # Column for each regime
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
    W_orig <- stvar$model$structural_pars$W
    n_zero <- sum(W_orig == 0, na.rm=TRUE)
    n_free <- sum(is.na(W_orig))
    n_sign <- d^2 - n_zero - n_free
    cat("The B-matrix (or equally W) is subject to", n_zero, "zero constraints and", n_sign, "sign constraints.\n")
    cat("The eigenvalues lambda_{mi} are", ifelse(is.null(stvar$model$structural_pars$C_lambda), "not subject to linear constraints.",
                                                  "subject to linear constraints."))
    cat("\n")
  }

  if(summary_print) {
    cat("Print approximate standard errors with the function 'print_std_errors'.\n")
  }
  invisible(stvar)
}
