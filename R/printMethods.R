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
  weightfun_pars <- check_weightfun_pars(p=p, d=d, weight_function=weight_function, weightfun_pars=stvar$model$weightfun_pars)
  cond_dist <- stvar$model$cond_dist
  parametrization <- stvar$model$parametrization
  identification <- stvar$model$identification
  AR_constraints <- stvar$model$AR_constraints
  mean_constraints <- stvar$model$mean_constraints
  weight_constraints <- stvar$model$weight_constraints
  B_constraints <- stvar$model$B_constraints
  IC <- stvar$IC
  var_names <- colnames(stvar$data)
  if(is.null(var_names)) var_names <- paste0("Var.", 1:d)
  all_mu <- stvar$uncond_moments$regime_means
  all_sd <- sqrt(stvar$uncond_moments$regime_vars)
  npars <- length(params)
  T_obs <- ifelse(is.null(stvar$data), NA, nrow(stvar$data) - p)
  params <- reform_constrained_pars(p=p, M=M, d=d, params=params,
                                    weight_function=weight_function, weightfun_pars=weightfun_pars,
                                    cond_dist=cond_dist, identification=identification, AR_constraints=AR_constraints,
                                    mean_constraints=mean_constraints, weight_constraints=weight_constraints,
                                    B_constraints=B_constraints)
  if(stvar$model$parametrization == "mean") {
    params <- change_parametrization(p=p, M=M, d=d, params=params, AR_constraints=NULL,
                                     mean_constraints=NULL, change_to="intercept")
  }
  all_phi0 <- pick_phi0(M=M, d=d, params=params)
  all_A <- pick_allA(p=p, M=M, d=d, params=params)
  all_Omega <- pick_Omegas(p=p, M=M, d=d, params=params, identification=identification)
  weightpars <- pick_weightpars(p=p, M=M, d=d, params=params, weight_function=weight_function, weightfun_pars=weightfun_pars,
                                cond_dist=cond_dist)
  distpars <- pick_distpars(params=params, cond_dist=cond_dist)

  if(weight_function == "mlogit") {
    all_gamma_m <- cbind(matrix(weightpars, ncol=M-1), 0) # Column per gamma_m, m=1,...,M-1, gamma_M=0.
  }

  cat(weight_function, cond_dist, "STVAR model,",
      ifelse(identification == "reduced_form", "reduced form model",
             ifelse(identification == "recursive", "recursive identification,", paste0("identified by ", identification, ","))),
      ifelse(is.null(AR_constraints), "no AR_constraints,", "AR_constraints used,"),
      ifelse(is.null(mean_constraints), paste0("no mean_constraints,", ifelse(is.null(B_constraints), "", ",")),
             paste0("mean_constraints used,", ifelse(is.null(B_constraints), "", ","))),
      ifelse(identification %in% c("reduced_form", "recursive"), "",
             ifelse(is.null(B_constraints), "no B_constraints,", "B_constraints used,")))
  cat("\n", paste0(" p = ", p, ", "))
  cat(paste0("M = ", M, ", "))
  cat(paste0("d = ", d, ", #parameters = " , npars, ","),
      ifelse(is.na(T_obs), "\n", paste0("#observations = ", T_obs, " x ", d, "")))
  if(weight_function == "mlogit") {
    cat("\n ", paste0("Switching variables: ", paste0(var_names[weightfun_pars[[1]]], collapse=", "), " with ",
                     weightfun_pars[[2]], ifelse(weightfun_pars[[2]] == 1, " lag.", " lags.")))
  } else if(weight_function %in% c("logistic", "exponential", "threshold")) {
    cat("\n ", paste0("Switching variable: ", paste0(var_names[weightfun_pars[1]], collapse=", "), " with lag ",
                      weightfun_pars[2], "."))
  }
  cat("\n\n")

  if(summary_print) {
    all_boldA_eigens <- get_boldA_eigens(stvar)
    all_omega_eigens <- get_omega_eigens(stvar)
    form_val2 <- function(txt, val) paste(txt, format_value(val))
    cat(paste(form_val2("loglik/T:", stvar$loglik/T_obs),
              form_val2("AIC:", IC$AIC),
              form_val2("HQIC:", IC$HQIC),
              form_val2("BIC:", IC$BIC),
              sep=", "), "\n\n")
  }

  plus <- c("+", rep(" ", times=d-1))
  round_lbrackets <- rep("(", times=d)
  round_rbrackets <- rep(")", times=d)
  Y <- paste0("y", 1:d)
  tmp_names <- paste0("tmp", 1:(p*(d + 2) + d + 2))

  for(m in seq_len(M)) {
    count <- 1
    cat(paste("Regime", m), "\n")
    if(cond_dist == "Student") {
      if(m == 1) {
        cat(paste0("Degrees of freedom: ", format_value(distpars), " (for all regimes)"), "\n")
      }
    }
    if(summary_print) {
      cat(paste("Moduli of 'bold A' eigenvalues: ", paste0(format_value(all_boldA_eigens[,m]), collapse=", ")),"\n")
      cat(paste("Cov. matrix 'Omega' eigenvalues:", paste0(format_value(all_omega_eigens[,m]), collapse=", ")),"\n")
    }
    if(weight_function == "relative_dens") {
      cat(paste("Weight param:", format_value(weightpars[m])), "\n")
    } else if(weight_function == "mlogit") {
      if(m < M) {
        cat(paste("Weight params:", paste0(format_value(all_gamma_m[,m]), collapse=", ")), "\n")
      }
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
    cat("Regime means:", paste0(format_value(all_mu[,m]), collapse=", "))
    if(summary_print) {
      cat("\nRegime sdevs:", paste0(format_value(all_sd[,m]), collapse=", "))
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
    df[, tmp_names[p*(d + 2) + 1]] <- left_brackets
    df[, c("Omega", tmp_names[(p*(d + 2) + 2):(p*(d + 2) + d)])] <- format_value(all_Omega[, , m])
    df[, tmp_names[p*(d + 2) + d + 1]] <- rep("]", times=d)
    df[, "1/2"] <- rep(" ", d)
    df[, tmp_names[p*(d + 2) + d + 2]] <- paste0("eps", 1:d)
    names_to_omit <- unlist(lapply(c("plus", "eq", "round_lbrackets", "round_rbrackets", tmp_names),
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
  if(!identification %in% c("reduced_form", "recursive")) { # No separate struct pars for recursively identified models
    cat("Structural parameters:\n")
    if(identification == "heteroskedasticity") {
      W <- format_value(pick_W(p=p, M=M, d=d, params=params, identification=identification))
      tmp <- c(rep(" ", times=d - 1), ",")
      df2 <- data.frame(left_brackets, W=W[,1])
      for(i1 in 2:d) {
        df2 <- cbind(df2, W[, i1])
        colnames(df2)[1 + i1] <- "tmp"
      }
      df2 <- cbind(df2, right_brackets)
      if(M > 1) {
        lambdas <- format_value(pick_lambdas(p=p, M=M, d=d, params=params, identification=identification))
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
      n_free <- ifelse(is.null(B_constraints), d^2, sum(is.na(B_constraints)))
      n_sign <- d^2 - n_zero - n_free
      cat("The impact matrix is subject to", n_zero, "zero constraints and", n_sign, "sign constraints.\n")
      cat("\n")
    }
  }

  if(summary_print) {
    cat("Print approximate standard errors with the function 'print_std_errors'.\n")
  }
  invisible(stvar)
}


#' @title Summary print method from objects of class 'stvarsum'
#'
#' @description \code{print.stvarsum} is a print method for object \code{'stvarsum'} generated
#'   by \code{summary.stvar}.
#'
#' @param x object of class 'stvarsum' generated by \code{summary.stvar}.
#' @param ... currently not used.
#' @param digits the number of digits to be printed.
#' @export

print.stvarsum <- function(x, ..., digits) {
  stvarsum <- x
  if(missing(digits)) digits <- stvarsum$digits
  print.stvar(stvarsum$stvar, ..., digits=digits, summary_print=TRUE)
  invisible(stvarsum)
}



#' @describeIn predict.stvar print method
#'
#' @param x object of class \code{'stvarpred'}
#' @param digits the number of decimals to print
#' @param ... currently not used.
#' @export

print.stvarpred <- function(x, ..., digits=2) {
  stvarpred <- x
  stopifnot(digits >= 0 & digits%%1 == 0)
  format_value <- format_valuef(digits)

  cat(paste0("Point forecast by ", stvarpred$pred_type, ", ",
             "two-sided prediction intervals with levels ", paste(stvarpred$pi, collapse=", "), "."), "\n")
  cat(paste0("Forecast ", stvarpred$nsteps, " steps ahead, based on ", stvarpred$nsim, " simulations.\n"))

  cat("\n")
  q <- stvarpred$q
  pred_ints <- stvarpred$pred_ints
  pred <- stvarpred$pred
  pred_type <- stvarpred$pred_type
  series_names <- colnames(stvarpred$pred)
  for(i1 in seq_len(stvarpred$stvar$model$d)) {
    cat(paste0(series_names[i1], ":"), "\n")
    df <- as.data.frame(lapply(1:length(stvarpred$q), function(i2) format_value(pred_ints[, i2, i1])))
    names(df) <- q
    df[, pred_type] <- format_value(pred[,i1])
    new_order <- as.character(c(q[1:(length(q)/2)], pred_type, q[(length(q)/2 + 1):length(q)]))
    print(df[, new_order])
    cat("\n")
  }
  cat(paste("Point forecasts and prediction intervals for transition weights can be obtained with $trans_pred",
            "and $trans_pred_ints, respectively.\n"))
  invisible(stvarpred)
}


#' @title Print method for the class hypotest
#'
#' @description \code{print.hypotest} is the print method for the class hypotest
#'  objects.
#' @param digits how many significant digits to print?
#' @param x object of class \code{'hypotest'} generated by the function \code{Wald_test}, \code{LR_test},
#'    \code{Rao_test}, or \code{Portmenteau_test}.
#' @param ... currently not in use.
#' @export

print.hypotest <- function(x, ..., digits=4) {
  stopifnot(digits >= 0 & digits%%1 == 0)
  format_value <- function(a) format(a, digits=digits)
  type <- x$type
  if(type == "Portmanteau test") {
    type <- paste0(x$type, " (", x$which_test, ")")
  }
  cat(paste0(type, ":"), "\n",
      paste0("test stat = ", format_value(x$test_stat),
             ", df = ", x$df,
             ", p-value = ", format_value(x$p_value)))
  invisible(x)
}


#' @describeIn GIRF print method
#' @inheritParams print.stvarpred
#' @param x object of class \code{'girf'} generated by the function \code{GIRF}.
#' @param N_to_print an integer specifying the horizon how far to print the estimates and
#'   confidence intervals. The default is that all the values are printed.
#' @export

print.girf <- function(x, ..., digits=2, N_to_print) {
  girf <- x
  girf_res <- girf$girf_res
  stopifnot(digits >= 0 & digits%%1 == 0)
  format_value <- format_valuef(digits)
  if(missing(N_to_print)) {
    N_to_print <- nrow(girf_res[[1]]$point_est)
  } else {
    stopifnot(N_to_print %in% 1:nrow(girf_res[[1]]$point_est))
  }
  if(length(girf$which_cumulative) > 0) {
    cat(paste0("The responses of the variables ",
               paste0(dimnames(girf_res[[1]]$point_est)[[2]][girf$which_cumulative], collapse=", "),
               " were cumulated."), "\n\n")
  }

  for(i1 in 1:length(girf_res)) {
    if(i1 > 1) cat("------------------------\n")
    cat(paste0("The GIRF of shock ", girf$shocks[i1], ":"), "\n")
    girf_i1 <- girf_res[[i1]]
    for(i2 in 1:dim(girf_i1$conf_ints)[3]) {
      cat(paste0("The response of ", dimnames(girf_i1$conf_ints)[[3]][i2], ":"), "\n")
      df <- as.data.frame(lapply(1:ncol(girf_i1$conf_ints[, , i2]), function(i3) format_value(girf_i1$conf_ints[, i3, i2])))
      q <- dimnames(girf_i1$conf_ints)[[2]]
      names(df) <- q

      df[, "mean"] <- format_value(girf_i1$point_est[, i2])
      new_order <- as.character(c(q[1:(length(q)/2)], "mean", q[(length(q)/2 + 1):length(q)]))
      print(utils::head(df[, new_order], n=N_to_print + 1))
      cat("\n")
    }
  }
  invisible(girf)
}


#' @describeIn GFEVD print method
#' @inheritParams print.stvarpred
#' @param x object of class \code{'gfevd'} generated by the function \code{GFEVD}.
#' @param N_to_print an integer specifying the horizon how far to print the estimates.
#'   The default is that all the values are printed.
#' @export

print.gfevd <- function(x, ..., digits=2, N_to_print) {
  gfevd <- x
  gfevd_res <- gfevd$gfevd_res
  stopifnot(digits >= 0 & digits%%1 == 0)
  format_value <- format_valuef(digits)
  if(missing(N_to_print)) {
    N_to_print <- nrow(gfevd_res[, , 1]) - 1
  } else {
    stopifnot(N_to_print %in% 1:nrow(gfevd_res[, , 1]))
  }
  if(length(gfevd$which_cumulative) > 0) {
    cat(paste0("The responses of the variables ",
               paste0(dimnames(gfevd_res)[[3]][gfevd$which_cumulative], collapse=", "),
               " were cumulated."), "\n\n")
  }

  for(i1 in 1:dim(gfevd_res)[3]) { # Go through GFEVDs of each variable and possibly mixing weights
    if(i1 > 1) cat("------------------------\n")
    cat(paste0("The GFEVD for ", dimnames(gfevd_res)[[3]][i1], ":"), "\n")
    print(round(gfevd_res[1:(N_to_print  + 1), , i1], digits=digits))
    cat("\n")
  }
  invisible(gfevd)
}


#' @describeIn linear_IRF print method
#' @inheritParams print.girf
#' @param x object of class \code{'irf'} generated by the function \code{linear_IRF}.
#' @param N_to_print an integer specifying the horizon how far to print the estimates and
#'   confidence intervals. The default is that all the values are printed.
#' @param shocks_to_print the responses to which should should be printed?
#'   A numeric vector with elements in \code{1,...,d}. The default is that responses to
#'   all the shocks are printed.
#' @param ... currently not used.
#' @export

print.irf <- function(x, ..., digits=2, N_to_print, shocks_to_print) {
  irf <- x
  point_est <- irf$point_est
  stopifnot(digits >= 0 & digits%%1 == 0)
  format_value <- format_valuef(digits)
  if(missing(N_to_print)) {
    N_to_print <- dim(point_est)[3] # n slices
  } else {
    stopifnot(N_to_print %in% 1:irf$N)
  }
  if(missing(shocks_to_print)) {
    shocks_to_print <- 1:irf$stvar$model$d
  } else {
    stopifnot(all(shocks_to_print %in% irf$stvar$model$d))
  }
  var_names <- colnames(irf$stvar$data)
  if(length(irf$which_cumulative) > 0) {
    cat(paste0("The responses of the variables ",
               paste0(var_names[irf$which_cumulative], collapse=", "),
               " were cumulated."), "\n\n")
  }

  # Loop through the shocks
  for(i1 in shocks_to_print) {
    if(i1 > 1) cat("------------------------\n")
    cat(paste0("The IRF of Shock ", i1, ":"), "\n")
    irf_i1 <- as.data.frame(format_value(t(point_est[, i1, ]))) # irf[variable, shock, horizon] -> irf_i1[variable, horizon]
    colnames(irf_i1) <- var_names
    rownames(irf_i1) <- 0:(nrow(irf_i1) - 1)
    print(utils::head(irf_i1, n=N_to_print + 1))
    cat("\n")
  }

  # If conf_ints not null: let know where conf_ints are found
  if(!is.null(irf$conf_ints)) {
    cat("Print confidence intervals from $conf_ints")
  }
  invisible(irf)
}
