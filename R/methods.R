#' @export
print.diagref_fit <- function(x, ...) {
  cat("Diagonal Reference Model\n")
  cat("Formula: "); print(x$formula)
  cat("\n")

  if (length(x$theta_hat) == 1) {
    cat(sprintf("Weight (w): %.4f", plogis(x$theta_hat)))
    if (!is.null(x$ci)) {
      cat(sprintf("  [%.4f, %.4f] (%d%% profile CI)\n",
                  x$ci$ci_w[1], x$ci$ci_w[2], round(x$level * 100)))
    } else cat("\n")
  } else {
    cat("Weights (w) per group:\n")
    w_vals <- plogis(x$theta_hat)
    for (i in seq_along(w_vals)) {
      cat(sprintf("  theta[%d]: w = %.4f\n", i, w_vals[i]))
    }
  }

  cat("\nScores (", x$parsed$var1_name, "/", x$parsed$var2_name, "):\n", sep = "")
  print(round(x$scores, 4))

  if (length(x$cov_coef) > 0) {
    cat("\nCovariate coefficients:\n")
    print(round(x$cov_coef, 4))
  }

  n       <- x$parsed$n
  K       <- x$K
  p_cov   <- length(x$cov_coef)
  p_total <- (K - 1) + length(x$theta_hat) + p_cov
  cat(sprintf("\nn = %d  |  df.resid = %d  |  sigma = %.4f\n",
              n, x$inner$df_resid, sqrt(x$inner$sigma2)))
  invisible(x)
}

#' @export
tidy.diagref_fit <- function(x, ...) {
  parsed  <- x$parsed
  K       <- x$K
  n_drm   <- K - 1
  inner   <- x$inner
  n       <- parsed$n

  # Variance-covariance for linear params from OLS
  sigma2  <- inner$sigma2
  XtX_inv <- tryCatch(
    chol2inv(chol(crossprod(inner$X))),
    error = function(e) solve(crossprod(inner$X))
  )
  se_coef <- sqrt(diag(XtX_inv) * sigma2)

  # --- Score parameters ---
  C       <- contr.sum(K)
  delta   <- inner$coef[seq_len(n_drm)]
  se_delta <- se_coef[seq_len(n_drm)]

  # Propagate delta SE to mu via C (K x K-1 matrix)
  # Var(mu) = C Var(delta) C'
  var_delta <- diag(se_delta^2, nrow = n_drm)
  var_mu    <- C %*% var_delta %*% t(C)
  se_mu     <- sqrt(pmax(diag(var_mu), 0))

  df_r <- inner$df_resid
  score_rows <- tibble::tibble(
    term      = paste0("score[", parsed$all_levels, "]"),
    type      = "score",
    estimate  = x$scores,
    std.error = se_mu,
    statistic = x$scores / se_mu,
    p.value   = 2 * pt(abs(x$scores / se_mu), df = df_r, lower.tail = FALSE)
  )

  # --- Weight parameter(s) ---
  p_w <- length(x$theta_hat)
  w_names <- if (!is.null(x$weight_formula)) {
    colnames(x$Z_w)
  } else {
    "w"
  }

  if (p_w == 1) {
    ci_w <- x$ci$ci_w
    w_rows <- tibble::tibble(
      term      = "w",
      type      = "weight",
      estimate  = plogis(x$theta_hat),
      std.error = NA_real_,
      statistic = NA_real_,
      p.value   = NA_real_,
      conf.low  = ci_w[1],
      conf.high = ci_w[2]
    )
  } else {
    ci_list <- x$ci
    w_rows <- tibble::tibble(
      term      = paste0("w[", w_names, "]"),
      type      = "weight",
      estimate  = plogis(x$theta_hat),
      std.error = NA_real_,
      statistic = NA_real_,
      p.value   = NA_real_,
      conf.low  = vapply(ci_list, function(ci) ci$ci_w[1], numeric(1)),
      conf.high = vapply(ci_list, function(ci) ci$ci_w[2], numeric(1))
    )
  }

  # --- Covariate coefficients ---
  cov_coef  <- x$cov_coef
  se_cov    <- se_coef[-(seq_len(n_drm))]

  cov_rows <- tibble::tibble(
    term      = names(cov_coef),
    type      = "covariate",
    estimate  = cov_coef,
    std.error = se_cov,
    statistic = cov_coef / se_cov,
    p.value   = 2 * pt(abs(cov_coef / se_cov), df = df_r, lower.tail = FALSE)
  )

  # Combine; fill missing columns with NA
  all_rows <- list(score_rows, w_rows, cov_rows)
  all_cols <- unique(unlist(lapply(all_rows, names)))
  for (i in seq_along(all_rows)) {
    missing <- setdiff(all_cols, names(all_rows[[i]]))
    for (col in missing) all_rows[[i]][[col]] <- NA_real_
    all_rows[[i]] <- all_rows[[i]][all_cols]
  }

  do.call(rbind, all_rows)
}

#' @export
glance.diagref_fit <- function(x, ...) {
  n      <- x$parsed$n
  K      <- x$K
  p_w    <- length(x$theta_hat)
  p_cov  <- length(x$cov_coef)
  p      <- (K - 1) + p_w + p_cov
  df_res <- x$inner$df_resid
  sigma  <- sqrt(x$inner$sigma2)
  rss    <- x$rss

  # Null model RSS (intercept only)
  y      <- x$parsed$y
  rss0   <- sum((y - mean(y))^2)
  r2     <- 1 - rss / rss0

  # F statistic
  p_extra <- p - 1  # vs intercept-only
  f_stat  <- ((rss0 - rss) / p_extra) / (rss / df_res)
  p_f     <- 1 - pf(f_stat, p_extra, df_res)

  tibble::tibble(
    r.squared     = r2,
    sigma         = sigma,
    statistic     = f_stat,
    p.value       = p_f,
    df            = p,
    df.residual   = df_res,
    nobs          = n
  )
}

#' @export
augment.diagref_fit <- function(x, ...) {
  mf <- x$parsed$mf
  tibble::as_tibble(mf) |>
    tibble::add_column(
      .fitted   = x$inner$fitted,
      .resid    = x$inner$residuals,
      .w        = if (length(x$w_hat) == 1) x$w_hat else x$w_hat,
      .drm_term = x$w_hat * x$scores[as.integer(x$parsed$var1)] +
                  (1 - x$w_hat) * x$scores[as.integer(x$parsed$var2)]
    )
}
