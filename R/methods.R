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
    cat("Weights (w) by group:\n")
    unique_wd <- unique(x$weight_data)
    unique_wd <- unique_wd[do.call(order, as.list(unique_wd)), , drop = FALSE]
    Z_uniq    <- build_weight_matrix(x$weight_formula, unique_wd)
    w_vals    <- plogis(drop(Z_uniq %*% x$theta_hat))
    labels    <- vapply(seq_len(nrow(unique_wd)), function(i) {
      row_vals <- lapply(unique_wd[i, , drop = FALSE], as.character)
      paste(paste(names(row_vals), unlist(row_vals), sep = "="), collapse = ", ")
    }, character(1))
    has_group_ci <- !is.null(x$ci) && length(x$ci) == length(w_vals)
    for (i in seq_along(w_vals)) {
      if (has_group_ci) {
        ci_i <- x$ci[[i]]$ci_w
        cat(sprintf("  %s: w = %.4f  [%.4f, %.4f] (%d%% profile CI)\n",
                    labels[i], w_vals[i], ci_i[1], ci_i[2],
                    round(x$level * 100)))
      } else {
        cat(sprintf("  %s: w = %.4f\n", labels[i], w_vals[i]))
      }
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
  parsed <- x$parsed
  K      <- x$K
  n_drm  <- K - 1
  inner  <- x$inner
  df_r   <- inner$df_resid
  level  <- x$level
  t_crit <- qt(1 - (1 - level) / 2, df = df_r)
  z_crit <- qnorm(1 - (1 - level) / 2)

  # OLS variance-covariance for linear parameters
  sigma2  <- inner$sigma2
  XtX_inv <- tryCatch(
    chol2inv(chol(crossprod(inner$X))),
    error = function(e) solve(crossprod(inner$X))
  )
  se_coef <- sqrt(diag(XtX_inv) * sigma2)

  # --- Score parameters ---
  C         <- contr.sum(K)
  se_delta  <- se_coef[seq_len(n_drm)]
  var_mu    <- C %*% diag(se_delta^2, nrow = n_drm) %*% t(C)
  se_mu     <- sqrt(pmax(diag(var_mu), 0))
  t_scores  <- x$scores / se_mu

  score_rows <- tibble::tibble(
    term      = paste0("score[", parsed$all_levels, "]"),
    type      = "score",
    estimate  = x$scores,
    std.error = se_mu,
    statistic = t_scores,
    p.value   = 2 * pt(abs(t_scores), df = df_r, lower.tail = FALSE),
    conf.low  = x$scores - t_crit * se_mu,
    conf.high = x$scores + t_crit * se_mu
  )

  # --- Weight parameter(s) ---
  # Always show predicted weights on the probability scale for unique
  # combinations of weight predictors.
  vcov_w   <- x$vcov_theta
  has_vcov <- !is.null(vcov_w) && !anyNA(vcov_w)

  if (is.null(x$weight_formula)) {
    # Scalar w: use profile CI; SE via delta method from Hessian
    theta <- x$theta_hat
    w_val <- plogis(theta)
    se_w  <- if (has_vcov) sqrt(vcov_w[1, 1]) * w_val * (1 - w_val) else NA_real_
    ci_w  <- x$ci$ci_w

    w_rows <- tibble::tibble(
      term      = "w",
      type      = "weight",
      estimate  = w_val,
      std.error = se_w,
      statistic = NA_real_,
      p.value   = NA_real_,
      conf.low  = ci_w[1],
      conf.high = ci_w[2]
    )
  } else {
    # Multi-param weight model. Choose display based on the number of unique
    # predictor combinations:
    #   few (<=20): per-group predicted weights on the probability scale
    #   many (>20): logit-scale coefficients, like logistic regression output
    unique_wd  <- unique(x$weight_data)
    n_patterns <- nrow(unique_wd)

    if (n_patterns <= 20) {
      # --- Per-group predicted weights (probability scale) ---
      unique_wd <- unique_wd[do.call(order, as.list(unique_wd)), , drop = FALSE]
      Z_uniq   <- build_weight_matrix(x$weight_formula, unique_wd)
      lin_pred <- drop(Z_uniq %*% x$theta_hat)
      w_pred   <- plogis(lin_pred)

      # SE: delta method (approximate, shown for reference)
      if (has_vcov) {
        se_lin <- vapply(seq_len(nrow(Z_uniq)), function(i) {
          l <- Z_uniq[i, ]
          sqrt(max(0, as.numeric(t(l) %*% vcov_w %*% l)))
        }, numeric(1))
        se_w <- se_lin * w_pred * (1 - w_pred)
      } else {
        se_w <- rep(NA_real_, nrow(unique_wd))
      }

      # CIs: profile likelihood (stored in fit object during diagref())
      if (!is.null(x$ci) && length(x$ci) == nrow(unique_wd)) {
        conf_low  <- vapply(x$ci, function(ci) ci$ci_w[1], numeric(1))
        conf_high <- vapply(x$ci, function(ci) ci$ci_w[2], numeric(1))
      } else {
        # Fallback to delta method if profile CI unavailable
        conf_low  <- if (has_vcov) plogis(lin_pred - z_crit * se_lin) else rep(NA_real_, nrow(unique_wd))
        conf_high <- if (has_vcov) plogis(lin_pred + z_crit * se_lin) else rep(NA_real_, nrow(unique_wd))
      }

      labels <- vapply(seq_len(nrow(unique_wd)), function(i) {
        row_vals <- lapply(unique_wd[i, , drop = FALSE], as.character)
        paste(paste(names(row_vals), unlist(row_vals), sep = "="), collapse = ", ")
      }, character(1))

      w_rows <- tibble::tibble(
        term      = paste0("w[", labels, "]"),
        type      = "weight",
        estimate  = w_pred,
        std.error = se_w,
        statistic = NA_real_,
        p.value   = NA_real_,
        conf.low  = conf_low,
        conf.high = conf_high
      )
    } else {
      # --- Logit-scale coefficients (continuous or mixed predictors) ---
      # estimate/conf.low/conf.high are on the logit scale (like logistic
      # regression). Use plogis(estimate) to get the probability at a
      # given covariate value.
      se_theta <- if (has_vcov) sqrt(pmax(0, diag(vcov_w))) else rep(NA_real_, length(x$theta_hat))
      z_stat   <- x$theta_hat / se_theta
      p_val    <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)

      w_rows <- tibble::tibble(
        term      = paste0("w[", colnames(x$Z_w), "]"),
        type      = "weight",
        estimate  = x$theta_hat,
        std.error = se_theta,
        statistic = z_stat,
        p.value   = p_val,
        conf.low  = x$theta_hat - z_crit * se_theta,
        conf.high = x$theta_hat + z_crit * se_theta
      )
    }
  }

  # --- Covariate coefficients ---
  cov_coef <- x$cov_coef
  se_cov   <- se_coef[-(seq_len(n_drm))]
  t_cov    <- cov_coef / se_cov

  cov_rows <- tibble::tibble(
    term      = names(cov_coef),
    type      = "covariate",
    estimate  = cov_coef,
    std.error = se_cov,
    statistic = t_cov,
    p.value   = 2 * pt(abs(t_cov), df = df_r, lower.tail = FALSE),
    conf.low  = cov_coef - t_crit * se_cov,
    conf.high = cov_coef + t_crit * se_cov
  )

  do.call(rbind, list(score_rows, w_rows, cov_rows))
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
logLik.diagref_fit <- function(object, ...) {
  n      <- object$parsed$n
  # MLE of sigma^2 uses n not n-p in the denominator
  sigma2 <- object$rss / n
  ll     <- -n / 2 * (log(2 * pi * sigma2) + 1)

  K    <- object$K
  p_w  <- length(object$theta_hat)
  p_cov <- length(object$cov_coef)
  # free params: (K-1) scores + p_w weight params + p_cov linear params + 1 for sigma^2
  df   <- (K - 1L) + p_w + p_cov + 1L

  structure(ll, df = df, nobs = n, class = "logLik")
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
