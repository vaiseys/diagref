# Profile likelihood confidence intervals for the weight parameters (theta).
# For the scalar case this is a 1D scan; for multi-parameter weights we
# profile each theta_j with the rest re-optimized at each grid point.

# Log-likelihood for Gaussian model (up to constant).
drm_loglik <- function(theta, Z_w, parsed) {
  w     <- compute_weights(theta, Z_w)
  inner <- drm_inner_fit(w, parsed)
  n     <- parsed$n
  -n / 2 * log(inner$rss / n)
}

# Profile CI for a single scalar theta (1D case).
# Returns list: ci_theta (on logit scale), ci_w (on [0,1] scale), profile_grid
profile_ci_scalar <- function(fit_obj, level = 0.95, n_grid = 200) {
  parsed    <- fit_obj$parsed
  Z_w       <- fit_obj$Z_w
  theta_hat <- fit_obj$theta_hat
  ll_hat    <- drm_loglik(theta_hat, Z_w, parsed)

  # Chi-squared cutoff (1 df)
  cutoff <- qchisq(level, df = 1) / 2

  theta_grid <- seq(
    max(theta_hat - 6, -10),
    min(theta_hat + 6,  10),
    length.out = n_grid
  )

  ll_grid <- vapply(theta_grid, function(th) {
    drm_loglik(th, Z_w, parsed)
  }, numeric(1))

  # Find crossing points (ll_hat - ll > cutoff)
  in_ci <- (ll_hat - ll_grid) <= cutoff

  if (!any(in_ci)) {
    warning("Profile likelihood CI not found within search range.", call. = FALSE)
    ci_theta <- c(NA_real_, NA_real_)
  } else {
    ci_theta <- c(theta_grid[which(in_ci)[1]], theta_grid[rev(which(in_ci))[1]])
  }

  ci_w <- plogis(ci_theta)

  list(
    ci_theta    = ci_theta,
    ci_w        = ci_w,
    profile_grid = tibble::tibble(
      theta = theta_grid,
      w     = plogis(theta_grid),
      loglik = ll_grid
    ),
    level = level
  )
}

# Profile CI for multi-parameter weight model.
# Profiles each theta_j independently (others re-optimized).
profile_ci_multi <- function(fit_obj, level = 0.95, n_grid = 100) {
  parsed    <- fit_obj$parsed
  Z_w       <- fit_obj$Z_w
  theta_hat <- fit_obj$theta_hat
  p_w       <- length(theta_hat)
  ll_hat    <- drm_loglik(theta_hat, Z_w, parsed)
  cutoff    <- qchisq(level, df = 1) / 2

  ci_list <- vector("list", p_w)

  for (j in seq_len(p_w)) {
    theta_j_grid <- seq(theta_hat[j] - 4, theta_hat[j] + 4, length.out = n_grid)

    ll_profile <- vapply(theta_j_grid, function(th_j) {
      # Re-optimize over all other theta with theta_j fixed
      if (p_w == 1) return(drm_loglik(th_j, Z_w, parsed))

      theta_init <- theta_hat[-j]
      opt <- optim(
        par    = theta_init,
        fn     = function(th_rest) {
          th_full <- numeric(p_w)
          th_full[j]  <- th_j
          th_full[-j] <- th_rest
          -drm_loglik(th_full, Z_w, parsed)
        },
        method = "BFGS"
      )
      -opt$value
    }, numeric(1))

    in_ci    <- (ll_hat - ll_profile) <= cutoff
    ci_theta <- if (any(in_ci)) {
      c(theta_j_grid[which(in_ci)[1]], theta_j_grid[rev(which(in_ci))[1]])
    } else {
      c(NA_real_, NA_real_)
    }

    ci_list[[j]] <- list(ci_theta = ci_theta, ci_w = plogis(ci_theta))
  }

  ci_list
}

# Dispatcher
compute_profile_ci <- function(fit_obj, level = 0.95, n_grid = 200) {
  if (length(fit_obj$theta_hat) == 1) {
    profile_ci_scalar(fit_obj, level = level, n_grid = n_grid)
  } else {
    profile_ci_multi(fit_obj, level = level, n_grid = n_grid)
  }
}
