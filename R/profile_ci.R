# Profile likelihood confidence intervals for weight parameters.

# Log-likelihood for Gaussian model (up to constant).
drm_loglik <- function(theta, Z_w, parsed) {
  w     <- compute_weights(theta, Z_w)
  inner <- drm_inner_fit(w, parsed)
  n     <- parsed$n
  -n / 2 * log(inner$rss / n)
}

# Profile CI for a single scalar theta (1D case).
profile_ci_scalar <- function(raw, level = 0.95, n_grid = 200) {
  parsed    <- raw$parsed
  Z_w       <- raw$Z_w
  theta_hat <- raw$theta_hat
  ll_hat    <- drm_loglik(theta_hat, Z_w, parsed)
  cutoff    <- qchisq(level, df = 1) / 2

  theta_grid <- seq(max(theta_hat - 6, -10), min(theta_hat + 6, 10),
                    length.out = n_grid)
  ll_grid    <- vapply(theta_grid,
                       function(th) drm_loglik(th, Z_w, parsed), numeric(1))

  in_ci <- (ll_hat - ll_grid) <= cutoff
  if (!any(in_ci)) {
    warning("Profile CI not found within search range.", call. = FALSE)
    ci_theta <- c(NA_real_, NA_real_)
  } else {
    ci_theta <- c(theta_grid[which(in_ci)[1]], theta_grid[rev(which(in_ci))[1]])
  }

  list(
    ci_theta     = ci_theta,
    ci_w         = plogis(ci_theta),
    profile_grid = tibble::tibble(theta  = theta_grid,
                                  w      = plogis(theta_grid),
                                  loglik = ll_grid),
    level        = level
  )
}

# Profile CI for a predicted weight corresponding to a linear combination
# eta = z' * theta.  For a group with weight-model row z, this gives a CI
# on the probability scale that is bounded in [0,1] and asymmetric when w
# is near the boundary — unlike the delta method approximation.
#
# Uses a null-space parameterisation: for fixed eta, theta lives on the
# (p_w-1)-dimensional affine subspace {theta : z'theta = eta}.  We optimise
# over that subspace with warm-starting from the centre of the grid outward.
profile_ci_lc <- function(raw, z, level = 0.95, n_grid = 200) {
  parsed    <- raw$parsed
  Z_w       <- raw$Z_w
  theta_hat <- raw$theta_hat
  p_w       <- length(theta_hat)

  eta_hat <- as.numeric(z %*% theta_hat)
  ll_hat  <- drm_loglik(theta_hat, Z_w, parsed)
  cutoff  <- qchisq(level, df = 1) / 2

  eta_grid <- seq(eta_hat - 6, eta_hat + 6, length.out = n_grid)
  ll_grid  <- numeric(n_grid)

  if (p_w == 1) {
    ll_grid <- vapply(eta_grid,
                      function(eta) drm_loglik(eta, Z_w, parsed), numeric(1))
  } else {
    z_sq   <- sum(z^2)
    Q      <- qr.Q(qr(matrix(z, ncol = 1)), complete = TRUE)
    N      <- Q[, -1, drop = FALSE]          # p_w x (p_w-1) null-space basis
    alpha0 <- drop(t(N) %*% theta_hat)       # starting alpha at MLE

    # Constrained -loglik for a given eta and free parameters alpha
    neg_ll_c <- function(eta, alpha) {
      theta <- (eta / z_sq) * z + drop(N %*% alpha)
      -drm_loglik(theta, Z_w, parsed)
    }

    opt_alpha <- function(eta, alpha_start) {
      if (length(alpha_start) == 1L) {
        opt <- optimize(function(a) neg_ll_c(eta, a),
                        interval = alpha_start + c(-5, 5))
        list(val = -opt$objective, par = opt$minimum)
      } else {
        opt <- optim(alpha_start, function(a) neg_ll_c(eta, a),
                     method = "BFGS", control = list(maxit = 100))
        list(val = -opt$value, par = opt$par)
      }
    }

    # Warm-start scanning outward from the centre in both directions
    centre <- which.min(abs(eta_grid - eta_hat))

    alpha_cur <- alpha0
    for (k in centre:n_grid) {
      res          <- opt_alpha(eta_grid[k], alpha_cur)
      ll_grid[k]   <- res$val
      alpha_cur    <- res$par
    }

    alpha_cur <- alpha0
    for (k in seq_len(centre - 1L)) {
      idx          <- centre - k
      res          <- opt_alpha(eta_grid[idx], alpha_cur)
      ll_grid[idx] <- res$val
      alpha_cur    <- res$par
    }
  }

  in_ci <- (ll_hat - ll_grid) <= cutoff
  if (!any(in_ci)) {
    warning("Profile CI not found within search range; try larger n_grid.",
            call. = FALSE)
    ci_eta <- c(NA_real_, NA_real_)
  } else {
    ci_eta <- c(eta_grid[which(in_ci)[1]], eta_grid[rev(which(in_ci))[1]])
  }

  list(ci_w = plogis(ci_eta), ci_eta = ci_eta, level = level)
}

# Dispatcher — chooses the right CI strategy based on the weight model.
# rep_data: representative display rows from make_representative_data(); NULL
# for purely continuous weight predictors (no categorical grouping exists).
compute_profile_ci <- function(raw, rep_data = NULL, weight_formula = NULL,
                                level = 0.95, n_grid = 200) {
  p_w <- length(raw$theta_hat)

  if (p_w == 1) {
    profile_ci_scalar(raw, level = level, n_grid = n_grid)
  } else if (!is.null(rep_data) && nrow(rep_data) <= 20) {
    # Profile CI for each representative point's predicted w = plogis(z'theta)
    Z_rep <- build_weight_matrix(weight_formula, rep_data)
    lapply(seq_len(nrow(Z_rep)), function(i) {
      profile_ci_lc(raw, Z_rep[i, ], level = level, n_grid = n_grid)
    })
  } else {
    # Purely continuous predictors: tidy() uses delta method on logit scale.
    NULL
  }
}
