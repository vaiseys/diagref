# Build the DRM design matrix columns for a given weight vector w_i.
# Uses sum-to-zero (contr.sum) coding: K-1 columns returned.
# var1, var2: integer indices (1-based) into all_levels
# w: numeric vector length n
# K: number of levels
# Returns n x (K-1) matrix
build_drm_matrix <- function(var1_int, var2_int, w, K) {
  C <- contr.sum(K)  # K x (K-1) sum-to-zero contrast matrix
  # DRM row for obs i: w_i * C[var1_i, ] + (1-w_i) * C[var2_i, ]
  w * C[var1_int, , drop = FALSE] + (1 - w) * C[var2_int, , drop = FALSE]
}

# Given theta (unconstrained), compute weight vector w (length n).
# Z_w: n x p weight model matrix; theta: p-vector
compute_weights <- function(theta, Z_w) {
  plogis(drop(Z_w %*% theta))
}

# Inner OLS fit for fixed weights w.
# Returns list: coef, rss, sigma2, fitted, residuals, df_resid
drm_inner_fit <- function(w, parsed) {
  K       <- length(parsed$all_levels)
  v1_int  <- as.integer(parsed$var1)
  v2_int  <- as.integer(parsed$var2)

  X_drm <- build_drm_matrix(v1_int, v2_int, w, K)
  X     <- cbind(X_drm, parsed$X_cov)

  fit <- lm.fit(X, parsed$y)

  n       <- parsed$n
  p       <- ncol(X)
  df_res  <- n - p
  rss     <- sum(fit$residuals^2)
  sigma2  <- rss / df_res

  list(
    coef      = fit$coefficients,
    rss       = rss,
    sigma2    = sigma2,
    fitted    = fit$fitted.values,
    residuals = fit$residuals,
    df_resid  = df_res,
    rank      = fit$rank,
    X         = X
  )
}

# Objective: RSS as a function of unconstrained theta (for optimizer).
drm_objective <- function(theta, Z_w, parsed) {
  w <- compute_weights(theta, Z_w)
  drm_inner_fit(w, parsed)$rss
}

# Recover all K sum-to-zero scores from the K-1 contrast coefficients.
# delta: (K-1)-vector of contrast coefficients
# Returns K-vector of scores (sums to zero)
recover_scores <- function(delta, K) {
  C <- contr.sum(K)
  drop(C %*% delta)
}

# Center non-intercept columns of Z_w for numerical stability.
# Returns list: Z_centered, col_means (for back-transforming theta later if needed)
center_weight_matrix <- function(Z_w) {
  col_means <- colMeans(Z_w)
  # Leave intercept column (all 1s) untouched
  is_intercept <- apply(Z_w, 2, function(col) all(col == 1))
  col_means[is_intercept] <- 0
  Z_centered <- sweep(Z_w, 2, col_means, "-")
  list(Z = Z_centered, means = col_means)
}

# Main fitting routine. Returns a raw fit object.
drm_fit <- function(parsed, Z_w) {
  p_w <- ncol(Z_w)

  # Center non-intercept columns for numerical stability. This transforms
  # the intercept to represent logit(w) at the mean of each predictor,
  # which avoids degeneracy when predictors are far from zero (e.g. raw age).
  centered   <- center_weight_matrix(Z_w)
  Z_c        <- centered$Z
  col_means  <- centered$means

  neg_ll <- function(theta) -drm_loglik(theta, Z_c, parsed)

  if (p_w == 1) {
    # Scalar w: 1D Brent with Hessian for SE
    opt <- optim(
      par     = 0,
      fn      = neg_ll,
      method  = "Brent",
      lower   = -10,
      upper   = 10,
      hessian = TRUE
    )
  } else {
    # Initialize intercept from scalar fit; other params start at 0
    scalar_opt <- optim(
      par    = 0,
      fn     = function(th) -drm_loglik(th, Z_c[, 1, drop = FALSE], parsed),
      method = "Brent",
      lower  = -10,
      upper  = 10
    )
    theta0 <- c(scalar_opt$par, rep(0, p_w - 1))
    opt <- optim(
      par     = theta0,
      fn      = neg_ll,
      method  = "BFGS",
      hessian = TRUE
    )
  }

  theta_c    <- opt$par
  vcov_c     <- tryCatch(solve(opt$hessian),
                         error = function(e) matrix(NA_real_, p_w, p_w))

  # Back-transform from centered to original scale so Z_w (uncentered)
  # can be used for all downstream predictions.
  # Slopes are unchanged; only intercept absorbs the centering shift.
  if (p_w > 1) {
    theta_hat <- theta_c
    theta_hat[1] <- theta_c[1] - sum(theta_c[-1] * col_means[-1])

    J           <- diag(p_w)
    J[1, -1]    <- -col_means[-1]
    vcov_theta  <- J %*% vcov_c %*% t(J)
  } else {
    theta_hat  <- theta_c
    vcov_theta <- vcov_c
  }

  w_hat <- compute_weights(theta_hat, Z_w)
  inner <- drm_inner_fit(w_hat, parsed)

  K     <- length(parsed$all_levels)
  n_drm <- K - 1
  delta  <- inner$coef[seq_len(n_drm)]
  scores <- recover_scores(delta, K)
  names(scores) <- parsed$all_levels

  cov_coef <- inner$coef[-(seq_len(n_drm))]

  list(
    theta_hat  = theta_hat,
    vcov_theta = vcov_theta,
    w_hat      = w_hat,
    scores     = scores,
    cov_coef   = cov_coef,
    inner      = inner,
    rss        = inner$rss,
    K          = K,
    parsed     = parsed,
    Z_w        = Z_w
  )
}
