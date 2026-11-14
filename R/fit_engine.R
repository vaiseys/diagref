# Fitting engine.
#
# For fixed weights w the DRM is linear in the scores and the covariate
# coefficients, so those are profiled out by OLS and only the weight-model
# coefficients theta are optimised numerically (variable projection). The
# objective is log(RSS), which is -2/n times the Gaussian log-likelihood with
# sigma^2 profiled out, plus a constant. Its gradient is analytic: by the
# envelope theorem,
#   dRSS / dtheta = -2 G' r,   G = d mu / d theta = w (1 - w) (s[a] - s[b]) Z,
# where s are the scores and r the residuals at the inner OLS solution.

# OLS fit of the linear part for fixed theta.
drm_linear_fit <- function(theta, Z, d) {
  w  <- plogis(drop(Z %*% theta))
  X  <- cbind(d$C2 + w * d$D, d$X_cov)
  qx <- qr(X)
  if (qx$rank < ncol(X)) {
    stop("The design matrix is rank-deficient; check for covariates that are ",
         "collinear with the dyadic variables.", call. = FALSE)
  }
  resid <- qr.resid(qx, d$y)
  list(w = w, X = X, coef = qr.coef(qx, d$y), resid = resid,
       rss = sum(resid^2))
}

# n x p Jacobian of the mean with respect to theta.
dmu_dtheta <- function(lin, Z, d) {
  delta <- lin$coef[seq_len(d$K - 1L)]
  (lin$w * (1 - lin$w) * drop(d$D %*% delta)) * Z
}

# Objective and gradient in theta for weight design Z. optim() calls fn and
# gr at the same point, so the last inner fit is cached.
drm_objective <- function(Z, d) {
  last_theta <- NULL
  last_fit   <- NULL
  at <- function(theta) {
    if (!identical(theta, last_theta)) {
      last_fit   <<- drm_linear_fit(theta, Z, d)
      last_theta <<- theta
    }
    last_fit
  }
  list(
    fn = function(theta) log(at(theta)$rss),
    gr = function(theta) {
      lin <- at(theta)
      -2 * drop(crossprod(dmu_dtheta(lin, Z, d), lin$resid)) / lin$rss
    }
  )
}

drm_fit <- function(d) {
  Z      <- d$Z
  p      <- ncol(Z)
  scalar <- is_scalar_weight(Z)

  # Single weight: coarse grid, then Brent on the best bracket. The grid
  # guards against local optima and is cheap (one OLS fit per point).
  one  <- drm_objective(intercept_matrix(d$n), d)
  grid <- seq(-10, 10, by = 0.5)
  i    <- which.min(vapply(grid, one$fn, numeric(1)))
  bracket  <- grid[c(max(i - 1L, 1L), min(i + 1L, length(grid)))]
  theta_1  <- optimize(one$fn, bracket, tol = 1e-10)$minimum

  A <- weight_reparam(Z)
  if (scalar) {
    theta <- theta_1
  } else {
    # Optimise in the standardised parameterisation, starting from the
    # single-weight solution (projected onto the weight design).
    Zs    <- Z %*% A
    obj   <- drm_objective(Zs, d)
    start <- qr.coef(qr(Zs), rep(theta_1, d$n))
    opt   <- optim(start, obj$fn, obj$gr, method = "BFGS",
                   control = list(maxit = 1000, reltol = 1e-12))
    if (opt$convergence != 0L) {
      warning("The optimiser did not converge (code ", opt$convergence, ").",
              call. = FALSE)
    }
    theta <- drop(A %*% opt$par)
  }
  names(theta) <- colnames(Z)

  # |logit(w)| > 8 means w is within 0.0003 of 0 or 1. The likelihood is then
  # flat in that direction, so the estimate is wherever the optimiser
  # stopped, its Wald SE is meaningless, and its profile CI is open-ended.
  if (any(abs(Z %*% theta) > 8)) {
    warning("Some fitted weights are at the boundary (w = 0 or 1). The ",
            "weight estimates are unstable there; use the profile ",
            "confidence intervals, which may be infinite.", call. = FALSE)
  }

  lin <- drm_linear_fit(theta, Z, d)
  inf <- drm_vcov(lin, Z, d)

  # Internal parameters are (delta, gamma, theta). Reported parameters are
  # (scores, weight, gamma), where weight is w = plogis(theta) for a single
  # weight and the logit-scale theta otherwise. Tr is the Jacobian of that
  # map, for the delta method.
  K1    <- d$K - 1L
  q     <- ncol(d$X_cov)
  delta <- lin$coef[seq_len(K1)]
  gamma <- lin$coef[K1 + seq_len(q)]
  names(gamma) <- colnames(d$X_cov)
  scores <- drop(d$C %*% delta)
  names(scores) <- d$levels

  w_est <- if (scalar) plogis(theta) else theta
  Tr <- matrix(0, d$K + p + q, K1 + q + p)
  Tr[seq_len(d$K), seq_len(K1)] <- d$C
  Tr[d$K + seq_len(p), K1 + q + seq_len(p)] <-
    if (scalar) w_est * (1 - w_est) else diag(p)
  Tr[d$K + p + seq_len(q), K1 + seq_len(q)] <- diag(q)

  V0 <- inf$vcov
  V0[is.na(V0)] <- 0
  V  <- Tr %*% V0 %*% t(Tr)
  w_rows <- d$K + seq_len(p)
  if (anyNA(inf$vcov)) V[w_rows, ] <- V[, w_rows] <- NA_real_

  est <- c(scores, w_est, gamma)
  names(est) <- c(paste0("score[", d$levels, "]"),
                  if (scalar) "w" else paste0("w[", colnames(Z), "]"),
                  colnames(d$X_cov))
  dimnames(V) <- list(names(est), names(est))

  theta_rows <- K1 + q + seq_len(p)
  list(
    coefficients  = est,
    vcov          = V,
    term_type     = rep(c("score", "weight", "covariate"), c(d$K, p, q)),
    scores        = scores,
    theta         = theta,
    vcov_theta    = inf$vcov[theta_rows, theta_rows, drop = FALSE],
    cov_coef      = gamma,
    w             = lin$w,
    fitted        = d$y - lin$resid,
    residuals     = lin$resid,
    rss           = lin$rss,
    sigma         = sqrt(inf$sigma2),
    df.residual   = inf$df.residual,
    A             = A,
    scalar_weight = scalar
  )
}

# Joint asymptotic covariance of (delta, gamma, theta): sigma^2 (J'J)^-1,
# with J = [X, G] the Jacobian of the mean (as for any nonlinear least
# squares fit; cf. nls()). If theta is not identified (J rank-deficient),
# fall back to the covariance of the linear part conditional on w.
drm_vcov <- function(lin, Z, d) {
  J      <- cbind(lin$X, dmu_dtheta(lin, Z, d))
  df     <- d$n - ncol(J)
  sigma2 <- lin$rss / df

  V    <- matrix(NA_real_, ncol(J), ncol(J))
  Jinv <- crossprod_inverse(J)
  if (!is.null(Jinv)) {
    V <- sigma2 * Jinv
  } else {
    warning("The weight is not identified (the dyadic variables never ",
            "differ, or w is at 0 or 1); its standard error is unavailable.",
            call. = FALSE)
    k <- ncol(lin$X)
    V[seq_len(k), seq_len(k)] <- sigma2 * crossprod_inverse(lin$X)
  }
  list(vcov = V, sigma2 = sigma2, df.residual = df)
}

# (M'M)^-1 via QR, or NULL if M is rank-deficient.
crossprod_inverse <- function(M) {
  qm <- qr(M)
  if (qm$rank < ncol(M)) return(NULL)
  out <- matrix(0, ncol(M), ncol(M))
  out[qm$pivot, qm$pivot] <- chol2inv(qr.R(qm))
  out
}
