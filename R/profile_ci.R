# Profile-likelihood confidence intervals for the weight model.
#
# For a linear combination t = a' theta, the CI is the set of t whose
# likelihood-ratio statistic n * log(RSS_profile(t) / RSS_hat) is below
# qchisq(level, 1). RSS_profile(t) minimises RSS over theta subject to
# a' theta = t. As in MASS::profile(), the profile is traced continuously
# outward from the estimate, so the interval is the region connected to it.
# (The weight model can have degenerate far-away optima, e.g. w becoming a
# step function of a continuous predictor; those are not followed.)
#
# Each endpoint is bracketed by stepping out from the estimate in doubling
# steps, then located with uniroot(). If the statistic has not crossed the
# cutoff once the move shifts logit(w) by more than 40 for some observation,
# their weights are saturated at 0 or 1, the likelihood is flat from there
# on, and the endpoint is infinite.
profile_ci_weight <- function(fit, a, level = 0.95) {
  prof <- weight_profile(fit, a)
  crit <- qchisq(level, df = 1)
  lr_minus_crit <- function(t) prof$lr(t) - crit

  # Step size for bracketing: one standard error, capped at 1 because the
  # Wald SE is huge (and useless) when a weight is near the boundary.
  se <- sqrt(drop(crossprod(a, fit$vcov_theta %*% a)))
  h  <- if (is.finite(se) && se > 0) min(se, 1) else 1
  # Largest move worth trying: shifting t by max_step moves logit(w) by 40
  # for the observation with the largest |z'a|.
  max_step <- 40 / max(abs(fit$data$Z %*% a))

  endpoint <- function(direction) {
    prof$restart()
    inner   <- prof$t_hat
    f_inner <- -crit
    step    <- h
    repeat {
      if (step > 2 * max_step) return(direction * Inf)
      outer   <- prof$t_hat + direction * min(step, max_step)
      f_outer <- lr_minus_crit(outer)
      if (f_outer > 0) {
        lo <- min(inner, outer)
        hi <- max(inner, outer)
        return(uniroot(lr_minus_crit, c(lo, hi),
                       f.lower = if (lo == inner) f_inner else f_outer,
                       f.upper = if (hi == inner) f_inner else f_outer,
                       tol = 1e-6 * h)$root)
      }
      inner   <- outer
      f_inner <- f_outer
      step    <- 2 * step
    }
  }

  ci <- c(endpoint(-1), endpoint(1))
  prof$check()
  ci
}

# Likelihood-ratio test of a' theta = t0 against the fitted model.
profile_lr_test <- function(fit, a, t0 = 0) {
  prof <- weight_profile(fit, a)
  lr   <- max(prof$lr(t0), 0)
  prof$check()
  c(statistic = sign(prof$t_hat - t0) * sqrt(lr),
    p.value   = pchisq(lr, df = 1, lower.tail = FALSE))
}

# Profile likelihood-ratio statistic for a' theta, as a function of t:
#   lr(t) = n * log(RSS_profile(t) / RSS_hat),
# where RSS_profile(t) minimises RSS subject to a' theta = t. With several
# weight parameters, theta_s = t * b / |b|^2 + N alpha (standardised scale),
# where N spans the null space of b, and alpha is optimised from (i) the
# previous solution, which follows the profile path, and (ii) whichever of
# the profile_starts() candidates starts lowest; the better result is kept.
weight_profile <- function(fit, a) {
  d     <- fit$data
  A     <- fit$A
  A_inv <- solve(A)
  Zs    <- d$Z %*% A
  b     <- drop(crossprod(A, a))            # a' theta = b' theta_s
  obj   <- drm_objective(Zs, d)
  alpha_prev     <- NULL
  better_optimum <- FALSE

  if (length(b) == 1L) {
    profile_log_rss <- function(t) obj$fn(t / b)
  } else {
    N      <- qr.Q(qr(b), complete = TRUE)[, -1L, drop = FALSE]
    base   <- b / sum(b^2)
    starts <- profile_starts(fit, a)
    profile_log_rss <- function(t) {
      theta_s <- function(al) t * base + drop(N %*% al)
      fn <- function(al) obj$fn(theta_s(al))
      gr <- function(al) drop(crossprod(N, obj$gr(theta_s(al))))
      cands <- lapply(starts(t), function(th) drop(crossprod(N, A_inv %*% th)))
      cands <- cands[which.min(vapply(cands, fn, numeric(1)))]
      if (!is.null(alpha_prev)) cands <- c(cands, list(alpha_prev))
      best <- NULL
      for (al in cands) {
        o <- optim(al, fn, gr, method = "BFGS",
                   control = list(maxit = 1000, reltol = 1e-10))
        if (is.null(best) || o$value < best$value) best <- o
      }
      alpha_prev <<- best$par
      best$value
    }
  }

  list(
    t_hat = sum(a * fit$theta),
    # A small negative LR is expected when the estimate sits on a flat ridge
    # towards w = 0 or 1 (the optimiser stops once the gradient vanishes);
    # only a meaningful improvement signals a different optimum.
    lr = function(t) {
      lr <- d$n * (profile_log_rss(t) - log(fit$rss))
      if (lr < -0.01) better_optimum <<- TRUE
      lr
    },
    restart = function() alpha_prev <<- NULL,
    check = function() {
      if (better_optimum) {
        warning("While profiling, a better fit was found than the reported ",
                "estimate, so the weight model may have several optima and ",
                "the estimates may be a local one.", call. = FALSE)
      }
    }
  )
}

# Starting values for the nuisance parameters when profiling a' theta = t.
# Returns a function of t giving a list of theta vectors (original scale),
# each equal to theta_hat moved along a direction v with a'v = 1:
#  - least squares: the move that changes the observations' logit(w) least
#    (in the sum-of-squares sense);
#  - for a single coefficient j of a model with an intercept: moves that
#    pivot logit(w) around the smallest and the largest observed z_j. For a
#    dummy this leaves one group's weight unchanged and moves only the
#    other's, which matters when that group's w is near 0 or 1: a move that
#    pushes everyone's w to the boundary flattens the likelihood and stalls
#    the optimiser.
profile_starts <- function(fit, a) {
  Z     <- fit$data$Z
  A     <- fit$A
  th    <- fit$theta
  t_hat <- sum(a * th)

  # Least-squares direction, computed in the well-conditioned scale
  b    <- drop(crossprod(A, a))
  Mb   <- solve(crossprod(Z %*% A), b)
  dirs <- list(drop(A %*% Mb) / sum(b * Mb))

  int <- which(colnames(Z) == "(Intercept)")
  j   <- which(a != 0)
  if (length(int) == 1L && length(j) == 1L && j != int) {
    for (pivot in range(Z[, j])) {
      v <- numeric(length(a))
      v[j]   <- 1 / a[j]
      v[int] <- -pivot / a[j]
      dirs <- c(dirs, list(v))
    }
  }
  function(t) lapply(dirs, function(v) th + (t - t_hat) * v)
}
