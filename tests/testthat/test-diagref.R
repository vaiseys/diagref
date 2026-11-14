sim_drm <- function(n = 800, K = 4, seed = 10, group_w = FALSE, beta_x = 0,
                    theta = c(0.4, 1.2)) {
  set.seed(seed)
  orig   <- sample(seq_len(K), n, replace = TRUE)
  dest   <- sample(seq_len(K), n, replace = TRUE)
  g      <- rbinom(n, 1, 0.5)
  age    <- round(runif(n, 18, 89))
  x      <- rnorm(n)
  scores <- seq(0, 3, length.out = K)
  w      <- if (group_w) plogis(theta[1] + theta[2] * g) else 0.6
  y      <- w * scores[orig] + (1 - w) * scores[dest] + beta_x * x + rnorm(n)
  data.frame(y, orig, dest, g, age, x)
}

score_names <- function(fit) grep("^score\\[", names(coef(fit)), value = TRUE)

# ---- Estimation -------------------------------------------------------------

test_that("diagref recovers the weight and the scores", {
  fit <- diagref(y ~ orig + dest, data = sim_drm(n = 3000))
  true_scores <- seq(0, 3, length.out = 4)
  expect_equal(unname(coef(fit)["w"]), 0.6, tolerance = 0.05)
  expect_equal(unname(fit$scores), true_scores - mean(true_scores),
               tolerance = 0.1)
})

test_that("scores sum to zero", {
  fit <- diagref(y ~ orig + dest, data = sim_drm(K = 3))
  expect_equal(sum(coef(fit)[score_names(fit)]), 0, tolerance = 1e-10)
})

test_that("covariates are estimated", {
  fit <- diagref(y ~ orig + dest + x, data = sim_drm(n = 1000, beta_x = 1.5))
  expect_equal(unname(coef(fit)["x"]), 1.5, tolerance = 0.1)
})

test_that("a weight model recovers group-specific weights", {
  fit <- diagref(y ~ orig + dest, data = sim_drm(n = 4000, group_w = TRUE),
                 weight_formula = ~g)
  expect_equal(unname(coef(fit)[c("w[(Intercept)]", "w[g]")]), c(0.4, 1.2),
               tolerance = 0.25)
  expect_equal(sort(unique(round(fit$w, 6))),
               sort(unname(round(plogis(cumsum(fit$theta)), 6))))
})

test_that("centring a weight predictor changes only the intercept", {
  df <- sim_drm()
  f1 <- diagref(y ~ orig + dest, data = df, weight_formula = ~age)
  f2 <- diagref(y ~ orig + dest, data = df, weight_formula = ~I(age - 50))
  expect_equal(unname(coef(f1)["w[age]"]),
               unname(coef(f2)["w[I(age - 50)]"]), tolerance = 1e-6)
  expect_equal(as.numeric(logLik(f1)), as.numeric(logLik(f2)),
               tolerance = 1e-10)
  expect_equal(unname(confint(f1, "w[age]")),
               unname(confint(f2, "w[I(age - 50)]")), tolerance = 1e-5)
})

test_that("the analytic gradient matches finite differences", {
  df  <- sim_drm(group_w = TRUE)
  fit <- diagref(y ~ orig + dest + x, data = df, weight_formula = ~g + age)
  d   <- fit$data
  obj <- drm_objective(d$Z, d)
  th  <- fit$theta + c(0.3, -0.2, 0.01)
  h   <- 1e-6
  num <- vapply(seq_along(th), function(j) {
    e <- replace(numeric(3), j, h)
    (obj$fn(th + e) - obj$fn(th - e)) / (2 * h)
  }, numeric(1))
  expect_equal(unname(obj$gr(th)), num, tolerance = 1e-5)
})

# ---- Inference --------------------------------------------------------------

test_that("estimates and SEs match nls() with a single weight", {
  df  <- sim_drm(n = 600)
  fit <- diagref(y ~ orig + dest + x, data = df)
  cf  <- coef(fit)
  ref <- nls(
    y ~ b0 + bx * x +
      plogis(th) * c(s1, s2, s3, -s1 - s2 - s3)[orig] +
      (1 - plogis(th)) * c(s1, s2, s3, -s1 - s2 - s3)[dest],
    data = df,
    start = list(s1 = cf[["score[1]"]], s2 = cf[["score[2]"]],
                 s3 = cf[["score[3]"]], th = qlogis(cf[["w"]]),
                 b0 = cf[["(Intercept)"]], bx = cf[["x"]])
  )
  V  <- vcov(ref)
  w  <- plogis(coef(ref)[["th"]])
  se_ref <- c(sqrt(diag(V)[c("s1", "s2", "s3")]),
              sqrt(sum(V[1:3, 1:3])),            # score[4] = -(s1 + s2 + s3)
              sqrt(V["th", "th"]) * w * (1 - w), # delta method for w
              sqrt(diag(V)[c("b0", "bx")]))
  td <- tidy(fit, conf.int = FALSE)

  expect_equal(unname(cf[["w"]]), w, tolerance = 1e-6)
  expect_equal(td$std.error, unname(se_ref), tolerance = 1e-4)
  expect_equal(fit$df.residual, df.residual(ref))
})

test_that("estimates and SEs match nls() with a weight model", {
  df  <- sim_drm(n = 800, group_w = TRUE)
  fit <- diagref(y ~ orig + dest, data = df, weight_formula = ~g)
  cf  <- coef(fit)
  ref <- nls(
    y ~ b0 +
      plogis(a0 + a1 * g) * c(s1, s2, s3, -s1 - s2 - s3)[orig] +
      (1 - plogis(a0 + a1 * g)) * c(s1, s2, s3, -s1 - s2 - s3)[dest],
    data = df,
    start = list(s1 = cf[["score[1]"]], s2 = cf[["score[2]"]],
                 s3 = cf[["score[3]"]], a0 = cf[["w[(Intercept)]"]],
                 a1 = cf[["w[g]"]], b0 = cf[["(Intercept)"]])
  )
  se_ref <- sqrt(diag(vcov(ref)))
  se     <- sqrt(diag(vcov(fit)))

  expect_equal(unname(cf[c("w[(Intercept)]", "w[g]")]),
               unname(coef(ref)[c("a0", "a1")]), tolerance = 1e-5)
  expect_equal(unname(se[c("score[1]", "score[2]", "score[3]",
                           "w[(Intercept)]", "w[g]", "(Intercept)")]),
               unname(se_ref[c("s1", "s2", "s3", "a0", "a1", "b0")]),
               tolerance = 1e-4)
})

test_that("score SEs do not depend on which level is last", {
  # The last score is minus the sum of the others; its SE needs the
  # covariances of the contrast coefficients. Reordering the levels moves a
  # different category into that position, so SEs per label must agree.
  df   <- sim_drm()
  df_a <- transform(df, orig = factor(orig, levels = 1:4),
                        dest = factor(dest, levels = 1:4))
  df_b <- transform(df, orig = factor(orig, levels = c(4, 1, 2, 3)),
                        dest = factor(dest, levels = c(4, 1, 2, 3)))
  td_a <- tidy(diagref(y ~ orig + dest, data = df_a), conf.int = FALSE)
  td_b <- tidy(diagref(y ~ orig + dest, data = df_b), conf.int = FALSE)
  sa <- td_a[td_a$type == "score", ]
  sb <- td_b[match(sa$term, td_b$term), ]
  expect_equal(sa$estimate,  sb$estimate,  tolerance = 1e-5)
  expect_equal(sa$std.error, sb$std.error, tolerance = 1e-5)
})

test_that("profile CI endpoints sit at the likelihood-ratio cutoff", {
  crit <- qchisq(0.95, 1)

  # Single weight
  fit <- diagref(y ~ orig + dest, data = sim_drm())
  d   <- fit$data
  ci  <- confint(fit, "w")
  lr  <- function(w) {
    d$n * log(drm_linear_fit(qlogis(w), intercept_matrix(d$n), d)$rss /
                fit$rss)
  }
  expect_equal(lr(ci[1]), crit, tolerance = 1e-4)
  expect_equal(lr(ci[2]), crit, tolerance = 1e-4)
  expect_true(ci[1] < 0.6 && 0.6 < ci[2])

  # Weight model: profile out the intercept by brute force
  fit <- diagref(y ~ orig + dest, data = sim_drm(group_w = TRUE),
                 weight_formula = ~g)
  d   <- fit$data
  ci  <- confint(fit, "w[g]")
  lr  <- function(b1) {
    rss <- optimize(function(b0) drm_linear_fit(c(b0, b1), d$Z, d)$rss,
                    c(-10, 10), tol = 1e-10)$objective
    d$n * log(rss / fit$rss)
  }
  expect_equal(lr(ci[1]), crit, tolerance = 1e-3)
  expect_equal(lr(ci[2]), crit, tolerance = 1e-3)
})

test_that("profile CIs are right when a group's weight is near 1", {
  # Group g = 1 has w ~ 1, so the likelihood is flat as w[g] -> Inf. A naive
  # inner start pushes both groups' weights to the boundary and stalls,
  # giving a spurious finite endpoint.
  crit <- qchisq(0.95, 1)
  expect_warning(
    fit <- diagref(y ~ orig + dest,
                   data = sim_drm(group_w = TRUE, theta = c(0.4, 8)),
                   weight_formula = ~g),
    "boundary"
  )
  d    <- fit$data
  lr   <- function(b1) {
    rss <- optimize(function(b0) drm_linear_fit(c(b0, b1), d$Z, d)$rss,
                    c(-15, 15), tol = 1e-10)$objective
    d$n * log(rss / fit$rss)
  }
  ci <- confint(fit, "w[g]")
  expect_true(is.finite(ci[1]))
  expect_equal(lr(ci[1]), crit, tolerance = 1e-3)
  expect_equal(ci[2], Inf)
  expect_lt(lr(30), crit)
})

test_that("weight-model tests are LR tests that agree with the profile CIs", {
  fit <- diagref(y ~ orig + dest, data = sim_drm(group_w = TRUE),
                 weight_formula = ~g)
  td  <- tidy(fit)
  wg  <- td[td$term == "w[g]", ]
  # LR test by hand: refit with the weight constant across g
  ref <- diagref(y ~ orig + dest, data = sim_drm(group_w = TRUE))
  lr  <- 2 * (as.numeric(logLik(fit)) - as.numeric(logLik(ref)))
  expect_equal(wg$statistic, sign(wg$estimate) * sqrt(lr), tolerance = 1e-5)
  expect_equal(wg$p.value, pchisq(lr, 1, lower.tail = FALSE),
               tolerance = 1e-5)
  # 0 lies outside the 95% CI exactly when p < 0.05
  for (r in seq_len(nrow(td))[td$type == "weight"]) {
    expect_equal(td$p.value[r] < 0.05,
                 td$conf.low[r] > 0 || td$conf.high[r] < 0)
  }
})

test_that("tidy() agrees with coef(), vcov() and confint()", {
  fit <- diagref(y ~ orig + dest + x, data = sim_drm(), weight_formula = ~g)
  td  <- tidy(fit)
  expect_named(td, c("term", "type", "estimate", "std.error", "statistic",
                     "p.value", "conf.low", "conf.high"))
  expect_equal(td$term, names(coef(fit)))
  expect_equal(td$estimate, unname(coef(fit)))
  expect_equal(td$std.error, unname(sqrt(diag(vcov(fit)))))
  expect_equal(td$conf.low, unname(confint(fit)[, 1]))
  expect_equal(unique(td$type), c("score", "weight", "covariate"))
  expect_false("conf.low" %in% names(tidy(fit, conf.int = FALSE)))
})

test_that("a single weight is reported on the probability scale", {
  td <- tidy(diagref(y ~ orig + dest, data = sim_drm()))
  w  <- td[td$term == "w", ]
  expect_equal(w$type, "weight")
  expect_true(is.na(w$statistic))
  expect_true(w$conf.low > 0 && w$conf.high < 1)
})

# ---- glance / augment / logLik ---------------------------------------------

test_that("glance() returns model statistics", {
  fit <- diagref(y ~ orig + dest, data = sim_drm(n = 300))
  gl  <- glance(fit)
  expect_equal(nrow(gl), 1L)
  expect_true(gl$r.squared >= 0 && gl$r.squared <= 1)
  expect_equal(gl$nobs, 300L)
  expect_equal(gl$AIC, AIC(fit))
  expect_equal(gl$BIC, BIC(fit))
})

test_that("logLik() counts all free parameters", {
  fit <- diagref(y ~ orig + dest + x, data = sim_drm(), weight_formula = ~g)
  # 3 scores + 2 weight + intercept + x + sigma
  expect_equal(attr(logLik(fit), "df"), 3L + 2L + 2L + 1L)
  ref <- lm(y ~ factor(orig) + factor(dest) + x, data = sim_drm())
  expect_true(is.finite(AIC(fit, ref)$AIC[1]))
})

test_that("augment() returns the model frame plus fitted columns", {
  df  <- sim_drm(n = 100)
  fit <- diagref(y ~ orig + dest, data = df)
  aug <- augment(fit)
  expect_true(all(c(".fitted", ".resid", ".w", ".drm_term") %in% names(aug)))
  expect_equal(nrow(aug), 100L)
  expect_equal(aug$.fitted + aug$.resid, df$y)
  expect_equal(aug$.fitted, aug$.drm_term + coef(fit)[["(Intercept)"]])
})

# ---- Data handling ----------------------------------------------------------

test_that("NAs in a weight predictor drop those rows from the whole model", {
  df <- sim_drm(group_w = TRUE)
  df$g[c(3, 50, 400)] <- NA
  expect_no_warning(fit <- diagref(y ~ orig + dest, data = df,
                                   weight_formula = ~g))
  ref <- diagref(y ~ orig + dest, data = df[!is.na(df$g), ],
                 weight_formula = ~g)
  expect_equal(nobs(fit), nrow(df) - 3L)
  expect_length(fit$w, nobs(fit))
  expect_equal(coef(fit), coef(ref), tolerance = 1e-8)
  expect_equal(nrow(augment(fit)), nrow(df) - 3L)
})

test_that("subsetted data with non-sequential row names fits correctly", {
  df    <- sim_drm(group_w = TRUE)
  sub   <- df[df$y > 0.5, ]                    # row names like "2", "5", ...
  reset <- sub; rownames(reset) <- NULL
  chr   <- sub; rownames(chr) <- paste0("id", seq_len(nrow(chr)))

  expect_no_warning(fit <- diagref(y ~ orig + dest, data = sub,
                                   weight_formula = ~g))
  ref  <- diagref(y ~ orig + dest, data = reset, weight_formula = ~g)
  fchr <- diagref(y ~ orig + dest, data = chr,   weight_formula = ~g)

  expect_equal(nobs(fit), nrow(sub))
  expect_equal(coef(fit),  coef(ref), tolerance = 1e-8)
  expect_equal(coef(fchr), coef(ref), tolerance = 1e-8)
})

test_that("factor level order is preserved", {
  df  <- sim_drm(K = 3)
  lab <- c("low", "mid", "high")
  df$orig <- factor(lab[df$orig], levels = lab)
  df$dest <- factor(lab[df$dest], levels = lab)
  fit <- diagref(y ~ orig + dest, data = df)
  expect_equal(names(fit$scores), lab)
  expect_true(all(diff(fit$scores) > 0))   # simulated as increasing
})

test_that("numeric codes are ordered numerically, not alphabetically", {
  fit <- diagref(y ~ orig + dest, data = sim_drm(n = 3000, K = 12))
  expect_equal(names(fit$scores), as.character(1:12))
})

test_that("unobserved factor levels are dropped", {
  df <- sim_drm(K = 3)
  df$orig <- factor(df$orig, levels = 1:5)
  df$dest <- factor(df$dest, levels = 1:5)
  fit <- diagref(y ~ orig + dest, data = df)
  expect_equal(names(fit$scores), as.character(1:3))
  expect_false(anyNA(coef(fit)))
})

test_that("transformed and interacted covariates work", {
  df  <- sim_drm()
  fit <- diagref(y ~ orig + dest + log(age) + factor(g):x, data = df)
  expect_true(all(c("log(age)", "factor(g)0:x", "factor(g)1:x") %in%
                    names(coef(fit))))
})

test_that("bad input gives clear errors", {
  df <- sim_drm()
  expect_error(diagref(y ~ orig, data = df), "two right-hand-side terms")
  expect_error(diagref(y ~ orig + dest, data = df, weight_formula = g ~ age),
               "one-sided")
  expect_error(diagref(y ~ 0 + orig + dest, data = df), "intercept")
  df$f <- factor(df$y > 0)
  expect_error(diagref(f ~ orig + dest, data = df), "numeric")
  expect_error(diagref(y ~ orig + dest, data = df, weight_formula = ~g + I(2 * g)),
               "collinear")
})

# ---- Example data -----------------------------------------------------------

test_that("gss2024 is consistent with its documentation", {
  g <- gss2024
  expect_equal(dim(g), c(3230L, 19L))
  expect_false(anyNA(g$degree) || anyNA(g$pardeg))
  expect_equal(g$pardeg, pmax(g$madeg, g$padeg, na.rm = TRUE))
  for (v in c("degree", "madeg", "padeg", "pardeg")) {
    expect_true(all(g[[v]] %in% c(0:4, NA)))
  }
  expect_true(all(g$female %in% c(0L, 1L, NA)))
  expect_true(all(g$partyid %in% c(0:6, NA)))
  expect_true(all(g$polviews %in% c(1:7, NA)))
})

test_that("diagref runs on gss2024", {
  fit <- diagref(polviews ~ degree + pardeg + female + age, data = gss2024)
  expect_equal(names(fit$scores), as.character(0:4))
  expect_true(coef(fit)[["w"]] > 0 && coef(fit)[["w"]] < 1)
})
