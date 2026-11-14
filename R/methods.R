#' Tidy, glance, augment and confint methods for diagonal reference models
#'
#' `tidy()` returns one row per parameter: the category scores, the weight,
#' and the covariate coefficients. With a single weight it is reported as
#' `w` on the probability scale. With a weight model it is reported as
#' logit-scale coefficients `w[...]`, so that `w[(Intercept)]` is logit(w)
#' when all weight predictors are zero.
#'
#' Standard errors come from the joint covariance of all parameters. For the
#' weight parameters, confidence intervals are profile-likelihood intervals,
#' which respect the 0-1 range of `w` and can be asymmetric; an infinite
#' endpoint means that the likelihood is flat in that direction (a weight
#' heading to 0 or 1). Likewise, the tests of weight-model coefficients are
#' likelihood-ratio tests, with `statistic` the signed square root of the LR
#' statistic, so they agree with the intervals. Near the boundary the Wald
#' standard error of a weight coefficient becomes very large and should not
#' be used for inference. A single weight `w` has no test, because `w = 0`
#' is not usually of interest. All other parameters have Wald intervals and
#' t tests.
#'
#' @param x,object A `diagref_fit` object.
#' @param conf.int Whether to include confidence intervals.
#' @param conf.level,level Confidence level.
#' @param parm Parameters to compute intervals for, by name or position
#'   (default: all).
#' @param ... Unused.
#'
#' @return `tidy()`: a tibble with columns `term`, `type` (`"score"`,
#'   `"weight"` or `"covariate"`), `estimate`, `std.error`, `statistic`,
#'   `p.value` and, optionally, `conf.low` and `conf.high`.
#'
#'   `glance()`: a one-row tibble of model statistics.
#'
#'   `augment()`: the model frame with columns `.fitted`, `.resid`, `.w`
#'   (each observation's weight) and `.drm_term` (its weighted score,
#'   \eqn{w_i \mu_{a_i} + (1 - w_i) \mu_{b_i}}).
#'
#'   `confint()`: a matrix with one row per parameter.
#'
#' @name diagref_tidiers
NULL

#' @rdname diagref_tidiers
#' @export
tidy.diagref_fit <- function(x, conf.int = TRUE, conf.level = 0.95, ...) {
  est  <- x$coefficients
  se   <- sqrt(diag(x$vcov))
  stat <- est / se
  pval <- 2 * pt(abs(stat), x$df.residual, lower.tail = FALSE)

  is_w <- x$term_type == "weight"
  if (x$scalar_weight) {
    stat[is_w] <- pval[is_w] <- NA_real_
  } else {
    # Likelihood-ratio tests, consistent with the profile CIs. Wald tests are
    # unreliable here: near w = 0 or 1 the SE explodes (Hauck-Donner effect).
    for (j in seq_len(sum(is_w))) {
      a  <- replace(numeric(sum(is_w)), j, 1)
      lr <- profile_lr_test(x, a)
      stat[is_w][j] <- lr[["statistic"]]
      pval[is_w][j] <- lr[["p.value"]]
    }
  }

  out <- tibble::tibble(
    term      = names(est),
    type      = x$term_type,
    estimate  = unname(est),
    std.error = unname(se),
    statistic = unname(stat),
    p.value   = unname(pval)
  )
  if (conf.int) {
    ci <- confint(x, level = conf.level)
    out$conf.low  <- unname(ci[, 1])
    out$conf.high <- unname(ci[, 2])
  }
  out
}

#' @rdname diagref_tidiers
#' @export
glance.diagref_fit <- function(x, ...) {
  y  <- x$data$y
  ll <- logLik(x)
  tibble::tibble(
    r.squared   = 1 - x$rss / sum((y - mean(y))^2),
    sigma       = x$sigma,
    logLik      = as.numeric(ll),
    AIC         = AIC(ll),
    BIC         = BIC(ll),
    deviance    = x$rss,
    df.residual = x$df.residual,
    nobs        = x$data$n
  )
}

#' @rdname diagref_tidiers
#' @export
augment.diagref_fit <- function(x, ...) {
  d   <- x$data
  s   <- unname(x$scores)
  out <- tibble::as_tibble(d$mf)
  out$.fitted   <- x$fitted
  out$.resid    <- x$residuals
  out$.w        <- x$w
  out$.drm_term <- x$w * s[d$a] + (1 - x$w) * s[d$b]
  out
}

#' @rdname diagref_tidiers
#' @export
confint.diagref_fit <- function(object, parm, level = 0.95, ...) {
  est <- object$coefficients
  if (missing(parm)) parm <- names(est)
  if (is.numeric(parm)) parm <- names(est)[parm]

  se <- sqrt(diag(object$vcov))
  tq <- qt((1 + level) / 2, object$df.residual)
  ci <- cbind(est - tq * se, est + tq * se)[parm, , drop = FALSE]

  w_names <- names(est)[object$term_type == "weight"]
  for (j in which(w_names %in% parm)) {
    a <- replace(numeric(length(w_names)), j, 1)
    bounds <- profile_ci_weight(object, a, level = level)
    ci[w_names[j], ] <- if (object$scalar_weight) plogis(bounds) else bounds
  }

  pct <- format(100 * c(1 - level, 1 + level) / 2, trim = TRUE,
                scientific = FALSE, digits = 3)
  colnames(ci) <- paste(pct, "%")
  ci
}

#' @export
coef.diagref_fit <- function(object, ...) object$coefficients

#' @export
vcov.diagref_fit <- function(object, ...) object$vcov

#' @export
nobs.diagref_fit <- function(object, ...) object$data$n

#' @export
logLik.diagref_fit <- function(object, ...) {
  n <- object$data$n
  # ML estimate of sigma^2 uses n, not the residual df
  ll <- -n / 2 * (log(2 * pi * object$rss / n) + 1)
  # Free parameters: mean-model parameters plus sigma^2
  df <- n - object$df.residual + 1L
  structure(ll, df = df, nobs = n, class = "logLik")
}

#' @export
print.diagref_fit <- function(x, digits = 4, ...) {
  fmt <- function(v) formatC(v, digits = digits, format = "f")

  cat("Diagonal reference model\n\n")
  cat("Formula:", paste(deparse(x$formula), collapse = " "), "\n")

  if (x$scalar_weight) {
    ci <- confint(x, "w")
    cat(sprintf("Weight:  w = %s  (95%% profile CI %s to %s)\n",
                fmt(plogis(x$theta)), fmt(ci[1]), fmt(ci[2])))
  } else {
    th    <- x$theta
    parts <- ifelse(names(th) == "(Intercept)", fmt(abs(th)),
                    paste(fmt(abs(th)), "*", names(th)))
    signs <- ifelse(th < 0, "- ", "+ ")
    eq    <- paste0(signs, parts, collapse = " ")
    eq    <- sub("^\\+ ", "", eq)
    cat("Weight:  logit(w) =", eq, "\n")
  }

  cat("\nScores (", x$data$dyad[1], " / ", x$data$dyad[2], "):\n", sep = "")
  print(round(x$scores, digits))

  cat("\nCovariates:\n")
  print(round(x$cov_coef, digits))

  cat(sprintf("\nn = %d, residual df = %d, sigma = %s\n",
              x$data$n, x$df.residual, fmt(x$sigma)))
  invisible(x)
}
