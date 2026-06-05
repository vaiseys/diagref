#' Fit a diagonal reference model
#'
#' @param formula A formula of the form `y ~ var1 + var2 + covariates`. The
#'   first two RHS terms are the dyadic variables (e.g., origin/destination
#'   class, partner1/partner2). Both are coerced to factors.
#' @param data A data frame.
#' @param weight_formula Optional one-sided formula specifying what predicts
#'   the salience weight `w`. Default (`NULL`) estimates a single scalar `w`.
#'   Example: `~gender` estimates a separate `w` per gender group.
#' @param level Confidence level for profile likelihood CIs. Default `0.95`.
#' @param n_grid Number of grid points for profile CI search. Default `200`.
#'
#' @return An object of class `diagref_fit`. Use [tidy()], [glance()], and
#'   [augment()] to extract results.
#'
#' @examples
#' # Simulated data
#' set.seed(1)
#' n <- 500
#' orig <- factor(sample(1:4, n, replace = TRUE))
#' dest <- factor(sample(1:4, n, replace = TRUE))
#' y <- 0.6 * as.numeric(orig) + 0.4 * as.numeric(dest) + rnorm(n)
#' df <- data.frame(y, orig, dest)
#' fit <- diagref(y ~ orig + dest, data = df)
#' tidy(fit)
#'
#' @export
diagref <- function(formula, data, weight_formula = NULL, level = 0.95,
                    n_grid = 200) {
  parsed      <- parse_diagref_formula(formula, data)
  data_subset <- data[parsed$kept_rows, , drop = FALSE]
  Z_w         <- build_weight_matrix(weight_formula, data_subset)

  # Keep original weight predictor columns for per-group display in tidy()
  weight_data <- if (!is.null(weight_formula)) {
    wvars <- intersect(all.vars(weight_formula), names(data_subset))
    data_subset[, wvars, drop = FALSE]
  } else {
    NULL
  }

  raw <- drm_fit(parsed, Z_w)
  ci  <- compute_profile_ci(raw,
                             weight_data    = weight_data,
                             weight_formula = weight_formula,
                             level          = level,
                             n_grid         = n_grid)

  structure(
    list(
      formula        = formula,
      weight_formula = weight_formula,
      weight_data    = weight_data,
      theta_hat      = raw$theta_hat,
      vcov_theta     = raw$vcov_theta,
      w_hat          = raw$w_hat,
      scores         = raw$scores,
      cov_coef       = raw$cov_coef,
      inner          = raw$inner,
      rss            = raw$rss,
      ci             = ci,
      level          = level,
      parsed         = parsed,
      Z_w            = Z_w,
      K              = raw$K
    ),
    class = "diagref_fit"
  )
}
