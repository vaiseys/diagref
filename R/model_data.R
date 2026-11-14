# Assemble everything the fitter needs from the user's formulas and data.
#
# The DRM mean is
#   mu_i = w_i * score[a_i] + (1 - w_i) * score[b_i] + x_i' gamma,
#   logit(w_i) = z_i' theta,
# with scores coded sum-to-zero: score = C %*% delta, C = contr.sum(K).
# For the DRM block we precompute C1 = C[a, ], C2 = C[b, ] and D = C1 - C2,
# so that the design for fixed w is  C2 + w * D  (n x (K-1)).
diagref_data <- function(formula, data, weight_formula = NULL) {
  if (!is.null(weight_formula) && length(weight_formula) != 2L) {
    stop("`weight_formula` must be one-sided, e.g. `~ female`.", call. = FALSE)
  }

  data <- complete_rows(data, formula, weight_formula)
  mf   <- model.frame(formula, data = data, drop.unused.levels = TRUE)
  tt   <- terms(mf)
  labs <- attr(tt, "term.labels")

  if (attr(tt, "response") == 0L) {
    stop("`formula` needs a response, e.g. `y ~ a + b`.", call. = FALSE)
  }
  if (length(labs) < 2L) {
    stop("`formula` must have at least two right-hand-side terms ",
         "(the dyadic pair).", call. = FALSE)
  }
  if (attr(tt, "intercept") == 0L) {
    stop("diagref models need an intercept; remove `0 +` or `- 1` ",
         "from `formula`.", call. = FALSE)
  }

  y <- model.response(mf)
  if (!is.numeric(y) && !is.logical(y)) {
    stop("The response must be numeric.", call. = FALSE)
  }
  y <- as.numeric(y)
  n <- length(y)

  a <- mf[[labs[1]]]
  b <- mf[[labs[2]]]
  if (is.null(a) || is.null(b) || !is.null(dim(a)) || !is.null(dim(b))) {
    stop("The first two terms of `formula` must be single variables ",
         "(the dyadic pair).", call. = FALSE)
  }

  lev <- shared_levels(a, b)
  K   <- length(lev)
  if (K < 2L) {
    stop("The dyadic variables must have at least two categories.",
         call. = FALSE)
  }
  C  <- contr.sum(K)
  C1 <- C[match(as.character(a), lev), , drop = FALSE]
  C2 <- C[match(as.character(b), lev), , drop = FALSE]

  X_cov <- if (length(labs) > 2L) {
    model.matrix(drop.terms(tt, 1:2, keep.response = FALSE), mf)
  } else {
    intercept_matrix(n)
  }
  if (qr(X_cov)$rank < ncol(X_cov)) {
    stop("The covariates are collinear.", call. = FALSE)
  }

  Z <- if (is.null(weight_formula)) {
    intercept_matrix(n)
  } else {
    model.matrix(weight_formula, data = data)
  }
  if (nrow(Z) != n) {
    stop("Internal error: weight model has ", nrow(Z), " rows but outcome ",
         "model has ", n, ".", call. = FALSE)
  }
  if (qr(Z)$rank < ncol(Z)) {
    stop("The predictors in `weight_formula` are collinear.", call. = FALSE)
  }

  list(
    y        = y,
    n        = n,
    K        = K,
    levels   = lev,
    a        = match(as.character(a), lev),
    b        = match(as.character(b), lev),
    C        = C,
    C2       = C2,
    D        = C1 - C2,
    X_cov    = X_cov,
    Z        = Z,
    dyad     = labs[1:2],
    response = deparse(attr(tt, "variables")[[2L]]),
    mf       = mf
  )
}

intercept_matrix <- function(n) {
  matrix(1, nrow = n, ncol = 1L, dimnames = list(NULL, "(Intercept)"))
}

# Subset `data` to rows that are complete for every term in `formula` and
# (if given) `weight_formula`. Rows are selected by position, so this works
# for subsetted data frames with non-sequential or character row names.
# Evaluating via model.frame() means transformed terms that yield NA/NaN
# (e.g. log() of a non-positive value) are dropped too.
complete_rows <- function(data, formula, weight_formula = NULL) {
  ok <- complete.cases(model.frame(formula, data = data, na.action = na.pass))
  if (!is.null(weight_formula)) {
    ok <- ok & complete.cases(model.frame(weight_formula, data = data,
                                          na.action = na.pass))
  }
  if (!any(ok)) {
    stop("No complete observations after removing missing values.",
         call. = FALSE)
  }
  data[ok, , drop = FALSE]
}

# Common level set for the two dyadic variables, in a meaningful order:
# - both numeric: sorted numerically (so 1, 2, ..., 10, not 1, 10, 2, ...)
# - otherwise: the first variable's level order (factor order is respected,
#   e.g. an ordered low/mid/high factor), followed by any levels that appear
#   only in the second. Plain character vectors are sorted, as factor() does.
# Levels not observed in either variable are dropped (as lm() does); keeping
# them would make the design rank-deficient.
shared_levels <- function(a, b) {
  if (is.numeric(a) && is.numeric(b)) {
    return(as.character(sort(unique(c(a, b)))))
  }
  lev_a <- levels(if (is.factor(a)) a else factor(a))
  lev_b <- levels(if (is.factor(b)) b else factor(b))
  lev   <- union(lev_a, lev_b)
  if (!is.factor(a) && !is.factor(b)) lev <- sort(lev)
  lev[lev %in% c(as.character(a), as.character(b))]
}

# Well-conditioned reparameterisation of the weight model for optimisation.
# Returns A such that Z %*% A has centred, unit-variance non-intercept
# columns; then Z theta = (Z A) theta_s, i.e. theta = A theta_s. Without an
# intercept column, centring would change the model, so A is the identity.
weight_reparam <- function(Z) {
  p   <- ncol(Z)
  A   <- diag(p)
  int <- which(colnames(Z) == "(Intercept)")
  if (length(int) == 1L && p > 1L) {
    j <- setdiff(seq_len(p), int)
    m <- colMeans(Z[, j, drop = FALSE])
    s <- apply(Z[, j, drop = FALSE], 2L, stats::sd)
    A[cbind(j, j)] <- 1 / s
    A[int, j]      <- -m / s
  }
  A
}

is_scalar_weight <- function(Z) {
  ncol(Z) == 1L && colnames(Z) == "(Intercept)"
}
