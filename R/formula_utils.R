# Parse the diagref formula into its components.
# Returns a list with: response, var1_name, var2_name, covar_terms, data
parse_diagref_formula <- function(formula, data) {
  tt <- terms(formula, data = data)
  vars <- attr(tt, "term.labels")

  if (length(vars) < 2) {
    stop("Formula must have at least two RHS variables (the dyadic pair).", call. = FALSE)
  }

  var1_name <- vars[1]
  var2_name <- vars[2]
  covar_names <- vars[-(1:2)]

  response_name <- as.character(attr(tt, "variables")[[2]])

  mf <- model.frame(formula, data = data, na.action = na.omit)
  kept_rows <- as.integer(rownames(mf))
  y <- as.numeric(model.response(mf))

  var1 <- mf[[var1_name]]
  var2 <- mf[[var2_name]]

  # Coerce dyadic vars to factor if not already
  if (!is.factor(var1)) var1 <- factor(var1)
  if (!is.factor(var2)) var2 <- factor(var2)

  # Enforce shared levels across both dyadic variables
  all_levels <- sort(union(levels(var1), levels(var2)))
  var1 <- factor(var1, levels = all_levels)
  var2 <- factor(var2, levels = all_levels)

  # Covariate design matrix (with intercept)
  if (length(covar_names) > 0) {
    covar_formula <- reformulate(covar_names, intercept = TRUE)
    X_cov <- model.matrix(covar_formula, data = mf)
  } else {
    X_cov <- matrix(1, nrow = length(y), ncol = 1,
                    dimnames = list(NULL, "(Intercept)"))
  }

  list(
    y         = y,
    var1      = var1,
    var2      = var2,
    all_levels = all_levels,
    X_cov     = X_cov,
    n         = length(y),
    response_name = response_name,
    var1_name = var1_name,
    var2_name = var2_name,
    covar_names = covar_names,
    mf        = mf,
    kept_rows = kept_rows
  )
}

# Build the weight model matrix from weight_formula and data.
# data_subset is the original data filtered to the rows used after na.omit.
# Returns an n x p matrix; default (NULL) is a single intercept column.
#
# Standard (with-intercept) model.matrix is used so that continuous and
# categorical predictors both work. The intercept represents the baseline
# weight; tidy() recovers per-group predicted weights via delta method.
build_weight_matrix <- function(weight_formula, data_subset) {
  n <- nrow(data_subset)
  if (is.null(weight_formula)) {
    return(matrix(1, nrow = n, ncol = 1,
                  dimnames = list(NULL, "(Intercept)")))
  }
  model.matrix(weight_formula, data = data_subset)
}

# Build a small representative data frame for display in print()/tidy().
# Categorical columns (factor, character, or integer ≤ 10 levels) are shown
# at all unique combinations; continuous columns are fixed at their sample
# mean. Returns NULL when all predictors are continuous (no natural grouping).
make_representative_data <- function(weight_data) {
  if (is.null(weight_data)) return(NULL)

  is_cat <- function(col) {
    is.factor(col) || is.character(col) ||
      (is.integer(col) && length(unique(col)) <= 10L)
  }
  cat_flag <- vapply(weight_data, is_cat, logical(1))

  if (!any(cat_flag)) return(NULL)

  cat_data <- weight_data[, cat_flag, drop = FALSE]
  u_cat    <- unique(cat_data)
  u_cat    <- u_cat[do.call(order, as.list(u_cat)), , drop = FALSE]
  rownames(u_cat) <- NULL

  if (all(cat_flag)) return(u_cat)

  cont_data <- weight_data[, !cat_flag, drop = FALSE]
  cont_rep  <- as.data.frame(lapply(cont_data, function(col) rep(mean(col), nrow(u_cat))))
  cbind(u_cat, cont_rep)[, names(weight_data), drop = FALSE]
}
