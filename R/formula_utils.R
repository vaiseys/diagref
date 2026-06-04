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
# The intercept is suppressed by default so that each column maps directly
# to one group/level: plogis(theta_j) = w for that group. Users who want a
# proper intercept for continuous predictors should write ~1 + continuous_var.
build_weight_matrix <- function(weight_formula, data_subset) {
  n <- nrow(data_subset)
  if (is.null(weight_formula)) {
    return(matrix(1, nrow = n, ncol = 1,
                  dimnames = list(NULL, "w")))
  }
  # Remove intercept so factor levels get their own columns
  no_int_formula <- update(weight_formula, ~ . - 1)
  model.matrix(no_int_formula, data = data_subset)
}
