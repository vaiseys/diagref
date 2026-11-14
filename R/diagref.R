#' Fit a diagonal reference model
#'
#' Fits Sobel's diagonal reference model (DRM) for an outcome that depends on
#' a pair of positions drawn from the same set of categories, such as origin
#' and destination class, or the education of two parents or partners:
#'
#' \deqn{y_i = w_i \mu_{a_i} + (1 - w_i) \mu_{b_i} + x_i'\gamma + \varepsilon_i}
#'
#' Here \eqn{a_i} and \eqn{b_i} are the categories of the two dyadic
#' variables, \eqn{\mu_k} is the score of category \eqn{k} (the scores sum to
#' zero, so the intercept is the mean score), and \eqn{w_i} is the weight, or
#' relative salience, of the first variable. The weight is either a single
#' constant or depends on predictors through
#' \eqn{\mathrm{logit}(w_i) = z_i'\theta}.
#'
#' @details
#' The model is fitted by maximum likelihood (equivalently, nonlinear least
#' squares). For fixed weights the model is linear, so the scores and
#' covariate coefficients are profiled out by OLS and only \eqn{\theta} is
#' optimised numerically, using an analytic gradient. Predictors in
#' `weight_formula` are standardised internally for the optimisation, and
#' estimates are reported on their original scale.
#'
#' Standard errors come from the joint asymptotic covariance of all
#' parameters, so uncertainty in the weight carries over to the scores and
#' covariates. Confidence intervals for weight parameters (from [tidy()] or
#' [confint()]) are profile-likelihood intervals; all others are Wald
#' intervals.
#'
#' Rows with a missing value in any variable of `formula` or `weight_formula`
#' are dropped.
#'
#' @param formula A formula `y ~ a + b + covariates`. The first two terms are
#'   the dyadic variables, and `w` is the weight of the first. They are
#'   treated as categorical and must share categories (numeric codes are
#'   fine); a warning is given if they have different numbers of categories.
#'   Further terms are covariates and enter linearly, as in [lm()].
#' @param data A data frame.
#' @param weight_formula An optional one-sided formula for the weight, such
#'   as `~ female` or `~ female + age`. The default, `NULL`, estimates a
#'   single weight `w`.
#'
#' @return An object of class `diagref_fit`, with methods for [tidy()],
#'   [glance()], [augment()], [coef()], [vcov()], [confint()], [logLik()]
#'   (and therefore [AIC()] and [BIC()]), [nobs()] and [print()].
#'
#' @references
#' Sobel, M. E. (1981). Diagonal mobility models: A substantively motivated
#' class of designs for the analysis of mobility effects. *American
#' Sociological Review*, 46(6), 893-906.
#'
#' @examples
#' # Political views by own and parents' education (GSS 2024)
#' fit <- diagref(polviews ~ degree + pardeg + female + age, data = gss2024)
#' fit
#' tidy(fit)
#'
#' # Let the weight of own education differ by gender
#' fit2 <- diagref(polviews ~ degree + pardeg + female + age, data = gss2024,
#'                 weight_formula = ~ female)
#' tidy(fit2)
#' AIC(fit, fit2)
#'
#' @export
diagref <- function(formula, data, weight_formula = NULL) {
  d   <- diagref_data(formula, data, weight_formula)
  fit <- drm_fit(d)
  structure(
    c(fit, list(call           = match.call(),
                formula        = formula,
                weight_formula = weight_formula,
                data           = d)),
    class = "diagref_fit"
  )
}
