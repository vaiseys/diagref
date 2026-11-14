# Future options

Ideas for later, not yet planned.

## Survey weights

`diagref()` currently treats all observations equally. (`weight_formula` is
the model's own salience weight `w`, not a survey weight.) Supporting survey
weights, e.g. `wtssps` in `gss2024`, would involve:

- **Estimates:** weighted least squares in the inner OLS step. The variable
  projection and analytic gradient carry over unchanged. Small change.
- **Standard errors:** a sandwich (robust) variance instead of
  sigma^2 (J'J)^-1, ideally design-based using the GSS clusters and strata
  (`vpsu`, `vstrat` in the full file), possibly via the `survey` package.
- **Inference:** with weights the likelihood is a pseudo-likelihood, so
  profile CIs and LR tests at the plain chi-square cutoff are not valid. Use
  Wald intervals from the sandwich variance, or Rao-Scott-adjusted LR tests.
  `logLik()`, `AIC()` and `BIC()` should be disabled or warn.

Likely interface: a `weights` argument, plus optional `cluster` / `strata`.

## Other candidates

- Binary outcomes: a logistic DRM (`family =` argument), for items like
  `abany`, `cappun` or `gunlaw`.
- `predict()` and `augment(newdata = )`.
- `anova()` likelihood-ratio tests between nested weight models
  (e.g. `~1` vs `~female`).
- Fill in the author field in `DESCRIPTION` (still the placeholder).
