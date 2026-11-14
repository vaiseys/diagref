
<!-- README.md is generated from README.Rmd. Please edit that file, then
     run devtools::build_readme(). -->

# diagref

diagref fits **diagonal reference models** (DRMs; Sobel 1981) in R with
a tidy, broom-style interface.

A DRM is for outcomes that depend on a pair of positions drawn from the
same set of categories: your own education and your parents’, your class
of origin and destination, your and your partner’s religion. Each
category $k$ has a score $\mu_k$, and the outcome is a weighted average
of the scores of the two positions,

$$y_i = w\,\mu_{a_i} + (1 - w)\,\mu_{b_i} + x_i'\gamma + \varepsilon_i .$$

The weight $w$ is the relative salience of the first position: with
$w = 0.8$, own education counts four times as much as parents’. The
weight can also depend on covariates,
$\mathrm{logit}(w_i) = z_i'\theta$.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("vaiseys/diagref")
```

## Example

`gss2024`, an extract of the 2024 General Social Survey, comes with the
package. Do gender-role attitudes follow people’s own education or their
parents’? (`fefam` runs from 1 = traditional to 4 = egalitarian.)

``` r
library(diagref)

fit <- diagref(fefam ~ degree + pardeg + female + age, data = gss2024)
fit
#> Diagonal reference model
#> 
#> Formula: fefam ~ degree + pardeg + female + age 
#> Weight:  w = 0.8169  (95% profile CI 0.6902 to 0.9539)
#> 
#> Scores (degree / pardeg):
#>       0       1       2       3       4 
#> -0.4762 -0.1524 -0.1140  0.2894  0.4532 
#> 
#> Covariates:
#> (Intercept)      female         age 
#>      2.8653      0.1841     -0.0012 
#> 
#> n = 2047, residual df = 2039, sigma = 0.8551
```

Own education carries most of the weight. `tidy()` gives every
parameter, with standard errors and confidence intervals:

``` r
tidy(fit)
#> # A tibble: 9 × 8
#>   term        type      estimate std.error statistic    p.value conf.low conf.high
#>   <chr>       <chr>        <dbl>     <dbl>     <dbl>      <dbl>    <dbl>     <dbl>
#> 1 score[0]    score     -0.476     0.0594      -8.01  1.91e- 15 -0.593   -0.360   
#> 2 score[1]    score     -0.152     0.0368      -4.14  3.59e-  5 -0.225   -0.0802  
#> 3 score[2]    score     -0.114     0.0650      -1.76  7.94e-  2 -0.241    0.0134  
#> 4 score[3]    score      0.289     0.0447       6.48  1.17e- 10  0.202    0.377   
#> 5 score[4]    score      0.453     0.0536       8.45  5.39e- 17  0.348    0.558   
#> 6 w           weight     0.817     0.0660      NA    NA          0.690    0.954   
#> 7 (Intercept) covariate  2.87      0.0633      45.3   3.63e-310  2.74     2.99    
#> 8 female      covariate  0.184     0.0383       4.81  1.64e-  6  0.109    0.259   
#> 9 age         covariate -0.00124   0.00109     -1.13  2.57e-  1 -0.00338  0.000903
```

Does the weight differ for women and men? `weight_formula` models it:

``` r
fit2 <- diagref(fefam ~ degree + pardeg + female + age, data = gss2024,
                weight_formula = ~ female)
subset(tidy(fit2), type == "weight")
#> # A tibble: 2 × 8
#>   term           type   estimate std.error statistic p.value conf.low conf.high
#>   <chr>          <chr>     <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#> 1 w[(Intercept)] weight    0.657     0.417      1.65  0.0999   -0.122      1.75
#> 2 w[female]      weight    2.00      1.49       2.20  0.0280    0.193    Inf

# Weights for men and women
th <- coef(fit2)[c("w[(Intercept)]", "w[female]")]
plogis(c(men = th[[1]], women = sum(th)))
#>       men     women 
#> 0.6586191 0.9341890

AIC(fit, fit2)
#>      df      AIC
#> fit   9 5178.132
#> fit2 10 5175.306
```

Women’s attitudes track their own education more closely. The interval
for `w[female]` is open-ended because the likelihood stays flat as
women’s weight approaches 1.

## Usage

- The **first two terms** on the right-hand side are the pair, and `w`
  is the weight of the first. Both are treated as categorical, so
  numeric codes are fine. A warning is given if they have different
  numbers of categories, which usually means a covariate was written
  first.
- **Further terms** are covariates and enter linearly, as in `lm()`.
- `weight_formula` is a one-sided formula for what predicts `w`. Leave
  it out to estimate a single weight.
- Rows with missing values in any model variable are dropped.

## Methods and inference

- `tidy()`, `glance()` and `augment()`, plus `coef()`, `vcov()`,
  `confint()`, `nobs()`, `logLik()`, `AIC()` and `BIC()`.
- Estimates are maximum likelihood. Standard errors come from the joint
  covariance of all parameters, so uncertainty in the weight carries
  over to the scores and covariates.
- Weight parameters get profile-likelihood confidence intervals and
  likelihood-ratio tests. These stay valid when a weight is near 0 or 1,
  where Wald standard errors break down.

## Status

diagref is in early development, and its interface may change. Planned
options, including survey weights and binary outcomes, are listed in
[TODO.md](TODO.md).

For DRMs within generalized nonlinear models, see also `Dref()` in the
[gnm](https://cran.r-project.org/package=gnm) package.

## Reference

Sobel, M. E. (1981). Diagonal mobility models: A substantively motivated
class of designs for the analysis of mobility effects. *American
Sociological Review*, 46(6), 893–906.
