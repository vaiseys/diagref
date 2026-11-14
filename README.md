
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
parents’? (`fefam` runs from 1 = traditional to 4 = egalitarian.) We
keep respondents aged 25 and over, since many younger ones have not
finished their education.

``` r
library(diagref)

gss25 <- subset(gss2024, age >= 25)
fit <- diagref(fefam ~ degree + pardeg + female + age, data = gss25)
fit
#> Diagonal reference model
#> 
#> Formula: fefam ~ degree + pardeg + female + age 
#> Weight:  w = 0.8803  (95% profile CI 0.7416 to 1.0000)
#> 
#> Scores (degree / pardeg):
#>       0       1       2       3       4 
#> -0.4629 -0.1464 -0.0979  0.2763  0.4310 
#> 
#> Covariates:
#> (Intercept)      female         age 
#>      2.8424      0.1759     -0.0009 
#> 
#> n = 1899, residual df = 1891, sigma = 0.8546
```

Own education carries most of the weight. `tidy()` gives every
parameter, with standard errors and confidence intervals:

``` r
tidy(fit)
#> # A tibble: 9 × 8
#>   term        type       estimate std.error statistic    p.value conf.low conf.high
#>   <chr>       <chr>         <dbl>     <dbl>     <dbl>      <dbl>    <dbl>     <dbl>
#> 1 score[0]    score     -0.463      0.0622     -7.44   1.47e- 13 -0.585    -0.341  
#> 2 score[1]    score     -0.146      0.0367     -3.99   6.78e-  5 -0.218    -0.0745 
#> 3 score[2]    score     -0.0979     0.0635     -1.54   1.23e-  1 -0.222     0.0266 
#> 4 score[3]    score      0.276      0.0439      6.30   3.72e- 10  0.190     0.362  
#> 5 score[4]    score      0.431      0.0531      8.12   8.39e- 16  0.327     0.535  
#> 6 w           weight     0.880      0.0735     NA     NA          0.742     1      
#> 7 (Intercept) covariate  2.84       0.0717     39.7    6.36e-251  2.70      2.98   
#> 8 female      covariate  0.176      0.0397      4.43   1.00e-  5  0.0980    0.254  
#> 9 age         covariate -0.000893   0.00122    -0.731  4.65e-  1 -0.00329   0.00150
```

Does the weight differ for women and men? `weight_formula` models it:

``` r
fit2 <- diagref(fefam ~ degree + pardeg + female + age, data = gss25,
                weight_formula = ~ female)
subset(tidy(fit2), type == "weight")
#> # A tibble: 2 × 8
#>   term           type   estimate std.error statistic p.value conf.low conf.high
#>   <chr>          <chr>     <dbl>     <dbl>     <dbl>   <dbl>    <dbl>     <dbl>
#> 1 w[(Intercept)] weight    0.936     0.509      2.05  0.0400   0.0399      2.63
#> 2 w[female]      weight    4.43     21.5        2.20  0.0278   0.0715    Inf

# Weights for men and women
th <- coef(fit2)[c("w[(Intercept)]", "w[female]")]
plogis(c(men = th[[1]], women = sum(th)))
#>       men     women 
#> 0.7183657 0.9953673

AIC(fit, fit2)
#>      df      AIC
#> fit   9 4802.345
#> fit2 10 4800.283
```

Women’s attitudes track their own education more closely. The interval
for `w[female]` is open-ended because the likelihood stays flat as
women’s weight approaches 1. For the same reason its Wald standard error
is huge and would suggest no difference at all, which is why weight
coefficients are tested with likelihood-ratio tests instead.

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
