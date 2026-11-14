# devtools::load_all(".")
library(tidyverse)
library(broom)

gss2024 <- readRDS("gss2024.rds") |> 
  haven::zap_labels()

d <- gss2024 |> 
  select(y = fefam, madeg, padeg, degree, sex, age) |> 
  mutate(pardeg = pmax(madeg, padeg, na.rm = TRUE),
         female = if_else(sex == 2, 1L, 0L),
         age45 = age-45) |> 
  select(y, pardeg, degree, female, age45) |> 
  drop_na()

fit <- diagref(y ~ degree + pardeg + female + age45,
               data = d)

tidy(fit)

fit2 <- lm(y ~ factor(degree) + factor(pardeg) + female + age45,
           data = d)

AIC(fit, fit2)
