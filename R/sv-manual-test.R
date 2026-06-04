# devtools::load_all(".")
library(tidyverse)

gss2024 <- readRDS("gss2024.rds") |> 
  haven::zap_labels()

d <- gss2024 |> 
  select(madeg, padeg, degree, polviews, sex) |> 
  mutate(pardeg = pmax(madeg, padeg, na.rm = TRUE),
         female = if_else(sex == 2, 1L, 0L)) |> 
  select(pardeg, degree, polviews, female) |> 
  drop_na()

fit <- diagref(polviews ~ degree + pardeg + female,
               data = d,
               weight_formula = ~female)

tidy(fit)
