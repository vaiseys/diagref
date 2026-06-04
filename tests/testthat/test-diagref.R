test_that("diagref fits a basic model and recovers w", {
  set.seed(42)
  n    <- 1000
  orig <- factor(sample(1:4, n, replace = TRUE))
  dest <- factor(sample(1:4, n, replace = TRUE))
  scores <- c(0, 0.5, 1.5, 3)
  w_true <- 0.7
  y <- w_true * scores[orig] + (1 - w_true) * scores[dest] + rnorm(n, sd = 0.5)
  df <- data.frame(y, orig, dest)

  fit <- diagref(y ~ orig + dest, data = df)

  # w should be recovered close to 0.7
  expect_equal(as.numeric(plogis(fit$theta_hat)), w_true, tolerance = 0.05)
})

test_that("scores sum to zero", {
  set.seed(1)
  df <- data.frame(
    y    = rnorm(200),
    orig = factor(sample(1:3, 200, replace = TRUE)),
    dest = factor(sample(1:3, 200, replace = TRUE))
  )
  fit <- diagref(y ~ orig + dest, data = df)
  expect_equal(sum(fit$scores), 0, tolerance = 1e-10)
})

test_that("tidy() returns a tibble with expected columns", {
  set.seed(2)
  df <- data.frame(
    y    = rnorm(200),
    orig = factor(sample(1:3, 200, replace = TRUE)),
    dest = factor(sample(1:3, 200, replace = TRUE))
  )
  fit <- diagref(y ~ orig + dest, data = df)
  td  <- tidy(fit)
  expect_true(inherits(td, "data.frame"))
  expect_true(all(c("term", "type", "estimate", "std.error") %in% names(td)))
  expect_true("w" %in% td$term)
})

test_that("glance() returns sensible model stats", {
  set.seed(3)
  df <- data.frame(
    y    = rnorm(300),
    orig = factor(sample(1:4, 300, replace = TRUE)),
    dest = factor(sample(1:4, 300, replace = TRUE))
  )
  fit <- diagref(y ~ orig + dest, data = df)
  gl  <- glance(fit)
  expect_equal(nrow(gl), 1)
  expect_true(gl$r.squared >= 0 && gl$r.squared <= 1)
  expect_equal(gl$nobs, 300L)
})

test_that("augment() returns original data plus fitted columns", {
  set.seed(4)
  df <- data.frame(
    y    = rnorm(100),
    orig = factor(sample(1:3, 100, replace = TRUE)),
    dest = factor(sample(1:3, 100, replace = TRUE))
  )
  fit <- diagref(y ~ orig + dest, data = df)
  aug <- augment(fit)
  expect_true(all(c(".fitted", ".resid", ".w", ".drm_term") %in% names(aug)))
  expect_equal(nrow(aug), 100L)
})

test_that("covariates are included correctly", {
  set.seed(5)
  n    <- 500
  orig <- factor(sample(1:3, n, replace = TRUE))
  dest <- factor(sample(1:3, n, replace = TRUE))
  age  <- rnorm(n)
  y    <- 0.5 * as.numeric(orig) + 0.5 * as.numeric(dest) + 1.5 * age + rnorm(n)
  df   <- data.frame(y, orig, dest, age)

  fit <- diagref(y ~ orig + dest + age, data = df)
  td  <- tidy(fit)

  expect_true("age" %in% td$term)
  age_est <- td$estimate[td$term == "age"]
  expect_equal(unname(age_est), 1.5, tolerance = 0.1)
})

test_that("profile CI contains the true w", {
  set.seed(6)
  n    <- 800
  orig <- factor(sample(1:4, n, replace = TRUE))
  dest <- factor(sample(1:4, n, replace = TRUE))
  scores <- c(0, 1, 2, 4)
  w_true <- 0.6
  y <- w_true * scores[orig] + (1 - w_true) * scores[dest] + rnorm(n)
  df <- data.frame(y, orig, dest)

  fit <- diagref(y ~ orig + dest, data = df)
  ci  <- fit$ci$ci_w
  expect_true(ci[1] <= w_true && w_true <= ci[2])
})
