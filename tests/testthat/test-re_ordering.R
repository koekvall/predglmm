# Tests that RE extraction is invariant to ordering of RE terms

test_that("pred_glmer invariant to RE term order", {
  skip_on_cran()
  set.seed(1)
  n <- 300
  g1 <- factor(sample(1:15, n, replace = TRUE))
  g2 <- factor(sample(1:10, n, replace = TRUE))
  b1 <- rnorm(15, sd = 0.5)[g1]
  b2 <- rnorm(10, sd = 0.8)[g2]
  y <- rpois(n, lambda = exp(1 + b1 + b2))
  d <- data.frame(y = y, g1 = g1, g2 = g2)

  fit_a <- lme4::glmer(y ~ 1 + (1 | g1) + (1 | g2),
                        data = d, family = poisson())
  fit_b <- lme4::glmer(y ~ 1 + (1 | g2) + (1 | g1),
                        data = d, family = poisson())
  expect_equal(pred_glmer(fit_a), pred_glmer(fit_b),
               tolerance = 1e-6)
})

test_that("pred_glmmtmb invariant to RE term order", {
  skip_on_cran()
  set.seed(1)
  n <- 300
  g1 <- factor(sample(1:15, n, replace = TRUE))
  g2 <- factor(sample(1:10, n, replace = TRUE))
  b1 <- rnorm(15, sd = 0.5)[g1]
  b2 <- rnorm(10, sd = 0.8)[g2]
  y <- rpois(n, lambda = exp(1 + b1 + b2))
  d <- data.frame(y = y, g1 = g1, g2 = g2)

  fit_a <- glmmTMB::glmmTMB(y ~ 1 + (1 | g1) + (1 | g2),
                              data = d, family = poisson())
  fit_b <- glmmTMB::glmmTMB(y ~ 1 + (1 | g2) + (1 | g1),
                              data = d, family = poisson())
  expect_equal(pred_glmmtmb(fit_a), pred_glmmtmb(fit_b),
               tolerance = 1e-6)
})

test_that("pred_glmer handles random slopes", {
  skip_on_cran()
  set.seed(2)
  n <- 300
  g <- factor(rep(1:30, each = 10))
  x <- rnorm(n)
  b0 <- rnorm(30, sd = 0.5)[g]
  b1 <- rnorm(30, sd = 0.3)[g]
  y <- rpois(n, lambda = exp(0.5 + 0.3 * x + b0 + b1 * x))
  d <- data.frame(y = y, x = x, g = g)

  fit <- lme4::glmer(y ~ x + (1 + x | g),
                      data = d, family = poisson())
  pred <- pred_glmer(fit)
  expect_length(pred, n)
  expect_true(all(pred > 0))
})

test_that("pred_glmmtmb handles random slopes", {
  skip_on_cran()
  set.seed(2)
  n <- 300
  g <- factor(rep(1:30, each = 10))
  x <- rnorm(n)
  b0 <- rnorm(30, sd = 0.5)[g]
  b1 <- rnorm(30, sd = 0.3)[g]
  y <- rpois(n, lambda = exp(0.5 + 0.3 * x + b0 + b1 * x))
  d <- data.frame(y = y, x = x, g = g)

  fit <- glmmTMB::glmmTMB(y ~ x + (1 + x | g),
                            data = d, family = poisson())
  pred <- pred_glmmtmb(fit)
  expect_length(pred, n)
  expect_true(all(pred > 0))
})

test_that("pred_glmer and pred_glmmtmb agree with random slopes", {
  skip_on_cran()
  set.seed(2)
  n <- 300
  g <- factor(rep(1:30, each = 10))
  x <- rnorm(n)
  b0 <- rnorm(30, sd = 0.5)[g]
  b1 <- rnorm(30, sd = 0.3)[g]
  y <- rpois(n, lambda = exp(0.5 + 0.3 * x + b0 + b1 * x))
  d <- data.frame(y = y, x = x, g = g)

  fit4 <- lme4::glmer(y ~ x + (1 + x | g),
                       data = d, family = poisson())
  fit_t <- glmmTMB::glmmTMB(y ~ x + (1 + x | g),
                              data = d, family = poisson())
  expect_equal(pred_glmer(fit4), pred_glmmtmb(fit_t),
               tolerance = 0.5)
})
