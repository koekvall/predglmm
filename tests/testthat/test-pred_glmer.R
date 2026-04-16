# Tests for pred_glmer

test_that("pred_glmer works on Poisson glmer", {
  skip_on_cran()
  d <- data.frame(
    y = rpois(200, lambda = exp(1 + rep(rnorm(20, sd = 0.5), each = 10))),
    group = factor(rep(1:20, each = 10))
  )
  fit <- lme4::glmer(y ~ 1 + (1 | group), data = d, family = poisson())
  pred <- pred_glmer(fit)

  expect_length(pred, 200)
  expect_true(all(pred > 0))
  # Marginal mean should exceed population mean for log link
  pop <- exp(as.numeric(lme4::getME(fit, "beta")))
  expect_true(mean(pred) > pop)
})

test_that("pred_glmer on lmerMod returns X %*% beta", {
  skip_on_cran()
  d <- data.frame(
    y = rnorm(100, mean = rep(rnorm(10), each = 10)),
    group = factor(rep(1:10, each = 10))
  )
  fit <- lme4::lmer(y ~ 1 + (1 | group), data = d)
  pred <- pred_glmer(fit)
  expected <- as.vector(
    lme4::getME(fit, "X") %*% lme4::getME(fit, "beta")
  )
  expect_equal(pred, expected)
})

test_that("pred_glmer rejects non-lme4 objects", {
  expect_error(pred_glmer(lm(1:10 ~ 1)), "lmerMod or glmerMod")
})
