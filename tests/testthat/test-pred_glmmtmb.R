# Tests for pred_glmmtmb

test_that("pred_glmmtmb works on Poisson glmmTMB", {
  skip_on_cran()
  d <- data.frame(
    y = rpois(200, lambda = exp(1 + rep(rnorm(20, sd = 0.5), each = 10))),
    group = factor(rep(1:20, each = 10))
  )
  fit <- glmmTMB::glmmTMB(y ~ 1 + (1 | group), data = d,
                           family = poisson())
  pred <- pred_glmmtmb(fit)

  expect_length(pred, 200)
  expect_true(all(pred > 0))
  pop <- exp(as.numeric(glmmTMB::getME(fit, "beta")))
  expect_true(mean(pred) > pop)
})

test_that("pred_glmmtmb agrees with pred_glmer", {
  skip_on_cran()
  set.seed(123)
  d <- data.frame(
    y = rpois(200, lambda = exp(1 + rep(rnorm(20, sd = 0.5), each = 10))),
    group = factor(rep(1:20, each = 10))
  )
  fit4 <- lme4::glmer(y ~ 1 + (1 | group), data = d,
                       family = poisson())
  fit_t <- glmmTMB::glmmTMB(y ~ 1 + (1 | group), data = d,
                             family = poisson())
  # Both should give similar marginal predictions
  expect_equal(mean(pred_glmer(fit4)), mean(pred_glmmtmb(fit_t)),
               tolerance = 0.1)
})

test_that("pred_glmmtmb rejects non-glmmTMB objects", {
  expect_error(pred_glmmtmb(lm(1:10 ~ 1)))
})
