# Offsets must enter the linear predictor of marginal predictions.
# For Poisson-log fits the marginal mean has the closed form
# exp(offset + eta + sigma^2 / 2), computed here from VarCorr directly.

make_offset_data <- function() {
  set.seed(31)
  n <- 200
  expo <- runif(n, 0.5, 5)
  g <- factor(rep(1:20, each = 10))
  b <- rnorm(20, sd = 0.5)[g]
  data.frame(y = rpois(n, expo * exp(0.2 + b)), expo = expo, g = g)
}

test_that("pred_glmer includes formula offsets", {
  skip_on_cran()
  d <- make_offset_data()
  fit <- lme4::glmer(y ~ 1 + offset(log(expo)) + (1 | g), data = d,
                     family = poisson())
  s2 <- as.numeric(lme4::VarCorr(fit)$g[1, 1])
  expected <- exp(log(d$expo) + unname(lme4::fixef(fit)[1]) + s2 / 2)
  expect_equal(pred_glmer(fit), expected, tolerance = 1e-10)
})

test_that("pred_glmer on lmerMod includes offsets", {
  skip_on_cran()
  d <- make_offset_data()
  fit <- lme4::lmer(y ~ 1 + offset(expo) + (1 | g), data = d)
  expected <- d$expo + unname(lme4::fixef(fit)[1])
  expect_equal(pred_glmer(fit), expected, tolerance = 1e-10)
})

test_that("pred_glmmtmb includes offsets given in the formula", {
  skip_on_cran()
  d <- make_offset_data()
  fit <- glmmTMB::glmmTMB(y ~ 1 + offset(log(expo)) + (1 | g), data = d,
                          family = poisson())
  s2 <- as.numeric(glmmTMB::VarCorr(fit)$cond$g[1, 1])
  expected <- exp(log(d$expo) + unname(glmmTMB::fixef(fit)$cond[1]) + s2 / 2)
  expect_equal(pred_glmmtmb(fit), expected, tolerance = 1e-10)
})

test_that("pred_glmmtmb includes offsets given as an argument", {
  skip_on_cran()
  d <- make_offset_data()
  fit <- glmmTMB::glmmTMB(y ~ 1 + (1 | g), offset = log(expo), data = d,
                          family = poisson())
  s2 <- as.numeric(glmmTMB::VarCorr(fit)$cond$g[1, 1])
  expected <- exp(log(d$expo) + unname(glmmTMB::fixef(fit)$cond[1]) + s2 / 2)
  expect_equal(pred_glmmtmb(fit), expected, tolerance = 1e-10)
})

test_that("boot_pred estimate matches offset-aware predictions", {
  skip_on_cran()
  d <- make_offset_data()
  fit <- glmmTMB::glmmTMB(y ~ 1 + offset(log(expo)) + (1 | g), data = d,
                          family = poisson())
  b <- boot_pred(fit, n_boot = 3, seed = 1)
  expect_equal(b$estimate, pred_glmmtmb(fit))
})
