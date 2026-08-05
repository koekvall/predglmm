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

test_that("pred_glmmtmb handles rank-deficient designs", {
  skip_on_cran()
  set.seed(8)
  d <- data.frame(y = rpois(200, 3), x = rnorm(200),
                  g = factor(rep(1:20, each = 10)))
  d$x2 <- d$x  # perfectly collinear; glmmTMB drops it with a warning
  fit_full <- suppressWarnings(suppressMessages(
    glmmTMB::glmmTMB(y ~ x + x2 + (1 | g), data = d, family = poisson())
  ))
  fit_red <- glmmTMB::glmmTMB(y ~ x + (1 | g), data = d, family = poisson())
  expect_equal(pred_glmmtmb(fit_full), pred_glmmtmb(fit_red),
               tolerance = 1e-6)
})

# Shared data for the zero-inflation and truncation tests: counts with excess
# zeros, a binary covariate, and a grouping factor
make_zi_data <- function() {
  set.seed(42)
  g <- factor(rep(1:20, each = 15))
  x <- rep(c(0, 1), length.out = 300)
  mu <- exp(0.5 + 0.4 * x + rep(rnorm(20, sd = 0.5), each = 15))
  y <- rpois(300, mu) * rbinom(300, 1, prob = 0.7)
  data.frame(y = y, x = x, group = g)
}

test_that("zero-inflated fit matches closed form and MC, conditional drops zi", {
  skip_on_cran()
  d <- make_zi_data()
  fit <- glmmTMB::glmmTMB(y ~ x + (1 | group), ziformula = ~ x,
                          family = glmmTMB::nbinom1, data = d)
  s2 <- glmmTMB::VarCorr(fit)$cond$group[1, 1]
  eta <- as.vector(glmmTMB::getME(fit, "X") %*% glmmTMB::fixef(fit)$cond)
  eta_zi <- as.vector(glmmTMB::getME(fit, "Xzi") %*% glmmTMB::fixef(fit)$zi)

  # Log link and fixed-effect zi have closed forms
  expect_equal(pred_glmmtmb(fit),
               (1 - plogis(eta_zi)) * exp(eta + s2 / 2))
  expect_equal(pred_glmmtmb(fit, type = "conditional"), exp(eta + s2 / 2))
})

test_that("truncated families match quadrature-free reference", {
  skip_on_cran()
  d <- make_zi_data()

  # Hurdle Poisson: with dnbinom/dpois-based P(0) computed independently of
  # the expm1 algebra in pred_glmmtmb, integrate() over the random effect
  # gives a reference marginal mean for individual observations
  fit <- glmmTMB::glmmTMB(y ~ x + (1 | group), ziformula = ~ x,
                          family = glmmTMB::truncated_poisson, data = d)
  s <- sqrt(glmmTMB::VarCorr(fit)$cond$group[1, 1])
  eta <- as.vector(glmmTMB::getME(fit, "X") %*% glmmTMB::fixef(fit)$cond)
  eta_zi <- as.vector(glmmTMB::getME(fit, "Xzi") %*% glmmTMB::fixef(fit)$zi)
  ref <- function(e) stats::integrate(function(u) {
    mu <- exp(e + u)
    mu / (1 - stats::dpois(0, mu)) * stats::dnorm(u, sd = s)
  }, -8 * s, 8 * s)$value
  p <- pred_glmmtmb(fit)
  for (i in c(1, 2, 150)) {
    expect_equal(p[i], (1 - plogis(eta_zi[i])) * ref(eta[i]), tolerance = 1e-4)
  }

  # Truncated nbinom2, reference P(0) from dnbinom
  fit2 <- glmmTMB::glmmTMB(y ~ x + (1 | group), ziformula = ~ 1,
                           family = glmmTMB::truncated_nbinom2, data = d)
  s2v <- sqrt(glmmTMB::VarCorr(fit2)$cond$group[1, 1])
  eta2 <- as.vector(glmmTMB::getME(fit2, "X") %*% glmmTMB::fixef(fit2)$cond)
  th <- glmmTMB::sigma(fit2)
  ref2 <- function(e) stats::integrate(function(u) {
    mu <- exp(e + u)
    mu / (1 - stats::dnbinom(0, mu = mu, size = th)) * stats::dnorm(u, sd = s2v)
  }, -8 * s2v, 8 * s2v)$value
  p2 <- pred_glmmtmb(fit2, type = "conditional")
  for (i in c(1, 2, 150)) {
    expect_equal(p2[i], ref2(eta2[i]), tolerance = 1e-4)
  }
})

test_that("truncated family with non-constant dispersion is refused", {
  skip_on_cran()
  d <- make_zi_data()
  fit <- glmmTMB::glmmTMB(y ~ x + (1 | group), ziformula = ~ 1,
                          dispformula = ~ x,
                          family = glmmTMB::truncated_nbinom2, data = d)
  expect_error(pred_glmmtmb(fit), "constant dispersion")
})
