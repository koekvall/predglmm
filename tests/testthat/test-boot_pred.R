# Tests for boot_pred. Fits are kept tiny and n_boot small so the full
# suite stays fast; these test plumbing, not statistical accuracy.

make_boot_data <- function() {
  set.seed(7)
  g <- factor(rep(1:8, each = 10))
  x <- factor(rep(c("a", "b"), times = 40))
  mu <- exp(1 + 0.5 * (x == "b") + rep(rnorm(8, sd = 0.3), each = 10))
  data.frame(y = rpois(80, mu), x = x, group = g)
}

test_that("default statistic gives correct dimensions and estimate (glmmTMB)", {
  skip_on_cran()
  d <- make_boot_data()
  fit <- glmmTMB::glmmTMB(y ~ x + (1 | group), data = d, family = poisson())
  b <- boot_pred(fit, n_boot = 20, seed = 1)

  expect_s3_class(b, "boot_pred")
  expect_equal(b$estimate, pred_glmmtmb(fit))
  expect_true(is.matrix(b$draws))
  expect_equal(nrow(b$draws), nrow(d))
  expect_equal(ncol(b$draws), 20 - b$n_fail)
  expect_equal(dim(b$ci), c(nrow(d), 2))
  expect_equal(colnames(b$ci), c("lower", "upper"))
  expect_equal(b$n_boot, 20L)
})

test_that("boot_pred supports glmer fits and matches pred_glmer", {
  skip_on_cran()
  d <- make_boot_data()
  fit <- lme4::glmer(y ~ x + (1 | group), data = d, family = poisson())
  b <- boot_pred(fit, n_boot = 20, seed = 1)

  expect_equal(b$estimate, pred_glmer(fit))
  expect_equal(nrow(b$draws), nrow(d))
  expect_equal(ncol(b$draws), 20 - b$n_fail)
})

test_that("custom named statistic propagates names to draws and ci", {
  skip_on_cran()
  d <- make_boot_data()
  fit <- glmmTMB::glmmTMB(y ~ x + (1 | group), data = d, family = poisson())
  contrast <- function(data) {
    m <- tapply(data$pred, data$x, mean)
    c(mean_a = unname(m["a"]), "a-b" = unname(m["a"] - m["b"]))
  }
  b <- boot_pred(fit, statistic = contrast, n_boot = 20, seed = 1)

  expect_equal(names(b$estimate), c("mean_a", "a-b"))
  expect_equal(rownames(b$ci), c("mean_a", "a-b"))
  expect_equal(rownames(b$draws), c("mean_a", "a-b"))
  expect_equal(nrow(b$draws), 2)
  # Print method runs and returns invisibly
  expect_invisible(print(b))
})

test_that("seeded serial calls are reproducible and restore the RNG state", {
  skip_on_cran()
  d <- make_boot_data()
  fit <- lme4::glmer(y ~ x + (1 | group), data = d, family = poisson())

  set.seed(99)
  before <- get(".Random.seed", envir = globalenv())
  b1 <- boot_pred(fit, n_boot = 20, seed = 2)
  after <- get(".Random.seed", envir = globalenv())
  expect_identical(before, after)

  b2 <- boot_pred(fit, n_boot = 20, seed = 2)
  expect_identical(b1$draws, b2$draws)
  expect_identical(b1$ci, b2$ci)
  expect_identical(b1$n_fail, b2$n_fail)
  expect_identical(b1$n_boundary, b2$n_boundary)
})

test_that("percentile and normal CIs use the same draws but differ", {
  skip_on_cran()
  d <- make_boot_data()
  fit <- lme4::glmer(y ~ x + (1 | group), data = d, family = poisson())
  contrast <- function(data) {
    m <- tapply(data$pred, data$x, mean)
    c(mean_a = unname(m["a"]), "a-b" = unname(m["a"] - m["b"]))
  }
  bp <- boot_pred(fit, statistic = contrast, n_boot = 30, seed = 3,
                  type = "percentile")
  bn <- boot_pred(fit, statistic = contrast, n_boot = 30, seed = 3,
                  type = "normal")

  expect_true(all(is.finite(bp$ci)))
  expect_true(all(is.finite(bn$ci)))
  # The normal interval is symmetric about the estimate by construction
  expect_true(all(bn$ci[, "lower"] <= bn$estimate))
  expect_true(all(bn$ci[, "upper"] >= bn$estimate))
  # Same seed, so the draws agree; only the interval type differs
  expect_identical(bp$draws, bn$draws)
  expect_false(identical(bp$ci, bn$ci))
})

test_that("boundary (singular) refits are kept and counted, not dropped", {
  skip_on_cran()
  # True group SD is zero, so many parametric refits estimate it at zero
  set.seed(11)
  d0 <- data.frame(y = rpois(60, exp(1)), group = factor(rep(1:6, each = 10)))
  fit <- suppressMessages(
    lme4::glmer(y ~ 1 + (1 | group), data = d0, family = poisson())
  )
  b <- boot_pred(fit, n_boot = 20, seed = 4)

  expect_equal(ncol(b$draws), 20 - b$n_fail)
  expect_true(b$n_boundary > 0)
  expect_true(b$n_boundary <= ncol(b$draws))
})

test_that("unsupported classes are rejected with a clear message", {
  expect_error(boot_pred(lm(1:10 ~ 1)), "glmmTMB or merMod")
})
