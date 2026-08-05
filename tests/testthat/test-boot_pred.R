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

# Reproduce the first bootstrap replicate outside boot_pred: same RNG setup,
# one simulate() draw, returned with the caller's RNG state restored.
first_sim <- function(fit, seed) {
  old_seed <- get(".Random.seed", envir = globalenv(), inherits = FALSE)
  on.exit(assign(".Random.seed", old_seed, envir = globalenv()))
  RNGkind("L'Ecuyer-CMRG")
  set.seed(seed)
  stats::simulate(fit)[[1]]
}

test_that("weighted binomial glmmTMB refits use proportions, not counts", {
  skip_on_cran()
  set.seed(12)
  n <- 120
  g <- factor(rep(1:12, each = 10))
  w <- sample(5:15, n, replace = TRUE)
  p <- plogis(0.3 + rnorm(12, sd = 0.6)[g])
  d <- data.frame(yprop = rbinom(n, w, p) / w, g = g, w = w)
  fit <- glmmTMB::glmmTMB(yprop ~ 1 + (1 | g), weights = w, data = d,
                          family = binomial())
  b <- boot_pred(fit, n_boot = 3, seed = 1)
  expect_equal(b$n_fail, 0L)

  # simulate() returns a (successes, failures) matrix; the refit must see
  # proportions with the original weights
  ysim <- first_sim(fit, seed = 1)
  expect_true(is.matrix(ysim))
  d_ref <- d
  d_ref$yprop <- ysim[, 1] / rowSums(ysim)
  refit <- suppressWarnings(
    glmmTMB::glmmTMB(yprop ~ 1 + (1 | g), weights = w, data = d_ref,
                     family = binomial())
  )
  expect_equal(unname(b$draws[, 1]), pred_glmmtmb(refit), tolerance = 1e-6)
})

test_that("glmmTMB refits align simulated responses with NA-dropped rows", {
  skip_on_cran()
  set.seed(13)
  d <- data.frame(y = rpois(80, 3), x = rnorm(80),
                  g = factor(rep(1:8, each = 10)))
  d$x[c(3, 17, 55)] <- NA  # incomplete rows dropped by the fit
  fit <- glmmTMB::glmmTMB(y ~ x + (1 | g), data = d, family = poisson())
  expect_equal(stats::nobs(fit), 77)
  b <- boot_pred(fit, n_boot = 3, seed = 2)
  expect_equal(nrow(b$draws), 77)
  expect_equal(b$n_fail, 0L)

  # First replicate reproduced by hand: simulated values go to the used rows
  ysim <- first_sim(fit, seed = 2)
  d_ref <- d
  d_ref$y <- NA_real_
  d_ref$y[stats::complete.cases(d[, c("y", "x", "g")])] <- ysim
  refit <- glmmTMB::glmmTMB(y ~ x + (1 | g), data = d_ref,
                            family = poisson())
  expect_equal(unname(b$draws[, 1]), pred_glmmtmb(refit), tolerance = 1e-6)
})

test_that("statistic sees data columns absent from the formula (glmmTMB)", {
  skip_on_cran()
  d <- make_boot_data()
  d$bin <- rep(1:10, times = 8)  # not in the model formula
  fit <- glmmTMB::glmmTMB(y ~ x + (1 | group), data = d, family = poisson())
  cellmean <- function(data) c(bin1 = mean(data$pred[data$bin == 1]))
  b <- boot_pred(fit, statistic = cellmean, n_boot = 3, seed = 1)
  expect_equal(unname(b$estimate), mean(pred_glmmtmb(fit)[d$bin == 1]))
  expect_equal(ncol(b$draws), 3 - b$n_fail)
})

test_that("statistic sees data columns absent from the formula (glmer)", {
  skip_on_cran()
  d <- make_boot_data()
  d$bin <- rep(1:10, times = 8)  # not in the model formula
  fit <- lme4::glmer(y ~ x + (1 | group), data = d, family = poisson())
  cellmean <- function(data) c(bin1 = mean(data$pred[data$bin == 1]))
  b <- boot_pred(fit, statistic = cellmean, n_boot = 3, seed = 1)
  expect_equal(unname(b$estimate), mean(pred_glmer(fit)[d$bin == 1]))
  expect_equal(ncol(b$draws), 3 - b$n_fail)
})

test_that("dropped rows carry NA predictions in the statistic's data", {
  skip_on_cran()
  set.seed(13)
  d <- data.frame(y = rpois(80, 3), x = rnorm(80),
                  g = factor(rep(1:8, each = 10)),
                  row_id = seq_len(80))  # not in the model formula
  na_rows <- c(3L, 17L, 55L)
  d$x[na_rows] <- NA  # incomplete rows dropped by the fit
  fit <- glmmTMB::glmmTMB(y ~ x + (1 | g), data = d, family = poisson())

  seen <- NULL
  probe <- function(data) {
    seen <<- data
    mean(data$pred, na.rm = TRUE)
  }
  b <- boot_pred(fit, statistic = probe, n_boot = 2, seed = 2)

  # The statistic sees the full 80-row data, with NA predictions exactly
  # in the dropped rows and all original columns intact
  expect_equal(nrow(seen), 80)
  expect_identical(which(is.na(seen$pred)), na_rows)
  expect_equal(seen$row_id, d$row_id)
  # The default statistic still aligns with the fitted observations
  b0 <- boot_pred(fit, n_boot = 2, seed = 2)
  expect_equal(length(b0$estimate), stats::nobs(fit))
  expect_equal(b0$estimate, pred_glmmtmb(fit))
})

test_that("glmer fits without a data frame fall back to the model frame", {
  skip_on_cran()
  d <- make_boot_data()
  # Variables in an environment, no data argument: nothing to recover, so
  # the statistic sees the model frame
  fit <- local({
    y <- d$y; x <- d$x; group <- d$group
    lme4::glmer(y ~ x + (1 | group), family = poisson())
  })
  b <- boot_pred(fit, n_boot = 2, seed = 1)
  expect_equal(b$estimate, pred_glmer(fit))
  expect_equal(ncol(b$draws), 2 - b$n_fail)
})

test_that("forked parallel replicates run and keep dimensions", {
  skip_on_cran()
  skip_on_os("windows")
  d <- make_boot_data()
  fit <- glmmTMB::glmmTMB(y ~ x + (1 | group), data = d, family = poisson())
  b <- boot_pred(fit, n_boot = 4, seed = 1, cores = 2)
  expect_equal(nrow(b$draws), nrow(d))
  expect_equal(ncol(b$draws), 4 - b$n_fail)
})

test_that("se = FALSE glmmTMB fits give NA n_boundary and still print", {
  skip_on_cran()
  d <- make_boot_data()
  fit <- glmmTMB::glmmTMB(y ~ x + (1 | group), data = d,
                          family = poisson(), se = FALSE)
  b <- boot_pred(fit, n_boot = 3, seed = 1)
  expect_true(is.na(b$n_boundary))
  expect_output(print(b), "boundary status unknown")
})
