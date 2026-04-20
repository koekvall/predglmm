# Tests for pred_nlme

test_that("pred_nlme on lme returns X %*% beta", {
  skip_on_cran()
  fit <- nlme::lme(distance ~ age,
                   random = ~ 1 | Subject,
                   data = nlme::Orthodont)
  pred <- pred_nlme(fit)
  X <- stats::model.matrix(fit, data = fit$data)
  beta <- nlme::fixef(fit)
  expected <- as.vector(X %*% beta)
  expect_equal(pred, expected)
})

test_that("pred_nlme on nlme returns positive values", {
  skip_on_cran()
  fit <- nlme::nlme(
    height ~ SSasymp(age, Asym, R0, lrc),
    fixed = Asym + R0 + lrc ~ 1,
    random = Asym ~ 1 | Seed,
    data = Loblolly
  )
  set.seed(1)
  pred <- pred_nlme(fit, nsim = 2000)
  expect_length(pred, nrow(Loblolly))
  expect_true(all(pred > 0))
})

test_that("pred_nlme MC differs from population for nonlinear model", {
  skip_on_cran()
  fit <- nlme::nlme(
    height ~ SSasymp(age, Asym, R0, lrc),
    fixed = Asym + R0 + lrc ~ 1,
    random = lrc ~ 1 | Seed,
    data = Loblolly
  )
  set.seed(1)
  pred_mc <- pred_nlme(fit, nsim = 5000)
  pred_pop <- as.numeric(predict(fit, level = 0))
  # Jensen's inequality: they should differ when RE is on
  # the nonlinear parameter lrc
  expect_false(isTRUE(all.equal(pred_mc, pred_pop,
                                tolerance = 0.01)))
})

test_that("pred_nlme rejects non-nlme objects", {
  expect_error(pred_nlme(lm(1:10 ~ 1)), "class lme or nlme")
})

test_that("pred_nlme is reproducible with seed and does not pollute RNG", {
  skip_on_cran()
  fit <- nlme::nlme(
    height ~ SSasymp(age, Asym, R0, lrc),
    fixed = Asym + R0 + lrc ~ 1,
    random = lrc ~ 1 | Seed,
    data = Loblolly
  )
  set.seed(99)
  before <- .Random.seed
  p1 <- pred_nlme(fit, nsim = 200, seed = 1)
  p2 <- pred_nlme(fit, nsim = 200, seed = 1)
  after <- .Random.seed
  expect_identical(p1, p2)
  expect_identical(before, after)
})

test_that("pred_nlme nested closed form: y = exp(b) gives sum-of-variances", {
  skip_on_cran()
  # Constant fixed mean with nested random effects on the same parameter:
  #   y_i = exp(beta + b1_{g1(i)} + b2_{g2(i)}) + epsilon_i
  # Marginal mean = exp(beta + (sigma1^2 + sigma2^2) / 2). All observations
  # share the same marginal mean. This isolates the multi-level summation
  # of Psi: if pred_nlme used only reStruct[[1]], the predicted mean would
  # be too small.
  set.seed(11)
  n_outer <- 12
  n_inner <- 5
  n_per <- 6
  d_local <- expand.grid(rep = seq_len(n_per),
                         inner = seq_len(n_inner),
                         outer = seq_len(n_outer))
  d_local$g1 <- factor(d_local$outer)
  d_local$g2 <- factor(paste(d_local$outer, d_local$inner, sep = ":"))
  beta_true <- 0.8
  sd1 <- 0.4
  sd2 <- 0.3
  b1 <- rnorm(nlevels(d_local$g1), sd = sd1)[d_local$g1]
  b2 <- rnorm(nlevels(d_local$g2), sd = sd2)[d_local$g2]
  d_local$y <- exp(beta_true + b1 + b2) + rnorm(nrow(d_local), sd = 0.05)

  # nlme::getData() evaluates mCall$data in the formula's environment; when
  # called from inside testthat the local frame is no longer reachable by the
  # time getData runs. Park the data on the global env for the duration of
  # the test so getData can find it.
  assign("d_nlme_test", d_local, envir = globalenv())
  on.exit(rm("d_nlme_test", envir = globalenv()), add = TRUE)

  fit <- nlme::nlme(
    y ~ exp(b),
    fixed = b ~ 1,
    random = b ~ 1 | g1 / g2,
    data = d_nlme_test,
    start = c(b = beta_true)
  )
  expect_length(fit$modelStruct$reStruct, 2)

  beta_hat <- nlme::fixef(fit)["b"]
  Psi1 <- as.numeric(as.matrix(fit$modelStruct$reStruct[[1]])) * fit$sigma^2
  Psi2 <- as.numeric(as.matrix(fit$modelStruct$reStruct[[2]])) * fit$sigma^2
  expected <- as.numeric(exp(beta_hat + (Psi1 + Psi2) / 2))

  pred <- pred_nlme(fit, nsim = 5000, seed = 17)
  expect_equal(mean(pred), expected, tolerance = 0.02)
  # Single-level fallback would predict exp(beta + max(Psi1, Psi2)/2),
  # noticeably smaller. Verify the gap.
  single_level <- as.numeric(exp(beta_hat + max(Psi1, Psi2) / 2))
  expect_gt(expected, single_level)
  expect_gt(mean(pred), single_level)
})

test_that("pred_nlme handles nested grouping factors", {
  skip_on_cran()
  # nlme::Pixel has Dog and Side (nested within Dog). Using SSlogis gives a
  # linear-in-Asym mean, so the marginal mean equals the population mean and
  # we can check correctness against a closed form. The test verifies that
  # (i) reStruct really has two levels for this fit and (ii) predictions are
  # finite and consistent with the population mean.
  fit_nested <- nlme::nlme(
    pixel ~ A + B * day + C * day^2,
    fixed = A + B + C ~ 1,
    random = A ~ 1 | Dog/Side,
    data = nlme::Pixel,
    start = c(A = 1100, B = 8, C = -0.5)
  )
  expect_length(fit_nested$modelStruct$reStruct, 2)

  pred_full <- pred_nlme(fit_nested, nsim = 1000, seed = 3)
  expect_length(pred_full, nrow(nlme::Pixel))
  expect_true(all(is.finite(pred_full)))

  beta <- nlme::fixef(fit_nested)
  pop <- with(nlme::Pixel,
              beta["A"] + beta["B"] * day + beta["C"] * day^2)
  expect_equal(pred_full, as.numeric(pop), tolerance = 1)
})
