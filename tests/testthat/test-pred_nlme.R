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
