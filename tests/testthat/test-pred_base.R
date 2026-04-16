# Tests for pred_base: analytical formulas and quadrature

test_that("identity link returns eta unchanged", {
  eta <- c(-1, 0, 2.5)
  sigma <- c(0.5, 1, 2)
  expect_equal(pred_base(eta, sigma, link_name = "identity"), eta)
})

test_that("log link returns exp(eta + sigma^2/2)", {
  eta <- c(0, 1, 2)
  sigma <- c(0.5, 1, 1.5)
  expected <- exp(eta + 0.5 * sigma^2)
  expect_equal(pred_base(eta, sigma, link_name = "log"), expected)
})

test_that("sqrt link returns eta^2 + sigma^2", {
  eta <- c(1, 2, 3)
  sigma <- c(0.5, 1, 2)
  expected <- eta^2 + sigma^2
  expect_equal(pred_base(eta, sigma, link_name = "sqrt"), expected)
})

test_that("sigma defaults to zero", {
  eta <- c(1, 2)
  expect_equal(pred_base(eta, link_name = "log"), exp(eta))
})

test_that("NA sigma values are replaced with zero", {
  eta <- c(1, 2)
  sigma <- c(NA, 0.5)
  expect_warning(
    result <- pred_base(eta, sigma, link_name = "log"),
    "NA values"
  )
  expect_equal(result, exp(c(1, 2 + 0.5 * 0.5^2)))
})

test_that("custom inv_link via quadrature matches log link", {
  eta <- c(0, 1, 2)
  sigma <- c(0.5, 0.8, 1.0)
  analytical <- pred_base(eta, sigma, link_name = "log")
  numerical <- pred_base(eta, sigma, inv_link = exp,
                         num_nodes = 30)
  expect_equal(numerical, analytical, tolerance = 1e-6)
})

test_that("error when neither link_name nor inv_link given", {
  expect_error(pred_base(c(1, 2)), "Exactly one")
})

test_that("unsupported link_name returns NA with warning", {
  expect_warning(
    result <- pred_base(c(1), c(0), link_name = "probit"),
    "not implemented"
  )
  expect_true(is.na(result))
})
