# Strong-correctness tests for pred_glmer and pred_glmmtmb across
# more complicated random-effect structures (crossed, nested, slopes,
# multiple grouping factors). For log and sqrt links the marginal mean
# has a closed form: we check that the implementation matches that form
# to machine tolerance, with sigma^2 computed via an independent helper
# (see helper-analytical.R).

# Each block builds a small simulated dataset, fits both lme4 and glmmTMB,
# and checks both adapters against the analytical formula.

simulate_poisson <- function(eta_total, link = "log") {
  mu <- if (link == "log") exp(eta_total) else pmax(eta_total, 0)^2
  rpois(length(mu), lambda = mu)
}

test_that("crossed grouping (1|g1) + (1|g2), Poisson-log", {
  skip_on_cran()
  set.seed(101)
  n <- 400
  g1 <- factor(sample(1:20, n, replace = TRUE))
  g2 <- factor(sample(1:15, n, replace = TRUE))
  b1 <- rnorm(20, sd = 0.4)[g1]
  b2 <- rnorm(15, sd = 0.3)[g2]
  d <- data.frame(y = simulate_poisson(0.5 + b1 + b2),
                  g1 = g1, g2 = g2)

  fit4 <- lme4::glmer(y ~ 1 + (1 | g1) + (1 | g2),
                       data = d, family = poisson())
  fit_t <- glmmTMB::glmmTMB(y ~ 1 + (1 | g1) + (1 | g2),
                              data = d, family = poisson())

  expect_equal(pred_glmer(fit4),
               analytical_pred(eta_glmer(fit4),
                               sigma_sq_glmer(fit4), "log"),
               tolerance = 1e-10)
  expect_equal(pred_glmmtmb(fit_t),
               analytical_pred(eta_glmmtmb(fit_t),
                               sigma_sq_glmmtmb(fit_t), "log"),
               tolerance = 1e-10)
})

test_that("nested grouping (1 | g1/g2), Poisson-log", {
  skip_on_cran()
  set.seed(102)
  n_outer <- 12
  n_inner <- 5
  n_per <- 8
  d <- expand.grid(rep = seq_len(n_per),
                   inner = seq_len(n_inner),
                   outer = seq_len(n_outer))
  d$g1 <- factor(d$outer)
  d$g2 <- factor(paste(d$outer, d$inner, sep = ":"))
  b1 <- rnorm(nlevels(d$g1), sd = 0.4)[d$g1]
  b2 <- rnorm(nlevels(d$g2), sd = 0.3)[d$g2]
  d$y <- simulate_poisson(0.3 + b1 + b2)

  fit4 <- lme4::glmer(y ~ 1 + (1 | g1 / g2),
                       data = d, family = poisson())
  fit_t <- glmmTMB::glmmTMB(y ~ 1 + (1 | g1 / g2),
                              data = d, family = poisson())

  expect_equal(pred_glmer(fit4),
               analytical_pred(eta_glmer(fit4),
                               sigma_sq_glmer(fit4), "log"),
               tolerance = 1e-10)
  expect_equal(pred_glmmtmb(fit_t),
               analytical_pred(eta_glmmtmb(fit_t),
                               sigma_sq_glmmtmb(fit_t), "log"),
               tolerance = 1e-10)
})

test_that("random slopes (1 + x | g), Poisson-log", {
  skip_on_cran()
  set.seed(103)
  n <- 400
  g <- factor(rep(1:25, each = 16))
  x <- rnorm(n)
  Psi <- matrix(c(0.4^2, 0.5 * 0.4 * 0.3,
                  0.5 * 0.4 * 0.3, 0.3^2), 2, 2)
  L <- t(chol(Psi))
  B <- L %*% matrix(rnorm(2 * 25), nrow = 2)
  b0 <- B[1, ][g]
  b1 <- B[2, ][g]
  d <- data.frame(y = simulate_poisson(0.2 + 0.3 * x + b0 + b1 * x),
                  x = x, g = g)

  fit4 <- lme4::glmer(y ~ x + (1 + x | g),
                       data = d, family = poisson())
  fit_t <- glmmTMB::glmmTMB(y ~ x + (1 + x | g),
                              data = d, family = poisson())

  expect_equal(pred_glmer(fit4),
               analytical_pred(eta_glmer(fit4),
                               sigma_sq_glmer(fit4), "log"),
               tolerance = 1e-10)
  expect_equal(pred_glmmtmb(fit_t),
               analytical_pred(eta_glmmtmb(fit_t),
                               sigma_sq_glmmtmb(fit_t), "log"),
               tolerance = 1e-10)
})

test_that("slopes plus extra grouping factor, Poisson-log", {
  skip_on_cran()
  set.seed(104)
  n <- 500
  g1 <- factor(rep(1:25, each = 20))
  g2 <- factor(sample(1:10, n, replace = TRUE))
  x <- rnorm(n)
  b0 <- rnorm(25, sd = 0.3)[g1]
  b1 <- rnorm(25, sd = 0.2)[g1]
  c2 <- rnorm(10, sd = 0.25)[g2]
  d <- data.frame(y = simulate_poisson(0.4 + 0.3 * x +
                                         b0 + b1 * x + c2),
                  x = x, g1 = g1, g2 = g2)

  fit4 <- lme4::glmer(y ~ x + (1 + x | g1) + (1 | g2),
                       data = d, family = poisson())
  fit_t <- glmmTMB::glmmTMB(y ~ x + (1 + x | g1) + (1 | g2),
                              data = d, family = poisson())

  expect_equal(pred_glmer(fit4),
               analytical_pred(eta_glmer(fit4),
                               sigma_sq_glmer(fit4), "log"),
               tolerance = 1e-10)
  expect_equal(pred_glmmtmb(fit_t),
               analytical_pred(eta_glmmtmb(fit_t),
                               sigma_sq_glmmtmb(fit_t), "log"),
               tolerance = 1e-10)
})

test_that("three crossed grouping factors, Poisson-log", {
  skip_on_cran()
  set.seed(105)
  n <- 600
  g1 <- factor(sample(1:15, n, replace = TRUE))
  g2 <- factor(sample(1:12, n, replace = TRUE))
  g3 <- factor(sample(1:10, n, replace = TRUE))
  b1 <- rnorm(15, sd = 0.3)[g1]
  b2 <- rnorm(12, sd = 0.25)[g2]
  b3 <- rnorm(10, sd = 0.2)[g3]
  d <- data.frame(y = simulate_poisson(0.3 + b1 + b2 + b3),
                  g1 = g1, g2 = g2, g3 = g3)

  fit4 <- lme4::glmer(y ~ 1 + (1 | g1) + (1 | g2) + (1 | g3),
                       data = d, family = poisson())
  fit_t <- glmmTMB::glmmTMB(y ~ 1 + (1 | g1) + (1 | g2) + (1 | g3),
                              data = d, family = poisson())

  expect_equal(pred_glmer(fit4),
               analytical_pred(eta_glmer(fit4),
                               sigma_sq_glmer(fit4), "log"),
               tolerance = 1e-10)
  expect_equal(pred_glmmtmb(fit_t),
               analytical_pred(eta_glmmtmb(fit_t),
                               sigma_sq_glmmtmb(fit_t), "log"),
               tolerance = 1e-10)
})

test_that("crossed grouping with sqrt link, Poisson-sqrt", {
  skip_on_cran()
  set.seed(106)
  n <- 400
  g1 <- factor(sample(1:20, n, replace = TRUE))
  g2 <- factor(sample(1:15, n, replace = TRUE))
  b1 <- rnorm(20, sd = 0.3)[g1]
  b2 <- rnorm(15, sd = 0.25)[g2]
  d <- data.frame(y = simulate_poisson(2 + b1 + b2, link = "sqrt"),
                  g1 = g1, g2 = g2)

  fit4 <- lme4::glmer(y ~ 1 + (1 | g1) + (1 | g2),
                       data = d, family = poisson(link = "sqrt"))
  fit_t <- glmmTMB::glmmTMB(y ~ 1 + (1 | g1) + (1 | g2),
                              data = d,
                              family = poisson(link = "sqrt"))

  expect_equal(pred_glmer(fit4),
               analytical_pred(eta_glmer(fit4),
                               sigma_sq_glmer(fit4), "sqrt"),
               tolerance = 1e-10)
  expect_equal(pred_glmmtmb(fit_t),
               analytical_pred(eta_glmmtmb(fit_t),
                               sigma_sq_glmmtmb(fit_t), "sqrt"),
               tolerance = 1e-10)
})
