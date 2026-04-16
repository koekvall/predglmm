# Simulation: MSE of fitted values vs true marginal mean
# as a function of random effect variance.
#
# Model: E[Y | b] = exp(beta + b), b ~ N(0, sigma^2)
# True marginal mean: mu = exp(beta + sigma^2 / 2)
#
# lme4 and glmmTMB: Poisson GLMM with log link
# nlme: Gaussian nonlinear model y = exp(beta + b) + e
#
# MSE = mean of (y_hat_i - mu_true)^2 over observations i
#
# Usage: Rscript sim_nlme.R
# Output: sim_mse.pdf

library(lme4)
library(glmmTMB)
library(nlme)
library(ggplot2)
library(patchwork)
devtools::load_all(".")

beta0 <- 1.5
sigma_vals <- seq(0.2, 1.2, by = 0.2)
sigma_e <- 1         # residual SD for nlme
n_groups <- 50
n_per_group <- 10
n_obs <- n_groups * n_per_group
n_sim <- 300

results <- data.frame()

set.seed(42)
for (sigma_re in sigma_vals) {
  mu_true <- exp(beta0 + sigma_re^2 / 2)

  for (i in seq_len(n_sim)) {
    group <- rep(1:n_groups, each = n_per_group)
    b <- rnorm(n_groups, sd = sigma_re)
    eta <- beta0 + b[group]

    y_pois <- rpois(n_obs, lambda = exp(eta))
    d_pois <- data.frame(y = y_pois, group = factor(group))

    y_gauss <- exp(eta) + rnorm(n_obs, sd = sigma_e)
    d_gauss <- data.frame(y = y_gauss, group = factor(group))

    # ── lme4 ──
    fit4 <- tryCatch(
      glmer(y ~ 1 + (1 | group), data = d_pois,
            family = poisson()),
      error = function(e) NULL)
    if (!is.null(fit4)) {
      bhat <- as.numeric(getME(fit4, "beta"))
      p_pop <- rep(exp(bhat), n_obs)
      p_eb  <- fitted(fit4)
      p_mar <- pred_glmer(fit4)
      results <- rbind(results, data.frame(
        sigma = sigma_re, package = "lme4",
        method = "Zero RE",
        mse = mean((p_pop - mu_true)^2)))
      results <- rbind(results, data.frame(
        sigma = sigma_re, package = "lme4",
        method = "Conditional",
        mse = mean((p_eb - mu_true)^2)))
      results <- rbind(results, data.frame(
        sigma = sigma_re, package = "lme4",
        method = "Marginal",
        mse = mean((p_mar - mu_true)^2)))
    }

    # ── glmmTMB ──
    fit_t <- tryCatch(
      glmmTMB(y ~ 1 + (1 | group), data = d_pois,
              family = poisson()),
      error = function(e) NULL)
    if (!is.null(fit_t)) {
      bhat <- as.numeric(glmmTMB::getME(fit_t, "beta"))
      p_pop <- rep(exp(bhat), n_obs)
      p_eb  <- predict(fit_t, type = "response")
      p_mar <- pred_glmmtmb(fit_t)
      results <- rbind(results, data.frame(
        sigma = sigma_re, package = "glmmTMB",
        method = "Zero RE",
        mse = mean((p_pop - mu_true)^2)))
      results <- rbind(results, data.frame(
        sigma = sigma_re, package = "glmmTMB",
        method = "Conditional",
        mse = mean((p_eb - mu_true)^2)))
      results <- rbind(results, data.frame(
        sigma = sigma_re, package = "glmmTMB",
        method = "Marginal",
        mse = mean((p_mar - mu_true)^2)))
    }

    # ── nlme ──
    fit_nl <- tryCatch(
      nlme(y ~ exp(b0),
           fixed = b0 ~ 1,
           random = b0 ~ 1 | group,
           data = d_gauss,
           start = c(b0 = 1.5)),
      error = function(e) NULL)
    if (!is.null(fit_nl)) {
      p_pop <- as.numeric(predict(fit_nl, level = 0))
      p_eb  <- as.numeric(predict(fit_nl, level = 1))
      p_mar <- pred_nlme(fit_nl, nsim = 5000)
      results <- rbind(results, data.frame(
        sigma = sigma_re, package = "nlme",
        method = "Zero RE",
        mse = mean((p_pop - mu_true)^2)))
      results <- rbind(results, data.frame(
        sigma = sigma_re, package = "nlme",
        method = "Conditional",
        mse = mean((p_eb - mu_true)^2)))
      results <- rbind(results, data.frame(
        sigma = sigma_re, package = "nlme",
        method = "Marginal",
        mse = mean((p_mar - mu_true)^2)))
    }
  }
  cat(sprintf("Done: sigma = %.1f\n", sigma_re))
}

# ── Average MSE over replications ──
library(dplyr)
avg <- results %>%
  group_by(sigma, package, method) %>%
  summarize(mse = mean(mse), .groups = "drop")

# ── Plot ──
col_vals <- c("Marginal"    = "#2166ac",
              "Zero RE"     = "#d6604d",
              "Conditional" = "#bababa")

common_theme <- list(
  theme_bw(base_size = 10),
  theme(legend.position  = "bottom",
        panel.grid.minor = element_blank(),
        plot.title       = element_text(size = 10))
)

make_panel <- function(data, title) {
  ggplot(data, aes(x = sigma, y = mse,
                   color = method, shape = method)) +
    geom_line(linewidth = 0.6) +
    geom_point(size = 2) +
    scale_y_log10() +
    scale_color_manual(values = col_vals) +
    labs(x = expression(sigma[b]),
         y = "MSE (log scale)",
         color = "Fitted value", shape = "Fitted value",
         title = title) +
    common_theme
}

p_a <- make_panel(avg[avg$package == "lme4", ],
                  "(a) lme4")
p_b <- make_panel(avg[avg$package == "glmmTMB", ],
                  "(b) glmmTMB")
p_c <- make_panel(avg[avg$package == "nlme", ],
                  "(c) nlme")

p <- (p_a | p_b | p_c) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave("sim_mse.pdf", p, width = 10, height = 3.5)
cat("Saved sim_mse.pdf\n")
