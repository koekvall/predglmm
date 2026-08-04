#' Calculate predictions (marginal expectations)
#'
#' @param eta An n-vector of linear predictor values (eta = X %*% beta)
#' @param sigma An n-vector of standard deviations for the random effects.
#'   Default is zeros (no random effects). NA values are replaced with zeros.
#' @param link_name Character string specifying the link function. Currently supports
#'   "identity", "log", or "sqrt". Use NA if providing a custom inv_link function.
#'   If both link_name and inv_link are provided, uses link_name if it's supported,
#'   otherwise falls back to inv_link with numerical integration.
#' @param num_nodes Number of nodes for Gaussian quadrature used to calculate predictions
#'   when using a custom inv_link function (default: 15)
#' @param inv_link Custom inverse link function. If provided, numerical integration
#'   via Gaussian quadrature is used. Use NA if providing link_name instead.
#'   Must be a vectorized function that accepts a numeric vector and returns
#'   a numeric vector of the same length.
#' @param link_warn Logical. If TRUE, shows warnings when both link_name and
#'   inv_link are provided. Default is FALSE to suppress warnings when called
#'   from adapter functions that intentionally supply both.
#' @return A vector of predicted (fitted) values on the response scale
#' @export
pred_base <- function(eta, sigma = rep(0, length(eta)), link_name = NA, num_nodes = 15,
                      inv_link = NA, link_warn = FALSE) {

  supported_links <- c("identity", "log", "sqrt")

  stopifnot(is.numeric(eta),
            is.numeric(sigma),
            length(sigma) == length(eta),
            is.numeric(num_nodes),
            length(num_nodes) == 1,
            length(link_name) == 1,
            is.function(inv_link) || is.na(inv_link),
            is.logical(link_warn) && length(link_warn) == 1)

  n <- length(eta)

  # At least one node
  num_nodes <- max(1, floor(num_nodes))

  sig_na <- is.na(sigma)
  if (any(sig_na)) {
    warning("NA values for sigma supplied; setting to zeros")
    sigma[sig_na] <- 0
  }

  # is.na() on a function warns, so test is.function() first
  if (is.na(link_name) && !is.function(inv_link)) {
    stop("Exactly one of link_name and inv_link is needed")
  }

  if (!is.na(link_name) && is.function(inv_link)) {
    if(link_name %in% supported_links) {
      inv_link <- NA
      if(link_warn)
        warning("link_name and inv_link both supplied; using supported link_name
              for speed")
    } else{
      link_name <- NA
      if(link_warn)
        warning("link_name and inv_link both supplied; using inv_link because
              link_name not implemented")
    }
  }

  if (is.na(link_name)) {
    grid_gauss <- mvQuad::createNIGrid(dim = 1, type = "GHN", level = num_nodes,
                                       level.trans = FALSE)
    nodes <- as.vector(mvQuad::getNodes(grid_gauss))
    weights <- as.vector(mvQuad::getWeights(grid_gauss))
    W <- matrix(rep(nodes, each = n), nrow = n, ncol = num_nodes,
                 byrow = FALSE)
    W <- W * sigma
    W <- W + eta
    W <- inv_link(W)
    W <- t(weights * t(W))
    eta <- rowSums(W)
  } else if (link_name == "identity") {
    # Do nothing; eta is prediction
  } else if (link_name == "log") {
    eta <- exp(eta + 0.5 * sigma^2)
  } else if (link_name == "sqrt") {
    eta <- eta^2 + sigma^2
  } else {
    warning("Requested link_name not implemented; returning NA")
    eta <- rep(NA, n)
  }

  # Return
  eta
}

#' Calculate predictions for lme4 fitted models
#'
#' Computes marginal predictions (expected values) for models fitted with
#' \code{lme4::glmer} or \code{lme4::lmer}. For linear mixed models (lmerMod),
#' returns the linear predictor. For generalized linear mixed models (glmerMod),
#' integrates over random effects to obtain predictions on the response scale.
#'
#' @param fit A fitted model object from \code{lme4::glmer} or \code{lme4::lmer}
#' @return A vector of predicted (fitted) values. For lmerMod, returns X %*% beta.
#'   For glmerMod, returns marginal expectations on the response scale.
#' @examples
#' \dontrun{
#' library(lme4)
#' fit <- glmer(y ~ x + (1|group), data = mydata, family = poisson())
#' predictions <- pred_glmer(fit)
#' }
#' @export
pred_glmer <- function(fit) {
  if(inherits(fit, "lmerMod")) {
    pred <- as.vector(lme4::getME(fit, "X") %*% lme4::getME(fit, "beta"))
  } else if (inherits(fit, "glmerMod")) {
    pred <- as.vector(lme4::getME(fit, "X") %*% lme4::getME(fit, "beta"))
    # H = Z %*% Lambda, so diag(H H^T) = diag(Z Sigma Z^T) gives Var(z_i^T b)
    # per observation. rowSums on a sparse matrix avoids the dense coercion
    # that apply(., 1, sum) would force.
    H <- lme4::getME(fit, "Z") %*% lme4::getME(fit, "Lambda")
    # unname() so predictions do not inherit observation names
    sigma <- unname(sqrt(Matrix::rowSums(H^2)))
    pred <- pred_base(eta = pred,
                      sigma = sigma,
                      link_name = fit@resp$family$link,
                      inv_link = fit@resp$family$linkinv
    )
  } else{
    stop("fit must be lmerMod or glmerMod")
  }
  pred
}

#' Calculate predictions for nlme fitted models
#'
#' Computes marginal predictions (expected values) for models fitted with
#' \code{nlme::lme} or \code{nlme::nlme}. For linear mixed models (\code{lme}),
#' the marginal prediction is simply the fixed-effects linear predictor
#' X \%*\% beta. For nonlinear mixed models (\code{nlme}), marginal predictions
#' are computed by Monte Carlo integration over the random effects distribution.
#'
#' For nonlinear fits with nested or crossed grouping factors (multiple entries
#' in \code{fit$modelStruct$reStruct}), the per-parameter random-effect
#' contribution at each observation is a sum of independent draws across
#' grouping levels, so the marginal distribution has covariance equal to the
#' sum of per-level covariances (indexed by parameter name). The Monte Carlo
#' draws use this combined covariance.
#'
#' @param fit A fitted model object from \code{nlme::lme} or \code{nlme::nlme}.
#' @param nsim Number of Monte Carlo samples for nonlinear models (default:
#'   1000). Ignored for linear models.
#' @param seed Optional integer seed. If supplied, the RNG state is restored
#'   on exit so the call does not pollute the global stream. Ignored for
#'   linear models.
#' @return A vector of predicted (fitted) values on the response scale.
#' @importFrom stats model.matrix formula
#' @examples
#' \dontrun{
#' library(nlme)
#' # Linear mixed model
#' fit_lme <- lme(distance ~ age, random = ~ 1 | Subject, data = Orthodont)
#' pred_nlme(fit_lme)
#'
#' # Nonlinear mixed model
#' fit_nlme <- nlme(height ~ SSasymp(age, Asym, R0, lrc),
#'                  fixed = Asym + R0 + lrc ~ 1,
#'                  random = Asym ~ 1 | Seed,
#'                  data = Loblolly)
#' pred_nlme(fit_nlme, nsim = 2000, seed = 1)
#' }
#' @export
pred_nlme <- function(fit, nsim = 1000, seed = NULL) {
  if (inherits(fit, "nlme")) {
    # Nonlinear mixed model: Monte Carlo integration over random effects.
    model_rhs <- formula(fit)[[3]]
    beta <- nlme::fixef(fit)
    d <- nlme::getData(fit)
    if (is.null(d)) {
      stop("Could not retrieve original data via nlme::getData(fit). ",
           "Refit with the data argument or save the data alongside the fit.")
    }

    # Combine RE covariances across grouping levels into a single Psi.
    # nlme parameterizes reStruct[[k]] in units of the residual SD, so multiply
    # by sigma^2 to recover the covariance on the response scale. Per-parameter
    # contributions across nested or crossed levels are independent, so the
    # combined marginal covariance is the name-indexed sum of per-level Psi.
    re_struct <- fit$modelStruct$reStruct
    re_names <- unique(unlist(lapply(re_struct,
                                     function(p) colnames(as.matrix(p)))))
    q <- length(re_names)
    Psi <- matrix(0, q, q, dimnames = list(re_names, re_names))
    for (k in seq_along(re_struct)) {
      Pk <- as.matrix(re_struct[[k]]) * fit$sigma^2
      nm <- colnames(Pk)
      Psi[nm, nm] <- Psi[nm, nm] + Pk
    }

    missing_re <- setdiff(re_names, names(beta))
    if (length(missing_re) > 0) {
      stop("Random-effect parameter names not found in fixef(): ",
           paste(missing_re, collapse = ", "))
    }

    L <- t(chol(Psi))
    n <- nrow(d)

    # Restore the global RNG state on exit so a user-supplied seed does not
    # change the stream seen by subsequent code.
    if (!is.null(seed)) {
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
        old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv),
                add = TRUE)
      } else {
        on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
      }
      set.seed(seed)
    }

    # Evaluate inside the formula's environment so user-defined model
    # functions (e.g., custom self-starters) are still in scope.
    enclos <- environment(formula(fit))
    if (is.null(enclos)) enclos <- parent.frame()
    data_list <- as.list(as.data.frame(d))

    pred_sum <- rep(0, n)
    for (s in seq_len(nsim)) {
      b <- as.vector(L %*% stats::rnorm(q))
      params <- beta
      params[re_names] <- params[re_names] + b
      env <- c(data_list, as.list(params))
      pred_sum <- pred_sum + eval(model_rhs, envir = env, enclos = enclos)
    }
    pred_sum / nsim

  } else if (inherits(fit, "lme")) {
    # Linear mixed model: marginal prediction is X %*% beta.
    X <- model.matrix(fit, data = fit$data)
    beta <- nlme::fixef(fit)
    as.vector(X %*% beta)
  } else {
    stop("fit must be of class lme or nlme")
  }
}

#' Per-observation SD of the random-effect contribution to a linear predictor
#'
#' Internal helper shared by the conditional and zero-inflation components of
#' \code{pred_glmmtmb}. Uses \code{reformulas::mkReTrms} only to reconstruct the
#' random-effects design matrix (Zt) and the indexing (Lind) describing where
#' parameters live in a sparse block-diagonal matrix; the nonzero entries of
#' that matrix are then overwritten with glmmTMB's estimated covariance blocks
#' (a VarCorr component), so after forceSymmetric() Psi is the covariance of
#' the stacked random effects U. Given Psi = Cov(U), the per-observation SD is
#' sigma_i = sqrt(Var(z_i^T U)) = sqrt((Z Psi Z^T)_{ii}).
#'
#' @param form Formula whose bars define the random effects (the conditional
#'   formula or the zi formula).
#' @param fr Model frame containing the grouping variables.
#' @param vc Corresponding component of \code{glmmTMB::VarCorr(fit)} (a list of
#'   covariance blocks).
#' @return An n-vector of SDs; zeros if the formula has no random effects.
#' @noRd
re_sd_glmmtmb <- function(form, fr, vc) {
  bars <- reformulas::findbars(form)
  if (is.null(bars)) return(rep(0, nrow(fr)))
  re_terms <- reformulas::mkReTrms(bars = bars,
                             fr = fr,
                             reorder.terms = FALSE # for compatibility w. glmmTMB
                             )

  # Iterate over grouping factors to fill in covariance matrix elements
  theta <- c()
  for (ii in seq_along(vc)) {
    theta <- c(theta, vc[[ii]][lower.tri(vc[[ii]], diag = TRUE)])
  }
  Psi <- re_terms$Lambdat

  # Fill in the upper triangular part of Psi with the extracted elements,
  # then symmetrize
  Psi@x <- theta[re_terms$Lind]
  Psi <- Matrix::forceSymmetric(Psi, uplo = "U")

  # sigma_i^2 = z_i^T Psi z_i with z_i = Zt[, i]. colSums(Zt * (Psi %*% Zt))
  # keeps everything sparse and avoids forming the dense n x n product
  # Z Psi Z^T. unname() so predictions do not inherit observation names.
  unname(sqrt(Matrix::colSums(re_terms$Zt * (Psi %*% re_terms$Zt))))
}

#' Calculate predictions for glmmTMB fitted models
#'
#' Computes marginal predictions (expected values) for models fitted with
#' \code{glmmTMB::glmmTMB}, integrating over the random effects in the
#' conditional model and, when present, in the zero-inflation model.
#'
#' For zero-inflated fits, \code{type = "response"} (the default) returns the
#' marginal mean of the response, (1 - pi_i) * mu_i, where pi_i is the
#' zero-inflation probability marginalized over any zi random effects and mu_i
#' the conditional-component mean marginalized over the conditional random
#' effects; \code{type = "conditional"} returns mu_i alone. For truncated
#' families (\code{truncated_poisson}, \code{truncated_nbinom1},
#' \code{truncated_nbinom2}) the conditional-component mean accounts for the
#' zero truncation, which requires the dispersion parameter, so those families
#' are supported only with constant dispersion (\code{dispformula = ~1}).
#'
#' @param fit A fitted model object from \code{glmmTMB::glmmTMB}
#' @param type \code{"response"} for the marginal mean of the response
#'   (includes the zero-inflation factor), \code{"conditional"} for the
#'   marginal mean of the conditional component only. Identical when the model
#'   has no zero inflation.
#' @param num_nodes Number of Gauss--Hermite nodes used whenever a component
#'   requires numerical integration (default: 15)
#' @return A vector of predicted (fitted) values on the response scale
#' @examples
#' \dontrun{
#' library(glmmTMB)
#' fit <- glmmTMB(y ~ x + (1|group), data = mydata, family = poisson())
#' predictions <- pred_glmmtmb(fit)
#' # Hurdle model: marginal mean including the zero component
#' fit2 <- glmmTMB(y ~ x + (1|group), ziformula = ~ x,
#'                 family = truncated_poisson, data = mydata)
#' predictions2 <- pred_glmmtmb(fit2)
#' }
#' @export
pred_glmmtmb <- function(fit, type = c("response", "conditional"),
                         num_nodes = 15) {
  stopifnot(inherits(fit, "glmmTMB"))
  type <- match.arg(type)
  VC <- glmmTMB::VarCorr(fit)
  fam <- fit$modelInfo$family

  # Conditional component: eta and per-observation RE standard deviation.
  # Use fixef()$cond rather than getME(fit, "beta"): for families with a
  # dispersion parameter (e.g. nbinom1/nbinom2), getME appends it to beta,
  # making the product with X non-conformable.
  eta <- as.vector(glmmTMB::getME(fit, "X") %*% glmmTMB::fixef(fit)$cond)
  sigma <- re_sd_glmmtmb(stats::formula(fit), fit$frame, VC$cond)

  trunc_families <- c("truncated_poisson", "truncated_nbinom1",
                      "truncated_nbinom2")
  if (fam$family %in% trunc_families) {
    # The zero-truncated mean is mu / (1 - P(0)), a per-family transform of
    # the untruncated mean mu = linkinv(eta + u); it has no closed-form
    # normal integral, so always integrate by quadrature. P(0) depends on the
    # dispersion parameter, which must be a scalar for the transform below.
    if (!identical(deparse(fit$modelInfo$allForm$dispformula), "~1")) {
      stop("truncated families are supported only with constant dispersion ",
           "(dispformula = ~1)")
    }
    disp <- glmmTMB::sigma(fit)
    trunc_mean <- switch(fam$family,
      # P(0) = exp(-mu); -expm1() keeps 1 - P(0) accurate for small mu
      truncated_poisson = function(x) {
        mu <- fam$linkinv(x)
        mu / (-expm1(-mu))
      },
      # nbinom1: Var = mu * (1 + phi), i.e. NB with size mu / phi, so
      # P(0) = (1 + phi)^(-mu / phi)
      truncated_nbinom1 = function(x) {
        mu <- fam$linkinv(x)
        mu / (-expm1(-(mu / disp) * log1p(disp)))
      },
      # nbinom2: Var = mu * (1 + mu / theta), size theta, so
      # P(0) = (theta / (theta + mu))^theta
      truncated_nbinom2 = function(x) {
        mu <- fam$linkinv(x)
        mu / (-expm1(disp * (log(disp) - log(disp + mu))))
      })
    pred <- pred_base(eta = eta, sigma = sigma, link_name = NA,
                      num_nodes = num_nodes, inv_link = trunc_mean)
  } else {
    # Supply both link name and inv link; pred_base uses the inv link with
    # quadrature if link_name is not supported analytically
    pred <- pred_base(eta = eta, sigma = sigma,
                      link_name = fam$link,
                      num_nodes = num_nodes,
                      inv_link = fam$linkinv)
  }

  # Zero-inflation factor: E[Y] = (1 - pi) * mu with pi and mu marginalized
  # over their (independent) random effects. glmmTMB's zi link is always
  # logit; with zi random effects the logit-normal mean has no closed form,
  # so integrate by quadrature.
  ziform <- fit$modelInfo$allForm$ziformula
  if (type == "response" && !identical(deparse(ziform), "~0")) {
    eta_zi <- as.vector(glmmTMB::getME(fit, "Xzi") %*% glmmTMB::fixef(fit)$zi)
    sigma_zi <- re_sd_glmmtmb(ziform, fit$frame, VC$zi)
    if (all(sigma_zi == 0)) {
      p_zi <- stats::plogis(eta_zi)
    } else {
      p_zi <- pred_base(eta = eta_zi, sigma = sigma_zi, link_name = NA,
                        num_nodes = num_nodes, inv_link = stats::plogis)
    }
    pred <- (1 - p_zi) * pred
  }

  pred
}
