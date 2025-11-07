#' Calculate predictions (marginal expectations)
#'
#' @param eta An n-vector of linear predictor values (eta = X %*% beta)
#' @param sigma An n-vector of standard deviations for the random effects.
#'   If NA (default), assumes no random effects (sigma = 0).
#' @param link_name Character string specifying the link function. Currently supports
#'   "identity", "log", or "sqrt". Use NA if providing a custom inv_link function.
#' @param num_nodes Number of nodes for Gaussian quadrature used to calculate predictions
#'   when using a custom inv_link function (default: 15)
#' @param inv_link Custom inverse link function. If provided, numerical integration
#'   via Gaussian quadrature is used. Use NA if providing link_name instead.
#' @return A vector of predicted (fitted) values on the response scale
#' @export
pred_base <- function(eta, sigma = NA, link_name = NA, num_nodes = 15,
                      inv_link = NA) {

  supported_links <- c("identity", "log", "sqrt")

  n <- length(eta)
  if (any(is.na(sigma))) {
    warning("NA values for sigma supplied; setting ALL to zeros (no REs)")
    sigma <- rep(0, n)
  }

  if (any(is.na(link_name)) && any(is.na(inv_link))) {
    stop("Exactly one of link_name and inv_link is needed")
  }

  if (!is.na(link_name) && !is.na(inv_link)) {
    if(link_name %in% supported_links) {
      inv_link <- NA
      warning("link_name and inv_link both supplied; trying to use using link_name")
    } else{
      link_name <- NA
      warning("link_name and inv_link both supplied; using inv_link because
              link_name not implemented")
    }
  }

  if (length(link_name) == 1 && is.na(link_name)) {
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
pred_glmer <- function(fit) {
  pred <- lme4::getME(fit, "X") %*% lme4::getME(fit, "beta")

  if (inherits(fit, "lmerMod")) {
    # Do nothing, eta = Xb is prediction
  } else { # Nonlinear
    H <- lme4::getME(fit, "Z") %*% lme4::getME(fit, "Lambda")
    sigma <- sqrt(rowSums(H^2))
    pred <- pred_base(eta = pred,
                      sigma = sigma,
                      link_name = fit@resp$family$link,
                      inv_link = fit@resp$family$invlink
    )
  }
  pred
}


#' Calculate predictions for glmmTMB fitted models
#'
#' Computes marginal predictions (expected values) for models fitted with
#' \code{glmmTMB::glmmTMB}. Integrates over random effects to obtain predictions
#' on the response scale. Uses \code{lme4::mkReTrms} to reconstruct the random
#' effects structure and manually fills the sparse precision matrix from variance
#' components.
#'
#' @param fit A fitted model object from \code{glmmTMB::glmmTMB}
#' @return A vector of predicted (fitted) values on the response scale
#' @examples
#' \dontrun{
#' library(glmmTMB)
#' fit <- glmmTMB(y ~ x + (1|group), data = mydata, family = poisson())
#' predictions <- pred_glmmtmb(fit)
#' }
pred_glmmtmb <- function(fit) {
  VC <- glmmTMB::VarCorr(fit)

  # Get lme4 structure
  re_terms <- lme4::mkReTrms(bars = lme4::findbars(formula(fit)),
                             fr = fit$frame,
                             reorder.terms = FALSE # for compatibility w. glmmTMB
                             )
  # Iterate over grouping factors to fill in covariance matrix elements
  theta <- c()
  for (ii in seq_along(VC$cond)) {
    theta <- c(theta, VC$cond[[ii]][lower.tri(VC$cond[[ii]], diag = TRUE)])
  }
  Psi <- re_terms$Lambdat

  # Fill in the upper triangular part of Psi with the extracted elements
  Psi@x <- theta[re_terms$Lind]

  # Make matrix symmetric
  Psi <- Matrix::forceSymmetric(Psi, uplo = "U")

  sigma <- sqrt(rowSums(crossprod(re_terms$Zt, Psi) * t(re_terms$Zt)))

  # Supply both link name and inv link; will use inv link if link_name
  # not supported
  pred_base(eta = glmmTMB::getME(fit, "X") %*% glmmTMB::getME(fit, "beta"),
            sigma = sigma,
            link_name = fit$modelInfo$family$link,
            inv_link = fit$modelInfo$family$linkinv)

}
