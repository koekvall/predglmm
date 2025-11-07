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
#' @return A vector of predicted (fitted) values on the natural scale
#' @export
pred_base <- function(eta, sigma = NA, link_name = NA, num_nodes = 15, inv_link = NA)
{
  n <- length(eta)
  if(is.na(sigma)){
    warning("sigma not supplied; setting to zeros (no REs)")
    sigma <- rep(0, n)
  }

  if(is.na(link_name) & is.na(inv_link)){
    stop("Exactly one of link_name and inv_link is needed")
  }

  if(!is.na(link_name) & !is.na(inv_link)){
    link_name <- NA
    warning("link_name and inv_link both supplied; using inv_link and ignoring
            link_name")
  }

  if(is.na(link_name)){
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
  } else if(link_name == "identity"){
    # Do nothing; eta is prediction
  } else if(link_name == "log"){
    eta <- exp(eta + 0.5 * sigma^2)
  } else if(link_name == "sqrt"){
    eta <- eta^2 + sigma^2
  } else{
    warning("Requested link not implemented; returning NA")
    eta <- rep(NA, n)
  }

  # Return
  eta
}

pred_glmer <- function(fit)
{
  X <- lme4::getME(fit, "X")
  Beta <- lme4::getME(fit, "beta")
  pred <- X %*% Beta

  if(class(fit) == "lmerMod"){
    # Do nothing, eta = Xb is prediciction
  } else{ # Nonlinear
    H <- lme4::getME(fit, "Z") %*% lme4::getME(fit, "Lambda")
    sigma <- sqrt(rowSums(H^2))
    pred <- pred_base(eta = pred, sigma = sigma, link = fit@resp$family$link)
  }
  pred
}



pred_glmmtmb <-function(fit)
{
  VC <- glmmTMB::VarCorr(fit)

  # Get lme4 structure
  re_terms <- lme4::mkReTrms(bars = lme4::findbars(formula(fit)),
                             fr = fit$frame,
                             reorder.terms = F # for compatibility w. glmmTMB
                             )
  # Iterate over grouping factors to fill in covariance matrix elements
  theta <- c()
  for(ii in 1:length(VC$cond)){
    theta <- c(theta, VC$cond[[ii]][lower.tri(VC$cond[[ii]], diag = T)])
  }
  Psi <- re_terms$Lambdat

  # Fill in the upper triangular part of Psi with the extracted elements
  Psi@x <- theta[re_terms$Lind]

  # Make matrix symmetric
  Psi <- Matrix::forceSymmetric(Psi, uplo = "U")

  sigma <- rowSums(crossprod(re_terms$Zt, Psi) * t(re_terms$Zt))

  pred_base(eta = getME(fit, "X") %*% getME(fit, "beta"),
            sigma = sigma,
            link = fit$modelInfo$family$link)

}
