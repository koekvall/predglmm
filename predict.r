#' Calculate predictions (marginal expectations)
#'
#' @param X An n x p matrix of predictors
#' @param Beta A p-vector of regression coefficients
#' @param sigma An n-vector of standard deviations for the latent variables,
#' @param type The inverse link function, currently one of
#' "Linear", "Exponential", or "Logistic"
#' @param num_nodes Number of nodes for Gaussian quadrature used to calculate predictions
#' @return A vector of predicted (fitted) values
#' @export
pred_base <- function(X, Beta, sigma, type = "Exponential",
                      num_nodes = 15)
{
  # Argument checking
  stopifnot(is.matrix(X),
            is.numeric(Beta), is.atomic(Beta),
            is.numeric(sigma), is.atomic(sigma), all(sigma >= 0),
            is.numeric(num_nodes), is.atomic(num_nodes),
            length(num_nodes) == 1, floor(num_nodes) == num_nodes,
            num_nodes >= 1)

  # Define constants
  p <- ncol(X)
  n <- nrow(X)
  stopifnot(floor(n) == n, length(Beta) == p, length(sigma) == n)
  Xb <-as.vector(X %*% Beta)

  if(type == "Linear"){
    return(Xb)
  }

  if(type == "Exponential"){
    Xb <- exp(Xb + 0.5 * sigma^2)
    return(Xb)
  }

  # Standard normal quadrature
  grid_gauss <- mvQuad::createNIGrid(dim = 1, type = "GHN", level = num_nodes,
                                     level.trans = FALSE)
  nodes <- as.vector(mvQuad::getNodes(grid_gauss))
  weights <- as.vector(mvQuad::getWeights(grid_gauss))
  W0 <- matrix(rep(nodes, each = n), nrow = n, ncol = num_nodes,
               byrow = FALSE)

  if(type == "Logistic"){
    W <- W0 * sigma
    W <- W + Xb
    W <- 1 / (1 + exp(-W))
    W <- t(weights * t(W))
    Xb <- rowSums(W)
  } else{
    warning("Requested type not implemented")
    Xb <- rep(NA, n)
  }

  return(Xb)
}


