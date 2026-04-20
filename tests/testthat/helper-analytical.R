# Helper: compute the per-observation random-effect variance
#   sigma_i^2 = z_i^T Psi z_i
# independently of the implementation in pred_glmer / pred_glmmtmb.
#
# For each grouping factor k with within-group design Z_k and covariance Psi_k,
# the contribution to sigma_i^2 is
#   z_{k,i}^T Psi_k z_{k,i}
# where z_{k,i} is the row of Z_k for observation i. Crucially, this does not
# depend on the grouping levels themselves — only on the within-group design
# (the LHS of the bar) and the per-group covariance from VarCorr. Summing over
# bars handles crossed / nested / multiple-factor structures uniformly.

sigma_sq_from_varcorr <- function(formula_obj, frame_data, vc_list) {
  re_terms <- lme4::findbars(formula_obj)
  stopifnot(length(re_terms) == length(vc_list))
  sig_sq <- numeric(nrow(frame_data))
  for (k in seq_along(re_terms)) {
    bar <- re_terms[[k]]
    re_form <- stats::as.formula(call("~", bar[[2]]))
    X_re <- stats::model.matrix(re_form, frame_data)
    Psi_g <- as.matrix(vc_list[[k]])
    attr(Psi_g, "stddev") <- NULL
    attr(Psi_g, "correlation") <- NULL
    sig_sq <- sig_sq + rowSums((X_re %*% Psi_g) * X_re)
  }
  sig_sq
}

sigma_sq_glmer <- function(fit) {
  sigma_sq_from_varcorr(stats::formula(fit), fit@frame, lme4::VarCorr(fit))
}

sigma_sq_glmmtmb <- function(fit) {
  sigma_sq_from_varcorr(stats::formula(fit), fit$frame,
                        glmmTMB::VarCorr(fit)$cond)
}

eta_glmer <- function(fit) {
  as.vector(lme4::getME(fit, "X") %*% lme4::getME(fit, "beta"))
}

eta_glmmtmb <- function(fit) {
  as.vector(glmmTMB::getME(fit, "X") %*% glmmTMB::getME(fit, "beta"))
}

# Closed-form marginal means under Gaussian random effects.
analytical_pred <- function(eta, sigma_sq, link) {
  switch(link,
         log      = exp(eta + sigma_sq / 2),
         sqrt     = eta^2 + sigma_sq,
         identity = eta,
         stop("no closed form for link ", link))
}
