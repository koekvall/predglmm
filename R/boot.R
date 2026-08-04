#' Parametric bootstrap confidence intervals for marginal predictions
#'
#' Simulates responses from a fitted model, refits the model to each simulated
#' response, computes marginal predictions on the refit (\code{pred_glmmtmb}
#' or \code{pred_glmer}, with their default arguments), and applies a
#' user-supplied statistic to the model frame augmented with those
#' predictions. Confidence intervals are formed from the bootstrap draws of
#' the statistic. For lme4 fits, refitting uses \code{lme4::refit}, which
#' reuses the fitted model structure and needs no access to the original
#' data. For glmmTMB fits, the original call is re-evaluated with the
#' simulated response substituted into the data, so the data must be
#' recoverable from the model call and the response must be a single column
#' of it (in particular, \code{cbind()} binomial responses are not
#' supported).
#'
#' A replicate is dropped only when refitting fails outright: the refit
#' throws an error or its optimizer reports non-convergence. Refits on the
#' boundary of the parameter space (a variance component estimated at zero;
#' flagged by \code{lme4::isSingular} for lme4 and by a non-positive-definite
#' Hessian for glmmTMB) are kept: zero is a legitimate value of the
#' estimator, and boundary refits arise exactly when a variance component is
#' small, so discarding them would compute the intervals from a selected
#' subset and bias them. The two counts are returned as \code{n_fail}
#' (dropped) and \code{n_boundary} (kept). A large \code{n_boundary} is not a
#' failure of the bootstrap but a sign the data are consistent with a zero
#' variance component; interval accuracy is reduced in that regime because
#' the parametric bootstrap is unreliable for parameters near a boundary.
#'
#' Two further cautions. First, if the model was selected using the same
#' data, the bootstrap conditions on that selection and understates
#' uncertainty: the fixed model is refit to every replicate, so variability
#' from the selection step is not propagated. Second, this is a parametric
#' bootstrap; responses are simulated from the fitted model, so the intervals
#' inherit all of the model's assumptions (distribution, link,
#' random-effects structure) and are not robust to their failure.
#'
#' @param fit A fitted model object from \code{glmmTMB::glmmTMB} or from
#'   \code{lme4::glmer} / \code{lme4::lmer}.
#' @param statistic Function taking the model frame with an added column
#'   \code{pred} (the marginal predictions) and returning a numeric vector,
#'   named if there is more than one element. Default \code{NULL} returns the
#'   prediction vector itself. The statistic sees the model frame (the
#'   variables as used in the formula), not the original data set.
#' @param n_boot Number of bootstrap replicates (default: 1000).
#' @param level Confidence level (default: 0.95).
#' @param cores Number of workers (default: 1, serial). On unix-alikes,
#'   parallel replicates use forking (\code{parallel::mclapply}); no sockets
#'   are opened, so no firewall prompts or blocked connections. On Windows a
#'   PSOCK cluster is used and stopped on exit.
#' @param seed Optional integer seed. If supplied, the RNG kind is switched
#'   to L'Ecuyer-CMRG so results are reproducible for a fixed
#'   \code{(seed, cores)}, and the caller's RNG state is restored on exit.
#' @param type \code{"percentile"} for quantile intervals from the draws,
#'   \code{"normal"} for \code{estimate +/- z * sd(draws)}.
#' @return An object of class \code{"boot_pred"}: a list with components
#'   \code{estimate} (the statistic on the original fit), \code{draws}
#'   (matrix, one column per kept replicate), \code{ci} (matrix with columns
#'   \code{lower} and \code{upper}), \code{n_boot}, \code{n_fail} (replicates
#'   dropped because the refit errored or did not converge),
#'   \code{n_boundary} (kept replicates whose refit was on the boundary),
#'   \code{level}, \code{type}, and \code{call}.
#' @examples
#' \dontrun{
#' library(glmmTMB)
#' fit <- glmmTMB(y ~ group + (1 | id), data = d, family = poisson())
#'
#' # CI for the difference in cell means of the predictions between groups
#' contrast <- function(data) {
#'   m <- tapply(data$pred, data$group, mean)
#'   c("A-B" = unname(m["A"] - m["B"]))
#' }
#' b <- boot_pred(fit, statistic = contrast, n_boot = 1000, seed = 1)
#' print(b)
#' }
#' @export
boot_pred <- function(fit, statistic = NULL, n_boot = 1000, level = 0.95,
                      cores = 1, seed = NULL, type = c("percentile", "normal")) {
  type <- match.arg(type)
  stopifnot(is.null(statistic) || is.function(statistic),
            is.numeric(n_boot), length(n_boot) == 1, n_boot >= 1,
            is.numeric(level), length(level) == 1, level > 0, level < 1,
            is.numeric(cores), length(cores) == 1, cores >= 1)
  n_boot <- as.integer(n_boot)
  cl_call <- match.call()

  # Per-class ingredients: the model frame the statistic sees, the prediction
  # and refit functions, and the two screens. "Failed" refits are dropped;
  # "boundary" refits are kept and counted (see Details).
  if (inherits(fit, "glmmTMB")) {
    frame <- fit$frame
    pred_fun <- pred_glmmtmb
    failed <- function(f) !isTRUE(f$fit$convergence == 0)
    boundary <- function(f) !isTRUE(f$sdr$pdHess)
    worker_pkg <- "glmmTMB"

    # glmmTMB has no structure-reusing refit (glmmTMB::refit re-evaluates the
    # original call and looks the data up in the caller's frame, which fails
    # inside package code), so recover the data once here and re-evaluate the
    # call per replicate with the simulated response substituted. Symbols in
    # the call (formula, family, ...) resolve in the formula's environment,
    # where the original call was evaluated.
    cc <- stats::getCall(fit)
    form_env <- environment(stats::formula(fit))
    if (is.null(form_env)) form_env <- parent.frame()
    dat <- tryCatch(eval(cc$data, envir = form_env), error = function(e) NULL)
    if (!is.data.frame(dat)) {
      stop("Could not retrieve the model data via getCall(fit)$data. ",
           "Refit with a data frame that is still available in the ",
           "environment of the model formula.")
    }
    resp <- names(fit$modelInfo$respCol)
    if (length(resp) != 1 || !(resp %in% names(dat))) {
      stop("the response must be a single column of the model data")
    }
    refit_fun <- function(y) {
      dat_boot <- dat
      dat_boot[[resp]] <- y
      cc_boot <- cc
      cc_boot$data <- quote(.predglmm_boot_data)
      eval(cc_boot, envir = list(.predglmm_boot_data = dat_boot),
           enclos = form_env)
    }
  } else if (inherits(fit, "merMod")) {
    frame <- fit@frame
    pred_fun <- pred_glmer
    refit_fun <- function(y) lme4::refit(fit, newresp = y)
    failed <- function(f) !isTRUE(f@optinfo$conv$opt == 0)
    boundary <- function(f) lme4::isSingular(f)
    worker_pkg <- "lme4"
  } else {
    stop("fit must be a glmmTMB or merMod (lme4) fit")
  }

  if ("pred" %in% names(frame)) {
    stop("the model frame already has a column named 'pred'; rename that ",
         "variable in the model")
  }

  if (is.null(statistic)) statistic <- function(data) data$pred

  frame_est <- frame
  frame_est$pred <- pred_fun(fit)
  estimate <- statistic(frame_est)
  if (!is.numeric(estimate) || length(estimate) < 1) {
    stop("statistic must return a numeric vector")
  }
  q <- length(estimate)

  # Restore the global RNG state on exit so a user-supplied seed does not
  # change the stream seen by subsequent code (same pattern as pred_nlme).
  # Restoring .Random.seed also restores the RNG kind, which is switched to
  # L'Ecuyer-CMRG below so that parallel streams are reproducible.
  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv),
              add = TRUE)
    } else {
      # Removing .Random.seed makes R reinitialize with the default kinds
      on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
    }
    RNGkind("L'Ecuyer-CMRG")
    set.seed(seed)
  }

  # One replicate returns the q statistic values plus a boundary flag, or
  # NULL if the refit fails; fitted objects and simulation matrices are never
  # accumulated, so memory use is O(q x n_boot) regardless of model size.
  # Warnings and messages from the refits (convergence chatter, singular-fit
  # notices) are suppressed because their information is aggregated into
  # n_fail and n_boundary; a thousand repetitions of them would drown the
  # console.
  one_rep <- function(ii) {
    y_boot <- stats::simulate(fit)[[1]]
    fit_boot <- tryCatch(
      suppressWarnings(suppressMessages(refit_fun(y_boot))),
      error = function(e) NULL)
    if (is.null(fit_boot) || failed(fit_boot)) return(NULL)
    frame_boot <- frame
    frame_boot$pred <- pred_fun(fit_boot)
    list(stat = statistic(frame_boot), boundary = boundary(fit_boot))
  }

  if (cores > 1 && .Platform$OS.type != "windows") {
    # Forking, not sockets: PSOCK workers can be blocked by firewalls and
    # would need the package installed; forked workers share this process
    res <- parallel::mclapply(seq_len(n_boot), one_rep, mc.cores = cores,
                              mc.set.seed = TRUE)
    err <- vapply(res, inherits, logical(1), what = "try-error")
    if (any(err)) {
      stop("a bootstrap replicate failed: ",
           conditionMessage(attr(res[[which(err)[1]]], "condition")))
    }
  } else if (cores > 1) {
    cl <- parallel::makePSOCKcluster(cores)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    if (!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
    # one_rep's environment chain reaches the predglmm namespace, which is
    # shipped by reference; the backend package must be loaded explicitly for
    # simulate() and refit() dispatch on the workers
    parallel::clusterCall(cl, function(p) requireNamespace(p, quietly = TRUE),
                          worker_pkg)
    res <- parallel::parLapply(cl, seq_len(n_boot), one_rep)
  } else {
    res <- lapply(seq_len(n_boot), one_rep)
  }

  keep <- !vapply(res, is.null, logical(1))
  n_fail <- sum(!keep)
  if (n_fail == n_boot) {
    stop("all ", n_boot, " bootstrap refits failed (error or non-convergence)")
  }
  res <- res[keep]
  n_boundary <- sum(vapply(res, function(r) isTRUE(r$boundary), logical(1)))
  # vapply enforces that every kept replicate has length q; matrix() keeps
  # the q x m shape also when q = 1, where vapply returns a plain vector
  draws <- matrix(vapply(res, function(r) as.vector(r$stat), numeric(q)),
                  nrow = q, dimnames = list(names(estimate), NULL))

  alpha <- (1 - level) / 2
  if (type == "percentile") {
    ci <- t(apply(draws, 1, stats::quantile, probs = c(alpha, 1 - alpha),
                  names = FALSE))
  } else {
    z <- stats::qnorm(1 - alpha)
    sds <- apply(draws, 1, stats::sd)
    ci <- cbind(estimate - z * sds, estimate + z * sds)
  }
  dimnames(ci) <- list(names(estimate), c("lower", "upper"))

  structure(list(estimate = estimate,
                 draws = draws,
                 ci = ci,
                 n_boot = n_boot,
                 n_fail = n_fail,
                 n_boundary = n_boundary,
                 level = level,
                 type = type,
                 call = cl_call),
            class = "boot_pred")
}

#' Print a boot_pred object
#'
#' @param x An object of class \code{"boot_pred"}.
#' @param digits Number of digits to display (default: 4).
#' @param max_rows Maximum number of statistic elements to display
#'   (default: 10).
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @export
print.boot_pred <- function(x, digits = 4, max_rows = 10, ...) {
  cat("Parametric bootstrap CIs (", x$type, ", level = ", x$level, ")\n",
      sep = "")
  tab <- cbind(estimate = x$estimate, x$ci)
  n_show <- min(nrow(tab), max_rows)
  print(round(tab[seq_len(n_show), , drop = FALSE], digits))
  if (nrow(tab) > n_show) {
    cat("... ", nrow(tab) - n_show, " more rows\n", sep = "")
  }
  cat(x$n_boot, " replicates: ", ncol(x$draws), " kept (", x$n_boundary,
      " on the boundary), ", x$n_fail, " dropped (refit error or ",
      "non-convergence)\n", sep = "")
  invisible(x)
}
