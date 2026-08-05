#' Parametric bootstrap confidence intervals for marginal predictions
#'
#' Computes confidence intervals for marginal predictions, or for functions
#' of them, by a parametric bootstrap. Responses are simulated from the
#' fitted model, the model is refit to each simulated response, marginal
#' predictions are computed on each refit with \code{\link{pred_glmmtmb}} or
#' \code{\link{pred_glmer}}, and the statistic is applied to the model data
#' augmented with those predictions. Intervals are formed from the bootstrap
#' draws of the statistic.
#'
#' @param fit A \code{glmmTMB} object from \code{glmmTMB::glmmTMB}, or a
#'   \code{merMod} object from \code{lme4::glmer} or \code{lme4::lmer}.
#' @param statistic Function of one argument, the data the model was fit to
#'   with an added column \code{pred} holding the marginal predictions
#'   (\code{NA} for rows dropped by the fit's NA handling), returning a
#'   numeric vector, named if longer than one. Columns of the data not
#'   appearing in the model formula are available. If \code{NULL} (default),
#'   the predictions for the fitted observations.
#' @param n_boot Positive integer. Number of bootstrap replicates. Default
#'   is 1000.
#' @param level Numeric confidence level in (0, 1). Default is \code{0.95}.
#' @param cores Positive integer. Number of workers. Default is 1 (serial).
#'   Parallel replicates use forking (\code{parallel::mclapply}) on
#'   unix-alikes and a PSOCK cluster on Windows.
#' @param seed Integer or \code{NULL}. If supplied, the RNG kind is set to
#'   L'Ecuyer-CMRG, making results reproducible for fixed
#'   \code{(seed, cores)}, and the caller's RNG state is restored on exit.
#' @param type Character, either \code{"percentile"} (default) or
#'   \code{"normal"}. Percentile intervals are empirical quantiles of the
#'   draws; normal intervals are the estimate plus and minus a standard
#'   normal quantile times the standard deviation of the draws.
#'
#' @return An object of class \code{boot_pred}: a list with components
#'   \code{estimate} (the statistic on the original fit), \code{draws}
#'   (matrix with one column per kept replicate), \code{ci} (matrix with
#'   columns \code{lower} and \code{upper}), \code{n_boot}, \code{n_fail},
#'   \code{n_boundary} (\code{NA} for glmmTMB fits made with
#'   \code{se = FALSE}, see Details), \code{level}, \code{type}, and
#'   \code{call}. Supports \code{print}.
#'
#' @details
#' A replicate is dropped only if the refit throws an error or its optimizer
#' reports non-convergence; the number dropped is returned as \code{n_fail}.
#' Refits with a variance component estimated at zero (singular fits for
#' lme4, flagged by a non-positive-definite Hessian for glmmTMB) are kept
#' and counted in \code{n_boundary}: such refits occur exactly when a
#' variance component is small, so dropping them would bias the intervals.
#' A large \code{n_boundary} indicates the data are compatible with a zero
#' variance component, a regime where the parametric bootstrap is
#' unreliable.
#'
#' The statistic receives the data frame passed as the model call's
#' \code{data} argument, so variables not appearing in the model formula are
#' available. For \code{merMod} fits with no recoverable data frame
#' (variables taken from the environment, or a \code{subset} argument), the
#' model frame is used instead and only formula variables are available.
#' Rows dropped by NA handling have \code{pred = NA}; a statistic that
#' aggregates over such rows should allow for this, for example with
#' \code{na.rm = TRUE}.
#'
#' glmmTMB fits are refit by re-evaluating the model call, so the data must
#' be recoverable from the call, the response must be a single column of the
#' data (for binomial responses, use proportions with \code{weights}, not
#' \code{cbind()}), and fits made with \code{subset} are not supported. For
#' fits made with \code{se = FALSE}, boundary status is unavailable and
#' \code{n_boundary} is \code{NA}.
#'
#' The intervals inherit the fitted model's assumptions, since responses are
#' simulated from it. If the model was selected using the same data, the
#' intervals condition on that selection and understate uncertainty.
#' @examples
#' \donttest{
#' set.seed(1)
#' d <- data.frame(
#'   x = factor(rep(c("a", "b"), times = 40)),
#'   group = factor(rep(1:8, each = 10))
#' )
#' d$y <- rpois(80, exp(0.5 + 0.4 * (d$x == "b") +
#'                        rep(rnorm(8, sd = 0.6), each = 10)))
#' fit <- lme4::glmer(y ~ x + (1 | group), data = d, family = poisson())
#'
#' # Percentile CIs for each observation's marginal prediction
#' boot_pred(fit, n_boot = 100, seed = 1)
#'
#' # CI for the difference in mean predictions between levels of x
#' contrast <- function(data) {
#'   m <- tapply(data$pred, data$x, mean)
#'   c("a-b" = unname(m["a"] - m["b"]))
#' }
#' boot_pred(fit, statistic = contrast, n_boot = 100, seed = 1)
#' }
#' @export
boot_pred <- function(fit, statistic = NULL, n_boot = 1000, level = 0.95,
                      cores = 1, seed = NULL,
                      type = c("percentile", "normal")) {
  type <- match.arg(type)
  stopifnot(is.null(statistic) || is.function(statistic),
            is.numeric(n_boot), length(n_boot) == 1, n_boot >= 1,
            is.numeric(level), length(level) == 1, level > 0, level < 1,
            is.numeric(cores), length(cores) == 1, cores >= 1)
  n_boot <- as.integer(n_boot)
  cl_call <- match.call()

  # Per-class ingredients: the data the statistic sees, the prediction
  # and refit functions, and the two screens. "Failed" refits are dropped;
  # "boundary" refits are kept and counted (see Details).
  if (inherits(fit, "glmmTMB")) {
    pred_fun <- pred_glmmtmb
    failed <- function(f) !isTRUE(f$fit$convergence == 0)
    # A fit made with se = FALSE has no sdreport, so the Hessian-based
    # boundary flag cannot be assessed; report NA rather than a bogus count
    if (is.null(fit$sdr)) {
      boundary <- function(f) NA
    } else {
      boundary <- function(f) !isTRUE(f$sdr$pdHess)
    }
    worker_pkg <- "glmmTMB"

    # glmmTMB has no structure-reusing refit (glmmTMB::refit re-evaluates the
    # original call and looks the data up in the caller's frame, which fails
    # inside package code), so recover the data once here and re-evaluate the
    # call per replicate with the simulated response substituted. Symbols in
    # the call (formula, family, ...) resolve in the formula's environment,
    # where the original call was evaluated.
    cc <- stats::getCall(fit)
    if (!is.null(cc$subset)) {
      stop("fits made with a 'subset' argument are not supported; subset ",
           "the data before fitting")
    }
    # Namespace the call head so refits also work on PSOCK workers, where
    # the backend package is loaded but not attached
    cc[[1]] <- quote(glmmTMB::glmmTMB)
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
    # simulate() returns one value per row used in the fit; rows dropped by
    # the fit's NA handling get NA responses in the refit data, so each
    # refit drops exactly the same rows
    used_rows <- setdiff(seq_len(nrow(dat)),
                         as.integer(attr(fit$frame, "na.action")))
    if (length(used_rows) != stats::nobs(fit)) {
      stop("the model data could not be aligned with the fitted ",
           "observations; refit so that every row of the data enters the ",
           "model frame or is dropped by NA handling")
    }
    # The statistic sees the recovered data, not the reduced model frame,
    # so columns not appearing in the formula remain available
    stat_data <- dat
    refit_fun <- function(y) {
      # simulate.glmmTMB returns a (successes, failures) matrix for
      # binomial-type responses fitted as proportions; convert back so the
      # response type matches the call, which still carries any weights
      if (is.matrix(y)) {
        n_trials <- rowSums(y)
        y <- ifelse(n_trials > 0, y[, 1] / n_trials, 0)
      }
      y_col <- rep(NA_real_, nrow(dat))
      y_col[used_rows] <- y
      dat_boot <- dat
      dat_boot[[resp]] <- y_col
      cc_boot <- cc
      cc_boot$data <- quote(.predglmm_boot_data)
      eval(cc_boot, envir = list(.predglmm_boot_data = dat_boot),
           enclos = form_env)
    }
  } else if (inherits(fit, "merMod")) {
    pred_fun <- pred_glmer
    refit_fun <- function(y) lme4::refit(fit, newresp = y)
    failed <- function(f) !isTRUE(f@optinfo$conv$opt == 0)
    boundary <- function(f) lme4::isSingular(f)
    worker_pkg <- "lme4"

    # lme4::refit does not need the original data, but the statistic should
    # still see it so that columns not appearing in the formula remain
    # available. Recover it from the call; if the fit was made without a
    # data argument, or the data cannot be aligned with the fitted
    # observations (e.g. a subset argument), fall back to the model frame,
    # which holds the formula variables only.
    form_env <- environment(stats::formula(fit))
    if (is.null(form_env)) form_env <- parent.frame()
    dat <- tryCatch(eval(stats::getCall(fit)$data, envir = form_env),
                    error = function(e) NULL)
    if (is.data.frame(dat)) {
      used_rows <- setdiff(seq_len(nrow(dat)),
                           as.integer(attr(fit@frame, "na.action")))
    }
    if (is.data.frame(dat) && length(used_rows) == stats::nobs(fit)) {
      stat_data <- dat
    } else {
      stat_data <- fit@frame
      used_rows <- seq_len(nrow(stat_data))
    }
  } else {
    stop("fit must be a glmmTMB or merMod (lme4) fit")
  }

  if ("pred" %in% names(stat_data)) {
    stop("the model data already has a column named 'pred'; rename that ",
         "variable before fitting")
  }

  # The pred column is NA for rows the fit dropped, so a statistic never
  # sees stale or misaligned values; the default statistic returns the
  # predictions for the rows actually used, aligning its output with the
  # fitted observations
  if (is.null(statistic)) statistic <- function(data) data$pred[used_rows]

  add_pred <- function(pred_vals) {
    p <- rep(NA_real_, nrow(stat_data))
    p[used_rows] <- pred_vals
    stat_data$pred <- p
    stat_data
  }

  estimate <- statistic(add_pred(pred_fun(fit)))
  if (!is.numeric(estimate) || length(estimate) < 1) {
    stop("statistic must return a numeric vector")
  }
  q <- length(estimate)

  # Restore the global RNG state on exit so a user-supplied seed does not
  # change the stream seen by subsequent code.
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
    list(stat = statistic(add_pred(pred_fun(fit_boot))),
         boundary = boundary(fit_boot))
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
  bnd <- vapply(res, function(r) as.logical(r$boundary), logical(1))
  n_boundary <- if (anyNA(bnd)) NA_integer_ else sum(bnd)
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

#' Print method for bootstrap confidence intervals
#'
#' @param x A \code{boot_pred} object from \code{\link{boot_pred}}.
#' @param digits Number of digits displayed. Default 4.
#' @param max_rows Maximum number of rows displayed. Default 10.
#' @param ... Unused.
#' @return \code{x}, invisibly.
#' @method print boot_pred
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
  bnd_txt <- if (is.na(x$n_boundary)) {
    "boundary status unknown"
  } else {
    paste0(x$n_boundary, " on the boundary")
  }
  cat(x$n_boot, " replicates: ", ncol(x$draws), " kept (", bnd_txt,
      "), ", x$n_fail, " dropped (refit error or ",
      "non-convergence)\n", sep = "")
  invisible(x)
}
