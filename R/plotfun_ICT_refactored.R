#' Plot a bayesICT object
#'
#' General plot method for bayesICT objects.  The method dispatches to one of
#' the specialized plotting functions while keeping the dispatch explicit and
#' checked.
#'
#' @param x A bayesICT object.
#' @param type Plot type: "prevalence", "trace", or "coef".
#' @param main Optional plot title.
#' @param ... Additional arguments passed to the selected plot method.
#'
#' @return The plot object when the selected backend returns one; otherwise the
#'   input object invisibly.
#' @method plot bayesICT
#' @export
plot.bayesICT <- function(x, type = c("prevalence", "trace", "coef"),
                          main = "", ...) {
  ## Check the requested plot type early so misspellings fail with a
  ## standard match.arg() message instead of failing later in dispatch.
  type <- match.arg(type)

  ## Use explicit switch-based dispatch rather than constructing function
  ## names as strings.  This keeps the S3 method easier to audit.
  plotter <- switch(type,
    prevalence = plot_prevalence.bayesICT,
    trace = plot_trace.bayesICT,
    coef = plot_coef.bayesICT
  )

  ## All specialized plotting functions use the same argument names for
  ## the fitted object and title, so the general method can simply pass
  ## through remaining graphical arguments via ....
  plotter(object = x, main = main, ...)
}

#' Logistic function
#'
#' @param x A scalar or vector.
#'
#' @return The inverse-logit transform of x.
#' @export
## Keep the historical exported helper, but delegate to R's stable
## implementation of the inverse-logit transform.
logistic <- function(x) stats::plogis(x)

#' Plot posterior sensitive-item prevalence
#'
#' Produces the base-R prevalence plot used by the original bayesICT plotting
#' code.  When covariates are present and only_prev is FALSE, a second panel
#' displays posterior summaries for the sensitive-item regression coefficients.
#'
#' @param object A bayesICT object.
#' @param covariate_names Optional labels for the sensitive-item coefficients.
#' @param only_prev If TRUE, plot only the prevalence posterior.
#' @param xlim Optional x-axis limits for prevalence.
#' @param inverse If TRUE, plot 1 - prevalence.
#' @param digit Number of digits for the displayed posterior variance.
#' @param main Optional plot title.
#' @param group Integer posterior prevalence group to plot when pas is vector-valued.
#' @param ci Length-two vector of posterior interval probabilities.
#' @param ... Additional graphical parameters passed to base plotting functions.
#'
#' @return The input object, invisibly.
#' @export
plot_prevalence.bayesICT <- function(object, covariate_names = NULL,
                                     only_prev = FALSE, xlim = NULL,
                                     inverse = FALSE, digit = 3,
                                     main = "", group = 1L,
                                     ci = c(0.025, 0.975), ...) {
  ## The prevalence plot uses base graphics and temporarily changes
  ## graphical parameters.  Save and restore the user's graphics state.
  old_par <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(old_par), add = TRUE)

  ## With an intercept-only sensitive-item regression there is no
  ## coefficient panel to draw, so reduce automatically to a one-panel plot.
  if (!is.null(object$npars) && object$npars <= 1L) {
    only_prev <- TRUE
    warning("No covariates found. Plotting the prevalence only.", call. = FALSE)
  }

  ## Extract posterior draws for the population prevalence parameter.
  ## pas can be vector-valued when auxiliary groups are present.
  pas <- bayesICT_plot_parameter_matrix(object, "pas")
  if (isTRUE(inverse)) {
    pas <- 1 - pas
  }

  ## Select the requested prevalence group after checking that it exists.
  group <- bayesICT_plot_validate_group(group, ncol(pas))
  pas_group <- pas[, group]

  if (is.null(xlim)) {
    xlim <- c(0, 1)
  }

  ## Posterior interval and variance displayed on the histogram.
  ci <- sort(ci)
  q <- stats::quantile(pas_group, probs = ci, na.rm = TRUE)
  posterior_var <- stats::var(pas_group, na.rm = TRUE)

  ## Compute the histogram first to choose a y-axis range that leaves room
  ## for the credible-interval annotation.
  hist_info <- graphics::hist(pas_group, probability = TRUE, plot = FALSE)
  ymax <- max(hist_info$density, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) {
    ymax <- 1
  }

  ## Draw the prevalence histogram and overlay posterior density, auxiliary
  ## reference value, and credible interval.  It is defined locally so it can
  ## be reused in one-panel and two-panel layouts.
  prevalence_panel <- function() {
    color <- grDevices::rgb(166, 165, 166, maxColorValue = 255, alpha = 100)

    graphics::plot(NA,
      xlim = xlim, ylim = c(0, ymax * 1.7), axes = FALSE,
      xlab = "", ylab = "", main = main, ...)

    graphics::hist(pas_group,
      probability = TRUE, add = TRUE, col = color,
      main = main, xlab = "Proportion affirmative on the sensitive item", ...)

    if (length(unique(stats::na.omit(pas_group))) > 1L) {
      graphics::lines(stats::density(pas_group, na.rm = TRUE), col = color, lwd = 4)
    }

    ## For auxiliary-prior fits, show the external estimate h as a vertical
    ## dashed reference line when it is available for this group.
    aux_h <- bayesICT_plot_aux_h(object, group)
    if (!is.null(aux_h)) {
      graphics::abline(v = aux_h, lty = 2)
    }

    graphics::axis(1, c(xlim[1], xlim[2]), labels = c(xlim[1], xlim[2]), las = 1)
    graphics::segments(q[1], ymax * 1.10, q[2], ymax * 1.10, lwd = 4, col = color)
    graphics::segments(q[1], ymax * 1.05, q[1], ymax * 1.15, lwd = 4, col = color)
    graphics::segments(q[2], ymax * 1.05, q[2], ymax * 1.15, lwd = 4, col = color)

    graphics::text(mean(q), ymax * 1.25, "95% Credible Interval", col = color)
    graphics::text(q[2] + 0.10, ymax * 1.10,
      paste("Var.:", round(posterior_var, digit)))
    graphics::mtext("Prevalence", side = 1, line = 1, at = mean(xlim))
  }

  ## If only the prevalence panel is requested, draw it and exit before
  ## setting up a two-row plotting layout.
  if (isTRUE(only_prev)) {
    prevalence_panel()
    return(invisible(object))
  }

  ## Otherwise stack the prevalence panel above the sensitive-item
  ## coefficient summaries.
  graphics::par(mfrow = c(2, 1), mai = c(0.6, 0.1, 0.1, 0.1), xpd = TRUE)
  prevalence_panel()
  bayesICT_plot_delta_panel(object, covariate_names = covariate_names,
                            main = main, ...)

  invisible(object)
}

#' Trace plots for a bayesICT object
#'
#' @param object A bayesICT object.
#' @param main Optional plot title.
#' @param pars Parameters to plot.  By default, all monitored parameters except
#'   log_lik and lp__ are plotted.
#' @param label If TRUE, try to replace Stan parameter names in the returned
#'   ggplot object by the same human-readable labels used by coef.bayesICT().
#' @param ... Additional arguments passed to rstan::traceplot().
#'
#' @return A ggplot object.
#' @export
plot_trace.bayesICT <- function(object, main = "", pars = NULL,
                                label = TRUE, ...) {
  ## Trace plots require ggplot2 for plotting and rstan for extracting the
  ## unpermuted iteration-by-chain arrays.
  bayesICT_plot_require_ggplot2()
  bayesICT_plot_require_rstan()

  ## Expand top-level requests such as "delta" into the concrete posterior
  ## draw columns, for example delta[1], delta[2], ....
  pars <- bayesICT_plot_default_pars(object, pars)
  par_info <- bayesICT_plot_expand_parameters(object, pars)

  ## Use unpermuted draws so each chain remains identifiable in the trace.
  arr <- rstan::extract(object$stanfit, pars = unique(par_info$top),
                        permuted = FALSE, inc_warmup = FALSE)
  dat <- bayesICT_plot_trace_data(arr, par_info, object, label = label)

  ## One facet per parameter, with chains overlaid within each facet.
  p <- ggplot2::ggplot(dat, ggplot2::aes(x = iteration, y = value, group = chain)) +
    ggplot2::geom_line(linewidth = 0.25, alpha = 0.75) +
    ggplot2::facet_wrap(stats::as.formula("~ parameter"), scales = "free_y") +
    ggplot2::labs(x = "Iteration", y = "Value") +
    ggplot2::theme_bw()

  p <- bayesICT_plot_add_title(p, main)
  p
}

#' Coefficient plots for a bayesICT object
#'
#' @param object A bayesICT object.
#' @param main Optional plot title.
#' @param pars Parameters to plot.  By default, all monitored parameters except
#'   log_lik and lp__ are plotted.
#' @param label If TRUE, try to replace Stan parameter names in the returned
#'   ggplot object by the same human-readable labels used by coef.bayesICT().
#' @param ... Additional arguments passed to rstan::plot().
#'
#' @return A ggplot object.
#' @export
plot_coef.bayesICT <- function(object, main = "", pars = NULL,
                               label = TRUE, ci_level = 0.8,
                               outer_level = 0.95, ...) {
  bayesICT_plot_require_ggplot2()

  ## Build interval summaries directly from object$sampledf rather than
  ## relying on rstan::plot(); this avoids backend-specific grid/PDF issues
  ## and gives us full control over labels.
  pars <- bayesICT_plot_default_pars(object, pars)
  par_info <- bayesICT_plot_expand_parameters(object, pars)
  dat <- bayesICT_plot_interval_data(object, par_info, label = label,
                                     ci_level = ci_level,
                                     outer_level = outer_level)

  ## Draw outer and inner posterior intervals horizontally, centered at the
  ## posterior median, with a zero reference line for regression parameters.
  p <- ggplot2::ggplot(dat, ggplot2::aes(y = parameter, x = median)) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = outer_low, xmax = outer_high),
                            linewidth = 0.4, alpha = 0.55) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = inner_low, xmax = inner_high),
                            linewidth = 1.2, alpha = 0.85) +
    ggplot2::geom_point(size = 1.7) +
    ggplot2::geom_vline(xintercept = 0, linetype = 2, linewidth = 0.25) +
    ggplot2::labs(x = "Posterior estimate", y = NULL) +
    ggplot2::theme_bw()

  p <- bayesICT_plot_add_title(p, main)
  p
}

## -------------------------------------------------------------------------
## Internal helpers
## -------------------------------------------------------------------------


#' Require ggplot2 for plotting
#'
#' Stops with a clear message if ggplot2 is unavailable.
#'
#' @return TRUE invisibly.
bayesICT_plot_require_ggplot2 <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package is required for this plot.", call. = FALSE)
  }
  invisible(TRUE)
}

#' Posterior draw names available for plotting
#'
#' Removes internal Stan quantities that are not normally useful for plots.
#'
#' @param object A bayesICT object.
#' @return Character vector of available posterior-draw column names.
bayesICT_plot_available_draw_names <- function(object) {
  nms <- names(object$sampledf)
  if (is.null(nms)) {
    character(0L)
  } else {
    bayesICT_plot_drop_internal_pars(nms)
  }
}

#' Extract the top-level Stan parameter name
#'
#' Converts names such as delta[1] or delta.1 to delta.
#'
#' @param x Character vector of Stan-style names.
#' @return Character vector of top-level names.
bayesICT_plot_top_name <- function(x) {
  sub("(\\[.*$|\\..*$)", "", x)
}

#' Normalize Stan parameter names
#'
#' Converts data-frame-safe names such as delta.1 back to Stan-like names
#' such as delta[1].
#'
#' @param x Character vector of parameter names.
#' @return Normalized character vector.
bayesICT_plot_normalize_name <- function(x) {
  x <- as.character(x)
  x <- sub("\\.([0-9]+)$", "[\\1]", x)
  x
}

#' Expand requested parameters to posterior-draw columns
#'
#' Matches user requests for top-level parameters against concrete columns in
#' object$sampledf.  For example, a request for delta matches delta[1],
#' delta[2], and so on.
#'
#' @param object A bayesICT object.
#' @param pars Character vector of requested parameters.
#' @return Data frame with sampledf column names, normalized Stan names, and
#'   top-level Stan names.
bayesICT_plot_expand_parameters <- function(object, pars) {
  if (is.null(object$sampledf)) {
    stop("object$sampledf is missing.")
  }

  ## Work only with posterior-draw columns that are actually present in the
  ## saved sample data frame.
  available <- bayesICT_plot_available_draw_names(object)
  if (length(available) == 0L) {
    stop("No posterior draws were found in object$sampledf.")
  }

  ## If pars is NULL, request all top-level parameters available after
  ## dropping internal quantities.
  requested <- if (is.null(pars)) {
    unique(bayesICT_plot_top_name(available))
  } else {
    as.character(pars)
  }

  cols <- character(0L)
  missing <- character(0L)

  ## Each requested parameter may be a scalar column or a vector parameter
  ## represented by several indexed columns.
  for (p in requested) {
    p <- as.character(p)
    exact <- available == p
    indexed <- startsWith(available, paste0(p, "[")) |
      startsWith(available, paste0(p, "."))
    hit <- available[exact | indexed]

    if (length(hit) == 0L) {
      missing <- c(missing, p)
    } else {
      cols <- c(cols, hit)
    }
  }

  cols <- unique(cols)

  ## Missing parameters are warned about rather than treated as fatal unless
  ## none of the requested parameters can be plotted.
  if (length(missing) > 0L) {
    warning("Ignoring parameter(s) not present in the fitted object: ",
            paste(unique(missing), collapse = ", "), call. = FALSE)
  }

  if (length(cols) == 0L) {
    stop("None of the requested parameters are present in the fitted object.")
  }

  ## Store both the concrete Stan name and the top-level name needed by
  ## rstan::extract().
  normalized <- bayesICT_plot_normalize_name(cols)
  top <- bayesICT_plot_top_name(normalized)

  data.frame(
    column = cols,
    stan_name = normalized,
    top = top,
    stringsAsFactors = FALSE
  )
}

#' Posterior interval summaries for coefficient plots
#'
#' Computes medians and inner/outer credible intervals from posterior draws.
#'
#' @param object A bayesICT object.
#' @param par_info Output from bayesICT_plot_expand_parameters().
#' @param label If TRUE, use coefficient labels when available.
#' @param ci_level Probability content of the inner interval.
#' @param outer_level Probability content of the outer interval.
#' @return Data frame used by plot_coef.bayesICT().
bayesICT_plot_interval_data <- function(object, par_info, label = TRUE,
                                        ci_level = 0.8,
                                        outer_level = 0.95) {
  if (!is.numeric(ci_level) || length(ci_level) != 1L ||
      ci_level <= 0 || ci_level >= 1) {
    stop("ci_level must be a number strictly between 0 and 1.")
  }
  if (!is.numeric(outer_level) || length(outer_level) != 1L ||
      outer_level <= 0 || outer_level >= 1) {
    stop("outer_level must be a number strictly between 0 and 1.")
  }

  ## The coefficient plot uses the already saved posterior draws in
  ## object$sampledf, avoiding a second extraction from the stanfit object.
  draws <- object$sampledf[, par_info$column, drop = FALSE]

  ## Convert interval widths into equal-tailed quantile probabilities.
  inner_p <- c((1 - ci_level) / 2, 1 - (1 - ci_level) / 2)
  outer_p <- c((1 - outer_level) / 2, 1 - (1 - outer_level) / 2)

  ## Summarize each requested parameter separately.
  stat <- lapply(seq_along(par_info$column), function(i) {
    x <- draws[[i]]
    c(
      mean = mean(x, na.rm = TRUE),
      median = stats::median(x, na.rm = TRUE),
      inner_low = stats::quantile(x, probs = inner_p[1], na.rm = TRUE, names = FALSE),
      inner_high = stats::quantile(x, probs = inner_p[2], na.rm = TRUE, names = FALSE),
      outer_low = stats::quantile(x, probs = outer_p[1], na.rm = TRUE, names = FALSE),
      outer_high = stats::quantile(x, probs = outer_p[2], na.rm = TRUE, names = FALSE)
    )
  })

  stat <- as.data.frame(do.call(rbind, stat))
  ## Reuse the print-file labeling helper when that file has been sourced.
  lab <- par_info$stan_name
  if (isTRUE(label) && exists("bayesICT_coef_labels", mode = "function")) {
    lab <- bayesICT_coef_labels(par_info$stan_name, object)
  }

  stat$parameter <- factor(lab, levels = rev(unique(lab)))
  stat$stan_name <- par_info$stan_name
  stat
}

#' Convert unpermuted Stan draws to trace-plot data
#'
#' Takes the iteration-by-chain-by-parameter array from rstan::extract() and
#' returns a long data frame suitable for ggplot2 faceting.
#'
#' @param arr Unpermuted posterior draw array from rstan::extract().
#' @param par_info Output from bayesICT_plot_expand_parameters().
#' @param object A bayesICT object.
#' @param label If TRUE, use coefficient labels when available.
#' @return Long-format trace-plot data frame.
bayesICT_plot_trace_data <- function(arr, par_info, object, label = TRUE) {
  ## The third dimension of the unpermuted rstan array contains parameter
  ## names; without those names we cannot align draws to requested parameters.
  dn <- dimnames(arr)
  if (is.null(dn) || length(dn) < 3L || is.null(dn[[3L]])) {
    stop("Could not extract named posterior draws from object$stanfit.")
  }

  available <- bayesICT_plot_normalize_name(dn[[3L]])
  want <- par_info$stan_name
  ## Match requested normalized Stan names to the parameter dimension of the
  ## extracted array.
  keep <- match(want, available)
  ok <- !is.na(keep)
  if (!any(ok)) {
    stop("None of the requested parameters could be extracted for trace plotting.")
  }

  keep <- keep[ok]
  want <- want[ok]

  lab <- want
  if (isTRUE(label) && exists("bayesICT_coef_labels", mode = "function")) {
    lab <- bayesICT_coef_labels(want, object)
  }

  ## Some extracted Stan arrays carry dimnames.  When a scalar column with a
  ## name is combined with longer unnamed vectors, data.frame() warns that row
  ## names from a short variable were discarded.  Strip names/dimensions from
  ## every component before assembling the plotting data.
  lab <- unname(as.character(lab))
  want <- unname(as.character(want))

  ## Reshape each parameter's iteration-by-chain matrix into long form.
  dims <- dim(arr)
  out <- vector("list", length(keep))
  for (j in seq_along(keep)) {
    mat <- arr[, , keep[j], drop = FALSE]
    val <- unname(as.numeric(mat))
    iter <- unname(rep(seq_len(dims[1L]), times = dims[2L]))
    chain <- unname(rep(seq_len(dims[2L]), each = dims[1L]))
    out[[j]] <- data.frame(
      iteration = iter,
      chain = factor(chain),
      value = val,
      parameter = rep.int(lab[j], length(val)),
      stan_name = rep.int(want[j], length(val)),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  }

  dat <- do.call(rbind, out)
  dat$parameter <- factor(dat$parameter, levels = unique(dat$parameter))
  dat
}

#' Extract posterior draws as a matrix
#'
#' Uses the shared print-function helper when available, with a local fallback
#' for standalone use of the plotting file.
#'
#' @param object A bayesICT object.
#' @param parameter Top-level Stan parameter name.
#' @return Matrix of posterior draws.
bayesICT_plot_parameter_matrix <- function(object, parameter) {
  ## Prefer the shared extraction helper from printfun_ICT_refactored.R so
  ## printing, prediction, and plotting agree on parameter-name handling.
  if (exists("bayesICT_parameter_draws", mode = "function")) {
    return(bayesICT_parameter_draws(object$sampledf, parameter))
  }

  if (is.null(object$sampledf)) {
    stop("object$sampledf is missing.")
  }

  ## Local fallback: match scalar, bracket-indexed, and data-frame-safe
  ## dot-indexed versions of the requested parameter.
  pattern <- paste0("^", parameter, "(\\[|\\.|$)")
  loc <- grep(pattern, names(object$sampledf))
  if (length(loc) == 0L) {
    stop("Could not find posterior draws for parameter '", parameter, "'.")
  }

  as.matrix(object$sampledf[, loc, drop = FALSE])
}

#' Validate a prevalence-group index
#'
#' @param group Requested group index.
#' @param n_groups Number of available prevalence groups.
#' @return Integer group index.
bayesICT_plot_validate_group <- function(group, n_groups) {
  if (length(group) != 1L || is.na(group)) {
    stop("group must be a single integer.")
  }

  group <- as.integer(group)
  if (group < 1L || group > n_groups) {
    stop("group must be between 1 and ", n_groups, ".")
  }

  group
}

#' Auxiliary reference value for prevalence plots
#'
#' @param object A bayesICT object.
#' @param group Requested auxiliary/prevalence group.
#' @return The auxiliary h value, or NULL if unavailable.
bayesICT_plot_aux_h <- function(object, group) {
  ## Auxiliary h values exist only for auxiliary-prior models.
  if (is.null(object$prior) || !(object$prior %in% c("auxiliary", "double_list_aux"))) {
    return(NULL)
  }

  h <- object$aux_info$h
  if (is.null(h) || length(h) < group || is.na(h[group])) {
    return(NULL)
  }

  h[group]
}

#' Draw the sensitive-item coefficient panel
#'
#' @param object A bayesICT object.
#' @param covariate_names Optional coefficient labels.
#' @param main Optional plot title.
#' @param ... Additional base-graphics arguments.
#' @return NULL invisibly.
bayesICT_plot_delta_panel <- function(object, covariate_names = NULL,
                                      main = "", ...) {
  ## The lower panel displays posterior summaries for the sensitive-item
  ## regression coefficients delta.
  modout <- bayesICT_plot_summary_rows(object, "delta")
  color <- grDevices::rgb(166, 165, 166, maxColorValue = 255, alpha = 100)

  if (nrow(modout) == 0L) {
    graphics::plot.new()
    graphics::text(0.5, 0.5, "No sensitive-item coefficients to plot")
    return(invisible(NULL))
  }

  ## Prefer named columns, but retain the original rstan::summary() positional
  ## fallback for older objects.
  keep_cols <- c("mean", "2.5%", "97.5%")
  if (!all(keep_cols %in% colnames(modout))) {
    ## Fall back to the column positions used by rstan::summary().
    modout <- modout[, c(1, 4, 8), drop = FALSE]
    colnames(modout) <- keep_cols
  } else {
    modout <- modout[, keep_cols, drop = FALSE]
  }

  ## Choose an x-axis range that covers all posterior means and intervals.
  xmin <- min(modout, na.rm = TRUE)
  xmax <- max(modout, na.rm = TRUE)
  if (!is.finite(xmin) || !is.finite(xmax) || xmin == xmax) {
    xmin <- xmin - 1
    xmax <- xmax + 1
  }

  k <- nrow(modout)
  graphics::plot(NA,
    xlim = c(xmin, xmax), ylim = c(0, k + 1), axes = FALSE,
    xlab = "", ylab = "", main = main, ...)

  ## Draw each coefficient as a horizontal 95% interval with a point at the
  ## posterior mean.
  for (i in seq_len(k)) {
    graphics::segments(modout[i, "2.5%"], i,
      modout[i, "97.5%"], i, lwd = 4, col = color)
    graphics::points(modout[i, "mean"], i, col = color, pch = 16)
  }

  ## Use caller-provided labels when given; otherwise infer labels from the
  ## fitted object.
  labels <- covariate_names
  if (is.null(labels)) {
    labels <- bayesICT_plot_delta_labels(object, rownames(modout))
  }

  labels <- rep_len(labels, k)
  graphics::text(modout[, "mean"], seq_len(k) + 0.1, labels)
  graphics::axis(side = 1, tick = TRUE)
  graphics::mtext("Coefficient Estimates", side = 1, line = 2)

  invisible(NULL)
}

#' Extract summary rows for a top-level Stan parameter
#'
#' @param object A bayesICT object.
#' @param parameter Top-level Stan parameter name.
#' @return Matrix/data frame subset of object$summaryout.
bayesICT_plot_summary_rows <- function(object, parameter) {
  if (is.null(object$summaryout)) {
    stop("object$summaryout is missing.")
  }

  ## The summary matrix uses Stan parameter names as row names.
  nms <- rownames(object$summaryout)
  if (is.null(nms)) {
    stop("object$summaryout must have row names.")
  }

  pattern <- paste0("^", parameter, "(\\[|$)")
  object$summaryout[grep(pattern, nms), , drop = FALSE]
}

#' Labels for sensitive-item coefficients
#'
#' Uses the shared coefficient-label helper when present, with a local fallback
#' based on object$xnames.
#'
#' @param object A bayesICT object.
#' @param par_names Stan parameter row names.
#' @return Character vector of labels.
bayesICT_plot_delta_labels <- function(object, par_names) {
  ## Share coefficient labels with the print methods whenever possible.
  if (exists("bayesICT_coef_labels", mode = "function")) {
    return(bayesICT_coef_labels(par_names, object))
  }

  xnames <- object$xnames
  if (is.null(xnames) || length(xnames) == 0L) {
    return(par_names)
  }

  idx <- vapply(par_names, function(x) {
    inside <- regmatches(x, regexpr("\\[.*\\]", x))
    if (length(inside) == 0L || is.na(inside) || inside == "") {
      return(NA_integer_)
    }
    as.integer(sub("^\\[([0-9]+).*$", "\\1", inside))
  }, integer(1L))

  ifelse(!is.na(idx) & idx >= 1L & idx <= length(xnames), xnames[idx], par_names)
}

#' Default parameters for trace and coefficient plots
#'
#' @param object A bayesICT object.
#' @param pars Optional user-supplied parameter vector.
#' @return Character vector of parameters to plot.
bayesICT_plot_default_pars <- function(object, pars = NULL) {
  if (!is.null(pars)) {
    return(pars)
  }

  ## Prefer the monitored-parameter list saved by rstan, falling back to the
  ## stanfit names if the internal slot is unavailable.
  pars <- tryCatch(object$stanfit@sim$pars_oi, error = function(e) NULL)
  if (is.null(pars)) {
    pars <- names(object$stanfit)
  }

  bayesICT_plot_drop_internal_pars(pars)
}

#' Remove internal Stan parameters from plotting requests
#'
#' Drops log_lik and lp__ so default plots focus on model parameters.
#'
#' @param pars Character vector of Stan parameter names.
#' @return Filtered character vector.
bayesICT_plot_drop_internal_pars <- function(pars) {
  ## Avoid bracket-sensitive regular expressions here.  The Stan names are
  ## either top-level names such as "log_lik" and "lp__", or indexed names
  ## such as "log_lik[1]".
  pars <- as.character(pars)
  is_log_lik <- pars == "log_lik" | startsWith(pars, "log_lik[")
  is_lp <- pars == "lp__"
  pars[!(is_log_lik | is_lp)]
}

#' Build a Stan-name to human-label map
#'
#' This helper is retained for compatibility with older relabeling paths.
#'
#' @param pars Character vector of parameter names.
#' @param object A bayesICT object.
#' @return Named character vector mapping Stan names to display labels.
bayesICT_plot_label_map <- function(pars, object) {
  ## Start from all Stan names so indexed parameters can be mapped to labels.
  all_names <- names(object$stanfit)
  if (is.null(all_names)) {
    return(stats::setNames(pars, pars))
  }

  keep <- bayesICT_plot_drop_internal_pars(all_names)
  labels <- if (exists("bayesICT_coef_labels", mode = "function")) {
    bayesICT_coef_labels(keep, object)
  } else {
    keep
  }

  label_map <- stats::setNames(labels, keep)

  ## If the user supplied top-level parameter names such as "delta", include
  ## identity mappings too.  rstan often expands these to indexed names in the
  ## returned ggplot data.
  for (p in pars) {
    if (!p %in% names(label_map)) {
      label_map[p] <- p
    }
  }

  label_map
}

#' Relabel parameter columns in an existing ggplot object
#'
#' This is useful for ggplot objects returned by external plotting backends.
#'
#' @param p A ggplot object.
#' @param label_map Named character vector.
#' @return The relabeled plot object.
bayesICT_plot_relabel_ggplot <- function(p, label_map) {
  if (is.null(p) || is.null(label_map) || length(label_map) == 0L) {
    return(p)
  }

  if (!inherits(p, "ggplot")) {
    return(p)
  }

  ## Replace parameter names only when they appear in the map; unmatched names
  ## are left unchanged.
  relabel_vec <- function(x) {
    x_chr <- as.character(x)
    y <- unname(label_map[x_chr])
    y[is.na(y)] <- x_chr[is.na(y)]
    if (is.factor(x)) {
      return(factor(y, levels = unique(y)))
    }
    y
  }

  ## Different plotting backends use different column names for facets or
  ## parameters, so check the common alternatives.
  candidate_cols <- c("parameter", "variable", "Parameter", "pars")
  for (nm in intersect(candidate_cols, names(p$data))) {
    p$data[[nm]] <- relabel_vec(p$data[[nm]])
  }

  p
}

#' Add a title to a ggplot object when requested
#'
#' @param p Plot object.
#' @param main Optional title.
#' @return Plot object, with title added when applicable.
bayesICT_plot_add_title <- function(p, main) {
  if (!is.character(main) || length(main) == 0L || identical(main, "")) {
    return(p)
  }

  if (requireNamespace("ggplot2", quietly = TRUE) && inherits(p, "ggplot")) {
    return(p + ggplot2::ggtitle(main))
  }

  p
}

#' Require rstan for trace extraction
#'
#' Stops with a clear message if rstan is unavailable.
#'
#' @return TRUE invisibly.
bayesICT_plot_require_rstan <- function() {
  if (!requireNamespace("rstan", quietly = TRUE)) {
    stop("The rstan package is required for this plot.", call. = FALSE)
  }
  invisible(TRUE)
}
