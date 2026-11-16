#' Extract coefficients from a bayesICT object
#'
#' Create a table of coefficient summaries from a bayesICT object.
#'
#' @param object A bayesICT object returned by bayesICT().
#' @param ... Further arguments; currently unused.
#'
#' @return A coefficient-summary matrix with human-readable row names.
#' @method coef bayesICT
#' @export
  # Start from the stored Stan summary table, then replace raw Stan
  # parameter names with labels suitable for user-facing output.
coef.bayesICT <- coefficients.bayesICT <- function(object, ...) {
  coefmat <- bayesICT_coef_matrix(object)
  # Replace raw Stan names such as delta[1] with interpretable labels.
  rownames(coefmat) <- bayesICT_coef_labels(rownames(coefmat), object)
  coefmat
}

#' Print a bayesICT object
#'
#' @param x A bayesICT object.
#' @param type Character string giving the print type: "text", "mcmc", or "coef".
#' @param ... Additional arguments passed to the selected print method.
#'
#' @return The input object, invisibly.
#' @export
  # Normalize the requested print mode and dispatch to the corresponding
  # small method.  This avoids constructing function names as strings.
print.bayesICT <- function(x, type = c("text", "mcmc", "coef"), ...) {
  type <- match.arg(type)
  # Select the corresponding print routine.
  printer <- switch(type,
    text = print_text.bayesICT,
    mcmc = print_mcmc.bayesICT,
    coef = print_coef.bayesICT
  )
  printer(object = x, ...)
  # Preserve S3 print convention by returning the object invisibly below.
  invisible(x)
}

#' Summarize a bayesICT object
#'
#' @param object A bayesICT object.
#' @param ... Additional arguments passed to print_text.bayesICT().
#'
#' @return The input object, invisibly.
#' @method summary bayesICT
#' @export
  # The summary method uses the same textual display as the default print
  # method, but returns the fitted object invisibly like standard summary methods.
summary.bayesICT <- function(object, ...) {
  print_text.bayesICT(object = object, ...)
  invisible(object)
}

#' Predicted prevalence from a bayesICT object
#'
#' @param object A bayesICT object.
#' @param type Posterior summary statistic: "median" or "mean".
#' @param list For double-list models, predict from the first list, second list,
#'   or their average. The default uses the average when both lists have
#'   separate prevalence parameters.
#' @param ... Further arguments; currently unused.
#'
#' @return A numeric vector of predicted sensitive-item prevalence values.
#'   For `list_variable_model = "diff_means"`, the values are treated-minus-
#'   control expected list-count differences under the fitted ordinal models.
#' @method predict bayesICT
#' @export
  # Prediction summarizes posterior draws of the sensitive-item prevalence
  # for each row of the model matrix.
predict.bayesICT <- function(object, type = c("median", "mean"),
                             list = c("auto", "first", "second", "average"),
                             ...) {
  type_missing <- missing(type)
  # Record whether the user supplied a prediction summary type; double-list
  # models historically defaulted to posterior means.
  type <- match.arg(type)
  # Select the corresponding print routine.
  list <- match.arg(list)

  has_second_delta <- bayesICT_has_parameter_draws(object$sampledf, "delta2")
  # Some double-list models share delta across lists, while others estimate delta2.
  has_second_X <- !is.null(object$X2)
  is_diff_means <- identical(object$list_variable_model, "diff_means")

  p1 <- if (is_diff_means) {
    bayesICT_diff_means_prevalence(
      object$X, object$sampledf, object$J,
      treated_coef = "delta", control_coef = "psi0",
      treated_cuts = "psiMT", control_cuts = "psiM"
    )
  } else {
    bayesICT_linear_prevalence(object$X, object$sampledf, "delta")
  }

  if (!bayesICT_has_double_list(object)) {
    # Single-list models predict from the primary sensitive-item regression.
    if (!identical(list, "auto") && !identical(list, "first")) {
      stop("Only list = 'auto' or list = 'first' is valid for single-list models.")
    }
    return(bayesICT_summarize_rows(p1, type))
  }

  ## Preserve the original convention that double-list predictions are reported
  ## as posterior means unless the user explicitly asks otherwise.  Some
  ## double-list Stan models estimate a common sensitive-item regression
  ## coefficient named `delta`; others estimate both `delta` and `delta2`.
  ## The original refactor incorrectly assumed that `delta2` was always
  ## present.
  if (type_missing) {
    type <- "mean"
  }

  if (identical(list, "first") ||
      (identical(list, "auto") && !(has_second_delta && has_second_X))) {
    draws <- p1
  } else {
    if (!(has_second_delta && has_second_X)) {
      stop("Second-list prediction requested, but posterior draws for 'delta2' or object$X2 are missing.")
    }

    p2 <- if (is_diff_means) {
      bayesICT_diff_means_prevalence(
        object$X2, object$sampledf, object$J,
        treated_coef = "delta2", control_coef = "psi02",
        treated_cuts = "psiMT2", control_cuts = "psiM2"
      )
    } else {
      bayesICT_linear_prevalence(object$X2, object$sampledf, "delta2")
    }
    # Second-list prevalence draws are only available in models with delta2.

    if (identical(list, "second")) {
      draws <- p2
    } else {
      draws <- 0.5 * p1 + 0.5 * p2
    }
  }

  bayesICT_summarize_rows(draws, type)
}

#' Print the main results from a bayesICT object
#'
#' @param object A bayesICT object.
#' @param digits Number of digits to display.
#'
#' @return The input object, invisibly.
#' @export
  # Print the model header, original call, MCMC run size, and coefficient table.
print_text.bayesICT <- function(object, digits = 5) {
  cat("Bayesian sensitive item model: ", object$type, "\n",
      "With number of control items = ", object$J, "\n", sep = "")

  cat("\nCall:\n",
      paste(deparse(object$Call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("MCMC run for", object$nsim, "iterations, with",
      bayesICT_warmup_used(object), "used. \n\n")

  cat("Coefficients:\n")
  print(round(coef(object), digits))
  # coef() performs the parameter filtering and row-label construction.
  cat("\n")

  invisible(object)
}

#' Print convergence diagnostics from a bayesICT object
#'
#' @param object A bayesICT object.
#' @param ... Additional arguments passed to print().
#'
#' @return The input object, invisibly.
#' @export
  # Delegate to the stanfit print method so that Stan diagnostics and
  # posterior summaries are shown in their usual format.
print_mcmc.bayesICT <- function(object, ...) {
  print(object$stanfit, ...)
  invisible(object)
}

#' Print coefficients of a bayesICT object
#'
#' @param object A bayesICT object.
#' @param digits Number of digits to display.
#'
#' @return The input object, invisibly.
#' @export
  # Print only the coefficient summary, using the same labels as coef().
print_coef.bayesICT <- function(object, digits = 5) {
  print(round(coef(object), digits))
  # coef() performs the parameter filtering and row-label construction.
  invisible(object)
}

## --------------------------------------------------------------------------
## Internal helpers
## --------------------------------------------------------------------------

#' Extract the non-internal coefficient-summary matrix
#'
#' Removes log-likelihood and Stan bookkeeping parameters from object$summaryout.
#'
#' @param object A bayesICT object.
#' @return A numeric matrix of posterior summaries for model parameters.
bayesICT_coef_matrix <- function(object) {
  if (is.null(object$summaryout)) {
    # summaryout is the central posterior-summary table stored by bayesICT().
    stop("object$summaryout is missing.")
  }

  coefmat <- object$summaryout
  if (is.null(rownames(coefmat))) {
    stop("object$summaryout must have row names containing Stan parameter names.")
  }

  keep <- !grepl("^log_lik(\\[|$)", rownames(coefmat)) &
  # Drop quantities used for model comparison or Stan bookkeeping.
    !grepl("^lp__$", rownames(coefmat))
  coefmat[keep, , drop = FALSE]
}

#' Convert Stan parameter names to bayesICT coefficient labels
#'
#' @param par_names Character vector of raw Stan parameter names.
#' @param object A bayesICT object used to recover covariate names and model type.
#' @return A character vector of display labels.
bayesICT_coef_labels <- function(par_names, object) {
  vapply(par_names, bayesICT_coef_label_one, character(1), object = object)
}

#' Convert one Stan parameter name to one display label
#'
#' @param par_name A single raw Stan parameter name.
#' @param object A bayesICT object.
#' @return A single character label.
bayesICT_coef_label_one <- function(par_name, object) {
  parsed <- bayesICT_parse_stan_name(par_name)
  # Split the raw name into a base parameter and optional Stan indices.
  base <- parsed$base
  index <- parsed$index

  switch(base,
  # Map each recognized Stan parameter family to the label used in summaries.
    psi0 = paste0("Control items:", bayesICT_xname(object, index[1])),
    delta = paste0("Sensitive item:", bayesICT_xname(object, index[1])),
    psi02 = paste0("Control items 2:", bayesICT_xname(object, index[1])),
    delta2 = paste0("Sensitive item 2:", bayesICT_xname(object, index[1])),
    psiM = bayesICT_control_level_label(object, index[1], list = 1),
    psiM2 = bayesICT_control_level_label(object, index[1], list = 2),
    rho0 = "Rho: intra-resp. corr. on controls",
    bias = bayesICT_bias_label(object, index[1], list = 1),
    bias2 = bayesICT_bias_label(object, index[1], list = 2),
    pas = bayesICT_pas_label(object, index[1]),
    par_name
  )
}

#' Parse an indexed Stan parameter name
#'
#' Splits strings such as "delta[2]" into base name "delta" and integer index 2.
#'
#' @param x A Stan parameter name.
#' @return A list with components `base` and `index`.
bayesICT_parse_stan_name <- function(x) {
  base <- sub("\\[.*$", "", x)
  # The part before '[' is the Stan parameter base name.
  bracket <- regmatches(x, regexpr("\\[.*\\]", x))
  # The bracketed part, when present, contains one or more integer indices.

  if (length(bracket) == 0L || is.na(bracket) || bracket == "") {
    return(list(base = base, index = integer(0)))
  }

  index <- substring(bracket, 2L, nchar(bracket) - 1L)
  # Remove the square brackets before converting comma-separated indices.
  index <- as.integer(strsplit(index, ",", fixed = TRUE)[[1]])

  list(base = base, index = index)
}

#' Safely recover a covariate name
#'
#' @param object A bayesICT object.
#' @param i Integer index of the covariate.
#' @return The covariate name, or the index as text if unavailable.
bayesICT_xname <- function(object, i) {
  xnames <- object$xnames
  # xnames are the column labels from the model matrix.
  if (length(xnames) == 0L || is.na(i) || i < 1L || i > length(xnames)) {
    return(if (is.na(i)) "" else as.character(i))
  }
  xnames[i]
}

#' Build a label for multinomial control-item level parameters
#'
#' @param object A bayesICT object.
#' @param i Control level index.
#' @param list List indicator for double-list models.
#' @return A display label.
bayesICT_control_level_label <- function(object, i, list = 1L) {
  if (is.na(i)) {
    i <- ""
  }

  if (bayesICT_has_double_list(object)) {
    # Distinguish first- and second-list control-level parameters when needed.
    paste0("Control levels for ", if (list == 1L) "1st" else "2nd", " controls:", i)
  } else {
    paste0("Control levels:", i)
  }
}

#' Build a label for auxiliary bias parameters
#'
#' @param object A bayesICT object.
#' @param i Bias-parameter index.
#' @param list List indicator for double-list models.
#' @return A display label.
bayesICT_bias_label <- function(object, i, list = 1L) {
  if (identical(object$bias, "separate") || list == 2L) {
    paste0("Bias sensitive item ", list)
  } else {
    "Bias of sensitive item"
  }
}

#' Build a label for sensitive-item prevalence parameters
#'
#' @param object A bayesICT object.
#' @param i Auxiliary group index, if present.
#' @return A display label.
bayesICT_pas_label <- function(object, i) {
  G <- bayesICT_aux_group_count(object)
  # Group-specific prevalence labels are only needed for auxiliary information.

  if (is.na(i) || G <= 1L) {
    "Proportion affirm. on sensitive"
  } else {
    paste0("Proportion affirm. on sensitive:", i)
  }
}

#' Recover the number of auxiliary groups
#'
#' @param object A bayesICT object.
#' @return An integer group count, defaulting to one when absent.
bayesICT_aux_group_count <- function(object) {
  G <- object$aux_info$G
  if (is.null(G) || length(G) == 0L || is.na(G)) {
    return(1L)
  }
  as.integer(G[1])
}

#' Test whether a bayesICT object uses a double-list model
#'
#' @param object A bayesICT object.
#' @return TRUE if the prior/model uses two lists, otherwise FALSE.
bayesICT_has_double_list <- function(object) {
  prior <- object$prior
  # The prior field records whether the object came from a double-list model.
  !is.null(prior) && prior %in% c("double_list", "double_list_aux")
}

#' Recover the number of warmup iterations used by Stan
#'
#' @param object A bayesICT object.
#' @return An integer warmup count, or NA if unavailable.
bayesICT_warmup_used <- function(object) {
  stanfit <- object$stanfit
  warmup2 <- tryCatch(stanfit@sim$warmup2[1], error = function(e) NULL)
  # Prefer Stan's stored warmup count when it is available.

  if (!is.null(warmup2) && length(warmup2) > 0L && !is.na(warmup2)) {
    return(warmup2)
  }

  if (!is.null(object$warmup)) {
    return(object$warmup)
  }

  NA_integer_
}

#' Compute posterior prevalence draws from a linear predictor
#'
#' @param X Design matrix for prediction.
#' @param sampledf Posterior draws as a data frame.
#' @param parameter Stan parameter base name, usually `delta` or `delta2`.
#' @return A matrix with observations in rows and posterior draws in columns.
bayesICT_linear_prevalence <- function(X, sampledf, parameter) {
  if (is.null(X)) {
    stop("The design matrix needed for prediction is missing.")
  }
  if (is.null(sampledf)) {
    stop("object$sampledf is missing.")
  }

  draws <- bayesICT_parameter_draws(sampledf, parameter, ncol(X))
  # Arrange draws as iterations by regression coefficient, then multiply by X.
  stats::plogis(X %*% t(draws))
  # Convert linear predictors to probabilities for each observation/draw pair.
}

#' Compute expected list counts under an ordered-logistic model
#'
#' Stan anchors the first cutpoint at zero and builds each remaining cutpoint
#' by adding exp(psiM[k]). For a count with `n_cutpoints` cutpoints, its
#' expectation is the sum of Pr(count >= k), k = 1, ..., n_cutpoints.
#'
#' @param X Design matrix for prediction.
#' @param sampledf Posterior draws as a data frame.
#' @param coefficient Coefficient parameter name.
#' @param cutpoints Cutpoint-spacing parameter name.
#' @param n_cutpoints Number of ordered-logistic cutpoints.
#' @return A matrix with observations in rows and posterior draws in columns.
bayesICT_ordered_count_draws <- function(X, sampledf, coefficient, cutpoints,
                                        n_cutpoints) {
  if (is.null(X)) {
    stop("The design matrix needed for prediction is missing.")
  }
  if (is.null(sampledf)) {
    stop("object$sampledf is missing.")
  }

  beta <- bayesICT_parameter_draws(sampledf, coefficient, ncol(X))
  eta <- X %*% t(beta)
  cuts <- matrix(0, nrow = ncol(eta), ncol = n_cutpoints)
  if (n_cutpoints > 1L) {
    spacings <- bayesICT_parameter_draws(sampledf, cutpoints, n_cutpoints - 1L)
    for (k in seq.int(2L, n_cutpoints)) {
      cuts[, k] <- cuts[, k - 1L] + exp(spacings[, k - 1L])
    }
  }

  expected_count <- matrix(0, nrow = nrow(X), ncol = ncol(eta))
  for (k in seq_len(n_cutpoints)) {
    expected_count <- expected_count +
      stats::plogis(sweep(eta, 2L, cuts[, k], FUN = "-"))
  }
  expected_count
}

#' Compute treated-minus-control expected count draws for difference-in-means
#'
#' @param X Design matrix for prediction.
#' @param sampledf Posterior draws as a data frame.
#' @param J Number of control items.
#' @param treated_coef,control_coef Stan coefficient parameter names.
#' @param treated_cuts,control_cuts Stan cutpoint-spacing parameter names.
#' @return A matrix with observations in rows and posterior draws in columns.
bayesICT_diff_means_prevalence <- function(X, sampledf, J, treated_coef,
                                           control_coef, treated_cuts,
                                           control_cuts) {
  bayesICT_ordered_count_draws(
    X, sampledf, treated_coef, treated_cuts, J + 1L
  ) - bayesICT_ordered_count_draws(
    X, sampledf, control_coef, control_cuts, J
  )
}


#' Test whether posterior draws for a parameter are available
#'
#' @param sampledf Posterior draws as a data frame.
#' @param parameter Stan parameter base name.
#' @return TRUE if matching columns exist.
bayesICT_has_parameter_draws <- function(sampledf, parameter) {
  if (is.null(sampledf)) {
    return(FALSE)
  }

  nms <- names(sampledf)
  if (is.null(nms)) {
    return(FALSE)
  }

  pattern <- paste0("^", parameter, "(\\[|\\.|$)")
  # Match both Stan-style names delta[1] and data-frame names such as delta.1.
  any(grepl(pattern, nms))
}

#' Extract posterior draws for one vector-valued parameter
#'
#' @param sampledf Posterior draws as a data frame.
#' @param parameter Stan parameter base name.
#' @param ncol_expected Optional expected number of columns.
#' @return A numeric matrix of posterior draws.
bayesICT_parameter_draws <- function(sampledf, parameter, ncol_expected = NULL) {
  nms <- names(sampledf)
  pattern <- paste0("^", parameter, "(\\[|\\.|$)")
  # Match both Stan-style names delta[1] and data-frame names such as delta.1.
  loc <- grep(pattern, nms)
  # Locate all columns belonging to this possibly vector-valued parameter.

  if (length(loc) == 0L) {
    stop("Could not find posterior draws for parameter '", parameter, "'.")
  }

  draws <- as.matrix(sampledf[, loc, drop = FALSE])
  # Keep a matrix even when the parameter is scalar.

  if (!is.null(ncol_expected) && ncol(draws) != ncol_expected) {
    stop(
      "Parameter '", parameter, "' has ", ncol(draws),
      " columns in object$sampledf, but the corresponding design matrix has ",
      ncol_expected, " columns."
    )
  }

  draws
}

#' Summarize each row of a posterior-draw matrix
#'
#' @param draws Matrix with observations or quantities in rows and posterior draws in columns.
#' @param type Summary type, either "mean" or "median".
#' @return A numeric vector of row summaries.
bayesICT_summarize_rows <- function(draws, type) {
  if (identical(type, "mean")) {
    # Prediction summaries are computed independently for each observation.
    apply(draws, 1L, mean, na.rm = TRUE)
  } else {
    apply(draws, 1L, stats::quantile, probs = 0.5, na.rm = TRUE)
  }
}
