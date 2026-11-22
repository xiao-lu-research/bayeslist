# Map predictions for unique covariate patterns back to respondents, then
# average across respondents within each posterior draw.
bayesICT_weighted_pas_draws <- function(p, idx, weights) {
  if (is.null(weights)) weights <- rep.int(1, length(idx))
  if (!is.numeric(weights) || length(weights) != length(idx) ||
      any(!is.finite(weights)) || any(weights < 0) || sum(weights) <= 0) {
    stop("weights must be finite, nonnegative, and match the number of respondents")
  }
  as.vector(crossprod(weights / sum(weights), p[idx, , drop = FALSE]))
}

bayesICT_pas_interval_table <- function(p, idx, X, weights,
                                       prob, center, group_by) {
  if (!is.numeric(prob) || length(prob) != 1L || !is.finite(prob) ||
      prob <= 0 || prob >= 1) stop("prob must be between 0 and 1")
  if (is.null(group_by)) group_by <- setdiff(colnames(X), "(Intercept)")
  if (!is.character(group_by) || anyNA(group_by) ||
      anyDuplicated(group_by) || !all(group_by %in% colnames(X)) ||
      "(Intercept)" %in% group_by) {
    stop("group_by must name non-intercept columns of object$X")
  }

  if (length(group_by)) {
    group_X <- X[, group_by, drop = FALSE]
    key <- do.call(paste, c(as.data.frame(group_X), sep = "\r"))
    keep <- !duplicated(key)
    group_id <- match(key, key[keep])
    group_X <- group_X[keep, , drop = FALSE]
  } else {
    group_id <- rep.int(1L, nrow(X))
    group_X <- X[1, FALSE, drop = FALSE]
  }
  n_groups <- nrow(group_X)
  group_draws <- t(matrix(vapply(seq_len(n_groups), function(g) {
    rows <- which(group_id == g)
    group_weights <- if (is.null(weights)) NULL else weights[rows]
    bayesICT_weighted_pas_draws(p, idx[rows], group_weights)
  }, numeric(ncol(p))), nrow = ncol(p), ncol = n_groups))

  alpha <- (1 - prob) / 2
  probs <- c(alpha, 0.5, 1 - alpha)
  q <- t(vapply(seq_len(n_groups), function(g) {
    stats::quantile(group_draws[g, ], probs = probs, names = FALSE)
  }, numeric(3)))
  out <- data.frame(group_X,
                    fit = if (center == "median") q[, 2] else rowMeans(group_draws),
                    lower = q[, 1], upper = q[, 3])
  interval_cols <- (ncol(out) - 1):ncol(out)
  names(out)[interval_cols] <- paste0(c("lower", "upper"), format(prob * 100))
  rownames(out) <- NULL
  out
}

posterior_interval_ord.bayesICT <- function(object, prob = 0.95,
                                  center = c("median", "mean"),
                                  which_list = c("average", "A", "B"),
                                  check = TRUE, group_by = NULL) {
  center <- match.arg(center)
  which_list <- match.arg(which_list)
  sd_ <- object$sampledf
  nms <- names(sd_)
  X   <- as.matrix(object$X)
  
  getm <- function(pat) {
    cols <- grep(pat, nms)
    if (!length(cols)) stop("No columns in sampledf match ", pat)
    as.matrix(sd_[, cols, drop = FALSE])
  }
  
  # list A (always present)
  psi0 <- getm("^psi0\\[");  delta <- getm("^delta\\[")
  psiM <- getm("^psiM\\[");  psiMT <- getm("^psiMT\\[")
  stopifnot(ncol(delta) == ncol(X), ncol(psi0) == ncol(X))
  
  # list B (only in double-list models)
  has_B <- all(vapply(c("^psi02\\[", "^delta2\\[", "^psiM2\\[", "^psiMT2\\["),
                      function(p) any(grepl(p, nms)), logical(1)))
  if (has_B) {
    X2 <- if (!is.null(object$X2)) as.matrix(object$X2) else X
    psi02 <- getm("^psi02\\["); delta2 <- getm("^delta2\\[")
    psiM2 <- getm("^psiM2\\["); psiMT2 <- getm("^psiMT2\\[")
    stopifnot(ncol(delta2) == ncol(X2), ncol(psi02) == ncol(X2))
  } else {
    X2 <- X
    if (which_list == "B") stop("This model has only one list.")
    message("Single-list model detected: using list A only.")
  }
  S <- nrow(delta)
  
  # one row per distinct covariate pattern
  key  <- apply(cbind(X, X2), 1, paste, collapse = "_")
  keep <- !duplicated(key)
  Xu   <- X[keep, , drop = FALSE]
  X2u  <- X2[keep, , drop = FALSE]
  
  exp_count <- function(eta, cuts) length(cuts) - rowSums(plogis(outer(-eta, cuts, "+")))
  make_cuts <- function(psi) c(0, cumsum(exp(psi)))
  
  pA <- matrix(NA_real_, nrow(Xu), S)
  pB <- if (has_B) matrix(NA_real_, nrow(Xu), S) else NULL
  for (s in seq_len(S)) {
    pA[, s] <- exp_count(Xu %*% delta[s, ], make_cuts(psiMT[s, ])) -
      exp_count(Xu %*% psi0[s, ],  make_cuts(psiM[s, ]))
    if (has_B) {
      pB[, s] <- exp_count(X2u %*% delta2[s, ], make_cuts(psiMT2[s, ])) -
        exp_count(X2u %*% psi02[s, ],  make_cuts(psiM2[s, ]))
    }
  }
  
  p_both <- if (has_B) 0.5 * (pA + pB) else pA
  p <- if (!has_B) pA else switch(which_list, average = p_both, A = pA, B = pB)
  idx <- match(key, key[keep])

  if (check && "pas" %in% nms) {
    pas_calc <- bayesICT_weighted_pas_draws(p_both, idx, object$weights)
    message(sprintf("Check: recomputed overall proportion = %.4f; model's pas = %.4f",
                    mean(pas_calc), mean(sd_$pas)))
  }
  bayesICT_pas_interval_table(p, idx, X, object$weights, prob, center, group_by)
}

#' Compute posterior uncertainty intervals from a \code{bayesICT} object
#'
#' Return survey-weighted PAS intervals for subgroups. Respondent predictions
#' are averaged with survey weights within each subgroup and posterior draw
#' before intervals are computed.
#'
#' @param object A \code{bayesICT} object from running the main function \code{\link{bayesICT}}.
#' @param prob Posterior interval probability, between 0 and 1.
#' @param center character indicating if the measure of central 
#' location should be the median (default) or the mean.
#' @param which_list character indicating if the probability of 
#' affirmative on the sensitive item should be the average of the 
#' two lists, list A alone or list B alone. 
#' @param check Logical; compare the reconstructed overall PAS with Stan's
#'   scalar \code{pas} when one is available.
#' @param group_by Character vector naming columns of \code{object$X} that define
#'   subgroups. By default, all non-intercept columns define groups.
#' @param ... Reserved for future use.
#'
#' @return A data frame with one survey-weighted PAS interval per subgroup.
#' @export
#' @importFrom stats median plogis
#'
bayesICT_posterior_interval <- function(object, prob = 0.95,
                                  center = c("median", "mean"),
                                  which_list = c("average", "A", "B"),
                                  check = TRUE, group_by = NULL, ...) {
  center <- match.arg(center)
  which_list <- match.arg(which_list)
  sd_  <- object$sampledf
  nms  <- names(sd_)
  
  # ---- detect the model type ------------------------------------------------
  model_code <- tryCatch(object$stanfit@stanmodel@model_code, error = function(e) NULL)
  is_ordinal <- any(grepl("^psiMT\\[", nms)) ||
    (!is.null(model_code) && any(grepl("ordered_logistic_expectation", model_code)))
  
  if (is_ordinal) {
    message("Ordered-logistic model detected: using predict_ord.bayeslist()")
    return(posterior_interval_ord.bayesICT(object, prob = prob, center = center,
                                 which_list = which_list, check = check,
                                 group_by = group_by))
  }
  message("No ordered-logistic structure found: assuming prevalence = logistic(X %*% delta)")
  
  # ---- logistic branch --------------------------------------------------------
  X  <- as.matrix(object$X)
  X2 <- if (!is.null(object$X2)) as.matrix(object$X2) else X
  
  colsA <- grep("^delta\\[", nms)
  colsB <- grep("^delta2\\[", nms)
  if (length(colsA) != ncol(X)) {
    stop("object$X has ", ncol(X), " columns but ", length(colsA),
         " columns match '^delta\\['. This model is not a plain logistic model;",
         " check the Stan code before using these predictions.")
  }
  has_B <- length(colsB) > 0
  if (has_B && length(colsB) != ncol(X2)) {
    stop("delta2 has ", length(colsB), " columns but X2 has ", ncol(X2))
  }
  if (which_list == "B" && !has_B) stop("This model has only one list (no delta2 columns).")
  
  # one row per distinct covariate pattern
  key  <- apply(cbind(X, X2), 1, paste, collapse = "_")
  keep <- !duplicated(key)
  Xu   <- X[keep, , drop = FALSE]
  X2u  <- X2[keep, , drop = FALSE]
  
  pA <- plogis(Xu %*% t(as.matrix(sd_[, colsA, drop = FALSE])))      # patterns x draws
  pB <- if (has_B) plogis(X2u %*% t(as.matrix(sd_[, colsB, drop = FALSE]))) else NULL
  
  p <- if (!has_B || which_list == "A") pA
  else if (which_list == "B") pB
  else (pA + pB) / 2
  idx <- match(key, key[keep])

  # ---- compare with the model's own overall estimate ------------------------
  if (check && "pas" %in% nms) {
    p_all <- if (has_B) (pA + pB) / 2 else pA
    pas_calc <- bayesICT_weighted_pas_draws(p_all, idx, object$weights)
    message(sprintf("Check: recomputed overall proportion = %.4f; model's pas = %.4f",
                    mean(pas_calc), mean(sd_$pas)))
    if (abs(mean(pas_calc) - mean(sd_$pas)) > 0.01) {
      warning("Recomputed proportion differs from the model's pas by more than 0.01;",
              " the logistic mapping may not match this model.")
    }
  }
  
  bayesICT_pas_interval_table(p, idx, X, object$weights, prob, center, group_by)
}
