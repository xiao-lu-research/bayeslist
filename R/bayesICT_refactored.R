#' Fitting Bayesian sensitive item models using the item count technique
#'
#' The main function for estimating Bayesian item count technique (ICT) models
#' to list questions with sensitive item.
#'
#' This represents the classical model for list experiments, where the answer to the sensitive item is regressed
#' on and potentially explained by covariates. 
#'
#' The function returns a \code{bayeslist} object that can be further investigated using standard functions such as \code{summary}, \code{plot}, \code{print}, \code{predict}, and \code{coef}. The model can be passed using a \code{formula} as in \code{lm()}. Convergence diagnostics can be performed using either \code{print(object, "mcmc")} or \code{plot(object, "trace")}.
#'
#' @param formula An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' The left-hand-side variable is the list variable. If it is a double or direct type of model, the left-hand-side variable
#' is the first list variable. 
#' @param data A data frame containing the variables in the model.
#' @param treat Variable name of the treatment.
#' @param J Number of control items.
#' @param type Type of the model. Options include "outcome" for the item count technique (ICT) models.
#' This represents the classical model for list experiments, where the 
#' answer to the sensitive item is regressed
#' on, and potentially explained by, covariates. 
#' @param nsim The number of iterations.
#' @param warmup The number of warmup iterations before samples are taken.
#' @param thin Thinning parameter.
#' @param CIsize The size of posterior confidence interval.
#' @param nchain The number of parallel chains.
#' @param control List. control parameters to pass to
#' the underlying `rstan` function. Defaults to empty.
#' @param seeds Random seeds to replicate the results.
#' @param prior Prior types. Options include NULL, "auxiliary", "double_list", and "double_list_aux".
#' If NULL, no informative priors will be used. Auxiliary priors are not
#' available with `list_variable_model = "diff_means"`; the single-list
#' auxiliary binomial model and beta-binomial model without bias are also
#' unavailable.
#' @param double_list Variable name of the second list.
#' @param double_list_treat Variable name of the treatment variable of the second list.
#' @param aux_info Auxiliary information for the informative priors. list(G,h,sigma_h,g), where: G (number of subgroups), h (auxiliary information for each subgroup), sigma_h (standard deviation for each subgroup), and g (subgroup indicator). If it is NULL, the following two parameters need to be specified when estimating the model with prior = "auxiliary".
#' @param aux_g Auxiliary information for the informative priors: name of the variable indicating the group of each observation.
#' @param aux_h Auxiliary information for the informative priors: name of the variable containing information of prevalence for each group
#' @param bias character If there is auxiliary information, should a parameter
#' for the bias of the sensitive item response be added to the model.
#' If "common" then a single bias parameter is incorporated for
#' double and single lists (the default). If
#' "separate" then a two bias parameters are incorporated for double lists
#' (one per list). If "no" then no bias terms are incorporated.
#' @param sigma_h Auxiliary information for the informative priors: name of the variable containing information of the standard deviation of prevalence for each group
#' @param mu_delta vector of means for the prior on delta
#' @param mu_psi0 vector of means for the prior on psi0
#' @param mu_psiM vector of means for the prior on psiM
#' @param mu_rho0 mean for the prior on variance parameter without treatment
#' @param sigma_delta vector of standard deviations for the prior on delta
#' @param sigma_psi0 vector of standard deviations for the prior on psi0
#' @param sigma_psiM vector of standard deviations for the prior on psiM
#' @param sigma_bias standard deviation for the prior on bias/bias2 of the sensitive item bias
#' @param sigma_rho0 standard deviation for the prior on variance parameter without treatment
#' @param list_variable_model Model for the control-list count: "beta-binomial",
#' "binomial", "multinomial", or "diff_means". The default is "beta-binomial".
#' The binomial option uses a binomial control-count likelihood, including with
#' two-list auxiliary priors; a single-list auxiliary binomial model is not
#' implemented. The difference-in-means option uses an ordered-logistic count
#' model and is available only without auxiliary priors.
#' @param weights Variable name of the survey weights. This should be a variable in the data data frame. By default it is NULL, meaning equal weights.
#' @param init_fun Optional. Set initial values by providing a function that returns a list for specifying the initial values of parameters  
#' @param parallel Logical or count. Indicating whether to do parallel computing. If it is TRUE then all 
#' processors are used. If it is a count then that number is used. The default is \code{nchain}.
#'
#' @return A \code{bayeslist} object. An object of class \code{bayeslist} contains the following elements
#'
#'   \describe{
#'
#'   \item{\code{Call}}{The matched call.}
#'   \item{\code{formula}}{Symbolic representation of the model.}
#'   \item{\code{type}}{Model type}
#'   \item{\code{nsim}}{Number of iterations.}
#'   \item{\code{Warmup}}{Number of warmup iterations.}
#'   \item{\code{thin}}{Thinning.}
#'   \item{\code{seeds}}{Random seeds for reproducibility. The default is 12345.}
#'   \item{\code{CIsize}}{Size of the posterior confidence interval.}
#'   \item{\code{data}}{Data used.}
#'   \item{\code{X}}{Independent variables.}
#'   \item{\code{Y}}{Dependent variables.}
#'   \item{\code{xnames}}{Names of the independent variables.}
#'   \item{\code{stanfit}}{Output from stan.}
#'   \item{\code{sampledf}}{Posterior samples.}
#'   \item{\code{summaryout}}{Summary of the stan-fit object.}
#'   \item{\code{npars}}{Number of control variables.}
#'   \item{\code{prior}}{Informative prior types.}
#'   \item{\code{double_list}}{The second list.}
#'   \item{\code{aux_info}}{Auxiliary information.}
#'   \item{\code{ulbs}}{Upper and lower bounds based on the specified confidence interval.}
#'   \item{\code{means}}{Mean estimates.}
#'   \item{\code{treat}}{Treatment.}
#'   \item{\code{outcome}}{Outcome to be predicted.}
#'
#' }
#'
#' @references Lu, X. and Traunmüller, R. (2021). Improving Studies of Sensitive Topics Using Prior Evidence: A Unified Bayesian Framework for List Experiments.
#'
#' @import Formula
#' @import rstan
#'
#' @importFrom stats coef model.frame model.matrix quantile as.formula binomial density glm sd runif rnorm update
#' @importFrom rlang new_function
#'
#' @export
#'
#' @examples
#'\dontrun{
#'# Estimate sensitive item outcome model using Sri Lanka data on male sexual violence
#'# Load Sri Lanka list experiment data
#'data(srilanka)
#'
#'# Model 1: intercept-only outcome model without prior information:
#'mod1 <- bayesICT(sexassault ~ 1, data = srilanka, treat = "treatment", J = 3, type = "outcome",
#'      nsim = 2000, thin = 1, CIsize = 0.95, nchain = 2, seeds = 342321, prior = NULL, parallel = 2)
#'summary(mod1) # summary of estimates
#'predict(mod1) # predicted prevalence for each observation
#'plot(mod1,"trace") # trace plot
#'plot(mod1,"coef") # coefficient plot
#'plot(mod1, only_prev = TRUE) # prevalence plot
#'
#'# Model 2: multivariate outcome model without prior information:
#'mod2 <- bayesICT(sexassault ~ age + edu, data = srilanka, treat = "treatment",
#'     J = 3, type = "outcome", 
#'     nsim = 2000, thin = 1, CIsize = 0.95, nchain = 2, seeds = 342321, prior = NULL, parallel = 2)
#'summary(mod2) # summary of estimates
#'predict(mod2) # predicted prevalence for each observation
#'plot(mod2,"trace") # trace plot
#'plot(mod2,"coef") # coefficient plot
#'plot(mod2) # prevalence + coefficient plot
#'}
#'
bayesICT <- function(formula,
                     data,
                     treat,
                     J = NULL,
                     type = c("outcome", "predict", "misreport"),
                     nsim = 1000,
                     warmup = NULL,
                     thin = 1,
                     CIsize = .95,
                     nchain = 1,
                     control = NULL,
                     seeds = 12345,
                     prior = NULL, # options: 1. "auxiliary", 2. "double_list" 3. "double_list_aux" 4. NULL
                     double_list = NULL, # variable name for the outcome of the second list
                     double_list_treat = NULL, # variable name of the treatment variable for the second list
                     aux_info = NULL, # list(G,h,g), where: G (number of subgroups), h (auxiliary information for each subgroup) sigma_h (the standard deviation of the auxiliary information for each subgroup) and g (subgroup indicator).
                     aux_g = NULL, # name of the variable indicating the group of each observation
                     aux_h = NULL, # name of the variable containing information of prevalence for each group
                     bias = NULL, # If there is auxiliary information, should a parameter for the bias of the sensitive item response be added to the model 
                     sigma_h= NULL, # name of the variable containing information of the uncertainty in prevalence for each group
                     mu_delta = 0, # vector of means for the prior on delta
                     mu_psi0 = 0, # vector of means for the prior on psi0
                     mu_psiM = 0.0, # vector of means for the prior on psi0
                     mu_rho0 = 0.2, # mean for the prior on rho0
                     sigma_delta = 10, # vector of standard deviations for the prior on delta
                     sigma_psi0 = 10, # vector of standard deviations for the prior on psi0
                     sigma_psiM = 10, # vector of standard deviations for the prior on psiM
                     sigma_bias = 10, # standard deviation for the prior on bias/bias2 of the sensitive item bias
                     sigma_rho0= 0.3, # standard deviation for the prior on variance parameter without treatment
                     list_variable_model = "beta-binomial", # model for the outcome variable: either "binary" or "multinomial". The default is binary.
                     weights = NULL, # survey weights, if any
                     init_fun = NULL, # Set initial values by providing a function that returns a list for specifying the initial values of parameters  
                     # direct = NULL, # data of direct item for the misreport model
                     parallel = nchain) {
  # Dispatch on the high-level model type before doing any outcome-specific
  # validation or data preparation.  This keeps the classical outcome model
  # separate from future predict/misreport implementations.
  type <- match.arg(type)

  switch(
    type,
    outcome = bayesICT_outcome(
      formula = formula,
      data = data,
      treat = treat,
      J = J,
      type = type,
      nsim = nsim,
      warmup = warmup,
      thin = thin,
      CIsize = CIsize,
      nchain = nchain,
      control = control,
      seeds = seeds,
      prior = prior,
      double_list = double_list,
      double_list_treat = double_list_treat,
      aux_info = aux_info,
      aux_g = aux_g,
      aux_h = aux_h,
      bias = bias,
      sigma_h = sigma_h,
      mu_delta = mu_delta,
      mu_psi0 = mu_psi0,
      mu_psiM = mu_psiM,
      mu_rho0 = mu_rho0,
      sigma_delta = sigma_delta,
      sigma_psi0 = sigma_psi0,
      sigma_psiM = sigma_psiM,
      sigma_bias = sigma_bias,
      sigma_rho0 = sigma_rho0,
      list_variable_model = list_variable_model,
      weights = weights,
      init_fun = init_fun,
      parallel = parallel,
      call = match.call()
    ),
    predict = stop("type = 'predict' is not yet implemented in bayesICT()."),
    misreport = stop("type = 'misreport' is not yet implemented in bayesICT().")
  )
}

#' Fit the classical ICT outcome model
#'
#' This internal worker performs validation, data construction, Stan model
#' selection, sampling, and assembly of the returned bayesICT object for
#' type = "outcome".
#'
#' @keywords internal
bayesICT_outcome <- function(formula,
                             data,
                             treat,
                             J = NULL,
                             type = "outcome",
                             nsim = 1000,
                             warmup = NULL,
                             thin = 1,
                             CIsize = .95,
                             nchain = 1,
                             control = NULL,
                             seeds = 12345,
                             prior = NULL,
                             double_list = NULL,
                             double_list_treat = NULL,
                             aux_info = NULL,
                             aux_g = NULL,
                             aux_h = NULL,
                             bias = NULL,
                             sigma_h = NULL,
                             mu_delta = 0,
                             mu_psi0 = 0,
                             mu_psiM = 0.0,
                             mu_rho0 = 0.2,
                             sigma_delta = 10,
                             sigma_psi0 = 10,
                             sigma_psiM = 10,
                             sigma_bias = 10,
                             sigma_rho0 = 0.3,
                             list_variable_model = "beta-binomial",
                             weights = NULL,
                             init_fun = NULL,
                             parallel = nchain,
                             call = match.call()) {
  if (missing(formula) || missing(data)) {
    stop("Formula and data should be given.")
  }

  # Normalize and validate the MCMC dimensions before any model-dependent work.
  timings <- bayesICT_validate_mcmc_args(
    J = J,
    nsim = nsim,
    warmup = warmup,
    thin = thin,
    CIsize = CIsize
  )
  J <- timings$J
  warmup <- timings$warmup

  # Convert NULL prior into the internal value "none" and validate the
  # list-variable likelihood family.
  prior_kind <- bayesICT_match_prior(prior)
  list_variable_model <- match.arg(
    list_variable_model,
    choices = c("beta-binomial", "binomial", "multinomial", "diff_means")
  )
  # Auxiliary priors are the only cases where bias parameters can enter.
  # If the user does not specify a bias convention, use the historical default.
  if (is.null(bias)) {
    bias <- if (prior_kind %in% c("auxiliary", "double_list_aux")) "common" else "no"
  }
  bias <- match.arg(bias, choices = c("no", "common", "separate"))

  ## The original bayesICT() ignored `bias` for non-auxiliary models.
  ## Preserve that behavior explicitly: no bias parameter is present unless
  ## prior is "auxiliary" or "double_list_aux".
  if (!(prior_kind %in% c("auxiliary", "double_list_aux"))) {
    bias <- "no"
  }

  if (list_variable_model == "diff_means" &&
      prior_kind %in% c("auxiliary", "double_list_aux")) {
    stop("Auxiliary difference-in-means models are not implemented; choose a different list_variable_model or omit the auxiliary prior.")
  }
  if (list_variable_model == "binomial" && prior_kind == "auxiliary") {
    stop("A single-list auxiliary binomial model is not implemented; choose beta-binomial or multinomial, or omit the auxiliary prior.")
  }
  if (list_variable_model == "beta-binomial" && prior_kind == "auxiliary" &&
      bias == "no") {
    stop("A single-list auxiliary beta-binomial model without bias is not implemented; use bias = 'common' or choose another model.")
  }

  # Bundle the discrete modeling choices into a compact specification object.
  # Downstream helpers use this object rather than repeating nested conditionals.
  spec <- bayesICT_make_spec(
    prior_kind = prior_kind,
    list_variable_model = list_variable_model,
    bias = bias
  )

  # Check that auxiliary and double-list inputs are present when required by
  # the specification.
  bayesICT_validate_model_args(
    spec = spec,
    aux_info = aux_info,
    aux_g = aux_g,
    aux_h = aux_h,
    sigma_h = sigma_h,
    double_list = double_list,
    double_list_treat = double_list_treat
  )

  # rstan consults the global mc.cores option.  Preserve the user's setting
  # and restore it when the function exits.
  old_mc_cores <- getOption("mc.cores")
  on.exit(options(mc.cores = old_mc_cores), add = TRUE)
  cores <- bayesICT_resolve_cores(parallel)
  options(mc.cores = cores)

  # Build a single complete model frame containing every variable needed by
  # this specification, then perform listwise deletion once.
  data_prepared <- bayesICT_prepare_model_frame(
    formula = formula,
    data = data,
    treat = treat,
    weights = weights,
    spec = spec,
    aux_info = aux_info,
    aux_g = aux_g,
    aux_h = aux_h,
    sigma_h = sigma_h,
    double_list = double_list,
    double_list_treat = double_list_treat
  )
  data0 <- data_prepared$data
  weights_vec <- data_prepared$weights

  if (nrow(data0) < 5) {
    warning("Less than 5 observations after list-wise deletion. Please check your data and missing values.")
  }

  # Extract the primary list response, treatment indicator, design matrix,
  # and normalized survey weights.
  data_info <- bayesICT_make_primary_data(
    formula = formula,
    data = data0,
    treat = treat,
    weights = weights_vec
  )

  # Expand scalar prior hyperparameters to vectors of the correct length and
  # reject ambiguous partially specified vectors.
  prior_info <- bayesICT_prepare_prior_vectors(
    K = data_info$K,
    J = J,
    mu_delta = mu_delta,
    mu_psi0 = mu_psi0,
    mu_psiM = mu_psiM,
    mu_rho0 = mu_rho0,
    sigma_delta = sigma_delta,
    sigma_psi0 = sigma_psi0,
    sigma_psiM = sigma_psiM,
    sigma_bias = sigma_bias,
    sigma_rho0 = sigma_rho0
  )

  # Auxiliary information is optional overall, but mandatory for auxiliary
  # prior specifications. Store it separately from the main model frame.
  aux_data <- NULL
  if (spec$has_aux) {
    aux_data <- bayesICT_make_aux_data(
      data = data0,
      aux_info = aux_info,
      aux_g = aux_g,
      aux_h = aux_h,
      sigma_h = sigma_h,
      N = data_info$N,
      weights = data_info$weights
    )
  }

  # Double-list models add a second response, design matrix, and treatment
  # indicator while retaining the primary-list design information.
  double_data <- NULL
  if (spec$has_double_list) {
    double_data <- bayesICT_make_double_list_data(
      formula = formula,
      data = data0,
      double_list = double_list,
      double_list_treat = double_list_treat,
      default_treat = treat
    )
  }

  # Select the compiled Stan model, build the Stan data list, choose monitored
  # parameters, and construct default initial values if needed.
  stan_spec <- bayesICT_make_stan_spec(
    spec = spec,
    J = J,
    data_info = data_info,
    prior_info = prior_info,
    aux_data = aux_data,
    double_data = double_data,
    init_fun = init_fun
  )

  # Assemble rstan::sampling() arguments once so every model branch uses the
  # same sampling path.
  sampling_args <- list(
    object = stan_spec$stanmodel,
    data = stan_spec$datlist,
    pars = stan_spec$pars,
    seed = seeds,
    iter = nsim,
    control = control,
    thin = thin,
    warmup = warmup,
    chains = nchain,
    cores = cores
  )
  if (!is.null(stan_spec$init)) {
    sampling_args$init <- stan_spec$init
  }

  # Fit the selected Stan model.
  stanout <- do.call(rstan::sampling, sampling_args)

  # Keep both the rstan summary matrix and a flat data-frame version of the
  # posterior draws for downstream print, plot, and predict methods.
  summaryout <- rstan::summary(stanout)$summary
  sampledf <- as.data.frame(stanout)

  # Warn if monitored parameters appear fixed at their initial values.  This
  # is a coarse diagnostic intended to catch obvious sampling failures.
  no_sampling <- apply(
    sampledf[, 1:min(length(stan_spec$pars), ncol(sampledf)), drop = FALSE],
    2,
    stats::sd,
    na.rm = TRUE
  ) < 1e-16
  if (any(no_sampling)) {
    warning("Some of the parameters did not vary from their initial values in the MCMC sampling.")
    print(no_sampling)
  }

  # Assemble the public S3 object.  The object stores both user-facing fields
  # and internal Stan metadata needed by summary(), plot(), predict(), and coef().
  out <- list()
  class(out) <- c("bayesICT", class(out))
  out$Call <- call
  out$formula <- formula
  out$J <- J
  out$type <- type
  out$nsim <- nsim
  out$warmup <- warmup
  out$thin <- thin
  out$seeds <- seeds
  out$CIsize <- CIsize
  out$data <- data0
  out$X <- data_info$X
  out$Y <- data_info$Y
  out$xnames <- colnames(data_info$X)
  out$stanfit <- stanout
  out$X2 <- if (!is.null(double_data)) double_data$X2 else NULL
  out$weights <- weights_vec
  out$sampledf <- sampledf
  out$summaryout <- summaryout
  out$npars <- data_info$K
  out$prior <- if (identical(prior_kind, "none")) NULL else prior_kind
  out$double_list <- if (!is.null(double_data)) double_data$Y2 else NULL
  out$double_list_var <- double_list
  out$double_list_treat <- double_list_treat
  out$aux_info <- aux_data
  out$ulbs <- apply(
    sampledf,
    2,
    stats::quantile,
    probs = c((1 - CIsize) / 2, 1 - (1 - CIsize) / 2)
  )
  out$means <- apply(sampledf, 2, mean)
  out$treat <- treat
  out$treat_list <- data_info$treat_list
  out$list_variable_model <- list_variable_model
  out$bias <- bias
  out$stanmodel <- stan_spec$stanmodel
  out$stan_data <- stan_spec$datlist
  out$stan_pars <- stan_spec$pars
  out$no_sampling <- no_sampling

  out
}

#' Normalize the prior argument
#'
#' Convert the public prior encoding, where NULL means no informative prior,
#' into the internal string representation used by the helper functions.
#'
#' @keywords internal
bayesICT_match_prior <- function(prior) {
  # The public API uses NULL for the no-prior model; internally strings are
  # simpler to match and combine into model-selection keys.
  if (is.null(prior)) {
    return("none")
  }
  match.arg(prior, choices = c("auxiliary", "double_list", "double_list_aux"))
}

#' Construct the internal model specification
#'
#' Collect the prior family, list-variable likelihood family, and auxiliary-bias
#' convention into a compact object used throughout the fitting code.
#'
#' @keywords internal
bayesICT_make_spec <- function(prior_kind, list_variable_model, bias) {
  # The logical flags are redundant with prior_kind, but make downstream tests
  # more readable and reduce repeated membership checks.
  list(
    prior_kind = prior_kind,
    list_variable_model = list_variable_model,
    bias = bias,
    has_aux = prior_kind %in% c("auxiliary", "double_list_aux"),
    has_double_list = prior_kind %in% c("double_list", "double_list_aux")
  )
}

#' Validate MCMC and interval arguments
#'
#' Check scalar MCMC controls and return normalized values used by the fitting
#' routine.
#'
#' @keywords internal
bayesICT_validate_mcmc_args <- function(J, nsim, warmup, thin, CIsize) {
  # J defines the support of the list response and the length of cutpoint
  # priors, so require a single positive integer.
  if (is.null(J) || !is.numeric(J) || length(J) != 1L || !is.finite(J) || J <= 0) {
    stop("The number of control items, J, must be specified and be a positive integer.")
  }
  if (J != as.integer(J)) {
    stop("The number of control items, J, must be an integer.")
  }
  J <- as.integer(J)

  # Preserve the original convention: use half the iterations as warmup when
  # the user does not provide a value.
  if (is.null(warmup)) {
    warmup <- floor(nsim / 2)
  }
  if (warmup < 0) {
    stop("Warm-up must be non-negative.")
  }
  if (thin < 1) {
    stop("Thinning factor must be positive.")
  }
  if (CIsize <= 0) {
    stop("Confidence interval size 'CIsize' must be positive.")
  }
  if (CIsize > 1) {
    stop("Confidence interval size 'CIsize' can not be larger than 1.")
  }

  list(J = J, warmup = warmup)
}

#' Validate model-dependent inputs
#'
#' Check that auxiliary and double-list variables are supplied when the selected
#' model specification requires them.
#'
#' @keywords internal
bayesICT_validate_model_args <- function(spec,
                                         aux_info,
                                         aux_g,
                                         aux_h,
                                         sigma_h,
                                         double_list,
                                         double_list_treat) {
  # If aux_info is absent, the auxiliary information must be reconstructible
  # from variables in the data frame.
  if (spec$has_aux && is.null(aux_info)) {
    missing_aux <- c(
      if (is.null(aux_g)) "aux_g",
      if (is.null(aux_h)) "aux_h",
      if (is.null(sigma_h)) "sigma_h"
    )
    if (length(missing_aux) > 0L) {
      stop(
        "For auxiliary priors, provide aux_info or all of aux_g, aux_h, and sigma_h. Missing: ",
        paste(missing_aux, collapse = ", "),
        "."
      )
    }
  }

  # Double-list specifications cannot be fit without the second list response.
  if (spec$has_double_list && is.null(double_list)) {
    stop("No double list item indicated.")
  }

  if (spec$has_double_list && is.null(double_list_treat)) {
    warning("double_list_treat is NULL; using treat as the treatment variable for the second list.")
  }

  invisible(TRUE)
}

#' Resolve the requested number of parallel cores
#'
#' Translate TRUE, FALSE, or an integer-like value into the number of cores to
#' pass to rstan.
#'
#' @keywords internal
bayesICT_resolve_cores <- function(parallel) {
  # TRUE means use all detected cores; FALSE means force serial sampling.
  if (is.logical(parallel)) {
    if (isTRUE(parallel)) {
      return(parallel::detectCores())
    }
    return(1L)
  }

  cores <- as.integer(round(parallel))
  if (!is.finite(cores) || cores < 1L) {
    stop("parallel must be TRUE, FALSE, or a positive integer count.")
  }
  cores
}

#' Prepare the complete model frame
#'
#' Gather all variables required by the selected model, normalize weights, and
#' perform listwise deletion once before constructing model matrices.
#'
#' @keywords internal
bayesICT_prepare_model_frame <- function(formula,
                                         data,
                                         treat,
                                         weights,
                                         spec,
                                         aux_info,
                                         aux_g,
                                         aux_h,
                                         sigma_h,
                                         double_list,
                                         double_list_treat) {
  # Work on a data-frame copy so we can add normalized weights without
  # modifying the user's original object.
  data <- as.data.frame(data)

  if (!is.null(weights)) {
    if (is.null(data[[weights]])) {
      stop("The weights variable '", weights, "' does not exist.")
    }
    weights_vec <- data[[weights]]
    # Treat non-finite weights as zero before normalization.
    weights_vec[!is.finite(weights_vec)] <- 0
    if (any(weights_vec < 0)) {
      stop("The weights variable must be non-negative.")
    }
    if (sum(weights_vec) <= 0) {
      stop("The weights variable must have positive finite total weight.")
    }
    data$data_weights <- weights_vec / sum(weights_vec)
  } else {
    data$data_weights <- rep(1 / nrow(data), nrow(data))
  }

  # Assemble the union of variables required by the formula, treatment,
  # weights, and any optional auxiliary or double-list components.
  vars_needed <- unique(c(
    all.vars(stats::as.formula(formula)),
    treat,
    "data_weights",
    if (spec$has_double_list) double_list,
    if (spec$has_double_list && !is.null(double_list_treat)) double_list_treat,
    if (spec$has_aux && is.null(aux_info)) c(aux_g, aux_h, sigma_h)
  ))
  vars_needed <- vars_needed[!is.na(vars_needed) & nzchar(vars_needed)]

  missing_vars <- setdiff(vars_needed, names(data))
  if (length(missing_vars) > 0L) {
    stop("The following variables are missing from data: ", paste(missing_vars, collapse = ", "), ".")
  }

  # model.frame() performs listwise deletion consistently across all variables
  # that enter the selected model.
  mf <- stats::model.frame(
    data[vars_needed],
    na.action = stats::na.omit
  )

  weights_vec <- mf$data_weights
  weight_total <- sum(weights_vec)
  if (!is.finite(weight_total) || weight_total <= 0) {
    stop("The retained observations must have positive total weight.")
  }
  weights_vec <- weights_vec / weight_total

  list(data = mf, weights = weights_vec)
}

#' Construct primary-list model data
#'
#' Extract the primary response, treatment indicator, design matrix, and weights
#' used by all outcome models.
#'
#' @keywords internal
bayesICT_make_primary_data <- function(formula, data, treat, weights) {
  # Formula::Formula preserves the standard formula interface while allowing
  # model.matrix() to construct the regression design.
  f <- Formula::Formula(formula)
  Y <- c(as.matrix(stats::model.frame(f, data)[, 1]))
  X <- stats::model.matrix(f, data=as.data.frame(unclass(data)))

  # Temporarily add the treatment variable to the formula so it is extracted
  # from the same row-aligned model frame as the response and covariates.
  fnew <- as.character(stats::as.formula(formula))
  fnew[3] <- paste(treat, "+", fnew[3])
  treat_list <- stats::model.frame(
    stats::as.formula(paste(fnew[2], fnew[1], fnew[3])),
    data
  )[, 2]

  list(
    N = length(Y),
    K = ncol(X),
    Y = Y,
    X = X,
    treat_list = treat_list,
    weights = weights
  )
}

#' Recycle scalar prior inputs safely
#'
#' Accept vectors of length one or the exact required length.  Other lengths are
#' rejected to avoid silently discarding user-supplied prior information.
#'
#' @keywords internal
bayesICT_recycle_or_stop <- function(x, n, name) {
  # Scalar hyperparameters are interpreted as common priors across dimensions.
  if (length(x) == 1L) {
    return(rep(x, n))
  }
  if (length(x) != n) {
    stop(name, " must have length 1 or length ", n, ".")
  }
  x
}

#' Prepare prior hyperparameter vectors
#'
#' Expand scalar prior hyperparameters to match the regression dimension and
#' control-item cutpoint dimension.
#'
#' @keywords internal
bayesICT_prepare_prior_vectors <- function(K,
                                           J,
                                           mu_delta,
                                           mu_psi0,
                                           mu_psiM,
                                           mu_rho0,
                                           sigma_delta,
                                           sigma_psi0,
                                           sigma_psiM,
                                           sigma_bias,
                                           sigma_rho0) {
  # delta and psi0 have one element per regression coefficient; psiM has one
  # element per non-baseline control-item category/cutpoint.
  list(
    mu_delta = bayesICT_recycle_or_stop(mu_delta, K, "mu_delta"),
    mu_psi0 = bayesICT_recycle_or_stop(mu_psi0, K, "mu_psi0"),
    mu_psiM = bayesICT_recycle_or_stop(mu_psiM, J - 1L, "mu_psiM"),
    mu_rho0 = mu_rho0,
    sigma_delta = bayesICT_recycle_or_stop(sigma_delta, K, "sigma_delta"),
    sigma_psi0 = bayesICT_recycle_or_stop(sigma_psi0, K, "sigma_psi0"),
    sigma_psiM = bayesICT_recycle_or_stop(sigma_psiM, J - 1L, "sigma_psiM"),
    sigma_bias = sigma_bias,
    sigma_rho0 = sigma_rho0
  )
}

#' Construct auxiliary-prior data
#'
#' Read auxiliary prevalence information either from a pre-built aux_info list
#' or from named variables in the model frame.
#'
#' @keywords internal
bayesICT_make_aux_data <- function(data, aux_info, aux_g, aux_h, sigma_h, N,
                                   weights) {
  # Prefer an explicit aux_info list when supplied; otherwise reconstruct the
  # same structure from named variables below.
  if (!is.null(aux_info)) {
    required <- c("G", "h", "sigma_h", "g")
    missing_aux <- setdiff(required, names(aux_info))
    if (length(missing_aux) > 0L) {
      stop("aux_info must contain G, h, sigma_h, and g. Missing: ",
           paste(missing_aux, collapse = ", "), ".")
    }

    G <- aux_info$G
    h <- aux_info$h
    sigma_h_vec <- aux_info$sigma_h
    g <- aux_info$g

    if (length(G) != 1L || !is.numeric(G) || is.na(G) || !is.finite(G) ||
        G < 1 || G != as.integer(G)) {
      stop("aux_info$G must be a positive integer.")
    }
    if (length(g) != N) {
      stop("The length of aux_info$g must match the retained data rows; check for missing values.")
    }
    if (length(h) != G || length(sigma_h_vec) != G) {
      stop("aux_info$h and aux_info$sigma_h must each contain one value per group.")
    }
  } else {
    # Each subgroup has one prevalence estimate and one uncertainty value.
    # Build the vectors by group rather than by unique values: distinct groups
    # can legitimately have the same h or sigma_h.
    g <- data[[aux_g]]
    h0 <- data[[aux_h]]
    sigma_h0 <- data[[sigma_h]]
  }

  if (is.factor(g)) {
    g <- as.character(g)
  }
  if (length(g) != N || anyNA(g)) {
    stop("Auxiliary group labels must be present for every retained observation.")
  }
  group_labels <- sort(unique(g))
  G_found <- length(group_labels)
  if (G_found < 1L || (!is.null(aux_info) && G_found != G)) {
    stop("Every auxiliary group must have at least one retained observation.")
  }
  g_index <- match(g, group_labels)

  if (is.null(aux_info)) {
    if (!is.numeric(h0) || !is.numeric(sigma_h0)) {
      stop("Auxiliary prevalence and standard deviation variables must be numeric.")
    }
    group_value <- function(values, field) {
      vapply(seq_len(G_found), function(gg) {
        one_group <- unique(values[g_index == gg])
        if (length(one_group) != 1L || is.na(one_group)) {
          stop("The auxiliary variable '", field,
               "' must have one value within each group.")
        }
        as.numeric(one_group)
      }, numeric(1))
    }
    h <- group_value(h0, aux_h)
    sigma_h_vec <- group_value(sigma_h0, sigma_h)
    G <- G_found
  }

  if (!is.numeric(h) || !is.numeric(sigma_h_vec) ||
      any(!is.finite(h)) || any(!is.finite(sigma_h_vec)) ||
      any(h <= 0 | h >= 1) || any(sigma_h_vec <= 0) ||
      any(sigma_h_vec^2 >= h * (1 - h))) {
    stop("Auxiliary h must be in (0, 1), with 0 < sigma_h^2 < h * (1 - h) for every group.")
  }
  if (length(weights) != N || any(!is.finite(weights)) ||
      any(weights < 0) ||
      any(vapply(seq_len(G), function(gg) sum(weights[g_index == gg]),
                 numeric(1)) <= 0)) {
    stop("Every auxiliary group must have positive total survey weight.")
  }

  list(G = as.integer(G), g = as.integer(g_index), h = unname(h),
       sigma_h = unname(sigma_h_vec))
}

#' Construct double-list model data
#'
#' Build the second-list response, design matrix, and treatment vector using the
#' same right-hand-side formula as the primary list.
#'
#' @keywords internal
bayesICT_make_double_list_data <- function(formula,
                                           data,
                                           double_list,
                                           double_list_treat,
                                           default_treat) {
  # If the second list has no separate treatment variable, use the primary
  # treatment variable as in the original implementation.
  double_treat_var <- if (is.null(double_list_treat)) default_treat else double_list_treat

  # Replace only the left-hand side of the formula, so the second list uses the
  # same covariate design as the first list.
  fnew <- as.character(stats::as.formula(formula))
  fnew[2] <- double_list
  f2 <- stats::as.formula(paste(fnew[2], fnew[1], fnew[3]))

  Y2 <- c(as.matrix(stats::model.frame(f2, data)[, 1]))
  X2 <- stats::model.matrix(f2, data=as.data.frame(unclass(data)))

  fnew2 <- as.character(f2)
  fnew2[3] <- paste(double_treat_var, "+", fnew2[3])
  treat2 <- stats::model.frame(
    stats::as.formula(paste(fnew2[2], fnew2[1], fnew2[3])),
    data
  )[, 2]

  list(
    N2 = length(Y2),
    Y2 = Y2,
    X2 = X2,
    treat2 = treat2,
    double_list_var = double_list,
    double_list_treat_var = double_treat_var
  )
}

#' Construct the Stan specification
#'
#' Select the compiled Stan model, build the data list, choose monitored
#' parameters, and provide default initial values if none are supplied.
#'
#' @keywords internal
bayesICT_make_stan_spec <- function(spec,
                                    J,
                                    data_info,
                                    prior_info,
                                    aux_data,
                                    double_data,
                                    init_fun) {
  # Keep all Stan-specific choices together so the main fitting function only
  # has to call rstan::sampling() once.
  stanmodel <- bayesICT_select_stanmodel(spec)
  datlist <- bayesICT_make_datlist(
    spec = spec,
    J = J,
    data_info = data_info,
    prior_info = prior_info,
    aux_data = aux_data,
    double_data = double_data
  )
  pars <- bayesICT_make_pars(spec)

  if (is.null(init_fun)) {
    init_fun <- bayesICT_make_default_init(
      spec = spec,
      J = J,
      K = data_info$K,
      prior_info = prior_info,
      pars = pars
    )
  }

  list(
    stanmodel = stanmodel,
    datlist = datlist,
    pars = pars,
    init = init_fun
  )
}

#' Select the compiled Stan model
#'
#' Map the internal model specification to the corresponding compiled Stan model
#' object stored in the package namespace.
#'
#' @keywords internal
bayesICT_select_stanmodel <- function(spec) {
  ## For non-auxiliary models, `bias` is irrelevant in the original code.
  ## Normalize here as a second line of defense so that older saved objects or
  ## direct calls to helpers with bias = "common" still select the correct model.
  bias_key <- if (spec$prior_kind %in% c("none", "double_list")) "no" else spec$bias
  key <- paste(spec$prior_kind, spec$list_variable_model, bias_key, sep = "/")

  if (spec$has_aux && spec$list_variable_model == "diff_means") {
    stop("Auxiliary difference-in-means models are not implemented; choose a different list_variable_model or omit the auxiliary prior.")
  }
  if (spec$prior_kind == "auxiliary" && spec$list_variable_model == "binomial") {
    stop("A single-list auxiliary binomial model is not implemented; choose beta-binomial or multinomial, or omit the auxiliary prior.")
  }
  if (spec$prior_kind == "auxiliary" &&
      spec$list_variable_model == "beta-binomial" && bias_key == "no") {
    stop("A single-list auxiliary beta-binomial model without bias is not implemented; use bias = 'common' or choose another model.")
  }

  # Each key has the form prior_kind/list_variable_model/bias_key.  This table
  # replaces the nested if/else tree in the original implementation.
  model_map <- list(
    "none/beta-binomial/no" = stanmodels$model_ICT_noaux,
    "none/binomial/no" = stanmodels$model_ICT_binomial_noaux,
    "none/multinomial/no" = stanmodels$model_ICT_multinomial_noaux,
    "none/diff_means/no" = stanmodels$model_ICT_multinomial_diff_means,

    "auxiliary/beta-binomial/common" = stanmodels$model_ICT_aux,
    "auxiliary/beta-binomial/separate" = stanmodels$model_ICT_aux,
    "auxiliary/multinomial/common" = stanmodels$model_ICT_multinomial_aux,
    "auxiliary/multinomial/separate" = stanmodels$model_ICT_multinomial_aux, # no seperate for single list
    "auxiliary/multinomial/no" = stanmodels$model_ICT_multinomial_aux_nobias,
    "double_list/beta-binomial/no" = stanmodels$model_ICT_doublelist_noaux,
    "double_list/binomial/no" = stanmodels$model_ICT_doublelist_binomial_noaux,
    "double_list/multinomial/no" = stanmodels$model_ICT_doublelist_multinomial_noaux,
    "double_list/diff_means/no" = stanmodels$model_ICT_doublelist_multinomial_diff_means,

    "double_list_aux/beta-binomial/common" = stanmodels$model_ICT_doublelist_aux,
    "double_list_aux/beta-binomial/separate" = stanmodels$model_ICT_doublelist_aux_bias_separate,
    "double_list_aux/beta-binomial/no" = stanmodels$model_ICT_doublelist_aux_nobias,
    "double_list_aux/binomial/common" = stanmodels$model_ICT_doublelist_binomial_aux_bias_common,
    "double_list_aux/binomial/separate" = stanmodels$model_ICT_doublelist_binomial_aux_bias_separate,
    "double_list_aux/binomial/no" = stanmodels$model_ICT_doublelist_binomial_aux_nobias,
    "double_list_aux/multinomial/common" = stanmodels$model_ICT_doublelist_multinomial_aux_bias_common,
    "double_list_aux/multinomial/separate" = stanmodels$model_ICT_doublelist_multinomial_aux_bias_separate, 
    "double_list_aux/multinomial/no" = stanmodels$model_ICT_doublelist_multinomial_aux_nobias
  )

  stanmodel <- model_map[[key]]
  if (is.null(stanmodel)) {
    stop("Unsupported combination of prior, list_variable_model, and bias: ", key, ".")
  }
  stanmodel
}

#' Build the Stan data list
#'
#' Combine primary-list, double-list, auxiliary-prior, and prior-hyperparameter
#' components into the list expected by the selected Stan program.
#'
#' @keywords internal
bayesICT_make_datlist <- function(spec, J, data_info, prior_info, aux_data, double_data) {
  # Start with the common fields for either single-list or double-list Stan
  # programs, then append family- and prior-specific fields below.
  if (spec$has_double_list) {
    datlist <- list(
      N = double_data$N2,
      J = J,
      Y1 = data_info$Y,
      Y2 = double_data$Y2,
      K = data_info$K,
      X = data_info$X,
      X2 = double_data$X2,
      treat = data_info$treat_list,
      treat2 = double_data$treat2,
      weights = data_info$weights,
      mu_psi0 = as.array(prior_info$mu_psi0),
      mu_delta = as.array(prior_info$mu_delta),
      sigma_psi0 = as.array(prior_info$sigma_psi0),
      sigma_delta = as.array(prior_info$sigma_delta)
    )
  } else {
    datlist <- list(
      N = data_info$N,
      J = J,
      Y = data_info$Y,
      K = data_info$K,
      X = data_info$X,
      treat = data_info$treat_list,
      weights = data_info$weights,
      sigma_psi0 = as.array(prior_info$sigma_psi0),
      sigma_delta = as.array(prior_info$sigma_delta)
    )
  }

  # No-prior single-list models need only the likelihood-family-specific prior
  # hyperparameters used by their Stan programs.
  if (spec$prior_kind == "none" && !spec$has_double_list) {
    if (spec$list_variable_model == "beta-binomial") {
      datlist$mu_rho0 <- prior_info$mu_rho0
      datlist$sigma_rho0 <- prior_info$sigma_rho0 /
        (prior_info$mu_rho0 * (1 - prior_info$mu_rho0))
    }
    if (spec$list_variable_model == "binomial") {
      datlist$mu_rho0 <- prior_info$mu_rho0
      datlist$sigma_rho0 <- prior_info$sigma_rho0 /
        (prior_info$mu_rho0 * (1 - prior_info$mu_rho0))
    }
    if (spec$list_variable_model %in% c("multinomial", "diff_means")) {
      datlist$mu_psiM <- prior_info$mu_psiM
      datlist$sigma_psiM <- as.array(prior_info$sigma_psiM)
    }
    if (spec$list_variable_model == "diff_means") {
      datlist$sigma_psiMT <- as.array(rep(prior_info$sigma_psiM[1], J))
    }
  }

  # Auxiliary models add group-level prevalence information and, when selected,
  # bias-prior scale information.
  if (spec$prior_kind == "auxiliary") {
    datlist$mu_delta <- as.array(prior_info$mu_delta)
    datlist$mu_psiM <- as.array(prior_info$mu_psiM)
    datlist$mu_psi0 <- as.array(prior_info$mu_psi0)
    datlist$sigma_psiM <- as.array(prior_info$sigma_psiM)
    datlist$sigma_bias <- prior_info$sigma_bias
    datlist$G <- aux_data$G
    datlist$h <- as.array(aux_data$h)
    datlist$sigma_h <- as.array(aux_data$sigma_h)
    datlist$g <- aux_data$g
  }

  # Double-list models without auxiliary information only need the second-list
  # data plus any multinomial cutpoint priors.
  if (spec$prior_kind == "double_list") {
    if (spec$list_variable_model %in% c("multinomial", "diff_means")) {
      datlist$sigma_psiM <- as.array(prior_info$sigma_psiM)
      datlist$mu_psiM <- prior_info$mu_psiM
    }
    if (spec$list_variable_model == "diff_means") {
      datlist$sigma_psiMT <- as.array(rep(prior_info$sigma_psiM[1], J))
    }
  }

  # Double-list auxiliary models combine the double-list fields with the
  # group-level auxiliary-prior fields.
  if (spec$prior_kind == "double_list_aux") {
    datlist$sigma_bias <- prior_info$sigma_bias
    datlist$G <- aux_data$G
    datlist$h <- as.array(aux_data$h)
    datlist$sigma_h <- as.array(aux_data$sigma_h)
    datlist$g <- aux_data$g

    if (spec$list_variable_model == "beta-binomial") {
      datlist$mu_rho0 <- prior_info$mu_rho0
      datlist$sigma_rho0 <- prior_info$sigma_rho0
    } else if (spec$list_variable_model == "multinomial") {
      datlist$mu_psiM <- prior_info$mu_psiM
      datlist$sigma_psiM <- prior_info$sigma_psiM
    }
  }

  datlist
}

#' Choose monitored Stan parameters
#'
#' Return the parameter names passed to rstan::sampling() for the selected model
#' specification.
#'
#' @keywords internal
bayesICT_make_pars <- function(spec) {
  # The monitored parameters must match both the selected Stan program and the
  # quantities expected by print(), plot(), and predict() methods.
  if (spec$prior_kind == "none") {
    if (spec$list_variable_model == "multinomial") {
      return(c("psi0", "delta", "psiM", "pas", "log_lik"))
    }
    if (spec$list_variable_model == "diff_means") {
      return(c("psi0", "delta", "psiM", "psiMT", "pas", "log_lik"))
    }
    if (spec$list_variable_model == "beta-binomial") {
      return(c("psi0", "delta", "rho0", "pas", "log_lik"))
    }
    return(c("psi0", "delta", "pas", "log_lik"))
  }

  if (spec$prior_kind == "auxiliary") {
    return(c(
      "psi0", "delta",
      if (spec$list_variable_model == "beta-binomial") "rho0",
      if (spec$list_variable_model == "multinomial") "psiM",
      if (spec$bias %in% c("common", "separate")) "bias",
      "pas", "log_lik"
    ))
  }

  if (spec$prior_kind == "double_list") {
    if (spec$list_variable_model == "multinomial") {
      return(c("psi0", "delta", "psi02", "delta2", "psiM", "psiM2", "pas", "log_lik"))
    }
    if (spec$list_variable_model == "binomial") {
      return(c("psi0", "delta", "psi02", "delta2", "pas", "log_lik"))
    }
    if (spec$list_variable_model == "diff_means") {
      return(c("psi0", "delta", "psiM", "psiMT", "psi02", "delta2",
               "psiM2", "psiMT2", "pas", "log_lik"))
    }
    return(c("psi0", "delta", "psi02", "delta2", "rho0", "rho02", "pas", "log_lik"))
  }

  if (spec$prior_kind == "double_list_aux") {
    return(c(
      "psi0", "delta", "psi02", "delta2",
      if (spec$list_variable_model == "beta-binomial") c("rho0", "rho02"),
      if (spec$list_variable_model == "multinomial") "psiM",
      if (!(spec$list_variable_model %in% c("beta-binomial", "binomial"))) "psiM2",
      if (spec$bias %in% c("common", "separate")) "bias",
      if (spec$bias == "separate") "bias2",
      "pas", "log_lik"
    ))
  }

  stop("Unsupported prior kind: ", spec$prior_kind, ".")
}

#' Construct default initial values
#'
#' Return a function that initializes only the parameters monitored for the
#' selected model.
#'
#' @keywords internal
bayesICT_make_default_init <- function(spec, J, K, prior_info, pars) {
  # Capture prior centers in the closure returned below.  rstan will call the
  # returned function once per chain.
  mu_psi0 <- prior_info$mu_psi0
  mu_delta <- prior_info$mu_delta
  mu_psiM <- prior_info$mu_psiM

  function() {
    # Initialize only parameters present in the selected Stan program.
    init <- list()

    if ("psi0" %in% pars) {
      init$psi0 <- as.array(mu_psi0 + stats::rnorm(n = K, mean = 0, sd = 0.1))
    }
    if ("delta" %in% pars) {
      init$delta <- as.array(mu_delta)
    }
    if ("psi02" %in% pars) {
      init$psi02 <- as.array(mu_psi0 + stats::rnorm(n = K, mean = 0, sd = 0.1))
    }
    if ("delta2" %in% pars) {
      init$delta2 <- as.array(mu_delta)
    }

    if ("rho0" %in% pars) {
      ## The original single-list no-prior beta-binomial code initialized
      ## lrho0 while reporting rho0.  Preserve that convention for that
      ## case; initialize rho0 directly elsewhere.
      if (spec$prior_kind == "none" && !spec$has_double_list) {
        init$lrho0 <- 0
      } else {
        init$rho0 <- stats::runif(n = 1)
      }
    }
    if ("rho02" %in% pars) {
      init$rho02 <- stats::runif(n = 1)
    }

    if ("psiM" %in% pars) {
      init$psiM <- as.array(sort(mu_psiM + stats::rnorm(n = J - 1L, mean = 0, sd = 0.1)))
    }
    if ("psiM2" %in% pars) {
      init$psiM2 <- as.array(sort(mu_psiM + stats::rnorm(n = J - 1L, mean = 0, sd = 0.1)))
    }
    if ("psiMT" %in% pars) {
      init$psiMT <- as.array(sort(stats::rnorm(n = J, mean = 0, sd = 0.1)))
    }
    if ("psiMT2" %in% pars) {
      init$psiMT2 <- as.array(sort(stats::rnorm(n = J, mean = 0, sd = 0.1)))
    }

    init
  }
}
