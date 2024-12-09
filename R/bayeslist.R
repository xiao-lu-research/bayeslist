#' Fitting Bayesian sensitive item models
#'
#' The main function for estimating Bayesian sensitive item models.
#' The function returns a \code{bayeslist} object that can be further investigated using standard functions such as \code{summary}, \code{plot}, \code{print}, \code{predict}, and \code{coef}. The model can be passed using a \code{formula} as in \code{lm()}. Convergence diagnotics can be performed using either \code{print(object, "mcmc")} or \code{plot(object, "trace")}.
#'
#' @param formula An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted.
#' @param data A data frame containing the variables in the model.
#' @param treat Variable name of the treatment.
#' @param J Number of control items.
#' @param type Type of the model. Options include "outcome", "predict", "misreport", for the sensitive item outcome model, predictor model and misreport model, respectively.
#' @param nsim The number of iterations.
#' @param burnin The number of burnin iterations.
#' @param thin Thinning parameter.
#' @param CIsize The size of posterior confidence interval.
#' @param nchain The number of parallel chains.
#' @param seeds Random seeds to replicate the results.
#' @param vb Logic. If TRUE, variational approximation will be used to supply initial values. The default is FALSE.
#' @param only_vb Logic. If TRUE, only variational approximation will be calculated. The default is FALSE.
#' @param prior Prior types. Options include "auxiliary", "double_list", "direct_item", and "BL" for beta-logistic prior. If NULL, no informative priors will be used.
#' @param direct_item Variable name of the direct item.
#' @param direct_item_misreport Variable name of the direct item for the misreporting model.
#' @param double_list Variable name of the second list.
#' @param double_list_treat Treatment variable of the second list.
#' @param aux_info Auxiliary information for the informative priors. list(G,h,g), where: G (number of subgroups), h (auxiliary information for each subgroup), and g (subgroup indicator). If is.NULL, the following two parameters need to be specified when estimating the model with prior = "auxiliary".
#' @param aux_g Auxiliary information for the informative priors: name of the variable indicating the group of each observation.
#' @param aux_h Auxiliary information for the informative priors: name of the variable containing information of prevalence for each group
#' @param BL_a The first shape hyperparameter for the beta-logistic prior, indicating the prior number of affirmative answers to the sensitive item.
#' @param BL_b The second shape hyperparameter for the beta-logistic prior, indicating the prior number of non-affirmative answers to the sensitive item.
#' @param conjugate_distance Logic. Indicating whether conjugate distance prior should be used. The default is FALSE.
#' @param conjugate_k Degrees of freedom to be scaled by conjugate distance prior. The default is NULL.
#' @param predictvar Variable name of the outcome to be predicted.
#' @param predictvar_type The type of the outcome variable to be predicted. Options include "linear" and "binary". The default is "binary".
#' @param parallel Logic. Indicating whether to do paralell computing. The default is TRUE.
#' @param robust Logic. Indicating whether to impose robust constraints on the intercept-only model. The default is FALSE.
#'
#' @return A \code{bayeslist} object. An object of class \code{bayeslist} contains the following elements
#'
#'   \describe{
#'
#'   \item{\code{Call}}{The matched call.}
#'   \item{\code{formula}}{Symbolic representation of the model.}
#'   \item{\code{type}}{Model type}
#'   \item{\code{nsim}}{Number of iterations.}
#'   \item{\code{Burnin}}{Number of burnin iterations.}
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
#'   \item{\code{only_vb}}{Whether only viariational approximation is used.}
#'   \item{\code{prior}}{Informative prior types.}
#'   \item{\code{direct_item}}{Direct item.}
#'   \item{\code{double_list}}{The second list.}
#'   \item{\code{aux_info}}{Auxiliary information.}
#'   \item{\code{ulbs}}{Upper and lower bounds based on the specified confidence interval.}
#'   \item{\code{means}}{Mean estimates.}
#'   \item{\code{treat}}{Treatment.}
#'   \item{\code{outcome}}{Outcome to be predicted.}
#'   \item{\code{direct}}{Direct item for the misreport model.}
#'   \item{\code{robust}}{Robust indicator.}
#'
#' }
#'
#' @references Lu, X. and Traunmüller, R. (2021). Improving Studies of Sensitive Topics Using Prior Evidence: A Unified Bayesian Framework for List Experiments, SSRN, \doi{10.2139/ssrn.3871089}.
#'
#' @import Formula
#' @import rstan
#' @import rstantools
#'
#' @importFrom stats coef model.frame model.matrix quantile as.formula binomial density glm sd
#'
#' @export
#'
#' @examples
#'# Estimate sensitive item outcome model using Sri Lanka data on male sexual violence
#'# Load Sri Lanka list experiment data
#'data(srilanka)
#'
#'# Model 1: intercept-only outcome model without prior information:
#'mod1 <- bayeslist(sexaussault ~ 1, data = srilanka, treat = "treatment", J = 3,
#'type = "outcome", nsim = 200, thin = 1, CIsize = 0.95, nchain = 1,
#'seeds = 342321, prior = NULL, parallel = TRUE)
#'summary(mod1) # summary of estimates
#'predict(mod1) # predicted prevalence for each observation
#'plot(mod1,"trace") # trace plot
#'plot(mod1,"coef") # coefficient plot
#'plot(mod1, only_prev = TRUE) # prevalence plot
#'
#'\donttest{
#'# Model 2: multivariate outcome model without prior information:
#'mod2 <- bayeslist(sexaussault ~ age + edu, data = srilanka, treat = "treatment", J = 3,
#'type = "outcome", nsim = 200, thin = 1, CIsize = 0.95, nchain = 1,
#'seeds = 342321, prior = NULL, parallel = TRUE)
#'summary(mod2) # summary of estimates
#'predict(mod2) # predicted prevalence for each observation
#'plot(mod2,"trace") # trace plot
#'plot(mod2,"coef") # coefficient plot
#'plot(mod2) # prevalence + coefficient plot
#'
#'# Model 3: intercept-only outcome model with prior information from medicolegal reports, i.e.,
#'# with a prior beta-logistic distribution BL(38, 146).
#' a <- 38; b <-146
#'mod3 <- bayeslist(sexaussault ~ 1, data = srilanka, treat = "treatment", J = 3,
#'type = "outcome", nsim = 200, thin = 1, CIsize = 0.95, nchain = 1,
#'seeds = 342321, prior = "BL", BL_a = a, BL_b = b,, parallel = TRUE)
#'summary(mod3)
#'predict(mod3)
#'plot(mod3,"trace")
#'plot(mod3,"coef")
#'plot(mod3, only_prev = TRUE)
#'
#'# Model 4: multivariate outcome model with prior information from a direct item.
#'# Load London list experiment data
#'data(london)
#'mod4 <- bayeslist(listCount ~ agegrp + gender + social_grade + qual,data = london, J = 4,
#'treat = "listTreat", seeds = 4597, nsim = 200, nchain = 1,
#'prior = "direct_item", direct_item = "baselineTurnout")
#'summary(mod4)
#'predict(mod4)
#'plot(mod4,"trace")
#'plot(mod4,"coef")
#'plot(mod4)
#'}
#'
bayeslist <- function(formula,
                      data,
                      treat,
                      J,
                      type = "outcome", # "outcome", "predict", or "misreport"
                      nsim = 1000,
                      burnin = NULL,
                      thin = 1,
                      CIsize = .95,
                      nchain = 1,
                      seeds = 12345,
                      vb = FALSE, # if true, variation inference will be conducted to initialize the chains
                      only_vb = FALSE, # if true, only do variational approximation without direct sampling
                      prior = NULL, # options: 1. "auxiliary", 2. "double_list" 3. "direct_item" 4. "BL" 5. NULL
                      direct_item = NULL, # variable name of the direct item
                      direct_item_misreport = NULL, # variable name of the direct item for the misreport model
                      double_list = NULL, # variable name for the outcome of the second list
                      double_list_treat = NULL, # treatment variable for the second list
                      aux_info = NULL, # list(G,h,g), where: G (number of subgroups), h (auxiliary information for each subgroup), and g (subgroup indicator).
                      aux_g = NULL, # name of the variable indicating the group of each observation
                      aux_h = NULL, # name of the variable containing information of prevalence for each group
                      BL_a = NULL, # first parameter for beta-logistic prior: prior number of affirmative answers
                      BL_b = NULL, # second parameter for beta-logistic prior: prior number of non-affirmative answers
                      conjugate_distance = FALSE, # Indicating whether conjugate distance prior will be used for the direct item. The default is FALSE.
                      conjugate_k = NULL, # degrees of freedom to be scaled by conjugate distance prior. The default is NULL: no scaling.
                      predictvar = NULL, # outcome variable to be predicted in the predict model
                      predictvar_type = "binary", # outcome type for the predictor model: either "binary" or "linear". The default is binary.
                      # direct = NULL, # data of direct item for the misreport model
                      parallel = TRUE,
                      robust = FALSE
) {
  if (is.null(burnin))
    burnin <- floor(nsim / 2)
  if (burnin < 0)
    stop("Burn-in must be non-negative.")
  if (thin < 1)
    stop("Thinning factor must be positive.")
  if (CIsize <= 0)
    stop("Confidence interval size 'CIsize' must be positive.")
  if (CIsize > 1)
    stop(paste0("Confidence interval size 'CIsize' ",
                "can not be larger than 1."))
  if (missing(formula) | missing(data)) {
    stop(paste0("Formula and data should be given."))
  }
  # Check prior information
  if (is.null(prior) == FALSE){
    if (prior == "auxiliary" & is.null(aux_info)){
      if (is.null(aux_g)){
        stop("No auxiliary information specified!")
      } else if (is.null(aux_h)){
        stop("No auxiliary information specified!")
      }
    }
    if (prior == "direct_item" && is.null(direct_item)){
      stop("No direct item indicated!")
    }
    if (prior == "double_list" && is.null(double_list)){
      stop("No double list item indicated!")
    }
  }
  # list-wise deletion for missing values
  if (type == "outcome") {
    if (is.null(prior)){
      data = data
    } else if (prior == "direct_item") {
      fnew <- as.character(as.formula(formula))
      fnew <- as.formula(paste(fnew[2],fnew[1],fnew[3],"+",paste(direct_item,collapse = "+"),"+",treat))
      data <- model.frame(fnew,data)
    } else if (prior == "auxiliary") {
      if (!is.null(aux_g) & !is.null(aux_h)){
        fnew <- as.character(as.formula(formula))
        fnew <- as.formula(paste(fnew[2],fnew[1],fnew[3],"+",aux_g, "+", aux_h,"+",treat))
        data <- model.frame(fnew,data)
      }
    } else if (prior == "double_list") {
      fnew <- as.character(as.formula(formula))
      fnew <- as.formula(paste(fnew[2],fnew[1],fnew[3],"+",double_list,"+",double_list_treat,"+",treat))
      data <- model.frame(fnew,data)
    } else {
      data = data
    }
  } else if (type == "predict") {
    if (is.null(prior)){
      fnew <- as.character(as.formula(formula))
      fnew <- as.formula(paste(fnew[2],fnew[1],fnew[3],"+",predictvar,"+",treat))
      data <- model.frame(fnew,data)
    } else if (prior == "direct_item") {
      fnew <- as.character(as.formula(formula))
      fnew <- as.formula(paste(fnew[2],fnew[1],fnew[3],"+",predictvar,"+",paste(direct_item,collapse = "+"),"+",treat))
      data <- model.frame(fnew,data)
    } else if (prior == "auxiliary") {
      if (!is.null(aux_g) & !is.null(aux_h)){
        fnew <- as.character(as.formula(formula))
        fnew <- as.formula(paste(fnew[2],fnew[1],fnew[3],"+",predictvar,"+",aux_g, "+", aux_h,"+",treat))
        data <- model.frame(fnew,data)
      }
    } else if (prior == "double_list") {
      fnew <- as.character(as.formula(formula))
      fnew <- as.formula(paste(fnew[2],fnew[1],fnew[3],"+",predictvar,"+",double_list,"+",double_list_treat,"+",treat))
      data <- model.frame(fnew,data)
    } else {
      data = data
    }
  } else if (type == "misreport") {
    if (is.null(prior)){
      fnew <- as.character(as.formula(formula))
      fnew <- as.formula(paste(fnew[2],fnew[1],fnew[3],"+",direct_item_misreport,"+",treat))
      data <- model.frame(fnew,data)
    } else if (prior == "direct_item") {
      fnew <- as.character(as.formula(formula))
      fnew <- as.formula(paste(fnew[2],fnew[1],fnew[3],"+",direct_item_misreport,"+",paste(direct_item,collapse = "+"),"+",treat))
      data <- model.frame(fnew,data)
    } else if (prior == "auxiliary") {
      if (!is.null(aux_g) & !is.null(aux_h)){
        fnew <- as.character(as.formula(formula))
        fnew <- as.formula(paste(fnew[2],fnew[1],fnew[3],"+",direct_item_misreport,"+",aux_g, "+", aux_h,"+",treat))
        data <- model.frame(fnew,data)
      }
    } else if (prior == "double_list") {
      fnew <- as.character(as.formula(formula))
      fnew <- as.formula(paste(fnew[2],fnew[1],fnew[3],"+",direct_item_misreport,"+",double_list,"+",double_list_treat,"+",treat))
      data <- model.frame(fnew,data)
    } else if (prior == "BL") {
      fnew <- as.character(as.formula(formula))
      fnew <- as.formula(paste(fnew[2],fnew[1],fnew[3],"+",direct_item_misreport,"+",treat))
      data <- model.frame(fnew,data)
    } else {
      data = data
    }
  } else {
    stop("Unknown type of sensitive item models!")
  }
  # check number of observations
  if (dim(data)[1] < 5) {
    warning("Less than 5 observations after list-wise deletion. Please check your data and missing values.")
  }
  # parallel computing
  if (parallel == TRUE){
    options(mc.cores = parallel::detectCores())
  }
  ##########################################################
  # Data info
  data0 <- data
  treat0 <- ifelse(is.null(double_list_treat),treat,double_list_treat)
  f <- Formula::Formula(formula)
  data <- model.frame(f, data)
  Y <- c(as.matrix(model.frame(f, data)[, 1]))
  X <- model.matrix(f, data)
  fnew <- as.character(as.formula(formula))
  fnew[3] <- paste(treat,"+",fnew[3])
  treat <- model.frame(as.formula(paste(fnew[2],fnew[1],fnew[3])),data0)[,2]


  n_covariate <- dim(X)[2]
  N <- length(Y)
  K <- dim(X)[2]
  # aux_info <- NULL

  ########################################################################
  # Bayesian estimation with informative priors
  if (is.null(prior) == FALSE) {
    if (prior == "BL"){
      ###############################################
      # beta-logistic prior
      ###############################################
      if (type == "outcome") {
        # model_outcome_noaux = stan_model("list_outcome_constrained1.0.stan")
        stanmodel <- stanmodels$model_outcome_cmp
        datlist <- list(
          N = N, # obs
          J = J, # number of control items
          Y = Y, # number of affirmative answers
          K = K, # number of covariates
          X = X, # covariate matrix
          treat = treat, # treatment indicator
          a = as.array(BL_a),
          b = as.array(BL_b)
        )
      } else if (type == "predict") {
        # model_predict_noaux = stan_model("list_predictor_constrained_logit1.0.stan")
        if (predictvar_type == "binary"){
          stanmodel <- stanmodels$model_predict_cmp
        } else {
          stanmodel <- stanmodels$model_predict_cmp_linear
        }
        fnew <- as.character(as.formula(formula))
        fnew[2] <- predictvar
        fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
        outcome <- c(as.matrix(model.frame(fnew, data0)[, 1]))

        datlist <- list(
          N = N, # obs
          J = J, # number of control items
          Y = Y, # number of affirmative answers
          K = K, # number of covariates
          X = X, # covariate matrix
          treat = treat, # treatment indicator
          outcome = outcome, # outcome to be predicted
          a = as.array(BL_a),
          b = as.array(BL_b)
        )
      } else {
        # model_misreport_noaux = stan_model("list_misreport_constrained1.0.stan")
        stanmodel <- stanmodels$model_misreport_cmp

        fnew <- as.character(as.formula(formula))
        fnew[2] <- direct_item_misreport
        fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
        direct <- c(as.matrix(model.frame(fnew, data0)[, 1]))

        datlist <- list(
          N = N, # obs
          J = J, # number of control items
          Y = Y, # number of affirmative answers
          K = K, # number of covariates
          X = X, # covariate matrix
          treat = treat, # treatment indicator
          direct = direct, # direct item
          a = as.array(BL_a),
          b = as.array(BL_b)
        )
      }

      if (type == "outcome") {
        pars = c("psi0","delta","rho0")
      } else if (type == "predict") {
        pars = c("psi0","delta","rho0","psi2","gamma0","phi")
      } else {
        pars = c("psi0", "delta", "rho0","gamma0","treat_e","U_e","Z_e")
      }

      if (vb == TRUE & only_vb == FALSE){
        stanvb <- try(vb(
          stanmodel,
          data = datlist,
          pars = pars,
          seed = seeds
        ),silent = T)

        if (inherits(stanvb, "try-error")){
          warning("Variational inference fails! Proceed to direct sampling.")
        } else{
          initvalues = summary(stanvb)$summary[,1]
        }

        if (type == "outcome") {
          init_fun = function(){
            list(psi0 = initvalues[1:K],
                 delta = initvalues[(K+1):(2*K)],
                 rho0 = initvalues[2*K+1])
          }
        } else if (type == "predict") {
          init_fun = function(){
            list(psi0 = initvalues[1:K],
                 delta = initvalues[(K+1):(2*K)],
                 rho0 = initvalues[2*K+1],
                 psi2 = initvalues[(2*K+2):(3*K+1)],
                 gamma0 = initvalues[(3*K+2)],
                 phi = initvalues[(3*K+3)]
            )
          }
        } else {
          init_fun = function(){
            list(psi0 = initvalues[1:K],
                 delta = initvalues[(K+1):(2*K)],
                 rho0 = initvalues[2*K+1],
                 gamma0 = initvalues[(2*K+2):(3*K+1)],
                 treat_e = initvalues[3*K+2],
                 U_e = initvalues[3*K+3],
                 Z_e = initvalues[3*K+4]
            )
          }
        }

        if (inherits(stanvb, "try-error")){
          stanout <- sampling(
            stanmodel,
            data = datlist,
            pars = pars,
            seed = seeds,
            iter = nsim,
            thin = thin,
            warmup = burnin,
            chains = nchain
          )
        } else {
          stanout <- sampling(
            stanmodel,
            data = datlist,
            pars = pars,
            init = init_fun,
            seed = seeds,
            iter = nsim,
            thin = thin,
            warmup = burnin,
            chains = nchain
          )
        }
      } else if (only_vb == TRUE) {
        stanvb <- try(vb(
          stanmodel,
          data = datlist,
          pars = pars,
          seed = seeds
        ),silent = T)

        if (inherits(stanvb, "try-error")){
          stop("Variational inference failed!")
        }
        stanout <- stanvb
      } else {
        stanout <- sampling(
          stanmodel,
          data = datlist,
          pars = pars,
          seed = seeds,
          iter = nsim,
          thin = thin,
          warmup = burnin,
          chains = nchain
        )
      }
    } else if (prior == "auxiliary") {
      ##############################################
      # I. Auxiliary information
      ##############################################
      if (is.null(aux_info)){
        fnew <- as.character(as.formula(formula))
        fnew[3] <- paste(aux_g,"+",aux_h,"+",fnew[3])
        g <- model.frame(as.formula(paste(fnew[2],fnew[1],fnew[3])),data0)[,2]
        h <- model.frame(as.formula(paste(fnew[2],fnew[1],fnew[3])),data0)[,3]
        h <- h[order(g)]
        h <- unique(h)
        h <- h[which(!is.na(h))]
        g0 <- unique(g)
        G <- length(g0[which(!is.na(g0))])
        aux_info <- list(G=G,g=g,h=h)
      } else {
        G = aux_info$G
        h = aux_info$h
        g = aux_info$g
        if (length(g) != N){
          stop("Dimension found in g does not match that of the data. The data may contain missing values.")
        }
      }

      if (type == "outcome") {
        # model_outcome_aux = stan_model("list_outcome_constrained_aux3.0_autosigma.stan")
        stanmodel <- stanmodels$model_outcome_aux
        datlist <- list(
          N = N, # obs
          J = J, # number of control items
          Y = Y, # number of affirmative answers
          K = K, # number of covariates
          X = X, # covariate matrix
          treat = treat, # treatment indicator
          G = G,
          h = as.array(h),
          g = g
        )
      } else if (type == "predict") {
        # model_predict_aux = stan_model("list_predictor_constrained_logit1.0_aux.stan")
        if (predictvar_type == "binary"){
          stanmodel <- stanmodels$model_predict_aux
        } else {
          stanmodel <- stanmodels$model_predict_aux_linear
        }
        fnew <- as.character(as.formula(formula))
        fnew[2] <- predictvar
        fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
        outcome <- c(as.matrix(model.frame(fnew, data0)[, 1]))

        datlist <- list(
          N = N, # obs
          J = J, # number of control items
          Y = Y, # number of affirmative answers
          K = K, # number of covariates
          X = X, # covariate matrix
          treat = treat, # treatment indicator
          outcome = outcome, # outcome to be predicted
          G = G,
          h = as.array(h),
          g = g
        )
      } else {
        # model_misreport_aux = stan_model("list_misreport_constrained1.0_aux.stan")
        stanmodel <- stanmodels$model_misreport_aux
        fnew <- as.character(as.formula(formula))
        fnew[2] <- direct_item_misreport
        fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
        direct <- c(as.matrix(model.frame(fnew, data0)[, 1]))

        datlist <- list(
          N = N, # obs
          J = J, # number of control items
          Y = Y, # number of affirmative answers
          K = K, # number of covariates
          X = X, # covariate matrix
          treat = treat, # treatment indicator
          direct = direct, # direct item
          G = G,
          h = as.array(h),
          g = g
        )
      }

      #######################################
      # variables to be stored
      if (type == "outcome") {
        pars = c("psi0","delta","rho0","sigma")
      } else if (type == "predict") {
        pars = c("psi0","delta","rho0","psi2","gamma0","phi","sigma")
      } else {
        pars = c("psi0", "delta", "rho0","gamma0","treat_e","U_e","Z_e","sigma")
      }

      ################################################
      # Variational inference to get start values
      if (vb == TRUE & only_vb == FALSE){
        stanvb <- try(vb(
          stanmodel,
          data = datlist,
          pars = pars,
          seed = seeds
        ),silent = T)

        if (inherits(stanvb, "try-error")){
          warning("Variational inference fails! Proceed to direct sampling.")
        } else{
          initvalues = summary(stanvb)$summary[,1]
        }

        if (type == "outcome") {
          init_fun = function(){
            list(psi0 = initvalues[1:K],
                 delta = initvalues[(K+1):(2*K)],
                 rho0 = initvalues[2*K+1],
                 sigma = initvalues[2*K+2])
          }
        } else if (type == "predict") {
          init_fun = function(){
            list(psi0 = initvalues[1:K],
                 delta = initvalues[(K+1):(2*K)],
                 rho0 = initvalues[2*K+1],
                 psi2 = initvalues[(2*K+2):(3*K+1)],
                 gamma0 = initvalues[(3*K+2)],
                 phi = initvalues[(3*K+3)],
                 sigma = initvalues[(3*K+4)]
            )
          }
        } else {
          init_fun = function(){
            list(psi0 = initvalues[1:K],
                 delta = initvalues[(K+1):(2*K)],
                 rho0 = initvalues[2*K+1],
                 gamma0 = initvalues[(2*K+2):(3*K+1)],
                 treat_e = initvalues[3*K+2],
                 U_e = initvalues[3*K+3],
                 Z_e = initvalues[3*K+4],
                 sigma = initvalues[3*K+5]
            )
          }
        }

        ####################################################
        # direct sampling when variational inference fails
        if (inherits(stanvb, "try-error")){
          stanout <- sampling(
            stanmodel,
            data = datlist,
            pars = pars,
            seed = seeds,
            iter = nsim,
            thin = thin,
            warmup = burnin,
            chains = nchain
          )
        } else {
          stanout <- sampling(
            stanmodel,
            data = datlist,
            pars = pars,
            init = init_fun,
            seed = seeds,
            iter = nsim,
            thin = thin,
            warmup = burnin,
            chains = nchain
          )
        }
      } else if (only_vb == TRUE) {
        stanvb <- try(vb(
          stanmodel,
          data = datlist,
          pars = pars,
          seed = seeds
        ),silent = T)

        if (inherits(stanvb, "try-error")){
          stop("Variational inference failed!")
        }
        stanout <- stanvb
      } else {
        stanout <- sampling(
          stanmodel,
          data = datlist,
          pars = pars,
          seed = seeds,
          iter = nsim,
          thin = thin,
          warmup = burnin,
          chains = nchain
        )
      }
      ##############################################
      # II. Direct item
    } else if (prior == "direct_item"){
      ivs = as.character(as.formula(formula))[3]# strsplit(formula,split = "~")[[1]][2]
      if (length(direct_item) == 1){
        fnew = paste(direct_item, "~", ivs)
        logitout = glm(fnew , data = data0, family = binomial(link = "logit"))
        logitoutsum = summary(logitout)
        mu_delta = logitoutsum$coefficients[,1]
        sigma_delta = logitoutsum$coefficients[,2]
        mu_delta00 = mu_delta
        sigma_delta00 = sigma_delta
      } else {
        direct_vec = NA
        for (i in 1:length(direct_item)){
          direct_vec = c(direct_vec, model.matrix(as.formula(paste("~",direct_item[i])),data0)[,2])
        }
        direct_vec = direct_vec[-1]
        data0_stack = data0[rep(seq_len(nrow(data0)), length(direct_item)),]
        data0_stack$direct_vec00000000 = direct_vec
        fnew = paste("direct_vec00000000", "~", ivs)
        logitout = glm(fnew , data = data0_stack, family = binomial(link = "logit"))
        logitoutsum = summary(logitout)
        mu_delta = logitoutsum$coefficients[,1]
        sigma_delta = logitoutsum$coefficients[,2]
        mu_delta00 = mu_delta
        sigma_delta00 = sigma_delta
      }

      if (conjugate_distance == TRUE){
        ###############################################
        # direct item with conjugate distance
        ###############################################
        # First get unbiased estimates
        if (type == "outcome") {
          # model_outcome_noaux = stan_model("list_outcome_constrained1.0.stan")
          stanmodel <- stanmodels$model_outcome_noaux
          datlist <- list(
            N = N, # obs
            J = J, # number of control items
            Y = Y, # number of affirmative answers
            K = K, # number of covariates
            X = X, # covariate matrix
            treat = treat # treatment indicator
          )
        } else if (type == "predict") {
          # model_predict_noaux = stan_model("list_predictor_constrained_logit1.0.stan")
          if (predictvar_type == "binary"){
            stanmodel <- stanmodels$model_predict_noaux
          } else {
            stanmodel <- stanmodels$model_predict_noaux_linear
          }
          fnew <- as.character(as.formula(formula))
          fnew[2] <- predictvar
          fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
          outcome <- c(as.matrix(model.frame(fnew, data0)[, 1]))

          datlist <- list(
            N = N, # obs
            J = J, # number of control items
            Y = Y, # number of affirmative answers
            K = K, # number of covariates
            X = X, # covariate matrix
            treat = treat, # treatment indicator
            outcome = outcome # outcome to be predicted
          )
        } else {
          # model_misreport_noaux = stan_model("list_misreport_constrained1.0.stan")
          stanmodel <- stanmodels$model_misreport_noaux
          fnew <- as.character(as.formula(formula))
          fnew[2] <- direct_item_misreport
          fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
          direct <- c(as.matrix(model.frame(fnew, data0)[, 1]))

          datlist <- list(
            N = N, # obs
            J = J, # number of control items
            Y = Y, # number of affirmative answers
            K = K, # number of covariates
            X = X, # covariate matrix
            treat = treat, # treatment indicator
            direct = direct # direct item
          )
        }

        if (type == "outcome") {
          pars = c("psi0","delta","rho0")
        } else if (type == "predict") {
          pars = c("psi0","delta","rho0","psi2","gamma0","phi")
        } else {
          pars = c("psi0", "delta", "rho0","gamma0","treat_e","U_e","Z_e")
        }

        if (vb == TRUE & only_vb == FALSE){
          stanvb <- try(vb(
            stanmodel,
            data = datlist,
            pars = pars,
            seed = seeds
          ),silent = T)

          if (inherits(stanvb, "try-error")){
            warning("Variational inference failed! Proceed to direct sampling.")
          } else{
            initvalues = summary(stanvb)$summary[,1]
          }

          if (type == "outcome") {
            init_fun = function(){
              list(psi0 = initvalues[1:K],
                   delta = initvalues[(K+1):(2*K)],
                   rho0 = initvalues[2*K+1])
            }
          } else if (type == "predict") {
            init_fun = function(){
              list(psi0 = initvalues[1:K],
                   delta = initvalues[(K+1):(2*K)],
                   rho0 = initvalues[2*K+1],
                   psi2 = initvalues[(2*K+2):(3*K+1)],
                   gamma0 = initvalues[(3*K+2)],
                   phi = initvalues[(3*K+3)]
              )
            }
          } else {
            init_fun = function(){
              list(psi0 = initvalues[1:K],
                   delta = initvalues[(K+1):(2*K)],
                   rho0 = initvalues[2*K+1],
                   gamma0 = initvalues[(2*K+2):(3*K+1)],
                   treat_e = initvalues[3*K+2],
                   U_e = initvalues[3*K+3],
                   Z_e = initvalues[3*K+4]
              )
            }
          }
          if (inherits(stanvb, "try-error")){
            stanout <- sampling(
              stanmodel,
              data = datlist,
              pars = pars,
              seed = seeds,
              iter = nsim,
              thin = thin,
              warmup = burnin,
              chains = nchain
            )
          } else {
            stanout <- sampling(
              stanmodel,
              data = datlist,
              pars = pars,
              init = init_fun,
              seed = seeds,
              iter = nsim,
              thin = thin,
              warmup = burnin,
              chains = nchain
            )
          }
        } else if (only_vb == TRUE) {
          stanvb <- try(vb(
            stanmodel,
            data = datlist,
            pars = pars,
            seed = seeds
          ),silent = T)

          if (inherits(stanvb, "try-error")){
            stop("Variational inference failed!")
          }
          stanout <- stanvb
        } else {
          stanout <- sampling(
            stanmodel,
            data = datlist,
            pars = pars,
            seed = seeds,
            iter = nsim,
            thin = thin,
            warmup = burnin,
            chains = nchain
          )
        }
        prior_est_mu = rstan::summary(stanout)$summary[,1]
        prior_est_sigma = rstan::summary(stanout)$summary[,3]
        if (type == "outcome") {
          mu_delta0 = prior_est_mu[(K+1):(2*K)]
          sigma_delta0 = prior_est_sigma[(K+1):(2*K)]
        } else if (type == "predict") {
          mu_delta0 = prior_est_mu[(K+1):(2*K)]
          sigma_delta0 = prior_est_sigma[(K+1):(2*K)]
        } else {
          mu_delta0 = prior_est_mu[(K+1):(2*K)]
          sigma_delta0 = prior_est_sigma[(K+1):(2*K)]
        }

        if (is.null(conjugate_k) == FALSE){
          #############################
          # conjugate distance scaled by df.
          difsq = (mu_delta - mu_delta0)^2
          varprior = NA
          for (i in 1:K){
            varprior[i] = max((sigma_delta^2)[i], difsq[i])
          }
          sigma_delta = sqrt(1/log(conjugate_k) * varprior)

          if (type == "outcome") {
            # model_outcome_infonorm = stan_model("list_outcome_constrained1.0_infonorm.stan")
            stanmodel <- stanmodels$model_outcome_infonorm
            datlist <- list(
              N = N, # obs
              J = J, # number of control items
              Y = Y, # number of affirmative answers
              K = K, # number of covariates
              X = X, # covariate matrix
              treat = treat, # treatment indicator
              mu_delta = as.array(mu_delta),
              sigma_delta = as.array(sigma_delta)
            )
          } else if (type == "predict") {
            # model_predict_infonorm = stan_model("list_predictor_constrained_logit1.0_infonorm.stan")
            if (predictvar_type == "binary"){
              stanmodel <- stanmodels$model_predict_infonorm
            } else {
              stanmodel <- stanmodels$model_predict_infonorm_linear
            }
            fnew <- as.character(as.formula(formula))
            fnew[2] <- predictvar
            fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
            outcome <- c(as.matrix(model.frame(fnew, data0)[, 1]))

            datlist <- list(
              N = N, # obs
              J = J, # number of control items
              Y = Y, # number of affirmative answers
              K = K, # number of covariates
              X = X, # covariate matrix
              treat = treat, # treatment indicator
              outcome = outcome, # outcome to be predicted
              mu_delta = as.array(mu_delta),
              sigma_delta = as.array(sigma_delta)
            )
          } else {
            # model_misreport_infonorm = stan_model("list_misreport_constrained1.0_infonorm.stan")
            stanmodel <- stanmodels$model_misreport_infonorm
            fnew <- as.character(as.formula(formula))
            fnew[2] <- direct_item_misreport
            fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
            direct <- c(as.matrix(model.frame(fnew, data0)[, 1]))

            datlist <- list(
              N = N, # obs
              J = J, # number of control items
              Y = Y, # number of affirmative answers
              K = K, # number of covariates
              X = X, # covariate matrix
              treat = treat, # treatment indicator
              direct = direct, # direct item
              mu_delta = as.array(mu_delta),
              sigma_delta = as.array(sigma_delta),
              mu_gamma = as.array(mu_delta00),
              sigma_gamma = as.array(sigma_delta00)
            )
          }

          if (type == "outcome") {
            pars = c("psi0","delta","rho0")
          } else if (type == "predict") {
            pars = c("psi0","delta","rho0","psi2","gamma0","phi")
          } else {
            pars = c("psi0", "delta", "rho0","gamma0","treat_e","U_e","Z_e")
          }

          if (vb == TRUE & only_vb == FALSE){
            stanvb <- try(vb(
              stanmodel,
              data = datlist,
              pars = pars,
              seed = seeds
            ),silent = T)

            if (inherits(stanvb, "try-error")){
              warning("Variational inference fails! Proceed to direct sampling.")
            } else{
              initvalues = summary(stanvb)$summary[,1]
            }

            if (type == "outcome") {
              init_fun = function(){
                list(psi0 = initvalues[1:K],
                     delta = initvalues[(K+1):(2*K)],
                     rho0 = initvalues[2*K+1])
              }
            } else if (type == "predict") {
              init_fun = function(){
                list(psi0 = initvalues[1:K],
                     delta = initvalues[(K+1):(2*K)],
                     rho0 = initvalues[2*K+1],
                     psi2 = initvalues[(2*K+2):(3*K+1)],
                     gamma0 = initvalues[(3*K+2)],
                     phi = initvalues[(3*K+3)]
                )
              }
            } else {
              init_fun = function(){
                list(psi0 = initvalues[1:K],
                     delta = initvalues[(K+1):(2*K)],
                     rho0 = initvalues[2*K+1],
                     gamma0 = initvalues[(2*K+2):(3*K+1)],
                     treat_e = initvalues[3*K+2],
                     U_e = initvalues[3*K+3],
                     Z_e = initvalues[3*K+4]
                )
              }
            }

            if (inherits(stanvb, "try-error")){
              stanout <- sampling(
                stanmodel,
                data = datlist,
                pars = pars,
                seed = seeds,
                iter = nsim,
                thin = thin,
                warmup = burnin,
                chains = nchain
              )
            } else {
              stanout <- sampling(
                stanmodel,
                data = datlist,
                pars = pars,
                init = init_fun,
                seed = seeds,
                iter = nsim,
                thin = thin,
                warmup = burnin,
                chains = nchain
              )
            }
          } else if (only_vb == TRUE) {
            stanvb <- try(vb(
              stanmodel,
              data = datlist,
              pars = pars,
              seed = seeds
            ),silent = T)

            if (inherits(stanvb, "try-error")){
              stop("Variational inference failed!")
            }
            stanout <- stanvb
          } else {
            stanout <- sampling(
              stanmodel,
              data = datlist,
              pars = pars,
              seed = seeds,
              iter = nsim,
              thin = thin,
              warmup = burnin,
              chains = nchain
            )
          }

        } else {
          ###########################################
          # conjugate distance without scaling by df
          difsq = (mu_delta - mu_delta0)^2
          varprior = NA
          for (i in 1:K){
            varprior[i] = max((sigma_delta^2)[i], difsq[i])
          }
          sigma_delta = sqrt(varprior)

          if (type == "outcome") {
            # model_outcome_infonorm = stan_model("list_outcome_constrained1.0_infonorm.stan")
            stanmodel <- stanmodels$model_outcome_infonorm
            datlist <- list(
              N = N, # obs
              J = J, # number of control items
              Y = Y, # number of affirmative answers
              K = K, # number of covariates
              X = X, # covariate matrix
              treat = treat, # treatment indicator
              mu_delta = as.array(mu_delta),
              sigma_delta = as.array(sigma_delta)
            )
          } else if (type == "predict") {
            # model_predict_infonorm = stan_model("list_predictor_constrained_logit1.0_infonorm.stan")
            if (predictvar_type == "binary"){
              stanmodel <- stanmodels$model_predict_infonorm
            } else {
              stanmodel <- stanmodels$model_predict_infonorm_linear
            }
            fnew <- as.character(as.formula(formula))
            fnew[2] <- predictvar
            fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
            outcome <- c(as.matrix(model.frame(fnew, data0)[, 1]))

            datlist <- list(
              N = N, # obs
              J = J, # number of control items
              Y = Y, # number of affirmative answers
              K = K, # number of covariates
              X = X, # covariate matrix
              treat = treat, # treatment indicator
              outcome = outcome, # outcome to be predicted
              mu_delta = as.array(mu_delta),
              sigma_delta = as.array(sigma_delta)
            )
          } else {
            # model_misreport_infonorm = stan_model("list_misreport_constrained1.0_infonorm.stan")
            stanmodel <- stanmodels$model_misreport_infonorm
            fnew <- as.character(as.formula(formula))
            fnew[2] <- direct_item_misreport
            fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
            direct <- c(as.matrix(model.frame(fnew, data0)[, 1]))

            datlist <- list(
              N = N, # obs
              J = J, # number of control items
              Y = Y, # number of affirmative answers
              K = K, # number of covariates
              X = X, # covariate matrix
              treat = treat, # treatment indicator
              direct = direct, # direct item
              mu_delta = as.array(mu_delta),
              sigma_delta = as.array(sigma_delta),
              mu_gamma = as.array(mu_delta00),
              sigma_gamma = as.array(sigma_delta00)
            )
          }

          if (type == "outcome") {
            pars = c("psi0","delta","rho0")
          } else if (type == "predict") {
            pars = c("psi0","delta","rho0","psi2","gamma0","phi")
          } else {
            pars = c("psi0", "delta", "rho0","gamma0","treat_e","U_e","Z_e")
          }

          if (vb == TRUE & only_vb == FALSE){
            stanvb <- try(vb(
              stanmodel,
              data = datlist,
              pars = pars,
              seed = seeds
            ),silent = T)

            if (inherits(stanvb, "try-error")){
              warning("Variational inference fails! Proceed to direct sampling.")
            } else{
              initvalues = summary(stanvb)$summary[,1]
            }

            if (type == "outcome") {
              init_fun = function(){
                list(psi0 = initvalues[1:K],
                     delta = initvalues[(K+1):(2*K)],
                     rho0 = initvalues[2*K+1])
              }
            } else if (type == "predict") {
              init_fun = function(){
                list(psi0 = initvalues[1:K],
                     delta = initvalues[(K+1):(2*K)],
                     rho0 = initvalues[2*K+1],
                     psi2 = initvalues[(2*K+2):(3*K+1)],
                     gamma0 = initvalues[(3*K+2)],
                     phi = initvalues[(3*K+3)]
                )
              }
            } else {
              init_fun = function(){
                list(psi0 = initvalues[1:K],
                     delta = initvalues[(K+1):(2*K)],
                     rho0 = initvalues[2*K+1],
                     gamma0 = initvalues[(2*K+2):(3*K+1)],
                     treat_e = initvalues[3*K+2],
                     U_e = initvalues[3*K+3],
                     Z_e = initvalues[3*K+4]
                )
              }
            }

            if (inherits(stanvb, "try-error")){
              stanout <- sampling(
                stanmodel,
                data = datlist,
                pars = pars,
                seed = seeds,
                iter = nsim,
                thin = thin,
                warmup = burnin,
                chains = nchain
              )
            } else {
              stanout <- sampling(
                stanmodel,
                data = datlist,
                pars = pars,
                init = init_fun,
                seed = seeds,
                iter = nsim,
                thin = thin,
                warmup = burnin,
                chains = nchain
              )
            }
          } else if (only_vb == TRUE) {
            stanvb <- try(vb(
              stanmodel,
              data = datlist,
              pars = pars,
              seed = seeds
            ),silent = T)

            if (inherits(stanvb, "try-error")){
              stop("Variational inference failed!")
            }
            stanout <- stanvb
          } else {
            stanout <- sampling(
              stanmodel,
              data = datlist,
              pars = pars,
              seed = seeds,
              iter = nsim,
              thin = thin,
              warmup = burnin,
              chains = nchain
            )
          }


        }

      } else {
        ###############################################
        # direct item without conjugate distance
        if (type == "outcome") {
          # model_outcome_infonorm = stan_model("list_outcome_constrained1.0_infonorm.stan")
          stanmodel <- stanmodels$model_outcome_infonorm
          datlist <- list(
            N = N, # obs
            J = J, # number of control items
            Y = Y, # number of affirmative answers
            K = K, # number of covariates
            X = X, # covariate matrix
            treat = treat, # treatment indicator
            mu_delta = as.array(mu_delta),
            sigma_delta = as.array(sigma_delta)
          )
        } else if (type == "predict") {
          # model_predict_infonorm = stan_model("list_predictor_constrained_logit1.0_infonorm.stan")
          if (predictvar_type == "binary"){
            stanmodel <- stanmodels$model_predict_infonorm
          } else {
            stanmodel <- stanmodels$model_predict_infonorm_linear
          }
          fnew <- as.character(as.formula(formula))
          fnew[2] <- predictvar
          fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
          outcome <- c(as.matrix(model.frame(fnew, data0)[, 1]))

          datlist <- list(
            N = N, # obs
            J = J, # number of control items
            Y = Y, # number of affirmative answers
            K = K, # number of covariates
            X = X, # covariate matrix
            treat = treat, # treatment indicator
            outcome = outcome, # outcome to be predicted
            mu_delta = as.array(mu_delta),
            sigma_delta = as.array(sigma_delta)
          )
        } else {
          # model_misreport_infonorm = stan_model("list_misreport_constrained1.0_infonorm.stan")
          stanmodel <- stanmodels$model_misreport_infonorm
          fnew <- as.character(as.formula(formula))
          fnew[2] <- direct_item_misreport
          fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
          direct <- c(as.matrix(model.frame(fnew, data0)[, 1]))

          datlist <- list(
            N = N, # obs
            J = J, # number of control items
            Y = Y, # number of affirmative answers
            K = K, # number of covariates
            X = X, # covariate matrix
            treat = treat, # treatment indicator
            direct = direct, # direct item
            mu_delta = as.array(mu_delta),
            sigma_delta = as.array(sigma_delta),
            mu_gamma = as.array(mu_delta00),
            sigma_gamma = as.array(sigma_delta00)
          )
        }

        if (type == "outcome") {
          pars = c("psi0","delta","rho0")
        } else if (type == "predict") {
          pars = c("psi0","delta","rho0","psi2","gamma0","phi")
        } else {
          pars = c("psi0", "delta", "rho0","gamma0","treat_e","U_e","Z_e")
        }

        if (vb == TRUE & only_vb == FALSE){
          stanvb <- try(vb(
            stanmodel,
            data = datlist,
            pars = pars,
            seed = seeds
          ),silent = T)

          if (inherits(stanvb, "try-error")){
            warning("Variational inference fails! Proceed to direct sampling.")
          } else{
            initvalues = summary(stanvb)$summary[,1]
          }

          if (type == "outcome") {
            init_fun = function(){
              list(psi0 = initvalues[1:K],
                   delta = initvalues[(K+1):(2*K)],
                   rho0 = initvalues[2*K+1])
            }
          } else if (type == "predict") {
            init_fun = function(){
              list(psi0 = initvalues[1:K],
                   delta = initvalues[(K+1):(2*K)],
                   rho0 = initvalues[2*K+1],
                   psi2 = initvalues[(2*K+2):(3*K+1)],
                   gamma0 = initvalues[(3*K+2)],
                   phi = initvalues[(3*K+3)]
              )
            }
          } else {
            init_fun = function(){
              list(psi0 = initvalues[1:K],
                   delta = initvalues[(K+1):(2*K)],
                   rho0 = initvalues[2*K+1],
                   gamma0 = initvalues[(2*K+2):(3*K+1)],
                   treat_e = initvalues[3*K+2],
                   U_e = initvalues[3*K+3],
                   Z_e = initvalues[3*K+4]
              )
            }
          }

          if (inherits(stanvb, "try-error")){
            stanout <- sampling(
              stanmodel,
              data = datlist,
              pars = pars,
              seed = seeds,
              iter = nsim,
              thin = thin,
              warmup = burnin,
              chains = nchain
            )
          } else {
            stanout <- sampling(
              stanmodel,
              data = datlist,
              pars = pars,
              init = init_fun,
              seed = seeds,
              iter = nsim,
              thin = thin,
              warmup = burnin,
              chains = nchain
            )
          }
        } else if (only_vb == TRUE) {
          stanvb <- try(vb(
            stanmodel,
            data = datlist,
            pars = pars,
            seed = seeds
          ),silent = T)

          if (inherits(stanvb, "try-error")){
            stop("Variational inference failed!")
          }
          stanout <- stanvb
        } else {
          stanout <- sampling(
            stanmodel,
            data = datlist,
            pars = pars,
            seed = seeds,
            iter = nsim,
            thin = thin,
            warmup = burnin,
            chains = nchain
          )
        }
      }


      ##############################################
      # III. Double list
    } else if (prior == "double_list"){

      if (type == "outcome") {
        # model_outcome_noaux = stan_model("list_outcome_constrained1.0.stan")
        stanmodel <- stanmodels$model_outcome_noaux
        datlist <- list(
          N = N, # obs
          J = J, # number of control items
          Y = Y, # number of affirmative answers
          K = K, # number of covariates
          X = X, # covariate matrix
          treat = treat # treatment indicator
        )
      } else if (type == "predict") {
        # model_predict_noaux = stan_model("list_predictor_constrained_logit1.0.stan")
        if (predictvar_type == "binary"){
          stanmodel <- stanmodels$model_predict_noaux
        } else {
          stanmodel <- stanmodels$model_predict_noaux_linear
        }
        fnew <- as.character(as.formula(formula))
        fnew[2] <- predictvar
        fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
        outcome <- c(as.matrix(model.frame(fnew, data0)[, 1]))

        datlist <- list(
          N = N, # obs
          J = J, # number of control items
          Y = Y, # number of affirmative answers
          K = K, # number of covariates
          X = X, # covariate matrix
          treat = treat, # treatment indicator
          outcome = outcome # outcome to be predicted
        )
      } else {
        # model_misreport_noaux = stan_model("list_misreport_constrained1.0.stan")
        stanmodel <- stanmodels$model_misreport_noaux
        fnew <- as.character(as.formula(formula))
        fnew[2] <- direct_item_misreport
        fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
        direct <- c(as.matrix(model.frame(fnew, data0)[, 1]))

        datlist <- list(
          N = N, # obs
          J = J, # number of control items
          Y = Y, # number of affirmative answers
          K = K, # number of covariates
          X = X, # covariate matrix
          treat = treat, # treatment indicator
          direct = direct # direct item
        )
      }

      if (type == "outcome") {
        pars = c("psi0","delta","rho0")
      } else if (type == "predict") {
        pars = c("psi0","delta","rho0","psi2","gamma0","phi")
      } else {
        pars = c("psi0", "delta", "rho0","gamma0","treat_e","U_e","Z_e")
      }

      if (vb == TRUE & only_vb == FALSE){
        stanvb <- try(vb(
          stanmodel,
          data = datlist,
          pars = pars,
          seed = seeds
        ),silent = T)

        if (inherits(stanvb, "try-error")){
          warning("Variational inference fails! Proceed to direct sampling.")
        } else{
          initvalues = summary(stanvb)$summary[,1]
        }

        if (type == "outcome") {
          init_fun = function(){
            list(psi0 = initvalues[1:K],
                 delta = initvalues[(K+1):(2*K)],
                 rho0 = initvalues[2*K+1])
          }
        } else if (type == "predict") {
          init_fun = function(){
            list(psi0 = initvalues[1:K],
                 delta = initvalues[(K+1):(2*K)],
                 rho0 = initvalues[2*K+1],
                 psi2 = initvalues[(2*K+2):(3*K+1)],
                 gamma0 = initvalues[(3*K+2)],
                 phi = initvalues[(3*K+3)]
            )
          }
        } else {
          init_fun = function(){
            list(psi0 = initvalues[1:K],
                 delta = initvalues[(K+1):(2*K)],
                 rho0 = initvalues[2*K+1],
                 gamma0 = initvalues[(2*K+2):(3*K+1)],
                 treat_e = initvalues[3*K+2],
                 U_e = initvalues[3*K+3],
                 Z_e = initvalues[3*K+4]
            )
          }
        }

        if (inherits(stanvb, "try-error")){
          stanout <- sampling(
            stanmodel,
            data = datlist,
            pars = pars,
            seed = seeds,
            iter = nsim,
            thin = thin,
            warmup = burnin,
            chains = nchain
          )
        } else {
          stanout <- sampling(
            stanmodel,
            data = datlist,
            pars = pars,
            init = init_fun,
            seed = seeds,
            iter = nsim,
            thin = thin,
            warmup = burnin,
            chains = nchain
          )
        }
      } else if (only_vb == TRUE) {
        stanvb <- try(vb(
          stanmodel,
          data = datlist,
          pars = pars,
          seed = seeds
        ),silent = T)

        if (inherits(stanvb, "try-error")){
          stop("Variational inference failed!")
        }
        stanout <- stanvb
      } else {
        stanout <- sampling(
          stanmodel,
          data = datlist,
          pars = pars,
          seed = seeds,
          iter = nsim,
          thin = thin,
          warmup = burnin,
          chains = nchain
        )
      }
      ######################################
      #### double list II.
      ######################################
      # priors from the first list experiment
      prior_est_mu = rstan::summary(stanout)$summary[,1]
      prior_est_sigma = rstan::summary(stanout)$summary[,3]

      fnew <- as.character(as.formula(formula))
      fnew[2] <- double_list
      fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))

      double_list <- c(as.matrix(model.frame(fnew, data0)[, 1]))
      X_doublelist <- model.matrix(fnew, data0)
      fnew <- as.character(fnew)
      fnew[3] <- paste(treat0,"+",fnew[3])
      treat_doublelist <- model.frame(as.formula(paste(fnew[2],fnew[1],fnew[3])),data0)[,2]
      N_doublelist <- length(double_list)

      if (type == "outcome") {
        mu_psi0 = prior_est_mu[1:K]
        sigma_psi0 = prior_est_sigma[1:K]
        mu_delta = prior_est_mu[(K+1):(2*K)]
        sigma_delta = prior_est_sigma[(K+1):(2*K)]
        mu_rho0 = prior_est_mu[2*K+1]
        sigma_rho0 = prior_est_sigma[2*K+1]

      } else if (type == "predict") {
        mu_psi0 = prior_est_mu[1:K]
        sigma_psi0 = prior_est_sigma[1:K]
        mu_delta = prior_est_mu[(K+1):(2*K)]
        sigma_delta = prior_est_sigma[(K+1):(2*K)]
        mu_rho0 = prior_est_mu[2*K+1]
        sigma_rho0 = prior_est_sigma[2*K+1]
        mu_psi2 = prior_est_mu[(2*K+2):(3*K+1)]
        sigma_psi2 = prior_est_sigma[(2*K+2):(3*K+1)]
        mu_gamma0 = prior_est_mu[(3*K+2)]
        sigma_gamma0 = prior_est_sigma[(3*K+2)]
        mu_phi = prior_est_mu[(3*K+3)]
        sigma_phi = prior_est_sigma[(3*K+3)]

      } else {
        mu_psi0 = prior_est_mu[1:K]
        sigma_psi0 = prior_est_sigma[1:K]
        mu_delta = prior_est_mu[(K+1):(2*K)]
        sigma_delta = prior_est_sigma[(K+1):(2*K)]
        mu_rho0 = prior_est_mu[2*K+1]
        sigma_rho0 = prior_est_sigma[2*K+1]
        mu_gamma0 = prior_est_mu[(2*K+2):(3*K+1)]
        sigma_gamma0 = prior_est_sigma[(2*K+2):(3*K+1)]
        mu_treate = prior_est_mu[3*K+2]
        sigma_treate = prior_est_sigma[3*K+2]
        mu_ue = prior_est_mu[3*K+3]
        sigma_ue = prior_est_sigma[3*K+3]
        mu_ze = prior_est_mu[3*K+4]
        sigma_ze = prior_est_sigma[3*K+4]
      }

      if (type == "outcome") {
        # model_outcome_doublelist = stan_model("list_outcome_constrained3.0_normal_prior.stan")
        stanmodel <- stanmodels$model_outcome_doublelist
        datlist <- list(
          N = N_doublelist, # obs
          J = J, # number of control items
          Y = double_list, # number of affirmative answers
          K = K, # number of covariates
          X = X_doublelist, # covariate matrix
          treat = treat_doublelist, # treatment indicator
          mu_psi0 = as.array(mu_psi0),
          sigma_psi0 = as.array(sigma_psi0),
          mu_delta = as.array(mu_delta),
          sigma_delta = as.array(sigma_delta),
          mu_rho0 = mu_rho0,
          sigma_rho0 = sigma_rho0
        )
      } else if (type == "predict") {
        # model_predict_doublelist = stan_model("list_predictor_constrained_logit1.0_doublelist.stan")
        if (predictvar_type == "binary"){
          stanmodel <- stanmodels$model_predict_doublelist
        } else {
          stanmodel <- stanmodels$model_predict_doublelist_linear
        }
        fnew <- as.character(as.formula(formula))
        fnew[2] <- predictvar
        fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
        outcome <- c(as.matrix(model.frame(fnew, data0)[, 1]))

        datlist <- list(
          N = N_doublelist, # obs
          J = J, # number of control items
          Y = double_list, # number of affirmative answers
          K = K, # number of covariates
          X = X_doublelist, # covariate matrix
          treat = treat_doublelist, # treatment indicator
          outcome = outcome, # outcome to be predicted
          mu_psi0 = as.array(mu_psi0),
          sigma_psi0 = as.array(sigma_psi0),
          mu_rho0 = mu_rho0,
          sigma_rho0 = sigma_rho0,
          mu_delta = as.array(mu_delta),
          sigma_delta = as.array(sigma_delta),
          mu_psi2 = as.array(mu_psi2),
          sigma_psi2 = as.array(sigma_psi2),
          mu_phi = mu_phi,
          sigma_phi = sigma_phi,
          mu_gamma0 = mu_gamma0,
          sigma_gamma0 = sigma_gamma0
        )
      } else {
        # model_misreport_doublelist = stan_model("list_misreport_constrained_normal_prior2.0.stan")
        stanmodel <- stanmodels$model_misreport_doublelist
        fnew <- as.character(as.formula(formula))
        fnew[2] <- direct_item_misreport
        fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
        direct <- c(as.matrix(model.frame(fnew, data0)[, 1]))

        datlist <- list(
          N = N_doublelist, # obs
          J = J, # number of control items
          Y = double_list, # number of affirmative answers
          K = K, # number of covariates
          X = X_doublelist, # covariate matrix
          treat = treat_doublelist, # treatment indicator
          direct = direct, # direct item
          mu_psi0 = as.array(mu_psi0),
          sigma_psi0 = as.array(sigma_psi0),
          mu_rho0 = mu_rho0,
          sigma_rho0 = sigma_rho0,
          mu_delta = as.array(mu_delta),
          sigma_delta = as.array(sigma_delta),
          mu_gamma0 = as.array(mu_gamma0),
          sigma_gamma0 = as.array(sigma_gamma0),
          mu_treate = mu_treate,
          sigma_treate = sigma_treate,
          mu_ue = mu_ue,
          sigma_ue = sigma_ue,
          mu_ze = mu_ze,
          sigma_ze = sigma_ze
        )
      }

      if (type == "outcome") {
        pars = c("psi0","delta","rho0")
      } else if (type == "predict") {
        pars = c("psi0","delta","rho0","psi2","gamma0","phi")
      } else {
        pars = c("psi0", "delta", "rho0","gamma0","treat_e","U_e","Z_e")
      }

      if (vb == TRUE & only_vb == FALSE){
        stanvb <- try(vb(
          stanmodel,
          data = datlist,
          pars = pars,
          seed = seeds
        ),silent = T)

        if (inherits(stanvb, "try-error")){
          warning("Variational inference fails! Proceed to direct sampling.")
        } else{
          initvalues = summary(stanvb)$summary[,1]
        }

        if (type == "outcome") {
          init_fun = function(){
            list(psi0 = initvalues[1:K],
                 delta = initvalues[(K+1):(2*K)],
                 rho0 = initvalues[2*K+1])
          }
        } else if (type == "predict") {
          init_fun = function(){
            list(psi0 = initvalues[1:K],
                 delta = initvalues[(K+1):(2*K)],
                 rho0 = initvalues[2*K+1],
                 psi2 = initvalues[(2*K+2):(3*K+1)],
                 gamma0 = initvalues[(3*K+2)],
                 phi = initvalues[(3*K+3)]
            )
          }
        } else {
          init_fun = function(){
            list(psi0 = initvalues[1:K],
                 delta = initvalues[(K+1):(2*K)],
                 rho0 = initvalues[2*K+1],
                 gamma0 = initvalues[(2*K+2):(3*K+1)],
                 treat_e = initvalues[3*K+2],
                 U_e = initvalues[3*K+3],
                 Z_e = initvalues[3*K+4]
            )
          }
        }

        if (inherits(stanvb, "try-error")){
          stanout <- sampling(
            stanmodel,
            data = datlist,
            pars = pars,
            seed = seeds,
            iter = nsim,
            thin = thin,
            warmup = burnin,
            chains = nchain
          )
        } else {
          stanout <- sampling(
            stanmodel,
            data = datlist,
            pars = pars,
            init = init_fun,
            seed = seeds,
            iter = nsim,
            thin = thin,
            warmup = burnin,
            chains = nchain
          )
        }
      } else if (only_vb == TRUE) {
        stanvb <- try(vb(
          stanmodel,
          data = datlist,
          pars = pars,
          seed = seeds
        ),silent = T)

        if (inherits(stanvb, "try-error")){
          stop("Variational inference failed!")
        }
        stanout <- stanvb
      } else {
        stanout <- sampling(
          stanmodel,
          data = datlist,
          pars = pars,
          seed = seeds,
          iter = nsim,
          thin = thin,
          warmup = burnin,
          chains = nchain
        )
      }
      ##############################################
      # Otherwise, error
    } else {
      stop("Please choose prior types among 'auxiliary', 'direct_item', 'double_list', and 'BL'!")
    }
  } else {
    ###############################################
    # Without informative priors
    ###############################################
    if (robust == TRUE){
      if (K != 1){
        warning("Robust estimation is only applicable to intercept-only models. The model will be estimated without robust constraints.")
      } else {
        # Use difference-in-means estimate to constrain the parameter space
        DiM_est = mean(Y[which(treat == 1)]) - mean(Y[which(treat == 0)])
        if (DiM_est > 1) {DiM_est = 1}
        if (DiM_est < 0) {DiM_est = 0}
        # hyper parameter for BL prior
        prior_a = DiM_est * N
        prior_b = (1 - DiM_est) * N
        if (prior_a == 0) {prior_a = prior_a + 1}
        if (prior_b == 0) {prior_b = prior_b + 1}

        if (type == "outcome") {
          # model_outcome_noaux = stan_model("list_outcome_constrained1.0.stan")
          stanmodel <- stanmodels$model_outcome_cmp
          datlist <- list(
            N = N, # obs
            J = J, # number of control items
            Y = Y, # number of affirmative answers
            K = K, # number of covariates
            X = X, # covariate matrix
            treat = treat, # treatment indicator
            a = as.array(prior_a),
            b = as.array(prior_b)
          )
        } else if (type == "predict") {
          # model_predict_noaux = stan_model("list_predictor_constrained_logit1.0.stan")
          if (predictvar_type == "binary"){
            stanmodel <- stanmodels$model_predict_cmp
          } else {
            stanmodel <- stanmodels$model_predict_cmp_linear
          }
          fnew <- as.character(as.formula(formula))
          fnew[2] <- predictvar
          fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
          outcome <- c(as.matrix(model.frame(fnew, data0)[, 1]))

          datlist <- list(
            N = N, # obs
            J = J, # number of control items
            Y = Y, # number of affirmative answers
            K = K, # number of covariates
            X = X, # covariate matrix
            treat = treat, # treatment indicator
            outcome = outcome, # outcome to be predicted
            a = as.array(prior_a),
            b = as.array(prior_b)
          )
        } else {
          # model_misreport_noaux = stan_model("list_misreport_constrained1.0.stan")
          stanmodel <- stanmodels$model_misreport_cmp
          fnew <- as.character(as.formula(formula))
          fnew[2] <- direct_item_misreport
          fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
          direct <- c(as.matrix(model.frame(fnew, data0)[, 1]))

          datlist <- list(
            N = N, # obs
            J = J, # number of control items
            Y = Y, # number of affirmative answers
            K = K, # number of covariates
            X = X, # covariate matrix
            treat = treat, # treatment indicator
            direct = direct, # direct item
            a = as.array(prior_a),
            b = as.array(prior_b)
          )
        }
      }
    } else {
      if (type == "outcome") {
        # model_outcome_noaux = stan_model("list_outcome_constrained1.0.stan")
        stanmodel <- stanmodels$model_outcome_noaux
        datlist <- list(
          N = N, # obs
          J = J, # number of control items
          Y = Y, # number of affirmative answers
          K = K, # number of covariates
          X = X, # covariate matrix
          treat = treat # treatment indicator
        )
      } else if (type == "predict") {
        # model_predict_noaux = stan_model("list_predictor_constrained_logit1.0.stan")
        if (predictvar_type == "binary"){
          stanmodel <- stanmodels$model_predict_noaux
        } else {
          stanmodel <- stanmodels$model_predict_noaux_linear
        }
        fnew <- as.character(as.formula(formula))
        fnew[2] <- predictvar
        fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
        outcome <- c(as.matrix(model.frame(fnew, data0)[, 1]))

        datlist <- list(
          N = N, # obs
          J = J, # number of control items
          Y = Y, # number of affirmative answers
          K = K, # number of covariates
          X = X, # covariate matrix
          treat = treat, # treatment indicator
          outcome = outcome # outcome to be predicted
        )
      } else {
        # model_misreport_noaux = stan_model("list_misreport_constrained1.0.stan")
        stanmodel <- stanmodels$model_misreport_noaux
        fnew <- as.character(as.formula(formula))
        fnew[2] <- direct_item_misreport
        fnew = as.formula(paste(fnew[2],fnew[1],fnew[3]))
        direct <- c(as.matrix(model.frame(fnew, data0)[, 1]))

        datlist <- list(
          N = N, # obs
          J = J, # number of control items
          Y = Y, # number of affirmative answers
          K = K, # number of covariates
          X = X, # covariate matrix
          treat = treat, # treatment indicator
          direct = direct # direct item
        )
      }
    }


    if (type == "outcome") {
      pars = c("psi0","delta","rho0")
    } else if (type == "predict") {
      pars = c("psi0","delta","rho0","psi2","gamma0","phi")
    } else {
      pars = c("psi0", "delta", "rho0","gamma0","treat_e","U_e","Z_e")
    }

    if (vb == TRUE & only_vb == FALSE){
      stanvb <- try(vb(
        stanmodel,
        data = datlist,
        pars = pars,
        seed = seeds
      ),silent = T)

      if (inherits(stanvb, "try-error")){
        warning("Variational inference fails! Proceed to direct sampling.")
      } else{
        initvalues = summary(stanvb)$summary[,1]
      }

      if (type == "outcome") {
        init_fun = function(){
          list(psi0 = initvalues[1:K],
               delta = initvalues[(K+1):(2*K)],
               rho0 = initvalues[2*K+1])
        }
      } else if (type == "predict") {
        init_fun = function(){
          list(psi0 = initvalues[1:K],
               delta = initvalues[(K+1):(2*K)],
               rho0 = initvalues[2*K+1],
               psi2 = initvalues[(2*K+2):(3*K+1)],
               gamma0 = initvalues[(3*K+2)],
               phi = initvalues[(3*K+3)]
          )
        }
      } else {
        init_fun = function(){
          list(psi0 = initvalues[1:K],
               delta = initvalues[(K+1):(2*K)],
               rho0 = initvalues[2*K+1],
               gamma0 = initvalues[(2*K+2):(3*K+1)],
               treat_e = initvalues[3*K+2],
               U_e = initvalues[3*K+3],
               Z_e = initvalues[3*K+4]
          )
        }
      }

      if (inherits(stanvb, "try-error")){
        stanout <- sampling(
          stanmodel,
          data = datlist,
          pars = pars,
          seed = seeds,
          iter = nsim,
          thin = thin,
          warmup = burnin,
          chains = nchain
        )
      } else {
        stanout <- sampling(
          stanmodel,
          data = datlist,
          pars = pars,
          init = init_fun,
          seed = seeds,
          iter = nsim,
          thin = thin,
          warmup = burnin,
          chains = nchain
        )
      }
    } else if (only_vb == TRUE) {
      stanvb <- try(vb(
        stanmodel,
        data = datlist,
        pars = pars,
        seed = seeds
      ),silent = T)

      if (inherits(stanvb, "try-error")){
        stop("Variational inference failed!")
      }
      stanout <- stanvb
    } else {
      stanout <- sampling(
        stanmodel,
        data = datlist,
        pars = pars,
        seed = seeds,
        iter = nsim,
        thin = thin,
        warmup = burnin,
        chains = nchain
      )
    }

  }



  summaryout <- rstan::summary(stanout)$summary
  sampledf <- as.data.frame(stanout)

  # output
  out <- list()
  class(out) <- c("bayeslist", class(out))
  out$Call <- match.call()
  out$formula <- formula
  out$J <- J
  out$type <- type
  out$nsim <- nsim
  out$burnin <- burnin
  out$thin <- thin
  out$seeds <- seeds
  out$CIsize  <- CIsize
  out$data   <- data
  out$X    <- X
  out$Y <- Y
  out$xnames <- colnames(X)
  out$stanfit <- stanout
  out$sampledf <- sampledf
  out$summaryout <- summaryout
  out$npars <- K
  out$only_vb <- only_vb
  out$prior = prior
  out$direct_item = direct_item
  out$double_list = double_list
  out$double_list_treat = double_list_treat
  out$aux_info = aux_info
  out$ulbs <-
    apply(sampledf, 2, quantile, probs = c((1 - CIsize) / 2, 1 - (1 - CIsize) / 2))
  out$means <- apply(sampledf, 2, mean)
  out$treat = treat
  out$predictvar = predictvar
  out$robust = robust
  # out$direct = direct

  return(out)

}
