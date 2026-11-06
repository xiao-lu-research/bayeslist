#' Extract posterior sample of the pointwise log-likelihood from a bayesICT object.
#'
#' Convenience function for extracting the pointwise log-likelihood
#' matrix or array from a fitted Stan model.
#'
#' @param object A bayesICT object of the bayesICT fitted model.
#' @param permuted A logical scalar indicating whether the draws after
#' the \code{warmup} period in each chain should be permuted and merged.
#' If FALSE, the original order is kept. For each \code{stanfit} object,
#' the permutation is fixed (i.e., extracting samples a second time
#' will give the same sequence of iterations).
#' @param ... additional values need in log_lik methods
#'
#' @return
#' Usually, an S x N matrix containing the pointwise log-likelihood
#' samples, where S is the number of samples and N is the number
#' of observations in the data. If \code{permuted} is \code{FALSE},
#' an S x N x R array is returned, where R is the number of fitted
#' chains.
#'
#' @aliases log_lik
#'
#' @references
#' Vehtari, A., Gelman, A., & Gabry J. (2016). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. \emph{In Statistics
#' and Computing}, \doi{10.1007/s11222-016-9696-4}.
#'
#' Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive
#' information criteria for Bayesian models. \emph{Statistics and Computing}.
#'  24, 997-1016.
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
#' and widely applicable information criterion in singular learning theory.
#' \emph{The Journal of Machine Learning Research}. 11, 3571-3594.
#'
#' @importFrom rstantools log_lik
#' @importFrom loo  extract_log_lik
#' @method log_lik bayesICT
#' @export
#' @export log_lik
#'
#' @examples
#'
#' \donttest{
#'  mod1 <- bayesICT(sexassault ~ 1, data = srilanka, treat = "treatment", J = 3, type = "outcome",
#'      nsim = 200, thin = 1, CIsize = 0.95, nchain = 2, seeds = 342321, prior = NULL)
#'
#'  log1 = log_lik(mod1)
#'  log1
#' }
#'
log_lik.bayesICT = function(object,permuted = TRUE,...){
  if(!is.bayesICT(object))
    stop("The current object is not of bayesICT class")

  log_lik = loo::extract_log_lik(object$stanfit,merge_chains = permuted)
  return(log_lik)
}
#' Extract posterior sample of the accumulated log-likelihood from a bayesICT object
#'
#' Convenience function for extracting the posterior sample of the accumulated
#' log-likelihood array from a fitted bayesICT object.
#'
#' @param object A bayesICT object of the time series fitted model.
#' @param permuted A logical scalar indicating whether the draws after
#' the warmup period in each chain should be permuted and merged.
#' If FALSE, the original order is kept. For each stanfit object,
#' the permutation is fixed (i.e., extracting samples a second time
#' will give the same sequence of iterations).
#'
#' @return
#' A real value with the accumulated log likelihood.
#'
#' @export
#'
#' @references
#' Vehtari, A., Gelman, A., & Gabry J. (2016). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. \emph{In Statistics
#' and Computing}, \doi{10.1007/s11222-016-9696-4}.
#'
#' Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive
#' information criteria for Bayesian models. \emph{Statistics and Computing}.
#'  24, 997-1016.
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
#' and widely applicable information criterion in singular learning theory.
#' \emph{The Journal of Machine Learning Research}. 11, 3571-3594.
#'
#' @examples
#' \donttest{
#'# Estimate sensitive item outcome model using Sri Lanka data on male sexual violence
#'# Load Sri Lanka list experiment data
#'data(srilanka)
#'
#'# Model 1: intercept-only outcome model without prior information:
#'  mod1 <- bayesICT(sexassault ~ 1, data = srilanka, treat = "treatment", J = 3, type = "outcome",
#'      nsim = 200, thin = 1, CIsize = 0.95, nchain = 2, seeds = 342321, prior = NULL)
#'
#'  log1 = loglik(mod1)
#'  log1
#' }
#'
loglik = function(object,permuted = TRUE){
  if(!is.bayesICT(object))
    stop("The current object is not of bayesICT class")

  loglik = data.frame(extract_stan(object = object,pars = "log_lik",permuted = permuted))

  if(permuted)
    loglik = as.numeric(loglik$loglik)
  else{
    colnames(loglik) = paste0("loglik.",1:ncol(loglik))
    loglik = as.matrix(loglik)
  }
  return(loglik)
}
#' Leave-one-out cross-validation
#'
#' The \code{loo} method for bayesICT objects. Computes approximate
#' leave-one-out cross-validation using Pareto smoothed importance
#' sampling (PSIS-LOO CV).
#'
#' @aliases loo
#'
#' @param x A bayesICT object
#' @param ... additional values need in loo methods
#'
#' @return an object from the loo class with the results of the Pareto-Smooth
#' Importance Sampling, leave one out cross validation for model selection.
#'
#' @seealso
#'  * The **loo** package [vignettes](https://mc-stan.org/loo/articles/index.html)
#'    for demonstrations.
#'  * [psis()] for the underlying Pareto Smoothed Importance Sampling (PSIS)
#'    procedure used in the LOO-CV approximation.
#'  * [pareto-k-diagnostic] for convenience functions for looking at diagnostics.
#'  * [loo_compare()] for model comparison.
#'
#' @references
#' Vehtari, A., Gelman, A., & Gabry J. (2016). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. \emph{In Statistics
#' and Computing}, \doi{10.1007/s11222-016-9696-4}.
#'
#' Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive
#' information criteria for Bayesian models. \emph{Statistics and Computing}.
#'  24, 997-1016.
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
#' and widely applicable information criterion in singular learning theory.
#' \emph{The Journal of Machine Learning Research}. 11, 3571-3594.
#'
#' @importFrom rstan loo
#' @method loo bayesICT
#' @export loo
#' @export
#'
#'
#' @examples
#' \donttest{
#'# Estimate sensitive item outcome model using Sri Lanka data on male sexual violence
#'# Load Sri Lanka list experiment data
#'data(srilanka)
#'
#'# Model 1: intercept-only outcome model without prior information:
#'  mod1 <- bayesICT(sexassault ~ 1, data = srilanka, treat = "treatment", J = 3, type = "outcome",
#'      nsim = 200, thin = 1, CIsize = 0.95, nchain = 2, seeds = 342321, prior = NULL)
#'
#'  loo1 = loo(mod1)
#'  loo1
#' }
#'
loo.bayesICT = function(x,...){
  if (!is.bayesICT(x))
    stop("The current object is not a bayesICT class")

  return(rstan::loo(x$stanfit) )
}
#' Widely Applicable Information Criterion (WAIC)
#'
#' Compute the widely applicable information criterion (WAIC)
#' based on the posterior likelihood using the \pkg{loo} package.
#' For more details see \code{\link[loo:waic]{waic}}.
#'
#' @param x A bayesICT object
#' @param ... additional values need in waic methods
#'
#' @aliases  waic
#'
#' @details See the \code{loo_compare} function of the \pkg{loo} package
#' for more details on model comparisons.
#'
#' @return An object of class \code{loo}. With the estimates of the
#' Watanabe-Akaike Information criteria.
#'
#' @references
#' Vehtari, A., Gelman, A., & Gabry J. (2016). Practical Bayesian model
#' evaluation using leave-one-out cross-validation and WAIC. \emph{In Statistics
#' and Computing}, \doi{10.1007/s11222-016-9696-4}.
#'
#' Gelman, A., Hwang, J., & Vehtari, A. (2014). Understanding predictive
#' information criteria for Bayesian models. \emph{Statistics and Computing}.
#'  24, 997-1016.
#'
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation
#' and widely applicable information criterion in singular learning theory.
#' \emph{The Journal of Machine Learning Research}. 11, 3571-3594.
#'
#' @importFrom loo waic
#' @method waic bayesICT
#' @export waic
#' @export
#'
#' @examples
#' \donttest{
#'# Estimate sensitive item outcome model using Sri Lanka data on male sexual violence
#'# Load Sri Lanka list experiment data
#'data(srilanka)
#'
#'# Model 1: intercept-only outcome model without prior information:
#'  mod1 <- bayesICT(sexassault ~ 1, data = srilanka, treat = "treatment", J = 3, type = "outcome",
#'      nsim = 200, thin = 1, CIsize = 0.95, nchain = 2, seeds = 342321, prior = NULL)
#'
#'  waic1 = waic(mod1)
#'  waic1
#' }
#'
waic.bayesICT = function(x,...){
  if (!is.bayesICT(x))
    stop("The current object is not a bayesICT class")

  return(loo::waic(log_lik.bayesICT(x)) )
}
#' Checks if is a bayesICT object
#'
#' @method is bayesICT
#'
#' @param object a bayesICT object
#'
#' @noRd
#'
is.bayesICT = function(object){
  y = FALSE
  if(is(object,"bayesICT") ) y = TRUE
  return (y)
}
#' Extract chains of an stanfit object implemented in rstan package
#'
#' @usage extract_stan(object,pars,permuted = TRUE,inc_warmup = FALSE,include = TRUE, ...)
#'
#' @param object a bayesICT object
#' @param pars n optional character vector providing the parameter names
#' (or other quantity names) of interest. If not specified, all parameters
#' and other quantities are used. The log-posterior with name lp__ is also
#' included by default.
#' @param permuted A logical scalar indicating whether the draws
#' after the warmup period in each chain should be permuted and merged.
#' If FALSE, the original order is kept. For each stanfit object,
#' the permutation is fixed (i.e., extracting samples a second time will
#' give the same sequence of iterations).
#' @param inc_warmup A logical scalar indicating whether to include the
#' warmup draws. This argument is only relevant if permuted is FALSE.
#' @param include A logical scalar indicating whether the parameters named
#' in pars should be included (TRUE) or excluded (FALSE).
#' @param ... Further arguments passed to  \code{extract} function.
#'
#' @author Asael Alonzo Matamoros
#'
#' @return a list with the posterior samples of the provided parameters.
#'
#' @importFrom rstan extract
#' @export
#'
#' @examples
#' \donttest{
#'  mod1 <- bayesICT(sexassault ~ 1, data = srilanka, treat = "treatment", J = 3, type = "outcome",
#'      nsim = 2000, thin = 1, CIsize = 0.95, nchain = 2, seeds = 342321, prior = NULL)
#'
#'  # Extracting the Proportion of formative on the sensitive item.
#'  PAS = extract_stan(mod1, pars = "pas")
#' }
#'
extract_stan = function(object,pars,permuted = TRUE,
                        inc_warmup = FALSE,include = TRUE,...){

  if(!is.bayesICT(object) )
    stop("The current object is not a bayesICT class")

  chains = rstan::extract(object$stanfit,pars,permuted,inc_warmup,include,...)

  return(chains)
}
