#' Extract coefficients from a \code{bayeslist} object
#'
#' Create a table of coefficient results from a \code{bayeslist} object.
#'
#' @param object A \code{bayeslist} object from running the main function \code{\link{bayeslist}}.
#' @param ... Further arguments to be passed according to \code{\link[stats]{coef}}.
#'
#' @return A table of coefficients with their corresponding lower and upper bounds.
#' @export
#'
#' @method coef bayeslist
#'
#'
coef.bayeslist <- coefficients.bayeslist <- function(object, ...) {
  coefmat <- object$summaryout
  N <- dim(coefmat)[1]
  coefmat = coefmat[-N,]
  if (object$type == "outcome") {
    if (! is.null(object$prior)){
      if (object$prior == "auxiliary"){
        vnames = c(paste0("Control items:",object$xnames),paste0("Sensitive item:",object$xnames),"Rho","Sigma")
      } else {
        vnames = c(paste0("Control items:",object$xnames),paste0("Sensitive item:",object$xnames),"Rho")
      }
    } else {
      vnames = c(paste0("Control items:",object$xnames),paste0("Sensitive item:",object$xnames),"Rho")
    }
  } else if (object$type == "predict"){
    if (! is.null(object$prior)){
      if (object$prior == "auxiliary"){
        vnames = c(paste0("Control items:",object$xnames),paste0("Sensitive item:",object$xnames),"Rho",paste0("Predict:",object$xnames),"Gamma","Phi","Sigma")
      } else {
        vnames = c(paste0("Control items:",object$xnames),paste0("Sensitive item:",object$xnames),"Rho",paste0("Predict:",object$xnames),"Gamma","Phi")
      }
    } else {
      vnames = c(paste0("Control items:",object$xnames),paste0("Sensitive item:",object$xnames),"Rho",paste0("Predict:",object$xnames),"Gamma","Phi")
    }
  } else if (object$type == "misreport") {
    if (! is.null(object$prior)){
      if (object$prior == "auxiliary"){
        vnames = c(paste0("Control items:",object$xnames),paste0("Sensitive item:",object$xnames),"Rho",paste0("Misreport:",object$xnames),"Treatment effect on misreport","U_e","Z_e","Sigma")
      } else {
        vnames = c(paste0("Control items:",object$xnames),paste0("Sensitive item:",object$xnames),"Rho",paste0("Misreport:",object$xnames),"Treatment effect on misreport","U_e","Z_e")
      }
    } else {
      vnames = c(paste0("Control items:",object$xnames),paste0("Sensitive item:",object$xnames),"Rho",paste0("Misreport:",object$xnames),"Treatment effect on misreport","U_e","Z_e")
    }
  } else {
    stop("Unknown type of sensitive item models!")
  }
  rownames(coefmat) = vnames
  return(coefmat)
}


#' Print returns from a \code{bayeslist} object
#'
#' General print function for \code{bayeslist} objects, which dispatches the chosen type
#' of printing to the corresponding function.
#'
#' @param x A \code{bayeslist} object to be printed.
#' @param type Character string giving the type of printing, such as
#'   \code{"text"}, \code{"mcmc"}, \code{"coef"}.
#' @param ... Additional arguments to be passed to print functions (check the See Also section).
#'
#' @return None.
#' @seealso \code{\link{print_text.bayeslist}}, \code{\link{print_mcmc.bayeslist}}, \code{\link{print_coef.bayeslist}}.
#' @export
#'
#'
#'
print.bayeslist <- function(x, type = "text", ...) {
  printFunName <- paste0("print_", type, ".bayeslist")
  do.call(printFunName, args = c(list(object = x), list(...)))
}


#' Summary of a \code{bayeslist} object
#'
#' General summary function for \code{bayeslist} objects.
#'
#' @param object A \code{bayeslist} object to be summarized.
#' @param ... Additional arguments to be passed to summary function.
#'
#' @return None.
#' @seealso \code{\link{print_text.bayeslist}}, \code{\link{print_mcmc.bayeslist}}, \code{\link{print_coef.bayeslist}}.
#' @method summary bayeslist
#' @S3method summary bayeslist
#' @export
#'
#'
#'
summary.bayeslist <- function(object, ...) {
  printFunName <- "print_text.bayeslist"
  do.call(printFunName, args = c(list(object = object), list(...)))
}


#' Predicted prevalence from a \code{bayeslist} object
#'
#' Prediction function for \code{bayeslist} objects.
#'
#' @param object A \code{bayeslist} object to be summarized.
#' @param ... Additional arguments to be passed to summary function.
#'
#' @return None.
#' @seealso \code{\link{print_text.bayeslist}}, \code{\link{print_mcmc.bayeslist}}, \code{\link{print_coef.bayeslist}}.
#' @export
#'
#'
#'
predict.bayeslist <- function(object, ...) {
  tmp = object$X %*% t(object$sampledf[,grep("delta",names(object$sampledf))])
  tmp = apply(tmp,1,quantile,probs = c(0.5))
  return(logistic(c(tmp)))
}

#' Print the main results from a \code{bayeslist} object.
#'
#' @param object A \code{bayeslist} object.
#' @param digits Number of digits to display.
#'
#' @return None.
#' @export
#'
#'
#'
print_text.bayeslist <- function(object, digits = 3) {
  cat("Bayesian sensitive item model: ",
      object$type, "With number of control items =", object$J,
      "\n")
  cat("\nCall:\n",
      paste(deparse(object$Call), sep = "\n", collapse = "\n"),
      "\n\n",
      sep = "")
  cat(
    "MCMC run for",
    object$nsim,
    "iterations, with",
    object$stanfit@sim$warmup2[1],
    "used. \n\n"
  )
  cat("Coefficients:\n")
  print(round(coef(object), digits))
  cat("\n")
}


#' Print convergence diagnostics from a \code{bayeslist} object
#'
#' \code{print_mcmc.bayeslist} prints a number of diagnostics about the convergence of a \code{bayeslist} objects.
#'
#'
#' @param object A \code{bayeslist} object.
#' @param ... Additional arguments to be passed to the \code{print} function.
#'
#' @return None.
#' @export
#'
#'
#'
print_mcmc.bayeslist <- function(object, ...) {
  print(object$stanfit, ...)
}



#' Print coefficients of a \code{bayeslist} object
#'
#' \code{print_coef.bayeslist} prints out coefficients from a \code{bayeslist} object from running the main function \code{\link{bayeslist}}.
#'
#' @param object A \code{bayeslist} object.
#' @param digits Number of digits to display.
#'
#' @return None.
#' @export
#'
#'
#'
print_coef.bayeslist <- function(object, digits = 3) {
  print(round(coef(object), digits))
}
