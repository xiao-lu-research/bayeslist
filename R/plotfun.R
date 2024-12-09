#' Plot bayeslist object
#'
#' General plot function for \code{bayeslist} objects, which dispatches the chosen
#' type of plotting to the corresponding function.
#'
#' @param x A \code{bayeslist} object to be plotted.
#' @param type Character string giving the type of plotting. The options are
#'   \code{"trace"} for trace plots, \code{"prevalence"} for prevalence plots. The default is "prevalence".
#' @param ... Additional arguments to be passed to subsequent plot functions (check the See Also section).
#'
#' @import ggplot2
#'
#' @return None.
#'
#' @seealso \code{\link{plot_trace.bayeslist}} and \code{\link{plot_coef.bayeslist}}.
#' @method plot bayeslist
#' @S3method plot bayeslist
#' @export
#'
#'
#'
plot.bayeslist <- function(x, type = "prevalence", ...) {
  printFunName <- paste0("plot_", type, ".bayeslist")
  do.call(printFunName, args = c(list(object = x), list(...)))
}

#' Logistic function
#'
#' Standard logistic function.
#'
#' @param x A scalar or vector to be logit-transformed.
#'
#' @return logit-transformed value
#' @export
#'
logistic <- function (x) exp(x)/(1 + exp(x))

#' Plots of prevalence for bayeslist
#'
#' \code{plot_prevalence.bayeslist} is used to produce plots of prevalence from a \code{bayeslist} object from the main function \code{\link{bayeslist}}.
#'
#' @param object A \code{bayeslist} object from running the main function \code{\link{bayeslist}}.
#' @param covariate_names Names of covariates.
#' @param only_prev Indicating whether only prevalence will be plotted. The default is FALSE.
#' @param xlim Limits of x-axis.
#' @param inverse Indicating whether prevalence should be calculated in the reverse order. The default is FALSE.
#' @param digit Digit number to be displayed.
#' @param ... Additional parameters to be passed.
#'
#' @importFrom grDevices rgb
#' @importFrom graphics axis hist lines mtext par points segments text
#'
#' @return None.
#' @export
#'
#'
plot_prevalence.bayeslist <- function(object,covariate_names = NULL, only_prev = FALSE, xlim = NULL, inverse = FALSE, digit = 3,...){
  old.par = par(no.readonly = TRUE)
  on.exit(par(old.par))
  if (object$npars == 1){
    only_prev = TRUE
    warning("No covariates found. Plot prevalence only.")
  }
  # Data preparation
  color1 = rgb(166,165,166,maxColorValue = 255, alpha = 100)
  # color2 = rgb(235,165,65,max = 255, alpha = 100)
  if (inverse == FALSE){
    tmp = object$X %*% t(object$sampledf[,grep("delta",names(object$sampledf))])
    tmp1 = logistic(c(tmp))
  } else {
    tmp = object$X %*% t(object$sampledf[,grep("delta",names(object$sampledf))])
    tmp1 = 1 - logistic(c(tmp))
  }
  hist1 = hist(tmp1,prob=TRUE)
  max_count = max(hist1$density)
  q1 = quantile(tmp1,probs = c(0.025,0.975))
  # variance reduction
  var1 = sd(tmp1)^2
  # Plot
  if (only_prev == TRUE){
    if (is.null(xlim)){
      xlim = c(0,1)
      plot(NA, xlim = c(0,1), ylim = c(0,max_count*1.7),axes=FALSE,xlab="",ylab="")
    } else {
      plot(NA, xlim = xlim, ylim = c(0,max_count*1.7),axes=FALSE,xlab="",ylab="")
    }
    hist(tmp1,prob=TRUE,add = T,col = color1)
    lines(density(tmp1),col = color1,lwd = 4)
    axis(1,c(0,(xlim[2])),labels = c(0,xlim[2]),las=1)
    segments(x0 = q1[1],y0 = max_count*1.1,x1 = q1[2],y1 = max_count*1.1,lwd = 4, col = color1)
    # End Segments
    segments(x0 = q1[1],y0 = max_count*1.05,x1 = q1[1],y1 = max_count*1.15,lwd = 4, col = color1)
    segments(x0 = q1[2],y0 = max_count*1.05,x1 = q1[2],y1 = max_count*1.15,lwd = 4, col = color1)
    # Add statistics: sd, efficiency gains in terms of precisions
    text(x = mean(c(q1,q1)),y = max_count*1.25,"Uninformative Prior\n95% Credible Interval",col = color1)
    text(x = q1[2]+0.1,y = max_count*1.1,paste("Var.:",round(var1,digit)))
    mtext("Prevalence",1,line=1,at=xlim[2]/2)
  } else {
    # plot
    par(mfrow=c(2,1),mai = c(0.6,0.1,0.1,0.1),xpd=TRUE)
    if (is.null(xlim)){
      xlim = c(0,1)
      plot(NA, xlim = c(0,1), ylim = c(0,max_count*1.7),axes=FALSE,xlab="",ylab="")
    } else {
      plot(NA, xlim = xlim, ylim = c(0,max_count*1.7),axes=FALSE,xlab="",ylab="")
    }
    hist(tmp1,prob=TRUE,add = T,col = color1)
    lines(density(tmp1),col = color1,lwd = 4)
    axis(1,0:1,las=1)
    segments(x0 = q1[1],y0 = max_count*1.1,x1 = q1[2],y1 = max_count*1.1,lwd = 4, col = color1)
    # End Segments
    segments(x0 = q1[1],y0 = max_count*1.05,x1 = q1[1],y1 = max_count*1.15,lwd = 4, col = color1)
    segments(x0 = q1[2],y0 = max_count*1.05,x1 = q1[2],y1 = max_count*1.15,lwd = 4, col = color1)
    # Add statistics: sd, efficiency gains in terms of precisions
    text(x = mean(c(q1,q1)),y = max_count*1.25,"95% Credible Interval",col = color1)
    text(x = q1[2]+0.1,y = max_count*1.1,paste("Var.:",round(var1,digit)))
    mtext("Prevalence",1,line=1,at=xlim[2]/2)
    ### Coefficient plot
    modout1 = object$summaryout[grep("delta",row.names(object$summaryout)),c(1,4,8)]
    k = dim(modout1)[1]
    xmin = min(c(modout1))
    xmax = max(c(modout1))
    xmean = mean(xmin,xmax)
    plot(NA, xlim = c(xmin,xmax), ylim = c(0,k),axes=FALSE,xlab="",ylab="")
    for (i in 1:k){
      segments(x0 = modout1[i,2],y0 = i,x1 = modout1[i,3],y1 = i,lwd = 4, col = color1)
      points(modout1[i,1],i,col = color1,pch = 16)
    }
    if (is.null(covariate_names)){
      for (i in 1:k){
        text(modout1[i,1],i+0.1,object$xnames[i])
      }
    } else {
      for (i in 1:k){
        text(modout1[i,1],i+0.1,covariate_names[i])
      }
    }
    axis(side = 1,tick = T)
    mtext("Coefficient Estimates", 1,line = 2)
  }
}
#' Trace plots for bayeslist
#'
#' \code{plot_trace.bayeslist} is used to produce trace plots from a \code{bayeslist} object from the main function \code{\link{bayeslist}}.
#'
#' @param object A \code{bayeslist} object from running the main function \code{\link{bayeslist}}.
#' @param ... Additional parameters to be passed to \code{\link[rstan]{traceplot}}.
#'
#' @return None.
#' @export
#'
#'
plot_trace.bayeslist <- function(object, ...) {
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
  names(object$stanfit)[-length(names(object$stanfit))] <- vnames
  rstan::traceplot(object$stanfit, ...)
}

#' Make coefficient plots for a \code{bayeslist} object
#'
#' \code{plot_coef.bayeslist} is used to produce coefficient plots from a \code{bayeslist} object.
#'
#' @param object A \code{bayeslist} object from running the main function \code{\link{bayeslist}}.
#' @param ... Additional parameters to be passed to \code{\link[rstan]{stan_plot}}.
#'
#' @return None.
#' @export
#'
#'
plot_coef.bayeslist <- function(object, ...) {
  if (object$type == "outcome") {
    vnames = c(paste0("Control items:",object$xnames),paste0("Sensitive item:",object$xnames),"Rho")
  } else if (object$type == "predict"){
    vnames = c(paste0("Control items:",object$xnames),paste0("Sensitive item:",object$xnames),"Rho",paste0("Predict:",object$xnames),"Gamma","Phi")
  } else if (object$type == "misreport") {
    vnames = c(paste0("Control items:",object$xnames),paste0("Sensitive item:",object$xnames),"Rho",paste0("Misreport:",object$xnames),"Treatment effect on misreport","U_e","Z_e")
  } else {
    stop("Unknown type of sensitive item models!")
  }
  names(object$stanfit)[-length(names(object$stanfit))] <- vnames
  rstan::plot(object$stanfit, ...)
}
