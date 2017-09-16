#' Summarizing RUF Model Fits
#' 
#' These functions are all \code{\link{methods}} for class \code{ruf} or
#' \code{summary.ruf} objects.
#' 
#' \code{summary.ruf} tries to be smart about formatting the coefficients,
#' standard errors, etc. and additionally gives \dQuote{significance stars} if
#' \code{signif.stars} is \code{TRUE}.
#' 
#' @aliases summary.ruf print.summary.ruf
#' @param object an object of class \code{"ruf"}, usually, a result of a call
#' to \code{\link{ruf.fit}}.
#' @param results name of the file to store the formatted results. The defualt
#' is to the console.
#' @param \dots further arguments passed to or from other methods.
#' @return \code{summary.ruf} returns an nothing.
#' @seealso \code{\link{ruf.fit}}, \code{\link{summary}}.
#' @references ``Resource utilization by an avian nest predator: relating
#' resources to a probabilistic measure of animal space use," by John M.
#' Marzluff, J. J. Millspaugh, P. Hurvitz, and Mark S. Handcock.
#' \emph{Ecology}, 2004, 85:1411-1427.
#' @keywords models regression
#' @examples
#' 
#' ## --- Continuing the Example from  '?ruf.fit':%\code{\link{ruf.fit}}:
#' \dontrun{
#' summary(d412.fit)
#' }
#' 
#' @method summary ruf
summary.ruf <- function(object, results = "", ...){
   if(!missing(results)){sink(file=results,append=TRUE)}
   cllik <-  -2*(object$lsfit$pllik-object$fit$pllik)
   if(object$standardized){
    cat(paste("\nStandardized Coefficients for name:",object$name,"\n"),file=results,append=TRUE)
  }else{
    cat(paste("\nUnstandardized Coefficients for name:",object$name,"\n"),file=results,append=TRUE)
  }
   cat(paste("\nMatern Log-Lik =",format(object$pllik,2),
                   "LS Log-Lik =",format(object$lsfit$pllik,2),"\n"),
       file=results,append=TRUE)
   cat(paste("\nChange in Log-Lik",format(0.5*cllik,2),
      "p-value =",format.pval(1-pchisq(cllik,2),eps=0.0001),"\n"),
       file=results,append=TRUE)
   cat(paste("\n"),file=results,append=TRUE)
#
   if(!object$fixrange){
    summarylik <- cbind(c(NA,NA,object$lsfit$beta),c(object$theta,object$beta),
                       object$asyse,object$lsasyse)
   }else{
    summarylik <- cbind(c(NA,NA,object$lsfit$beta),c(object$theta,object$beta),
                       object$asyse,object$lsasyse)
    summarylik[1,3] <- NA
   }
   dimnames(summarylik) <- list(c("range","smoothness",dimnames(object$x)[[2]]),
                                c("LS estimate","MLE","s.e.","LS s.e."))
   summarylik <- summarylik[,c(2,3,1,4)]
   print(round(summarylik,6),file=file,append=TRUE)
   if(!missing(results)){sink()}
}
