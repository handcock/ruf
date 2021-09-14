#' Estimates of the Resource Utilization Function
#' 
#' Function to calculate maximum likelihood fits of the Resource Utilization
#' Function using on a Matern covariance function.
#' 
#' 
#' @param formula An formula object, of the form \code{RU ~ <covariate terms>},
#' where \code{RU} is the name of the RU values.  The \code{<model terms>} are
#' a list of variables that are covariates.  All variables must be in the
#' data.frame specified by the \code{data} argument.
#' @param space An formula object, of the form \code{~ x + y}, where \code{x}
#' and \code{y} are the names of the variables containing the two coordinates
#' of the RU variable.
#' @param data a data.frame containing columns for the spatial coordinates,
#' covariates, and RUF values as described. The variables can have any names
#' and are specified in the formula and 'space' formula in the \code{ruf.fit}
#' call.
#' 
#' coordinates: Two variables where each row has the 2-D coordinates of the n
#' RU locations.
#' 
#' RUF: a vector of the RU values at the n locations given by the coordinates.
#' covariates: a set of p vectors of covariates to be fit
#' @param subset an optional vector specifying a subset of observations to be
#' used in the fitting process.
#' @param na.action a function which indicates what should happen when the data
#' contain \code{NA}s.  The default is set by the \code{na.action} setting of
#' \code{\link{options}}, and is \code{\link{na.fail}} if that is unset.  The
#' \dQuote{factory-fresh} default is \code{\link{na.omit}}.
#' @param theta 2-vector of Matern correlation parameters: spatial range and
#' smoothness.
#' @param standardized logical: Compute standardized coefficients?
#' @param algorithm.control An optional list of control parameters. See the
#' \bold{Algorithm Tuning} section below for the options.
#' @param name Name of the bird to report in the formatted results.
#' @param fixrange logical: if TRUE the range of the model is fixed at the
#' starting value rather than estimated.
#' @param fixsmoothness logical: if TRUE the smoothness of the model is fixed
#' at the starting value rather than estimated.
#' @param \dots further arguments passed to or from other methods.
#' @return \item{par}{MLE of the spatial range and smoothness}
#' \item{value}{Value of the log-likelihood at the MLE.} \item{counts}{number
#' of log-likelihood and gradient evaluations taken.} \item{convergence}{0:
#' converged. See `optim' for details.} \item{message}{Message associated with
#' the type of convergence.} \item{hessian}{Hessian matrix of the ML
#' estimates.} \item{pplik}{MLE of the spatial range and smoothness}
#' \item{beta}{MLE of the regression coefficients} \item{value}{Value of the
#' log-likelihood at the MLE.} \item{asycor}{asymptotic correlation matrix of
#' the MLE.} \item{asyse}{asymptotic standard errors of the MLE.}
#' @note The code uses the \code{optim} function to maximize the
#' log-likelihood.  If the \code{fixrange} or \code{fixsmoothness} options are
#' used the one-dimensional \code{optimize} is used instead.
#' @section Algorithm Tuning:
#' 
#' There are a large number of parameters that can modified the computational
#' aspects. These are specified via the \code{algorithm.control} argument. The
#' \code{algorithm.control} argument is a list that can supply any of the
#' following components: \describe{ \item{list("maxit")}{count; The maximum
#' number of iterations in the Newton-Raphson optimization.  Defaults to
#' \code{15}.  \code{maxit} gives the total number of likelihood function
#' evaluations.} \item{sample}{The number of values of the data to subsample
#' for the re-sampled MLE method for large data sets. Default is no
#' sub-sampling, i.e., use the full data unless the data set has over 2000
#' points in it.} \item{nresamples}{The number of re-samples of the full
#' dataset to take. This is for use with the re-sampled MLE method for large
#' data sets. Default is no sub-sampling, i.e., use the full data unless the
#' data set has over 2000 points in it.} \item{reltol}{Relative convergence
#' tolerance.  The algorithm stops if it is unable to reduce the value by a
#' factor of `reltol * (abs(val) + reltol)' at a step.  Defaults to
#' `sqrt(.Machine\$double.eps)', typically about `1e-8'.} \item{factr}{controls
#' the convergence of the `"L-BFGS-B"' method.  Convergence occurs when the
#' reduction in the objective is within this factor of the machine tolerance.
#' Default is `1e7', that is a tolerance of about `1e-8'.}
#' \item{list("trace")}{non-negative integer; If positive, tracing information
#' on the progress of the optimization is produced. Higher values may produce
#' more tracing information: for method \code{"L-BFGS-B"} there are six levels
#' of tracing.  (To understand exactly what these do see the source code for
#' \code{\link[stats]{optim}}: higher levels give more detail.)}
#' \item{list("method")}{character; The name of the optimization method to use
#' for the maximum likelihood estimation. See \code{\link[stats]{optim}} for
#' the options.  The default method for the maximixation of the smoothness and
#' range is \code{"BFGS"}, a quasi-Newton method (also known as a variable
#' metric algorithm). It is attributed to Broyden, Fletcher, Goldfarb and
#' Shanno. This uses function values and gradients to build up a picture of the
#' surface to be optimized.  The default method for one-dimensional searchs is
#' that of Nelder and Mead \code{Nelder-Mead}.} }
#' @references ``Resource utilization by an avian nest predator: relating
#' resources to a probabilistic measure of animal space use," by John M.
#' Marzluff, J. J. Millspaugh, P. Hurvitz, and Mark S. Handcock.
#' \emph{Ecology}, 2004, 85:1411-1427.
#' @keywords models regression
#' @examples
#' 
#' #
#' # attach the small test data within the library
#' #
#' data(d412)
#' #
#' # Set initial estimates at the spatial range and smoothness
#' #
#' hval <- c(0.2, 1.5)
#' #
#' # Estimate the maximum likelihood values
#' # with unstandardized coefficients
#' #
#' d412.fit <- ruf.fit(ruf ~ CWED + IJI + NP + MSI,
#'          space= ~ x + y,
#'          data=d412, theta=hval,
#'          name="Bird 412",
#'          standardized=FALSE)
#' #
#' # Show the details of the results
#' #
#' summary(d412.fit)
#' #
#' # Estimate the maximum likelihood values
#' # with standardized coefficients
#' #
#' d412.fit <- ruf.fit(ruf ~ CWED + IJI + NP + MSI,
#'          space= ~ x + y,
#'          data=d412, theta=hval,
#'          name="Bird 412 standardized",
#'          standardized=TRUE)
#' summary(d412.fit)
#' 
#' @export ruf.fit
"ruf.fit" <- function (formula, space, data, subset, na.action,
    theta=NULL,
    standardized=FALSE,
    algorithm.control=list(),
    name="", fixrange=FALSE, fixsmoothness=TRUE, ...) 
{
    cl <- match.call()
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action"), 
               names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
##  Defaults :
    con <- list(maxit=20, method="Nelder-Mead", #   method="BFGS",
                trace=0, reltol=1e-8, factr=1e-8,
                sample=NULL, nresamples=10
             )
    con[(namc <- names(algorithm.control))] <- algorithm.control
#   trms <- ruf.getterms(formula)
#   termnames <- ruf.gettermnames(trms)
#   y <- try(as.vector(eval(trms[[2]],data)))
#   if(inherits(y,"try-error")){
#    stop("Invalid ruf. Is the left-hand-side of the formula correct?")
#   }
#   m <- ruf.getmodel(trms, y)
    mf <- eval(mf, parent.frame())
    mt <- attr(mf, "terms")
    y <- model.response(mf, "numeric")
    w <- model.weights(mf)
    offset <- model.offset(mf)
#
    coords <- model.matrix(object=space, data=data)[,-1]
#
    if (!is.null(offset) && length(offset) != NROW(y)) 
        stop("Number of offsets is ", length(offset), ", should equal ", 
            NROW(y), " (number of observations)")
    if (is.empty.model(mt)) {
        x <- NULL
        z <- list(coefficients = numeric(0), residuals = y, fitted.values = 0 * 
            y, weights = w, rank = 0, df.residual = length(y))
        if (!is.null(offset)) 
            z$fitted.values <- offset
    }
    else {
        x <- model.matrix(object=mt, data=mf)
    }
#
# Back to the code
#
  if(length(y) > 5000 | is.numeric(con$sample)){
   if(!is.numeric(con$sample)){con$sample <- 2000}
   warning("Estimating using a resampled MLE. The standard errors include the variation due to the resampling.")
   betas <- matrix(0,ncol=ncol(x),nrow=con$nresamples)
   thetas <- matrix(0,ncol=length(theta),nrow=con$nresamples)
   colnames(thetas) <- c("range","smoothness")
   colnames(betas) <- colnames(x)
   asycovtheta <- diag(length(theta))
   asycovbeta <- diag(ncol(x))
   for(i in 1:con$nresamples){
     samplesome <- sample((1:nrow(coords)),size=con$sample, replace=FALSE)
     y <- y[samplesome]
     x <- x[samplesome, , drop = FALSE]
     coords <- coords[samplesome,]
     av <- apply(coords,2,min)
     coords <- sweep(coords, 2, av,"-")
#
     z <- ruf.fit.mainloop(coords, x, y, theta, standardized, con, fixrange, fixsmoothness)
     if(i ==1){
       asycovtheta <- robust.inverse(-z$hessian)
       asycovbeta <- z$fit$covbeta
     }else{
       asycovtheta <- asycovtheta + robust.inverse(-z$hessian)
       asycovbeta <- asycovbeta + z$fit$covbeta
     }
     betas[i,] <- z$fit$beta
     thetas[i,] <- z$fit$theta
     cat(paste("Finished ",round(100*i/con$nresamples),"% of the fitting.\n",sep=""))
   }
#
   z$fit$beta <- apply(betas,2,mean)
   z$fit$theta <- apply(thetas,2,mean)
   if(!fixrange & !fixsmoothness){
     z$fit$hessian <- -robust.inverse(var(thetas) + asycovtheta/con$nresamples)
   }else{
     if(!fixsmoothness){
      asycovtheta <- asycovtheta/con$nresamples + var(thetas)[2:length(z$fit$theta),2:length(z$fit$theta)]
      z$fit$hessian <- -robust.inverse(asycovtheta)
     }else{
      asycovtheta <- asycovtheta/con$nresamples + var(thetas)[(1:length(z$fit$theta))[-2],(1:length(z$fit$theta))[-2]]
      z$fit$hessian <- -robust.inverse(asycovtheta)
     }
   }
   z$fit$covbeta <- var(betas) + asycovbeta/con$nresamples
  }else{
    z <- ruf.fit.mainloop(coords, x, y, theta, standardized, con, fixrange, fixsmoothness)
  }
#
   ffit <- z$fit
   z$pllik <- ffit$pllik
   asycov <- diag(length(ffit$theta)+length(ffit$beta))
   diag(asycov) <- NA
   acov <- try(-robust.inverse(z$hessian),silent=TRUE)
   if(!inherits(acov, "try-error")){
    if(!fixrange & !fixsmoothness){
     asycov[1:length(ffit$theta),1:length(ffit$theta)] <- -robust.inverse(z$hessian)
    }else{
     if(!fixsmoothness){
      asycov[2:length(ffit$theta),2:length(ffit$theta)] <- -robust.inverse(z$hessian)
     }else{
      asycov[(1:length(ffit$theta))[-2],(1:length(ffit$theta))[-2]] <- -robust.inverse(z$hessian)
     }
    }
   }else{
    asycov[1:length(ffit$theta),1:length(ffit$theta)] <- 0
   }
   asycov[length(ffit$theta)+(1:length(ffit$beta)),
          length(ffit$theta)+(1:length(ffit$beta))] <- ffit$covbeta
   dimnames(asycov) <- list(c(names(ffit$theta),names(ffit$beta)),
                            c(names(ffit$theta),names(ffit$beta)))
   asyse <- sqrt(diag(asycov))
#  asycor <- cov2cor(asycov)
   asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
   dimnames(asycor) <- dimnames(asycov)
   z$beta <- ffit$beta
   z$asycor <- asycor
   z$asycov <- asycov
   z$asyse <- asyse
#
   ntheta <- theta-theta-1
   ffit <- ruf.llik(ntheta, y, x,coords, llik.only=FALSE)
   asycov <- diag(length(ffit$theta)+length(ffit$beta))
   asycov[1:length(ffit$theta),1:length(ffit$theta)] <- NA
   asycov[length(ffit$theta)+(1:length(ffit$beta)),
          length(ffit$theta)+(1:length(ffit$beta))] <- ffit$covbeta
   dimnames(asycov) <- list(c(names(ffit$theta),names(ffit$beta)),
                            c(names(ffit$theta),names(ffit$beta)))
   asyse <- sqrt(diag(asycov))
   asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
#  asycor <- cov2cor(asycov)
   dimnames(asycor) <- dimnames(asycov)

    class(z) <- "ruf"
    z$na.action <- attr(mf, "na.action")
    z$x <- x
    z$lsfit <- ffit
    z$lsasyse <- asyse
    z$offset <- offset
    z$standardized <- standardized
    z$name <- name
    z$fixrange <- fixrange
    z$fixsmoothness <- fixsmoothness
    z$contrasts <- attr(x, "contrasts")
    z$call <- cl
    z$terms <- mt
    cat(paste("Fitting completed successfully!\n"))
    z
}
"ruf.fit.mainloop" <- function (coords, x, y,
    theta=NULL,
    standardized=FALSE,
    con,
    fixrange=FALSE, fixsmoothness=TRUE) 
{
  av <- apply(coords,2,min)
  coords <- sweep(coords, 2, av,"-")
#
# degcov <- apply(x,2,var) > 1e-6
  if(standardized){
   not.intercept <- is.na(match(colnames(x),"(Intercept)"))
   av <- apply(x[,not.intercept,drop=FALSE],2,mean)
   x[,not.intercept] <- sweep(x[,not.intercept,drop=FALSE], 2, av, "-")
   sd <- sqrt(apply(x[,not.intercept,drop=FALSE],2,var))
   x[,not.intercept] <- sweep(x[,not.intercept,drop=FALSE], 2, sd, "/")
  }
#
#
  if(is.null(theta)){
   theta <- c(0.5*sqrt(max(apply(coords,2,var))),2)
  }
  drange <- range(dist(coords))
  z <- list(theta=theta)
  if(!fixrange & !fixsmoothness){
   z <- stats::optim(par=theta,
    fn=ruf.llik,
    hessian=TRUE,
#   method=method,
#   method="BFGS",
#   control=list(trace=con$trace,fnscale=-1, reltol=con$reltol),
    method="L-BFGS-B",
    lower=c(drange[1],0.5),upper=c(drange[2]/2,10),
    control=list(maxit=con$maxit,trace=con$trace,fnscale=-1,factr=1e15*con$factr),
    y=y,x=x,coords=coords,llik.only=TRUE
   )
   z$theta <- z$par
   names(z$theta) <- c("range","smoothness")
  }else{
   if(!fixrange){
    options(warn=-1)
    z <- stats::optim(par=theta[-2],
     fn=ruf.vllik,
     hessian=TRUE,
     method=con$method,
     lower=drange[1],upper=drange[2]/3,
     control=list(maxit=con$maxit,trace=con$trace,fnscale=-1, reltol=con$reltol),
      y=y,x=x,coords=coords,llik.only=TRUE, v=theta[2]
    )
    options(warn=0)
    z$theta <- c(z$par,theta[2])
    names(z$theta) <- c("range","smoothness")
#  z <- stats::optimize(#par=theta[-2],
#   f=ruf.vllik,
#   interval=c(drange[1],drange[2]/3),
#   maximum=TRUE, tol=sqrt(con$reltol),
#   y=y,x=x,coords=coords,llik.only=TRUE, v=theta[2]
#   )
#   z$theta <- c(z$maximum,theta[2])
   }
   if(!fixsmoothness){
    z <- stats::optim(par=theta[-1],
     fn=ruf.hllik,
     hessian=TRUE,
     method=con$method,
##  method="BFGS",
##  method="L-BFGS-B",
##  lower=c(1,0.5),upper=c(10000,10),
     lower=0.5,upper=10,
##  control=list(trace=con$trace,fnscale=-1,factr=1e15*con$factr),
     control=list(maxit=con$maxit,trace=con$trace,fnscale=-1, reltol=con$reltol),
     y=y,x=x,coords=coords,llik.only=TRUE, h=theta[1]
    )
    z$theta <- c(theta[1],z$par)
    names(z$theta) <- c("range","smoothness")
#   z <- stats::optimize(#par=theta[-1],
#    f=ruf.hllik,
#    interval=c(0.5,10),
#    maximum=TRUE, tol=sqrt(con$reltol),
#    y=y,x=x,coords=coords,llik.only=TRUE, h=theta[1]
#    )
#    z$theta <- c(z$maximum,theta[2])
   }
   }
#
   ffit <- ruf.llik(z$theta, y,x,coords,llik.only=FALSE)
   z$fit <- ffit
   z
}
#
#   Likelihood functions
#
  ruf.llik <- function(theta,y,x,coords,llik.only=TRUE){
   bad <- FALSE
   if(theta[1]< 0 | theta[2] < 0.1 | theta[2] > 11){bad <- TRUE}
   if(mean(abs(theta- -1))<1e-8){bad <- FALSE}
   if(bad){
    pllik <- -10000000
   }else{
    f <- x
    zn <- y
    num <- length(y)
    vcon <- num*(1+log(2*pi))
    if(mean(abs(theta- -1))<1e-8){
     v <- diag(length(zn))
     vdet <- 0
     xm <- f
     Vzn <- zn
     beta <- qr.coef(qr(f),zn)
     w <- zn - f %*% beta
    }else{
     v <- vmatcovform(y,coords,theta)
     v <- try(chol(v,pivot=FALSE),silent=TRUE)
     if(!inherits(v, "try-error")){
      vdet <- 2*sum(log(diag(v)))
      xm <- try(backsolve(v,f,transpose=FALSE),silent=TRUE)
#                   -1
#     calculate y = v  zno
#     by back-substitution
#
      Vzn <- try(backsolve(v,zn,transpose=FALSE),silent=TRUE)
      xmqr <- qr(xm)
      beta <- qr.coef(xmqr,Vzn)
      znres <- zn - f %*% beta
#
      w <- try(backsolve(v,znres,transpose=FALSE),silent=TRUE)
     }
    }
    if(!inherits(v, "try-error")){
     pllik <- -10000000
     slope <- NA
     if(!inherits(w, "try-error")){
      slope <- mean(w*w)
      pllik <- -0.5 * ( vcon + num*log(slope) + vdet )
      if(is.na(pllik)){pllik <- -10000000}
     }
    }else{
     pllik <- -10000000
     cat(paste("Note: non-invertible: range = ",theta[1],"smoothness = ",theta[2],". This should not be a problem. Continuing the fit...\n"),sep="")
    }
   }
#
# Penalty
#
#  pllik <- pllik - max(c(0,theta[2] - 20))
#  pllik <- pllik - max(c(0,theta[2] - 5))
#  pllik <- pllik - 2*log(1+theta[2])
#  pllik <- pllik - 2*(log(1+theta[2])-log(1+10))*(theta[2]>10)
   if(mean(abs(theta- -1))>1e-8){
     pllik <- pllik - 2*(log(1+theta[2])-log(1+2))*(theta[2]>2)
   }
#
   if(llik.only | bad){
    pllik
   }else{
    if(!exists("Vzn")){
     Vzn <- zn
     xm <- f
     slope <- 0
     beta <- rep(0,ncol(xm))
    }
    names(beta) <- colnames(x)
    names(theta) <- c("range","smoothness")
    covbeta <- summary(stats::lm(Vzn ~ -1 + xm))$cov.unscaled
    list(theta=theta, pllik=pllik, beta=beta,slope=slope, covbeta=covbeta)
   }
  }
  ruf.hllik <- function(theta,y,x,coords,llik.only=TRUE,h=1){
   ruf.llik(theta=c(h,theta),y=y,x=x,coords=coords,llik.only=llik.only)
  }
  ruf.vllik <- function(theta,y,x,coords,llik.only=TRUE,v=2){
   ruf.llik(theta=c(theta,v),y=y,x=x,coords=coords,llik.only=llik.only)
  }
