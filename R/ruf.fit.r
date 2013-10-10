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
    coords <- model.matrix(space, data)[,-1]
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
        x <- model.matrix(mt, mf, contrasts)
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
   z <- optim(par=theta,
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
    z <- optim(par=theta[-2],
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
#  z <- optimize(#par=theta[-2],
#   f=ruf.vllik,
#   interval=c(drange[1],drange[2]/3),
#   maximum=TRUE, tol=sqrt(con$reltol),
#   y=y,x=x,coords=coords,llik.only=TRUE, v=theta[2]
#   )
#   z$theta <- c(z$maximum,theta[2])
   }
   if(!fixsmoothness){
    z <- optim(par=theta[-1],
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
#   z <- optimize(#par=theta[-1],
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
    covbeta <- summary(lm(Vzn ~ -1 + xm))$cov.unscaled
    list(theta=theta, pllik=pllik, beta=beta,slope=slope, covbeta=covbeta)
   }
  }
  ruf.hllik <- function(theta,y,x,coords,llik.only=TRUE,h=1){
   ruf.llik(theta=c(h,theta),y=y,x=x,coords=coords,llik.only=llik.only)
  }
  ruf.vllik <- function(theta,y,x,coords,llik.only=TRUE,v=2){
   ruf.llik(theta=c(theta,v),y=y,x=x,coords=coords,llik.only=llik.only)
  }
