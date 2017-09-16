#' Older Style Interface to the Resource Utilization Function
#' 
#' This is an older version of the \code{ruf.fit} function retained for legacy
#' purposes.  It is a function to calculate maximum likelihood fits of the
#' Resource Utilization Function using on a Matern covariance function.
#' 
#' 
#' @param geodata a list containing elements `coords', `covariates', and
#' `data'. The componets are:
#' 
#' coords: an n x 2 matrix where each row has the 2-D coordinates of the n RU
#' locations.
#' 
#' data: a vector of the RU values at the n locations given by `coords'.
#' 
#' covariates: an n x p matrix of covariates to be fit
#' @param theta 2-vector of Matern correlation parameters: spatial range and
#' smoothness.
#' @param cov indices of the columns of `covariates' to include in the model.
#' @param standardized logical: Compute standardized coeffients?
#' @param sample If numerical, the number of values of the data to subsample.
#' Default is no sub-sampling, i.e., use the full data.
#' @param trace level of diagnostics to print during the search.
#' @param reltol Relative convergence tolerance.  The algorithm stops if it is
#' unable to reduce the value by a factor of `reltol * (abs(val) + reltol)' at
#' a step.  Defaults to `sqrt(.Machine\$double.eps)', typically about `1e-8'.
#' @param factr controls the convergence of the `"L-BFGS-B"' method.
#' Convergence occurs when the reduction in the objective is within this factor
#' of the machine tolerance. Default is `1e7', that is a tolerance of about
#' `1e-8'.
#' @param results name of the file to store the formatted results.
#' @param birdname Name of the bird to report in the formatted results.
#' @param fixrange logical: if TRUE the range of the model is fixed at the
#' starting value rather than estimated.
#' @return \item{par}{MLE of the spatial range and smoothness}
#' \item{value}{Value of the log-likelihood at the MLE.} \item{counts}{number
#' of log-likelihood and gradient evaluations taken.} \item{convergence}{0:
#' converged. See `optim' for details.} \item{message}{Message associated with
#' the type of convergence.} \item{hessian}{Hessian matrix of the ML
#' estimates.} \item{pplik}{MLE of the spatial range and smoothness}
#' \item{beta}{MLE of the regression coefficients} \item{value}{Value of the
#' log-likelihood at the MLE.} \item{asycor}{asymptotic correlation matrix of
#' the MLE.} \item{asyse}{asymptotic standard errors of the MLE.}
#' @note The code uses the `optim' function to maximize the log-likelihood.
#' @seealso ruf.fit
#' @references ``Resource utilization by an avian nest predator: relating
#' resources to a probabilistic measure of animal space use," by John M.
#' Marzluff, J. J. Millspaugh, P. Hurvitz, and Mark S. Handcock.
#' \emph{Ecology}, 2004, 85:1411-1427.
#' @examples
#' 
#' #
#' # attach the small test data within the library
#' #
#' data(s412)
#' #
#' # Set initial estimates at the spatial range and smoothness
#' #
#' hval <- c(0.2, 1.5)
#' #
#' # Estimate the maximum likelihood values
#' # with unstandardized coefficients
#' #
#' s412.fit <- maternmle(s412,theta=hval,
#'          cov=-c(1,2), birdname="412",
#'          standardized=FALSE,
#'          results = "s412.out")
#' #
#' # Show the details of the results
#' # The formatted output has been sent to the file "s412.out"
#' #
#' s412.fit
#' #
#' # Estimate the maximum likelihood values
#' # with standardized coefficients
#' #
#' s412.fit <- maternmle(s412,theta=hval,
#'          cov=-c(1,2), birdname="412",
#'          standardized=TRUE,
#'          results = "s412.out")
#' s412.fit
#' 
#' @export maternmle
"maternmle"<- function(geodata, theta=NULL, 
        cov=1:ncol(geodata$covariate), standardized=FALSE,
        sample=500, trace=0, reltol=1e-5,factr=1e-5,
        results = "maternresults.out", birdname=NULL, fixrange=FALSE){
#
  geodata$coords <- as.matrix(geodata$coords)
  geodata$data <- unlist(geodata$data)
  dimnames(geodata$covariate)[[2]][1] <- "intercept"
#
  if(length(geodata$data) > 1000 & !missing(sample)){
   rrr <- sample((1:nrow(geodata$coords)),size=sample, replace=FALSE)
   print(summary(lm(geodata$data ~ -1 + as.matrix(geodata$covariate))))
   geodata <- list(covariate=geodata$covariate[rrr,],
                   coords=geodata$coords[rrr,],
                   data=geodata$data[rrr])
  }
  av <- apply(geodata$coords,2,min)
  geodata$coords <- sweep(geodata$coords, 2, av,"-")
#
  geodata$covariate <- geodata$covariate[,cov]
  degcov <- apply(geodata$covariate,2,var) > 1e-6
  covariates <- geodata$covariate
  if(standardized){
   av <- apply(covariates,2,mean)
   covariates <- sweep(covariates, 2, av, "-")
   sd <- sqrt(apply(covariates,2,var))
   covariates <- sweep(covariates, 2, sd, "/")
  }
  geodata$covariate <- cbind(1,covariates[,degcov])
#
  llik <- function(x,geodata,llik.only=TRUE){
   bad <- FALSE
   if(x[1]<0.0001 | x[2] < 0.1 | x[2] > 60){bad <- TRUE}
   if(mean(abs(x- -1))<1e-8){bad <- FALSE}
   if(bad){
    pllik <- -10000000
   }else{
    f <- as.matrix(geodata$covariate)
    zn <- geodata$data
    num <- length(geodata$data)
    vcon <- num*(1+log(2*pi))
    if(mean(abs(x- -1))<1e-8){
     v <- diag(length(zn))
    }else{
     v <- vmatcov(geodata,x)
    }
    v <- try(chol(v,pivot=FALSE))
    if(!inherits(v, "try-error")){
     vdet <- 2*sum(log(diag(v)))
     xm <- try(backsolve(v,f,transpose=TRUE))
#                   -1
#    calculate y = v  zno
#    by back-substitution
#
     Vzn <- try(backsolve(v,zn,transpose=TRUE))
     xmqr <- qr(xm)
     beta <- qr.coef(xmqr,Vzn)
     znres <- zn - f %*% beta
#
     w <- try(backsolve(v,znres,transpose=TRUE))
     pllik <- -10000000
     slope <- NA
     if(!inherits(w, "try-error")){
      slope <- mean(w*w)
      pllik <- -0.5 * ( vcon + num*log(slope) + vdet )
      if(is.na(pllik)){pllik <- -10000000}
     }
    }else{
     cat(paste("non-invertible: thta1 = ",x[1],"thta2 = ",x[2],"\n"))
     pllik <- -10000000
    }
   }
   if(llik.only){
    pllik
   }else{
    covbeta <- summary(lm(Vzn ~ -1 + xm))$cov.unscaled
    list(theta=x, pllik=pllik, beta=beta,slope=slope, covbeta=covbeta)
   }
  }
  hllik <- function(x,geodata,llik.only=TRUE,h=1){
   llik(x=c(h,x),geodata=geodata,llik.only=llik.only)
  }
#
  if(!fixrange){
  Lout <- optim(par=theta,
   fn=llik,
   hessian=TRUE,
   method="Nelder-Mead",
#  method="BFGS",
#  method="L-BFGS-B",
#  lower=c(1,0.25),upper=c(10000,60),
#  control=list(trace=trace,fnscale=-1,factr=1e15*factr),
   control=list(trace=trace,fnscale=-1, reltol=reltol),
   geodata=geodata,llik.only=TRUE
   )
   Lout$theta <- Lout$par
   ffit <- llik(Lout$theta, geodata,llik.only=FALSE)
  }else{
  Lout <- optim(par=theta[-1],
   fn=hllik,
   hessian=TRUE,
   method="Nelder-Mead",
#  method="BFGS",
#  method="L-BFGS-B",
#  lower=c(1,0.25),upper=c(10000,60),
#  control=list(trace=trace,fnscale=-1,factr=1e15*factr),
   control=list(trace=trace,fnscale=-1, reltol=reltol),
   geodata=geodata,llik.only=TRUE, h=theta[1]
   )
   Lout$theta <- c(theta[1],Lout$par)
   ffit <- llik(x=Lout$theta, geodata,llik.only=FALSE)
  }
#
   Lout$pllik <- ffit$pllik
   asycov <- diag(length(ffit$theta)+length(ffit$beta))
   acov <- try(-robust.inverse(Lout$hessian))
   if(!inherits(acov, "try-error")){
    if(!fixrange){
     asycov[1:length(ffit$theta),1:length(ffit$theta)] <- -robust.inverse(Lout$hessian)
    }else{
     asycov[2:length(ffit$theta),2:length(ffit$theta)] <- -robust.inverse(Lout$hessian)
    }
   }else{
    asycov[1:length(ffit$theta),1:length(ffit$theta)] <- 0
   }
   asycov[length(ffit$theta)+(1:length(ffit$beta)),
          length(ffit$theta)+(1:length(ffit$beta))] <- ffit$covbeta
   dimnames(asycov) <- list(c(names(ffit$theta),names(ffit$beta)),
                            c(names(ffit$theta),names(ffit$beta)))
   asyse <- sqrt(diag(asycov))
   asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
   dimnames(asycor) <- dimnames(asycov)
   Lout$beta <- ffit$beta
   Lout$asycor <- asycor
   Lout$asycor <- asycor
   Lout$asyse <- asyse
#
   ntheta <- Lout$theta-Lout$theta-1
   ffit <- llik(ntheta, geodata, llik.only=FALSE)
   asycov <- diag(length(ffit$theta)+length(ffit$beta))
   asycov[1:length(ffit$theta),1:length(ffit$theta)] <- NA
   asycov[length(ffit$theta)+(1:length(ffit$beta)),
          length(ffit$theta)+(1:length(ffit$beta))] <- ffit$covbeta
   dimnames(asycov) <- list(c(names(ffit$theta),names(ffit$beta)),
                            c(names(ffit$theta),names(ffit$beta)))
   asyse <- sqrt(diag(asycov))
   asycor <- diag(1/asyse) %*% asycov %*% diag(1/asyse)
   dimnames(asycor) <- dimnames(asycov)
#
   cllik <-  -2*(ffit$pllik-Lout$pllik)
   sink(file=results,append=TRUE)
   if(standardized){
    cat(paste("\nStandardized Coefficients for bird name:",birdname,"\n"),file=results,append=TRUE)
  }else{
    cat(paste("\nUnstandardized Coefficients for bird name:",birdname,"\n"),file=results,append=TRUE)
  }
   cat(paste("\nMatern Log-Lik =",Lout$pllik,"LS Log-Lik =",ffit$pllik,"\n"),file=results,append=TRUE)
   cat(paste("\nChange in Log-Lik",0.5*cllik,"p-value =",1-pchisq(cllik,2),"\n"),file=results,append=TRUE)
   cat(paste("\n"),file=results,append=TRUE)
#
   if(!fixrange){
    summarylik <- cbind(c(NA,NA,ffit$beta),c(Lout$theta,Lout$beta),
                       Lout$asyse,asyse)
   }else{
    summarylik <- cbind(c(NA,NA,ffit$beta),c(Lout$theta,Lout$beta),
                       Lout$asyse,asyse)
    summarylik[1,3] <- NA
   }
   dimnames(summarylik) <- list(c("range","smoothness",
                                   dimnames(geodata$covariate)[[2]]),
                        c("LS estimate","MLE","s.e.","LS s.e."))
   print(round(summarylik,6),file=file,append=TRUE)
   sink()
#
   Lout
}
#
# get the Matern covariance
#


#' Matern Correlation for the Resource Utilization Function
#' 
#' Function to calculate Matern correlation matrix for Resource Utilization
#' Function using on a Matern covariance function.
#' 
#' 
#' @param geodata a list containing elements `coords', `covariates', and `data'
#' as described
#' 
#' coords: an n x 2 matrix where each row has the 2-D coordinates of the n RU
#' locations.
#' 
#' data: a vector of the RU values at the n locations given by `coords'.
#' covariates: an n x p matrix of covariates to be fit
#' @param theta 2-vector of Matern correlation parameters: spatial range and
#' smoothness.
#' @return \item{value}{n x n Matern correlation matrix for theta}
#' @note The code uses the FORTRAN code for speed.
#' @export vmatcov
"vmatcov"<- function(geodata,theta) {
  N <- length(geodata$data)
  coords <- geodata$coords
  cov <- array(0, dim = c(N,N))
  storage.mode(cov) <- "double"
  storage.mode(theta) <- "double"
  storage.mode(coords) <- "double"
  .Fortran("matcov",
        cov,
        coords,
        theta,
        as.integer(N),
   PACKAGE="ruf")[[1]]
}
#' @keywords internal
"vmatcovform"<- function(y,coords,theta) {
  N <- length(y)
  cov <- array(0, dim = c(N,N))
  storage.mode(cov) <- "double"
  storage.mode(theta) <- "double"
  storage.mode(coords) <- "double"
  .Fortran("matcov",
        cov,
        coords,
        theta,
        as.integer(N),
   PACKAGE="ruf")[[1]]
}
