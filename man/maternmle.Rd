% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/maternmle.r
\name{maternmle}
\alias{maternmle}
\title{Older Style Interface to the Resource Utilization Function}
\usage{
maternmle(
  geodata,
  theta = NULL,
  cov = 1:ncol(geodata$covariate),
  standardized = FALSE,
  sample = 500,
  trace = 0,
  reltol = 1e-05,
  factr = 1e-05,
  results = "maternresults.out",
  birdname = NULL,
  fixrange = FALSE
)
}
\arguments{
\item{geodata}{a list containing elements `coords', `covariates', and
`data'. The componets are:

coords: an n x 2 matrix where each row has the 2-D coordinates of the n RU
locations.

data: a vector of the RU values at the n locations given by `coords'.

covariates: an n x p matrix of covariates to be fit}

\item{theta}{2-vector of Matern correlation parameters: spatial range and
smoothness.}

\item{cov}{indices of the columns of `covariates' to include in the model.}

\item{standardized}{logical: Compute standardized coeffients?}

\item{sample}{If numerical, the number of values of the data to subsample.
Default is no sub-sampling, i.e., use the full data.}

\item{trace}{level of diagnostics to print during the search.}

\item{reltol}{Relative convergence tolerance.  The algorithm stops if it is
unable to reduce the value by a factor of `reltol * (abs(val) + reltol)' at
a step.  Defaults to `sqrt(.Machine\$double.eps)', typically about `1e-8'.}

\item{factr}{controls the convergence of the `"L-BFGS-B"' method.
Convergence occurs when the reduction in the objective is within this factor
of the machine tolerance. Default is `1e7', that is a tolerance of about
`1e-8'.}

\item{results}{name of the file to store the formatted results.}

\item{birdname}{Name of the bird to report in the formatted results.}

\item{fixrange}{logical: if TRUE the range of the model is fixed at the
starting value rather than estimated.}
}
\value{
\item{par}{MLE of the spatial range and smoothness}
\item{value}{Value of the log-likelihood at the MLE.} \item{counts}{number
of log-likelihood and gradient evaluations taken.} \item{convergence}{0:
converged. See `optim' for details.} \item{message}{Message associated with
the type of convergence.} \item{hessian}{Hessian matrix of the ML
estimates.} \item{pplik}{MLE of the spatial range and smoothness}
\item{beta}{MLE of the regression coefficients} \item{value}{Value of the
log-likelihood at the MLE.} \item{asycor}{asymptotic correlation matrix of
the MLE.} \item{asyse}{asymptotic standard errors of the MLE.}
}
\description{
This is an older version of the \code{ruf.fit} function retained for legacy
purposes.  It is a function to calculate maximum likelihood fits of the
Resource Utilization Function using on a Matern covariance function.
}
\note{
The code uses the `optim' function to maximize the log-likelihood.
}
\examples{

#
# attach the small test data within the library
#
data(s412)
#
# Set initial estimates at the spatial range and smoothness
#
hval <- c(0.2, 1.5)
#
# Estimate the maximum likelihood values
# with unstandardized coefficients
#
s412.fit <- maternmle(s412,theta=hval,
         cov=-c(1,2), birdname="412",
         standardized=FALSE,
         results = "s412.out")
#
# Show the details of the results
# The formatted output has been sent to the file "s412.out"
#
s412.fit
#
# Estimate the maximum likelihood values
# with standardized coefficients
#
s412.fit <- maternmle(s412,theta=hval,
         cov=-c(1,2), birdname="412",
         standardized=TRUE,
         results = "s412.out")
s412.fit

}
\references{
``Resource utilization by an avian nest predator: relating
resources to a probabilistic measure of animal space use," by John M.
Marzluff, J. J. Millspaugh, P. Hurvitz, and Mark S. Handcock.
\emph{Ecology}, 2004, 85:1411-1427.
}
\seealso{
ruf.fit
}
