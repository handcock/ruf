\name{d412}
\alias{d412}
\title{Data for a test Jay in `ruf'}
\description{
  Data set in standard data.frame form for a Jay to test 
  maximum likelihood fits of the Resource 
  Utilization Function using on a Matern covariance function.
}
\arguments{
  \item{data.frame}{a data.frame containing colums for the spatial coordinates, covariates, 
          and RUF values as described. The variables can have any names and are
specified in the formula and 'space' formula in the 'ruffit' call.

          coordinatess: Two variables where each row has the 2-D coordinates 
          of the n RU locations. 

          RUF: a vector of the RU values at the n locations given by the 
 coordinatess.
          covariates: a set of p vectors of covariates to be fit}
}
\note{
  See the `ruf' library for details.}
\examples{
#
# attach the "ruf" library
#
library(ruf)
#
# attach the small test data within the library
#
data(d412)
}
\references{
``Resource utilization by an avian nest predator: relating resources to a probabilistic
measure of animal space use,"
by John M. Marzluff, J. J. Millspaugh, P. Hurvitz, and Mark S. Handcock.
\emph{Ecology}, 2004, 85:1411-1427.
}
\keyword{models}
\keyword{regression}
