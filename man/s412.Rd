\name{s412}
\alias{s412}
\title{Data for a test Jay in `ruf'}
\description{
  Data set in `geodata' form for a Jay to test 
  maximum likelihood fits of the Resource 
  Utilization Function using on a Matern covariance function.
}
\arguments{
  \item{geodata}{a list containing elements `coords', `covariates', 
          and `data' as described

          coords: an n x 2 matrix where each row has the 2-D coordinates 
          of the n RU locations. 

          data: a vector of the RU values at the n locations given by `coords'.
          covariates: an n x p matrix of covariates to be fit}
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
data(s412)
}
\references{
``Resource utilization by an avian nest predator: relating resources to a probabilistic
measure of animal space use,"
by John M. Marzluff, J. J. Millspaugh, P. Hurvitz, and Mark S. Handcock.
\emph{Ecology}, 2004, 85:1411-1427.
}
\keyword{models}
\keyword{regression}
