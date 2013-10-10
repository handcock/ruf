\name{vmatcov}
\alias{vmatcov}
\title{Matern Correlation for the Resource Utilization Function}
\description{
  Function to calculate Matern correlation matrix for Resource
  Utilization Function using on a Matern covariance function.
}
\usage{
vmatcov(
        geodata,
        theta)
}
\arguments{
  \item{geodata}{a list containing elements `coords', `covariates', 
          and `data' as described

          coords: an n x 2 matrix where each row has the 2-D coordinates 
          of the n RU locations. 

          data: a vector of the RU values at the n locations given by `coords'.
          covariates: an n x p matrix of covariates to be fit}
  \item{theta}{2-vector of Matern correlation parameters: spatial range and
smoothness.}
}
\value{
  \item{value}{n x n Matern correlation matrix for theta}
}
\note{
  The code uses the FORTRAN code for speed.}
