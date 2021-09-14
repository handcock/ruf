#' Software to Implement Resource Utilization Function Estimation
#'
#' This package provides functionality 
#' to determine maximum
#' likelihood fits of the Resource Utilization Function based on a
#' Matern covariance function.
#' @docType package
#' @name ruf
#' @useDynLib ruf, .registration = TRUE
#' @importFrom utils packageDescription
#' @importFrom stats contrasts dist is.empty.model lm model.matrix model.offset model.response model.weights optim pchisq var
NULL

#' Data for a test Jay
#' @name d412
#' @docType data
#' @keywords datasets
#' @format An data.frame object
#' 
#' @description Data set in standard \code{data.frame} form for a Jay to test maximum likelihood
#' fits of the Resource Utilization Function using on a Matern covariance
#' function.
#' 
#' data.frame: a data.frame containing colums for the spatial
#' coordinates, covariates, and RUF values as described. The variables can have
#' any names and are specified in the formula and 'space' formula in the
#' 'ruffit' call.
#' 
#' coordinates: Two variables where each row has the 2-D coordinates of the n
#' RU locations.
#' 
#' RUF: a vector of the RU values at the n locations given by the coordinatess.
#' covariates: a set of p vectors of covariates to be fit
#' @note See the `ruf' library for details.
#' @references ``Resource utilization by an avian nest predator: relating
#' resources to a probabilistic measure of animal space use," by John M.
#' Marzluff, J. J. Millspaugh, P. Hurvitz, and Mark S. Handcock.
#' \emph{Ecology}, 2004, 85:1411-1427.
#' @keywords models regression
#' @examples
#' #
#' # attach the small test data within the library
#' #
#' data(d412)
#' 
NULL

#' Data for a test Jay in `ruf'
#' @name s412
#' @docType data
#' @keywords datasets
#' @format An data.frame object
#' 
#' @description Data set in `geodata' form for a Jay to test maximum likelihood fits of the
#' Resource Utilization Function using on a Matern covariance function.
#' 
#' geodata: a list containing elements `coords', `covariates', and `data'
#' as described
#' 
#' coords: an n x 2 matrix where each row has the 2-D coordinates of the n RUF
#' locations.
#' 
#' data: a vector of the RU values at the n locations given by `coords'.
#' covariates: an n x p matrix of covariates to be fit
#' @note See the `ruf' library for details.
#' @references ``Resource utilization by an avian nest predator: relating
#' resources to a probabilistic measure of animal space use," by John M.
#' Marzluff, J. J. Millspaugh, P. Hurvitz, and Mark S. Handcock.
#' \emph{Ecology}, 2004, 85:1411-1427.
#' @keywords models regression
#' @examples
#' #
#' # attach the small test data within the library
#' #
#' data(s412)
#' 
NULL
