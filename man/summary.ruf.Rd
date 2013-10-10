\name{summary.ruf}
\alias{summary.ruf}
\alias{print.summary.ruf}
\title{Summarizing RUF Model Fits}
\usage{
\method{summary}{ruf}(object, results = "", \dots)
}
\arguments{
  \item{object}{an object of class \code{"ruf"}, usually, a result of a
    call to \code{\link{ruf.fit}}.}
  \item{results}{name of the file to store the formatted results. The defualt
is to the console.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\description{
  These functions are all \code{\link{methods}} for class \code{ruf} or
  \code{summary.ruf} objects.
}
\details{
  \code{summary.ruf} tries to be smart about formatting the
  coefficients, standard errors, etc. and additionally gives
  \dQuote{significance stars} if \code{signif.stars} is \code{TRUE}.
}
\value{
  \code{summary.ruf} returns an nothing.
}
\seealso{
  \code{\link{ruf.fit}}, \code{\link{summary}}.
}
\references{
``Resource utilization by an avian nest predator: relating resources to a probabilistic
measure of animal space use,"
by John M. Marzluff, J. J. Millspaugh, P. Hurvitz, and Mark S. Handcock.
\emph{Ecology}, 2004, 85:1411-1427.
}
\examples{
## --- Continuing the Example from  '?ruf.fit':%\code{\link{ruf.fit}}:
%\dontshow{example("ruf", echo = FALSE)}
\dontrun{
summary(d412.fit)
}
}
\keyword{models}
\keyword{regression}
