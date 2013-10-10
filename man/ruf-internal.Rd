\name{ruf-internal}
\title{Internal ruf Objects}

\alias{robust.inverse}
\alias{vmatcovform}
\alias{ruf.fit.mainloop}
\alias{ruf.hllik}
\alias{ruf.llik}
\alias{ruf.vllik}

\description{
  Internal ruf functions.
}
\usage{
robust.inverse(H, tol = sqrt(.Machine$double.eps))
}
\details{
  Most of these are not to be called by the user.
}
\arguments{
\item{H}{Square matrix to be inverted.}
\item{tol}{Tolerance for the inversion.}
%\item{\dots}{further arguments passed to or used by methods.}
}
\seealso{ruf}
\references{
``Resource utilization by an avian nest predator: relating resources to a probabilistic
measure of animal space use,"
by John M. Marzluff, J. J. Millspaugh, P. Hurvitz, and Mark S. Handcock.
\emph{Ecology}, 2004, 85:1411-1427.
}
\keyword{internal}
