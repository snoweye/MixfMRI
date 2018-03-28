\name{Smoothing}
\alias{gcv.smooth2d}
\title{Generate datasets with smoothing for MixfMRI simulations}
\description{
  Generate datasets with smoothing for MixfMRI simulations
}
\usage{
  gcv.smooth2d(y, interval)
}
\arguments{
  \item{y}{a set of p-values in 2d phantom}
  \item{interval}{an interval for \code{optimize} function.}
}
\details{
  The function is used to smooth for Dr. Maitra's 2d phantom simulation.
  The smoothing method is based on Garcia (2010), CSDA.
  
}
\value{
  Return a list containing two elements \code{im.smooth} and
  \code{par.val}.
}
\references{
  \url{http://maitra.public.iastate.edu/}
}
\author{
  Ranjan Maitra.
}
\keyword{programming}

