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
  Chen, W.-C. and Maitra, R. (2021)
  \dQuote{A Practical Model-based Segmentation Approach for Accurate
	  Activation Detection in Single-Subject functional Magnetic
	  Resonance Imaging Studies},
  \emph{arXiv:2102.03639}.
}
\author{
  Ranjan Maitra.
}
\keyword{programming}

