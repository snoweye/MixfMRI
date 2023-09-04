\name{Threshold.RF}
\alias{Threshold.RF}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Random Field Theory Thersholds.}
\description{Calculates the Random Field theory threshold to give that
  results in a specified p-value.
}
\usage{
Threshold.RF(p.val, sigma, voxdim = c(1, 1, 1), num.vox,
             type = c("Normal", "t"), df = NULL)
}
\arguments{
  \item{p.val}{The required p-value.}
  \item{sigma}{The 3D covariance matrix of the random field.}
  \item{voxdim}{The dimesnions of a voxel.}
  \item{num.vox}{The number of voxels that constitute the random field.}
  \item{type}{The type of random field, "Normal" or "t".}
  \item{df}{The degrees of the t distributed field.}
}
\details{Calculates the threshold that produces an expected Euler
  characteristic equal to the required p-value.
}
\value{Returns the Random Field threshold.
}
\author{J. L. Marchini}

\seealso{\code{\link{EC.3D}}}

\examples{
a <- Threshold.RF(p.val = 0.05, sigma = diag(1, 3), voxdim = c(1, 1, 1), num.vox = 10000)

EC.3D(a, sigma = diag(1, 3), voxdim = c(1, 1, 1), num.vox = 10000)

}
\keyword{utilities}
