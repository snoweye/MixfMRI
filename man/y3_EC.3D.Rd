\name{EC.3D}
\alias{EC.3D}

\title{Expected Euler Characteristic for a 3D Random Field}
\description{Calculates the Expected Euler Characteristic for a 3D
  Random Field thesholded a level u.}
\usage{
EC.3D(u, sigma, voxdim = c(1, 1, 1), num.vox, type = c("Normal", "t"), df = NULL)
}
\arguments{
  \item{u}{The threshold for the field.}
  \item{sigma}{The spatial covariance matrix of the field.}
  \item{voxdim}{The dimensions of the cuboid 'voxels' upon which the
    discretized field is observed.}
  \item{num.vox}{The number of voxels that make up the field.}
  \item{type}{The marginal distribution of the Random Field (only Normal
  and t at present).}
  \item{df}{The degrees of freedom of the t field.}
}
\details{The Euler Characteristic \eqn{\chi_u} (Adler, 1981) is a
  topological measure that essentially counts the number of isolated
  regions of the random field above the threshold \eqn{u} minus the
  number of 'holes'. As \eqn{u} increases the holes disappear and
  \eqn{\chi_u} counts the number of local maxima. So when \eqn{u}
  becomes close to the maximum of the random field
  \eqn{Z_{\textrm{max}}} we have that
  \deqn{P( \textrm{reject}  H_0 | H_0  \textrm{true}) =
  P(Z_{\textrm{max}}) = P(\chi_u > 0) \approx E(\chi_u)}
  where \eqn{H_0} is the null hypothesis that there is no signicant
  positive actiavtion/signal present in the field. Thus the Type I error
  of the test can be controlled through knowledge of the Expected Euler characteristic.
  }
\value{The value of the expected Euler Characteristic.
}
\references{
  Adler, R. (1981) \emph{The Geometry of Random Fields.}. New York: Wiley.
  Worlsey, K. J. (1994) Local maxima and the expected euler
  characteristic of excursion sets of \eqn{\chi^2}, \eqn{f} and \eqn{t}
    fields. \emph{Advances in Applied Probability}, \bold{26}, 13-42.}
\author{J. L. Marchini}


\seealso{\code{\link{Threshold.RF}} }

\examples{
EC.3D(4.6, sigma = diag(1, 3), voxdim = c(1, 1, 1), num.vox = 10000)

EC.3D(4.6, sigma = diag(1, 3), voxdim = c(1, 1, 1), num.vox = 10000, type = "t", df = 100)

}
\keyword{utilities}
