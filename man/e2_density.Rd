\name{Density function of p-values}
\alias{dpval}
\alias{dmixpval}
\title{Density function of p-values}
\description{
  These functions based on normal assumption and transformation to
  derive a (mixture) density function of p-values.
}
\usage{
  dpval(x, mu = 0, log = FALSE)
  dmixpval(x, eta, mu)
}
\arguments{
  \item{x}{support of p-values which should be between 0 and 1.}
  \item{mu}{hypothetical mean of testing statistics (in normal distribution)
            for producing p-values.}
  \item{log}{if return log of density.}
  \item{eta}{mixing proportion of \code{K} components if a mixture is assumed.}
}
\details{

  Note that \code{eta} and \code{mu} in \code{dmixpval()} are of length
  \code{K} for \code{K} component mixtures.
}
\value{
  Corresponding density values (to the input \code{x}) are returned.
}
\references{
  \url{http://maitra.public.iastate.edu/}
}
\author{
  Wei-Chen Chen and Ranjan Maitra.
}
\seealso{
  \code{\link{gendataset}()}, \code{\link{qvalue}()}.
}
\examples{
library(MixfMRI, quietly = TRUE)
set.seed(1234)
da <- gendataset(phantom = shepp1fMRI, overlap = 0.01)
x <- da$pval[!is.na(da$pval)][1:100]
dpval(x)
dmixpval(x, mu = da$mu, eta = da$eta)
}
\keyword{programming}

