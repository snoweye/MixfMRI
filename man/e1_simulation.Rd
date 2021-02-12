\name{Simulations}
\alias{gendataset}
\title{Generate datasets for MixfMRI simulations}
\description{
  Generate datasets for MixfMRI simulations
}
\usage{
  gendataset(phantom, overlap, smooth = FALSE)
}
\arguments{
  \item{phantom}{a phantom dataset.}
  \item{overlap}{a desired overlap level.}
  \item{smooth}{if \code{gcv.smooth2d()} be applied to the data.}
}
\details{
  This is a function to generate simulated fMRI data based on the input
  \code{phantom} and the desired \code{overlap} level for the fMRI p-value.
}
\value{
  Return a list contains \code{eta} for mixing proportion,
  \code{overlap} for the desired level, \code{mu} for center of p-values,
  \code{class.id} for the true classifications where p-values belong to,
  \code{tval} for the testing statistics, and \code{pval} for the p-values
  of interesting in simulations.
}
\references{
  Chen, W.-C. and Maitra, R. (2021)
  \dQuote{A Practical Model-based Segmentation Approach for Accurate
	  Activation Detection in Single-Subject functional Magnetic
	  Resonance Imaging Studies},
  \emph{arXiv:2102.03639}.
}
\author{
  Wei-Chen Chen and Ranjan Maitra.
}
\seealso{
  \code{\link{set.global}()}.
}
\examples{
\donttest{
.rem <- function(){

  library(MixfMRI, quietly = TRUE)
  set.seed(1234)
  da <- gendataset(phantom = shepp1fMRI, overlap = 0.01)$pval
  da2 <- gendataset(phantom = shepp2fMRI, overlap = 0.01)$pval

  par(mfrow = c(2, 2), mar = rep(0.05, 4))
  image(shepp1fMRI[50:210, 50:210], axes = FALSE)
  image(shepp2fMRI[50:210, 50:210], axes = FALSE)
  image(da[50:210, 50:210], axes = FALSE)
  image(da2[50:210, 50:210], axes = FALSE)

}
}
}
\keyword{programming}

