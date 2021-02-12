\name{MixfMRI-package}
\alias{MixfMRI-package}
\alias{MixfMRI}
\docType{package}
\title{fMRI Clustering Analysis}
\description{
  Utilizing model-based clustering (unsupervised)
  for fMRI data especially in a distributed manner. The methods
  includes 2D and 3D clustering analyses and segmentation analyses for
  fMRI signals where p-values are significant levels of active voxels
  which respond to stimulate of interesting. The analyses are
  mainly identifying active voxels/signals from normal brain behaviors.
  Workflows are also implemented utilizing high performance techniques.
}
\details{
\tabular{ll}{
  Package: \tab MixfMRI\cr
  Type: \tab Package\cr
  License: \tab GPL (>= 2)\cr
  LazyLoad: \tab yes\cr
}

  The main function of this package is \code{\link{fclust}()} that implements
  model-based clustering algorithm for fMRI signal data and provides
  unsupervised clustering results for the data. Several workflows implemented
  with high-performance computing techniques are also built in for automatically
  process clustering, hypothesis, cluster merging, and visualizations.
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
  \code{\link{fclust}()}, \code{\link{set.global}()}.
}
\examples{
library(MixfMRI, quietly = TRUE)
\donttest{
.rem <- function(){

  demo(fclust3d,'MixfMRI',ask=FALSE,echo=FALSE)
  demo(fclust2d,'MixfMRI',ask=FALSE,echo=FALSE)

}
}
}
\keyword{package}

