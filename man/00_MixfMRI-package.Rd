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
  The main function of this package is \code{\link{fclust}()} that implements
  model-based clustering algorithm for fMRI signal data and provides
  unsupervised clustering results for the data. Several workflows implemented
  with high-performance computing techniques are also built in for automatically
  process clustering, hypothesis, cluster merging, and visualizations.
}
\references{
  Chen, W.-C. and Maitra, R. (2023)
  \dQuote{A practical model-based segmentation approach for improved
          activation detection in single-subject functional magnetic
          resonance imaging studies},
  \emph{Human Brain Mapping}, \bold{44}(16), 5309--5335.
  (\emph{doi:10.1002/hbm.26425})
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

