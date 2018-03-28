\name{False Discovery Rates for Spatial Signals}
\alias{fdr.bh.p1}
\alias{fdr.bh.p2}
\title{False Discovery Rates for Spatial Signals using Benjamini and Heller (2007)}
\description{
  Compute q-values Benjamini and Heller's (2007) approach for controlling FDR
  for spatial signals.
}
\usage{
  fdr.bh.p1(p, w = rep(1, length(p)), q = 0.05)
  fdr.bh.p2(p, w = rep(1, length(p)), q = 0.05)
}
\arguments{
  \item{p}{a p-value vector. No NA is allowed and all values are in [0, 1].}
  \item{w}{a weight vector for p-values.}
  \item{q}{a desired cutoff for adjusting p-values.}
}
\details{
  These functions implement first two procedures in Benjamini and Heller (2007)
  for controlling FDR for spatial signals.
}
\value{
  Return the number of rejected hypotheses and all
  corresponding q-values for the input p-values.
}
\references{
  \url{http://maitra.public.iastate.edu/}
}
\author{
  Wei-Chen Chen.
}
\seealso{
  \code{\link{qvalue}()}.
}
\examples{
library(MixfMRI, quietly = TRUE)
set.seed(1234)
da <- gendataset(phantom = shepp1fMRI, overlap = 0.01)
p <- da$pval[!is.na(da$pval)][1:100]
fdr.bh.p1(p)
fdr.bh.p2(p)
}
\keyword{programming}

