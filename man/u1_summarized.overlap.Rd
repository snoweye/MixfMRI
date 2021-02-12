\name{Summarized Overlap}
\alias{summarized.overlap}
\title{Summarized Overlap}
\description{
  Compute summarized overlap on a given overlap (symmetric) matrix.
}
\usage{
  summarized.overlap(overlap.mat)
}
\arguments{
  \item{overlap.mat}{an overlap (symmetric) matrix.}
}
\details{
  \code{overlap.mat} is a \code{p * p} matrix containing pair wised overlaps
  of \code{p} experiments. \code{overlap.mat} is assumed a symmetric matrix.
  This function returns a summarized overlap based on the input
  \code{overlap.mat} that charactorizes the overall overlap behavior of the
  \code{p} experiments.
}
\value{
  A single value is returned.
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
\examples{
library(MixfMRI, quietly = TRUE)
set.seed(1234)
p <- 10  # 10 experiments.
overlap.mat <- diag(1, p)
overlap.mat[lower.tri(overlap.mat)] <- runif(p * (p - 1) / 2)
overlap.mat[upper.tri(overlap.mat)] <- t(overlap.mat)[upper.tri(overlap.mat)]
summarized.overlap(overlap.mat)
}
\keyword{programming}

