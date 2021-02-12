\name{Compute Q values}
\alias{qvalue}
\title{Q-values using Benjamini and Hochberg (1995)}
\description{
  Compute q-values Benjamini and Hochberg's (1995) approach for controlling FDR.
}
\usage{
  qvalue(p, method = c("BH1995", "BY2001"))
}
\arguments{
  \item{p}{a p-value vector.}
  \item{method}{using method by either BH1995 or BY2001}
}
\details{
  This function compute q-values using Benjamini and Hochberg's (1995)
  approach for controlling FDR. The function \code{bh.fdr} is originally
  written by Dr. Dan Nettleton.

  The Benjamini and Yeekutieli's (2001) approach for controlling FDR using
  the function \code{by.fdr} is coded by Wei-Chen Chen.
}
\value{
  Return corresponding q-values for the input p-values.
}
\references{
  Chen, W.-C. and Maitra, R. (2021)
  \dQuote{A Practical Model-based Segmentation Approach for Accurate
	  Activation Detection in Single-Subject functional Magnetic
	  Resonance Imaging Studies},
  \emph{arXiv:2102.03639}.
}
\author{
  Dan Nettleton.

  Modified by Wei-Chen Chen.
}
\seealso{
  \code{\link{dpval}()}, \code{\link{dmixpval}()}.
}
\examples{
library(MixfMRI, quietly = TRUE)
set.seed(1234)
da <- gendataset(phantom = shepp1fMRI, overlap = 0.01)
p <- da$pval[!is.na(da$pval)][1:100]
qvalue(p)
}
\keyword{programming}

