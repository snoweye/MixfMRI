\name{Threshold.FDR}
\alias{Threshold.FDR}
\title{False Discovery Rate (FDR) Threshold}
\description{Calculates the False Discovery Rate (FDR) threshold for a
  given vector of statistic values.
}
\usage{
Threshold.FDR(x, q, cV.type = 2, type = c("Normal", "t", "F"), df1 = NULL, df2 = NULL)
}
\arguments{
  \item{x}{A vector of test statistic values.}
  \item{q}{The desired False Discovery Rate threshold.}
  \item{cV.type}{A flag that specifies the assumptions about the joint
    distribution of p-values. Choose cV.type = 2 for fMRI data (see
    Genovese et al (2001)}
  \item{type}{The distribution of the statistic values. Either "Normal",
    "t" or "F".}
  \item{df1}{The degrees of freedom of the t-distribution or the first
    degrees of freedom parameter for the F distribution.}
  \item{df2}{The second degrees of freedom parameter for the F distribution.}
}

\value{Returns the FDR threshold.
}
\references{Genovese et al. (2001) Thresholding of Statistical Maps in
  Functional NeuroImaging Using the False Discovery Rate.}
\author{J. L. Marchini}

\examples{
x <- c(rnorm(1000), rnorm(100, mean = 3))
Threshold.FDR(x = x, q = 0.20, cV.type = 2) 
}
\keyword{utilities}
