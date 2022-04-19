\name{Likelihood Mixture Tests with Identity Cov Matrix or Only p-values}
\alias{lmt.I}
\alias{lmt.pv}
\title{Likelihood Mixture Tests with Identity Cov Matrix or Only p-values}
\description{
  These functions test two mixture Gaussian fMRI models with diagonal
  covariance matrices and different numbers of clusters.
  These functions are similar to the \code{EMCluster::lmt(}), but is coded
  for fMRI models in \pkg{MixfMRI}.
}
\usage{
  lmt.I(fcobj.0, fcobj.a, X.gbd, PV.gbd, tau = 0.5, n.mc.E.delta = 1000,
        n.mc.E.chi2 = 1000, verbose = FALSE)
  lmt.pv(fcobj.0, fcobj.a, X.gbd, PV.gbd, tau = 0.5, n.mc.E.delta = 1000,
        n.mc.E.chi2 = 1000, verbose = FALSE)
}
\arguments{
  \item{fcobj.0}{a \code{fclust} object for the null hypothesis.}
  \item{fcobj.a}{a \code{fclust} object for the alternative hypothesis.}
  \item{X.gbd}{a data matrix of \code{N} voxel locations.
               \code{dim(X.gbd) = N x 3} for 3D data and
               \code{N x 2} for 2D data.}
  \item{PV.gbd}{a p-value vector of signals associated with voxels.
                \code{length(PV.gbd) = N}.}
  \item{tau}{proportion of null and alternative hypotheses.}
  \item{n.mc.E.delta}{number of Monte Carlo simulations for expectation
                      of delta (difference of logL).}
  \item{n.mc.E.chi2}{number of Monte Carlo simulations for expectation
                     of chisquare statistics.}
  \item{verbose}{if verbose.}
}
\details{
  This function calls several subroutines to compute information,
  likelihood ratio statistics, degrees of freedom, non-centrality
  of chi-squared distributions \dots etc. Based on Monte Carlo methods
  to estimate parameters of likelihood mixture tests, this function
  return a p-value for testing H0: \code{fcobj.0} v.s. Ha: \code{fcobj.a}.

  \code{lmt.pv()} only uses \code{PV.gbd}.
}
\value{
  A list of class \code{lmt.I} are returned.
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
  \code{EMCluster::lmt()}.
}
\examples{
library(MixfMRI, quietly = TRUE)
library(EMCluster, quietly = TRUE)
.FC.CT$model.X <- "I"
.FC.CT$check.X.unit <- FALSE
.FC.CT$CONTROL$debug <- 0
  
\donttest{
.rem <- function(){

  ### Fit toy1.
  set.seed(1234)
  X.gbd <- toy1$X.gbd
  PV.gbd <- toy1$PV.gbd
  ret.2 <- fclust(X.gbd, PV.gbd, K = 2)
  ret.3 <- fclust(X.gbd, PV.gbd, K = 3)
  ret.4 <- fclust(X.gbd, PV.gbd, K = 4)
  ret.5 <- fclust(X.gbd, PV.gbd, K = 5)
  
  ### ARI
  RRand(toy1$CLASS.gbd, ret.2$class)
  RRand(toy1$CLASS.gbd, ret.3$class)
  RRand(toy1$CLASS.gbd, ret.4$class)
  RRand(toy1$CLASS.gbd, ret.5$class)
  
  ### Test toy1.
  (lmt.23 <- lmt.I(ret.2, ret.3, X.gbd, PV.gbd))
  (lmt.24 <- lmt.I(ret.2, ret.4, X.gbd, PV.gbd))
  (lmt.25 <- lmt.I(ret.2, ret.5, X.gbd, PV.gbd))
  (lmt.34 <- lmt.I(ret.3, ret.4, X.gbd, PV.gbd))
  (lmt.35 <- lmt.I(ret.3, ret.5, X.gbd, PV.gbd))
  (lmt.45 <- lmt.I(ret.4, ret.5, X.gbd, PV.gbd))
  
  ### Test toy1 using p-values only.
  (lmt.pv.23 <- lmt.pv(ret.2, ret.3, X.gbd, PV.gbd))
  (lmt.pv.24 <- lmt.pv(ret.2, ret.4, X.gbd, PV.gbd))
  (lmt.pv.25 <- lmt.pv(ret.2, ret.5, X.gbd, PV.gbd))
  (lmt.pv.34 <- lmt.pv(ret.3, ret.4, X.gbd, PV.gbd))
  (lmt.pv.35 <- lmt.pv(ret.3, ret.5, X.gbd, PV.gbd))
  (lmt.pv.45 <- lmt.pv(ret.4, ret.5, X.gbd, PV.gbd))

}
}
}
\keyword{programming}
