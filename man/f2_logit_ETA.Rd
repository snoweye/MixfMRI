\name{Covariance Matrices of Logit ETA}
\alias{cov.logit.ETA}
\title{Covariance Matrices of Logit ETA}
\description{
  These functions computes covariance matrix of logit ETA.
}
\usage{
  cov.logit.ETA(x, fcobj, cov.param = NULL)
}
\arguments{
  \item{x}{an input list of two elements \code{X.gbd} and \code{PV.gbd}.}
  \item{fcobj}{a \code{fclust} object.}
  \item{cov.param}{a covariance matrix of \code{dim = d * d} for parameters,
                   which is also a return of \code{cov.param()}.
                   \code{d} is total number of parameters which is dependent
                   on data and models.}
}
\details{
  These functions are required to compute covariance matrices of
  logit ETA.

  Input the returns of \code{cov.param()}
  to \code{cov.logit.ETA()} to obtain the cov matrix for logit ETA
  by the multivariate delta method on the cov matrix for parameters.
}
\value{
  A matrix.
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
  \code{EMCluster::lmt()}, \code{lmt.I()}.
}
\examples{
library(MixfMRI, quietly = TRUE)
.FC.CT$model.X <- "I"
.FC.CT$CONTROL$debug <- 0
K <- 3
  
\donttest{
.rem <- function(){

  ### Fit toy1.
  set.seed(1234)
  X.gbd <- toy1$X.gbd
  X.range <- apply(X.gbd, 2, range)
  X.gbd <- t((t(X.gbd) - X.range[1,]) / (X.range[2,] - X.range[1,]))
  PV.gbd <- toy1$PV.gbd
  fcobj <- fclust(X.gbd, PV.gbd, K = K, min.1st.prop = 0.5)
  
  ### Test cov matrix of posterior z.
  x <- list(X.gbd = X.gbd, PV.gbd = PV.gbd)
  post.z <- post.prob(x, fcobj)
  cov.param <- cov.param(x, fcobj, post.z)
  cov.logit.ETA <- cov.logit.ETA(x, fcobj, cov.param = cov.param$cov)
  
  ### Compute cov matrxi of eta_k - eta_1 for all k > 1.
  A <- cbind(rep(-1, K - 1), diag(1, K - 1))
  ETA <- fcobj$param$ETA
  log.or <- log(ETA / (1 - ETA)) \%*\% t(A)
  cov.log.or <- A \%*\% cov.logit.ETA \%*\% t(A)

}
}
}
\keyword{programming}
