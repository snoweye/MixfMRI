\name{Covariance Matrices}
\alias{post.prob}
\alias{cov.param}
\alias{cov.post.z}
\alias{cov.logit.z}
\title{Covariance Matrices}
\description{
  These functions compute posterior probabilities,
  Fisher information with covariance matrix of parameters,
  covariance matrix of posterior probabilities, and covariance matrix of
  logit posterior probabilities.
}
\usage{
  post.prob(x, fcobj)
  cov.param(x, fcobj, post.z, drop.ETA1 = FALSE)
  cov.post.z(x, fcobj, post.z, cov.param = NULL, all.x = FALSE,
                 drop.ETA1 = FALSE)
  cov.logit.z(x, fcobj, post.z, cov.param = NULL, cov.post.z = NULL,
                  all.x = FALSE, drop.ETA1 = FALSE)
}
\arguments{
  \item{x}{an input list of two elements \code{X.gbd} and \code{PV.gbd}.}
  \item{fcobj}{a \code{fclust} object.}
  \item{post.z}{a matrix of \code{dim = N * K} for posterior probabilities,
                which is also the return value of \code{post.prob()}.}
  \item{cov.param}{a covariance matrix of \code{dim = d * d} for parameters,
                   which is also a return of \code{cov.param()}.
                   \code{d} is total number of parameters which is dependent
                   on data and models.}
  \item{cov.post.z}{a covariance list of length equal to number of active
                    voxels, which is also a return of \code{cov.post.z()}.}
  \item{all.x}{all cov matrices for all observations are returned if TRUE,
               while for only active observations (those of class ids are
               greater than 1) if FALSE.}
  \item{drop.ETA1}{if drop the \code{ETA[1]} from the cov matrix due to the
                   \code{min.1st.prop} constrain.}
}
\details{
  These functions are required to compute covariance matrices of
  parameters and posterior probabilities.

  Use \code{post.prob()} to get the posterior probabilities.

  Input the returns of \code{post.prob()}
  to \code{cov.param()} to obtain the cov matrix for parameters
  (inversed Fisher information obtained from inner product of gradient
   of log observed data likelihood). A list is returned with \code{I}
  for Fisher information, and \code{cov} for the covariance matrix which
  is inverted by \code{ginv()}.

  Input the returns of \code{post.prob()} and \code{cov.param()}
  to \code{cov.post.z()} to obtain the cov matrix for posterior
  probabilities by the multivariate delta method on the cov matrix for
  parameters.

  Input the returns of \code{post.prob()}, \code{cov.param()}, and
  \code{cov.post.z()} to \code{cov.logit.z()} to obtain cov matrix
  for logit posterior probabilities by the multivariate delta method on
  cov matrix of posterior probabilities.
}
\value{
  A matrix or a list is returned.

  The \code{cov.param()} will return a list containing two elements
  \code{I} for the Fisher information, and \code{cov} for the covariance matrix
  by generalized inversed of the Fisher information. The dimension of both
  elements are \code{d * d} where \code{d = K * 7 - 4} for 2D data and
  \code{d = K * 9 - 4} for 3D data if \code{drop.ETA1 = TRUE}, otherwise
  they are \code{d = K * 7 - 3} and \code{d = K * 9 -4}, respectively.

  The \code{cov.post.z()} will return a list containing cov matrices of
  posterior probabilities for each valid/selected voxel.

  The \code{cov.logit.z()} will return a list containing cov matrices of
  logit posterior probabilities for each valid/selected voxel.
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
library(EMCluster, quietly = TRUE)
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
  
  ### Test cov matrix of posterior z and logit posterior z.
  x <- list(X.gbd = X.gbd, PV.gbd = PV.gbd)
  post.z <- post.prob(x, fcobj)
  cov.param <- cov.param(x, fcobj, post.z = post.z)
  cov.post.z <- cov.post.z(x, fcobj, post.z = post.z,
                               cov.param = cov.param$cov)
  cov.logit.z <- cov.logit.z(x, fcobj, post.z = post.z,
                                 cov.param = cov.param$cov,
                                 cov.post.z = cov.post.z)
  
  ### Compute cov matrix of log odds ratio for all k > 1.
  A <- cbind(rep(-1, K - 1), diag(1, K - 1))
  logit.p <- log(post.z[fcobj$class != 1,] / (1 - post.z[fcobj$class != 1,]))
  log.or <- logit.p \%*\% t(A)
  cov.log.or <- lapply(cov.logit.z, function(x) A \%*\% x \%*\% t(A))
  
  ### Check if 0 vector covered by 95\% confidence ellipsoid.
  id <- 1
  plot(log.or[id,],
       xlim = log.or[id, 1] + c(-5, 5),
       ylim = log.or[id, 2] + c(-5, 5),
       main = "1st observation", xlab = "x", ylab = "y")
  plotBN(log.or[id,], cov.log.or[[id]])
  points(0, 0, col = 2)

}
}
}
\keyword{programming}
