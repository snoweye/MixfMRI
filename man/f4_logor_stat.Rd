\name{Compute Statistics for Log Odds Ratio of Posterior Probability}
\alias{logor.stat}
\title{Compute Statistics for Log Odds Ratio of Posterior Probability}
\description{
  The function computes statistics for log odds ratio of posterior probability.
}
\usage{
  logor.stat(x, fcobj, post.z, cov.param = NULL, cov.post.z = NULL,
                 cov.logit.z = NULL, all.x = FALSE, drop.ETA1 = FALSE)
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
  \item{cov.logit.z}{a covariance list of length equal to number of active
                     voxels, which is also a return of \code{cov.logit.z()}.}
  \item{all.x}{all cov matrices for all observations are returned if TRUE,
               while for only active observations (those of class ids are
               greater than 1) if FALSE.}
  \item{drop.ETA1}{if drop the \code{ETA[1]} from the cov matrix due to the
                   \code{min.1st.prop} constrain.}
}
\details{
  For posterior probability,
  this function compute log odd ratio, cov matrix of log odd ratio,
  degrees of freedom, and testing statistics.
}
\value{
  A list is returned with four elements: \code{log.or}, \code{cov.log.or},
  \code{df}, and \code{test.stat}.
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
  \code{post.prob()}, \code{cov.param()}, \code{cov.post.z()},
  \code{cov.logit.z()}.
}
\examples{
library(MixfMRI, quietly = TRUE)
.FC.CT$model.X <- "I"
.FC.CT$CONTROL$debug <- 0
K <- 3

### Fit toy1.
set.seed(1234)
X.gbd <- toy1$X.gbd
X.range <- apply(X.gbd, 2, range)
X.gbd <- t((t(X.gbd) - X.range[1,]) / (X.range[2,] - X.range[1,]))
PV.gbd <- toy1$PV.gbd
fcobj <- fclust(X.gbd, PV.gbd, K = K, min.1st.prop = 0.5)

### Test log odds ratio.
x <- list(X.gbd = X.gbd, PV.gbd = PV.gbd)
post.z <- post.prob(x, fcobj)
lor <- logor.stat(x, fcobj, post.z)

### Check if 95\% CE covers log odd ratio = 1.
id <- !is.na(lor$df)
id.cover.0 <- which(lor$test.stat[id] < pchisq(0.95, lor$df[id]))

### Get voxels needed for merging.
id.active <- which(fcobj$class != 1)
id.merge <- id.active[id][id.cover.0]

### Check results.
post.z[id.merge,]
cbind(toy1$X.gbd[id.merge,], toy1$PV.gbd[id.merge])
}
\keyword{programming}
