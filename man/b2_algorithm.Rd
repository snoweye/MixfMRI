\name{algorithm}
\alias{ecm.step.gbd}
\alias{apecma.step.gbd}
\alias{em.step.gbd}

\title{Main algorithms implemented in fclust}
\description{
  Main algorithms implemented in fclust.
}
\usage{
  ecm.step.gbd(PARAM.org)

  apecma.step.gbd(PARAM.org)

  em.step.gbd(PARAM.org)
}
\arguments{
  \item{PARAM.org}{an initialized \code{PARAM}, usually returned by
                   \code{\link{set.global}()}, \code{\link{initial.em.gbd}()},
                   and \code{\link{initial.RndEM.gbd}()}.}
}
\details{
  These are main algorithms implemented in \code{\link{fclust}()}.
}
\value{
  Return an optimized \code{PARAM}.
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
  \code{\link{set.global}()}, \code{\link{fclust}()}, \code{\link{PARAM}},
  \code{\link{PARAM.org}}.
}
\examples{
library(MixfMRI, quietly = TRUE)
library(EMCluster, quietly = TRUE)
# .FC.CT$algorithm <- "em"
# .FC.CT$model.X <- "V"
# .FC.CT$ignore.X <- TRUE
.FC.CT$check.X.unit <- FALSE

### Test toy1.
set.seed(1234)
X.gbd <- toy1$X.gbd
PV.gbd <- toy1$PV.gbd
PARAM <- set.global(X.gbd, PV.gbd, K = 2)
PARAM.new <- initial.em.gbd(PARAM)
PARAM.toy1 <- em.step.gbd(PARAM.new)
id.toy1 <- .MixfMRIEnv$CLASS.gbd
print(PARAM.toy1$ETA)
RRand(toy1$CLASS.gbd, id.toy1)

\donttest{
.rem <- function(){

  ### Test toy2.
  set.seed(1234)
  X.gbd <- toy2$X.gbd
  PV.gbd <- toy2$PV.gbd
  PARAM <- set.global(X.gbd, PV.gbd, K = 3)
  PARAM.new <- initial.em.gbd(PARAM)
  PARAM.toy2 <- em.step.gbd(PARAM.new)
  id.toy2 <- .MixfMRIEnv$CLASS.gbd
  print(PARAM.toy2$ETA)
  RRand(toy2$CLASS.gbd, id.toy2)

}
}
}
\keyword{programming}

