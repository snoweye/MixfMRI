\name{initial}
\alias{initial.em.gbd}
\alias{initial.RndEM.gbd}
\title{Main initialization functions}
\description{
  Main initialization functions.
}
\usage{
  initial.em.gbd(PARAM)

  initial.RndEM.gbd(PARAM)
}
\arguments{
  \item{PARAM}{a list of uninitialized parameters, as usual, the returned
               values of \code{\link{set.global}()}, to be initialized
               according to data (inside \code{PARAM}).}
}
\details{
  \code{initial.em.gbd()} takes in a template of \code{PARAM} (uninitialized),
  and usually is available by calling \code{\link{set.global}()}, then
  return an initialized \code{PARAM} which is ready for EM runs.

  Internally, there are six different initializations implemented for
  the function \code{initial.em.gbd()} including \code{prob.extend},
  \code{prob.simple}, \code{qnorm.extend}, \code{qnorm.simple}, \code{extend},
  and \code{simple}. These methods are mainly based on transformation of
  original space of data (p-values and voxel locations) into more linear
  space such that the Euclidean distance more makes sense (fairly) to
  classify data in groups.

  \code{initial.RndEM.gbd()} implement RndEM initialization algorithm
  based on repeated calling \code{initial.em.gbd()}.

  Note that all configurations are included in \code{PARAM} set by
  \code{\link{set.global}()}.
}
\value{
  These functions return an initialized \code{PARAM} for EM runs based on
  pre-stored configuration within the input uninitialized \code{PARAM}.
}
\references{
  \url{http://maitra.public.iastate.edu/}
}
\author{
  Wei-Chen Chen and Ranjan Maitra.
}
\seealso{
  \code{\link{set.global}()}, \code{\link{fclust}()}, \code{\link{PARAM}}.
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
\keyword{programming}

