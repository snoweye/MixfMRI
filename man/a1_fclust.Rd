\name{Main functions}
\alias{fclust}
\alias{set.global}
\alias{PARAM}
\alias{PARAM.org}
\title{Main MixfMRI function}
\description{
  Main MixfMRI functions.
}
\usage{
  fclust(X.gbd, PV.gbd, K = 2,
    PARAM.init = NULL,
    min.1st.prop = .FC.CT$INIT$min.1st.prop,
    max.PV = .FC.CT$INIT$max.PV,
    class.method = .FC.CT$INIT$class.method[1],
    RndEM.iter = .FC.CT$CONTROL$RndEM.iter,
    algorithm = .FC.CT$algorithm[1],
    model.X = .FC.CT$model.X[1],
    ignore.X = .FC.CT$ignore.X,
    stop.unstable = TRUE,
    MPI.gbd = .FC.CT$MPI.gbd, common.gbd = .FC.CT$common.gbd)

  set.global(X.gbd, PV.gbd, K = 2,
    min.1st.prop = .FC.CT$INIT$min.1st.prop,
    max.PV = .FC.CT$INIT$max.PV,
    class.method = .FC.CT$INIT$class.method[1],
    RndEM.iter = .FC.CT$CONTROL$RndEM.iter,
    algorithm = .FC.CT$algorithm[1],
    model.X = .FC.CT$model.X[1],
    ignore.X = .FC.CT$ignore.X,
    check.X.unit = .FC.CT$check.X.unit,
    MPI.gbd = .FC.CT$MPI.gbd, common.gbd = .FC.CT$common.gbd)
}
\arguments{
  \item{X.gbd}{a data matrix of \code{N} voxel locations.
               \code{dim(X.gbd) = N x 3} for 3D data and
               \code{N x 2} for 2D data.}
  \item{PV.gbd}{a p-value vector of signals associated with voxels.
                \code{length(PV.gbd) = N}.}
  \item{K}{number of clusters to be estimated.}
  \item{PARAM.init}{initial parameters.}
  \item{min.1st.prop}{lower bound of mixing proportion (ETA) of the 1st cluster
                      (uniform).}
  \item{max.PV}{upper bound of p-values where initializations pick from.}
  \item{class.method}{classification method for initializations.}
  \item{RndEM.iter}{number of RndEM iterations.}
  \item{algorithm}{either ``ecm'' (ECM), ``apecma'' (APECMa) or ``em'' (EM) algorithm.}
  \item{model.X}{either ``I'' or ``V'' for covariance matrix.}
  \item{ignore.X}{if \code{X.gbd} used in model, TRUE for PV.gbd only.}
  \item{check.X.unit}{if \code{X.gbd} are all in \code{[0, 1]}.}
  \item{stop.unstable}{if \code{fclust} stops if unstable results occur.}
  \item{MPI.gbd}{if MPI (``EGM'' algorithm) is used.}
  \item{common.gbd}{if \code{X.gbd} and \code{PV.gbd} are in common across all
                    ranks when \code{MPI.gbd = TRUE}.}
}
\details{
  The \code{fclust()} contains initialization and EM algorithms for clustering
  fMRI signal data which have two parts: \code{X.gbd} for voxel information
  either 2D or 3D, \code{PV.gbd} for p-value of signals associated with
  voxels. Each signal is assumed as a mixture distribution with \code{K}
  components with mixing proportion \code{ETA}, and each component has
  two independent coordinates with density functions: Beta and multivariate
  Normal distributions.

  Beta density:
  The 1st component is restricted by \code{min.1st.prop} and Beta(1, 1)
  distribution. The other \code{K - 1} components have Beta(alpha, beta)
  distribution with \code{alpha < 1 < beta}.

  Multivariate Normal density:
  \code{model.X = "I"} is for identity cov matrix of multivariate Normal
  distribution, and \code{"V"} for unstructured cov matrix.
  \code{ignore.X = TRUE} is to ignore \code{X.gbd} and normal density,
  i.e. only Beta density is used.

  Currently, APECMa and EM algorithms are implemented with EGM algorithm
  to speed up convergence if MPI is available.
  RndEM initialization is also implemented for better chance of good
  initial values for convergence.

  The \code{set.global()} has purposes: create a template/storage of
  parameters, save configurations, and called by \code{fclust()} to initial
  the parameters, such as \code{\link{initial.em.gbd}()} or
  \code{\link{initial.RndEM.gbd}()}.
}
\value{
  A list with class \code{fclust} by \code{fclust()} is returned
  which can be summarized by \code{\link{print.fclust}()}.

  A list \code{PARAM} or \code{PARAM.org} is returned by \code{set.global()}:
  \item{N.gbd}{number of observations (within the rank), and should be
               equal to \code{N.all} if \code{MPI.gbd = FALSE}.}
  \item{N.all}{numbers of observations (of all ranks
               if \code{MPI.gbd = TRUE}).}
  \item{N}{total number of observations (\code{sum(N.all)}).}
  \item{p}{dimension of an observation (3 for 2D signals, 4 for 3D signals),
           equivalent to total number of coordinates.}
  \item{p.X}{dimension of \code{X.gbd} (2 for 2D signals, 3 for 3D signals,
             0 when \code{ignore.X = TRUE}, number of voxel coordinates.}
  \item{K}{number of clusters.}
  \item{ETA}{mixing proportion, length \code{K}.}
  \item{log.ETA}{\code{log(ETA)}.}
  \item{BETA}{a list of length \code{K} containing parameters (alpha, beta)
              of Beta density.} 
  \item{MU}{a matrix of dimension \code{p.X by K}.}
  \item{SIGMA}{a list of length \code{K}, and each is of dimension
               \code{K x K}.}
  \item{logL}{log likelihood value.}
  \item{min.1st.prop}{carried from input.}
  \item{max.PV}{carried from input.}
  \item{class.method}{classification method of initializations.}
  \item{min.N.CLASS}{\code{p + 1}.}
  \item{model.X}{carried from input.}
}
\references{
  \url{http://maitra.public.iastate.edu/}
}
\author{
  Wei-Chen Chen and Ranjan Maitra.
}
\seealso{
  \code{\link{print.fclust}()}.
}
\examples{
library(MixfMRI, quietly = TRUE)
library(EMCluster, quietly = TRUE)
# .FC.CT$algorithm <- "em"
# .FC.CT$model.X <- "V"
# .FC.CT$ignore.X <- TRUE
.FC.CT$check.X.unit <- FALSE
set.seed(1234)

### Test toy1.
X.gbd <- toy1$X.gbd[, -3]
PV.gbd <- toy1$PV.gbd
PARAM <- fclust(X.gbd, PV.gbd, K = 2)
print(PARAM)
id.toy1 <- .MixfMRIEnv$CLASS.gbd
print(RRand(toy1$CLASS.gbd, id.toy1))

\donttest{
  ### Test toy2.
  X.gbd <- toy2$X.gbd[, -3]
  PV.gbd <- toy2$PV.gbd
  PARAM <- fclust(X.gbd, PV.gbd, K = 3)
  print(PARAM)
  id.toy2 <- .MixfMRIEnv$CLASS.gbd
  print(RRand(toy2$CLASS.gbd, id.toy2))
}
}
\keyword{programming}

