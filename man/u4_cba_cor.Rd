\name{Generalized Cluster-Based Analysis (CBA) Method}
\alias{cba.cor}
\alias{cba.cor.2d}
\alias{cba.cor.3d}
\title{Generalized Cluster-Based Analysis (CBA) Method}
\description{
  Find clusters in 2D or 3D based on a generalized CBA method. The CBA
  method is originally proposed by Heller, et.al. (2006) using the correlation
  of two time series as the similarity of two spatial locations.
}
\usage{
  cba.cor(da.ts, da.m = NULL, adj.dist = TRUE, fun.sim = stats::cor)
  cba.cor.2d(da.ts, da.m = NULL, adj.dist = TRUE, fun.sim = stats::cor)
  cba.cor.3d(da.ts, da.m = NULL, adj.dist = TRUE, fun.sim = stats::cor)
}
\arguments{
  \item{da.ts}{a time series array of dimensions \code{x * y * z * t}.}
  \item{da.m}{a mask determining inside of brain or not.}
  \item{adj.dist}{if adjust correlations by distance.}
  \item{fun.sim}{a function computing simility of two locations.}
}
\details{
  These functions implement the 2D and 3D versions of CBA proposed by
  Heller, et.al. (2006).

  \code{da.ts} should have dimensions \code{x * y * z * t} for 3D data
  and \code{x * y * time} for 2D data. Similarly, \code{da.m} would have
  \code{x * y * z} and \code{x * y} correspondingly.

  \code{da.m} has values 0 or 1 indicating outside or inside a brain,
  respectively.

  \code{fun.sim(a, B)} is a function return similarity between a location
  \code{a} and N neighboring locations \code{B} where \code{a} is of dimension
  \code{t * 1} and \code{B} is of dimension \code{t * N}.
  Ideally, \code{fun.sim}()
  should return values of similarity which take values
  between 0 and 1 where 0 means totally different and
  1 means completely identical of two spatial locations.
  By default, \code{stats::cor} is used.
  See the example section next for user defined functions for
  \code{fun.sim}().
}
\value{
  Return the cluster ids for each voxel. NA for outside of brain if
  \code{da.m} is provided.
}
\references{
  Heller, et.al. (2006) ``Cluster-based analysis of FMRI data'',
  NeuroImage, 33(2), 599-608.

  Chen, W.-C. and Maitra, R. (2021)
  \dQuote{A Practical Model-based Segmentation Approach for Accurate
	  Activation Detection in Single-Subject functional Magnetic
	  Resonance Imaging Studies},
  \emph{arXiv:2102.03639}.
}
\author{
  Wei-Chen Chen.
}
\seealso{
  \code{\link{fdr.bh.p1}()}, \code{\link{fdr.bh.p2}()}.
}
\examples{
### Simulated data
library(MixfMRI, quietly = TRUE)
dim <- c(4, 5, 4, 10)
set.seed(123)
da.ts <- array(rnorm(prod(dim)), dim = dim)
id.class <- suppressWarnings(cba.cor(da.ts))
table(id.class)

fun.tanh <- function(a, B){
  d <- 1 / apply(B, 2, function(b){ dist(rbind(as.vector(a), b)) })
  tanh(d)
}
id.class.tanh <- suppressWarnings(cba.cor(da.ts, fun.sim = fun.tanh))
table(id.class.tanh)

fun.logit <- function(a, B){
  d <- dist(t(cbind(a, B)))[1:ncol(B)]
  (1 / (1 + exp(-d))) * 2 - 1
}
id.class.logit <- suppressWarnings(cba.cor(da.ts, fun.sim = fun.logit))
table(id.class.logit)

\donttest{
.rem <- function(){

  ### Real data
  # library(AnalyzeFMRI, quietly = TRUE)
  # library(oro.nifti, quietly = TRUE)

  # fn <- "pb02_volreg_tlrc.nii"
  # da <- readNIfTI(fn)
  # da.ts <- da@.Data

  # fn <- "mask_anat.nii"
  # da <- readNIfTI(fn)
  # da.m <- da@.Data

  # id.class <- suppressWarnings(cba.cor(da.ts, da.m))
  # dim(id.class) <- dim(da.m)
  # length(table(id.class))

}
}
}
\keyword{programming}

