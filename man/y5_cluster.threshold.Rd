\name{cluster.threshold}
\alias{cluster.threshold}
\alias{cluster_mass}

\title{Cluster threshold an array.}
\description{Calculate contiguous clusters of locations in a 3D array that are
  above some threshold and with some minimum size.
}
\usage{
cluster.threshold(x, nmat = NULL, level.thr = 0.5, size.thr)
}
\arguments{
  \item{x}{A 3D array}
  \item{nmat}{A matrix with 3 columns specifying the neighbourhood
    system. Default is 6 nearest neighbours in 3D.}
  \item{level.thr}{The level at which to threshold the array
    values. Default is 0.5 and is designed to cluster 0-1 arrays.}
  \item{size.thr}{The cluster size threshold.}
}
\details{
  Note: This function is directly copied from "AnalyzeFMRI".
}
\value{Returns an array of the same size as x with a 1 at all locations
  which have a value above level.thr and are in a cluster of similiar
  locations with size greater than size.thr. 
}
\author{J. L. Marchini}

\examples{

x <- array(0, dim = c(64, 64, 21))
x[10:20, 10:20, 1:5] <- 1
x[30:40, 30:40, 6:7] <- 1
x[50, 50, 8:9] <- 1

a <- cluster.threshold(x, size.thr = 400)
sum(x) ## should be 849
sum(a) ## should be 605

}
\keyword{utilities}
