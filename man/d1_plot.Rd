\name{Plotting}
\alias{plotfclust}
\alias{plotfclustpv}
\alias{plotpv}
\alias{plotpvlegend}
\title{Main plotting function}
\description{
  Main plotting function in \pkg{MixfMRI}.
}
\usage{
  plotfclust(da, posterior, main = NULL, xlim = NULL, ylim = NULL)
  plotfclustpv(da, posterior, main = NULL, xlim = NULL, ylim = NULL)

  plotpv(da, posterior, PARAM, zlim = c(0, 0.01), plot.mean = TRUE,
         xlab = "", ylab = "", main = NULL, xlim = NULL, ylim = NULL,
         col = my.YlOrRd(), ignore.bg = FALSE)
  plotpvlegend(zlim = c(0, 0.01), n.level = 20, main = NULL,
               col = my.YlOrRd())
}
\arguments{
  \item{da}{a data set to be plotted.}
  \item{posterior}{a posterior data set to be plotted.}
  \item{PARAM}{a returning parameter object from \code{\link{fclust}()}.}
  \item{main}{title of the plot.}
  \item{xlim}{limits of x-axis.}
  \item{ylim}{limits of y-axis.}
  \item{zlim}{limits of z-axis.}
  \item{xlab}{labels of x-axis.}
  \item{ylab}{labels of y-axis.}
  \item{plot.mean}{if plotting mean values of each cluster.}
  \item{col}{colors to be drawn.}
  \item{ignore.bg}{if ignoring the background.}
  \item{n.level}{number of levels to be plotted.}
}
\details{
  These are example functions to plot results, simulations, and datasets.
}
\value{
  Return plots.
}
\references{
  \url{http://maitra.public.iastate.edu/}
}
\author{
  Wei-Chen Chen and Ranjan Maitra.
}
\seealso{
  \code{\link{set.global}()}.
}
\examples{
library(MixfMRI, quietly = TRUE)
set.seed(1234)
  
\donttest{
  ### Check 2d data.
  da <- pval.2d.complex
  id <- !is.na(da)
  PV.gbd <- da[id]
  hist(PV.gbd, nclass = 100, main = "p-value")
  
  ### Test 2d data.
  id.loc <- which(id, arr.ind = TRUE)
  X.gbd <- t(t(id.loc) / dim(da))
  ret <- fclust(X.gbd, PV.gbd, K = 3)
  print(ret)
  
  ### p-values of rest clusters.
  ret.lrt <- lrt(PV.gbd, ret$class, K = 3)
  print(ret.lrt)
  ret.lrt2 <- lrt2(PV.gbd, ret$class, K = 3)
  print(ret.lrt2)
  
  ### Plotting.
  par(mfrow = c(2, 2), mar = c(0, 0, 2, 0))
  plotpv(da, ret$posterior, ret$param,
         zlim = c(0.005, 0.008), main = "Mean of Beta Distribution")
  plotpv(da, ret$posterior, ret$param,
         plot.mean = FALSE, main = "p-value")
  par(mar = c(5.1, 4.1, 4.1, 2.1))
  plotpvlegend(zlim = c(0.005, 0.008), main = "Mean of Beta Distribution")
  plotpvlegend(zlim = c(0, 0.01), main = "p-value")
}
}
\keyword{programming}
