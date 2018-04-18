\name{Print Objects}
\alias{print.fclust}
\title{Print fclust related outputs}
\description{
  Print flcust related outputs.
}
\usage{
  \method{print}{fclust}(x, ...)
}
\arguments{
  \item{x}{an object with the class attributes.}
  \item{...}{other arguments to the \code{print} function.}
}
\details{
  \code{x} is the return result from \code{\link{fclust}()}.
}
\value{
  A summary of \code{\link{fclust}} object is printed.
}
\references{
  \url{http://maitra.public.iastate.edu/}
}
\author{
  Wei-Chen Chen and Ranjan Maitra.
}
\seealso{
  \code{\link{set.global}()}, \code{\link{fclust}()}.
}
\examples{
library(MixfMRI, quietly = TRUE)
set.seed(1234)

### Check 2d data.
da <- pval.2d.complex
id <- !is.na(da)
PV.gbd <- da[id]
# hist(PV.gbd, nclass = 100, main = "p-value")

### Test 2d data.
id.loc <- which(id, arr.ind = TRUE)
X.gbd <- t(t(id.loc) / dim(da))
ret <- fclust(X.gbd, PV.gbd, K = 2)
print(ret)
}
\keyword{programming}

