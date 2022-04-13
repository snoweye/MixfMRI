\name{LRT}
\alias{lrt}
\alias{lrt2}
\alias{lrt.betamean}
\alias{lrt.betaab}
\title{Likelihood ratio tests}
\description{
  Likelihood ratio tests for merging clusters.
}
\usage{
  lrt(PV.gbd, CLASS.gbd, K, H0.alpha = .FC.CT$LRT$H0.alpha,
      H0.beta = .FC.CT$LRT$H0.beta)

  lrt2(PV.gbd, CLASS.gbd, K, H0.mean = .FC.CT$LRT$H0.mean,
    upper.beta = .FC.CT$INIT$BETA.beta.max, proc = c("1", "2", "weight"))

  lrt.betamean(PV.gbd, CLASS.gbd, K, proc = c("1", "2"))

  lrt.betaab(PV.gbd, CLASS.gbd, K, proc = c("1", "2"))
}
\arguments{
  \item{PV.gbd}{a p-value vector of signals associated with voxels.
                \code{length(PV.gbd) = N}.}
  \item{CLASS.gbd}{a classification vector of signals associated with voxels.
                   \code{length(CLASS.gbd) = N}.}
  \item{K}{number of clusters.}
  \item{H0.alpha}{null hypothesis for the alpha parameter of Beta distribution.}
  \item{H0.beta}{null hypothesis for the beta parameter of Beta distribution.}
  \item{H0.mean}{null hypothesis for the mean of Beta distribution.}
  \item{upper.beta}{BETA.beta.max, maximum value of beta parameter of
                    Beta distribution.}
  \item{proc}{q-value procedure for adjusting p-values.}
}
\details{
  These functions perform likelihood ratio tests for merging clusters.
  Only p-values coordinates (Beta density) are tested, while
  voxel location coordinates (multivariate Normal density) are not involved
  in testing.

  \code{lrt.betamean} tests if means of any two pairs of mixture
  (p-value) component were the same.
  The chi-square distribution with 1 degree of freedom is used.

  \code{lrt.betaab} tests if alpha and beta of any two pairs of mixture
  (p-value) components were the same.
  The chi-square distribution with 2 degrees of freedom is used.

  Procedure to adjust/select plausible p-values,
  \code{proc = "1"} uses q-value \code{qvalue()},
  \code{proc = "2"} uses \code{fdr.bh.p2()}, and
  \code{proc = "weight"} uses a weighted version of \code{fdr.bh.p2()}.
}
\value{
  A matrix contains MLEs of parameters of Beta distribution under the null
  hypothesis and the union of null and alternative hypotheses. The matrix
  also contains testing statistics and p-values.
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
  \code{\link{PARAM}}.
}
\examples{
library(MixfMRI, quietly = TRUE)
set.seed(1234)

### Test 2d data.
da <- pval.2d.mag
id <- !is.na(da)
PV.gbd <- da[id]
id.loc <- which(id, arr.ind = TRUE)
X.gbd <- t(t(id.loc) / dim(da))
ret <- fclust(X.gbd, PV.gbd, K = 2, min.1st.prop = 0.95)
# print(ret)

### p-values of rest clusters.
ret.lrt <- lrt(PV.gbd, ret$class, K = 2)
print(ret.lrt)
\donttest{
.rem <- function(){

  ret.lrt2 <- lrt2(PV.gbd, ret$class, K = 3)
  print(ret.lrt2)

}
}
}
\keyword{programming}

