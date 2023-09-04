\name{Threshold.Bonferroni}
\alias{Threshold.Bonferroni}
%- Also NEED an `\alias' for EACH other topic documented here.
\title{Calculates Bonferroni Threshold}
\description{Calculate the Bonferroni threshold for n iid tests that
  results in an overall p-value of p.val. The tests can be distributed
  as Normal, t or F.
}
\usage{
Threshold.Bonferroni(p.val, n, type = c("Normal", "t", "F"), df1 = NULL, df2 = NULL)
}
\arguments{
  \item{p.val}{The required overall p-value.}
  \item{n}{The number of tests.}
  \item{type}{The distribution of the tests. One of "Normal", "t" or "F"}
  \item{df1}{The degrees of freedom of the t-distribution or the first
    degrees of freedom parameter for the F distribution.}
  \item{df2}{The second degrees of freedom parameter for the F distribution.}
}
\details{
  Note: This function is directly copied from "AnalyzeFMRI".
}
\value{Returns the Bonferroni threshold.
}
\author{Pierre Lafaye De Micheaux and J. L. Marchini.}

\examples{
Threshold.Bonferroni(0.05, 1000)

Threshold.Bonferroni(0.05, 1000, type = c("t"), df1 = 20)

Threshold.Bonferroni(0.05, 1000, type = c("F"), df1 = 3, df2 = 100)

}
\keyword{utilities}
