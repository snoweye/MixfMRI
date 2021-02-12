\name{MixfMRI Control}
\alias{.FC.CT}
\docType{data}
\title{Sets of controls in MixfMRI}
\description{
  These sets of controls are used to provide default values in this package.
}
\format{
  Objects contain several parameters for methods.
}
\details{
  The elements of \code{.FC.CT} are default values for main controls of
  \pkg{MixfMRI} including
\tabular{lcl}{
  Elements             \tab Default   \tab Usage \cr
  \code{algorithm}     \tab "apecma"  \tab implemented algorithm \cr
  \code{optim.method}  \tab "BFGS"    \tab optimization method \cr
  \code{model.X}       \tab "I"       \tab cov matrix structure \cr
  \code{ignore.X}      \tab FALSE     \tab if using voxel information \cr
  \code{check.X.unit}  \tab TRUE      \tab if checking X in \code{[0, 1]} \cr
  \code{CONTROL}       \tab a list    \tab see CONTROL next for details \cr
  \code{INIT}          \tab a list    \tab see INIT next for details \cr
  \code{LRT}           \tab a list    \tab see LRT next for details \cr
  \code{MPI.gbd}       \tab FALSE     \tab if MPI speedup available \cr
  \code{common.gbd}    \tab TRUE      \tab if X in common gbd format \cr
}

  The elements of \code{CONTROL} are default values for optimization controls
  of implemented EM algorithm including
\tabular{lcl}{
  Elements            \tab Default    \tab Usage \cr
  \code{max.iter}     \tab 1000       \tab maximum number of EM iterations \cr
  \code{abs.err}      \tab 1e-4       \tab absolute error of convergence \cr
  \code{rel.err}      \tab 1e-6       \tab relative error of convergence \cr
  \code{debug}        \tab 1          \tab debugging level \cr
  \code{RndEM.iter}   \tab 10         \tab RndEM iterations \cr
  \code{exp.min}      \tab log(.Machine$double.xmin)         \tab minimum exponential power \cr
  \code{exp.max}      \tab log(.Machine$double.xmax)         \tab maximum exponential power \cr
  \code{sigma.ill}    \tab 1e-6       \tab ill condition limit \cr
  \code{DS.max}       \tab 1e+4       \tab maximum \code{chol()} cov matrix \cr
  \code{DS.min}       \tab 1e-6       \tab minimum \code{chol()} cov matrix \cr
}


  The elements of \code{INIT} are default values or limitations for initial
  parameters implemented for EM algorithm including
\tabular{lcl}{
  Elements              \tab Default       \tab Usage \cr
  \code{min.1st.prop}   \tab 0.8           \tab minimum proportion of 1st cluster \cr
  \code{max.PV}         \tab 0.1           \tab maximum p-value for initialization \cr
  \code{BETA.alpha.min} \tab 0 + 1e-6      \tab minimum value of alpha parameter of Beta distribution \cr
  \code{BETA.alpha.max} \tab 1 - 1e-6      \tab maximum value of alpha parameter of Beta distribution \cr
  \code{BETA.beta.min}  \tab 1 + 1e-6      \tab minimum value of beta parameter of Beta distribution \cr
  \code{BETA.beta.max}  \tab 1e+6          \tab maximum value of beta parameter of Beta distribution \cr
  \code{max.try.iter}   \tab 10            \tab maximum retry iterations if result is unstable \cr
  \code{class.method}   \tab "prob.extned" \tab classification method at initializations \cr
}

  The elements of \code{LRT} are default values or limitations for likelihood
  ratio tests including
\tabular{lcl}{
  Elements              \tab Default   \tab Usage \cr
  \code{H0.alpha}       \tab 1         \tab null hypothesis alpha parameter of Beta distribution \cr
  \code{H0.beta}        \tab 1         \tab null hypothesis beta parameter of Beta distribution \cr
  \code{H0.mean}        \tab 0.05      \tab null hypothesis mean of Beta distribution \cr
}
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
  \code{\link{set.global}()}, \code{\link{fclust}()}.
}
\keyword{global variables}
