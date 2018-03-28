### Get the Cov(logit eta_k) and 95% CE.

### Partial logL. Return a matrix with dimension M * N.
# partial.logL.I <- function(x, ETA, BETA, MU, SIGMA, post.z)

### Paritial all logit(eta_k). Return a K * (K - 1)
partial.logit.p <- function(ETA){
  K <- length(ETA)
  d.dETA <- 1 / ETA + 1 / (1 - ETA)
  ret <- diag(d.dETA[-K])
  ret <- rbind(ret, rep(-d.dETA[K], K - 1))
  ret
} # End of partial.logit.p().

cov.logit.ETA <- function(x, fcobj, cov.param = NULL){
  K <- fcobj$param$K
  ETA <- fcobj$param$ETA

  nabla.logit.ETA <- partial.logit.p(ETA)
  if(is.null(cov.param)){
    cov.param <- cov.param(x, fcobj, drop.ETA1 = FALSE)$cov
  }

  cov.logit.ETA <- nabla.logit.ETA %*% cov.param[1:(K-1), 1:(K-1)] %*%
                   t(nabla.logit.ETA)
  cov.logit.ETA
} # End of cov.logit.ETA().

