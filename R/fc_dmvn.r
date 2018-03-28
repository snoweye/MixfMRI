### Density of mulitvariate normal distribution.
# X: N.gbd by p.X
# MU: length p.X
# full: p.X by p.X

dmvn <- function(X, MU, full, log = FALSE){
  N <- nrow(X)
  p <- ncol(X)

  tmp.S <- full
  logdet <- log(abs(det(tmp.S)))

  B <- W.plus.y(X, -MU, N, p)
  tmp.S <- solve(tmp.S)
  distval <- rowSums((B %*% tmp.S) * B)

  ret <- -(.MixfMRIEnv$p.times.logtwopi + logdet + distval) * 0.5
  if(!log){
    ret <- exp(ret)
  }
  ret
} # End of dmvn().


# SIGMA: list(full, U, U.check)
dmvn.chol <- function(X, MU, U, log = FALSE){
  N <- nrow(X)
  p <- ncol(X)

  logdet <- sum(log(abs(diag(U)))) * 2

  B <- W.plus.y(X, -MU, N, p)
  B <- B %*% backsolve(U, diag(1, p))
  distval <- rowSums(B * B)

  ret <- -(.MixfMRIEnv$p.times.logtwopi + logdet + distval) * 0.5
  if(!log){
    ret <- exp(ret)
  }
  ret
} # End of dmvn.chol().

