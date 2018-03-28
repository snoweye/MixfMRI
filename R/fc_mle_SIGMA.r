### This file contains major functions for solving MLE of SIGMA.

mle.SIGMA.I <- function(PARAM, i.k){
  N.gbd <- PARAM$N.gbd
  p.X <- PARAM$p.X

  B <- W.plus.y(.MixfMRIEnv$X.gbd, -PARAM$MU[, i.k], N.gbd, p.X) *
       sqrt(.MixfMRIEnv$Z.gbd[, i.k] / .MixfMRIEnv$Z.colSums[i.k])
  tmp.SIGMA <- colSums.gbd(B^2)
  tmp.SIGMA
} # End of mle.SIGMA.I().

mle.SIGMA.V <- function(PARAM, i.k){
  N.gbd <- PARAM$N.gbd
  p.X <- PARAM$p.X

  B <- W.plus.y(.MixfMRIEnv$X.gbd, -PARAM$MU[, i.k], N.gbd, p.X) *
       sqrt(.MixfMRIEnv$Z.gbd[, i.k] / .MixfMRIEnv$Z.colSums[i.k])
  tmp.SIGMA <- crossprod.gbd(B)
  tmp.SIGMA
} # End of mle.SIGMA.V().

