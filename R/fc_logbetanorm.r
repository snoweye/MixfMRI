### This file contains functions for log densities of beta and normal.
### These will majorly update .MixfMRIEnv$W.gbd.

logbetanorm.I <- function(PARAM, i.k){
  # log of beta density.
  ret <- dbeta(.MixfMRIEnv$PV.gbd,
               PARAM$BETA[[i.k]][1], PARAM$BETA[[i.k]][2], log = TRUE)

  # log of normal density.
  if(PARAM$p.X > 0){
    for(i.p in 1:PARAM$p.X){
      ret <- ret + dnorm(.MixfMRIEnv$X.gbd[, i.p],
                         mean = PARAM$MU[i.p, i.k],
                         sd = PARAM$SIGMA[[i.k]]$U[i.p], log = TRUE)
    }
  }

  # update global matrix.
  .MixfMRIEnv$W.gbd[, i.k] <- ret

  invisible()
} # End of logbetanorm.I().

logbetanorm.V <- function(PARAM, i.k){
  # log of beta density.
  ret <- dbeta(.MixfMRIEnv$PV.gbd,
               PARAM$BETA[[i.k]][1], PARAM$BETA[[i.k]][2], log = TRUE)

  # log of normal density.
  if(PARAM$p.X > 0){
    ret <- ret + dmvn.chol(.MixfMRIEnv$X.gbd, PARAM$MU[, i.k],
                           PARAM$SIGMA[[i.k]]$U, log = TRUE)
  }

  # update global matrix.
  .MixfMRIEnv$W.gbd[, i.k] <- ret

  invisible()
} # End of logbetanorm.V().

