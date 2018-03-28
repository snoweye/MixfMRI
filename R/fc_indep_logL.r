### This is an independent function to provide logL.
### Note that Inf, -Inf, NA, NaN is drop from the summation.

indep.logL.I <- function(PARAM){
  nrow <- nrow(.MixfMRIEnv$X.gbd)
  p.X <- PARAM$p.X

  ret <- matrix(0, nrow = nrow, ncol = PARAM$K)
  for(i.k in 1:PARAM$K){
    log.a <- dbeta(.MixfMRIEnv$PV.gbd, PARAM$BETA[[i.k]][1],
                   PARAM$BETA[[i.k]][2], log = TRUE)
    log.b <- 0.0
    if(p.X > 0){
      for(i.p in 1:p.X){
        log.b <- log.b + dnorm(.MixfMRIEnv$X.gbd[, i.p], PARAM$MU[i.p, i.k],
                               PARAM$SIGMA[[i.k]]$U[i.p], log = TRUE)
      }
    }
    ret[, i.k] <- log.b + log.a + PARAM$log.ETA[i.k]
  }

  ret <- rowSums(exp(ret))

  if(.MixfMRIEnv$CONTROL$debug > 10){
    .MixfMRIEnv$cat("  >>Not finite: ", sep = "", quiet = TRUE)
    base::cat(.MixfMRIEnv$COMM.RANK, ":", sum(!is.finite(ret)), " ", sep = "")
    .MixfMRIEnv$cat("\n", sep = "", quiet = TRUE)
  }

  ret <- sum.gbd(log(ret[is.finite(ret)]))
  ret
} # End of indep.logL.I().

indep.logL.V <- function(PARAM){
  nrow <- nrow(.MixfMRIEnv$X.gbd)
  p.X <- PARAM$p.X

  ret <- matrix(0, nrow = nrow, ncol = PARAM$K)
  for(i.k in 1:PARAM$K){
    log.a <- dbeta(.MixfMRIEnv$PV.gbd, PARAM$BETA[[i.k]][1],
                   PARAM$BETA[[i.k]][2], log = TRUE)
    log.b <- 0.0
    if(p.X > 0){
      log.b <- dmvn(.MixfMRIEnv$X.gbd, PARAM$MU[, i.k],
                    PARAM$SIGMA[[i.k]]$full, log = TRUE)
    }
    ret[, i.k] <- log.b + log.a + PARAM$log.ETA[i.k]
  }

  ret <- rowSums(exp(ret))

  if(.MixfMRIEnv$CONTROL$debug > 10){
    .MixfMRIEnv$cat("  >>Not finite: ", sep = "", quiet = TRUE)
    base::cat(.MixfMRIEnv$COMM.RANK, ":", sum(!is.finite(ret)), " ", sep = "")
    .MixfMRIEnv$cat("\n", sep = "", quiet = TRUE)
  }

  ret <- sum.gbd(log(ret[is.finite(ret)]))
  ret
} # End of indep.logL.V().

indep.logL.pv <- indep.logL.I
