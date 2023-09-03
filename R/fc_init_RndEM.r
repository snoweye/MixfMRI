### This implements RndEM algorithm.

initial.RndEM.gbd <- function(PARAM){
  logL.save <- -Inf
  i.iter <- 1

  PARAM.org <- PARAM
  PARAM.save <- PARAM
  repeat{
    if(i.iter > .MixfMRIEnv$CONTROL$RndEM.iter){
      break
    }

    PARAM <- try(initial.em.gbd(PARAM.org), silent = TRUE)
    if(.MixfMRIEnv$any(inherits(PARAM, "try-error"))){
      .MixfMRIEnv$cat("The initial may be unstable. Skip to next initial.\n",
                        quiet = TRUE)
      if(.MixfMRIEnv$CONTROL$debug > 0){
        .MixfMRIEnv$cat(PARAM, "\n", quiet = TRUE)
      }
      i.iter <- i.iter + 1
      next
    }

    ### Check # of class.
    N.CLASS <- get.N.CLASS(PARAM$K)
    if(any(N.CLASS < PARAM$min.N.CLASS)){
      if(.MixfMRIEnv$CONTROL$debug > 0){
        .MixfMRIEnv$cat("N.CLASS: ", N.CLASS, "\n", quiet = TRUE)
      }
      i.iter <- i.iter + 1
      next
    }

    ### Check SIGMA
    if(PARAM$p.X > 0){
      tmp.check <- lapply(PARAM$SIGMA, function(x){ x$U.check })
      if(!all(do.call("c", tmp.check))){
        if(.MixfMRIEnv$CONTROL$debug > 0){
          .MixfMRIEnv$cat("SIGMA: some U.check fails.\n", quiet = TRUE)
        }
        i.iter <- i.iter + 1
        next
      }
    }

    if(.MixfMRIEnv$CONTROL$debug > 0){
      .MixfMRIEnv$cat("Initial: ", format(Sys.time(), "%H:%M:%S"),
                        ", iter: ", i.iter, ", logL: ",
                                    sprintf("%-20.10f", PARAM$logL), "\n",
                        sep = "", quiet = TRUE)
    }

    if(logL.save < PARAM$logL){
      logL.save <- PARAM$logL
      PARAM.save <- PARAM
      PARAM.save$initial.i.iter <- i.iter
    }

    i.iter <- i.iter + 1
  }

  if(.MixfMRIEnv$CONTROL$debug > 0){
    .MixfMRIEnv$cat("Using initial iter: ", PARAM.save$initial.i.iter, "\n",
                      sep = "", quiet = TRUE)
  }

  PARAM <- PARAM.save
  e.step.gbd(PARAM)
  PARAM <- em.onestep.gbd(PARAM)
  PARAM$logL <- logL.step.gbd()
  em.update.class.gbd()

  PARAM
} # End of initial.RndEM.gbd().

