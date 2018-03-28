### This file contains major functions for EM iterations.

### E-step.
e.step.gbd <- function(PARAM, update.logL = TRUE){
  for(i.k in 1:PARAM$K){
    .MixfMRIEnv$logbetanorm(PARAM, i.k)
  }

  compute.expectation(PARAM, update.logL = update.logL)
  invisible()
} # End of e.step.gbd().


### z_nk / sum_k z_n might have numerical problems if z_nk all underflowed.
compute.expectation <- function(PARAM, update.logL = TRUE){
  N.gbd <- PARAM$N.gbd
  K <- PARAM$K

  .MixfMRIEnv$U.gbd <- W.plus.y(.MixfMRIEnv$W.gbd, PARAM$log.ETA, N.gbd, K)
  .MixfMRIEnv$Z.gbd <- exp(.MixfMRIEnv$U.gbd)

  tmp.id <- rowSums(.MixfMRIEnv$U.gbd < .MixfMRIEnv$CONTROL$exp.min) == K |
            rowSums(.MixfMRIEnv$U.gbd > .MixfMRIEnv$CONTROL$exp.max) > 0

  tmp.flag <- sum(tmp.id)
  if(tmp.flag > 0){
    tmp.gbd <- .MixfMRIEnv$U.gbd[tmp.id,]

    if(tmp.flag == 1){
      tmp.scale <- max(tmp.gbd) - .MixfMRIEnv$CONTROL$exp.max / K
    } else{
      tmp.scale <- unlist(apply(tmp.gbd, 1, max)) -
                   .MixfMRIEnv$CONTROL$exp.max / K
    }
    .MixfMRIEnv$Z.gbd[tmp.id,] <- exp(tmp.gbd - tmp.scale)
  }

  .MixfMRIEnv$W.gbd.rowSums <- rowSums(.MixfMRIEnv$Z.gbd)
  .MixfMRIEnv$Z.gbd <- .MixfMRIEnv$Z.gbd / .MixfMRIEnv$W.gbd.rowSums

  .MixfMRIEnv$Z.colSums <- colSums.gbd(.MixfMRIEnv$Z.gbd)

  if(update.logL){
    .MixfMRIEnv$W.gbd.rowSums <- log(.MixfMRIEnv$W.gbd.rowSums)
    if(tmp.flag > 0){
      .MixfMRIEnv$W.gbd.rowSums[tmp.id] <-
        .MixfMRIEnv$W.gbd.rowSums[tmp.id] + tmp.scale
    }
  }

  invisible()
} # End of compute.expectation().


### M-step.
m.step.gbd <- function(PARAM){
  ### MLE For ETA
  PARAM$ETA <- .MixfMRIEnv$Z.colSums / sum(.MixfMRIEnv$Z.colSums)
  if(PARAM$ETA[1] <= PARAM$min.1st.prop){
    PARAM$ETA <- c(PARAM$min.1st.prop,
                   (1 - PARAM$min.1st.prop) * .MixfMRIEnv$Z.colSums[-1] /
                   sum(.MixfMRIEnv$Z.colSums[-1]))
  }
  PARAM$log.ETA <- log(PARAM$ETA)

  ### MLE for BETA, MU, and SIGMA
  p.X <- PARAM$p.X
  for(i.k in 1:PARAM$K){
    #### MLE for BETA
    PARAM$BETA[[i.k]] <- cm.step.gbd.BETA.k(PARAM, i.k)

    ### MLE for MU and SIGMA
    if(PARAM$p.X > 0){
      PARAM$MU[, i.k] <- cm.step.gbd.MU.k(PARAM, i.k)
      PARAM$SIGMA[[i.k]] <- cm.step.gbd.SIGMA.k(PARAM, i.k)
    }
  }

  PARAM
} # End of m.step.gbd().


### log likelihood.
logL.step.gbd <- function(){
  sum.gbd(.MixfMRIEnv$W.gbd.rowSums)
} # End of logL.step.gbd().

### entropy.
entropy.step.gbd <- function(){
  tmp.gbd <- .MixfMRIEnv$Z.gbd * log(.MixfMRIEnv$Z.gbd)
  tmp.gbd[! is.finite(tmp.gbd)] <- 0
  -sum.gbd(tmp.gbd)
} # End of entropy.step.gbd().


### Check log likelihood convergence.
check.em.convergence <- function(PARAM.org, PARAM.new, i.iter){
  abs.err <- PARAM.new$logL - PARAM.org$logL
  rel.err <- abs.err / abs(PARAM.org$logL)
  convergence <- 0

  if(abs.err < 0){
    convergence <- 4
  } else if(i.iter > .MixfMRIEnv$CONTROL$max.iter){
    convergence <- 2
  } else if(rel.err < .MixfMRIEnv$CONTROL$rel.err){
    convergence <- 1
  }

  if(.MixfMRIEnv$CONTROL$debug > 1){
    .MixfMRIEnv$cat("  check.em.convergence:",
                      " abs: ", abs.err,
                      ", rel: ", rel.err,
                      ", conv: ", convergence, "\n",
                      sep = "", quiet = TRUE)
  }

  list(algorithm = .MixfMRIEnv$CHECK$algorithm,
       iter = i.iter, abs.err = abs.err, rel.err = rel.err,
       convergence = convergence)
} # End of check.em.convergence().


### EM-step.
em.step.gbd <- function(PARAM.org){
  .MixfMRIEnv$CHECK <- list(algorithm = "em", i.iter = 0, abs.err = Inf,
                              rel.err = Inf, convergence = 0)
  i.iter <- 1
  PARAM.org$logL <- -.Machine$double.xmax

  ### For debugging.
  if((!is.null(.MixfMRIEnv$CONTROL$save.log)) &&
     .MixfMRIEnv$CONTROL$save.log){
    if(! exists("SAVE.iter", envir = .MixfMRIEnv)){
      .MixfMRIEnv$SAVE.param <- NULL
      .MixfMRIEnv$SAVE.iter <- NULL
      .MixfMRIEnv$CLASS.iter.org <- unlist(apply(.MixfMRIEnv$Z.gbd, 1,
                                                   which.max))
    }
  }

  repeat{
    ### For debugging.
    if((!is.null(.MixfMRIEnv$CONTROL$save.log)) &&
        .MixfMRIEnv$CONTROL$save.log){
      time.start <- proc.time()
    }

    PARAM.new <- try(em.onestep.gbd(PARAM.org))
    if(.MixfMRIEnv$any(class(PARAM.new) == "try-error")){
      .MixfMRIEnv$cat("Results of previous iterations are returned.\n",
                        quiet = TRUE)
      .MixfMRIEnv$CHECK$convergence <- 99
      PARAM.new <- PARAM.org
      break
    }

    .MixfMRIEnv$CHECK <- check.em.convergence(PARAM.org, PARAM.new, i.iter)
    if(.MixfMRIEnv$CHECK$convergence > 0){
      break
    }

    ### For debugging.
    if((!is.null(.MixfMRIEnv$CONTROL$save.log)) &&
        .MixfMRIEnv$CONTROL$save.log){
      tmp.time <- proc.time() - time.start

      .MixfMRIEnv$SAVE.param <- c(.MixfMRIEnv$SAVE.param, PARAM.new)
      CLASS.iter.new <- unlist(apply(.MixfMRIEnv$Z.gbd, 1, which.max))
      tmp <- sum.gbd(CLASS.iter.new != .MixfMRIEnv$CLASS.iter.org)

      tmp.all <- c(tmp / PARAM.new$N, PARAM.new$logL,
                   PARAM.new$logL - PARAM.org$logL,
                   (PARAM.new$logL - PARAM.org$logL) / PARAM.org$logL)
      .MixfMRIEnv$SAVE.iter <- rbind(.MixfMRIEnv$SAVE.iter,
                                     c(tmp, tmp.all, tmp.time))
      .MixfMRIEnv$CLASS.iter.org <- CLASS.iter.new
    }

    PARAM.org <- PARAM.new
    i.iter <- i.iter + 1
  }

  PARAM.new
} # End of em.step.gbd().


em.onestep.gbd <- function(PARAM){
  PARAM <- m.step.gbd(PARAM)
  e.step.gbd(PARAM)

  PARAM$logL <- logL.step.gbd()

  if(.MixfMRIEnv$CONTROL$debug > 0){
    .MixfMRIEnv$cat(">>em.onestep: ", format(Sys.time(), "%H:%M:%S"),
                      ", iter: ", .MixfMRIEnv$CHECK$iter, ", logL: ",
                                  sprintf("%-30.15f", PARAM$logL), "\n",
                      sep = "", quiet = TRUE)
    if(.MixfMRIEnv$CONTROL$debug > 4){
      logL <- .MixfMRIEnv$indep.logL(PARAM)
      .MixfMRIEnv$cat("  >>indep.logL: ", sprintf("%-30.15f", logL), "\n",
                        sep = "", quiet = TRUE)
    }
    if(.MixfMRIEnv$CONTROL$debug > 20){
      mb.print(PARAM, .MixfMRIEnv$CHECK)
    }
  }

  PARAM
} # End of em.onestep.gbd().


### Obtain classifications.
em.update.class.gbd <- function(){
  .MixfMRIEnv$CLASS.gbd <- unlist(apply(.MixfMRIEnv$Z.gbd, 1, which.max))
  invisible()
} # End of em.update.class.gbd().

