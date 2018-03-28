### This file contains major functions for ECM iterations.

### CM-step.
cm.step.gbd.ETA.1st <- function(PARAM){
  ### MLE For ETA[c(1, K)].
  w <- 1 - sum(PARAM$ETA[-c(1, PARAM$K)])
  tmp <- w * .MixfMRIEnv$Z.colSums[c(1, PARAM$K)] /
         sum(.MixfMRIEnv$Z.colSums[c(1, PARAM$K)])
  if(tmp[1] <= PARAM$min.1st.prop){
    PARAM$ETA[c(1, PARAM$K)] <- c(PARAM$min.1st.prop,
                                  w - PARAM$min.1st.prop)
  } else{
    PARAM$ETA[c(1, PARAM$K)] <- tmp
  }
  PARAM$log.ETA <- log(PARAM$ETA)

  PARAM
} # End of cm.step.gbd.ETA.1st().

cm.step.gbd.ETA.rest <- function(PARAM){
  ### MLE For ETA[-1]
  PARAM$ETA[-1] <- (1 - PARAM$ETA[1]) * .MixfMRIEnv$Z.colSums[-1] /
                   sum(.MixfMRIEnv$Z.colSums[-1])
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
} # End of cm.step.gbd.ETA.rest().


### ECM-step.
ecm.step.gbd <- function(PARAM.org){
  .MixfMRIEnv$CHECK <- list(algorithm = "ecm", i.iter = 0, abs.err = Inf,
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

    PARAM.new <- try(ecm.onestep.gbd(PARAM.org))
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
} # End of ecm.step.gbd().


ecm.onestep.gbd <- function(PARAM){
  PARAM <- cm.step.gbd.ETA.1st(PARAM)
  e.step.gbd(PARAM)
  PARAM <- cm.step.gbd.ETA.rest(PARAM)
  e.step.gbd(PARAM)

  PARAM$logL <- logL.step.gbd()

  if(.MixfMRIEnv$CONTROL$debug > 0){
    .MixfMRIEnv$cat(">>ecm.onestep: ", format(Sys.time(), "%H:%M:%S"),
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
} # End of ecm.onestep.gbd().

