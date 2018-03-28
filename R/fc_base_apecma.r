### This file contains major functions for EM iterations.

### E-step.
apea.step.gbd.k <- function(PARAM, i.k, update.logL = TRUE){
  .MixfMRIEnv$logbetanorm(PARAM, i.k)
  compute.expectation(PARAM, update.logL = update.logL)
  invisible()
} # End of apea.step.gbd.k().

### CM-step
cm.step.gbd.ETA.BETA.MU.SIGMA.k <- function(PARAM, i.k){
  ### MLE For ETA
  PARAM$ETA <- .MixfMRIEnv$Z.colSums / sum(.MixfMRIEnv$Z.colSums)
  if(PARAM$ETA[1] <= PARAM$min.1st.prop){
    PARAM$ETA <- c(PARAM$min.1st.prop,
                   (1 - PARAM$min.1st.prop) * .MixfMRIEnv$Z.colSums[-1] /
                   sum(.MixfMRIEnv$Z.colSums[-1]))
  }
  PARAM$log.ETA <- log(PARAM$ETA)

  ### MLE for BETA
  PARAM$BETA[[i.k]] <- cm.step.gbd.BETA.k(PARAM, i.k)

  ### MLE for MU and SIGMA
  if(PARAM$p.X > 0){
    PARAM$MU[, i.k] <- cm.step.gbd.MU.k(PARAM, i.k)
    PARAM$SIGMA[[i.k]] <- cm.step.gbd.SIGMA.k(PARAM, i.k)
  }

  PARAM
} # End of cm.step.gbd.ETA.BETA.MU.SIGMA.k().

cm.step.gbd.BETA.k <- function(PARAM, i.k){
  #### MLE for BETA
  if(i.k != 1){
    ret <- mle.beta.constr(PARAM, i.k)
  } else{
    ret <- c(1.0, 1.0)
  }

  ret
} # End of cm.step.gbd.BETA.k().

cm.step.gbd.MU.k <- function(PARAM, i.k){
  ### MLE for MU
  ret <- colSums.gbd(.MixfMRIEnv$X.gbd * .MixfMRIEnv$Z.gbd[, i.k] /
                     .MixfMRIEnv$Z.colSums[i.k])
  ret
} # End of cm.step.gbd.MU.k().

cm.step.gbd.SIGMA.k <- function(PARAM, i.k){
  U <- PARAM$SIGMA[[i.k]]$U
  full <- PARAM$SIGMA[[i.k]]$full

  ### MLE for SIGMA
  tmp.full <- .MixfMRIEnv$mle.SIGMA(PARAM, i.k)
  ### If SIGMA is degenerated, then there is no update.
  tmp.U <- .MixfMRIEnv$decompsigma(tmp.full)
  U.check <- tmp.U$check
  if(tmp.U$check){
    U <- tmp.U$value
    full <- tmp.full
  }

  ret <- list(full = full, U = U, U.check = U.check)
  ret
} # End of cm.step.gbd.SIGMA.k().


### APECMa-step.
apecma.step.gbd <- function(PARAM.org){
  .MixfMRIEnv$CHECK <- list(algorithm = "apecma", i.iter = 0, abs.err = Inf,
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

    PARAM.new <- try(apecma.onestep.gbd(PARAM.org))
    if(.MixfMRIEnv$any(class(PARAM.new) == "try-error")){
      .MixfMRIEnv$cat("Results of previous iterations are returned.\n",
                        quiet =TRUE)
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
} # End of apecma.step.gbd().

apecma.onestep.gbd <- function(PARAM){
  for(i.k in 1:PARAM$K){
    PARAM <- cm.step.gbd.ETA.BETA.MU.SIGMA.k(PARAM, i.k)
    apea.step.gbd.k(PARAM, i.k,
                    update.logL = ifelse(i.k == PARAM$K, TRUE, FALSE))
  }

  PARAM$logL <- logL.step.gbd()

  if(.MixfMRIEnv$CONTROL$debug > 0){
    .MixfMRIEnv$cat(">>apecma.onestep: ", format(Sys.time(), "%H:%M:%S"),
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
} # End of apecma.onestep.gbd().

