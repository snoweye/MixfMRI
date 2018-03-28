### This file gives a simple initialization.

initial.em.gbd.simple <- function(PARAM){
  N.all <- PARAM$N.all
  N.gbd <- PARAM$N.gbd
  K <- PARAM$K
  p.X <- PARAM$p.X

  ### For the 1st cluster.
  gcenter.X <- colMeans.gbd(.MixfMRIEnv$X.gbd)
  gcenter.PV <- 0.5

  ### For the rests.
  id.PV.gbd <- which(.MixfMRIEnv$PV.gbd < PARAM$max.PV)
  if(length.gbd(id.PV.gbd) < PARAM$min.N.CLASS * (K - 1)){
    ### Avoid too tiny subset or too large K.
    stop("K may be too large.")
  }
  tmp.X.gbd <- .MixfMRIEnv$X.gbd[id.PV.gbd,]
  tmp.PV.gbd <- .MixfMRIEnv$PV.gbd[id.PV.gbd]

  ### Get comm size and rank if MPI is loaded.
  if(.MixfMRIEnv$MPI.gbd){
    COMM.SIZE <- pbdMPI::spmd.comm.size()
    COMM.RANK <- pbdMPI::spmd.comm.rank()
  }

  ### Random centering K - 1 points.
  iter <- 0
  repeat{
    if(.MixfMRIEnv$MPI.gbd){
      tmp.count <- as.integer(sample(1:COMM.SIZE, K - 1, replace = TRUE,
                                     prob = N.all) - 1)
      tmp.count <- pbdMPI::spmd.bcast.integer(tmp.count)
      n.center <- sum(tmp.count == COMM.RANK)

      ### Sample by rank.
      if(n.center > 0){
        id.center <- sample.int(length(id.PV.gbd), n.center, replace = FALSE)
        center.X <- matrix(tmp.X.gbd[id.center,], ncol = p.X)
        center.PV <- tmp.PV.gbd[id.center]
      } else{
        center.X <- NULL
        center.PV <- NULL
      }

      ### Gather all centers and pv's.
      center.X <- do.call("rbind", pbdMPI::allgather(center.X))
      center.X <- matrix(rbind(gcenter.X, center.X), ncol = p.X)
      center.PV <- do.call("c", pbdMPI::allgather(center.PV))
      center.PV <- c(gcenter.PV, center.PV)
    } else{
      id.center <- sample.int(length(id.PV.gbd), K - 1)
      center.X <- matrix(rbind(gcenter.X, tmp.X.gbd[id.center,]), ncol = p.X)
      center.PV <- c(gcenter.PV, tmp.PV.gbd[id.center])
    }

    ### Classifying to the closest center.
    CLASS.gbd <- rep(1, N.gbd)
    if(PARAM$max.PV == 1){
      K.start <- 1
    } else{
      K.start <- 2
    }
    tmp.dist <- NULL
    for(i.k in K.start:K){
      if(PARAM$p.X > 0){
        tmp <- rowMeans(W.plus.y(tmp.X.gbd, -center.X[i.k,],
                                 length(id.PV.gbd), p.X)^2)
      } else{
        tmp <- 0
      }
      tmp <- tmp + (tmp.PV.gbd - center.PV[i.k])^2
      tmp.dist <- cbind(tmp.dist, tmp)
    }
    CLASS.gbd[id.PV.gbd] <- apply(tmp.dist, 1, which.min) + (K.start - 1)
    check <- all(tabulate.gbd(CLASS.gbd, nbins = K) > PARAM$min.N.CLASS)

    if(.MixfMRIEnv$CONTROL$debug > 2){
      tmp <- tabulate.gbd(CLASS.gbd, nbins = K)
      .MixfMRIEnv$cat("n.class: ",
                        paste(tmp, collapse = " "),
                        "\n", sep = "", quiet = TRUE)
    }

    if(check){
      break
    } else{
      iter <- iter + 1
    }
    if(iter > .FC.CT$INIT$max.try.iter){
      .MixfMRIEnv$cat("max.try.iter reached.\n")
      break
    }
  }

  initial.em.gbd.update(PARAM, CLASS.gbd, K)
} # End of initial.em.gbd.simple().
