### This file gives a prob extend initialization.

initial.em.gbd.prob.extend <- function(PARAM){
  N.all <- PARAM$N.all
  N.gbd <- PARAM$N.gbd
  K <- PARAM$K
  p.X <- PARAM$p.X

  ### For the 1st cluster.
  gcenter.X <- rep(0, p.X)
  gcenter.PV <- 0

  ### For the rests.
  id.PV.gbd <- which(.MixfMRIEnv$PV.gbd < PARAM$max.PV)
  if(length_gbd(id.PV.gbd) < PARAM$min.N.CLASS * (K - 1)){
    ### Avoid too tiny subset or too large K.
    stop("K may be too large.")
  }
  tmp.X.gbd <- .MixfMRIEnv$qnorm.X.gbd[id.PV.gbd,]
  tmp.PV.gbd <- .MixfMRIEnv$qnorm.PV.gbd[id.PV.gbd]

  ### Get comm size and rank if MPI is loaded.
  if(.MixfMRIEnv$MPI.gbd){
    COMM.SIZE <- pbdMPI::spmd.comm.size()
    COMM.RANK <- pbdMPI::spmd.comm.rank()

    ### Stupid but work.
    all.tmp.X.gbd <- pbdMPI::gather(tmp.X.gbd)
    all.tmp.PV.gbd <- pbdMPI::gather(tmp.PV.gbd)
    if(COMM.RANK == 0){
      all.tmp.X.gbd <- do.call("rbind", all.tmp.X.gbd)
      all.tmp.PV.gbd <- do.call("c", all.tmp.PV.gbd)
      prob.center <- all.tmp.PV.gbd^2
      prob.center <- prob.center / sum(prob.center)
    }
  } else{
    prob.center <- tmp.PV.gbd^2
    prob.center <- prob.center / sum(prob.center)
  }

  ### Random centering K - 1 points.
  iter <- 0
  repeat{
    if(.MixfMRIEnv$MPI.gbd){
      if(COMM.RANK == 0){
        id.center <- sample(1:length(prob.center), K - 1, prob = prob.center)

        center.X <- matrix(all.tmp.X.gbd[id.center,], ncol = p.X)
        center.X <- pbdMPI::spmd.bcast.double(center.X)
        center.PV <- all.tmp.PV.gbd[id.center]
        center.PV <- pbdMPI::spmd.bcast.double(center.PV)
      } else{
        center.X <- pbdMPI::spmd.bcast.double(rep(0.0, (K - 1) * p.X))
        center.X <- matrix(center.X, ncol = p.X)
        center.PV <- pbdMPI::spmd.bcast.double(rep(0.0, K - 1))
      }

      ### Gather all centers and pv's.
      center.X <- matrix(rbind(gcenter.X, center.X), ncol = p.X)
      center.PV <- c(gcenter.PV, center.PV)
    } else{
      id.center <- sample(1:length(prob.center), K - 1, prob = prob.center)
      center.X <- matrix(rbind(gcenter.X, tmp.X.gbd[id.center,]), ncol = p.X)
      center.PV <- c(gcenter.PV, tmp.PV.gbd[id.center])
    }

    ### Classifying to the closest center.
    CLASS.gbd <- rep(1, N.gbd)
    tmp.dist <- NULL
    for(i.k in 1:K){
      if(PARAM$p.X > 0){
        tmp <- rowMeans(W.plus.y(tmp.X.gbd, -center.X[i.k,],
                                 length(id.PV.gbd), p.X)^2)
      } else{
        tmp <- 0
      }
      tmp <- tmp + (tmp.PV.gbd - center.PV[i.k])^2
      tmp.dist <- cbind(tmp.dist, tmp)
    }
    CLASS.gbd[id.PV.gbd] <- apply(tmp.dist, 1, which.min)
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
} # End of initial.em.gbd.prob.extend().
