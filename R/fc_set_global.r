### This function initializes global variables.

set.global <- function(X.gbd, PV.gbd, K = 2,
    min.1st.prop = .FC.CT$INIT$min.1st.prop,
    max.PV = .FC.CT$INIT$max.PV,
    class.method = .FC.CT$INIT$class.method[1],
    RndEM.iter = .FC.CT$CONTROL$RndEM.iter,
    algorithm = .FC.CT$algorithm[1],
    model.X = .FC.CT$model.X[1],
    ignore.X = .FC.CT$ignore.X,
    check.X.unit = .FC.CT$check.X.unit,
    MPI.gbd = .FC.CT$MPI.gbd, common.gbd = .FC.CT$common.gbd){
  ### Check options.
  if(K > 1){
    if(min.1st.prop < 0 || min.1st.prop >= 1){
      stop(paste("min.1st.prop is out of range [0, 1) for K > 1.", sep = ""))
    }
  } else{
    min.1st.prop <- 1
  }
  if(max.PV <= 0 || max.PV > 1){
    stop(paste("max.PV is out of range (0, 1].", sep = ""))
  }

  ### check MPI.
  if(MPI.gbd && !is.loaded("pkg_initialize", PACKAGE = "pbdMPI")){
    warning("MPI is required for MPI.gbd = TRUE.")
    MPI.gbd <- FALSE
  }
  .MixfMRIEnv$MPI.gbd <- MPI.gbd

  ### Rearrange data.
  if(.MixfMRIEnv$MPI.gbd){
    if(common.gbd){
      if(pbdMPI::spmd.comm.rank() != 0){
        X.gbd <- matrix(0, nrow = 0, ncol = ncol(X.gbd))
        PV.gbd <- matrix(0, nrow = 0, ncol = 1)
      }
    }
    X.gbd <- pbdMPI::comm.load.balance(X.gbd)
    PV.gbd <- as.vector(pbdMPI::comm.load.balance(matrix(PV.gbd, ncol = 1)))
  }

  ### Check data.
  if(! all(PV.gbd >= 0 & PV.gbd <= 1)){
    stop("PV.gbd should be in [0, 1].")
  }
  if(check.X.unit){
    if(! all(X.gbd >= 0 & X.gbd <= 1)){
      stop("X.gbd should be in [0, 1].")
    }
  }
  if(ignore.X){
    .MixfMRIEnv$X.gbd <- matrix(0, nrow = nrow(X.gbd), ncol = 0)
  } else{
    .MixfMRIEnv$X.gbd <- X.gbd
  }
  .MixfMRIEnv$PV.gbd <- PV.gbd
  .MixfMRIEnv$log.PV.gbd <- log(PV.gbd)
  .MixfMRIEnv$log.1.PV.gbd <- log(1 - PV.gbd)

  ### Get data information.
  N.gbd <- nrow(.MixfMRIEnv$X.gbd)
  if(.MixfMRIEnv$MPI.gbd){
    N.all <- pbdMPI::spmd.allgather.integer(as.integer(N.gbd),
                                            integer(pbdMPI::spmd.comm.size()))
    N <- sum(N.all)
  } else{
    N.all <- N.gbd
    N <- N.gbd
  }
  p.X <- ncol(.MixfMRIEnv$X.gbd)
  p.PV <- 1
  p <- p.X + p.PV

  ### Duplicate and rescale data for initialization only.
  ### 1. put pv and x to (0.001, 0.999).
  ### 2. take qnorm with mean 0.5 to get new scaled data.
  if(class.method %in%
     c("qnorm.simple", "qnorm.extend", "prob.simple", "prob.extend")){
    PV.max <- max.gbd(PV.gbd)
    PV.min <- min.gbd(PV.gbd)
    qnorm.PV.gbd <- (PV.gbd - PV.min + 0.001) / (PV.max - PV.min) * 0.998
    qnorm.PV.gbd <- qnorm(qnorm.PV.gbd, 0.5)

    max.X <- apply(X.gbd, 2, function(x){ max.gbd(x) })
    min.X <- apply(X.gbd, 2, function(x){ min.gbd(x) })
    qnorm.X.gbd <- W.plus.y(X.gbd, -min.X + 0.001,
                            nrow(X.gbd), ncol(X.gbd))
    qnorm.X.gbd <- W.divided.by.y(qnorm.X.gbd, max.X - min.X,
                                  nrow(X.gbd), ncol(X.gbd)) * 0.998
    qnorm.X.gbd <- apply(qnorm.X.gbd, 2, function(x){ qnorm(x, 0.5) })

    .MixfMRIEnv$qnorm.PV.gbd <- qnorm.PV.gbd
    .MixfMRIEnv$qnorm.X.gbd <- qnorm.X.gbd
  }

  ### Set global storages.
  .MixfMRIEnv$CONTROL <- .FC.CT$CONTROL
  .MixfMRIEnv$CONTROL$RndEM.iter <- RndEM.iter

  .MixfMRIEnv$Z.gbd <- matrix(0.0, nrow = N.gbd, ncol = K)
  .MixfMRIEnv$Z.colSums <- rep(0.0, K)

  .MixfMRIEnv$W.gbd <- matrix(0.0, nrow = N.gbd, ncol = K)
  .MixfMRIEnv$W.gbd.rowSums <- rep(0.0, N.gbd)

  .MixfMRIEnv$U.gbd <- matrix(0.0, nrow = N.gbd, ncol = K)
  .MixfMRIEnv$CLASS.gbd <- rep(0, N.gbd)
  .MixfMRIEnv$CHECK <- list(algorithm = algorithm[1], i.iter = 0,
                              abs.err = Inf, rel.err = Inf, convergence = 0)

  ### Set global constants.
  .MixfMRIEnv$p.times.logtwopi <- p.X * log(2 * pi)

  ### Set verbose functions.
  if(.MixfMRIEnv$MPI.gbd){
    .MixfMRIEnv$cat <- pbdMPI::comm.cat
    .MixfMRIEnv$print <- pbdMPI::comm.print
    .MixfMRIEnv$any <- pbdMPI::comm.any
  } else{
    .MixfMRIEnv$cat <- function(..., quiet = TRUE){
                           base::cat(...)
                         }
    .MixfMRIEnv$print <- function(x, ..., quiet = TRUE){
                             base::print(x, ...)
                           }
    .MixfMRIEnv$any <- function(..., na.rm = FALSE){
                           base::any(..., na.rm = na.rm)
                         }
  }

  ### Set functions for variance covariance matrix.
  if(model.X[1] == "I"){
    .MixfMRIEnv$initial.SIGMA <- initial.SIGMA.I
    .MixfMRIEnv$var.gbd <- var.gbd.I
    .MixfMRIEnv$mle.SIGMA <- mle.SIGMA.I
    .MixfMRIEnv$logbetanorm <- logbetanorm.I
    .MixfMRIEnv$indep.logL <- indep.logL.I
    .MixfMRIEnv$decompsigma <- decompsigma.I
  } else if(model.X[1] == "V"){
    .MixfMRIEnv$initial.SIGMA <- initial.SIGMA.V
    .MixfMRIEnv$var.gbd <- var.gbd.V
    .MixfMRIEnv$mle.SIGMA <- mle.SIGMA.V
    .MixfMRIEnv$logbetanorm <- logbetanorm.V
    .MixfMRIEnv$indep.logL <- indep.logL.V
    .MixfMRIEnv$decompsigma <- decompsigma.V
  } else{
    stop("model.X is not implemented.")
  }

  ### Set SIGMA.
  if(ncol(.MixfMRIEnv$X.gbd) > 0){
    full <- .MixfMRIEnv$var.gbd(X.gbd)
    tmp.U <- .MixfMRIEnv$decompsigma(full)
    SIGMA.1 <- list(full = full, U = tmp.U$value, U.check = tmp.U$check)
    if(K > 1){
      full <- full / (K - 1)
    }
    tmp.U <- .MixfMRIEnv$decompsigma(full)
    SIGMA.rest <- list(full = full, U = tmp.U$value, U.check = tmp.U$check)
  } else{
    full <- matrix(0, nrow = 0, ncol = 0)
    SIGMA.1 <- list(full = full, U = NULL, U.check = NULL)
    SIGMA.rest <- SIGMA.1
  }
  SIGMA <- rep(c(list(SIGMA.1), list(SIGMA.rest)), c(1, (K - 1)))

  ### Set parameters.
  if(K > 1){
    ETA <- c(min.1st.prop,
             rep((1 - min.1st.prop) / (K - 1), K - 1))
  } else{
    ETA <- 1
  }
  PARAM <- list(N.gbd = N.gbd, N.all = N.all,
                N = N, p = p, p.X = p.X, p.PV = p.PV,
                K = K,
                ETA = ETA,
                log.ETA = NULL,
                BETA = rep(c(list(c(1, 1)),
                             list(c(0.5, 2))),
                           c(1, (K - 1))),
                MU = matrix(0, nrow = p.X, ncol = K),
                SIGMA = SIGMA,
                logL = NULL,
                min.1st.prop = min.1st.prop,
                max.PV = max.PV,
                class.method = class.method,
                min.N.CLASS = p + 1,
                model.X = model.X[1])
  PARAM$log.ETA <- log(PARAM$ETA)
  invisible(PARAM)
} # End of set.global().
