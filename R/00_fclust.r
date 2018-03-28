### For general methods.

fclust <- function(X.gbd, PV.gbd, K = 2,
    PARAM.init = NULL,
    min.1st.prop = .FC.CT$INIT$min.1st.prop,
    max.PV = .FC.CT$INIT$max.PV,
    class.method = .FC.CT$INIT$class.method[1],
    RndEM.iter = .FC.CT$CONTROL$RndEM.iter,
    algorithm = .FC.CT$algorithm[1],
    model.X = .FC.CT$model.X[1],
    ignore.X = .FC.CT$ignore.X,
    stop.unstable = TRUE,
    MPI.gbd = .FC.CT$MPI.gbd, common.gbd = .FC.CT$common.gbd){
  if(any(PV.gbd == 0 | PV.gbd == 1)){
    stop("PV.gbd should be in (0, 1) and exclude 0 and 1.")
  }

  if(algorithm[1] %in% .FC.CT$algorithm){
    if(is.null(PARAM.init)){
      ### guess initial values.
      PARAM.org <- set.global(X.gbd, PV.gbd, K = K,
                              min.1st.prop = min.1st.prop,
                              max.PV = max.PV,
                              class.method = class.method[1],
                              RndEM.iter = RndEM.iter,
                              algorithm = algorithm[1],
                              model.X = model.X[1],
                              ignore.X = ignore.X,
                              MPI.gbd = MPI.gbd, common.gbd = common.gbd)
      PARAM <- initial.RndEM.gbd(PARAM.org)
    } else{
      ### use initial values.
      PARAM.org <- set.global(X.gbd, PV.gbd, K = PARAM.init$K,
                              min.1st.prop = PARAM.init$min.1st.prop,
                              max.PV = PARAM.init$max.PV,
                              class.method = class.method[1],
                              RndEM.iter = 0,
                              algorithm = algorithm[1],
                              model.X = model.X[1],
                              ignore.X = ignore.X,
                              MPI.gbd = MPI.gbd, common.gbd = common.gbd)
      PARAM <- PARAM.org
      PARAM$ETA <- PARAM.init$ETA
      PARAM$BETA <- PARAM.init$BETA
      PARAM$MU <- PARAM.init$MU
      PARAM$SIGMA <- PARAM.init$SIGMA
      e.step.gbd(PARAM)

      em.update.class.gbd()
    }

    ### Check # of class.
    N.CLASS <- get.N.CLASS(PARAM$K)
    if(any(N.CLASS < PARAM$min.N.CLASS)){
      if(stop.unstable){
        stop("Initialization may be unstable.")
      }
    }

    ### Check SIGMA
    if(PARAM$p.X > 0){
      tmp.check <- lapply(PARAM$SIGMA, function(x){ x$U.check })
      if(!all(do.call("c", tmp.check))){
        if(stop.unstable){
          stop("Initialization may be unstable.")
        }
      }
    }

    ### Run from the stable initialization.
    method.step <- switch(algorithm[1],
                          "em" = em.step.gbd,
                          "apecma" = apecma.step.gbd,
                          "ecm" = ecm.step.gbd,
                          NULL)
    PARAM.new <- method.step(PARAM)
    PARAM.new$entropy <- entropy.step.gbd()
    em.update.class.gbd()

    # Get class numbers.
    N.CLASS <- get.N.CLASS(K)

    # Get IC's.
    AIC <- aic(PARAM.new)
    BIC <- bic(PARAM.new)
    ICL.BIC <- icl.bic(PARAM.new)

    # For return.
    ret <- list(algorithm = algorithm[1],
                model.X = model.X[1],
                ignore.X = ignore.X,
                param = PARAM.new,
                check = .MixfMRIEnv$CHECK,
                n.class = N.CLASS,
                aic = AIC,
                bic = BIC,
                icl.bic = ICL.BIC,
                class = .MixfMRIEnv$CLASS.gbd,
                posterior = .MixfMRIEnv$Z.gbd)
  } else{
    stop("The algorithm is not found.")
  }

  class(ret) <- "fclust"
  ret
} # end of fclust().

