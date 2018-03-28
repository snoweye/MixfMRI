### This file is modified from EMcluster.


### partial logL. Return a matrix with dimension M * N.
partial.logL.pv <- function(x, ETA, BETA, MU, SIGMA, post.z){
  K <- length(ETA)
  N <- dim(x$X.gbd)[1]

  BETA.digamma <- lapply(BETA, function(b){ digamma(b) - digamma(sum(b)) })
  # SIGMA.inv <- lapply(SIGMA, function(s){ 1 / s$full })

  SS <- lapply(1:N, function(i){
    Sbuf.ETA <- NULL
    if(K != 1){
      tmp <- post.z[i, ] / ETA
      Sbuf.ETA <- tmp[-K] - tmp[K]
    }

    Sbuf.BETA <- list()
    for(j in 2:K){
      Sbuf.BETA[[j]] <- post.z[i, j] *
                        ((BETA[[j]] - 1) / x$PV.gbd[i] - BETA.digamma[[j]])
    }

    c(Sbuf.ETA, do.call("c", Sbuf.BETA))
  })

  do.call("cbind", SS)
} # End of partial.logL.pv().


### Generate dataset.
GenDataSet.pv <- function(N, ETA, BETA, MU, SIGMA, X.gbd){
  K <- ncol(MU)
  n <- nrow(X.gbd)
  SD <- do.call("cbind", lapply(SIGMA, function(s) sqrt(s$full)))

  p.k <- matrix(0, nrow = n, ncol = K)
  for(i.k in 1:K){
    p.k[, i.k] <- apply(X.gbd, 1,
                        function(x){
                          sum(dnorm(x, MU[, i.k], SD[, i.k], log = TRUE))
                        }) +
                  log(ETA[i.k])
  }
  p.k <- exp(p.k)
  id <- apply(p.k, 1, function(prob) sample(1:K, 1, prob = prob))

  ### For PV.gbd
  PV.gbd <- rep(0, length(id))
  for(i.k in 1:K){
    n <- sum(id == i.k)
    PV.gbd[id == i.k] <- rbeta(n, BETA[[i.k]][1], BETA[[i.k]][2])
  }

  list(x = list(X.gbd = X.gbd, PV.gbd = PV.gbd), id = id)
} # End of GenDataSet.pv().


### Obtain parameters.
get.E.chi2.pv <- function(x, fcobj.0, fcobj.a, given = c("0", "a"), tau = 0.5,
    n.mc = 1000, verbose = TRUE){
  K.0 <- fcobj.0$param$K
  ETA.0 <- fcobj.0$param$ETA
  BETA.0 <- fcobj.0$param$BETA
  MU.0 <- fcobj.0$param$MU
  SIGMA.0 <- fcobj.0$param$SIGMA

  K.a <- fcobj.a$param$K
  ETA.a <- fcobj.a$param$ETA
  BETA.a <- fcobj.a$param$BETA
  MU.a <- fcobj.a$param$MU
  SIGMA.a <- fcobj.a$param$SIGMA

  if(given[1] == "0"){
    ETA <- ETA.0
    BETA <- BETA.0
    MU <- MU.0
    SIGMA <- SIGMA.0
  } else if(given[1] == "a"){
    ETA <- ETA.a
    BETA <- BETA.a
    MU <- MU.a
    SIGMA <- SIGMA.a
  } else{
    stop("given should be '0' or 'a'.")
  }

  ### Obtain nabla logL via Monte Carlo.
  x.new <- GenDataSet.pv(n.mc, ETA, BETA, MU, SIGMA, x$X.gbd)$x
  post.z0 <- post.prob(x.new, fcobj.0)
  post.za <- post.prob(x.new, fcobj.a)
  pl0 <- partial.logL.pv(x.new, ETA.0, BETA.0, MU.0, SIGMA.0, post.z0)
  pla <- partial.logL.pv(x.new, ETA.a, BETA.a, MU.a, SIGMA.a, post.za)

  ### Drop eta_1 if on the boundary.
  if(fcobj.0$param$ETA[1] == fcobj.0$param$min.1st.prop){
    pl0 <- pl0[-1,]
  }
  if(fcobj.a$param$ETA[1] == fcobj.a$param$min.1st.prop){
    pla <- pla[-1,]
  }

  ### Expected degrees of freedom.
  nl <- rbind(pl0, pla)
  mu <- rowMeans(nl)
  nl <- nl - mu
  J <- nl %*% t(nl) / n.mc
  par.df <- eigen(J, TRUE, only.values = TRUE)$values
  par.df <- par.df[par.df > 0]
  par.df <- sum((par.df / par.df[1] > 1e-10) |
                (cumsum(par.df) / sum(par.df) < 0.90))

  ### Expected noncenteriality
  M0 <- nrow(pl0)
  Ma <- nrow(pla)
  if(given == "0"){
    id <- M0 + (1:Ma)
  } else{
    id <- 1:M0
  }
  J.nc <- matrix(J[id, id], nrow = length(id))
  nu <- matrix(mu[id], nrow = length(id))
  ### Numerical unstable.
  # par.nc <- t(nu) %*% solve(J.nc) %*% nu
  tmp <- eigen(J.nc, TRUE)
  tmp.nc <- tmp$values[tmp$values > 0]
  tmp.nc <- sum((tmp.nc / tmp.nc[1] > 1e-8) |
                (cumsum(tmp.nc) / sum(tmp.nc) < 0.90))
  nu <- t(nu) %*% tmp$vectors[, 1:tmp.nc]
  par.nc <- nu %*% diag(1 / tmp$values[1:tmp.nc], tmp.nc, tmp.nc) %*% t(nu)

  ### For returns.
  ret <- c(par.df, par.nc)
  if(verbose){
    cat("K.0=", K.0, ", M0=", M0, " v.s. K.a=", K.a, ", Ma=", Ma,
        " | given=", given, " : df=", ret[1], ", nc=", ret[2],
        ".\n", sep = "")
  }
  ret
} # End of get.E.chi2.pv().


get.E.delta.pv <- function(x, fcobj.0, fcobj.a, tau = 0.5, n.mc = 1000){
  N <- nrow(x$X.gbd)

  ETA.0 <- fcobj.0$param$ETA
  BETA.0 <- fcobj.0$param$BETA
  MU.0 <- fcobj.0$param$MU
  SIGMA.0 <- fcobj.0$param$SIGMA

  ETA.a <- fcobj.a$param$ETA
  BETA.a <- fcobj.a$param$BETA
  MU.a <- fcobj.a$param$MU
  SIGMA.a <- fcobj.a$param$SIGMA

  n.mc.0 <- rbinom(1, n.mc, c(tau, 1-tau))
  n.mc.a <- n.mc - n.mc.0

  E.delta <- NULL
  if(n.mc.0 > 0){
    tmp <- lapply(1:n.mc.0, function(i){
      x.new <- GenDataSet.pv(N, ETA.0, BETA.0, MU.0, SIGMA.0, x$X.gbd)$x
      logL.I(x.new, fcobj.a) - logL.I(x.new, fcobj.0)
    })
    E.delta <- c(E.delta, tmp)
  }
  if(n.mc.a > 0){
    tmp <- lapply(1:n.mc.a, function(i){
      x.new <- GenDataSet.pv(N, ETA.a, BETA.a, MU.a, SIGMA.a, x$X.gbd)$x
      logL.I(x.new, fcobj.a) - logL.I(x.new, fcobj.0)
    })
    E.delta <- c(E.delta, tmp)
  }

  do.call("sum", E.delta) / n.mc
} # End of get.E.delta.pv().
