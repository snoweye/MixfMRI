### This file is for Cov(z_nk) and Cov(z_nk - z_n1).

### Partial logL. Return a matrix with dimension M * N.
# partial.logL.I <- function(x, ETA, BETA, MU, SIGMA, post.z)

### Compute posterior. Return a matrix with dimension N * K.
# post.prob <- function(x, fcobj)

### Partial all z_nk. Return a list of length = length(index.N).
partial.post.z <- function(x, ETA, BETA, MU, SIGMA, post.z, index.N){
  BETA.digamma <- lapply(BETA, function(b){ digamma(b) - digamma(sum(b)) })
  SIGMA.inv <- lapply(SIGMA, function(s){ 1 / s$full })

  ret <- list()
  for(i.n in 1:length(index.N)){
    ret[[i.n]] <- partial.post.z.1(x, ETA, BETA, MU, SIGMA, post.z,
                                   BETA.digamma, SIGMA.inv, index.N[i.n])
  }

  ret
} # End of partial.post.z().

### Partial one z_nk. Return a matrix with dimension K * M for index i.n.
partial.post.z.1 <- function(x, ETA, BETA, MU, SIGMA, post.z,
    BETA.digamma, SIGMA.inv, i.n){
  K <- length(ETA)

  ### Make a K * K weight matrix.
  Zbuf.w <- list()
  for(i.k in 1:K){
    tmp <- rep(-post.z[i.n, i.k], K)
    tmp[i.k] <- tmp[i.k] + 1
    Zbuf.w[[i.k]] <- post.z[i.n,] / ETA * tmp
  }

  ### For ETA.
  Zbuf.ETA <- list()
  if(K > 1){
    for(i.k in 1:(K - 1)){
      Zbuf.ETA[[i.k]] <- Zbuf.w[[i.k]][-K] +
                         post.z[i.n, K] / ETA[K] * post.z[i.n, i.k]
    }
    Zbuf.ETA[[K]] <- -colSums(do.call("rbind", Zbuf.ETA))
  } else{
    Zbuf.ETA[[1]] <- vector(mode = "numeric", length = 0)
  }

  ### For BETA.
  Zbuf.BETA <- list()
  for(i.k in 1:K){
    tmp <- list()
    for(j.k in 2:K){
      tmp[[j.k]] <- Zbuf.w[[i.k]][j.k] *
                    partial.f.beta(x, BETA, BETA.digamma, i.n, j.k)
    }
    Zbuf.BETA[[i.k]] <- do.call("c", tmp)
  }
        
  ### For MU.
  Zbuf.MU <- list()
  Zbuf.SIGMA <- list()
  for(i.k in 1:K){
    tmp.mu <- list()
    tmp.sigma <- list()
    for(j.k in 1:K){
      tmp <- partial.f.mu.sigma(x, MU, SIGMA, SIGMA.inv, i.n, j.k)
      tmp.mu[[j.k]] <- Zbuf.w[[i.k]][j.k] * tmp$mu
      tmp.sigma[[j.k]] <- Zbuf.w[[i.k]][j.k] * tmp$sigma
    }
    Zbuf.MU[[i.k]] <- do.call("c", tmp.mu)
    Zbuf.SIGMA[[i.k]] <- do.call("c", tmp.sigma)
  }

  ### For Zbuf. Combine all ETA, BETA, MU, and SIGMA.
  Zbuf <- list()
  for(i.k in 1:K){
    Zbuf[[i.k]] <- c(Zbuf.ETA[[i.k]], Zbuf.BETA[[i.k]], Zbuf.MU[[i.k]],
                     Zbuf.SIGMA[[i.k]])
  }
  Zbuf <- do.call("rbind", Zbuf)

  Zbuf
} # End of partial.post.z.1().

partial.f.beta <- function(x, BETA, BETA.digamma, i.n, j.k){
  dbeta(x$PV.gbd[i.n], BETA[[j.k]][1], BETA[[j.k]][2]) *
  ((BETA[[j.k]] - 1) / x$PV.gbd[i.n] - BETA.digamma[[j.k]])
} # End of partial.f.beta().

partial.f.mu.sigma <- function(x, MU, SIGMA, SIGMA.inv, i.n, j.k){
  f <- dnorm(x$X.gbd[i.n,], MU[, j.k], sd = sqrt(SIGMA[[j.k]]$full))
  tmp.mu <- (x$X.gbd[i.n,] - MU[, j.k]) * SIGMA.inv[[j.k]]
  tmp.sigma <- (tmp.mu^2 - SIGMA.inv[[j.k]]) / 2
  list(mu = f * tmp.mu, sigma = f * tmp.sigma)
} # End of partial.f.mu.sigma().

cov.param <- function(x, fcobj, post.z, drop.ETA1 = FALSE){
  K <- fcobj$param$K
  ETA <- fcobj$param$ETA
  BETA <- fcobj$param$BETA
  MU <- fcobj$param$MU
  SIGMA <- fcobj$param$SIGMA

  ### cov of param = {ETA, BETA, MU, SIGMA}.
  pl <- partial.logL.I(x, ETA, BETA, MU, SIGMA, post.z)
  ### Drop eta_1 if on the boundary.
  if(fcobj$param$ETA[1] == fcobj$param$min.1st.prop && drop.ETA1){
    pl <- pl[-1,]
  }

  ### Get cov of param
  I <- pl %*% t(pl)
  cov <- ginv(I)

  list(I = I, cov = cov)
} # End of cov.param().

cov.post.z <- function(x, fcobj, post.z, cov.param = NULL, all.x = FALSE,
    drop.ETA1 = FALSE){
  K <- fcobj$param$K
  ETA <- fcobj$param$ETA
  BETA <- fcobj$param$BETA
  MU <- fcobj$param$MU
  SIGMA <- fcobj$param$SIGMA

  ### nabla of post.z.
  if(all.x){
    index.N <- 1:nrow(x$X.gbd)
  } else{
    index.N <- which(fcobj$class != 1)
  }
  nabla.post.z <- partial.post.z(x, ETA, BETA, MU, SIGMA, post.z, index.N)
  ### Drop eta_1 if on the boundary.
  if(fcobj$param$ETA[1] == fcobj$param$min.1st.prop && drop.ETA1){
    nabla.post.z <- lapply(nabla.post.z, function(x) x[, -1])
  }

  ### Get cov of post.z
  if(is.null(cov.param)){
    cov.param <- cov.param(x, fcobj, drop.ETA1 = drop.ETA1)$cov
  }
  cov.post.z <- lapply(nabla.post.z, function(x) x %*% cov.param %*% t(x))
  cov.post.z
} # End of cov.post.z().

cov.logit.z <- function(x, fcobj, post.z, cov.param = NULL,
    cov.post.z = NULL, all.x = FALSE, drop.ETA1 = FALSE){
  K <- fcobj$param$K

  ### nabla of logit.z.
  if(all.x){
    index.N <- 1:nrow(x$X.gbd)
  } else{
    index.N <- which(fcobj$class != 1)
  }
  nabla.logit.z <- lapply(index.N, function(i){ partial.logit.p(post.z[i,]) })

  ### Get cov of logit.z
  if(is.null(cov.param)){
    cov.param <- cov.param(x, fcobj, post.z, drop.ETA1 = drop.ETA1)$cov
  }
  if(is.null(cov.post.z)){
    cov.post.z <- cov.post.z(x, fcobj, post.z, cov.param = cov.param,
                                 all.x = all.x, drop.ETA1 = drop.ETA1)
  }
  if(length(nabla.logit.z) != length(cov.post.z)){
    stop("nabla.logit.z and cov.post.z are not of the same length.")
  }
  cov.logit.z <- lapply(1:length(nabla.logit.z),
                        function(i){
                          nabla.logit.z[[i]] %*%
                          cov.post.z[[i]][1:(K-1), 1:(K-1)] %*%
                          t(nabla.logit.z[[i]])
                        })
  cov.logit.z
} # End of cov.logit.z().

logor.stat <- function(x, fcobj, post.z, cov.param = NULL,
    cov.post.z = NULL, cov.logit.z = NULL, all.x = FALSE, drop.ETA1 = FALSE){
  K <- fcobj$param$K

  ### Get cov of logor
  if(is.null(cov.param)){
    cov.param <- cov.param(x, fcobj, post.z, drop.ETA1 = drop.ETA1)$cov
  }
  if(is.null(cov.post.z)){
    cov.post.z <- cov.post.z(x, fcobj, post.z, cov.param = cov.param,
                                 all.x = all.x, drop.ETA1 = drop.ETA1)
  }
  if(is.null(cov.logit.z)){
    cov.logit.z <- cov.logit.z(x, fcobj, post.z, cov.param = cov.param,
                                   cov.post.z = cov.post.z, all.x = all.x,
                                   drop.ETA1 = drop.ETA1)
  }

  ### Get log odds ratio and its cov matrix.
  if(all.x){
    index.N <- 1:nrow(x$X.gbd)
  } else{
    index.N <- which(fcobj$class != 1)
  }
  logit.p <- log(post.z[index.N,] / (1 - post.z[index.N,]))
  A <- cbind(rep(-1, K - 1), diag(1, K - 1))
  log.or <- logit.p %*% t(A)
  cov.log.or <- lapply(cov.logit.z, function(x) A %*% x %*% t(A))

  ### Get statistics and p-values.
  df <- rep(NA, nrow(log.or))
  test.stat <- rep(NA, nrow(log.or))
  for(i in 1:nrow(log.or)){
    ts <- matrix(log.or[i,] - 1, nrow = 1)    # check if log.or == 1
    if(all(is.finite(ts)) && all(is.finite(cov.log.or[[i]]))){
      # df[i] <- sum(eigen(cov.log.or[[i]], only.values = TRUE)$values > 0)
      df[i] <- as.integer(rankMatrix(cov.log.or[[i]]))
      test.stat[i] <- ts %*% ginv(cov.log.or[[i]]) %*% t(ts)
    }
  }

  ### Return.
  list(log.or = log.or, cov.log.or = cov.log.or, df = df,
       test.stat = test.stat)
} # End of logor.stat().

