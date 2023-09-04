### This file contains files for estimating parameters emperically.
# x <- rbeta(1000, 0.5, 2)
# (m <- mean(x))
# (v <- var(x))
# (alpha <- 1/v * (m^2 - m^3 - m * v))
# (beta <- (1 - m) / m * alpha)
# alpha / (alpha + beta)
# alpha * beta / (alpha + beta)^2 / (alpha + beta + 1)
# M1 <- mean(x)
# M2 <- mean(x^2)
# (alpha <- M1 * (M1 - M2) / (M2 - M1^2))
# (beta <- (1 - M1) * (M1 - M2) / (M2 - M1^2))


### Estimate BETA based on classified data.
initial.BETA <- function(CLASS.gbd, K){
  BETA <- list()

  ### For the first cluster.
  BETA[[1]] <- c(1.0, 1.0)

  ### For the rest cluster.
  for(i.k in 2:K){
    tmp <- .MixfMRIEnv$PV.gbd[CLASS.gbd == i.k]
    tmp.M <- c(sum_gbd(tmp), sum_gbd(tmp^2))

    ### First two moments, M1 and M2.
    tmp.N <- length_gbd(tmp)
    tmp.M <- tmp.M / tmp.N

    ### Get alpha and beta.
    alpha <- tmp.M[1] * (tmp.M[1] - tmp.M[2]) / (tmp.M[2] - tmp.M[1]^2)
    beta <- (1 - tmp.M[1]) * (tmp.M[1] - tmp.M[2]) / (tmp.M[2] - tmp.M[1]^2)

    if(alpha > .FC.CT$INIT$BETA.alpha.max ||
       alpha < .FC.CT$INIT$BETA.alpha.min){
      alpha <- runif.bcast(1, min = .FC.CT$INIT$BETA.alpha.min,
                              max = .FC.CT$INIT$BETA.alpha.max)
    }
    if(beta < .FC.CT$INIT$BETA.beta.min){
      beta <- rexp.bcast(1) + .FC.CT$INIT$BETA.beta.min
    }

    BETA[[i.k]] <- c(alpha, beta)
  }

  BETA
} # End of initial.BETA().


### Estimate MU based on classified data.
initial.MU <- function(CLASS.gbd, K){
  p.X <- ncol(.MixfMRIEnv$X.gbd)
  MU.CLASS <- matrix(0, nrow = p.X, ncol = K)

  for(i.k in 1:K){
    tmp.id <- which(CLASS.gbd == i.k)
    tmp.MU <- colMeans.gbd(matrix(.MixfMRIEnv$X.gbd[tmp.id,], ncol = p.X))
    MU.CLASS[, i.k] <- tmp.MU
  }

  MU.CLASS
} # End of initial.MU().


### Estimate SIGMA based on classified data.
initial.SIGMA.I <- function(MU, CLASS.gbd, K){
  p.X <- ncol(.MixfMRIEnv$X.gbd)
  SIGMA.CLASS <- list()

  for(i.k in 1:K){  
    tmp.id <- which(CLASS.gbd == i.k)
    tmp.N <- length_gbd(tmp.id)
    B <- W.plus.y(.MixfMRIEnv$X.gbd[tmp.id,], -MU[, i.k],
                  length(tmp.id), p.X) /
         sqrt(tmp.N - 1)
    tmp.full <- colSums.gbd(B^2)
    tmp.U <- decompsigma.I(tmp.full)
    SIGMA.CLASS[[i.k]] <- list(full = tmp.full, U = tmp.U$value,
                               U.check = tmp.U$check)
  }

  SIGMA.CLASS
} # End of initial.SIGMA.I().

initial.SIGMA.V <- function(MU, CLASS.gbd, K){
  p.X <- ncol(.MixfMRIEnv$X.gbd)
  SIGMA.CLASS <- list()

  for(i.k in 1:K){  
    tmp.id <- which(CLASS.gbd == i.k)
    tmp.N <- length_gbd(tmp.id)
    B <- W.plus.y(.MixfMRIEnv$X.gbd[tmp.id,], -MU[, i.k],
                  length(tmp.id), p.X) /
         sqrt(tmp.N - 1)
    tmp.full <- crossprod.gbd(B)
    tmp.U <- decompsigma.V(tmp.full)
    SIGMA.CLASS[[i.k]] <- list(full = tmp.full, U = tmp.U$value,
                               U.check = tmp.U$check)
  }

  SIGMA.CLASS
} # End of initial.SIGMA.V().


### This function collects N.CLASS
get.N.CLASS <- function(K){
  tabulate.gbd(.MixfMRIEnv$CLASS.gbd, nbins = K)
} # End of get.N.CLASS().

