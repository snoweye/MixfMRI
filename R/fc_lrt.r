### Perform H0: alpha = 1 & beta = 1 v.s. Ha: alpha < 1 & beta > 1.

lrt <- function(PV.gbd, CLASS.gbd, K, H0.alpha = .FC.CT$LRT$H0.alpha,
    H0.beta = .FC.CT$LRT$H0.beta){
  if((H0.alpha <= 0) || (H0.alpha > H0.beta)){
    stop("It should be 0 < H0.alpha <= H0.beta.")
  }

  K.CLASS <- as.integer(max.gbd(CLASS.gbd))
  if(K.CLASS > K){
    stop("CLASS.gbd and K are not matched.")
  }

  ### For constrOptim().
  fn <- function(theta, x.gbd){
    -sum.gbd(dbeta(x.gbd, theta[1], theta[2], log = TRUE))
  }
  ui <- rbind(c(1, 0), c(-1, 0), c(0, 1))
  ci <- c(0, -H0.alpha, H0.beta)

  ret <- NULL
  for(i.k in 1:K){
    N.class.gbd <- sum.gbd(CLASS.gbd == i.k)

    if(N.class.gbd > 0){    # in case of empty cluster.
      tmp.gbd <- PV.gbd[CLASS.gbd == i.k]

      ### logL under H0.
      logL.0 <- sum.gbd(dbeta(tmp.gbd, H0.alpha, H0.beta, log = TRUE))

      ### logL under Ha.
      tmp <- constrOptim(c(H0.alpha - 0.01, H0.beta + 0.01), fn, NULL, ui, ci,
                         method = "Nelder-Mead", x.gbd = tmp.gbd)
      Ha.alpha <- tmp$par[1] 
      Ha.beta <- tmp$par[2] 
      logL.a <- -tmp$value

      ### LRT.
      lrt <- -2 * (logL.0 - logL.a)
      if(lrt < 0){
        lrt <- 0
      }

      ### p-value.
      p.value <- pchisq(lrt, 2, lower.tail = FALSE)
      ret <- rbind(ret, c(i.k, Ha.alpha, Ha.beta,
                          logL.0, logL.a, lrt, p.value))
    } else{
      ret <- rbind(ret, c(i.k, rep(NA, 6)))
    }
  }
  ret <- cbind(ret, qvalue(ret[, ncol(ret)]))
  colnames(ret) <- c("i.k", "mle.alpha", "mle.beta",
                     "logL.0", "logL.a", "lrt", "p.value", "q.value")

  ret
} # End of lrt().
