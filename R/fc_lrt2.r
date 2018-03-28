### Perform H0: alpha/(alpha + beta) >= 0.05 vs Ha: alpha/(alpha + beta) < 0.05

lrt2 <- function(PV.gbd, CLASS.gbd, K, H0.mean = .FC.CT$LRT$H0.mean,
    upper.beta = .FC.CT$INIT$BETA.beta.max,
    proc = c("1", "2", "weight")){
  if((H0.mean <= 0) || (H0.mean >= 1.0)){
    stop("It should be 0 < H0.mean < 1.0.")
  }

  K.CLASS <- as.integer(max.gbd(CLASS.gbd))
  if(K.CLASS > K){
    stop("CLASS.gbd and K are not matched.")
  }
  beta.scale <- H0.mean / (1 - H0.mean)

  ### For optim() and H0.
  fn.H0 <- function(theta, x.gbd){
    -sum.gbd(dbeta(x.gbd, theta[1], theta[2], log = TRUE))
  }
  ui.H0 <- rbind(c(1, 0), c(1 - H0.mean, -H0.mean), c(0, 1))
  ci.H0 <- c(0, 0, 0)

  ### For constrOptim() and Ha.
  fn.Ha <- function(theta, x.gbd){
    -sum.gbd(dbeta(x.gbd, theta[1], theta[2], log = TRUE))
  }
  ui.Ha <- rbind(c(1, 0), c(0, 1))
  ci.Ha <- c(0, 0)

  ret <- NULL
  for(i.k in 1:K){
    N.class.gbd <- sum.gbd(CLASS.gbd == i.k)

    if(N.class.gbd > 0){    # in case of empty cluster.
      tmp.gbd <- PV.gbd[CLASS.gbd == i.k]

      ### logL under H0.
      tmp <- constrOptim(c(beta.scale * 1.02, 1.01), fn.H0, NULL, ui.H0, ci.H0,
                         method = "Nelder-Mead", x.gbd = tmp.gbd)
      H0.alpha <- tmp$par[1]
      H0.beta <- tmp$par[2]
      logL.0 <- -tmp$value

      ### logL under Ha.
      tmp <- constrOptim(c(beta.scale * 0.99, 1.01), fn.Ha, NULL, ui.Ha, ci.Ha,
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
      ret <- rbind(ret, c(i.k, H0.alpha, H0.beta, Ha.alpha, Ha.beta,
                          logL.0, logL.a, lrt, p.value))
    } else{
      ret <- rbind(ret, c(i.k, rep(NA, 8)))
    }
  }

  if(proc[1] == "1"){
    ret <- cbind(ret, qvalue(ret[, ncol(ret)]))
  } else if(proc[1] == "2"){
    qvalue.p2 <- fdr.bh.p2(ret[, ncol(ret)], q = 0.05)$adjp
    ret <- cbind(ret, qvalue.p2)
  } else if(proc[1] == "weight"){
    w.size <- rep(0, K)
    for(i.k in 1:K){
      w.size[i.k] <- sum.gbd(CLASS.gbd == i.k)
    }
    w.size <- w.size / sum(w.size) * K 
    qvalue.w <- fdr.bh.p2(ret[, ncol(ret)], w = w.size, q = 0.05)$adjp
    ret <- cbind(ret, qvalue.w)
  } else{
    stop("proc is not found.")
  }
  colnames(ret) <- c("i.k", "mle.0.alpha", "mle.0.beta",
                     "mle.a.alpha", "mle.a.beta",
                     "logL.0", "logL.a", "lrt", "p.value", "q.value")

  ret
} # End of lrt2().
