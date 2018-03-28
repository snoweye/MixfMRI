### Perform H0: alpha_1/(alpha_1 + beta_1) == alpha_2/(alpha_2 + beta_2) vs
###         Ha: alpha_1/(alpha_1 + beta_1) != alpha_2/(alpha_2 + beta_2)

lrt.betamean <- function(PV.gbd, CLASS.gbd, K, proc = c("1", "2")){
  K.CLASS <- as.integer(max.gbd(CLASS.gbd))
  if(K.CLASS > K){
    stop("CLASS.gbd and K are not matched.")
  }

  ### For constrOptim() and H0.
  fn.H0 <- function(theta, x.1, x.2){
    -sum(dbeta(x.1, theta[1], theta[2], log = TRUE)) -
    sum(dbeta(x.2, theta[1] * theta[3], theta[2] * theta[3], log = TRUE))
  }

  ### For constrOptim() and Ha.
  fn.Ha <- function(theta, x.1, x.2){
    -sum(dbeta(x.1, theta[1], theta[2], log = TRUE)) -
    sum(dbeta(x.2, theta[3], theta[4], log = TRUE))
  }

  ### For initial Beta parameters.
  fn.init <- function(x.gbd){
    x.bar <- mean(x.gbd)
    v.bar <- var(x.gbd) 
    a.hat <- x.bar * (x.bar * (1 - x.bar) / v.bar - 1)
    b.hat <- (1 - x.bar) * (x.bar * (1 - x.bar) / v.bar - 1)
    ret <- c(a.hat, b.hat)
    if(v.bar > x.bar * (1 - x.bar))
      ret <- -ret
    ret
  }

  ### Run for all (K-1)*K/2 pairs.
  ret <- NULL
  for(i.k.1 in 2:(K - 1)){
    N.class.gbd.1 <- sum(CLASS.gbd == i.k.1)

    if(N.class.gbd.1 > 1){    # in case of empty cluster.
      x.1.gbd <- PV.gbd[CLASS.gbd == i.k.1]

      for(i.k.2 in (i.k.1 + 1):K){
        N.class.gbd.2 <- sum(CLASS.gbd == i.k.2)

        if(N.class.gbd.2 > 1){    # in case of empty cluster.
          x.2.gbd <- PV.gbd[CLASS.gbd == i.k.2]

          ### logL under H0.
          x.all <- c(x.1.gbd, x.2.gbd)
          theta <- fn.init(x.all)
          theta <- c(theta, 1.0)
          tmp <- constrOptim(theta, fn.H0, NULL, diag(3), rep(1e-16, 3),
                             method = "Nelder-Mead", x.1 = x.1.gbd, x.2 = x.2.gbd)
          H0.alpha <- tmp$par[1]
          H0.beta <- tmp$par[2]
          H0.scale <- tmp$par[3]
          logL.0 <- -tmp$value

          ### logL under Ha.
          theta <- c(fn.init(x.1.gbd), fn.init(x.2.gbd))
          tmp <- constrOptim(theta, fn.Ha, NULL, diag(4), rep(1e-16, 4),
                             method = "Nelder-Mead", x.1 = x.1.gbd, x.2 = x.2.gbd)
          Ha.alpha.1 <- tmp$par[1]
          Ha.beta.1 <- tmp$par[2]
          Ha.alpha.2 <- tmp$par[3]
          Ha.beta.2 <- tmp$par[4]
          logL.a <- -tmp$value

          ### LRT.
          lrt <- -2 * (logL.0 - logL.a)
          if(lrt < 0){
            lrt <- 0
          }

          ### p-value.
          p.value <- pchisq(lrt, 1, lower.tail = FALSE)
          ret <- rbind(ret, c(i.k.1, i.k.2, logL.0, logL.a, lrt, p.value))
        }
      }
    }
  }


  if(proc[1] == "1"){
    ret <- cbind(ret, qvalue(ret[, ncol(ret)]))
  } else if(proc[1] == "2"){
    qvalue.p2 <- fdr.bh.p2(ret[, ncol(ret)], q = 0.05)$adjp
    ret <- cbind(ret, qvalue.p2)
  } else{
    stop("proc is not found.")
  }
  colnames(ret) <- c("i.k.1", "i.k.2",
                     "logL.0", "logL.a", "lrt", "p.value", "q.value")

  ret
} # End of lrt.betamean().

