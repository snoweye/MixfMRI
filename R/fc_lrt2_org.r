### Perform H0: alpha/(alpha + beta) = 0.05 vs. Ha: alpha/(alpha + beta) < 0.05.

lrt2.org <- function(PV.gbd, CLASS.gbd, K, H0.mean = 0.05,
    upper.beta = .FC.CT$INIT$BETA.beta.max){
  if((H0.mean <= 0) || (H0.mean >= 0.5)){
    stop("It should be 0 < H0.mean < 0.5.")
  }

  K.CLASS <- as.integer(max.gbd(CLASS.gbd))
  if(K.CLASS > K){
    stop("CLASS.gbd and K are not matched.")
  }
  beta.scale <- H0.mean / (1 - H0.mean)

  ### For optim() and H0.
  fn.H0 <- function(theta, x.gbd){
    -sum_gbd(dbeta(x.gbd, beta.scale * theta[1], theta[1], log = TRUE))
  }

  ### For constrOptim() and Ha.
  fn.Ha <- function(theta, x.gbd){
    -sum_gbd(dbeta(x.gbd, theta[1], theta[2], log = TRUE))
  }
  ui.Ha <- rbind(c(1, 0), c(H0.mean - 1, H0.mean))
  ci.Ha <- c(0, 0)

  ret <- NULL
  for(i.k in 1:K){
    N.class.gbd <- sum_gbd(CLASS.gbd == i.k)

    if(N.class.gbd > 0){    # in case of empty cluster.
      tmp.gbd <- PV.gbd[CLASS.gbd == i.k]

      ### logL under H0.
      tmp <- optim(c(1.01), fn.H0, NULL, x.gbd = tmp.gbd,
                   method = "Brent", lower = 0, upper = upper.beta)
      H0.beta <- tmp$par[1] 
      H0.alpha <- beta.scale * H0.beta
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
      p.value <- pchisq(lrt, 1, lower.tail = FALSE)
      ret <- rbind(ret, c(i.k, H0.beta, Ha.alpha, Ha.beta,
                          logL.0, logL.a, lrt, p.value))
    } else{
      ret <- rbind(ret, c(i.k, rep(NA, 7)))
    }
  }
  colnames(ret) <- c("i.k", "mle.0.beta",
                     "mle.a.alpha", "mle.a.beta",
                     "logL.0", "logL.a", "lrt", "p.value")

  ret
} # End of lrt2.org().
