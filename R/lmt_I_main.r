### For LMT.
### Copy from EMCluster.
### For identity cov matrix only.

lmt.I <- function(fcobj.0, fcobj.a, X.gbd, PV.gbd, tau = 0.5,
    n.mc.E.delta = 1000, n.mc.E.chi2 = 1000, verbose = FALSE){
  if(class(fcobj.0) != "fclust" || class(fcobj.a) != "fclust"){
    stop("fcobj.0 and fcobj.a should be both in \"fclust\" class.")
  }
  if(fcobj.0$param$K == fcobj.a$param$K){
    stop("fcobj.0 and fcobj.a should have different numbers of clusters.")
  }

  ### Reset global
  set.global(X.gbd, PV.gbd)

  ### x
  x <- list(X.gbd = X.gbd, PV.gbd = PV.gbd)

  ### K
  k.0 <- fcobj.0$param$K
  k.a <- fcobj.a$param$K

  ### logL
  ll.0 <- fcobj.0$param$logL
  ll.a <- fcobj.a$param$logL

  ### Likelihood ratio statistics
  delta.hat <- fcobj.a$param$logL - fcobj.0$param$logL
  E.delta <- get.E.delta.I(x, fcobj.0, fcobj.a, tau = tau, n.mc = n.mc.E.delta)

  ### Chi-squared statistics
  E.chi2.0 <- get.E.chi2.I(x, fcobj.0, fcobj.a, "0", tau = tau,
                           n.mc = n.mc.E.chi2, verbose = verbose)
  E.chi2.a <- get.E.chi2.I(x, fcobj.0, fcobj.a, "a", tau = tau,
                           n.mc = n.mc.E.chi2, verbose = verbose)

  ### Testing statistics.
  T <- 2 * (delta.hat - E.delta)

  ### p-values.
  pv.0 <- pchisq.my(T, E.chi2.0[1], E.chi2.0[2], lower.tail = FALSE)
  pv.a <- pchisq.my(T, E.chi2.a[1], E.chi2.a[2], lower.tail = FALSE)
  pv <- pv.0 * tau + pv.a * (1 - tau)

  ### Return.
  ret <- list(K.0 = k.0, K.a = k.a, ll.0 = ll.0, ll.a = ll.a,
              delta.hat = delta.hat,
              E.delta = E.delta, E.chi2.0 = E.chi2.0, E.chi2.a = E.chi2.a,
              T = T, pv.0 = pv.0, pv.a = pv.a, pv = pv)
  class(ret) <- "lmt.I"
  ret
} # End of lmt.I().

print.lmt.I <- function(x, digits = max(4, getOption("digits") - 3), ...){
  K.0 <- x$K.0
  K.a <- x$K.a
  ll.0 <- x$ll.0
  ll.a <- x$ll.a
  delta.hat <- x$delta.hat
  E.delta <- x$E.delta
  E.chi2.0 <- x$E.chi2.0
  E.chi2.a <- x$E.chi2.a
  T <- x$T
  pv.0 <- x$pv.0
  pv.a <- x$pv.a
  pv <- x$pv

  cat("- H.0: K = ", K.0, " vs H.a: K = ", K.a, "\n",
      "    ll.0 = ", ll.0, ", ll.a = ", ll.a, "\n",
      "    df.0 = ", E.chi2.0[1], ", nc.0 = ", E.chi2.0[2],
      ", df.a = ", E.chi2.a[1], ", nc.a = ", E.chi2.a[2], "\n",
      "    delta.hat = ", delta.hat, ", E.delta = ", E.delta,
      ", T = ", T, "\n",
      "    pv.0 = ", pv.0, ", pv.a = ", pv.a, " pv = ", pv,
      "\n", sep = "")

  invisible()
} # End of print.lmt.I().
