
mle.beta.constr <- function(PARAM, i.k, method = .FC.CT$optim.method[1]){
  beta.theta <- PARAM$BETA[[i.k]]
  if(method[1] == "BFGS"){
    ret <- constrOptim(theta = beta.theta,
                       f = beta.constrOptim.f,
                       grad = beta.constrOptim.grad,
                       ui = rbind(c(1, 0),
                                  c(-1, 0),
                                  c(0, 1)),
                       ci = c(.FC.CT$INIT$BETA.alpha.min,
                              -.FC.CT$INIT$BETA.alpha.max,
                              .FC.CT$INIT$BETA.beta.min),
                       method = "BFGS",
                       i.k = i.k)
    beta.theta <- ret$par
  }

  ret <- constrOptim(theta = beta.theta,
                     f = beta.constrOptim.f,
                     grad = beta.constrOptim.grad,
                     ui = rbind(c(1, 0),
                                c(-1, 0),
                                c(0, 1)),
                     ci = c(.FC.CT$INIT$BETA.alpha.min,
                            -.FC.CT$INIT$BETA.alpha.max,
                            .FC.CT$INIT$BETA.beta.min),
                     method = "Nelder-Mead",
                     i.k = i.k)

  if(.MixfMRIEnv$CONTROL$debug > 2){
    .MixfMRIEnv$cat("  mle.beta.constr:",
                      "\n    par: ", paste(formatC(ret$par, 6), collapse = " "),
                      ", value: ", ret$value,
                      ", counts: ", paste(ret$counts, collapse = " "),
                      ", conv: ", ret$convergence, "\n",
                      sep = "", quiet = TRUE)
  }
  
  ret$par
} # End of mle.beta.constr().

beta.constrOptim.f <- function(x, i.k){
  tmp <- dbeta(.MixfMRIEnv$PV.gbd, x[1], x[2], log = TRUE)

  ret <- -sum_gbd(tmp * .MixfMRIEnv$Z.gbd[, i.k])
  ret
} # End of beta.constrOptim.f().

# d/dx B(x, y) = B(x, y) * (Psi(x) - Psi(x + y))
# where Psi(.) is a digamma function.
beta.constrOptim.grad <- function(x, i.k){
  z <- sum_gbd(.MixfMRIEnv$Z.gbd[, i.k])
  A <- digamma(x[1]) * z
  B <- digamma(x[2]) * z
  AB <- digamma(x[1] + x[2]) * z

  tmp.PV.Z <- sum_gbd(.MixfMRIEnv$log.PV.gbd * .MixfMRIEnv$Z.gbd[, i.k])
  tmp.1.PV.Z <- sum_gbd(.MixfMRIEnv$log.1.PV.gbd * .MixfMRIEnv$Z.gbd[, i.k])

  ret <- -c(AB - A + tmp.PV.Z, AB - B + tmp.1.PV.Z)
  ret
} # End of beta.constrOptim.grad().

