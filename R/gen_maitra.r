
fbarmax <- function(mu, barom, maxom,
    eta = c(0.97675624, 0.01917240, 0.00407136)){
  # require(MixSim)
  sigma <- array(c(1, 1, 1), dim = c(1, 1, 3))
  # eta <- c(0.97675624, 0.01917240, 0.00407136)
  mu <- matrix(c(0, mu[1], mu[2]), nrow = 3, ncol = 1)
  ommap <- overlap(Pi = eta, Mu = mu, S = sigma)
  ((barom  - ommap$BarOmega) / barom)^2 + ((maxom - ommap$MaxOmega) / maxom)^2
} # End of fbarmax().

gendataset <- function(phantom, overlap, smooth = FALSE){
  eta <- table(phantom) / sum(table(phantom))
  mu <- c(0, optim(runif(2) + 3, fn = fbarmax,
                   barom = overlap, maxom = overlap, eta = eta)$par)
  tval <- matrix(rnorm(256^2), ncol = 256) + mu[phantom+1]

  if(smooth == TRUE){
    ystat <- tval
    ystat[is.na(ystat)] <- 0
    ## smooth by Garcia method
    ystat <-  gcv.smooth2d(ystat, interval = c(0,10))$im.smooth
    ystat[is.na(tval)] <- NA
    tval <- ystat
  }

  pval <- pnorm(tval, lower.tail = F)
  ret <- list(eta = eta, overlap = overlap,
              mu = mu, class.id = phantom + 1, tval = tval, pval = pval)
  ret
} # End of gendataset().

dpval <- function(x, mu = 0, log = FALSE){
  f <- function(p, mu = 0){
    dnorm(x = qnorm(p = 1 - p) - mu) / dnorm(x = qnorm(p = 1 - p))
  }
  f.log <- function(p, mu = 0){
    dnorm(x = qnorm(p = 1 - p) - mu, log = TRUE) -
    dnorm(x = qnorm(p = 1 - p), log = TRUE)
  }

  if(log){
    ret <- sapply(x, f.log, mu = mu)
  } else{
    ret <- sapply(x, f, mu = mu)
  }
  ret
} # End of dpval().

dmixpval <- function(x, eta, mu){
  K <- length(eta)
  if(K != length(mu)){
    stop("length(eta) != length(mu)")
  }

  ret <- 0
  for(i.k in 1:K){
    ret <- ret + exp(dpval(x, mu = mu[i.k], log = TRUE) + log(eta[i.k]))
  }

  ret
} # End of dmixpval().
