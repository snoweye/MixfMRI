library(MixfMRI, quietly = TRUE)
set.seed(1234)

N <- 1000
K <- 3
p.X <- 3
eta <- c(0.60, 0.30, 0.10)
beta <- matrix(c(1, 1, 0.5, 6, 0.8, 4), ncol = K)          # 2 * K
mu <- matrix(c(0, 0, 0, -2, -2, -2, 4, 4, 4), ncol = K)    # p.X * K
std <- matrix(c(1, 1, 1, 2, 1, 2, 1, 2, 1), ncol = K)      # p.X * K

### For toy1
N.K <- tabulate(sample(1:K, N, prob = eta, replace = TRUE))
CLASS.gbd <- NULL
PV.gbd <- NULL
X.gbd <- NULL
for(i.k in 1:K){
  CLASS.gbd <- c(CLASS.gbd, rep(i.k, N.K[i.k]))
  PV.gbd <- c(PV.gbd, rbeta(N.K[i.k], beta[1, i.k], beta[2, i.k]))
  tmp <- NULL
  for(i.p in 1:p.X){
    tmp <- cbind(tmp, rnorm(N.K[i.k], mu[i.p, i.k], std[i.p, i.k]))
  }
  X.gbd <- rbind(X.gbd, tmp)
}
toy1 <- list(N = N, K = K, p.X = p.X, eta = eta, beta = beta, mu = mu,
             std = std, N.K = N.K, CLASS.gbd = CLASS.gbd, PV.gbd = PV.gbd,
             X.gbd = X.gbd)

### For toy2
eta <- c(0.8, 0.10, 0.10)
N.K <- tabulate(sample(1:K, N, prob = eta, replace = TRUE))
CLASS.gbd <- NULL
PV.gbd <- NULL
X.gbd <- NULL
for(i.k in 1:K){
  CLASS.gbd <- c(CLASS.gbd, rep(i.k, N.K[i.k]))
  PV.gbd <- c(PV.gbd, rbeta(N.K[i.k], beta[1, i.k], beta[2, i.k]))
  tmp <- NULL
  for(i.p in 1:p.X){
    tmp <- cbind(tmp, rnorm(N.K[i.k], mu[i.p, i.k], std[i.p, i.k]))
  }
  X.gbd <- rbind(X.gbd, tmp)
}
toy2 <- list(N = N, K = K, p.X = p.X, eta = eta, beta = beta, mu = mu,
             std = std, N.K = N.K, CLASS.gbd = CLASS.gbd, PV.gbd = PV.gbd,
             X.gbd = X.gbd)

# save(toy1, toy2, file = "../data/toy.rda")

