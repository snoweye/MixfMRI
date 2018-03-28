library(MixfMRI, quietly = TRUE)
library(MASS, quietly = TRUE)
.FC.CT$model.X <- "I"
.FC.CT$CONTROL$debug <- 0
K <- 3

### Fit toy1.
set.seed(1234)
X.gbd <- toy1$X.gbd
X.range <- apply(X.gbd, 2, range)
X.gbd <- t((t(X.gbd) - X.range[1,]) / (X.range[2,] - X.range[1,]))
PV.gbd <- toy1$PV.gbd
fcobj <- fclust(X.gbd, PV.gbd, K = K, min.1st.prop = 0.5)

### Test cov matrix of posterior z.
x <- list(X.gbd = X.gbd, PV.gbd = PV.gbd)
post.z <- post.prob(x, fcobj)
cov.param <- cov.param(x, fcobj, post.z)
cov.logit.ETA <- cov.logit.ETA(x, fcobj, cov.param = cov.param$cov)

### Compute cov matrxi of z_nk - z_n1 for all k > 1.
A <- cbind(rep(-1, K - 1), diag(1, K - 1))
ETA <- fcobj$param$ETA
log.or <- log(ETA / (1 - ETA)) %*% t(A)
cov.log.or <- A %*% cov.logit.ETA %*% t(A)

### Check if 0 were within 95% confidence ellipsoid.
log.or %*% ginv(cov.log.or) %*% t(log.or)
partial <- MixfMRI:::partial.logit.p(fcobj$param$ETA)
log.or %*% A %*% partial %*%
  cov.param$I[1:(K-1), 1:(K-1)] %*%
  t(partial) %*% t(A) %*% t(log.or)

