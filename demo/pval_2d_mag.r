library(MixfMRI, quietly = TRUE)
set.seed(1234)

### Test 2d data.
da <- pval.2d.mag
id <- !is.na(da)
PV.gbd <- da[id]
id.loc <- which(id, arr.ind = TRUE)
X.gbd <- t(t(id.loc) / dim(da))
ret <- fclust(X.gbd, PV.gbd, K = 3)
print(ret)

### p-values of rest clusters.
ret.lrt <- lrt(PV.gbd, ret$class, K = 3)
print(ret.lrt)
ret.lrt2 <- lrt2(PV.gbd, ret$class, K = 3)
print(ret.lrt2)

### Plotting.
par(mfrow = c(2, 2), mar = c(0, 0, 2, 0))
plotfclust(da, ret$posterior, main = "Posterior")
plotfclustpv(da, ret$posterior, main = "p-value")
plotpv(da, ret$posterior, ret$param, main = "Mean of Beta")
par(mar = c(5.1, 4.1, 4.1, 2.1))
hist(PV.gbd, nclass = 40, main = "p-value")
