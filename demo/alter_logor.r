library(MixfMRI, quietly = TRUE)
set.seed(1234)
phantom <- shepp2fMRI

### Generate 2d data.
da <- gendataset(phantom = phantom, overlap = 0.01)$pval
id <- !is.na(da)
PV.gbd <- da[id]
id.loc <- which(id, arr.ind = TRUE)
X.gbd <- t(t(id.loc) / dim(da))

### Test 2d data.
i.k <- 4
set.seed(1234 + i.k)
fcobj <- fclust(X.gbd, PV.gbd, K = i.k)

### Test log odds ratio.
x <- list(X.gbd = X.gbd, PV.gbd = PV.gbd)
post.z <- post.prob(x, fcobj)
lor <- logor.stat(x, fcobj, post.z)

### Check if 95% CE covers log odd ratio = 1.
id.notna <- !is.na(lor$df)
id.cover <- which(lor$test.stat[id.notna] < pchisq(0.95, lor$df[id.notna]))

### Get voxels needed for merging.
id.active <- which(fcobj$class != 1)
id.merge <- id.active[id.notna][id.cover]

### Get new class.
ret.class <- da
ret.class[id] <- fcobj$class
logor.class <- ret.class
logor.class[id][id.merge] <- 1

### Plot
xlim <- ylim <- c(0.19, 0.82)
par(mfrow = c(2, 2), mar = c(0, 0, 3, 0))
image(da, xlim = xlim, ylim = ylim, axes = FALSE,
      col = rev(my.YlOrRd(10)), main = "Simulated True")
image(ret.class, xlim = xlim, ylim = ylim, axes = FALSE,
      col = rev(my.YlGnBu(10)), main = "Estimated (K = 4)")
image(logor.class, xlim = xlim, ylim = ylim, axes = FALSE,
      col = rev(my.YlGnBu(10)), main = "Log Odds Ratio")

