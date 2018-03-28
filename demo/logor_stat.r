library(MixfMRI, quietly = TRUE)
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

### Check results.
post.z[id.merge,]
cbind(toy1$X.gbd[id.merge,], toy1$PV.gbd[id.merge])
