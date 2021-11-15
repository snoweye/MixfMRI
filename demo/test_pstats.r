library(MixfMRI, quietly = TRUE)
library(EMCluster, quietly = TRUE)
# .FC.CT$algorithm <- "apecma"
# .FC.CT$model.X <- "V"
# .FC.CT$ignore.X <- TRUE
set.seed(1234)

### Test pstats
id <- which(!is.na(pstats))
PV.gbd <- pstats[id]

X.gbd <- which(!is.na(pstats), arr.ind = TRUE)
X.range <- apply(X.gbd, 2, range)
X.gbd <- t((t(X.gbd) - X.range[1,]) / (X.range[2,] - X.range[1,]))

PARAM <- fclust(X.gbd, PV.gbd, K = 3)
print(PARAM)

id.pstats <- .MixfMRIEnv$CLASS.gbd
table(id.pstats)

