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

PARAM <- fclust(X.gbd, PV.gbd, K = 4)
print(PARAM)

id.pstats <- .MixfMRIEnv$CLASS.gbd
table(id.pstats)

### Section 3.3: eta = 0.05 with chi-square df=1
ret.lrt2 <- lrt2(PV.gbd, id.pstats, K = 4, H0.mean = 0.05, proc = "1")
print(ret.lrt2)

### Section 5.1: for paired mu's with chi-square df=2
ret.lrt.ab <- lrt.betaab(PV.gbd, id.pstats, K = 4, proc = "1")
print(ret.lrt.ab)

