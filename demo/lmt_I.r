library(MixfMRI, quietly = TRUE)
library(EMCluster, quietly = TRUE)
.FC.CT$model.X <- "I"
.FC.CT$CONTROL$debug <- 0
.FC.CT$check.X.unit <- FALSE

### Fit toy1.
set.seed(1234)
X.gbd <- toy1$X.gbd
PV.gbd <- toy1$PV.gbd
ret.2 <- fclust(X.gbd, PV.gbd, K = 2)
ret.3 <- fclust(X.gbd, PV.gbd, K = 3)
ret.4 <- fclust(X.gbd, PV.gbd, K = 4)
ret.5 <- fclust(X.gbd, PV.gbd, K = 5)

### ARI
RRand(toy1$CLASS.gbd, ret.2$class)
RRand(toy1$CLASS.gbd, ret.3$class)
RRand(toy1$CLASS.gbd, ret.4$class)
RRand(toy1$CLASS.gbd, ret.5$class)

### Test toy1.
(lmt.23 <- lmt.I(ret.2, ret.3, X.gbd, PV.gbd))
(lmt.24 <- lmt.I(ret.2, ret.4, X.gbd, PV.gbd))
(lmt.25 <- lmt.I(ret.2, ret.5, X.gbd, PV.gbd))
(lmt.34 <- lmt.I(ret.3, ret.4, X.gbd, PV.gbd))
(lmt.35 <- lmt.I(ret.3, ret.5, X.gbd, PV.gbd))
(lmt.45 <- lmt.I(ret.4, ret.5, X.gbd, PV.gbd))
