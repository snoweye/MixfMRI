### SHELL> mpiexec -np 4 Rscript mpi_fclust2d.r

library(MixfMRI, quietly = TRUE)
library(EMCluster, quietly = TRUE)
library(pbdMPI, quietly = TRUE)
init()
.FC.CT$MPI.gbd <- TRUE
# .FC.CT$algorithm <- "em"
# .FC.CT$model.X <- "V"
# .FC.CT$ignore.X <- TRUE
set.seed(1234)

### Test toy1.
X.gbd <- toy1$X.gbd[, -3]
X.range <- apply(X.gbd, 2, range)
X.gbd <- t((t(X.gbd) - X.range[1,]) / (X.range[2,] - X.range[1,]))
PV.gbd <- toy1$PV.gbd
PARAM <- fclust(X.gbd, PV.gbd, K = 3)
print(PARAM)
id.toy1 <- do.call("c", allgather(.MixfMRIEnv$CLASS.gbd))
comm.print(RRand(toy1$CLASS.gbd, id.toy1))

### Test toy2.
X.gbd <- toy2$X.gbd[, -3]
X.range <- apply(X.gbd, 2, range)
X.gbd <- t((t(X.gbd) - X.range[1,]) / (X.range[2,] - X.range[1,]))
PV.gbd <- toy2$PV.gbd
PARAM <- fclust(X.gbd, PV.gbd, K = 3)
print(PARAM)
id.toy2 <- do.call("c", allgather(.MixfMRIEnv$CLASS.gbd))
comm.print(RRand(toy2$CLASS.gbd, id.toy2))

finalize()
