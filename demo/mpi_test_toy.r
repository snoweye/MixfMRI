### SHELL> mpiexec -np 2 Rscript mpi_toy_test.r

library(MixfMRI, quietly = TRUE)
library(EMCluster, quietly = TRUE)
library(pbdMPI, quietly = TRUE)
init()
.FC.CT$MPI.gbd <- TRUE
# .FC.CT$algorithm <- "em"
# .FC.CT$model.X <- "V"
# .FC.CT$ignore.X <- TRUE

### Test toy1.
set.seed(1234)
X.gbd <- toy1$X.gbd
X.range <- apply(X.gbd, 2, range)
X.gbd <- t((t(X.gbd) - X.range[1,]) / (X.range[2,] - X.range[1,]))
PV.gbd <- toy1$PV.gbd
PARAM <- set.global(X.gbd, PV.gbd, K = 3)
PARAM.new <- initial.em.gbd(PARAM)
PARAM.toy1 <- em.step.gbd(PARAM.new)
id.toy1 <- do.call("c", allgather(.MixfMRIEnv$CLASS.gbd))
comm.print(PARAM.toy1$ETA)
comm.print(RRand(toy1$CLASS.gbd, id.toy1))

### Test toy2.
set.seed(1234)
X.gbd <- toy2$X.gbd
X.range <- apply(X.gbd, 2, range)
X.gbd <- t((t(X.gbd) - X.range[1,]) / (X.range[2,] - X.range[1,]))
PV.gbd <- toy2$PV.gbd
PARAM <- set.global(X.gbd, PV.gbd, K = 3)
PARAM.new <- initial.em.gbd(PARAM)
PARAM.toy2 <- em.step.gbd(PARAM.new)
id.toy2 <- do.call("c", allgather(.MixfMRIEnv$CLASS.gbd))
comm.print(PARAM.toy2$ETA)
comm.print(RRand(toy2$CLASS.gbd, id.toy2))

finalize()
