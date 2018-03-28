library(pbdMPI, quietly = TRUE)
init(set.seed = FALSE)
library(MixfMRI, quietly = TRUE)

### Configuration
case <- "ds.pval"
.FC.CT$model.X <- "I"
.FC.CT$ignore.X <- FALSE
.FC.CT$CONTROL$debug <- 1
.FC.CT$CONTROL$RndEM.iter <- 50
.FC.CT$INIT$min.1st.prop <- 0.99
.FC.CT$INIT$max.PV <- 0.05
seed <- 1234
.FC.CT$MPI.gbd <- TRUE

### Test 3d data
da <- eval(parse(text = paste(case, sep = "")))
id <- !is.na(da)
PV.gbd <- da[id]
id.loc <- which(id, arr.ind = TRUE)
X.gbd <- t(t(id.loc) / dim(da))

### Set default
min.1st.prop <- .FC.CT$INIT$min.1st.prop
max.PV <- .FC.CT$INIT$max.PV
RndEM.iter <- .FC.CT$CONTROL$RndEM.iter
algorithm <- .FC.CT$algorithm[1]
model.X <- .FC.CT$model.X[1]
ignore.X <- .FC.CT$ignore.X
stop.unstable <- TRUE
MPI.gbd <- .FC.CT$MPI.gbd
common.gbd <- .FC.CT$common.gbd
K <- 5
set.seed(1234 + K)

### Initialization
if(comm.size() == 4){
  time <- system.time({
    PARAM.org <- set.global(X.gbd, PV.gbd, K = K,
                            min.1st.prop = min.1st.prop,
                            max.PV = max.PV,
                            RndEM.iter = RndEM.iter,
                            algorithm = algorithm[1],
                            model.X = model.X[1],
                            ignore.X = ignore.X,
                            MPI.gbd = MPI.gbd, common.gbd = common.gbd)
    PARAM.new <- initial.RndEM.gbd(PARAM.org)
  })
  comm.print(time)
  comm.print(PARAM.new)
  if(comm.rank() == 0){
    save(PARAM.new, file = "PARAM_new.rda")
  }
} else{
  load("PARAM_new.rda")
}

### APECMa
.FC.CT$algorithm <- "apecma"
time <- system.time({
  PARAM <- try(fclust(X.gbd, PV.gbd, K = K, PARAM.init = PARAM.new))
})
comm.print(time)
comm.print(PARAM)

### ECM
.FC.CT$algorithm <- "ecm"
time <- system.time({
  PARAM <- try(fclust(X.gbd, PV.gbd, K = K, PARAM.init = PARAM.new))
})
comm.print(time)
comm.print(PARAM)

### EM
.FC.CT$algorithm <- "em"
time <- system.time({
  PARAM <- try(fclust(X.gbd, PV.gbd, K = K, PARAM.init = PARAM.new))
})
comm.print(time)
comm.print(PARAM)

finalize()

