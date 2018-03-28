library(pbdMPI, quietly = TRUE)
init(set.seed = FALSE)
library(MixfMRI, quietly = TRUE)
source("00-set_condition.r")
.FC.CT$MPI.gbd <- TRUE

### Test 3d data
da <- eval(parse(text = paste(case, sep = "")))
id <- !is.na(da)
PV.gbd <- da[id]
id.loc <- which(id, arr.ind = TRUE)
X.gbd <- t(t(id.loc) / dim(da))

### Load results
fn.in <- paste("./output/ret.", case, ".tp.rda", sep = "")
load(fn.in)

new.ret.PARAM <- ret.PARAM
for(i.k in K.min:K.max){
  if(length(ret.PARAM) >= i.k &&
     class(ret.PARAM[[i.k]]) != "try-error"){
    ret <- lrt2(PV.gbd, ret.PARAM[[i.k]]$class, ret.PARAM[[i.k]]$param$K)

    ### Check q-value.
    merge.id <- ret[, ncol(ret)] > 0.05 | is.na(ret[, ncol(ret)])
    merge.id[1] <- TRUE    # in case the first cluster is rejected.
    new.K <- i.k - sum(merge.id[-1])

    ### Check if merging is required.
    if(new.K != i.k){
      comm.cat("\nmerging K.org = ", i.k, " to K = ", new.K, "\n", sep = "",
               quiet = TRUE)

      ### Rebuild PARAM.init from a right dimension.
      PARAM.init <- ret.PARAM[[i.k]]$param    # copy from a right dimension.
      PARAM.init$N.gbd <- PARAM.init$N
      PARAM.init$N.all <- PARAM.init$N
      PARAM.init$K <- new.K
      PARAM.init$ETA <- c(sum(ret.PARAM[[i.k]]$param$ETA[merge.id]),
                          ret.PARAM[[i.k]]$param$ETA[!merge.id])
      PARAM.init$log.ETA <- log(PARAM.init$ETA)

      merge.id[1] <- FALSE    # retain the first cluster is necessary.
      PARAM.init$BETA <- ret.PARAM[[i.k]]$param$BETA[!merge.id]
      PARAM.init$MU <- matrix(ret.PARAM[[i.k]]$param$MU[, !merge.id],
                              ncol = new.K)
      PARAM.init$SIGMA <- ret.PARAM[[i.k]]$param$SIGMA[!merge.id]
      PARAM.init$logL <- NULL
      PARAM.init$initial.i.iter <- NULL

      ### Start a rerun from rebuilt PARAM.init.
      time <- system.time({
        PARAM.merge <- try(fclust(X.gbd, PV.gbd, K = new.K,
                                  PARAM.init = PARAM.init,
                                  stop.unstable = FALSE))
      })
      if(.MixfMRIEnv$any(class(PARAM.merge) == "try-error")){
        new.ret.PARAM[[i.k]] <- PARAM.merge
        next
      } else{
        PARAM.merge$time <- time
        new.ret.PARAM[[i.k]] <- PARAM.merge
      }

      new.ret.PARAM[[i.k]]$class <-
        unlist(allgather(new.ret.PARAM[[i.k]]$class))
      new.ret.PARAM[[i.k]]$posterior <-
        do.call("rbind", allgather(new.ret.PARAM[[i.k]]$posterior))

      if(comm.rank() == 0){
        print(new.ret.PARAM[[i.k]])
      }
    }
  }
}

if(comm.rank() == 0){
  ### Save to file
  ret.PARAM <- new.ret.PARAM
  fn.out <- paste("./output/ret.", case, ".tp.new.rda", sep = "")
  save(ret.PARAM, file = fn.out)
  print(proc.time())
}

finalize()
