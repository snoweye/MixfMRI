library(pbdMPI, quietly = TRUE)
init(set.seed = FALSE)
library(MixfMRI, quietly = TRUE)
source("00-set_condition.r")
.FC.CT$MPI.gbd <- TRUE

### Test 2d data.
da <- eval(parse(text = paste(case, sep = "")))
id <- !is.na(da)
PV.gbd <- da[id]
id.loc <- which(id, arr.ind = TRUE)
X.gbd <- t(t(id.loc) / dim(da))

# ret.PARAM <- NULL
for(i.k in K.min:K.max){
  set.seed(seed + i.k)

  i.try <- 1
  repeat{
    time <- system.time({
      PARAM <- try(fclust(X.gbd, PV.gbd, K = i.k))
    })
    if(.MixfMRIEnv$any(class(PARAM) != "try-error")){
      break
    } else{
      if(i.try > max.try.fclust){
        break
      }
      i.try <- i.try + 1
    }
  }

  if(.MixfMRIEnv$any(class(PARAM) == "try-error")){
    next
  }

  PARAM$time <- time
  PARAM$class <- unlist(allgather(PARAM$class))
  PARAM$posterior <- do.call("rbind", allgather(PARAM$posterior))

  if(comm.rank() == 0){
    print(PARAM)

    fn.out <- paste("./output/ret.", case, ".tp.K_", i.k, ".rda", sep = "")
    save(PARAM, file = fn.out)

    # ret.PARAM[[i.k]] <- PARAM
  }
}

if(comm.rank() == 0){
  # fn.out <- paste("./output/ret.", case, ".tp.rda", sep = "")
  # save(ret.PARAM, file = fn.out)
  print(proc.time())
}

finalize()
