library(MixfMRI, quietly = TRUE)
source("00-set_condition.r")

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
new.ret <- NULL
for(i.k in K.min:K.max){
  if(length(ret.PARAM) >= i.k &&
     class(ret.PARAM[[i.k]]) != "try-error"){
    ret <- lrt2(PV.gbd, ret.PARAM[[i.k]]$class, ret.PARAM[[i.k]]$param$K,
                proc = "1")
    cat("K = ", i.k, sep = "")
    print(ret)

    ### Check q-value.
    merge.id <- ret[, ncol(ret)] > 0.05 | is.na(ret[, ncol(ret)])
    merge.id[1] <- TRUE    # in case the first cluster is rejected.
    new.K <- i.k - sum(merge.id[-1])

    ### Check if merging is required.
    if(new.K != i.k){
      cat("\nmerging K.org = ", i.k, " to K = ", new.K, "\n", sep = "",
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

      ### Do merge fails to the 1st cluster.
      tmp.posterior <- ret.PARAM[[i.k]]$posterior
      tmp.posterior[, 1] <- rowSums(tmp.posterior[, c(1, which(merge.id))])
      tmp.posterior <- matrix(tmp.posterior[, !merge.id], ncol = new.K)

      new.ret.PARAM[[i.k]]$param <- PARAM.init
      new.ret.PARAM[[i.k]]$posterior <- tmp.posterior
      new.ret.PARAM[[i.k]]$class <- apply(tmp.posterior, 1, which.max)
      new.ret.PARAM[[i.k]]$n.class <- tabulate(new.ret.PARAM[[i.k]]$class,
                                               nbins = new.K) 
    }

    new.ret <- rbind(new.ret, c(i.k, new.K))
  }
}
colnames(new.ret) <- c("K.org", "K")
print(new.ret)

### Save to file
ret.PARAM <- new.ret.PARAM
fn.out <- paste("./output/ret.", case, ".tp.new3.rda", sep = "")
save(ret.PARAM, file = fn.out)
