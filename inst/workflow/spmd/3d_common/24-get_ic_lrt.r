library(MixfMRI, quietly = TRUE)
source("00-set_condition.r")

### Test 3d data
da <- eval(parse(text = paste(case, sep = "")))
id <- !is.na(da)
PV.gbd <- da[id]
id.loc <- which(id, arr.ind = TRUE)
X.gbd <- t(t(id.loc) / dim(da))

### Load results
fn.in <- paste("./output/ret.", case, ".tp.new2.rda", sep = "")
load(fn.in)

ret <- NULL
for(i.k in K.min:K.max){
  if(length(ret.PARAM) >= i.k &&
     class(ret.PARAM[[i.k]]) != "try-error"){
    ret <- rbind(ret, c(i.k, ret.PARAM[[i.k]]$param$K,
                        ret.PARAM[[i.k]]$param$logL,
                        ret.PARAM[[i.k]]$aic, ret.PARAM[[i.k]]$bic,
                        ret.PARAM[[i.k]]$icl.bic))
  }
}
colnames(ret) <- c("K.org", "K", "logL", "AIC", "BIC", "ICL-BIC")
print(ret)

# cat("\nLRT: H0: alpha = beta = 1 v.s. Ha: alpha < 1, beta > 1\n")
# for(i.k in K.min:K.max){
#   if(length(ret.PARAM) >= i.k &&
#      class(ret.PARAM[[i.k]]) != "try-error"){
#     ret <- lrt(PV.gbd, ret.PARAM[[i.k]]$class, ret.PARAM[[i.k]]$param$K)
#     cat("K.org = ", i.k, ", K = ", ret.PARAM[[i.k]]$param$K, "\n", sep = "")
#     print(ret)
#   }
# }

cat("\nLRT2: H0: mean >= 0.05 v.s. Ha: mean < 0.05\n")
for(i.k in K.min:K.max){
  if(length(ret.PARAM) >= i.k &&
     class(ret.PARAM[[i.k]]) != "try-error"){
    ret <- lrt2(PV.gbd, ret.PARAM[[i.k]]$class, ret.PARAM[[i.k]]$param$K)
    cat("K.org = ", i.k, ", K = ", ret.PARAM[[i.k]]$param$K, "\n", sep = "")
    print(ret)
  }
}
