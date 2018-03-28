library(MixfMRI, quietly = TRUE)
source("00-set_condition.r")

### Test 3d data
da <- eval(parse(text = paste(case, sep = "")))
id <- !is.na(da)
PV.gbd <- da[id]
id.loc <- which(id, arr.ind = TRUE)
X.gbd <- t(t(id.loc) / dim(da))

### Load results
ret <- NULL
ret.PARAM <- NULL
for(i.k in K.min:K.max){
  fn.in <- paste("./output/ret.", case, ".tp.K_", i.k, ".rda", sep = "")
  load(fn.in)
  ret.PARAM[[i.k]] <- PARAM
  ret <- rbind(ret, c(i.k, ret.PARAM[[i.k]]$param$logL,
                      ret.PARAM[[i.k]]$aic, ret.PARAM[[i.k]]$bic))
}
colnames(ret) <- c("K", "logL", "AIC", "BIC")
ret <- data.frame(ret)

### Plot
fn.out <- paste("./plot/", case, "_ic.pdf", sep = "")
pdf(fn.out, width = 10, height = 5)
  par(mfrow = c(1, 2))
  plot(NULL, NULL, type = "b",
       xlim = range(ret$K), ylim = range(ret$logL),
       xlab = "K", ylab = "logL",
       main = case)
  points(ret$K, ret$logL, pch = 1, col = 1)
  lines(ret$K, ret$logL, lty = 1, col = 1)

  xlim <- range(ret$K)
  ylim <- range(c(ret$AIC, ret$BIC))
  plot(NULL, NULL, type = "n",
       xlim = xlim, ylim = ylim,
       xlab = "K", ylab = "IC",
       main = case)
  points(ret$K, ret$AIC, pch = 2, col = 2)
  lines(ret$K, ret$AIC, lty = 2, col = 2)
  points(ret$K, ret$BIC, pch = 3, col = 3)
  lines(ret$K, ret$BIC, lty = 3, col = 3)
  legend(xlim[2] - (xlim[2] - xlim[1]) * 0.5,
         ylim[2] - (ylim[2] - ylim[1]) * 0.1,
         c("AIC", "BIC"),
         pch = 2:3, lty = 2:3, col = 2:3)
dev.off()
