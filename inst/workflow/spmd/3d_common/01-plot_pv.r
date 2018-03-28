library(MixfMRI, quietly = TRUE)
source("00-set_condition.r")

### Test 2d data.
da <- eval(parse(text = paste(case, sep = "")))
id <- !is.na(da)
PV.gbd <- da[id]
# PV.gbd[PV.gbd > 0.3] <- runif(sum(PV.gbd > 0.3), 0.3, 1)

### plot hist
fn.plot <- paste("./plot/", case, "_hist.pdf", sep = "")
pdf(fn.plot, height = 4, width = 6)
  hist(PV.gbd, nclass = 100,
       xlab = "p-value", main = case)
dev.off()
