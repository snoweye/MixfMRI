library(MixfMRI, quietly = TRUE)
source("00-set_condition.r")

### Test 2d data
da <- eval(parse(text = paste(case, sep = "")))
id <- !is.na(da)
PV.gbd <- da[id]
id.loc <- which(id, arr.ind = TRUE)
X.gbd <- t(t(id.loc) / dim(da))

### Load results
fn.in <- paste("./output/ret.", case, ".tp.rda", sep = "")
load(fn.in)

mean.BETA <- list()
for(i.K in 2:length(ret.PARAM)){
  tmp <- lapply(ret.PARAM[[i.K]]$param$BETA[-1],
                function(x) x[1] / sum(x))
  mean.BETA[[i.K]] <- unlist(tmp) 
}

### Plot
zlim <- NULL
zlim <- range(unlist(mean.BETA))
# for(i.K in 2:length(ret.PARAM)){
#   if(is.null(zlim)){
#     zlim <- range(mean.BETA[[i.K]])
#   }
# 
#   fn.out <- paste("./plot/pv.", case, "_K_", i.K, ".pdf", sep = "")
#   pdf(fn.out, width = 6, height = 5)
#     nf <- layout(matrix(c(1, 1, 2, 3), ncol = 2, byrow = TRUE),
#                  widths = c(4, 2), heights = c(1, 4), respect = TRUE)
#     par(mar = c(0, 0, 0, 0))
#     plot(NULL, NULL, type = "n", xlim = c(0, 1), ylim = c(0, 1),
#          xlab = "", ylab = "", main = "", axes = FALSE)
#     text(0.5, 0.5, fn.out)
#     par(mar = c(5.1, 4.1, 4.1, 2.1))
# 
#     plotpv(da, ret.PARAM[[i.K]]$posterior, ret.PARAM[[i.K]]$param,
#            zlim = zlim)
# 
#     plotpvlegend(zlim = zlim)
#   dev.off()
# }

### All together
fn.out <- paste("./plot/pv_", case, ".pdf", sep = "")
pdf(fn.out, width = 8, height = 9)
  nf <- layout(matrix(c(rep(1, 4), 2:17), ncol = 4, byrow = TRUE),
               widths = rep(1, 4), heights = c(1, rep(4, 4)), respect = FALSE)
  par(mar = c(0, 0, 0, 0))
  plot(NULL, NULL, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "", ylab = "", main = "", axes = FALSE)
  text(0.5, 0.5, fn.out)
  par(mar = c(0, 0, 2, 0))

  for(i.K in 2:length(ret.PARAM)){
    plotpv(da, ret.PARAM[[i.K]]$posterior, ret.PARAM[[i.K]]$param,
           zlim = zlim)
  }

  par(mar = c(2.1, 5.1, 2.1, 5.1))
  plotpvlegend(zlim = zlim)
dev.off()

