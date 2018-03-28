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

slices <- dim(da)[3]
n.row <- floor(sqrt(slices))
n.col <- ceiling(slices / n.row) 
id <- rep(1:slices, each = prod(dim(da)[1:2]))[! is.na(da)]

### plotfclust
for(i.k in K.min:K.max){
  if(!is.null(ret.PARAM[[i.k]]) && length(ret.PARAM) >= i.k){
    fn.plot <- paste("./plot/", case, "_fclust.K_", i.k, ".pdf", sep = "")
    pdf(fn.plot, height = n.row * 2, width = n.col * 2)
      par(mfrow = c(n.row, n.col), mar = c(0, 0, 0, 0))
      for(i.slice in 1:slices){
        plotfclust(da[,, i.slice], ret.PARAM[[i.k]]$posterior[id == i.slice,],
                   main = "")
      }
    dev.off()
  }
}

### plotfclustpv
for(i.k in K.min:K.max){
  if(!is.null(ret.PARAM[[i.k]]) && length(ret.PARAM) >= i.k){
    fn.plot <- paste("./plot/", case, "_fclustpv.K_", i.k, ".pdf", sep = "")
    pdf(fn.plot, height = n.row * 2, width = n.col * 2)
      par(mfrow = c(n.row, n.col), mar = c(0, 0, 0, 0))
      for(i.slice in 1:slices){
        plotfclustpv(da[,, i.slice], ret.PARAM[[i.k]]$posterior[id == i.slice,],
                     main = "")
      }
    dev.off()
  }
}

### plotpv
for(i.k in K.min:K.max){
  if(!is.null(ret.PARAM[[i.k]]) && length(ret.PARAM) >= i.k){
    fn.plot <- paste("./plot/", case, "_fclustpv.K_", i.k, ".pdf", sep = "")
    pdf(fn.plot, height = n.row * 2, width = n.col * 2)
      par(mfrow = c(n.row, n.col), mar = c(0, 0, 0, 0))
      for(i.slice in 1:slices){
        plotpv(da[,, i.slice], ret.PARAM[[i.k]]$posterior[id == i.slice,],
               ret.PARAM[[i.k]]$param, main = "")
      }
    dev.off()
  }
}

