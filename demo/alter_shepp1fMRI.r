library(MixfMRI, quietly = TRUE)
set.seed(1234)
phantom <- shepp1fMRI
da.all <- gendataset(phantom = phantom, overlap = 0.01)
da.pval <- da.all$pval
da.tval <- da.all$tval[!is.na(da.all$tval)]
pval.cutoff <- 0.05
qval.cutoff <- 0.05

### Bonferroni
cf.Bonferroni <- Threshold.Bonferroni(pval.cutoff, length(da.tval),
                                      type = "Normal") 
Bonferroni.class <- da.pval
Bonferroni.class[!is.na(da.pval)][da.tval >= cf.Bonferroni] <- 1
Bonferroni.class[!is.na(da.pval)][da.tval < cf.Bonferroni] <- 2

### FDR
cf.FDR <- Threshold.FDR(x = da.tval, q = qval.cutoff, type = "Normal") 
FDR.class <- da.pval
FDR.class[!is.na(da.pval)][da.tval >= cf.FDR] <- 1
FDR.class[!is.na(da.pval)][da.tval < cf.FDR] <- 2

### RF
cf.RF <- Threshold.RF(pval.cutoff, diag(1, 2), voxdim = c(1, 1),
                      num.vox = length(da.tval), type = "Normal") 
RF.class <- da.pval
RF.class[!is.na(da.pval)][da.tval >= cf.RF] <- 1
RF.class[!is.na(da.pval)][da.tval < cf.RF] <- 2

### CT 1st order neighborhood
x <- c(rep(NA, length(da.pval)), da.pval, rep(NA, length(da.pval)))
x[is.na(x)] <- 1
dim(x) <- c(dim(da.pval), 3)
nmat <- expand.grid(-1:1, -1:1, -1:1)[10:18,][-5,]
xx <- cluster.threshold(1 - x, nmat = nmat, level.thr = 0.999, 4)
CT1st.class <- da.pval
CT1st.class[xx[,, 2] != 0 & !is.na(da.pval)] <- 1
CT1st.class[xx[,, 2] == 0 & !is.na(da.pval)] <- 2

### CT 2nd order neighborhood
nmat <- expand.grid(-1:1, -1:1, -1:1)[c(11, 13, 15, 17),]
xx <- cluster.threshold(1 - x, nmat = nmat, level.thr = 0.999, 4)
CT2nd.class <- da.pval
CT2nd.class[xx[,, 2] != 0 & !is.na(da.pval)] <- 1
CT2nd.class[xx[,, 2] == 0 & !is.na(da.pval)] <- 2

### Plot
xlim <- ylim <- c(0.19, 0.82)
par(mfrow = c(2, 3), mar = c(0, 0, 3, 0))
image(Bonferroni.class, xlim = xlim, ylim = ylim, axes = FALSE,
      col = rev(my.YlOrRd(10)), main = "Boneferroni")
image(FDR.class, xlim = xlim, ylim = ylim, axes = FALSE,
      col = rev(my.YlOrRd(10)), main = "FDR")
image(RF.class, xlim = xlim, ylim = ylim, axes = FALSE,
      col = rev(my.YlOrRd(10)), main = "RF")
image(CT1st.class, xlim = xlim, ylim = ylim, axes = FALSE,
      col = rev(my.YlOrRd(10)), main = "cluster.threshold 1st")
image(CT2nd.class, xlim = xlim, ylim = ylim, axes = FALSE,
      col = rev(my.YlOrRd(10)), main = "cluster.threshold 2nd")
