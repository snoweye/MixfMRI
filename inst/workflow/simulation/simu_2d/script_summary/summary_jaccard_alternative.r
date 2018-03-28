options("width"=240)

library(MixfMRI, quietly = TRUE)
library(AnalyzeFMRI, quietly = TRUE)

phantom <- c("shepp0fMRI", "shepp1fMRI", "shepp2fMRI")
overlap <- c("0.01", "0.1", "0.25", "0.50", "0.75", "0.95")

if(!file.exists("new/summary_alternative")){
  dir.creat("new/summary_alternative")
}

pval.cutoff <- 0.05
qval.cutoff <- 0.05

i.x <- 71
i.y.1 <- 205
i.y.2 <- 195
i.cex.1 <- 0.5
i.cex.2 <- 0.5
i.adj <- c(0, 0)
i.xlim <- c(70, 190)
i.ylim <- c(50, 210)

### Table and Plot
ret.alt.table <- NULL
for(i.phantom in 1:length(phantom)){

  file.plot <- paste("./new/summary_alternative/", phantom[i.phantom],
                     "_jaccard.pdf", sep = "")

  pdf(file.plot, height = 12, width = 6)
    par(mfrow = c(6, 3), mar = c(0, 0, 0, 0))

  for(i.overlap in 1:length(overlap)){
    tmp.phantom <- eval(parse(text = phantom[i.phantom]))
    tmp.overlap <- eval(parse(text = overlap[i.overlap]))
    set.seed(1234)
    da.all <- gendataset(phantom = tmp.phantom, overlap = tmp.overlap)
    da.pval <- da.all$pval
    da.tval <- da.all$tval[!is.na(da.all$tval)]

    ### True class
    tmp.class <- tmp.phantom[!is.na(tmp.phantom)] + 1
    tmp.class[tmp.class == 1] <- 0
    tmp.class[tmp.class > 1] <- 1

    ### Bonferroni
    cf.Bonferroni <- Threshold.Bonferroni(pval.cutoff, length(da.tval),
                                          type = "Normal") 
    Bonferroni.class <- rep(0, length(da.tval))
    Bonferroni.class[da.tval > cf.Bonferroni] <- 1
    Bonferroni.R <- Jaccard.Index(Bonferroni.class, tmp.class)
    ### plot
    posterior <- matrix(0, nrow = length(tmp.class), ncol = 2)
    posterior[Bonferroni.class == 0, 1] <- 1
    posterior[Bonferroni.class == 1, 2] <- 1
    plotfclust(da.pval, posterior, main = "",
               xlim = i.xlim, ylim = i.ylim)
    text(i.x, i.y.1, labels = "Bonferroni",
         cex = i.cex.1, adj = i.adj)

    ### FDR
    cf.FDR <- Threshold.FDR(x = da.tval, q = qval.cutoff, type = "Normal") 
    FDR.class <- rep(0, length(da.tval))
    FDR.class[da.tval > cf.FDR] <- 1
    FDR.R <- Jaccard.Index(FDR.class, tmp.class)
    ### plot
    posterior <- matrix(0, nrow = length(tmp.class), ncol = 2)
    posterior[FDR.class == 0, 1] <- 1
    posterior[FDR.class == 1, 2] <- 1
    plotfclust(da.pval, posterior, main = "",
               xlim = i.xlim, ylim = i.ylim)
    text(i.x, i.y.1, labels = "FDR",
         cex = i.cex.1, adj = i.adj)

    ### RF
    cf.RF <- Threshold.RF(pval.cutoff, diag(1, 2), voxdim = c(1, 1),
                          num.vox = length(da.tval), type = "Normal") 
    RF.class <- rep(0, length(da.tval))
    RF.class[da.tval > cf.RF] <- 1
    RF.R <- Jaccard.Index(RF.class, tmp.class)
    ### plot
    posterior <- matrix(0, nrow = length(tmp.class), ncol = 2)
    posterior[RF.class == 0, 1] <- 1
    posterior[RF.class == 1, 2] <- 1
    plotfclust(da.pval, posterior, main = "",
               xlim = i.xlim, ylim = i.ylim)
    text(i.x, i.y.1, labels = "RF",
         cex = i.cex.1, adj = i.adj)

    ### CT 1st order neighborhood
    x <- c(rep(NA, length(da.pval)), da.pval, rep(NA, length(da.pval)))
    x[is.na(x)] <- 1
    dim(x) <- c(dim(da.pval), 3)
    nmat <- expand.grid(-1:1, -1:1, -1:1)[10:18,][-5,]
    xx <- cluster.threshold(1 - x, nmat = nmat, level.thr = 0.999, 4)
    CT1st.class <- da.pval
    ### When 0 and 1 are wrong, this can affect JI but not adjR. 
    # CT1st.class[xx[,, 2] != 0 & !is.na(da.pval)] <- 0
    # CT1st.class[xx[,, 2] == 0 & !is.na(da.pval)] <- 1
    CT1st.class[xx[,, 2] != 0 & !is.na(da.pval)] <- 1
    CT1st.class[xx[,, 2] == 0 & !is.na(da.pval)] <- 0
    CT1st.class <- CT1st.class[!is.na(da.pval)]
    CT1st.R <- Jaccard.Index(CT1st.class, tmp.class)
    ### plot
    posterior <- matrix(0, nrow = length(tmp.class), ncol = 2)
    posterior[CT1st.class == 0, 1] <- 1
    posterior[CT1st.class == 1, 2] <- 1
    plotfclust(da.pval, posterior, main = "",
               xlim = i.xlim, ylim = i.ylim)
    text(i.x, i.y.1, labels = "CT1st",
         cex = i.cex.1, adj = i.adj)

    ### CT 2nd order neighborhood
    nmat <- expand.grid(-1:1, -1:1, -1:1)[c(11, 13, 15, 17),]
    xx <- cluster.threshold(1 - x, nmat = nmat, level.thr = 0.999, 4)
    CT2nd.class <- da.pval
    ### When 0 and 1 are wrong, this can affect JI but not adjR. 
    # CT2nd.class[xx[,, 2] != 0 & !is.na(da.pval)] <- 0
    # CT2nd.class[xx[,, 2] == 0 & !is.na(da.pval)] <- 1
    CT2nd.class[xx[,, 2] != 0 & !is.na(da.pval)] <- 1
    CT2nd.class[xx[,, 2] == 0 & !is.na(da.pval)] <- 0
    CT2nd.class <- CT2nd.class[!is.na(da.pval)]
    CT2nd.R <- Jaccard.Index(CT2nd.class, tmp.class)
    ### plot
    posterior <- matrix(0, nrow = length(tmp.class), ncol = 2)
    posterior[CT2nd.class == 0, 1] <- 1
    posterior[CT2nd.class == 1, 2] <- 1
    plotfclust(da.pval, posterior, main = "",
               xlim = i.xlim, ylim = i.ylim)
    text(i.x, i.y.1, labels = "CT2nd",
         cex = i.cex.1, adj = i.adj)

    ### FDR (HY2001)
    pvalue <- pnorm(da.tval, lower.tail = FALSE)
    qvalue <- qvalue(pvalue, method = "BY2001")
    HY.class <- rep(0, length(da.tval))
    HY.class[qvalue < qval.cutoff] <- 1
    HY.R <- Jaccard.Index(HY.class, tmp.class)
    ### plot
    posterior <- matrix(0, nrow = length(tmp.class), ncol = 2)
    posterior[HY.class == 0, 1] <- 1
    posterior[HY.class == 1, 2] <- 1
    plotfclust(da.pval, posterior, main = "",
               xlim = i.xlim, ylim = i.ylim)
    text(i.x, i.y.1, labels = "HY",
         cex = i.cex.1, adj = i.adj)

    ### Return
    tmp <- c(phantom[i.phantom], overlap[i.overlap], Bonferroni.R, FDR.R, RF.R,
             CT1st.R, CT2nd.R, HY.R)
    ret.alt.table <- rbind(ret.alt.table, tmp)
  }

  dev.off()
}

colnames(ret.alt.table) <-
  c("phantom", "overlap", "Bonferroni.JI", "FDR.JI", "RF.JI", "CT1st.JI", "CT2nd.JI",
    "HY.JI")
rownames(ret.alt.table) <- NULL
print(ret.alt.table)

save(ret.alt.table, file = "new/summary_alternative/summary_jaccard_alternative.rda")
