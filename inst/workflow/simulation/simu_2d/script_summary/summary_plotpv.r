library(MixfMRI, quietly = TRUE)
prefix <- "simu"
load("new/summary/summary.rda")
if(!file.exists("new/summary_plotpv")){
  dir.create("new/summary_plotpv")
}

### Plot
i.counter <- 0
for(i.case in 1:nrow(ret.summary.org)){
  i.counter <- i.counter + 1

  tmp.phantom <- eval(parse(text = ret.summary.org[i.case, 1]))
  tmp.overlap <- eval(parse(text = ret.summary.org[i.case, 2]))
  set.seed(1234)
  da <- gendataset(phantom = tmp.phantom, overlap = tmp.overlap)$pval

  tmp.prefix <- ret.summary.org[i.case, 1:5]

  file <- paste("./", prefix, "/",
                paste(tmp.prefix, collapse = "_"),
                "/output/ret.simu.2d.tp", sep = "")

  file.rda <- paste(file, ".rda", sep = "")
  if(file.exists(file.rda)){
    load(file.rda)
    ret.PARAM.org <- ret.PARAM
  } else{
    ret.PARAM.org <- NULL
  }

  file.rda <- paste(file, ".new3.rda", sep = "")
  if(file.exists(file.rda)){
    load(file.rda)
    ret.PARAM.new <- ret.PARAM
  } else{
    ret.PARAM.new <- NULL
  }

  file.plot <- paste("./new/summary_plotpv/",
                     paste(tmp.prefix[-5], collapse = "_"),
                     ".pdf", sep = "")
  i.x <- 71
  i.y.1 <- 205
  i.y.2 <- 195
  i.cex.1 <- 0.5
  i.cex.2 <- 0.5
  i.adj <- c(0, 0)
  i.xlim <- c(70, 190)
  i.ylim <- c(50, 210)

  if(i.counter == 1){
    pdf(file.plot, height = 6, width = 18)
      par(mfrow = c(3, 9), mar = c(0, 0, 0, 0))
  }

    ### K.AIC
    i.k <- as.integer(ret.summary[i.case, 6])
    i.k.org <- as.integer(ret.summary.org[i.case, 6])
    if(is.na(i.k) || is.na(i.k.org)){
      plot(NULL, NULL, xlim = i.xlim, ylim = i.ylim,
           main = "", xlab = "", ylab = "", axes = FALSE)
    } else{ 
      plotpv(da, ret.PARAM.org[[i.k.org]]$posterior,
             ret.PARAM.org[[i.k.org]]$param,
             main = "", xlim = i.xlim, ylim = i.ylim, zlim = c(0, 0.4))
      text(i.x, i.y.1, labels = paste(tmp.prefix, collapse = "_"),
           cex = i.cex.1, adj = i.adj, col = "#FFFFFF")
      text(i.x, i.y.2, labels = paste("K.AIC=", i.k, sep = ""),
           cex = i.cex.2, adj = i.adj, col = "#FFFFFF")
    }

    ### K.BIC
    i.k <- as.integer(ret.summary[i.case, 7])
    i.k.org <- as.integer(ret.summary.org[i.case, 6])
    if(is.na(i.k) || is.na(i.k.org)){
      plot(NULL, NULL, xlim = i.xlim, ylim = i.ylim,
           main = "", xlab = "", ylab = "", axes = FALSE)
    } else{ 
      plotpv(da, ret.PARAM.org[[i.k.org]]$posterior,
             ret.PARAM.org[[i.k.org]]$param,
             main = "", xlim = i.xlim, ylim = i.ylim, zlim = c(0, 0.4))
      text(i.x, i.y.1, labels = paste(tmp.prefix, collapse = "_"),
           cex = i.cex.1, adj = i.adj, col = "#FFFFFF")
      text(i.x, i.y.2, labels = paste("K.BIC=", i.k, sep = ""),
           cex = i.cex.2, adj = i.adj, col = "#FFFFFF")
    }

    ### K.ICL.BIC
    i.k <- as.integer(ret.summary[i.case, 8])
    i.k.org <- as.integer(ret.summary.org[i.case, 8])
    if(is.na(i.k) || is.na(i.k.org)){
      plot(NULL, NULL, xlim = i.xlim, ylim = i.ylim,
           main = "", xlab = "", ylab = "", axes = FALSE)
    } else{ 
      plotpv(da, ret.PARAM.org[[i.k.org]]$posterior,
             ret.PARAM.org[[i.k.org]]$param,
             main = "", xlim = i.xlim, ylim = i.ylim, zlim = c(0, 0.4))
      text(i.x, i.y.1, labels = paste(tmp.prefix, collapse = "_"),
           cex = i.cex.1, adj = i.adj, col = "#FFFFFF")
      text(i.x, i.y.2, labels = paste("K.ICL.BIC=", i.k, sep = ""),
           cex = i.cex.2, adj = i.adj, col = "#FFFFFF")
    }

    ### K.AIC.merge
    i.k <- as.integer(ret.summary[i.case, 9])
    i.k.org <- as.integer(ret.summary.org[i.case, 8])
    if(is.na(i.k) || is.na(i.k.org)){
      plot(NULL, NULL, xlim = i.xlim, ylim = i.ylim,
           main = "", xlab = "", ylab = "", axes = FALSE)
    } else{ 
      plotpv(da, ret.PARAM.new[[i.k.org]]$posterior,
             ret.PARAM.new[[i.k.org]]$param,
             main = "", xlim = i.xlim, ylim = i.ylim, zlim = c(0, 0.4))
      text(i.x, i.y.1, labels = paste(tmp.prefix, collapse = "_"),
           cex = i.cex.1, adj = i.adj, col = "#FFFFFF")
      text(i.x, i.y.2, labels = paste("K.AIC.merge=", i.k, sep = ""),
           cex = i.cex.2, adj = i.adj, col = "#FFFFFF")
    }
    abline(v = i.x - 1)

    ### K.BIC.merge
    i.k <- as.integer(ret.summary[i.case, 10])
    i.k.org <- as.integer(ret.summary.org[i.case, 10])
    if(is.na(i.k) || is.na(i.k.org)){
      plot(NULL, NULL, xlim = i.xlim, ylim = i.ylim,
           main = "", xlab = "", ylab = "", axes = FALSE)
    } else{ 
      plotpv(da, ret.PARAM.new[[i.k.org]]$posterior,
             ret.PARAM.new[[i.k.org]]$param,
             main = "", xlim = i.xlim, ylim = i.ylim, zlim = c(0, 0.4))
      text(i.x, i.y.1, labels = paste(tmp.prefix, collapse = "_"),
           cex = i.cex.1, adj = i.adj, col = "#FFFFFF")
      text(i.x, i.y.2, labels = paste("K.BIC.merge=", i.k, sep = ""),
           cex = i.cex.2, adj = i.adj, col = "#FFFFFF")
    }

    ### K.ICL.BIC.merge
    i.k <- as.integer(ret.summary[i.case, 11])
    i.k.org <- as.integer(ret.summary.org[i.case, 10])
    if(is.na(i.k) || is.na(i.k.org)){
      plot(NULL, NULL, xlim = i.xlim, ylim = i.ylim,
           main = "", xlab = "", ylab = "", axes = FALSE)
    } else{ 
      plotpv(da, ret.PARAM.new[[i.k.org]]$posterior,
             ret.PARAM.new[[i.k.org]]$param,
             main = "", xlim = i.xlim, ylim = i.ylim, zlim = c(0, 0.4))
      text(i.x, i.y.1, labels = paste(tmp.prefix, collapse = "_"),
           cex = i.cex.1, adj = i.adj, col = "#FFFFFF")
      text(i.x, i.y.2, labels = paste("K.ICL.BIC.merge=", i.k, sep = ""),
           cex = i.cex.2, adj = i.adj, col = "#FFFFFF")
    }

    ### K.AIC.maitra
    i.k <- as.integer(ret.summary[i.case, 12])
    i.k.org <- as.integer(ret.summary.org[i.case, 12])
    if(is.na(i.k) || is.na(i.k.org)){
      plot(NULL, NULL, xlim = i.xlim, ylim = i.ylim,
           main = "", xlab = "", ylab = "", axes = FALSE)
    } else{ 
      plotpv(da, ret.PARAM.new[[i.k.org]]$posterior,
             ret.PARAM.new[[i.k.org]]$param,
             main = "", xlim = i.xlim, ylim = i.ylim, zlim = c(0, 0.4))
      text(i.x, i.y.1, labels = paste(tmp.prefix, collapse = "_"),
           cex = i.cex.1, adj = i.adj, col = "#FFFFFF")
      text(i.x, i.y.2, labels = paste("K.AIC.maitra=", i.k, sep = ""),
           cex = i.cex.2, adj = i.adj, col = "#FFFFFF")
    }
    abline(v = i.x - 1)

    ### K.BIC.maitra
    i.k <- as.integer(ret.summary[i.case, 13])
    i.k.org <- as.integer(ret.summary.org[i.case, 13])
    if(is.na(i.k) || is.na(i.k.org)){
      plot(NULL, NULL, xlim = i.xlim, ylim = i.ylim,
           main = "", xlab = "", ylab = "", axes = FALSE)
    } else{ 
      plotpv(da, ret.PARAM.new[[i.k.org]]$posterior,
             ret.PARAM.new[[i.k.org]]$param,
             main = "", xlim = i.xlim, ylim = i.ylim, zlim = c(0, 0.4))
      text(i.x, i.y.1, labels = paste(tmp.prefix, collapse = "_"),
           cex = i.cex.1, adj = i.adj, col = "#FFFFFF")
      text(i.x, i.y.2, labels = paste("K.BIC.maitra=", i.k, sep = ""),
           cex = i.cex.2, adj = i.adj, col = "#FFFFFF")
    }

    ### K.ICL.BIC.maitra
    i.k <- as.integer(ret.summary[i.case, 14])
    i.k.org <- as.integer(ret.summary.org[i.case, 14])
    if(is.na(i.k) || is.na(i.k.org)){
      plot(NULL, NULL, xlim = i.xlim, ylim = i.ylim,
           main = "", xlab = "", ylab = "", axes = FALSE)
    } else{ 
      plotpv(da, ret.PARAM.new[[i.k.org]]$posterior,
             ret.PARAM.new[[i.k.org]]$param,
             main = "", xlim = i.xlim, ylim = i.ylim, zlim = c(0, 0.4))
      text(i.x, i.y.1, labels = paste(tmp.prefix, collapse = "_"),
           cex = i.cex.1, adj = i.adj, col = "#FFFFFF")
      text(i.x, i.y.2, labels = paste("K.ICL.BIC.maitra=", i.k, sep = ""),
           cex = i.cex.2, adj = i.adj, col = "#FFFFFF")
    }

  if(i.counter == 3){
    dev.off()
    i.counter <- 0
  }
}
