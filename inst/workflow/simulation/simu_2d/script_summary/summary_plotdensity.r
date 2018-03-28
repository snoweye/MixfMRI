phantom <- c("shepp0fMRI", "shepp1fMRI", "shepp2fMRI")
overlap <- c("0.01", "0.1", "0.25", "0.50", "0.75", "0.95")

library(MixfMRI, quietly = TRUE)
prefix <- "simu"

col.1 <- rgb(1, 0, 0)
col.2 <- rgb(0, 1, 0)
col.3 <- rgb(0, 0, 1)
x <- exp(seq(log(1e-8), log(0.999), length = 50))
pv.xlim = c(50, 206)
pv.ylim = c(50, 206)
pv.zlim = c(0, 0.05)
pv.col = rev(my.YlOrRd())

### Get density and data.
for(i.phantom in 1:length(phantom)){
  tmp.phantom <- eval(parse(text = phantom[i.phantom]))
  tmp.class <- tmp.phantom[!is.na(tmp.phantom)] + 1

  ret.da <- list()
  ret.y.1 <- list()
  ret.y.2 <- list()
  ret.y.3 <- list()
  ret.y.m <- list()
  ret.pval <- list()
  for(i.overlap in 1:length(overlap)){
    tmp.overlap <- eval(parse(text = overlap[i.overlap]))

    set.seed(1234)
    da.all <- gendataset(phantom = tmp.phantom, overlap = tmp.overlap)
    eta <- da.all$eta
    mu <- da.all$mu

    ret.da[[length(ret.da) + 1]] <- da.all$pval
    ret.y.1[[length(ret.y.1) + 1]] <- dpval(x, mu[1])
    ret.y.2[[length(ret.y.2) + 1]] <- dpval(x, mu[2])
    ret.y.3[[length(ret.y.3) + 1]] <- dpval(x, mu[3])
    ret.y.m[[length(ret.y.m) + 1]] <- dmixpval(x, eta, mu)
    ret.pval[[length(ret.pval) + 1]] <- da.all$pval[!is.na(da.all$pval)]
  }

  ret.all.y.1 <- do.call("cbind", ret.y.1)
  ret.all.y.2 <- do.call("cbind", ret.y.2)
  ret.all.y.3 <- do.call("cbind", ret.y.3)
  ret.all.y.m <- do.call("cbind", ret.y.m)
  ret.all.pval <- do.call("cbind", ret.pval)
  xlim <- range(x)
  ylim <- range(c(ret.all.y.1, ret.all.y.2, ret.all.y.3, ret.all.y.m))

  ### Get histogram
  tmp.h <- list()
  ylim.h <- NULL
  for(i.overlap in 1:length(overlap)){
    tmp.h[[i.overlap]] <- hist(ret.all.pval[, i.overlap], nclass = 50,
                               plot = FALSE)
    ylim.h <- c(ylim.h, max(tmp.h[[i.overlap]]$counts))
  }
  ylim.h <- c(0, max(ylim.h))
  
  ### Start to plot.
  fn.out <- paste("./new/summary_plotdensity/", phantom[i.phantom], "_0.975.pdf", sep = "")
  pdf(fn.out, height = 3 * length(overlap), width = 3 * 5)
    par(mfcol = c(length(overlap), 5))
  
    ### Density plot.
    for(i.overlap in 1:length(overlap)){
      par(mar = c(5.1, 4.1, 4.1, 2.1))
      plot(NULL, NULL, xlim = xlim, ylim = ylim, log = "xy",
           main = "Distribution", xlab = "p-values", ylab = "Density")
      lines(list(x = x, y = ret.all.y.m[, i.overlap]),
            lty = 1, col = 1, lwd = 2)
      lines(list(x = x, y = ret.all.y.1[, i.overlap]),
            lty = 2, col = col.1, lwd = 2)
      lines(list(x = x, y = ret.all.y.2[, i.overlap]),
            lty = 3, col = col.2, lwd = 2)
      lines(list(x = x, y = ret.all.y.3[, i.overlap]),
            lty = 4, col = col.3, lwd = 2)
      legend(1e-8, 1e-3,
             c("Mixture", "Component 1", "Component 2", "Component 3"),
             lty = 1:4, col = c(1, col.1, col.2, col.3), lwd = 2)
    }
  
    ### Historgram.
    for(i.overlap in 1:length(overlap)){
      plot(tmp.h[[i.overlap]], ylim = ylim.h,
           main = "Samples", xlab = "p-values", ylab = "Frequence")
      abline(v = 0.05, lty = 3, col = 2, lwd = 1)
    }
  
    ### noX, ICL-BIC
    par(mar = c(0, 0, 1.5, 0.5))
    for(i.overlap in 1:length(overlap)){
      fn.in <- paste("./", prefix, "/", phantom[i.phantom],
                     "_", overlap[i.overlap],
                     "_0.975_0.05_noX/output/ret.simu.2d.tp.rda", sep = "")
      if(!file.exists(fn.in)){
        plot(NULL, NULL, type = "n", xlim = c(0, 1), ylim = c(0, 1),
             axes = FALSE, main = "", xlab = "", ylab = "")
        next
      }
      load(fn.in)
      tmp.ic <- Inf
      for(i.k in 2:length(ret.PARAM)){
        if(is.null(ret.PARAM[[i.k]]$icl.bic)){
          tmp.ic <- c(tmp.ic, NA)
        } else{
          tmp.ic <- c(tmp.ic, ret.PARAM[[i.k]]$icl.bic)
        }
      }
      id.k <- which.min(tmp.ic)
  
      fn.in <- paste("./", prefix, "/", phantom[i.phantom],
                     "_", overlap[i.overlap],
                     "_0.975_0.05_noX/output/ret.simu.2d.tp.new3.rda", sep = "")
      if(!file.exists(fn.in)){
        plot(NULL, NULL, type = "n", xlim = c(0, 1), ylim = c(0, 1),
             axes = FALSE, main = "", xlab = "", ylab = "")
        next
      }
      load(fn.in)
      tmp.posterior <- ret.PARAM[[id.k]]$posterior
      tmp.PARAM <- ret.PARAM[[id.k]]$param
  
      i.R <- RRand(ret.PARAM[[id.k]]$class, tmp.class)$adjRand
      i.R <- sprintf("%.4f", i.R)

      plotpv(ret.da[[i.overlap]], tmp.posterior, tmp.PARAM,
             main = "", xlab = "", ylab = "",
             xlim = pv.xlim, ylim = pv.ylim, zlim = pv.zlim, col = pv.col)
      title(main = paste("noX, hatK=", tmp.PARAM$K, ", adjR=", i.R, sep = ""),
            cex.main = 1.0)
    }
  
    ### Case I, ICL-BIC
    par(mar = c(0, 0, 1.5, 0.5))
    for(i.overlap in 1:length(overlap)){
      fn.in <- paste("./", prefix, "/", phantom[i.phantom],
                     "_", overlap[i.overlap],
                     "_0.975_0.05_I/output/ret.simu.2d.tp.rda", sep = "")
      if(!file.exists(fn.in)){
        plot(NULL, NULL, type = "n", xlim = c(0, 1), ylim = c(0, 1),
             axes = FALSE, main = "", xlab = "", ylab = "")
        next
      }
      load(fn.in)
      tmp.ic <- Inf
      for(i.k in 2:length(ret.PARAM)){
        if(is.null(ret.PARAM[[i.k]]$icl.bic)){
          tmp.ic <- c(tmp.ic, NA)
        } else{
          tmp.ic <- c(tmp.ic, ret.PARAM[[i.k]]$icl.bic)
        }
      }
      id.k <- which.min(tmp.ic)
  
      fn.in <- paste("./", prefix, "/", phantom[i.phantom],
                     "_", overlap[i.overlap],
                     "_0.975_0.05_I/output/ret.simu.2d.tp.new3.rda", sep = "")
      if(!file.exists(fn.in)){
        plot(NULL, NULL, type = "n", xlim = c(0, 1), ylim = c(0, 1),
             axes = FALSE, main = "", xlab = "", ylab = "")
        next
      }
      load(fn.in)
      tmp.posterior <- ret.PARAM[[id.k]]$posterior
      tmp.PARAM <- ret.PARAM[[id.k]]$param
  
      i.R <- RRand(ret.PARAM[[id.k]]$class, tmp.class)$adjRand
      i.R <- sprintf("%.4f", i.R)

      plotpv(ret.da[[i.overlap]], tmp.posterior, tmp.PARAM,
             main = "", xlab = "", ylab = "",
             xlim = pv.xlim, ylim = pv.ylim, zlim = pv.zlim, col = pv.col)
      title(main = paste("I, hatK=", tmp.PARAM$K, ", adjR=", i.R, sep = ""),
            cex.main = 1.0)
    }
  
    ### Case V, ICL-BIC
    par(mar = c(0, 0, 1.5, 0.5))
    for(i.overlap in 1:length(overlap)){
      fn.in <- paste("./", prefix, "/", phantom[i.phantom],
                     "_", overlap[i.overlap],
                     "_0.975_0.05_V/output/ret.simu.2d.tp.rda", sep = "")
      if(!file.exists(fn.in)){
        plot(NULL, NULL, type = "n", xlim = c(0, 1), ylim = c(0, 1),
             axes = FALSE, main = "", xlab = "", ylab = "")
        next
      }
      load(fn.in)
      tmp.ic <- Inf
      for(i.k in 2:length(ret.PARAM)){
        if(is.null(ret.PARAM[[i.k]]$icl.bic)){
          tmp.ic <- c(tmp.ic, NA)
        } else{
          tmp.ic <- c(tmp.ic, ret.PARAM[[i.k]]$icl.bic)
        }
      }
      id.k <- which.min(tmp.ic)
  
      fn.in <- paste("./", prefix, "/", phantom[i.phantom],
                     "_", overlap[i.overlap],
                     "_0.975_0.05_V/output/ret.simu.2d.tp.new3.rda", sep = "")
      if(!file.exists(fn.in)){
        plot(NULL, NULL, type = "n", xlim = c(0, 1), ylim = c(0, 1),
             axes = FALSE, main = "", xlab = "", ylab = "")
        next
      }
      load(fn.in)
      tmp.posterior <- ret.PARAM[[id.k]]$posterior
      tmp.PARAM <- ret.PARAM[[id.k]]$param

      i.R <- RRand(ret.PARAM[[id.k]]$class, tmp.class)$adjRand
      i.R <- sprintf("%.4f", i.R)

      plotpv(ret.da[[i.overlap]], tmp.posterior, tmp.PARAM,
             main = "", xlab = "", ylab = "",
             xlim = pv.xlim, ylim = pv.ylim, zlim = pv.zlim, col = pv.col)
      title(main = paste("V, hatK=", tmp.PARAM$K, ", adjR=", i.R, sep = ""),
            cex.main = 1.0)
    }
  dev.off()
}
