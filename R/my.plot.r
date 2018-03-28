### Plot fclust image.

plotpv <- function(da, posterior, PARAM, zlim = c(0, 0.01), plot.mean = TRUE,
    xlab = "", ylab = "", main = NULL, xlim = NULL, ylim = NULL,
    col = my.YlOrRd(), ignore.bg = FALSE){
  da.clust <- da
  id <- !is.na(da)
  K <- ncol(posterior)
  class <- apply(posterior, 1, which.max)

  if(is.null(main)){
    main <- paste("K = ", K, sep = "")
  }

  x <- 0:nrow(da)
  y <- 0:ncol(da)
  if(is.null(xlim)){
    xlim <- range(x)
  }
  if(is.null(ylim)){
    ylim <- range(y)
  }

  ### For the clusters.
  if(plot.mean){
    for(i.k in 1:K){
      alpha <- PARAM$BETA[[i.k]][1]
      beta <- PARAM$BETA[[i.k]][2]
      tmp <- which(id)[class == i.k]
      if(length(tmp) > 0){
        da.clust[tmp] <- alpha / (alpha + beta)
      }
    }
  }

  ### Plot.
  # image(da.clust, col = gray(0.2),
  #       axes = FALSE, main = main, xlim = xlim, ylim = ylim)
  # image(da.clust, zlim = zlim, col = my.YlOrRd(), add = TRUE,
  #       xlim = xlim, ylim = ylim)

  ### Background image.
  if(!ignore.bg){
    da.bg <- da
    da.bg[is.na(da)] <- 0
    da.bg[!is.na(da)] <- 1
    col.bg <- c("#000000", "#AAAAAA")
    n.x <- dim(da.bg)[1]
    n.y <- dim(da.bg)[2]
    image(1:n.x, 1:n.y, da.bg, col = col.bg, axes = FALSE,
          main = main, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
  }
  image(1:n.x, 1:n.y, da.clust, zlim = zlim, col = col, add = TRUE,
        xlim = xlim, ylim = ylim)

  invisible()
} # End of plotpv().

plotpvlegend <- function(zlim = c(0, 0.01), n.level = 20, main = NULL,
    col = my.YlOrRd()){
  z <- matrix(seq(zlim[1], zlim[2], length = n.level), nrow = 1)
  x <- 0
  y <- seq(zlim[1], zlim[2], length = n.level)
  image(x, y, z, zlim = zlim, col = col, axes = FALSE,
        xaxs = "r", yaxs = "r", xlab = "", ylab = "", main = main)
  axis(2)
  invisible()
} # End of plotpvlegend().
