plotfclustpv <- function(da, posterior, main = NULL, xlim = NULL, ylim = NULL){
  ### Force the col.image here.
  col.image <- c("#000000",
                 "#F9070E", "#0063FF", "#01B901", "#C000C5", "#FF7F00",
                 "#FFFF33", "#A65628", "#F781BF",
                 "#AF0000", "#009600", "#0000FF", "#FF4500",
                 "#00AFAF", "#8A2BE2", "#458B74")

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

  image(x, y, da, zlim = c(0, 1), col = NULL, axes = FALSE, main = main,
        xlim = xlim, ylim = ylim)
  for(i.k in 1:K){
    id.clust <- id
    id.clust[id.clust == TRUE][class != i.k] <- FALSE
    da.clust <- da
    da.clust[!id.clust] <- NA

    col <- col.image[(i.k - 1) %% length(col.image) + 1]
    col.clust <- rev(my.alpha.col(col))
    image(x, y, da.clust, zlim = c(0, 1), col = col.clust, add = TRUE,
          xlim = xlim, ylim = ylim)
  }

  invisible()
} # End of plotfclustpv().
