my.alpha.col <- function(x, tl.level = 64, from = "#333333", to = "#FFFFFF"){
  # alpha <- colorRampPalette(c(from, to))(tl.level + 1)[-1]
  # alpha <- gsub(".....", "", alpha)
  col.fill <- colorRampPalette(c("#EEEEEE", x))(tl.level + 1)[-1]
  # col.fill <- paste(col.fill, alpha, sep = "")
  col.fill
} # End of my.alpha.col().

my.alpha.append <- function(cols, from = "#333333", to = "#CCCCCC"){
  tl.level <- length(cols)
  alpha.fill <- colorRampPalette(c(from, to))(tl.level)
  alpha <- gsub(".....", "", alpha.fill)
  cols.fill <- paste(cols, alpha, sep = "")
  cols.fill
} # End of my.alpha.append().

my.YlOrRd <- function(level = 10){
  # col <- c("#FFFFCC", "#FFEDA0", "#FED976", "#FEB24C", "#FD8D3C", "#FC4E2A",
  #          "#E31A1C", "#BD0026", "#800026")
  col <- brewer.pal(9, "YlOrRd")
  col.fill <- NULL
  for(i.col in 1:(length(col) - 1)){
    tmp <- colorRampPalette(c(col[i.col], col[i.col + 1]))(level + 1)
    col.fill <- c(col.fill, tmp[-(level + 1)])
  }
  col.fill
} # End of my.YlOrRd().

my.YlGnBu <- function(level = 10){
  # col <- c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0",
  #          "#225EA8", "#253494", "#081D58")
  col <- brewer.pal(9, "YlGnBu")
  col.fill <- NULL
  for(i.col in 1:(length(col) - 1)){
    tmp <- colorRampPalette(c(col[i.col], col[i.col + 1]))(level + 1)
    col.fill <- c(col.fill, tmp[-(level + 1)])
  }
  col.fill
} # End of my.YlGnBu().

my.Reds <- function(level = 10){
  col <- brewer.pal(9, "Reds")
  col.fill <- NULL
  for(i.col in 1:(length(col) - 1)){
    tmp <- colorRampPalette(c(col[i.col], col[i.col + 1]))(level + 1)
    col.fill <- c(col.fill, tmp[-(level + 1)])
  }
  col.fill
} # End of my.Reds().

