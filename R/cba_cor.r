### compute correlation

cba.cor <- function(da.ts, da.m = NULL, adj.dist = TRUE,
    fun.sim = stats::cor){
  dim.da.ts <- dim(da.ts)

  if(!is.null(da.m)){
    dim.da.m <- dim(da.m)

    if(length(dim.da.ts) != (length(dim.da.m) + 1)){
      stop("Dimension lengths of da.ts and da.m do not match.")
    }
    if(any(dim.da.ts[-length(dim.da.ts)] != dim.da.m)){
      stop("Dimensions of da.ts and da.m are not consistent.")
    }
  }

  if(length(dim.da.ts) == 4){
    cba.cor.3d(da.ts, da.m = da.m, adj.dist = adj.dist, fun.sim = fun.sim)
  } else if(length(dim.da.ts) == 3){
    cba.cor.2d(da.ts, da.m = da.m, adj.dist = adj.dist, fun.sim = fun.sim)
  } else{
    stop("Dimension of da.ts should be for 2D or 3D.")
  }
} # End of cba.cor().

