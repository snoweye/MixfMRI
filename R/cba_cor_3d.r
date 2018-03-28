### compute correlation
# dim(da.ts) = 46 x 55 x 43 x 105

cba.cor.3d <- function(da.ts, da.m = NULL, adj.dist = TRUE,
    fun.sim = stats::cor){
  dim.da.ts <- dim(da.ts)
  dim.da <- dim.da.ts[-length(dim.da.ts)]
  dim(da.ts) <- c(prod(dim.da.ts[-4]), dim.da.ts[4])

  ### Build inside brain ids
  id.inside <- 1:prod(dim.da)
  if(!is.null(da.m)){
    id.inside <- id.inside[da.m == 1]
  }

  ### Build a mask for neighbor
  mask.neighbor <- NULL
  for(i.3 in c(-1, 0, 1)){
    for(i.2 in c(-1, 0, 1)){
      for(i.1 in c(-1, 0, 1)){
        if(i.1 == 0 && i.2 == 0 && i.3 == 0){
          next
        }
        mask.neighbor <- cbind(mask.neighbor, c(i.1, i.2, i.3))
      }
    }
  }
  mask.neighbor.dist <- colSums(mask.neighbor^2)

  ### Set conversion functions
  loc2id <- function(x, dim.da){
    x[1,] + (x[2,] - 1) * dim.da[1] + (x[3,] - 1) * prod(dim.da[1:2])
  }
  id2loc <- function(id, dim.da){
    flag <- id %% dim.da[1] == 0
    id[flag] <- id[flag] - 1
    x.3 <- id %/% prod(dim.da[1:2]) + 1
    x.2 <- (id - (x.3 - 1) * prod(dim.da[1:2])) %/% dim.da[1] + 1
    x.1 <- (id - (x.3 - 1) * prod(dim.da[1:2]) - (x.2 - 1) * dim.da[1])
    x.1[flag] <- x.1[flag] + 1
    rbind(x.1, x.2, x.3)
  }

  ### Obtain the neighbor with max (adj) cor
  id.pair <- rep(0, prod(dim.da))
  for(i.0 in 1:prod(dim.da)){
    id.current <- id2loc(i.0, dim.da)
    id.neighbor <- mask.neighbor + as.vector(id.current)
    id.within.cube <- colSums(id.neighbor > 0 & id.neighbor <= dim.da) == 3
    id.neighbor <- id.neighbor[, id.within.cube]

    id.x <- i.0
    id.Y <- loc2id(id.neighbor, dim.da)
    x <- matrix(da.ts[id.x,], ncol = 1)
    Y <- matrix(da.ts[id.Y,], ncol = length(id.Y))
    tmp.cor <- fun.sim(x, Y)
    tmp.cor[!(id.Y %in% id.inside)] <- NA

    ### Adjust by distance if TRUE
    if(adj.dist){
      tmp.dist <- mask.neighbor.dist[id.within.cube]
      m.1 <- median(tmp.cor[tmp.dist == 1], na.rm = TRUE)
      m.2 <- median(tmp.cor[tmp.dist == 2], na.rm = TRUE)
      m.3 <- median(tmp.cor[tmp.dist == 3], na.rm = TRUE)
      tmp.cor[tmp.dist == 2] <- sqrt(tmp.cor[tmp.dist == 2] * m.1 * m.1 / m.2)
      tmp.cor[tmp.dist == 3] <- sqrt(tmp.cor[tmp.dist == 3] * m.1 * m.1 / m.3)
      tmp.cor[is.nan(tmp.cor)] <- NA
    }

    if(!all(is.na(tmp.cor))){
      id.pair[i.0] <- id.Y[which.max(tmp.cor)]
    }
  }

  ### Get clusters
  id.class <- rep(0, prod(dim.da))
  for(i.pair in 1:prod(dim.da)){
    if(id.pair[i.pair] == 0){
      next
    }
    tmp <- c(i.pair, id.pair[i.pair])
    tmp.class <- id.class[tmp]
    if(tmp.class[1] != 0 && tmp.class[2] != 0){
      id.class[id.class == max(tmp.class)] <- min(tmp.class)
    } else if(tmp.class[1] == 0 && tmp.class[2] != 0){
      id.class[tmp[1]] <- tmp.class[2]
    } else if(tmp.class[1] != 0 && tmp.class[2] == 0){
      id.class[tmp[2]] <- tmp.class[1]
    } else{
      id.class[tmp] <- max(id.class) + 1
    }
  }

  ### Trim cluster ids
  org.id <- sort(unique(id.class))
  if(length(org.id) < max(org.id)){
    if(0 %in% org.id){
      table.id <- cbind(org.id, c(0, 1:(length(org.id) - 1)))
    } else{
      table.id <- cbind(org.id, 1:length(org.id))
    }
    for(i in 1:nrow(table.id)){
      if(table.id[i, 1] != table.id[i, 2]){
        id.class[id.class == table.id[i, 1]] <- table.id[i, 2]
      }
    }
  }
  id.class[id.class == 0] <- NA

  ### Return
  ret <- id.class
  ret
} # End of cba.cor.3d().

# .rem <- function(){
#   dim <- c(46, 55, 43, 105)
#   set.seed(123)
#   da.ts <- array(rnorm(prod(dim)), dim = dim)
#   cba.cor(da.ts)
# }

