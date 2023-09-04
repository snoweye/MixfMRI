### Copy some functions from "AnalyzeFMRI/R/spatial.mixture.R" since
### "AnalyzeFMRI" is archived by CRAN.

cluster.threshold <- function(x, nmat = NULL, level.thr = 0.5, size.thr) {

  ## thresholds an array at level.thr
  ## calculates the number of contiguous clusters and their sizes
  ## answer is an array in which all voxels that are contained clusters of size greater
  ## than or equal to size.thr are 1, otherwise 0.
  ## nmat is a (Kx3) matrix specifying the neighbourhood system
  ## i.e if a row of nmat is (0, 1, -1) then x[10, 10, 10] and x[10, 11, 9] are neighbours

  if(is.null(nmat)) { ## default is 6 adjacent neighbours
    nmat <- expand.grid(-1:1, -1:1, -1:1)
    nmat <- nmat[c(5, 11, 13, 15, 17, 23), ]
  }   
  
  res <- .C("cluster_mass",
            mat = as.single(aperm(x, c(3, 2, 1))),
            as.integer(dim(x)),
            as.integer(t(nmat)),
            as.integer(dim(nmat)),
            as.single(level.thr),
            num.c = integer(1),
            res.c = single(1000 * 6),
            PACKAGE = "MixfMRI")

  res.c <- matrix(res$res.c, 1000, 6, byrow = TRUE)[1:res$num.c, ]

  mat1 <- array(res$mat, dim = dim(x)[3:1])
  mat1  <-  aperm(mat1, c(3, 2, 1))

  m <- (res.c[, 5] < size.thr) * (1:res$num.c)
  m <- m[m != 0]

  for(i in 1:length(m)) 
    mat1[mat1 == m[i]] <- 0
  
  mat1 <- 1 * (mat1 > 0)
  
  return(mat1)
}
