runif.bcast <- function(n, min = 0, max = 1){
  ret <- as.double(runif(n, min, max))
  if(.MixfMRIEnv$MPI.gbd){
    ret <- pbdMPI::spmd.bcast.double(ret, double(1))
  }
  ret
} # End of runif.bcast().

rexp.bcast <- function(n, rate = 1){
  ret <- as.double(rexp(n, rate))
  if(.MixfMRIEnv$MPI.gbd){
    ret <- pbdMPI::spmd.bcast.double(ret, double(1))
  }
  ret
} # End of rexp.bcast().

length_gbd <- function(x.gbd){
  ret <- as.integer(length(x.gbd))
  if(.MixfMRIEnv$MPI.gbd){
    ret <- pbdMPI::spmd.allreduce.integer(ret, integer(1))
  }
  ret
} # End of length_gbd().

nrow.gbd <- function(x.gbd){
  ret <- as.integer(nrow(x.gbd))
  if(.MixfMRIEnv$MPI.gbd){
    ret <- pbdMPI::spmd.allreduce.integer(ret, integer(1))
  }
  ret
} # End of nrow.gbd().

sum_gbd <- function(x.gbd){
  ret <- as.double(sum(x.gbd))
  if(.MixfMRIEnv$MPI.gbd){
    ret <- pbdMPI::spmd.allreduce.double(ret, double(1))
  }
  ret
} # End of sum_gbd().

mean_gbd <- function(x.gbd){
  if(.MixfMRIEnv$MPI.gbd){
    n.gbd <- length(x.gbd) 
    n <- pbdMPI::spmd.allreduce.integer(n.gbd, integer(1))
    mean.x.gbd <- as.double(sum(x.gbd / n))
    ret <- pbdMPI::spmd.allreduce.double(mean.x.gbd, double(1))
  } else{
    ret <- as.double(mean(x.gbd))
  }
  ret
} # End of mean_gbd().

var.gbd.I <- function(X.gbd){
  p.X <- ncol(X.gbd)

  if(.MixfMRIEnv$MPI.gbd){
    n.gbd <- as.integer(nrow(X.gbd))
    n <- pbdMPI::spmd.allreduce.integer(n.gbd, integer(1))
    mean.X.gbd <- as.double(colSums(X.gbd / n))
    mean.X <- pbdMPI::spmd.allreduce.double(mean.X.gbd, double(p.X))

    tmp.X.gbd <- W.plus.y(X.gbd, -mean.X, nrow(X.gbd), p.X)
    var.X.gbd <- as.double(colSums(tmp.X.gbd^2) / (n - 1))
    ret <- pbdMPI::spmd.allreduce.double(var.X.gbd, double(p.X))
  } else{
    ret <- as.double(apply(X.gbd, 2, var))
  }
  ret
} # End of var.gbd.I().

var.gbd.V <- function(X.gbd){
  p.X <- ncol(X.gbd)

  if(.MixfMRIEnv$MPI.gbd){
    n.gbd <- as.integer(nrow(X.gbd))
    n <- pbdMPI::spmd.allreduce.integer(n.gbd, integer(1))
    mean.X.gbd <- as.double(colSums(X.gbd / n))
    mean.X <- pbdMPI::spmd.allreduce.double(mean.X.gbd, double(p.X))

    tmp.X.gbd <- W.plus.y(X.gbd, -mean.X, nrow(X.gbd), p.X)
    tmp.X.gbd <- tmp.X.gbd / sqrt(n - 1)
    var.X.gbd <- as.double(crossprod(tmp.X.gbd))
    ret <- pbdMPI::spmd.allreduce.double(var.X.gbd, double(length(var.X.gbd)))
  } else{
    ret <- as.double(var(X.gbd))
  }
  dim(ret) <- c(p.X, p.X)
  ret
} # End of var.gbd.V().

colSums.gbd <- function(x, na.rm = FALSE, dims = 1L){
  ret <- as.double(colSums(x, na.rm , dims))
  if(.MixfMRIEnv$MPI.gbd){
    ret <- pbdMPI::spmd.allreduce.double(ret, double(ncol(x)))
  }
  ret
} # End of colSums.gbd().

colMeans.gbd <- function(x, na.rm = FALSE, dims = 1L){
  if(.MixfMRIEnv$MPI.gbd){
    n <- nrow.gbd(x)
    ret <- colSums.gbd(matrix(x / n, ncol = ncol(x)), na.rm, dims)
  } else{
    ret <- colMeans(x, na.rm, dims)
  }
  ret
} # End of colSums.gbd().

tabulate.gbd <- function(bin, nbins = max(1, bin, na.rm = TRUE)){
  ret <- as.integer(tabulate(bin, nbins = nbins))
  if(.MixfMRIEnv$MPI.gbd){
    ret <- pbdMPI::spmd.allreduce.integer(ret, integer(nbins))
  }
  ret
} # End of tabulate().

max.gbd <- function(..., na.rm = FALSE){
  ret <- as.double(max(..., na.rm))
  if(.MixfMRIEnv$MPI.gbd){
    ret <- pbdMPI::spmd.allreduce.double(ret, double(1), op = "max")
  }
  ret
} # End of max.gbd().

min.gbd <- function(..., na.rm = FALSE){
  ret <- as.double(min(..., na.rm))
  if(.MixfMRIEnv$MPI.gbd){
    ret <- pbdMPI::spmd.allreduce.double(ret, double(1), op = "min")
  }
  ret
} # End of min.gbd().

crossprod.gbd <- function(x, y = NULL, ...){
  nrow <- ncol(x)
  ncol <- nrow
  if(!is.null(y)){
    ncol <- ncol(y)
  }
  ret <- as.double(crossprod(x, y = y, ...))
  if(.MixfMRIEnv$MPI.gbd){
    ret <- pbdMPI::spmd.allreduce.double(ret, double(length(ret)), op = "sum")
  }
  dim(ret) <- c(nrow, ncol)
  ret
} # End of crossprod.gbd().

