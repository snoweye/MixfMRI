summarized.overlap <- function(overlap.mat){
  x <- eigen(overlap.mat, symmetric = T, only.values = T)
  p <- nrow(overlap.mat)
  if (p == 1) 1 else
  max(min((p * x$values[1] / sum(x$values) - 1) / (p - 1), 1), 0)
} # End of summarized.overlap().
