#### q-value
bh.fdr <- function(p){
  # This function computes q-values using Benjamini and Hochberg's (1995)
  # approach for controlling FDR.
  # Author: Dan Nettleton
  m <- length(p)
  # k <- 1:m
  ord <- order(p)
  p[ord] <- (p[ord] * m) / (1:m)
  qval <- p
  for(i in (m - 1):1) {
    qval[ord[i]] <- min(c(qval[ord[i]], qval[ord[i + 1]]))
  }
  return(qval)
} # End of bh.fdr().

by.fdr <- function(p){
  # This function computes q-values using Benjamini and Yekutieli's (2001)
  # approach for controlling FDR.
  # Author: Wei-Chen Chen 
  m <- length(p)
  hs <- sum(1 / (1:m))   # harmonic sum
  ord <- order(p)
  p[ord] <- (p[ord] * m * hs) / (1:m)
  qval <- p
  for(i in (m - 1):1) {
    qval[ord[i]] <- min(c(qval[ord[i]], qval[ord[i + 1]]))
  }
  qval[qval > 1] <- 1
  qval[qval < 0] <- 0
  return(qval)
} # End of by.fdr().

qvalue <- function(p, method = c("BH1995", "BY2001")){
  ret <- p
  id.ok <- (!is.na(ret)) & (ret >= 0) & (ret <= 1)

  if(method[1] == "BH1995"){
    ret[id.ok] <- bh.fdr(ret[id.ok])
  } else{
    ret[id.ok] <- by.fdr(ret[id.ok])
  }

  ret[!id.ok] <- NaN
  ret
} # End of qvalue().
