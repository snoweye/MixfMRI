### Based on the proecdures 1, 2, and 3 in Benjamini & Heller (2007)
### `False discovery rates for spatial signals', JASA 102(480), 1272-1281
### Author: Wei-Chen Chen Oct. 24, 2014.

### Input:
###   p: p-values of length m.
###   w: weights for the p-values, default 1/m.
###   q: a desired cutoff for adjusting p-values.
### Output:
###   A list of two elements k and adjp.
###     k: the number of rejected small p-values.
###     adjp: the adjuated p-values.

### Procedure 1:
fdr.bh.p1 <- function(p, w = rep(1, length(p)), q = 0.05){
  ### Backup.
  ret.p <- p
  ret.w <- w
  id.ok <- (!is.na(ret.p)) & (ret.p >= 0) & (ret.p <= 1)
  p <- p[id.ok]
  w <- w[id.ok]

  ### Check again.
  if(any(is.na(p) | (p < 0) | (p > 1))){
    stop("p should have no NA and all are in [0, 1].")
  }

  m <- length(p)
  if(length(w) != m){
    stop("length(w) != m where m is the length of p.")
  }
  if(sum(w) != m){
    w <- w / sum(w) * m
  }

  order <- order(p)
  p.order <- p[order]
  w.order <- (cumsum(w) / m)
  k <- sum(p.order <= w.order * q)
  adjp <- p
  adjp[order] <- p.order / w.order

  ret <- list(k = k, adjp = adjp)

  tmp <- ret$adjp
  ret$adjp <- ret.p
  ret$adjp[id.ok] <- tmp
  ret$adjp[!id.ok] <- NaN

  ret
} # End of fdr.bh.p1()


### Procedure 2:
fdr.bh.p2 <- function(p, w = rep(1, length(p)), q = 0.05){
  ### Backup.
  ret.p <- p
  ret.w <- w
  id.ok <- (!is.na(ret.p)) & (ret.p >= 0) & (ret.p <= 1)
  p <- p[id.ok]
  w <- w[id.ok]

  ### Check again.
  if(any(is.na(p) | (p < 0) | (p > 1))){
    stop("p should have no NA and all are in [0, 1].")
  }

  m <- length(p)
  if(length(w) != m){
    stop("length(w) != m where m is the length of p.")
  }
  if(sum(w) != m){
    w <- w / sum(w) * m
  }

  ### Stage 1.
  q.new <- q / (q + 1)
  ret <- fdr.bh.p1(p, w = w, q = q.new)

  ### Stage 2.
  if(ret$k > 0){
    hat.m0.w <- m - sum(w[1:ret$k])
    order <- order(p)
    p.order <- p[order]
    w.order <- (cumsum(w) / hat.m0.w)
    k <- sum(p.order <= w.order * q.new)
    adjp <- p
    adjp[order] <- p.order / w.order
    ret <- list(k = k, adjp = adjp)
  }

  tmp <- ret$adjp
  ret$adjp <- ret.p
  ret$adjp[id.ok] <- tmp
  ret$adjp[!id.ok] <- NaN

  ret
} # End of fdr.bh.p2().
