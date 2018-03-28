##
## Garcia (CSDA, 2010) 
##
##

gcv.score <- function(s, Lambda, DCT2y, N)
  {
    Gamma <- 1/(1 + s * Lambda^2)
    RSS <- sum((DCT2y * (Gamma - 1))^2)
    TrH <- sum(Gamma);
    GCVs <- RSS / prod(N) / (1 - TrH/prod(N))^2
##    print(c(s,GCVs))
    GCVs
  }


gcv.smooth2d <- function(y, interval)
  {
    if (is.null(dim(y)))
      n <- length(y)
    else
      n <- dim(y)

    if (length(n) == 1) 
      lambda <- -2 + 2 * cos((0:(n - 1)) * pi / n)
    else 
      lambda <- kronecker( X = -2 + 2 * cos((0:(n[2] - 1)) * pi / n[2]),
                          Y = -2 + 2 * t(cos((0:(n[1] - 1)) * pi / n[1])), 
                          FUN = "+")
    
    dct2y <- DCT2(y, inverse = FALSE)
    
    par.val <- optimize(gcv.score, interval = interval, Lambda = lambda,
                     DCT2y = dct2y, N = n)

    shat <- par.val$minimum
    
    gamma <- 1/(1 + shat * lambda^2)

    c(list(im.smooth = DCT2(gamma * dct2y, inverse = TRUE), par.val = par.val))
  }
    
