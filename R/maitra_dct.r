## from
## https://stackoverflow.com/questions/11215162/how-to-perform-a-fast-dct-discrete-cosine-transform-in-r
##
##  4 down vote favorite
##  Using Rprof revealed the dct in the dtt package to be the main offender in a piece of R code that was running quite slowly. 
##
## Take your vector and extend it to a vector twice as long as follows: If your vector is v=(1,2,3) then double the entries to w=(1,2,3,3,2,1). Note the ordering. If your vector is v=(1,2,4,9) then double the entries to w=(1,2,4,9,9,4,2,1)

## Let N be the length of your ORIGINAL vector (before you doubled its length).

## Then the first N coefficients of .5 * fft(w)/exp(complex(imaginary=pi / 2 / N)*(seq(2*N)-1)) should agree with computing dct(v) except it should be dramatically faster in almost all cases. 

## DCT <- function(x, inverse = FALSE) {
##    w <- c(x, rev(x))
##    N <- length(x)
##    Re(.5 * fft(w, inverse)/exp(complex(imaginary=pi / 2 / N)*(seq(2*N)-1)))[1:N]
##  }

## no longer needed


## compared with dct in dtt package

##y <- rnorm(50000)
##> system.time(DCT(y))
##   user  system elapsed 
##  0.029   0.002   0.031 
##> system.time(dct(y))
##   user  system elapsed 
##157.696   0.008 157.971 

## but DCT in fftw is faster


DCT2 <- function(x, inverse = FALSE, type=2)
  {
    # require(fftw)
    if (is.vector(x)) {
      if (inverse)
        IDCT(x, type = type)
      else 
        DCT(x, type = type)
    }
    else {
      y <- x
      if (inverse) {
        y <- t(apply(X = y, MARGIN = 1, FUN = IDCT, type = type))
        y <- apply(X = y, MARGIN = 2, FUN = IDCT, type = type)
      }
      else {
        y <- t(apply(X = y, MARGIN = 1, FUN = DCT, type = type))
        y <- apply(X = y, MARGIN = 2, FUN = DCT, type = type)
      } 
      y      
    }
  }

DCT3 <- function(x, inverse = FALSE, type = 2)
  {
    if (length(dim(x)) <=2) 
      DCT2(x, inverse = inverse, type = type)
    else {
      if (length(dim(x)) > 3)
        cat("DCT not implemented for arrays of dimension more than 4\n")
      else {
        y <- x
        y <- array(apply(X = y, MARGIN = 3, FUN = DCT2, inverse = inverse, type = type), dim = dim(x))
        y <- aperm(apply(X = y, MARGIN = c(1, 2), FUN = DCT2, inverse = inverse, type = type), perm = c(2, 3, 1))
        y
      } 
    }
  }
