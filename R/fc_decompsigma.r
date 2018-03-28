### This file contains functions for decomposition of SIGMA of MVN.
### These will be majorly used in m.step() to update PARAM$U and PARAM$U.check.

decompsigma.I <- function(full, sigma.ill = .FC.CT$CONTROL$sigma.ill){
  DS <- sqrt(full)
  DS.check <- TRUE
  if(min(DS) / max(DS) < sigma.ill){
    DS.check <- FALSE
  }

  list(value = DS, check = DS.check)
} # End of decompsigma.I().

decompsigma.V <- function(full, sigma.ill = .FC.CT$CONTROL$sigma.ill){
  DS <- try(chol(full))
  if(class(DS) == "try-error"){
    DS.check <- FALSE
  } else{
    ds <- diag(DS)
    DS.check <- TRUE
    if(any(ds < .MixfMRIEnv$CONTROL$DS.min |
           ds > .MixfMRIEnv$CONTROL$DS.max)){
      DS.check <- FALSE
    }
  }

  list(value = DS, check = DS.check)
} # End of decompsigma.V().

