# Compute AIC and BIC for given PARAM.

aic <- function(PARAM){
  d <- compute.d(PARAM)
  -2 * PARAM$logL + 2 * d
} # End of aic().

bic <- function(PARAM){
  d <- compute.d(PARAM)
  -2 * PARAM$logL + d * log(PARAM$N)
} # End of bic().

icl.bic <- function(PARAM){
  d <- compute.d(PARAM)
  -2 * PARAM$logL + d * log(PARAM$N) + 2 * PARAM$entropy
} # End of icl.bic().

compute.d <- function(PARAM){
  if(PARAM$model.X == "I"){
    d.normal <- PARAM$p.X * 2
  } else if(PARAM$model.X == "V"){
    d.normal <- PARAM$p.X + (PARAM$p.X + 1) * PARAM$p.X / 2
  } else{
    stop("model.X is not implemented.")
  }

  d <- (PARAM$K - 1) +                       # ETA
       d.normal +                            # The first cluster
       (2 + d.normal) * (PARAM$K - 1)        # The rest clusters.

  if(PARAM$ETA[1] == PARAM$min.1st.prop){
    d <- d - 1
  }

  d
} # End of compute.d()
