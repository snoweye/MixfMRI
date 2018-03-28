### This file gives initializations.

# x <- rbeta(1000, 0.5, 2)
# (m <- mean(x))
# (v <- var(x))
# (alpha <- 1/v * (m^2 - m^3 - m * v))
# (beta <- (1 - m) / m * alpha)
# alpha / (alpha + beta)
# alpha * beta / (alpha + beta)^2 / (alpha + beta + 1)

initial.em.gbd <- function(PARAM){
  if(PARAM$class.method == "prob.extend"){
    PARAM <- initial.em.gbd.prob.extend(PARAM)
  } else if(PARAM$class.method == "prob.simple"){
    PARAM <- initial.em.gbd.prob.simple(PARAM)
  } else if(PARAM$class.method == "qnorm.extend"){
    PARAM <- initial.em.gbd.qnorm.extend(PARAM)
  } else if(PARAM$class.method == "qnorm.simple"){
    PARAM <- initial.em.gbd.qnorm.simple(PARAM)
  } else if(PARAM$class.method == "extend"){
    PARAM <- initial.em.gbd.extend(PARAM)
  } else if(PARAM$class.method == "simple"){
    PARAM <- initial.em.gbd.simple(PARAM)
  } else{
    stop("class.method is not found.")
  }

  PARAM
} # End of initial.em.gbd().


### initial.em.gbd.update called by the last step of each initial method to
### update the new parameters according to new initial classifications.
initial.em.gbd.update <- function(PARAM, CLASS.gbd, K){
  ### Update parameters.
  PARAM$BETA <- initial.BETA(CLASS.gbd, K)
  if(PARAM$p.X > 0){
    PARAM$MU <- initial.MU(CLASS.gbd, K)
    PARAM$SIGMA <- .MixfMRIEnv$initial.SIGMA(PARAM$MU, CLASS.gbd, K)
  }

  ### Do one step em to update Z and logL.
  e.step.gbd(PARAM)
  PARAM <- em.onestep.gbd(PARAM)
  PARAM$logL <- logL.step.gbd()
  em.update.class.gbd()

  PARAM
} # End of initial.em.gbd.update().

