### This files provides functions for output.

mb.print <- function(PARAM, CHECK){
  ETA <- PARAM$ETA
  BETA <- matrix(do.call("c", PARAM$BETA), ncol = PARAM$K)
  if(PARAM$p.X > 0){
    MU <- PARAM$MU
    # SIGMA <- matrix(do.call("c",
    #                         lapply(PARAM$SIGMA,
    #                                function(x){
    #                                  if(is.matrix(x$full)){
    #                                    x$full[lower.tri(x$full, diag = TRUE)]
    #                                  } else{
    #                                    x$full
    #                                  })),
    #                 ncol = PARAM$K)
    SIGMA <- matrix(do.call("c", lapply(PARAM$SIGMA, function(x){ x$full })),
                    ncol = PARAM$K)
  }

  cat("\n")
  cat("Algorithm: ", CHECK$algorithm,
      "  Model.X: ", PARAM$model.X,
      "  Ignore.X: ", PARAM$ignore.X, "\n", sep = "")
  cat("Convergence: ", CHECK$convergence,
      "  iter: ", CHECK$iter,
      "  abs.err: ", CHECK$abs.err,
      "  rel.err: ", CHECK$rel.err, "\n", sep = "")
  cat("logL: ", PARAM$logL, "\n", sep = "")
  cat("K: ", PARAM$K, "\n", sep = "")
  cat("\nETA:\n")
  print(ETA)
  cat("\nBETA:\n")
  print(BETA)
  if(PARAM$p.X > 0){
    cat("\nMU:\n")
    print(MU)
    cat("\nSIGMA:\n")
    print(SIGMA)
  }
  cat("\n")

  invisible()
} # End of mb.print().
