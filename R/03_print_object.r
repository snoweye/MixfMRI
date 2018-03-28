### This files provides functions for output.

print.fclust <- function(x, ...){
  CHECK <- x$check
  PARAM <- x$param
  N.CLASS <- x$n.class
  Ignore.X <- x$ignore.X
  AIC <- x$aic
  BIC <- x$bic
  ICL.BIC <- x$icl.bic

  ETA <- PARAM$ETA
  BETA <- do.call("cbind", PARAM$BETA)
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

  cat("Algorithm: ", CHECK$algorithm,
      "  Model.X: ", PARAM$model.X,
      "  Ignore.X: ", Ignore.X, "\n",
      "- Convergence: ", CHECK$convergence,
      "  iter: ", CHECK$iter,
      "  abs.err: ", CHECK$abs.err,
      "  rel.err: ", CHECK$rel.err, "\n", sep = "")
  cat("- N: ", PARAM$N, "  p.X: ", PARAM$p.X, "  K: ", PARAM$K,
      "  logL: ", PARAM$logL, "\n",
      "- AIC: ", AIC, "  BIC: ", BIC, " ICL-BIC: ", ICL.BIC,
      "\n", sep = "")
  cat("- n.class: ", paste(N.CLASS, collapse = " "), "\n",
      sep = "")
  cat("- init.class.method: ", PARAM$init.class.method, "\n", sep = "")
  cat("- ETA: (min.1st.prop: ", PARAM$min.1st.prop,
      "  max.PV: ", PARAM$max.PV, ")\n", sep = "")
  print(ETA)
  cat("- BETA: (2 by K)\n")
  print(BETA)
  if(PARAM$p.X > 0){
    cat("- MU: (p.X by K)\n")
    print(MU)
    cat("- SIGMA: (d.normal by K)\n")
    print(SIGMA)
  }
  cat("\n")

  invisible()
} # End of print.fclust().
