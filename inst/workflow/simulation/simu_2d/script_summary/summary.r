options("width"=240)

phantom <- c("shepp0fMRI", "shepp1fMRI", "shepp2fMRI")
overlap <- c("0.01", "0.1", "0.25", "0.50", "0.75", "0.95")
min.1st.prop <- c("0.975")
max.PV <- c("0.05")
model <- c("I")

prefix <- "simu"
ret.summary <- NULL
ret.summary.org <- NULL
K.min <- 2
K.max <- 8
for(i.phantom in 1:length(phantom)){
  for(i.overlap in 1:length(overlap)){
    for(i.min.1st.prop in 1:length(min.1st.prop)){
      for(i.max.PV in 1:length(max.PV)){
        for(i.model in 1:length(model)){
          tmp.prefix <- c(phantom[i.phantom],
                          overlap[i.overlap],
                          min.1st.prop[i.min.1st.prop],
                          max.PV[i.max.PV],
                          model[i.model])
          file <- paste("./", prefix, "/",
                        paste(tmp.prefix, collapse = "_"),
                        "/output/ret.simu.2d.tp", sep = "")

          ### Load original first run.
          file.rda <- paste(file, ".rda", sep = "")
          if(! file.exists(file.rda)){
            tmp.summary <- rep(NA, 9)
            tmp.summary.org <- rep(NA, 9)
          } else{
            load(file.rda)
            tmp.ic <- lapply(ret.PARAM,
                        function(x){
                          if(class(x) == "try-error"){
                            tmp <- rep(NA, 4)
                          } else{
                            tmp <- c(x$param$K, x$aic, x$bic, x$icl.bic)
                          }
                          tmp
                        })
            tmp.ic <- do.call("rbind", tmp.ic)
            tmp.id.org <- apply(tmp.ic[, 2:4], 2, which.min)
            tmp.K.org <- tmp.ic[tmp.id.org, 1]

            ### Load new merged run.
            file.new.rda <- paste(file, ".new3.rda", sep = "")
            if(! file.exists(file.new.rda)){
              tmp.K.new <- rep(NA, 6)
              tmp.K.new.org <- rep(NA, 6)
            } else{
              load(file.new.rda)
              tmp.ic <- lapply(ret.PARAM,
                          function(x){
                            if(class(x) == "try-error"){
                              tmp <- rep(NA, 4)
                            } else{
                              tmp <- c(x$param$K, x$aic, x$bic, x$icl.bic)
                            }
                            tmp
                          })
              tmp.ic <- do.call("rbind", tmp.ic)
              tmp.id.new <- apply(tmp.ic[, 2:4], 2, which.min)
              tmp.K.new <- c(tmp.ic[tmp.id.new, 1], tmp.ic[tmp.id.org, 1])

              tmp.K.id <- NULL
              for(i.k in 1:length(ret.PARAM)){
                if(! is.null(ret.PARAM[[i.k]]) &&
                   class(ret.PARAM[[i.k]]) != "try-error"){
                  tmp.K.id <- c(tmp.K.id, i.k)
                }
              }
              tmp.K.new.org <- c(tmp.K.id[tmp.id.new], tmp.K.id[tmp.id.org])
            }

            tmp.summary <- c(tmp.K.org, tmp.K.new)
            tmp.summary.org <- c(tmp.K.org, tmp.K.new.org)
          }

          ### ret.summary is the best K for reporting.
          ### ret.summary.org is the original run where the best K comes from.
          ### Use ret.summary.org to get right run for ploting and adjR.
          ret.summary <- rbind(ret.summary, c(tmp.prefix, tmp.summary))
          ret.summary.org <- rbind(ret.summary.org,
                                   c(tmp.prefix, tmp.summary.org))
        }
      }
    }
  }
}
colnames(ret.summary) <-
  c("phantom", "overlap", "min.1st.prop", "max.PV", "model",
    "K.AIC", "K.BIC", "K.ICL.BIC",
    "K.AIC.merge", "K.BIC.merge", "K.ICL.BIC.merge",
    "K.AIC.maitra", "K.BIC.maitra", "K.ICL.BIC.maitra")
cat("K:\n")
print(ret.summary)
colnames(ret.summary.org) <-
  c("phantom", "overlap", "min.1st.prop", "max.PV", "model",
    "K.AIC", "K.BIC", "K.ICL.BIC",
    "K.AIC.merge", "K.BIC.merge", "K.ICL.BIC.merge",
    "K.AIC.maitra", "K.BIC.maitra", "K.ICL.BIC.maitra")
cat("K.org:\n")
print(ret.summary.org)

save(ret.summary, ret.summary.org, file = "new/summary/summary.rda")
