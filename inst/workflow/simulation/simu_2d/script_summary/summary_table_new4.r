library(MixfMRI, quietly = TRUE)
prefix <- "simu"
load("new/summary/summary.rda")

### Plot
ret.table <- NULL
for(i.case in 1:nrow(ret.summary.org)){
  tmp.phantom <- eval(parse(text = ret.summary.org[i.case, 1]))
  tmp.overlap <- eval(parse(text = ret.summary.org[i.case, 2]))
  set.seed(1235)
  da <- gendataset(phantom = tmp.phantom, overlap = tmp.overlap)$pval

  tmp.prefix <- ret.summary.org[i.case, 1:5]

  file <- paste("./", prefix, "/",
                paste(tmp.prefix, collapse = "_"),
                "/output/ret.simu.2d.tp", sep = "")

  file.rda <- paste(file, ".rda", sep = "")
  if(file.exists(file.rda)){
    load(file.rda)
    ret.PARAM.org <- ret.PARAM
  } else{
    ret.PARAM.org <- NULL
  }

  file.rda <- paste(file, ".new4.rda", sep = "")
  if(file.exists(file.rda)){
    load(file.rda)
    ret.PARAM.new <- ret.PARAM
  } else{
    ret.PARAM.new <- NULL
  }

  tmp <- ret.summary[i.case, c(1:8, 12:14)]
  tmp.class <- tmp.phantom[!is.na(tmp.phantom)] + 1
  tmp.class[tmp.class != 1] <- 2

  ### K.AIC.maitra
  i.k.org <- as.integer(ret.summary.org[i.case, 12])
  if(is.na(i.k.org)){
    tmp <- c(tmp, NA)
  } else{
    new.class <- ret.PARAM.new[[i.k.org]]$class
    new.class[new.class != 1] <- 2
    i.R <- RRand(new.class, tmp.class)$adjRand
    tmp <- c(tmp, i.R)
  }

  ### K.BIC.maitra
  i.k.org <- as.integer(ret.summary.org[i.case, 13])
  if(is.na(i.k.org)){
    tmp <- c(tmp, NA)
  } else{
    new.class <- ret.PARAM.new[[i.k.org]]$class
    new.class[new.class != 1] <- 2
    i.R <- RRand(new.class, tmp.class)$adjRand
    tmp <- c(tmp, i.R)
  }

  ### K.ICL.BIC.maitra
  i.k.org <- as.integer(ret.summary.org[i.case, 14])
  if(is.na(i.k.org)){
    tmp <- c(tmp, NA)
  } else{
    new.class <- ret.PARAM.new[[i.k.org]]$class
    new.class[new.class != 1] <- 2
    i.R <- RRand(new.class, tmp.class)$adjRand
    tmp <- c(tmp, i.R)
  }

  ret.table <- rbind(ret.table, tmp)
}

colnames(ret.table)[12:14] <-
   c("AIC.maitra.R", "BIC.maitra.R", "ICL.BIC.maitra.R")
rownames(ret.table) <- NULL
print(ret.table)

save(ret.table, file = "new/summary_table/summary_table_new4.rda")
