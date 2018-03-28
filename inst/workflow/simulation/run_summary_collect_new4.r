ret.summary <- list()
ret.jaccard.summary <- list()
ret.alternative <- list()
ret.jaccard.alternative <- list()
for(i.rep in 1:25){
  prefix <- paste("simu_2d_", sprintf("%02d", i.rep), "/new/",
                  sep = "")

  fn <- paste(prefix, "summary_table/summary_table_new4.rda", sep = "")
  load(fn)
  ret.summary[[i.rep]] <- ret.table

  fn <- paste(prefix, "summary_table/summary_jaccard_table.rda", sep = "")
  load(fn)
  ret.jaccard.summary[[i.rep]] <- ret.table

  fn <- paste(prefix, "summary_alternative/summary_alternative.rda", sep = "")
  load(fn)
  ret.alternative[[i.rep]] <- ret.alt.table

  fn <- paste(prefix, "summary_alternative/summary_jaccard_alternative.rda", sep = "")
  load(fn)
  ret.jaccard.alternative[[i.rep]] <- ret.alt.table
}

save(ret.summary, ret.alternative, ret.jaccard.summary, ret.jaccard.alternative,
     file = "summary_I_0.975_new4.rda")
