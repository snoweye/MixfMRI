library(MixfMRI, quietly = TRUE)
set.seed(1234)

### Test maitra's simulation.
da <- gendataset(phantom = shepp1fMRI, overlap = 0.01)$pval
id <- !is.na(da)
PV.gbd <- da[id]
id.loc <- which(id, arr.ind = TRUE)
X.gbd <- t(t(id.loc) / dim(da))

### Run noX and 2 clusters.
.FC.CT$ignore.X <- TRUE
.FC.CT$INIT$class.method <- "extend"
PARAM <- fclust(X.gbd, PV.gbd, K = 2)

### Determine min.1st.prop and max.PV.
PARAM$param$ETA
quantile(PV.gbd, prob = PARAM$param$ETA[2], na.rm = TRUE)
quantile(PV.gbd, prob = c(0.03, 0.04, 0.05), na.rm = TRUE)
table(PV.gbd < 0.03) / sum(table(PV.gbd < 0.03))
